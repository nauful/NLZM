#include "ppmd.h"

u32 UnitsToBytes(u32 num_units) {
  return PPMD::UNIT_SIZE * num_units;
}

void CopyUnits(void *dst, void *src, u32 num_units) {
  memcpy(dst, src, UnitsToBytes(num_units));
}

PPMD::Alloc::Alloc() {
  size = 0;

  heap_start = 0;
  text_area = 0;
  units_start = 0;
  units_lo = 0;
  units_hi = 0;

  block_list = 0;
}

bool PPMD::Alloc::Start(u32 suballoc_size) {
  const u32 t = suballoc_size << 20U;

  heap_start = new u8[t];
  if (!heap_start) {
    return false;
  }

  memset(heap_start, 0, t);
  size = t;
  return true;
}

void PPMD::Alloc::BlockNode::Link(PPMD::Alloc *alloc, PPMD::Alloc::BlockNode *p) {
  u32 i = alloc->PtrToIndex(p);
  p->SetFreeStamp(1);

  if (HasNext()) {
    GetNext(alloc)->idx_prev = i;
  }

  p->idx_next = idx_next;
  p->idx_prev = alloc->PtrToIndex(this);

  idx_next = i;
}

void PPMD::Alloc::BlockNode::Unlink(PPMD::Alloc *alloc) {
  if (HasNext()) {
    GetNext(alloc)->idx_prev = idx_prev;
  }

  if (HasPrev()) {
    GetPrev(alloc)->idx_next = idx_next;
  }

  header = 0;
  idx_next = idx_prev = 0;
}

PPMD::Alloc::BlockNode *PPMD::Alloc::BlockNode::Remove(PPMD::Alloc *alloc) {
  BlockNode *p = GetNext(alloc);
  PLATFORM_ASSERT_DEBUG(p->GetFreeStamp());

  p->Unlink(alloc);

  return p;
}

void PPMD::Alloc::BlockNode::Insert(PPMD::Alloc *alloc, BlockNode *p, int num_units) {
  Link(alloc, p);

  p->SetBlockNumUnits(num_units);
}

void PPMD::Alloc::Init() {
  u32 diff = (7 * size) / 8;
  diff -= diff % UNIT_SIZE;

  text_area = heap_start;
  units_hi = heap_start + size;
  units_lo = units_hi - diff;
  num_units_available = 0;

  units_start = units_lo;

  block_list = (BlockNode *)(units_hi -= N_INDEXES * sizeof(BlockNode));
  memset(block_list, 0, N_INDEXES * sizeof(BlockNode));

#if 0
  u32 sz = (units_hi - units_lo) / UNIT_SIZE;
  BlockNode *p = (BlockNode *)units_lo;
  units_lo += UnitsToBytes(sz);

  for (; sz > 128; sz -= 128, p += 128) {
   block_list[N_INDEXES - 1].Insert(this, p, 128);
  }

  block_list[sz - 1].Insert(this, p, sz);
#endif
}

u32 PPMD::Alloc::GetUsedMemory() {
  u32 ret = size - (units_hi - units_lo) - (units_start - text_area) - (num_units_available * UNIT_SIZE);

#if 0
  u32 verify_available = 0;
  for (int index = 0; index < N_INDEXES; index++) {
    for (BlockNode *p = block_list[index].HasNext() ? block_list[index].GetNext(this) : 0;
      p; p = p->HasNext() ? p->GetNext(this) : 0) {
      verify_available += p->GetBlockNumUnits();
    }
  }
  PLATFORM_ASSERT(verify_available == num_units_available);
#endif

  return ret;
}

void PPMD::Alloc::Stop() {
  if (size) {
    printf("StopSubAllocator %d KB\n", GetUsedMemory() >> 10);

    size = 0;
    delete[] heap_start;
  }
}

void PPMD::Alloc::SplitBlock(BlockNode *pv, u32 old_index, u32 new_index) {
  u32 unit_diff = old_index - new_index;
  BlockNode *p = pv + (new_index + 1);

  num_units_available += unit_diff;
  block_list[unit_diff - 1].Insert(this, p, unit_diff);
}

void *PPMD::Alloc::AllocUnitsRare(u32 index) {
  u32 i = index;

  while (++i < N_INDEXES) {
    if (block_list[i].HasNext()) {
      BlockNode *ret = block_list[i].Remove(this);
      num_units_available -= i + 1;

      SplitBlock(ret, i, index);

      return ret;
    }
  }

  // Shrink text area
  u32 ub = UnitsToBytes(index + 1);
  if (text_area + ub < units_start) {
    units_start -= ub;
    return units_start;
  }
  else {
    return nullptr;
  }
}

PPMD::State *PPMD::Alloc::AllocUnits(u32 num_units) {
  u32 index = num_units - 1;

  if (block_list[index].HasNext()) {
    num_units_available -= num_units;
    return (State *)block_list[index].Remove(this);
  }

  void *ret = units_lo;
  int ub = UnitsToBytes(index + 1);
  if (units_lo + ub <= units_hi) {
    units_lo += ub;
    return (State *)ret;
  }

  return (State *)AllocUnitsRare(index);
}

PPMD::Context *PPMD::Alloc::AllocContext() {
  if (units_lo < units_hi) {
    units_hi -= UNIT_SIZE;
    return (Context *)units_hi;
  }

  if (block_list[0].HasNext()) {
    num_units_available -= 1;
    return (Context *)block_list[0].Remove(this);
  }

  return (Context *)AllocUnitsRare(0);
}

void PPMD::Alloc::FreeUnits(void *ptr, u32 num_units) {
  BlockNode *p = (BlockNode *)ptr;

  u32 sz = num_units;
  while (p + sz < block_list && p[sz].GetFreeStamp()) {
    BlockNode *p1 = p + sz;
    sz += p1->GetBlockNumUnits();

    p1->Unlink(this);
  }

  if (ptr > units_start) {
    num_units_available += num_units;

    for (; sz > 128; sz -= 128, p += 128) {
      block_list[N_INDEXES - 1].Insert(this, p, 128);
    }

    block_list[sz - 1].Insert(this, p, sz);
  }
  else {
    // Expand text area
    num_units_available -= sz - num_units;
    units_start += UnitsToBytes(sz);
  }
}

void PPMD::Alloc::FreeContext(PPMD::Context *ptr) {
  FreeUnits(ptr, 1);
}

PPMD::State *PPMD::Alloc::ExpandUnits(PPMD::State *old_ptr, u32 num_units) {
  void *ptr = AllocUnits(num_units + 1);

  if (ptr) {
    CopyUnits(ptr, old_ptr, num_units);
    FreeUnits(old_ptr, num_units);
  }

  return (State *)ptr;
}

PPMD::State *PPMD::Alloc::ShrinkUnits(PPMD::State *old_ptr, u32 old_num_units, u32 new_num_units) {
  u32 i0 = old_num_units - 1;
  u32 i1 = new_num_units - 1;

  if (i0 == i1) {
    return old_ptr;
  }

  if (block_list[i1].HasNext()) {
    num_units_available -= new_num_units;
    void *ptr = block_list[i1].Remove(this);

    CopyUnits(ptr, old_ptr, new_num_units);
    FreeUnits(old_ptr, old_num_units);
    return (State *)ptr;
  }
  else {
    SplitBlock((BlockNode *)old_ptr, i0, i1);
    return (State *)old_ptr;
  }
}

PPMD::State *PPMD::Alloc::MoveUnitsUp(PPMD::State *old_ptr, u32 num_units) {
  u32 index = num_units - 1;

  if (!block_list[index].HasNext() || old_ptr > (void *)block_list[index].GetNext(this)) {
    return old_ptr;
  }
  else {
    num_units_available -= num_units;

    void *ptr = block_list[index].Remove(this);
    CopyUnits(ptr, old_ptr, num_units);
    FreeUnits(old_ptr, num_units);

    return (State *)ptr;
  }
}
