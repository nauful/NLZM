#pragma once

#include "platform.h"

/* Notes:
1. NStates & NMasked contain number of symbols minus 1
2. contexts example:
 MaxOrder:
  ABCD    context
   BCD    suffix
   BCDE   successor
 lower orders:
   BCD    context
    CD    suffix
   BCDE   successor
3. Code notation:
 0 - binary PPM-context
 1 - non-binary PPM-context without masked symbols
 2 - non-binary PPM-context with masked symbols
 LES - last encountered symbol in context

PPMd Memory Map:
{
  [ 0 ]           contains subset of original raw text, that is required to create context
                  records, Some symbols are not written, when max order context was reached
  [ Text ]        free area
  [ UnitsStart ]  Ppmd_State vectors and Ppmd7_Context records
  [ LoUnit ]      free area for Ppmd_State and Ppmd7_Context items
  [ HiUnit ]      Ppmd7_Context records
  [ Size ]        end of array
}

These addresses don't cross at any time.
And the following conditions is true for addresses:
  (0  <= Text < UnitsStart <= LoUnit <= HiUnit <= Size)

Raw text is BYTE--aligned.
the data in block [ UnitsStart ... Size ] contains 12-bytes aligned UNITs.

Last UNIT of array at offset (Size - 12) is root order-0 Ppmd7_Context record.
The code can free UNITs memory blocks that were allocated to store Ppmd_State vectors.
The code doesn't free UNITs allocated for Ppmd7_Context records.

The code calls Ppmd7_RestartModel(), when there is no free memory for allocation.
And Ppmd7_RestartModel() changes the state to start state, with full free block.

The code allocates UNITs with the following order:

Allocation of 1 UNIT for Context record
  - from free space (HiUnit) down to (LoUnit)
  - from FreeList[0]
  - Ppmd7_AllocUnitsRare()

Ppmd7_AllocUnits() for Ppmd_State vectors:
  - from FreeList[i]
  - from free space (LoUnit) up to (HiUnit)
  - Ppmd7_AllocUnitsRare()

Ppmd7_AllocUnitsRare()
  - if (GlueCount == 0)
       {  Glue lists, GlueCount = 255, allocate from FreeList[i]] }
  - loop for all higher sized FreeList[...] lists
  - from (UnitsStart - Text), GlueCount--
  - ERROR

Each Record with Context contains the Ppmd_State vector, where each
Ppmd_State contains the link to Successor.
There are 3 types of Successor:
  1) NULL-Successor   - NULL pointer. NULL-Successor links can be stored
                        only in 0-order Root Context Record.
                        We use 0 value as NULL-Successor
  2) RAW-Successor    - the link to position in raw text,
                        that "RAW-Successor" is being created after first
                        occurrence of new symbol for some existing context record.
                        (RAW-Successor > 0).
  3) RECORD-Successor - the link to Ppmd7_Context record of (Order+1),
                        that record is being created when we go via RAW-Successor again.

For any successors at any time: the following condtions are true for Successor links:
(NULL-Successor < RAW-Successor < UnitsStart <= RECORD-Successor)

---------- Symbol Frequency, SummFreq and Range in Range_Coder ----------

Ppmd7_Context::SummFreq = Sum(Stats[].Freq) + Escape_Freq

The PPMd code tries to fulfill the condition:
  (SummFreq <= (256 * 128 = RC::kBot))

We have (Sum(Stats[].Freq) <= 256 * 124), because of (MAX_FREQ = 124)
So (4 = 128 - 124) is average reserve for Escape_Freq for each symbol.
If (Ppmd_State::Freq) is not aligned for 4, the reserve can be 5, 6 or 7.
SummFreq and Escape_Freq can be changed in Ppmd7_Rescale() and *Update*() functions.
Ppmd7_Rescale() can remove symbols only from max-order contexts. So Escape_Freq can increase after multiple calls of Ppmd7_Rescale() for
max-order context.

When the PPMd code still break (Total <= RC::Range) condition in range coder,
we have two ways to resolve that problem:
  1) we can report error, if we want to keep compatibility with original PPMd code that has no fix for such cases.
  2) we can reduce (Total) value to (RC::Range) by reducing (Escape_Freq) part of (Total) value.
*/

#pragma pack(push, 1)
namespace PPMD {
  template<class T>
  T Clamp(const T &v, const T &lo, const T &hi) {
    return (v >= lo) ? ((v <= hi) ? v : hi) : lo;
  }

  template<class T>
  T Max(T a, T b) {
    return a >= b ? a : b;
  }

  template<class T>
  T Min(T a, T b) {
    return a <= b ? a : b;
  }

  template<class T>
  void Swap(T &a, T &b) {
    T tmp = a;
    a = b;
    b = tmp;
  }

  struct RangeCoder {
    enum {
      TOP = 1 << 24,
      BOT = 1 << 15
    };

    u32 low;
    u32 code;
    u32 range, range_bin;

    RangeCoder() {
      low = 0;
      code = 0;
      range = range_bin = 0;
    }

    void InitEncoder();
    void InitDecoder(const u8 **stream);
    void FlushEncoder(u8 **stream);
    void NormalizeEncoder(u8 **stream);
    void NormalizeDecoder(const u8 **stream);

    void Arrange(u32 Scale);
    void Interval(u32 rlow, u32 high);
    u32 GetCurrentCount();
    u32 BinStart(u32 f0, u32 shift);
    bool BinDecode(u32 tmp);
    void BinCorrect0(u32 tmp);
    void BinCorrect1(u32 tmp);
  };

  enum {
    UNIT_SIZE = 12,
    N_INDEXES = 128,

    UP_FREQ = 5,
    BIN_INT_BITS = 7,
    BIN_PERIOD_BITS = 7,
    BIN_TOT_BITS = BIN_INT_BITS + BIN_PERIOD_BITS,
    BIN_INTERVAL = 1 << BIN_INT_BITS,
    BIN_SCALE = 1 << BIN_TOT_BITS,
    BIN_ROUND = 16,
    MAX_FREQ = 124,
    O_BOUND = 8,
    MAX_ORDER = 12
  };

  struct Context;
  struct State;
  struct Model;

  struct Alloc {
    u32 size, num_units_available;

    u8 *heap_start;
    u8 *text_area, *units_start;
    u8 *units_lo, *units_hi;

    u32 PtrToIndex(void *p) const {
      PLATFORM_ASSERT_DEBUG(p >= heap_start);
      PLATFORM_ASSERT_DEBUG(p < heap_start + size);

      return (u32)((u8 *)p - heap_start);
    }

    u8 *IndexToPtr(u32 i) {
      PLATFORM_ASSERT_DEBUG(i < size);

      return heap_start + i;
    }

    struct BlockNode {
      u32 header, idx_next, idx_prev;

      // Free stamp overlaps with
      // State::idx_successor
      // Context::idx_suffix
      // Indices are limited to 0x7FFFFFFF, or < 2G.

      bool GetFreeStamp() const { return (header & 0x80000000u) != 0; }
      u32 GetBlockNumUnits() const { return header & 0x7FFFFFFFu; }

      void SetFreeStamp(u32 stamp) {
        header = (stamp << 31) | (header & 0x7FFFFFFFu);
      }

      void SetBlockNumUnits(u32 num_units) {
        PLATFORM_ASSERT_DEBUG(num_units < 0x80000000u);
        header = (header & 0x80000000u) | num_units;
      }

      BlockNode *GetNext(Alloc *alloc) const { return alloc->IndexToNode(idx_next); }
      bool HasNext() const { return idx_next != 0; }

      BlockNode *GetPrev(Alloc *alloc) const { return alloc->IndexToNode(idx_prev); }
      bool HasPrev() const { return idx_prev != 0; }

      void Link(Alloc *alloc, BlockNode *p);
      void Unlink(Alloc *alloc);

      void Insert(Alloc *alloc, BlockNode *p, int num_units);
      BlockNode *Remove(Alloc *alloc);
    };

    Alloc();

    Context *IndexToContext(u32 i) { return (Context *)IndexToPtr(i); }
    State *IndexToState(u32 i) { return (State *)IndexToPtr(i); }
    BlockNode *IndexToNode(u32 i) { return (BlockNode *)IndexToPtr(i); }

    bool Start(u32 suballoc_size);
    void Init();
    void Stop();

    u32 GetUsedMemory();
    
    State *AllocUnits(u32 num_units);
    void FreeUnits(void *ptr, u32 num_units);

    State *ExpandUnits(State *old_ptr, u32 num_units);
    State *ShrinkUnits(State *old_ptr, u32 old_num_units, u32 new_num_units);
    State *MoveUnitsUp(State *old_ptr, u32 num_units);

    Context *AllocContext();
    void FreeContext(Context *ptr);

  protected:
    BlockNode *block_list;

    void SplitBlock(BlockNode *pv, u32 old_index, u32 new_index);
    void *AllocUnitsRare(u32 indx);
  };

  // SEE-contexts for PPM-contexts with masked symbols
  struct SEE {
    u16 acc;
    u8 shift;
    u8 count;

    void Init(u32 v0);
    u32 Mean() const { return acc >> shift; }
    void UpdateSymbol();
    void SetShift();
  };

  struct State {
    u32 idx_successor;
    u8 symbol, freq;

    Context *GetSuccessor(Alloc *alloc) const;
  };

  struct Context {
    u32 idx_suffix;
    u8 num_stats, flags;

    union {
      struct {
        u16 sum_freq;
        u32 idx_stats;
      } transition;

      State one;
    };

    State *GetStats(Alloc *alloc) const;
    Context *GetSuffix(Alloc *alloc) const;
  };

  struct Model {
    Context *max_context;
    State *found_state;

    u8 ns_to_bs_index[256], tbl_q[260];
    int bsum, order_fall, run_len;
    int run_len_init, order_max;
    u8 prev_success_see;

    SEE tbl_see[23][32];
    u16 tbl_bin_see[25][64];

    u8 cutoff_mode;

    u8 escape_mask[256 >> 3];
    u8 num_masked;

    static int MaskIndex(int c) { return c >> 3; }
    static int MaskBit(int c) { return 1 << (c & 7); }
    bool EscMaskIsSet(int c) { return escape_mask[MaskIndex(c)] & MaskBit(c); }
    void EscMaskSet(int c) { escape_mask[MaskIndex(c)] |= MaskBit(c); }
    void EscMaskClear() { memset(escape_mask, 0, sizeof(escape_mask)); }
  };

  void Init(Model *model, Alloc *alloc);
  void StartModel(Model *model, Alloc *alloc, int order_max, int cutoff_mode);
  void RestoreModel(Model *model, Alloc *alloc, Context *pc);
  void AuxCutOff(Model *model, Alloc *alloc, State *p, int order);

  u32 CreateSuccessors(Alloc *alloc, State *fs, bool skip, State *p, Context *pc);
  u32 ReduceOrder(Model *model, Alloc *alloc, State *p, Context *pc);
  void UpdateModel(Model *model, Alloc *alloc, Context *min_ctx);
  void Rescale(Model *model, Alloc *alloc, Context *ctx);
  u32 CutOff(Model *model, Alloc *alloc, Context *ctx, int order);

  PPMD::SEE* GetSEE(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx, int num_stats, int flags);
  u16 *GetBinSEE(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx);

  void EncodeFile(PPMD::Model *model, PPMD::Alloc *alloc, u8 **DecodedFile, u8 **EncodedFile, int order_max, u8 cutoff_mode);
  void DecodeFile(PPMD::Model *model, PPMD::Alloc *alloc, u8 **DecodedFile, const u8 **EncodedFile, int order_max, u8 cutoff_mode);
}
#pragma pack(pop)
