#include "ppmd.h"

void PPMD::SEE::Init(u32 v0) {
  shift = BIN_PERIOD_BITS - 4;
  acc = v0 << shift;
  count = 7;
}

void PPMD::SEE::UpdateSymbol() {
  if (!--count) {
    SetShift();
  }
}

void PPMD::SEE::SetShift() {
  u32 i = acc >> shift;
  i = BIN_PERIOD_BITS - (i > 40) - (i > 280) - (i > 1020);

  if (i < shift) {
    acc >>= 1;
    shift--;
  }
  else if (i > shift) {
    acc <<= 1;
    shift++;
  }

  count = 5 << shift;
}

PPMD::Context *PPMD::State::GetSuccessor(PPMD::Alloc *alloc) const {
  return alloc->IndexToContext(idx_successor);
}

PPMD::State *PPMD::Context::GetStats(PPMD::Alloc *alloc) const {
  return alloc->IndexToState(transition.idx_stats);
}

PPMD::Context *PPMD::Context::GetSuffix(PPMD::Alloc *alloc) const {
  return alloc->IndexToContext(idx_suffix);
}

void PPMD::Init(PPMD::Model *model, PPMD::Alloc *alloc) {
  model->max_context = nullptr;

  int i, k, m, s;

  model->ns_to_bs_index[0] = 2 * 0;
  model->ns_to_bs_index[1] = 2 * 1;
  model->ns_to_bs_index[2] = 2 * 1;
  memset(model->ns_to_bs_index + 3, 2 * 2, 26);
  memset(model->ns_to_bs_index + 29, 2 * 3, 256 - 29);

  for (i = 0; i < UP_FREQ; i++) {
    model->tbl_q[i] = i;
  }

  for (m = i = UP_FREQ, k = s = 1; i < 260; i++) {
    model->tbl_q[i] = m;

    if (!--k) {
      k = ++s;
      m++;
    }
  }
}

void PPMD::StartModel(PPMD::Model *model, PPMD::Alloc *alloc, int order_max, int cutoff_mode) {
  int i, k, s;
  u8 i2f[25];

  alloc->Init();

  model->EscMaskClear();

  if (cutoff_mode >= 0) {
    model->cutoff_mode = cutoff_mode;
  }

  model->order_fall = model->order_max = order_max;
  model->run_len = model->run_len_init = -((order_max < 13) ? order_max : 13);
  model->max_context = (Context *)alloc->AllocContext();
  model->max_context->num_stats = 255;
  model->max_context->transition.sum_freq = model->max_context->num_stats + 2;
  model->max_context->transition.idx_stats = alloc->PtrToIndex(alloc->AllocUnits(256 >> 1));

  model->prev_success_see = 0;
  model->max_context->idx_suffix = 0;
  model->max_context->flags = 0;

  State *stats = model->max_context->GetStats(alloc);
  for (i = 0; i < 256; i++) {
    stats[i].symbol = i;
    stats[i].freq = 1;
    stats[i].idx_successor = 0;
  }

  if (cutoff_mode >= 0) {
    for (k = i = 0; i < 25; i2f[i++] = k + 1) {
      while (model->tbl_q[k] == i) {
        k++;
      }
    }

    const int esc_coef[12] = { 16, -10, 1, 51, 14, 89, 23, 35, 64, 26, -42, 43 };

    for (k = 0; k < 64; k++) {
      for (s = i = 0; i < 6; i++) {
        s += esc_coef[2 * i + ((k >> i) & 1)];
      }

      s = 128 * PPMD::Clamp(s, 32, 256 - 32);

      for (i = 0; i < 25; i++) {
        model->tbl_bin_see[i][k] = BIN_SCALE - s / i2f[i];
      }
    }

    for (i = 0; i < 23; i++) {
      for (k = 0; k < 32; k++) {
        model->tbl_see[i][k].Init(8 * i + 5);
      }
    }
  }
}

void PPMD::AuxCutOff(PPMD::Model *model, PPMD::Alloc *alloc, State *p, int order) {
  if (order < model->order_max) {
    p->idx_successor = CutOff(model, alloc, p->GetSuccessor(alloc), order + 1);
  }
  else {
    p->idx_successor = 0;
  }
}

void PPMD::RestoreModel(PPMD::Model *model, PPMD::Alloc *alloc, Context *pc) {
  printf("RestoreModel %d/%d KB -> ", alloc->GetUsedMemory() >> 10, alloc->size >> 10);

  if (!model->cutoff_mode || alloc->GetUsedMemory() < (alloc->size >> 1)) {
    StartModel(model, alloc, model->order_max, -1);
  }
  else {
    alloc->text_area = alloc->heap_start;

    while (model->max_context->idx_suffix) {
      model->max_context = model->max_context->GetSuffix(alloc);
    }

    model->order_fall = model->order_max;
    CutOff(model, alloc, model->max_context, 0);
  }

  printf("%d KB\n", alloc->GetUsedMemory() >> 10);
}

u32 PPMD::ReduceOrder(PPMD::Model *model, PPMD::Alloc *alloc, State *p, Context *pc) {
  Context *pc1 = pc;
  u32 idx_up = model->found_state->idx_successor = alloc->PtrToIndex(alloc->text_area);
  u8 sym = model->found_state->symbol;

  model->order_fall++;

  if (p) {
    pc = pc->GetSuffix(alloc);
    goto LOOP_ENTRY;
  }

  for (;;) {
    if (!pc->idx_suffix) {
      return alloc->PtrToIndex(pc);
    }

    pc = pc->GetSuffix(alloc);
    if (pc->num_stats) {
      p = pc->GetStats(alloc);
      while (p->symbol != sym) {
        ++p;
      }

      u32 cf = 2 * (p->freq < MAX_FREQ - 3);
      p->freq += cf;
      pc->transition.sum_freq += cf;
    }
    else {
      p = &(pc->one);
      p->freq += p->freq < 11;
    }

  LOOP_ENTRY:
    if (p->idx_successor) {
      break;
    }

    p->idx_successor = idx_up;
    model->order_fall++;
  }

  if (p->idx_successor <= idx_up) {
    p->idx_successor = CreateSuccessors(alloc, p, false, nullptr, pc);
  }

  if (model->order_fall == 1 && pc1 == model->max_context) {
    model->found_state->idx_successor = p->idx_successor;
    alloc->text_area--;
  }

  return p->idx_successor;
}

u32 PPMD::CreateSuccessors(PPMD::Alloc *alloc, State *fs, bool skip, State *p, Context *pc) {
  u32 idx_up_text = fs->idx_successor;
  State *pstack[MAX_ORDER + 1], **pcur = pstack;
  u8 sym = fs->symbol;

  if (!skip) {
    *pcur++ = fs;
    if (!pc->idx_suffix) {
      goto NO_LOOP;
    }
  }

  if (p) {
    pc = pc->GetSuffix(alloc);
    goto LOOP_ENTRY;
  }

  do {
    pc = pc->GetSuffix(alloc);
    if (pc->num_stats) {
      p = pc->GetStats(alloc);
      while (p->symbol != sym) {
        ++p;
      }

      if (p->freq < MAX_FREQ) {
        p->freq += 1;
        pc->transition.sum_freq += 1;
      }
    }
    else {
      p = &(pc->one);
      p->freq += p->freq < 11;
    }

  LOOP_ENTRY:
    if (p->idx_successor != idx_up_text) {
      pc = p->GetSuccessor(alloc);
      break;
    }

    *pcur++ = p;
  }
  while (pc->idx_suffix);

NO_LOOP:
  if (pcur == pstack) {
    return alloc->PtrToIndex(pc);
  }

  Context ct;
  ct.num_stats = 0;
  ct.flags = 0x10 * (sym >= 0x40);
  ct.one.symbol = sym = *alloc->IndexToPtr(idx_up_text);
  ct.one.idx_successor = alloc->PtrToIndex(alloc->IndexToPtr(idx_up_text) + 1);
  ct.flags |= 0x08 * (sym >= 0x40);

  if (pc->num_stats) {
    p = pc->GetStats(alloc);
    while (p->symbol != sym) {
      ++p;
    }

    u32 cf = p->freq - 1;
    u32 s0 = pc->transition.sum_freq - pc->num_stats - cf;
    cf = 1 + ((2 * cf <= s0) ? (12 * cf > s0) : ((cf + 2 * s0) / s0));
    ct.one.freq = cf < 7 ? cf : 7;
  }
  else {
    ct.one.freq = pc->one.freq;
  }

  do {
    Context *pc1 = (Context *)alloc->AllocContext();
    if (!pc1) {
      return 0;
    }

    *pc1 = ct;
    pc1->idx_suffix = alloc->PtrToIndex(pc);
    pc = pc1;

    --pcur;
    (*pcur)->idx_successor = alloc->PtrToIndex(pc);
  }
  while (pcur != pstack);

  return alloc->PtrToIndex(pc);
}

void PPMD::UpdateModel(PPMD::Model *model, PPMD::Alloc *alloc, Context *min_ctx) {
  // Tabulated escapes for exponential symbol distribution
  const u8 ExpEscape[16] = { 51, 43, 18, 12, 11, 9, 8, 7, 6, 5, 4, 3, 3, 2, 2, 2 };

  u8 assign_flag, found_sym = model->found_state->symbol;
  u32 flag, ns, cf, sf, s0, found_freq = model->found_state->freq;
  u32 idx_successor, idx_found_successor = model->found_state->idx_successor;

  Context *pc = nullptr;
  State *p = nullptr;

  if (min_ctx->idx_suffix) {
    pc = min_ctx->GetSuffix(alloc);

    if (pc->num_stats) {
      if ((p = pc->GetStats(alloc))->symbol != found_sym) {
        do {
          p++;
        }
        while (p->symbol != found_sym);

        if (p[0].freq >= p[-1].freq) {
          PPMD::Swap(p[0], p[-1]);
          p--;
        }
      }

      if (p->freq < MAX_FREQ) {
        cf = 1 + (found_freq < 32);
        p->freq += cf;
        pc->transition.sum_freq += cf;
      }
    }
    else {
      p = &(pc->one);
      p->freq += p->freq < 11;
    }
  }

  pc = model->max_context;
  if (!model->order_fall && idx_found_successor) {
    model->found_state->idx_successor = CreateSuccessors(alloc, model->found_state, true, p, min_ctx);
    if (!model->found_state->idx_successor) {
      goto RESTART_MODEL;
    }

    model->max_context = model->found_state->GetSuccessor(alloc);
    return;
  }

  *alloc->text_area++ = found_sym;
  idx_successor = alloc->PtrToIndex(alloc->text_area);
  if (alloc->text_area >= alloc->units_start) {
    goto RESTART_MODEL;
  }

  if (idx_found_successor) {
    // Successor is text pointer
    if (alloc->IndexToPtr(idx_found_successor) < alloc->units_start) {
      idx_found_successor = CreateSuccessors(alloc, model->found_state, false, p, min_ctx);
    }
  }
  else {
    idx_found_successor = ReduceOrder(model, alloc, p, min_ctx);
  }

  if (!idx_found_successor) {
    goto RESTART_MODEL;
  }

  if (!--model->order_fall) {
    idx_successor = idx_found_successor;
    alloc->text_area -= model->max_context != min_ctx;
  }

  s0 = min_ctx->transition.sum_freq - found_freq;
  ns = min_ctx->num_stats;
  assign_flag = 0x08 * (found_sym >= 0x40);

  for (; pc != min_ctx; pc = pc->GetSuffix(alloc)) {
    if ((flag = pc->num_stats) != 0) {
      if ((flag & 1) != 0) {
        p = alloc->ExpandUnits(pc->GetStats(alloc), (flag + 1) >> 1);
        if (!p) {
          goto RESTART_MODEL;
        }

        pc->transition.idx_stats = alloc->PtrToIndex(p);
      }

      pc->transition.sum_freq += model->tbl_q[ns + 4] >> 3;
    }
    else {
      p = alloc->AllocUnits(1);
      if (!p) {
        goto RESTART_MODEL;
      }

      *p = pc->one;
      pc->transition.idx_stats = alloc->PtrToIndex(p);
      p->freq = (p->freq <= MAX_FREQ / 3) ? (2 * p->freq - 1) : (MAX_FREQ - 15);
      pc->transition.sum_freq = p->freq + (ns > 1) + ExpEscape[model->tbl_q[model->bsum >> 8]];
    }

    cf = 2 * found_freq * (pc->transition.sum_freq + 4);
    sf = s0 + pc->transition.sum_freq;
    if (cf <= 6 * sf) {
      cf = 1 + (cf > sf) + (cf > 3 * sf);
      pc->transition.sum_freq += 4;
    }
    else {
      cf = 4 + (cf > 8 * sf) + (cf > 10 * sf) + (cf > 13 * sf);
      pc->transition.sum_freq += cf;
    }

    ++pc->num_stats;
    p = pc->GetStats(alloc) + pc->num_stats;
    p->idx_successor = idx_successor;
    p->symbol = found_sym;
    p->freq = cf;
    pc->flags |= assign_flag;
  }

  model->max_context = alloc->IndexToContext(idx_found_successor);
  return;

RESTART_MODEL:
  RestoreModel(model, alloc, pc);
}

void PPMD::Rescale(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx) {
  State *p = model->found_state;
  State *stats = ctx->GetStats(alloc);

  ctx->flags &= 0x14;

  // Current node to rank 0
  for (p = model->found_state; p != stats; p--) {
    PPMD::Swap(p[0], p[-1]);
  }

  bool cur_order_fall = model->order_fall != 0;
  int f0 = p->freq;
  int sum_freq0 = ctx->transition.sum_freq;
  int escape_freq = sum_freq0 - f0;

  ctx->transition.sum_freq = p->freq = (f0 + cur_order_fall) >> 1;

  // Sort symbols by freqs
  for (int i = 1; i <= ctx->num_stats; i++) {
    p++;

    escape_freq -= p->freq;
    p->freq = (p->freq + cur_order_fall) >> 1;
    ctx->transition.sum_freq += p->freq;

    if (p->freq) {
      ctx->flags |= 0x08 * (p->symbol >= 0x40);
    }

    if (p->freq > p[-1].freq) {
      State *p1 = p;
      State tmp = *p;
      do {
        p1[0] = p1[-1];
      }
      while (tmp.freq > (--p1)[-1].freq);

      *p1 = tmp;
    }
  }

  // Remove zero freq nodes
  if (!p->freq) {
    int num_removes = 0;
    do {
      ++num_removes;
      --p;
    }
    while (!p->freq);

    escape_freq += num_removes;
    int old_units = (ctx->num_stats + 2) >> 1;
    ctx->num_stats -= num_removes;

    if (!ctx->num_stats) {
      ctx->one = *stats;
      alloc->FreeUnits(stats, old_units);

      ctx->flags &= 0x18;
      ctx->one.freq = PPMD::Min(MAX_FREQ / 3, (2 * ctx->one.freq + escape_freq - 1) / escape_freq);

      model->found_state = &ctx->one;
      return;
    }

    int new_units = (ctx->num_stats + 2) >> 1;
    stats = alloc->ShrinkUnits(stats, old_units, new_units);
    ctx->transition.idx_stats = alloc->PtrToIndex(stats);
  }

  ctx->transition.sum_freq += (escape_freq + 1) >> 1;

  int cf;
  if (model->order_fall || (ctx->flags & 0x04) == 0) {
    sum_freq0 -= escape_freq;
    
    cf = sum_freq0 - f0;
    cf = PPMD::Clamp(u32((f0 * ctx->transition.sum_freq - sum_freq0 * stats->freq + cf - 1) / cf), 2U, MAX_FREQ / 2U - 18U);
  }
  else {
    cf = 2;
  }

  model->found_state = stats;
  model->found_state->freq += cf;
  ctx->transition.sum_freq += cf;
  ctx->flags |= 0x04;
}

u32 PPMD::CutOff(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx, int order) {
  int new_num_stats, num_units, esc_freq, scale;
  State *p;
  if (!ctx->num_stats) {
    p = &ctx->one;

    if ((u8 *)p->GetSuccessor(alloc) >= alloc->units_start) {
      AuxCutOff(model, alloc, p, order);

      if (p->idx_successor || order < O_BOUND) {
        return alloc->PtrToIndex(ctx);
      }
    }

    alloc->FreeContext(ctx);
    return 0;
  }

  num_units = (ctx->num_stats + 2) >> 1;
  State *p0 = ctx->GetStats(alloc);
  p0 = alloc->MoveUnitsUp(p0, num_units);
  ctx->transition.idx_stats = alloc->PtrToIndex(p0);

  new_num_stats = ctx->num_stats;
  for (p = p0 + new_num_stats; p >= p0; p--) {
    // Successor is text pointer
    if ((u8 *)p->GetSuccessor(alloc) < alloc->units_start) {
      p->idx_successor = 0;
      PPMD::Swap(*p, p0[new_num_stats--]);
    }
    else {
      AuxCutOff(model, alloc, p, order);
    }
  }

  if (new_num_stats != ctx->num_stats && order) {
    ctx->num_stats = new_num_stats;
    p = p0;

    if (new_num_stats < 0) {
      alloc->FreeUnits(p, num_units);
      alloc->FreeContext(ctx);

      return 0;
    }
    else if (new_num_stats == 0) {
      ctx->flags = (ctx->flags & 0x10) + 0x08 * (p->symbol >= 0x40);
      p->freq = 1 + (2 * (p->freq - 1)) / (ctx->transition.sum_freq - p->freq);
      ctx->one = *p;
      alloc->FreeUnits(p, num_units);
    }
    else {
      int new_units = (new_num_stats + 2) >> 1;
      p = alloc->ShrinkUnits(p0, num_units, new_units);
      ctx->transition.idx_stats = alloc->PtrToIndex(p);

      scale = ctx->transition.sum_freq > 16 * new_num_stats;
      esc_freq = ctx->transition.sum_freq - p->freq;
      ctx->flags = (ctx->flags & (0x10 + 0x04 * scale)) + 0x08 * (p->symbol >= 0x40);
      p->freq = (p->freq + scale) >> scale;
      ctx->transition.sum_freq = p->freq;

      do {
        ++p;

        esc_freq -= p->freq;
        p->freq = (p->freq + scale) >> scale;
        ctx->transition.sum_freq += p->freq;
        ctx->flags |= 0x08 * (p->symbol >= 0x40);
      }
      while (--new_num_stats);

      esc_freq = (esc_freq + scale) >> scale;
      ctx->transition.sum_freq += esc_freq;
    }
  }

  return alloc->PtrToIndex(ctx);
}

PPMD::SEE *PPMD::GetSEE(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx, int num_stats, int flags) {
  if (ctx->num_stats != 0xFF) {
    int i0 = model->tbl_q[ctx->num_stats + 3] - 4;
    int i1 = (ctx->transition.sum_freq > 10 * (ctx->num_stats + 1)) + 
      (2 * (2 * ctx->num_stats < ctx->GetSuffix(alloc)->num_stats + model->num_masked) + ctx->flags);
    return &model->tbl_see[i0][i1];
  }
  else {
    return nullptr;
  }
}

u16 *PPMD::GetBinSEE(PPMD::Model *model, PPMD::Alloc *alloc, PPMD::Context *ctx) {
  PPMD::State &rs = ctx->one;
  int c = model->ns_to_bs_index[ctx->GetSuffix(alloc)->num_stats] + model->prev_success_see + ctx->flags + ((model->run_len >> 26) & 0x20);
  return &model->tbl_bin_see[model->tbl_q[rs.freq - 1]][c];
}
