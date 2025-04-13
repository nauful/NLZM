#include "ppmd.h"

void PPMD::EncodeFile(PPMD::Model *model, PPMD::Alloc *alloc, u8 **DecodedFile, u8 **EncodedFile, int order_max, u8 cutoff_mode) {
  int symbol = 0;

  PPMD::RangeCoder rc;
  rc.InitEncoder();

  PPMD::StartModel(model, alloc, order_max, cutoff_mode);
  PPMD::Context *min_context = model->max_context;

  while (1) {
    u8 *ptr = *DecodedFile;
    symbol = *ptr++;
    *DecodedFile = ptr;

    if (!symbol) {
      symbol = -1;
    }

    if (!min_context->num_stats) {
      u16 *bs = GetBinSEE(model, alloc, min_context);
      model->bsum = *bs;

      u32 tmp = rc.BinStart(*bs, PPMD::BIN_TOT_BITS);
      *bs -= (*bs + PPMD::BIN_ROUND) >> PPMD::BIN_PERIOD_BITS;

      PPMD::State &rs = min_context->one;
      if (rs.symbol == symbol) {
        rc.BinCorrect0(tmp);

        *bs += PPMD::BIN_INTERVAL;

        rs.freq += rs.freq < 196;
        model->run_len++;
        model->prev_success_see = 1;
        model->found_state = &rs;
      }
      else {
        rc.BinCorrect1(tmp);

        model->EscMaskClear();
        model->EscMaskSet(rs.symbol);

        model->num_masked = 0;
        model->prev_success_see = 0;
        model->found_state = nullptr;
      }
    }
    else {
      PPMD::State *p = min_context->GetStats(alloc);
      int i = p->symbol;
      u32 lo_cnt = p->freq;

      rc.Arrange(min_context->transition.sum_freq);

      if (i == symbol) {
        rc.Interval(0, lo_cnt);
        model->prev_success_see = 2 * lo_cnt > min_context->transition.sum_freq;
        p->freq += 4;
        min_context->transition.sum_freq += 4;
      }
      else {
        model->prev_success_see = 0;

        for (i = 1; i <= min_context->num_stats; i++) {
          if (p[i].symbol == symbol) {
            break;
          }

          lo_cnt += p[i].freq;
        }

        if (i <= min_context->num_stats) {
          rc.Interval(lo_cnt, lo_cnt + p[i].freq);
          p[i].freq += 4;
          min_context->transition.sum_freq += 4;

          if (p[i].freq > p[i - 1].freq) {
            PPMD::Swap(p[i], p[i - 1]);
            i--;
          }

          p = &p[i];
        }
        else {
          rc.Interval(lo_cnt, min_context->transition.sum_freq);
          model->num_masked = min_context->num_stats;

          model->EscMaskClear();
          for (i = 0; i <= min_context->num_stats; i++) {
            model->EscMaskSet(p[i].symbol);
          }

          p = nullptr;
        }
      }

      model->found_state = p;
      if (p && p->freq > PPMD::MAX_FREQ) {
        Rescale(model, alloc, min_context);
      }
    }

    while (!model->found_state) {
      rc.NormalizeEncoder(EncodedFile);
      do {
        if (!min_context->idx_suffix) {
          rc.FlushEncoder(EncodedFile);
          return;
        }

        model->order_fall++;
        min_context = min_context->GetSuffix(alloc);
      }
      while (min_context->num_stats == model->num_masked);

      PPMD::State *p = min_context->GetStats(alloc);
      PPMD::SEE *see = GetSEE(model, alloc, min_context, min_context->num_stats, min_context->flags);
      int see_freq = see ? see->Mean() : 1;
    
      int i, found_idx = 0, flag_found_sym = 0;
      u32 lo_cnt = 0, sum_cnt = 0;
      for (i = 0; i <= min_context->num_stats; i++) {
        u8 sym = p[i].symbol;
        if (model->EscMaskIsSet(sym)) {
          continue;
        }
    
        if (sym == symbol) {
          flag_found_sym = 1;
          found_idx = i;
          lo_cnt = sum_cnt;
        }

        model->EscMaskSet(sym);
        sum_cnt += p[i].freq;
      }
    
      int total = PPMD::Max(see_freq, 1) + sum_cnt;
      rc.Arrange(total);
    
      if (flag_found_sym) {
        p += found_idx;
        rc.Interval(lo_cnt, lo_cnt + p->freq);
    
        if (see) {
          see->acc -= see_freq;
          see->UpdateSymbol();
        }
    
        model->found_state = p;
        p->freq += 4;
        min_context->transition.sum_freq += 4;
        if (p->freq > PPMD::MAX_FREQ) {
          Rescale(model, alloc, min_context);
        }
    
        model->run_len = model->run_len_init;
      }
      else {
        rc.Interval(sum_cnt, total);
        model->num_masked = min_context->num_stats;
    
        if (see) {
          see->acc += sum_cnt;
        }
      }
    }

    if (!model->order_fall && (u8 *)model->found_state->GetSuccessor(alloc) >= alloc->units_start) {
      model->max_context = model->found_state->GetSuccessor(alloc);
    }
    else {
      UpdateModel(model, alloc, min_context);
    }

    rc.NormalizeEncoder(EncodedFile);
    min_context = model->max_context;
  }

  rc.FlushEncoder(EncodedFile);
}

void PPMD::DecodeFile(PPMD::Model *model, PPMD::Alloc *alloc, u8 **DecodedFile, const u8 **EncodedFile, int order_max, u8 cutoff_mode) {
  PPMD::RangeCoder rc;
  rc.InitDecoder(EncodedFile);

  PPMD::StartModel(model, alloc, order_max, cutoff_mode);
  PPMD::Context *min_context = model->max_context;

  while (1) {
    if (!min_context->num_stats) {
      u16 *bs = GetBinSEE(model, alloc, min_context);
      model->bsum = *bs;

      u32 tmp = rc.BinStart(*bs, PPMD::BIN_TOT_BITS);
      *bs -= (*bs + PPMD::BIN_ROUND) >> PPMD::BIN_PERIOD_BITS;

      PPMD::State &rs = min_context->one;
      if (!rc.BinDecode(tmp)) {
        *bs += PPMD::BIN_INTERVAL;
        rc.BinCorrect0(tmp);

        model->found_state = &rs;
        rs.freq += rs.freq < 196;

        model->run_len++;
        model->prev_success_see = 1;
      }
      else {
        rc.BinCorrect1(tmp);

        model->EscMaskClear();
        model->EscMaskSet(rs.symbol);

        model->num_masked = 0;
        model->prev_success_see = 0;
        model->found_state = nullptr;
      }
    }
    else {
      PPMD::State *p = min_context->GetStats(alloc);
      u32 i = p->symbol;
      u32 lo_cnt = p->freq;

      rc.Arrange(min_context->transition.sum_freq);

      u32 rc_count = rc.GetCurrentCount();
      if (rc_count < lo_cnt) {
        rc.Interval(0, lo_cnt);
        model->prev_success_see = 2 * lo_cnt > min_context->transition.sum_freq;
        p->freq += 4;
        min_context->transition.sum_freq += 4;
      }
      else {
        model->prev_success_see = 0;

        for (i = 1; i <= min_context->num_stats; i++) {
          if (lo_cnt + p[i].freq > rc_count) {
            break;
          }

          lo_cnt += p[i].freq;
        }

        if (i <= min_context->num_stats) {
          rc.Interval(lo_cnt, lo_cnt + p[i].freq);
          p[i].freq += 4;
          min_context->transition.sum_freq += 4;

          if (p[i].freq > p[i - 1].freq) {
            PPMD::Swap(p[i], p[i - 1]);
            i--;
          }

          p = &p[i];
        }
        else {
          rc.Interval(lo_cnt, min_context->transition.sum_freq);
          model->num_masked = min_context->num_stats;

          model->EscMaskClear();
          for (i = 0; i <= min_context->num_stats; i++) {
            model->EscMaskSet(p[i].symbol);
          }

          p = nullptr;
        }
      }

      model->found_state = p;
      if (p && p->freq > PPMD::MAX_FREQ) {
        Rescale(model, alloc, min_context);
      }
    }

    while (!model->found_state) {
      rc.NormalizeDecoder(EncodedFile);
      do {
        if (!min_context->idx_suffix) {
          return;
        }

        model->order_fall++;
        min_context = min_context->GetSuffix(alloc);
      }
      while (min_context->num_stats == model->num_masked);

      PPMD::State *p = min_context->GetStats(alloc);
      PPMD::SEE *see = GetSEE(model, alloc, min_context, min_context->num_stats, min_context->flags);
      int see_freq = see ? see->Mean() : 1;
    
      u8 px[256];
      int i, found_idx = 0;
      u32 sum_cnt = 0;
      for (i = 0; i <= min_context->num_stats; i++) {
        u8 sym = p[i].symbol;
        if (model->EscMaskIsSet(sym)) {
          continue;
        }
    
        model->EscMaskSet(sym);
        sum_cnt += p[i].freq;
        px[found_idx++] = i;
      }
    
      int total = PPMD::Max(see_freq, 1) + sum_cnt;
      rc.Arrange(total);
    
      u32 count = rc.GetCurrentCount();
      if (count < sum_cnt) {
        u32 hi_cnt = 0;
        i = 0;
    
        do {
          found_idx = px[i];
          hi_cnt += p[found_idx].freq;
          ++i;
        }
        while (hi_cnt <= count);
    
        p += found_idx;
        rc.Interval(hi_cnt - p->freq, hi_cnt);
    
        if (see) {
          see->acc -= see_freq;
          see->UpdateSymbol();
        }
    
        model->found_state = p;
        p->freq += 4;
        min_context->transition.sum_freq += 4;
        if (p->freq > PPMD::MAX_FREQ) {
          Rescale(model, alloc, min_context);
        }
    
        model->run_len = model->run_len_init;
      }
      else {
        rc.Interval(sum_cnt, total);
        model->num_masked = min_context->num_stats;
    
        if (see) {
          see->acc += sum_cnt;
        }
      }
    }

    u8 *ptr = *DecodedFile;
    *ptr++ = model->found_state->symbol;
    *DecodedFile = ptr;

    if (!model->order_fall && (u8 *)model->found_state->GetSuccessor(alloc) >= alloc->units_start) {
      model->max_context = model->found_state->GetSuccessor(alloc);
    }
    else {
      UpdateModel(model, alloc, min_context);
    }

    rc.NormalizeDecoder(EncodedFile);
    min_context = model->max_context;
  }
}
