#include "platform.h"

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "ppmd.h"

struct BenchmarkBlock {
  const char *label;
  LARGE_INTEGER timer_freq, timer_start, timer_cur;

  BenchmarkBlock(const char *label) : label(label) {
    QueryPerformanceFrequency(&timer_freq);
    QueryPerformanceCounter(&timer_start);
    timer_cur = timer_start;
  }

  ~BenchmarkBlock() {
    QueryPerformanceCounter(&timer_cur);
    double dt = (double)(timer_cur.QuadPart - timer_start.QuadPart) / timer_freq.QuadPart;
    printf("%s: %.2fms\n", label, dt * 1000.0);
  }
};

u8 *read_file(const char *fname, u32 *fsize) {
  FILE *Freq = fopen(fname, "rb");
  PLATFORM_ASSERT(Freq);

  fseek(Freq, 0, SEEK_END);
  u32 _fsize = ftell(Freq);
  fseek(Freq, 0, SEEK_SET);

  u8 *res = new u8[_fsize + 1];
  _fsize = fread(res, 1, _fsize, Freq);
  res[_fsize] = 0;
  if (fsize) {
    *fsize = _fsize;
  }

  fclose(Freq);

  return res;
}

int main(int argc, const char **argv) {
  const int suballoc_size = 4;
  const int order_max_d = PPMD::Clamp<int>(6, 1, PPMD::MAX_ORDER);
  const int cutoff_mode = 1;

  u32 buf_len = 0;
  u8 *buf_src = read_file("book.txt", &buf_len);
  u8 *buf_packed = new u8[2 * buf_len];
  u8 *buf_unpacked = new u8[2 * buf_len];

  memset(buf_packed, 0, 2 * buf_len);
  memset(buf_unpacked, 0, 2 * buf_len);

  {
    BenchmarkBlock bench("pack");

    PPMD::Model model;
    PPMD::Alloc alloc;

    PPMD::Init(&model, &alloc);
    PLATFORM_ASSERT(alloc.Start(suballoc_size));
    u8 *tmp_in = buf_src, *tmp_out = buf_packed;
    PPMD::EncodeFile(&model, &alloc, &tmp_in, &tmp_out, order_max_d, cutoff_mode);
    alloc.Stop();

    printf("%d -> %d, %d KB\n", tmp_in - buf_src, tmp_out - buf_packed, (tmp_out - buf_packed) >> 10);
  }

  {
    BenchmarkBlock bench("unpack");

    PPMD::Model model;
    PPMD::Alloc alloc;

    PPMD::Init(&model, &alloc);
    PLATFORM_ASSERT(alloc.Start(suballoc_size));
    u8 *tmp_in = buf_unpacked, *tmp_out = buf_packed;
    PPMD::DecodeFile(&model, &alloc, &tmp_in, (const u8 **)&tmp_out, order_max_d, cutoff_mode);
    alloc.Stop();

    printf("%d -> %d\n", tmp_out - buf_packed, tmp_in - buf_unpacked);
  }

  PLATFORM_ASSERT(!memcmp(buf_src, buf_unpacked, buf_len));

  delete[] buf_src;
  delete[] buf_packed;
  delete[] buf_unpacked;

  _CrtDumpMemoryLeaks();
  return 0;
}
