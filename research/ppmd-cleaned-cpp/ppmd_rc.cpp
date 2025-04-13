#include "ppmd.h"

void PPMD::RangeCoder::InitEncoder() {
  low = 0;
  range = u32(-1);
}

void PPMD::RangeCoder::InitDecoder(const u8 **stream) {
  low = 0;
  range = u32(-1);

  code = 0;

  const u8 *ptr = *stream;
  for (int i = 0; i < 4; i++) {
    code = (code << 8) | *ptr++;
  }

  *stream = ptr;
}

void PPMD::RangeCoder::FlushEncoder(u8 **stream) {
  u8 *ptr = *stream;

  for (int i = 0; i < 4; i++) {
    *ptr++ = u8(low >> 24);
    low <<= 8;
  }

  *stream = ptr;
}

void PPMD::RangeCoder::NormalizeEncoder(u8 **stream) {
  u8 *ptr = *stream;

  while ((low ^ (low + range)) < TOP ||
    (range < BOT && ((range = -(int)low & (BOT - 1)), 1))) {
    *ptr++ = u8(low >> 24);

    range <<= 8;
    low <<= 8;
  }

  *stream = ptr;
}

void PPMD::RangeCoder::NormalizeDecoder(const u8 **stream) {
  const u8 *ptr = *stream;

  while ((low ^ (low + range)) < TOP ||
    (range < BOT && ((range = -(int)low & (BOT - 1)), 1))) {
    code = (code << 8) | *ptr++;

    range <<= 8;
    low <<= 8;
  }

  *stream = ptr;
}

void PPMD::RangeCoder::Arrange(u32 scale) {
  range /= scale;
}

void PPMD::RangeCoder::Interval(u32 rlow, u32 high) {
  low += rlow * range;
  range *= high - rlow;
}

u32 PPMD::RangeCoder::GetCurrentCount() {
  return (code - low) / range;
}

u32 PPMD::RangeCoder::BinStart(u32 f0, u32 shift) {
  range_bin = range;
  range >>= shift;
  return f0 * range;
}

bool PPMD::RangeCoder::BinDecode(u32 tmp) {
  return code - low >= tmp;
}

void PPMD::RangeCoder::BinCorrect0(u32 tmp) {
  range = tmp;
}

void PPMD::RangeCoder::BinCorrect1(u32 tmp) {
  low += tmp;
  range = range_bin - tmp;
}

