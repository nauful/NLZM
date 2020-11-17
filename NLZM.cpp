// NLZM 1.03 64-bit - Written by Nauful
// Released into the public domain.
// Please drop a comment if you find this useful.

#ifndef _CRT_DISABLE_PERFCRIT_LOCKS
# define _CRT_DISABLE_PERFCRIT_LOCKS
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
# define _CRT_SECURE_NO_WARNINGS
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef MSVC
# include <crtdbg.h>
# include <intrin.h>
# define _fpos64(f) _ftelli64(f)
# define ASSERT(x) { if (!(x)) { printf("Assert failed " #x "\n"); __debugbreak(); } }
#else
# define _fpos64(f) ftello64(f)
# define ASSERT(x) { if (!(x)) { printf("Assert failed " #x "\n"); exit(-1); } }
#endif

#ifdef DEBUG
# define ASSERT_DEBUG(x) ASSERT(x)
#else
# define ASSERT_DEBUG(x)
#endif

typedef unsigned char byte;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

typedef short int16;
typedef int int32;
typedef long long int64;

template<typename T, size_t N> constexpr size_t CountOf(T const (&)[N]) { return N; }
template<typename T> constexpr T Min(const T a, const T b) { return (a < b) ? a : b; }
template<typename T> constexpr T Max(const T a, const T b) { return (a < b) ? b : a; }

template<typename T> constexpr T Clamp(const T v, const T low, const T high) {
	if (v < low) {
		return low;
	}
	else if (v > high) {
		return high;
	}
	else {
		return v;
	}
}

uint32 clz32(uint32 x) {
#ifdef MSVC
	unsigned long v0;
	_BitScanReverse(&v0, x);
	return v0;
#else
	uint32 r, q;

	r = (x > 0xFFFF) << 4; x >>= r;
	q = (x > 0xFF) << 3; x >>= q; r |= q;
	q = (x > 0xF) << 2; x >>= q; r |= q;
	q = (x > 0x3) << 1; x >>= q; r |= q;
	r |= (x >> 1);

	return r;
#endif
}

uint32 popcnt32(uint32 x) {
	x -= ((x >> 1) & 0x55555555);
	x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
	x = (((x >> 4) + x) & 0x0F0F0F0F);
	x += (x >> 8);
	x += (x >> 16);

	return x & 0x0000003F;
}

uint32 ctz32(uint32 x) {
#ifdef MSVC
	unsigned long v0;
	_BitScanForward(&v0, x);
	return v0;
#else
	return popcnt32((x & uint32(-int32(x))) - 1);
#endif
}

const int LOG2_LUT_SIZE_BITS = 8;
const int LOG2_LUT_SCALE_BITS = 5;
const int LOG2_LUT_PRECISION = 16;

uint16 log2_lut[1 << LOG2_LUT_SIZE_BITS];

void log2_init() {
	const int table_size = 1 << LOG2_LUT_SIZE_BITS;
	const int half_table_size = 1 << (LOG2_LUT_SIZE_BITS - 1);
	const int scale = 1 << LOG2_LUT_SCALE_BITS;

	for (int i = 1; i < table_size; i++) {
		uint32 next = 1 << LOG2_LUT_PRECISION;
		uint16 acc = 0;

		for (int sum = 0; sum < scale; sum++) {
			const uint32 v = (i * next) >> LOG2_LUT_SIZE_BITS;
			const uint32 num_bits = LOG2_LUT_PRECISION - clz32(v);

			acc += num_bits - 1;
			next = v << (num_bits - 1);
		}

		log2_lut[i] = acc;
	}

	log2_lut[0] = log2_lut[1];
}

static uint32 crc32_table[16][256];

void crc32_init() {
	// Reversed bit order for CRC32
	const static uint32 POLY = 0xEDB88320;

	uint32 n, crc, k;

	for (n = 0; n < 256; n++) {
		crc = n;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
		crc32_table[0][n] = crc;
	}

	for (n = 0; n < 256; n++) {
		crc = crc32_table[0][n];

		for (k = 1; k < 16; k++) {
			crc = crc32_table[0][crc & 0xff] ^ (crc >> 8);
			crc32_table[k][n] = crc;
		}
	}
}

uint32 crc32_calc(const byte *window, int64 n, uint64 crci) {
	uint64 crc = crci ^ 0xFFFFFFFF;
	const byte *next = window;

	while (n && (byte(size_t(next)) & 15)) {
		crc = crc32_table[0][(crc ^ (*next++)) & 0xFF] ^ (crc >> 8);
		--n;
	}

	while (n >= 16) {
		crc ^= *reinterpret_cast<const uint64 *>(next);
		uint64 high = *reinterpret_cast<const uint64 *>(next + 8);

		crc =
			crc32_table[15][crc & 0xFF] ^
			crc32_table[14][(crc >> 8) & 0xFF] ^
			crc32_table[13][(crc >> 16) & 0xFF] ^
			crc32_table[12][(crc >> 24) & 0xFF] ^
			crc32_table[11][(crc >> 32) & 0xFF] ^
			crc32_table[10][(crc >> 40) & 0xFF] ^
			crc32_table[9][(crc >> 48) & 0xFF] ^
			crc32_table[8][crc >> 56] ^
			crc32_table[7][high & 0xFF] ^
			crc32_table[6][(high >> 8) & 0xFF] ^
			crc32_table[5][(high >> 16) & 0xFF] ^
			crc32_table[4][(high >> 24) & 0xFF] ^
			crc32_table[3][(high >> 32) & 0xFF] ^
			crc32_table[2][(high >> 40) & 0xFF] ^
			crc32_table[1][(high >> 48) & 0xFF] ^
			crc32_table[0][high >> 56];

		next += 16;
		n -= 16;
	}

	while (n) {
		crc = crc32_table[0][(crc ^ (*next++)) & 0xFF] ^ (crc >> 8);
		--n;
	}

	crci = crc ^ 0xFFFFFFFF;
	return crci;
}

uint32 crc32_file(FILE *f) {
	byte buf[0x4000];

	uint32 n = 0, crc = 0;
	while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
		crc = crc32_calc(buf, n, crc);
	}

	return crc;
}

const int16 CDF_ADAPT_BITS = 7;
const int16 CDF_SCALE_BITS = 14;

const int16 CDF_ADAPT_TOTAL = 1 << CDF_ADAPT_BITS;
const int16 CDF_SCALE_TOTAL = 1 << CDF_SCALE_BITS;
const int16 CDF_SCALE_TOTAL_MASK = CDF_SCALE_TOTAL - 1;

struct CDF1 {
	uint16 cell[3];
};

struct CDF2 {
	uint16 cell[5];
};

struct CDF3 {
#ifdef USE_SSE
	union {
		__m128i _cell[2];
#endif
		uint16 cell[9];
#ifdef USE_SSE
	};
#endif
};

struct CDF4 {
#ifdef USE_SSE
	union {
		__m128i _cell[3];
#endif
		uint16 cell[17];
#ifdef USE_SSE
	};
#endif
};

const int16 cdf_mixin1[2] = { CDF_SCALE_TOTAL - CDF_ADAPT_TOTAL - 1, CDF_ADAPT_TOTAL + 1 };

int16 cdf_mixin2[1 << 2][1 << 2];
int16 cdf_mixin3[1 << 3][1 << 3];
int16 cdf_mixin4[1 << 4][1 << 4];

#ifdef USE_SSE
__m128i _cdf_mixin3[8];
__m128i _cdf_mixin4[16][2];
#endif

const int16 cdf_initial1[3] = {
	0, CDF_SCALE_TOTAL / 2,
	CDF_SCALE_TOTAL
};

const int16 cdf_initial2[5] = {
	0, CDF_SCALE_TOTAL / 4, 2 * CDF_SCALE_TOTAL / 4, 3 * CDF_SCALE_TOTAL / 4,
	CDF_SCALE_TOTAL
};

const int16 cdf_initial3[9] = {
	0, CDF_SCALE_TOTAL / 8, 2 * CDF_SCALE_TOTAL / 8, 3 * CDF_SCALE_TOTAL / 8,
	4 * CDF_SCALE_TOTAL / 8, 5 * CDF_SCALE_TOTAL / 8, 6 * CDF_SCALE_TOTAL / 8, 7 * CDF_SCALE_TOTAL / 8,
	CDF_SCALE_TOTAL
};

const int16 cdf_initial4[17] = {
	0, CDF_SCALE_TOTAL / 16, 2 * CDF_SCALE_TOTAL / 16, 3 * CDF_SCALE_TOTAL / 16,
	4 * CDF_SCALE_TOTAL / 16, 5 * CDF_SCALE_TOTAL / 16, 6 * CDF_SCALE_TOTAL / 16, 7 * CDF_SCALE_TOTAL / 16,
	8 * CDF_SCALE_TOTAL / 16, 9 * CDF_SCALE_TOTAL / 16, 10 * CDF_SCALE_TOTAL / 16, 11 * CDF_SCALE_TOTAL / 16,
	12 * CDF_SCALE_TOTAL / 16, 13 * CDF_SCALE_TOTAL / 16, 14 * CDF_SCALE_TOTAL / 16, 15 * CDF_SCALE_TOTAL / 16,
	CDF_SCALE_TOTAL
};

template<int num_syms, int adapt_bits, int total>
void init_mixin_table(int16(&mixin)[num_syms][num_syms]) {
	const int mixin_bias = (1 << adapt_bits) - 1 - num_syms;
	ASSERT_DEBUG(mixin_bias > 0);

	for (int y = 0; y < num_syms; y++) {
		for (int x = 0; x <= y; x++) {
			mixin[y][x] = x;
		}

		for (int x = y + 1; x < num_syms; x++) {
			mixin[y][x] = total + x + mixin_bias;
		}
	}
}

void cdf_init_mixins() {
	init_mixin_table<1 << 2, CDF_ADAPT_BITS, CDF_SCALE_TOTAL>(cdf_mixin2);
	init_mixin_table<1 << 3, CDF_ADAPT_BITS, CDF_SCALE_TOTAL>(cdf_mixin3);
	init_mixin_table<1 << 4, CDF_ADAPT_BITS, CDF_SCALE_TOTAL>(cdf_mixin4);

#ifdef USE_SSE
	for (int i = 0; i < 8; i++) {
		_cdf_mixin3[i] = _mm_set_epi16(
			cdf_mixin3[i][7], cdf_mixin3[i][6], cdf_mixin3[i][5], cdf_mixin3[i][4],
			cdf_mixin3[i][3], cdf_mixin3[i][2], cdf_mixin3[i][1], cdf_mixin3[i][0]);
	}

	for (int i = 0; i < 16; i++) {
		_cdf_mixin4[i][0] = _mm_set_epi16(
			cdf_mixin4[i][7], cdf_mixin4[i][6], cdf_mixin4[i][5], cdf_mixin4[i][4],
			cdf_mixin4[i][3], cdf_mixin4[i][2], cdf_mixin4[i][1], cdf_mixin4[i][0]);

		_cdf_mixin4[i][1] = _mm_set_epi16(
			cdf_mixin4[i][15], cdf_mixin4[i][14], cdf_mixin4[i][13], cdf_mixin4[i][12],
			cdf_mixin4[i][11], cdf_mixin4[i][10], cdf_mixin4[i][9], cdf_mixin4[i][8]);
	}
#endif
}

void cdf_init(CDF1 &cdf) {
	for (int i = 0; i < 3; i++) {
		cdf.cell[i] = cdf_initial1[i];
	}
}

void cdf_init(CDF2 &cdf) {
	for (int i = 0; i < 5; i++) {
		cdf.cell[i] = cdf_initial2[i];
	}
}

void cdf_init(CDF3 &cdf) {
	for (int i = 0; i < 9; i++) {
		cdf.cell[i] = cdf_initial3[i];
	}
}

void cdf_init(CDF4 &cdf) {
	for (int i = 0; i < 17; i++) {
		cdf.cell[i] = cdf_initial4[i];
	}
}

void cdf_update(CDF1 &cdf, int y) {
	cdf.cell[1] += (cdf_mixin1[y] - cdf.cell[1]) >> CDF_ADAPT_BITS;
}

void cdf_update(CDF2 &cdf, int y) {
	for (int i = 0; i < 4; i++) {
		cdf.cell[i] += (cdf_mixin2[y][i] - cdf.cell[i]) >> CDF_ADAPT_BITS;
	}
}

void cdf_update(CDF3 &cdf, int y) {
#ifdef USE_SSE
	__m128i upd = _mm_srai_epi16(_mm_sub_epi16(_cdf_mixin3[y], cdf._cell[0]), CDF_ADAPT_BITS);

	cdf._cell[0] = _mm_add_epi16(cdf._cell[0], upd);
#else
	for (int i = 0; i < 8; i++) {
		cdf.cell[i] += (cdf_mixin3[y][i] - cdf.cell[i]) >> CDF_ADAPT_BITS;
	}
#endif
}

void cdf_update(CDF4 &cdf, int y) {
#ifdef USE_SSE
	__m128i upd0 = _mm_srai_epi16(_mm_sub_epi16(_cdf_mixin4[y][0], cdf._cell[0]), CDF_ADAPT_BITS);
	__m128i upd1 = _mm_srai_epi16(_mm_sub_epi16(_cdf_mixin4[y][1], cdf._cell[1]), CDF_ADAPT_BITS);

	cdf._cell[0] = _mm_add_epi16(cdf._cell[0], upd0);
	cdf._cell[1] = _mm_add_epi16(cdf._cell[1], upd1);
#else
	for (int i = 0; i < 16; i++) {
		cdf.cell[i] += (cdf_mixin4[y][i] - cdf.cell[i]) >> CDF_ADAPT_BITS;
	}
#endif
}

int cdf_lookup(const CDF1 &cdf, int f) {
	return f >= cdf.cell[1];
}

int cdf_lookup(const CDF2 &cdf, int f) {
	int r = 2 * (f >= cdf.cell[2]);
	r += f >= cdf.cell[1 + r];

	return r;
}

int cdf_lookup(const CDF3 &cdf, int f) {
#ifdef USE_SSE
	__m128i _f = _mm_cvtsi32_si128(f);
	_f = _mm_shuffle_epi32(_mm_unpacklo_epi16(_f, _f), 0);

	__m128i _cmp0 = _mm_cmpgt_epi16(cdf._cell[0], _f);
	int idx_mask = _mm_movemask_epi8(_cmp0) | 0x10000;

	int idx = ctz32(idx_mask) >> 1;
	return idx - 1;
#else
	int r = 4 * (f >= cdf.cell[4]);
	r += 2 * (f >= cdf.cell[2 + r]);
	r += f >= cdf.cell[1 + r];

	return r;
#endif
}

int cdf_lookup(const CDF4 &cdf, int f) {
#ifdef USE_SSE
	__m128i _f = _mm_cvtsi32_si128(f);
	_f = _mm_shuffle_epi32(_mm_unpacklo_epi16(_f, _f), 0);

	__m128i _cmp0 = _mm_cmpgt_epi16(cdf._cell[0], _f);
	__m128i _cmp1 = _mm_cmpgt_epi16(cdf._cell[1], _f);
	int idx_mask = _mm_movemask_epi8(_mm_packs_epi16(_cmp0, _cmp1)) | 0x10000;

	int idx = ctz32(idx_mask);
	return idx - 1;
#else
	int r = 8 * (f >= cdf.cell[8]);
	r += 4 * (f >= cdf.cell[4 + r]);
	r += 2 * (f >= cdf.cell[2 + r]);
	r += f >= cdf.cell[1 + r];

	return r;
#endif
}

uint16 cdf_cost(const CDF1 &cdf, int y) { return log2_lut[(cdf.cell[y + 1] - cdf.cell[y]) >> (CDF_SCALE_BITS - LOG2_LUT_SIZE_BITS)]; }
uint16 cdf_cost(const CDF2 &cdf, int y) { return log2_lut[(cdf.cell[y + 1] - cdf.cell[y]) >> (CDF_SCALE_BITS - LOG2_LUT_SIZE_BITS)]; }
uint16 cdf_cost(const CDF3 &cdf, int y) { return log2_lut[(cdf.cell[y + 1] - cdf.cell[y]) >> (CDF_SCALE_BITS - LOG2_LUT_SIZE_BITS)]; }
uint16 cdf_cost(const CDF4 &cdf, int y) { return log2_lut[(cdf.cell[y + 1] - cdf.cell[y]) >> (CDF_SCALE_BITS - LOG2_LUT_SIZE_BITS)]; }

typedef uint32 rans_t;

const rans_t RANS_MID = 1 << 16;

rans_t rans_enc_put(rans_t x, byte **pptr, uint32 start, uint32 freq) {
	uint32 x_max = ((RANS_MID >> CDF_SCALE_BITS) << 16);
	x_max *= freq;

	if (x >= x_max) {
		*--(pptr[0]) = (byte)x;
		*--(pptr[0]) = (byte)(x >> 8);
		x >>= 16;
	}

	return ((x / freq) << CDF_SCALE_BITS) + (x % freq) + start;
}

rans_t rans_dec_consume(rans_t x, uint32 start, uint32 freq) {
	return freq * (x >> CDF_SCALE_BITS) + (x & CDF_SCALE_TOTAL_MASK) - start;
}

void rans_enc_flush(rans_t r, byte **pptr) {
	uint32 x = r;

	pptr[0] -= 4;
	pptr[0][0] = (byte)(x);
	pptr[0][1] = (byte)(x >> 8);
	pptr[0][2] = (byte)(x >> 16);
	pptr[0][3] = (byte)(x >> 24);
}

rans_t rans_dec_init(byte **pptr) {
	uint32 x = pptr[0][0];
	x |= pptr[0][1] << 8;
	x |= pptr[0][2] << 16;
	x |= pptr[0][3] << 24;
	pptr[0] += 4;

	return x;
}

rans_t rans_dec_normalize(rans_t x, byte **pptr) {
	if (x < RANS_MID) {
		x = (x << 16) + (pptr[0][0] << 8) + pptr[0][1];
		*pptr += 2;
	}

	return x;
}

struct CodeFrame {
	byte *start, *end, *ptr_bits;

	uint32 word, word_bits;
	uint32 write_limit;
	uint32 num_ops;

	uint32 *buf_rans, num_rans_syms, max_rans_syms, est_rans_bits;

	void Init(byte *start, byte *end, uint32 *buf_rans, uint32 max_rans_syms);

	bool Size() const;
	bool NeedsFlush() const;

	void WriteRange(uint16 start, uint16 freq);
	void WriteBits(uint32 v, uint32 nb);

	void WriteCDF(const CDF1 &cdf, int y);
	void WriteCDF(const CDF2 &cdf, int y);
	void WriteCDF(const CDF3 &cdf, int y);
	void WriteCDF(const CDF4 &cdf, int y);

	uint32 Flush();
};

struct DecodeFrame {
	byte *start, *end, *ptr_bits, *ptr_rans;

	uint32 word, word_bits;
	uint32 num_ops;

	rans_t rans_state[4];
	byte rans_index;

	uint32 Init(byte *start, byte *end);

	int ReadCDF(const CDF1 &cdf);
	int ReadCDF(const CDF2 &cdf);
	int ReadCDF(const CDF3 &cdf);
	int ReadCDF(const CDF4 &cdf);

	uint32 ReadBits(uint32 nb);
};

void CodeFrame::Init(byte *start, byte *end, uint32 *buf_rans, uint32 max_rans_syms) {
	ASSERT_DEBUG(end - start >= 1024);

	this->start = start;
	this->end = end;
	this->ptr_bits = start + 12;

	this->word = 0;
	this->word_bits = 0;
	this->write_limit = (15 * (end - start)) / 16;
	this->num_ops = 0;

	this->buf_rans = buf_rans;
	this->max_rans_syms = max_rans_syms;
	this->num_rans_syms = 0;
	this->est_rans_bits = 0;
}

bool CodeFrame::Size() const { return end - start; }

bool CodeFrame::NeedsFlush() const {
	return num_rans_syms + 8 >= max_rans_syms ||
		(ptr_bits - start) + (est_rans_bits >> (8 + LOG2_LUT_SCALE_BITS)) + 64 >= write_limit;
}

void CodeFrame::WriteRange(uint16 start, uint16 freq) {
	ASSERT(num_rans_syms <= max_rans_syms);
	ASSERT(freq > 0);

	++num_ops;

	buf_rans[num_rans_syms++] = (freq << 16) + start;
	est_rans_bits += log2_lut[freq >> (CDF_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
}

void CodeFrame::WriteCDF(const CDF1 &cdf, int y) { WriteRange(cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]); }
void CodeFrame::WriteCDF(const CDF2 &cdf, int y) { WriteRange(cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]); }
void CodeFrame::WriteCDF(const CDF3 &cdf, int y) { WriteRange(cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]); }
void CodeFrame::WriteCDF(const CDF4 &cdf, int y) { WriteRange(cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]); }

void CodeFrame::WriteBits(uint32 v, uint32 nb) {
	ASSERT(v < (1u << nb));
	++num_ops;

	word |= v << (32 - word_bits - nb);
	word_bits += nb;

	while (word_bits >= 8) {
		ASSERT(ptr_bits < end);

		*ptr_bits++ = word >> 24;
		word <<= 8;
		word_bits -= 8;
	}
}

uint32 CodeFrame::Flush() {
	for (int i = 0; i < 4; i++) {
		ASSERT(ptr_bits < end);

		*ptr_bits++ = word >> 24;
		word <<= 8;
		word_bits -= Min(word_bits, 8u);
	}

	byte *wptr = end - 1;
	rans_t st[4] = { RANS_MID, RANS_MID, RANS_MID, RANS_MID };
	for (uint32 i = num_rans_syms - 1; i + 1 > 0; i--) {
		st[i & 3] = rans_enc_put(st[i & 3], &wptr, buf_rans[i] & 0xFFFF, buf_rans[i] >> 16);
	}

	rans_enc_flush(st[3], &wptr);
	rans_enc_flush(st[2], &wptr);
	rans_enc_flush(st[1], &wptr);
	rans_enc_flush(st[0], &wptr);

	ASSERT(wptr >= ptr_bits);
	uint32 num_rans_bytes = (end - 1) - wptr;
	memmove(ptr_bits, wptr, num_rans_bytes);

	uint32 num_bits_bytes = ptr_bits - start;
	start[0] = byte(num_ops >> 24);
	start[1] = byte(num_ops >> 16);
	start[2] = byte(num_ops >> 8);
	start[3] = byte(num_ops);

	start[4] = byte(num_bits_bytes >> 24);
	start[5] = byte(num_bits_bytes >> 16);
	start[6] = byte(num_bits_bytes >> 8);
	start[7] = byte(num_bits_bytes);

	start[8] = byte(num_rans_bytes >> 24);
	start[9] = byte(num_rans_bytes >> 16);
	start[10] = byte(num_rans_bytes >> 8);
	start[11] = byte(num_rans_bytes);

	ptr_bits = start + 12;

	word = 0;
	word_bits = 0;
	num_ops = 0;

	num_rans_syms = 0;
	est_rans_bits = 0;

	return num_bits_bytes + num_rans_bytes;
}

uint32 DecodeFrame::Init(byte *start, byte *end) {
	this->start = start;
	this->end = end;
	this->num_ops = (start[0] << 24) + (start[1] << 16) + (start[2] << 8) + start[3];
	if (!this->num_ops) {
		return -1;
	}

	this->ptr_bits = start + 12;
	uint32 num_bits_bytes = (start[4] << 24) + (start[5] << 16) + (start[6] << 8) + start[7];
	uint32 num_rans_bytes = (start[8] << 24) + (start[9] << 16) + (start[10] << 8) + start[11];

	this->ptr_rans = start + num_bits_bytes;
	this->word = 0;
	this->word_bits = 0;

	for (int i = 0; i < 4; i++) {
		this->rans_state[i] = rans_dec_init(&ptr_rans);
	}
	this->rans_index = 0;

	return num_bits_bytes + num_rans_bytes;
}

int DecodeFrame::ReadCDF(const CDF1 &cdf) {
	--num_ops;
	rans_t &rs = rans_state[rans_index++ & 3];

	uint16 f = rs & CDF_SCALE_TOTAL_MASK;
	int y = cdf_lookup(cdf, f);
	rs = rans_dec_consume(rs, cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]);
	rs = rans_dec_normalize(rs, &ptr_rans);

	return y;
}

int DecodeFrame::ReadCDF(const CDF2 &cdf) {
	--num_ops;
	rans_t &rs = rans_state[rans_index++ & 3];

	uint16 f = rs & CDF_SCALE_TOTAL_MASK;
	int y = cdf_lookup(cdf, f);
	rs = rans_dec_consume(rs, cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]);
	rs = rans_dec_normalize(rs, &ptr_rans);

	return y;
}

int DecodeFrame::ReadCDF(const CDF3 &cdf) {
	--num_ops;
	rans_t &rs = rans_state[rans_index++ & 3];

	uint16 f = rs & CDF_SCALE_TOTAL_MASK;
	int y = cdf_lookup(cdf, f);
	rs = rans_dec_consume(rs, cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]);
	rs = rans_dec_normalize(rs, &ptr_rans);

	return y;
}

int DecodeFrame::ReadCDF(const CDF4 &cdf) {
	--num_ops;
	rans_t &rs = rans_state[rans_index++ & 3];

	uint16 f = rs & CDF_SCALE_TOTAL_MASK;
	int y = cdf_lookup(cdf, f);
	rs = rans_dec_consume(rs, cdf.cell[y], cdf.cell[y + 1] - cdf.cell[y]);
	rs = rans_dec_normalize(rs, &ptr_rans);

	return y;
}

uint32 DecodeFrame::ReadBits(uint32 nb) {
	--num_ops;

	while (word_bits < 24) {
		ASSERT(ptr_bits < end);

		word |= *ptr_bits++ << (24 - word_bits);
		word_bits += 8;
	}

	ASSERT(nb <= word_bits);

	uint32 y = word >> (32 - nb);
	word <<= nb;
	word_bits -= nb;

	return y;
}

const static uint32 MATCH_MIN = 2;
const static uint32 MATCH_NICE_LENGTH = 64;
const static uint32 MATCH_SKIP_UPDATES_NICE_LENGTH = 7;
const static uint32 MATCH_NICE_RK_LENGTH = 256;
const static uint32 MATCH_MAX = MATCH_MIN + 255 + 7;

#define HASH4(x) ((x) * 987660757u)
#define VALUE4(p) (*(uint32 *)(p))
#define VALUE3(p) (VALUE4(p) & 0xFFFFFFu)
#define VALUE2(p) (*(uint16 *)(p))

uint32 get_match_min(uint32 dist);

struct MatchTable {
	uint16 max_len;
	uint32 delta[MATCH_MAX + 1];

	void Update(uint32 delta, uint16 len);
	void CarryFrom(MatchTable &prev, uint32 shift);
};

struct RingDictionary {
	byte *hist, *lookahead;
	uint32 hist_bits, hist_mask;
	uint32 hist_pos, lookahead_len;

	uint32 MatchLengthSigned(uint32 p0, uint32 p1, uint16 max_len, uint16 initial_len) const;
	uint32 MatchLength(uint32 p0, uint32 p1, uint16 max_len) const;
	int CharAtPosition(uint32 p) const;

	void Shift(uint32 shift);
};

struct MatchFinderHT {
	uint32 hash_bits, num_rows, window_bits, hash_shift, hash_mask, window_mask;
	uint32 *rows;

	uint32 Init(uint32 hash_bits, uint32 num_rows, uint32 window_bits);
	void Release();
	void FindAndUpdate(MatchTable &mt, uint32 h4, uint32 p, const RingDictionary &dict);
	void Shift(uint32 shift);
};

struct MatchFinderBT {
	const static int MAX_TESTS = 256;

	uint32 hash_bits, window_bits, hash_shift;
	uint32 *nodes, *heads, *tree;

	uint32 Init(uint32 hash_bits, uint32 window_bits);
	void Release();
	void FindAndUpdate(MatchTable &mt, uint32 h4, uint32 p, const RingDictionary &dict);
	void Shift(uint32 shift);
};

struct MatchFinderRK256 {
	const static uint32 BLOCK_BITS = 8;
	const static uint32 BLOCK_SIZE = 1 << BLOCK_BITS;
	const static uint32 BLOCK_MASK = BLOCK_SIZE - 1;

	const static uint32 ADDH = 0x2F0FD693u;

	//uint32 remh = 1; for (int i = 0; i < BLOCK_SIZE; i++) { remh *= addh; }
	const static uint32 REMH = 0x0E4EA401u;

	uint32 rolling_hash_add(uint32 p, int y) { return (y + p) * ADDH; }
	uint32 rolling_hash_add_remove(uint32 p, int yr, int yl) { return (yr + p - yl * REMH) * ADDH; }

	uint32 hash_shift, window_bits, hash_mask, window_mask;
	uint32 *table;

	uint32 carry_match_from, carry_match_to, carry_match_len;
	uint32 rh, rh_end;

	uint32 Init(uint32 hash_bits, uint32 window_bits);
	void Release();
	void FindAndUpdate(MatchTable &mt, uint32 p, const RingDictionary &dict);
	void Shift(uint32 shift);
};

uint32 get_match_min(uint32 dist) {
	uint32 match_min = MATCH_MIN;

	if (dist & ~0xFF) { ++match_min; }
	if (dist & ~0xFFF) { ++match_min; }
	if (dist & ~0xFFFFF) { ++match_min; }

	return match_min;
}

void MatchTable::CarryFrom(MatchTable &prev, uint32 shift) {
	if (prev.max_len <= shift) {
		max_len = 0;
	}
	else {
		max_len = prev.max_len - shift;
		for (uint16 i = 0; i <= max_len; i++) {
			delta[i] = prev.delta[i + shift];
		}
	}
}

void MatchTable::Update(uint32 mdelta, uint16 mlen) {
	ASSERT(mlen >= get_match_min(mdelta));
	ASSERT(mlen <= MATCH_MAX);
	ASSERT(mdelta > 0);

	int i = 0;
	while (i <= mlen && i <= max_len) {
		delta[i] = Min(mdelta, delta[i]);
		++i;
	}

	while (i <= mlen) {
		delta[i] = mdelta;
		++i;
	}

	max_len = Max(max_len, mlen);
}

uint32 RingDictionary::MatchLengthSigned(uint32 p0, uint32 p1, uint16 max_len, uint16 initial_len) const {
	ASSERT(p0 < p1);
	p0 += initial_len;
	p1 += initial_len;

	ASSERT((p0 >= hist_pos && p0 - hist_pos < lookahead_len) || hist_pos - p0 <= hist_mask);
	ASSERT((p1 >= hist_pos && p1 - hist_pos < lookahead_len) || hist_pos - p1 <= hist_mask);

	uint16 mlen = initial_len;
	while (mlen < max_len) {
		byte c0 = p0 >= hist_pos ? lookahead[p0 - hist_pos] : hist[p0 & hist_mask];
		byte c1 = p1 >= hist_pos ? lookahead[p1 - hist_pos] : hist[p1 & hist_mask];

		if (c0 != c1) {
			return mlen | ((c0 < c1) << 31);
		}

		++p0;
		++p1;
		++mlen;
	}

	return mlen;
}

uint32 RingDictionary::MatchLength(uint32 p0, uint32 p1, uint16 max_len) const {
	return RingDictionary::MatchLengthSigned(p0, p1, max_len, 0) & 0x7FFFFFFF;
}

int RingDictionary::CharAtPosition(uint32 p) const {
	ASSERT((p >= hist_pos && p - hist_pos < lookahead_len) || hist_pos - p <= hist_mask);
	return p >= hist_pos ? lookahead[p - hist_pos] : hist[p & hist_mask];
}

void RingDictionary::Shift(uint32 shift) {
	ASSERT(!(shift & hist_mask));
	hist_pos -= shift;
}

uint32 MatchFinderHT::Init(uint32 hash_bits, uint32 num_rows, uint32 window_bits) {
	this->hash_bits = hash_bits;
	this->hash_shift = 32 - hash_bits;
	this->window_bits = window_bits;
	this->num_rows = num_rows;
	this->window_mask = (1u << window_bits) - 1;
	this->hash_mask = (1u << (32 - window_bits)) - 1;

	rows = new uint32[num_rows << hash_bits];
	memset(rows, -1, (4U * num_rows) << hash_bits);
	return (4 * num_rows) << hash_bits;
}

void MatchFinderHT::Release() {
	delete[] rows;
}

void MatchFinderHT::FindAndUpdate(MatchTable &mt, uint32 h4, uint32 p, const RingDictionary &dict) {
	uint32 h_row = h4 & hash_mask;
	uint32 *ht = rows + (h4 >> hash_shift);
	uint32 carry = p | (h_row << window_bits);

	uint16 max_len = Min(dict.lookahead_len + dict.hist_pos - p, MATCH_MAX);

	uint32 best = MATCH_MIN - 1;
	uint32 best_p = -1;
	for (uint32 i = 0; i < num_rows; i++) {
		uint32 row = ht[i];

		if (best < max_len && (row >> window_bits) == h_row) {
			uint32 sp = row & window_mask;

			if (sp < p && p - sp <= dict.hist_mask) {
				uint16 mlen = dict.MatchLength(sp, p, max_len);

				if (mlen > best && mlen >= get_match_min(p - sp)) {
					mt.Update(p - sp, mlen);
					best = mlen;
				}
			}
		}

		ht[i] = carry;
		carry = row;
	}
}

void MatchFinderHT::Shift(uint32 shift) {
	uint32 *row = rows;
	uint32 *end = rows + (num_rows << hash_bits);
	while (row < end) {
		uint32 h = *rows >> window_bits;
		uint32 p = *rows & window_mask;

		if (p >= shift && *rows != uint32(-1)) {
			p -= shift;
			*rows = p | (h << window_bits);
		}
		else {
			*rows = -1;
		}

		++row;
	}
}

uint32 MatchFinderBT::Init(uint32 hash_bits, uint32 window_bits) {
	this->hash_bits = hash_bits;
	this->hash_shift = 32 - hash_bits;
	this->window_bits = window_bits;

	nodes = new uint32[(1u << hash_bits) + (2u << window_bits)];
	heads = nodes;
	tree = nodes + (1u << hash_bits);

	memset(heads, -1, 4u << hash_bits);
	memset(tree, -1, 8u << window_bits);

	return (4 << hash_bits) + (8 << window_bits);
}

void MatchFinderBT::Release() {
	delete[] nodes;
}

void MatchFinderBT::FindAndUpdate(MatchTable &mt, uint32 h4, uint32 p, const RingDictionary &dict) {
	uint32 *pending_left = tree + ((p & dict.hist_mask) << 1);
	uint32 *pending_right = pending_left + 1;
	uint32 left_len = 0, right_len = 0;

	uint32 sp = heads[h4 >> hash_shift];
	heads[h4 >> hash_shift] = p;
	ASSERT((h4 >> hash_shift) < (1u << (32 - hash_shift)));

	uint16 max_len = Min(dict.lookahead_len + dict.hist_pos - p, MATCH_MAX);
	uint16 tests_left = MAX_TESTS;
	while (sp != (uint32)-1 && p > sp && p - sp <= dict.hist_mask && tests_left-- > 0) {
		ASSERT(sp < p);

		uint32 *pair = tree + ((sp & dict.hist_mask) << 1);
		uint32 mlen_signed = dict.MatchLengthSigned(sp, p, max_len, Min(left_len, right_len));
		uint32 mlen = mlen_signed & 0x7FFFFFFF;

		if (mlen >= get_match_min(p - sp)) {
			mt.Update(p - sp, mlen);
		}

		if (mlen == max_len) {
			*pending_left = pair[0];
			*pending_right = pair[1];
			return;
		}

		if (mlen_signed >> 31) {
			*pending_left = sp;
			pending_left = pair + 1;
			sp = *pending_left;
			right_len = mlen;
		}
		else {
			*pending_right = sp;
			pending_right = pair;
			sp = *pending_right;
			left_len = mlen;
		}
	}

	*pending_right = -1;
	*pending_left = -1;
}

void MatchFinderBT::Shift(uint32 shift) {
	uint32 *n = heads;
	uint32 *end = heads + (1u << hash_bits) + (2u << window_bits);
	while (n < end) {
		*n = *n >= shift && *n != uint32(-1) ? *n - shift : -1;
		++n;
	}
}

uint32 MatchFinderRK256::Init(uint32 hash_bits, uint32 window_bits) {
	this->hash_shift = 32 - hash_bits;
	this->window_bits = window_bits;
	this->hash_mask = (1u << (32 - window_bits)) - 1;
	this->window_mask = (1u << window_bits) - 1;

	table = new uint32[1u << hash_bits];
	memset(table, -1, 4u << hash_bits);

	rh = 0;
	rh_end = 0;

	carry_match_len = 0;
	carry_match_from = carry_match_to = 0;

	return 4 << hash_bits;
}

void MatchFinderRK256::Release() {
	delete[] table;
}

void MatchFinderRK256::FindAndUpdate(MatchTable &mt, uint32 p, const RingDictionary &dict) {
	if (carry_match_len > 0) {
		if (p - carry_match_to < carry_match_len) {
			uint32 shift = p - carry_match_to;

			uint32 delta = carry_match_to - carry_match_from;
			uint32 mlen = carry_match_len - shift;
			if (mlen >= get_match_min(delta)) {
				mt.Update(delta, Min(mlen, MATCH_MAX));
			}
		}
		else {
			carry_match_len = 0;
		}
	}

	while (dict.lookahead_len >= (p - dict.hist_pos) + BLOCK_SIZE && rh_end < p + BLOCK_SIZE) {
		int c0 = dict.lookahead[rh_end - dict.hist_pos];

		if (rh_end >= BLOCK_SIZE) {
			int c1 = (rh_end - BLOCK_SIZE) >= dict.hist_pos ? dict.lookahead[rh_end - BLOCK_SIZE - dict.hist_pos] : dict.hist[(rh_end - BLOCK_SIZE) & dict.hist_mask];

			rh = rolling_hash_add_remove(rh, c0, c1);
		}
		else {
			rh = rolling_hash_add(rh, c0);
		}

		++rh_end;
		if (!(rh_end & BLOCK_MASK) && rh_end < p + BLOCK_SIZE) {
			uint32 h4 = rh >> hash_shift;
			table[h4] = p | (rh << window_bits);
		}
	}

	if (carry_match_len < MATCH_NICE_RK_LENGTH) {
		uint32 row_end = table[rh >> hash_shift];
		uint32 row_rh = row_end >> window_bits;
		uint32 sp = row_end & window_mask;

		if (row_rh == (rh & hash_mask) && sp < p && p - sp <= dict.hist_mask) {
			uint32 max_len = dict.lookahead_len + dict.hist_pos - p;
			uint32 mlen = dict.MatchLength(sp, p, max_len);

			if (mlen >= carry_match_len && mlen >= get_match_min(p - sp)) {
				mt.Update(p - sp, Min(mlen, MATCH_MAX));

				carry_match_from = sp;
				carry_match_to = p;
				carry_match_len = mlen;
			}
		}
	}

	if (!(rh_end & BLOCK_MASK) && rh_end == p + BLOCK_SIZE) {
		uint32 h4 = rh >> hash_shift;
		table[h4] = p | (rh << window_bits);
	}
}

void MatchFinderRK256::Shift(uint32 shift) {
	if (rh_end >= shift) {
		rh_end -= shift;
	}
	else {
		rh = 0;
		rh_end = 0;
	}
}

struct RepModel {
	uint32 table[4];

	void Init();
	void Add(uint32 delta);
	byte GetMatch(uint32 delta) const;
};

struct Model {
	const static int CMD_Literal = 0;
	const static int CMD_Dict = 1;
	const static int CMD_RepDelta = 2;
	//const static int CMD_RolzDict = 3;

	RepModel rep4;

	CDF2 cmd;

	CDF4 lit_hi, lit_lo[16];
	CDF4 len_ext_hi, len_ext_lo[16]; //, dist_lo;
	CDF3 len_direct, dist_slot_hi[4], dist_slot_lo[4][8];
};

struct ParseNode {
	uint16 link, len;
	uint32 cost, delta;
	byte cmd;
};

void RepModel::Init() {
	for (int i = 0; i < 4; i++) {
		table[i] = i + 1;
	}
}

void RepModel::Add(uint32 delta) {
	for (int i = 0; i < 4; i++) {
		if (table[i] == delta) {
			return;
		}
	}

	table[3] = table[2];
	table[2] = table[1];
	table[1] = table[0];
	table[0] = delta;
}

byte RepModel::GetMatch(uint32 delta) const {
	for (int i = 0; i < 4; i++) {
		if (table[i] == delta) {
			return i;
		}
	}

	return -1;
}

void model_init(Model &m) {
	m.rep4.Init();

	cdf_init(m.cmd);

	cdf_init(m.lit_hi);
	for (int i = 0; i < 16; i++) {
		cdf_init(m.lit_lo[i]);
	}

	cdf_init(m.len_direct);
	cdf_init(m.len_ext_hi);
	for (int i = 0; i < 16; i++) {
		cdf_init(m.len_ext_lo[i]);
	}

	//cdf_init(m.dist_lo);
	for (int c = 0; c < 4; c++) {
		cdf_init(m.dist_slot_hi[c]);
		for (int i = 0; i < 8; i++) {
			cdf_init(m.dist_slot_lo[c][i]);
		}
	}
}

uint32 model_cost_match(Model &m, uint32 delta, uint16 len) {
	ASSERT(len >= get_match_min(delta));
	ASSERT(delta > 0);

	uint32 cost = cdf_cost(m.cmd, Model::CMD_Dict);

	uint32 lv = len - get_match_min(delta);
	cost += cdf_cost(m.len_direct, Min(lv, 7u));

	uint32 lc = Min(lv, 3u);
	if (lv >= 7) {
		lv -= 7;
		int lo = lv & 0xF;
		int hi = lv >> 4;

		cost += cdf_cost(m.len_ext_hi, hi);
		cost += cdf_cost(m.len_ext_lo[hi], lo);
	}

	uint32 dv = delta - 1;
	if (dv >= 4) {
		int nb = clz32(dv) + 1;

		int add_bits = nb - 2;
		int top = dv >> add_bits;

		//if (add_bits > 4) {
		//	cost += (add_bits - 4) << LOG2_LUT_SCALE_BITS;
		//	cost += cdf_cost(m.dist_lo, dv & 0xF);
		//}
		//else {
		cost += add_bits << LOG2_LUT_SCALE_BITS;
		//}

		dv = ((nb - 1) << 1) + (top & 1);
	}

	uint32 dv_lo = dv & 0x7;
	uint32 dv_hi = dv >> 3;
	cost += cdf_cost(m.dist_slot_hi[lc], dv_hi);
	cost += cdf_cost(m.dist_slot_lo[lc][dv_hi], dv_lo);

	return cost;
}

uint32 model_cost_rep(Model &m, uint32 rep_idx, uint32 delta, uint16 len) {
	uint32 cost = cdf_cost(m.cmd, Model::CMD_RepDelta);

	uint32 lv = len - get_match_min(delta);
	cost += cdf_cost(m.len_direct, Min(lv, 7u));

	uint32 lc = Min(lv, 3u);
	if (lv >= 7) {
		lv -= 7;
		int lo = lv & 0xF;
		int hi = lv >> 4;

		cost += cdf_cost(m.len_ext_hi, hi);
		cost += cdf_cost(m.len_ext_lo[hi], lo);
	}

	cost += 2 << LOG2_LUT_SCALE_BITS;

	return cost;
}

void model_encode_match(CodeFrame &frame, Model &m, uint32 delta, uint16 len) {
	ASSERT(len >= get_match_min(delta));
	ASSERT(delta > 0);

	frame.WriteCDF(m.cmd, Model::CMD_Dict);
	cdf_update(m.cmd, Model::CMD_Dict);

	uint32 lv = len - get_match_min(delta);
	uint32 lc = Min(lv, 3u);

	frame.WriteCDF(m.len_direct, Min(lv, 7u));
	cdf_update(m.len_direct, Min(lv, 7u));

	if (lv >= 7) {
		lv -= 7;
		int lo = lv & 0xF;
		int hi = lv >> 4;

		frame.WriteCDF(m.len_ext_hi, hi);
		frame.WriteCDF(m.len_ext_lo[hi], lo);

		cdf_update(m.len_ext_hi, hi);
		cdf_update(m.len_ext_lo[hi], lo);
	}

	uint32 dv = delta - 1;
	if (dv < 4) {
		uint32 dv_lo = dv & 0x7;
		uint32 dv_hi = dv >> 3;

		frame.WriteCDF(m.dist_slot_hi[lc], dv_hi);
		frame.WriteCDF(m.dist_slot_lo[lc][dv_hi], dv_lo);

		cdf_update(m.dist_slot_hi[lc], dv_hi);
		cdf_update(m.dist_slot_lo[lc][dv_hi], dv_lo);
	}
	else {
		int nb = clz32(dv) + 1;

		int add_bits = nb - 2;
		int top = dv >> add_bits;
		uint32 add_mask = (1u << add_bits) - 1;
		uint32 dv_add = dv & add_mask;

		dv = ((nb - 1) << 1) + (top & 1);
		uint32 dv_lo = dv & 0x7;
		uint32 dv_hi = dv >> 3;

		frame.WriteCDF(m.dist_slot_hi[lc], dv_hi);
		frame.WriteCDF(m.dist_slot_lo[lc][dv_hi], dv_lo);

		cdf_update(m.dist_slot_hi[lc], dv_hi);
		cdf_update(m.dist_slot_lo[lc][dv_hi], dv_lo);

		if (add_bits < 4) {
			frame.WriteBits(dv_add, add_bits);
		}
		else {
			if (add_bits > 4) {
				frame.WriteBits(dv_add >> 4, add_bits - 4);
			}

			//dv_add &= 0xF;
			//frame.WriteCDF(m.dist_lo, dv_add);
			//cdf_update(m.dist_lo, dv_add);
			frame.WriteBits(dv_add & 0xF, 4);
		}
	}
}

void model_encode_rep(CodeFrame &frame, Model &m, byte rep_idx, uint16 len) {
	frame.WriteCDF(m.cmd, Model::CMD_RepDelta);
	cdf_update(m.cmd, Model::CMD_RepDelta);

	uint32 lv = len - get_match_min(m.rep4.table[rep_idx]);
	uint32 lc = Min(lv, 3u);

	frame.WriteCDF(m.len_direct, Min(lv, 7u));
	cdf_update(m.len_direct, Min(lv, 7u));

	if (lv >= 7) {
		lv -= 7;
		int lo = lv & 0xF;
		int hi = lv >> 4;

		frame.WriteCDF(m.len_ext_hi, hi);
		frame.WriteCDF(m.len_ext_lo[hi], lo);

		cdf_update(m.len_ext_hi, hi);
		cdf_update(m.len_ext_lo[hi], lo);
	}

	frame.WriteBits(rep_idx, 2);
}

uint32 model_decode_lv(DecodeFrame &frame, Model &m) {
	uint32 lv = frame.ReadCDF(m.len_direct);
	cdf_update(m.len_direct, lv);

	if (lv == 7) {
		int hi = frame.ReadCDF(m.len_ext_hi);
		int lo = frame.ReadCDF(m.len_ext_lo[hi]);

		cdf_update(m.len_ext_hi, hi);
		cdf_update(m.len_ext_lo[hi], lo);
		lv += (hi << 4) + lo;
	}

	return lv;
}

uint32 model_decode_dv(DecodeFrame &frame, Model &m, uint32 lv) {
	uint32 lc = Min(lv, 3u);

	uint32 dv_hi = frame.ReadCDF(m.dist_slot_hi[lc]);
	uint32 dv_lo = frame.ReadCDF(m.dist_slot_lo[lc][dv_hi]);
	uint32 dv = (dv_hi << 3) + dv_lo;

	cdf_update(m.dist_slot_hi[lc], dv_hi);
	cdf_update(m.dist_slot_lo[lc][dv_hi], dv_lo);

	if (dv >= 4) {
		uint32 add_bits = (dv >> 1) - 1;
		dv = (2 + (dv & 1)) << add_bits;

		if (add_bits < 4) {
			dv += frame.ReadBits(add_bits);
		}
		else {
			add_bits -= 4;
			if (add_bits > 0) {
				dv += frame.ReadBits(add_bits) << 4;
			}

			dv += frame.ReadBits(4);
			//uint32 dv_lo = frame.ReadCDF(m.dist_lo);
			//cdf_update(m.dist_lo, dv_lo);
			//dv += dv_lo;
		}
	}

	return dv;
}

uint32 model_cost_literal(Model &m, int y) {
	int lo = y & 0xF;
	int hi = y >> 4;

	return
		cdf_cost(m.cmd, Model::CMD_Literal) +
		cdf_cost(m.lit_hi, hi) +
		cdf_cost(m.lit_lo[hi], lo);
}

void model_encode_literal(CodeFrame &frame, Model &m, int y) {
	int lo = y & 0xF;
	int hi = y >> 4;

	frame.WriteCDF(m.cmd, 0);
	frame.WriteCDF(m.lit_hi, hi);
	frame.WriteCDF(m.lit_lo[hi], lo);

	cdf_update(m.cmd, Model::CMD_Literal);
	cdf_update(m.lit_hi, hi);
	cdf_update(m.lit_lo[hi], lo);
}

int model_decode_cmd(DecodeFrame &frame, Model &m) {
	int y = frame.ReadCDF(m.cmd);
	cdf_update(m.cmd, y);

	return y;
}

int model_decode_literal(DecodeFrame &frame, Model &m) {
	int hi = frame.ReadCDF(m.lit_hi);
	int lo = frame.ReadCDF(m.lit_lo[hi]);

	cdf_update(m.lit_hi, hi);
	cdf_update(m.lit_lo[hi], lo);

	return (hi << 4) + lo;
}

const uint32 PARSE_TABLE_SIZE = 1 << 12;

struct CarriedState {
	RepModel rep4;
};

uint32 parse_table(ParseNode(&table)[PARSE_TABLE_SIZE + 1], Model &m, RingDictionary &dict,
	MatchFinderHT &ht2, MatchFinderHT &ht3, MatchFinderBT &bt4, MatchFinderRK256 &rk, MatchTable &mt_carry,
	uint32 max_parse_len) {
	CarriedState carried_states[0x200];

	max_parse_len = Min(max_parse_len, PARSE_TABLE_SIZE);
	ASSERT(max_parse_len <= dict.lookahead_len);

	table[0].cost = 0;
	table[0].cmd = -1;
	table[0].len = 0;
	table[0].link = -1;
	carried_states[0].rep4 = m.rep4;

	table[1].cost = -1;
	table[1].cmd = Model::CMD_Literal;
	table[1].len = 0;
	table[1].link = 0;
	carried_states[1] = carried_states[0];

	MatchTable mt;
	uint32 p = 0, end_p = 1;
	while (p < end_p) {
		ASSERT(p < PARSE_TABLE_SIZE);
		uint32 np = 1 + p;

		int y = dict.lookahead[p];
		uint32 lit_cost = model_cost_literal(m, y);
		if (table[np].cost > table[p].cost + lit_cost) {
			table[np].cost = table[p].cost + lit_cost;
			table[np].cmd = Model::CMD_Literal;
			table[np].link = p;
			table[np].len = 0;

			carried_states[np & 0x1FF] = carried_states[p & 0x1FF];
		}

		mt.max_len = 0;
		mt.CarryFrom(mt_carry, 1);
		if (mt.max_len > 0 && dict.hist_pos + p >= mt.delta[mt.max_len]) {
			uint32 delta = mt.delta[mt.max_len];

			uint32 sp = dict.hist_pos + p - delta;
			while (mt.max_len < MATCH_MAX && dict.lookahead_len > mt.max_len + p &&
				dict.CharAtPosition(sp + mt.max_len) == dict.lookahead[p + mt.max_len]) {
				++mt.max_len;
				mt.delta[mt.max_len] = delta;
			}
		}

		if (mt.max_len < MATCH_NICE_LENGTH) {
			if (dict.lookahead_len >= 4 + p) {
				uint32 h2 = HASH4(VALUE2(dict.lookahead + p));
				uint32 h3 = HASH4(VALUE3(dict.lookahead + p));
				uint32 h4 = HASH4(VALUE4(dict.lookahead + p));

				ht2.FindAndUpdate(mt, h2, dict.hist_pos + p, dict);
				ht3.FindAndUpdate(mt, h3, dict.hist_pos + p, dict);
				bt4.FindAndUpdate(mt, h4, dict.hist_pos + p, dict);
			}

			if (dict.lookahead_len >= 256 + p) {
				rk.FindAndUpdate(mt, dict.hist_pos + p, dict);
			}
		}
		else if (!(p & MATCH_SKIP_UPDATES_NICE_LENGTH)) {
			if (dict.lookahead_len >= 4 + p) {
				uint32 h2 = HASH4(VALUE2(dict.lookahead + p));
				uint32 h3 = HASH4(VALUE3(dict.lookahead + p));

				ht2.FindAndUpdate(mt, h2, dict.hist_pos + p, dict);
				ht3.FindAndUpdate(mt, h3, dict.hist_pos + p, dict);
			}

			if (dict.lookahead_len >= 256 + p) {
				rk.FindAndUpdate(mt, dict.hist_pos + p, dict);
			}
		}

		mt_carry = mt;

		uint16 max_len = Min<uint32>(mt.max_len, max_parse_len - p);
		if (max_len < MATCH_MIN) {
			max_len = 0;
		}

		while (end_p < max_len + p) {
			++end_p;
			table[end_p].cost = -1;
			table[end_p].link = -1;
		}

		byte checked_rep4 = 0;

		uint16 tstep = (max_len - MATCH_MIN) >> 4;
		tstep += tstep == 0;
		for (uint16 tlen = max_len; tlen >= MATCH_MIN; tlen -= Min(tlen, tstep)) {
			uint16 min_len = get_match_min(mt.delta[tlen]);
			if (tlen < min_len) {
				continue;
			}

			uint32 np = tlen + p;
			uint32 match_cost = model_cost_match(m, mt.delta[tlen], tlen);
			if (table[np].cost > table[p].cost + match_cost) {
				table[np].cost = table[p].cost + match_cost;
				table[np].cmd = Model::CMD_Dict;
				table[np].link = p;
				table[np].len = tlen;
				table[np].delta = mt.delta[tlen];

				carried_states[np & 0x1FF] = carried_states[p & 0x1FF];
				carried_states[np & 0x1FF].rep4.Add(table[np].delta);
			}

			byte rep_idx = carried_states[p & 0x1FF].rep4.GetMatch(mt.delta[tlen]);
			if (rep_idx == byte(-1)) {
				continue;
			}
			checked_rep4 |= 1 << rep_idx;

			uint32 rep_cost = model_cost_rep(m, rep_idx, mt.delta[tlen], tlen);
			if (table[np].cost > table[p].cost + rep_cost) {
				table[np].cost = table[p].cost + rep_cost;
				table[np].cmd = Model::CMD_RepDelta;
				table[np].link = p;
				table[np].len = tlen;
				table[np].delta = rep_idx;

				carried_states[np & 0x1FF] = carried_states[p & 0x1FF];
				carried_states[np & 0x1FF].rep4.Add(mt.delta[tlen]);
			}
		}

		if (checked_rep4 != 15) {
			const RepModel &rep4 = carried_states[p & 0x1FF].rep4;
			for (byte rep_idx = 0; rep_idx < 4; ++rep_idx) {
				if (checked_rep4 & (1 << rep_idx) || rep4.table[rep_idx] >= dict.hist_pos + p) {
					continue;
				}

				uint16 mlen = dict.MatchLength(dict.hist_pos + p - rep4.table[rep_idx], dict.hist_pos + p, max_parse_len - p);
				mlen = Min<uint16>(mlen, MATCH_MAX);
				if (mlen >= get_match_min(rep4.table[rep_idx])) {
					while (end_p < mlen + p) {
						++end_p;
						table[end_p].cost = -1;
						table[end_p].link = -1;
					}

					uint32 rep_cost = model_cost_rep(m, rep_idx, rep4.table[rep_idx], mlen);
					uint32 np = mlen + p;
					if (table[np].cost > table[p].cost + rep_cost) {
						table[np].cost = table[p].cost + rep_cost;
						table[np].cmd = Model::CMD_RepDelta;
						table[np].link = p;
						table[np].len = mlen;
						table[np].delta = rep_idx;

						carried_states[np & 0x1FF] = carried_states[p & 0x1FF];
						carried_states[np & 0x1FF].rep4.Add(rep4.table[rep_idx]);
					}
				}
			}
		}

		++p;
	}

	ParseNode running_node, tmp_node;
	memset(&running_node, -1, sizeof(running_node));

	uint16 running_end = -1;
	uint16 cur = p;
	while (cur != (uint16)-1) {
		uint16 prev = table[cur].link;

		tmp_node = table[cur];
		table[cur] = running_node;
		table[cur].link = running_end;
		running_node = tmp_node;

		running_end = cur;
		cur = prev;
	}
	ASSERT(running_end == 0);
	return end_p;
}

//void RolzDict::Init() {
//	memset(cells, -1, sizeof(cells));
//	memset(head, 0, sizeof(head));
//}
//
//void RolzDict::Add(int ctx, uint32 pos, uint16 len) {
//	if (len < 3) {
//		return;
//	}
//
//	RolzEntry &e = cells[ctx][head[ctx]++ & TABLE_MASK];
//
//	e.pos = pos;
//	e.len = len;
//}
//
//RolzEntryMatch RolzDict::Lookup(int ctx, uint32 pos) {
//	RolzEntryMatch m;
//	m.depth = -1;
//
//	uint32 lo = 0;
//	uint32 hi = TABLE_SIZE;
//
//	while (lo < hi) {
//		uint32 mid = lo + (hi - lo) / 2;
//		const RolzEntry &e = cells[ctx][(head[ctx] - mid - 1) & TABLE_MASK];
//		if (e.pos == pos) {
//			m.depth = mid;
//			m.e = e;
//			return m;
//		}
//		else if (pos > e.pos || e.pos == uint32(-1)) {
//			hi = mid;
//		}
//		else {
//			lo = mid + 1;
//		}
//	}
//
//	return m;
//}

char print_fill_table[81];
int print_fill_max = 0;

void print_fill(int n) {
	print_fill_max = Max(print_fill_max, n);

	if (n < print_fill_max) {
		int fill = print_fill_max - n;
		fill = Min(fill, 80);
		printf("%s\r", print_fill_table + (80 - fill));
	}
	else {
		printf("\r");
	}
}

void encode_file(FILE *fin, FILE *fout, uint32 hist_bits) {
	fseek(fin, 0, SEEK_END);
	uint64 flen = _fpos64(fin);
	fseek(fin, 0, SEEK_SET);

	while (hist_bits > 10 && flen < (1ull << (hist_bits - 1))) {
		--hist_bits;
	}

	const uint32 window_size = 1u << hist_bits;

	const uint32 frame_bits = Clamp(hist_bits - 2, 14U, 17U);
	const uint32 frame_size = 1u << frame_bits;
	const uint32 chunk_size = ((frame_size * 15) / 16) - 0x200;
	const uint32 chunk_feed_size = chunk_size + MATCH_MAX + 1;

	uint32 *freqs_buf = new uint32[4 * frame_size];

	RingDictionary dict;
	dict.hist_bits = hist_bits;
	dict.hist_mask = (1u << dict.hist_bits) - 1;
	dict.hist = new byte[1u << dict.hist_bits];
	dict.lookahead = nullptr;
	dict.hist_pos = 0;
	dict.lookahead_len = 0;

	byte *chunk_mem = new byte[chunk_feed_size];
	byte *frame_mem = new byte[frame_size];

	Model model;
	model_init(model);

	ParseNode parse_nodes[PARSE_TABLE_SIZE + 1];

	MatchFinderHT ht2, ht3;
	MatchFinderBT bt4;
	MatchFinderRK256 rk;

	uint32 mf_mem = 0;
	mf_mem += ht2.Init(12, 1, dict.hist_bits);
	mf_mem += ht3.Init(12 + Clamp(dict.hist_bits, 15u, 17u) - 15u, 2, dict.hist_bits);
	mf_mem += bt4.Init(13 + Clamp(dict.hist_bits, 16u, 20u) - 16u, dict.hist_bits);
	mf_mem += rk.Init(15 + Clamp(dict.hist_bits, 16u, 22u) - 16u, dict.hist_bits);

	printf("Model: %d KB\n", (sizeof(model) + 1023) >> 10);
	printf("Parser: %d KB\n", (sizeof(parse_nodes) + 1023) >> 10);
	printf("Dictionary: %d KB\n", (dict.hist_mask + 1 + 1023) >> 10);
	printf("Frame: %d KB\n", ((1u << frame_bits) + 1023) >> 10);
	printf("Dictionary search: %d KB\n", (mf_mem + 1023) >> 10);
	printf("Working...\r");

	uint32 write_pos = 0;
	frame_mem[write_pos++] = byte(hist_bits >> 8);
	frame_mem[write_pos++] = byte(hist_bits);
	frame_mem[write_pos++] = byte(frame_bits >> 8);
	frame_mem[write_pos++] = byte(frame_bits);

	CodeFrame frame;
	MatchTable mt_carry;
	mt_carry.max_len = 0;

	uint64 file_pos_reader = 0, file_pos_writer = 0, print_pos = 0;

	uint32 chunk_read = fread(chunk_mem, 1, chunk_feed_size, fin);
	file_pos_reader += chunk_read;

	uint32 crc32 = 0;
	crc32 = crc32_calc(chunk_mem, chunk_read, crc32);

	clock_t t0 = clock();

	while (chunk_read > 0) {
		uint32 p_end = Min(chunk_size, chunk_read);
		frame.Init(frame_mem + write_pos, frame_mem + frame_size, freqs_buf, 4 * frame_size);

		if (dict.hist_pos >= 2 * window_size) {
			dict.Shift(window_size);
			ht2.Shift(window_size);
			ht3.Shift(window_size);
			bt4.Shift(window_size);
			rk.Shift(window_size);
		}

		uint32 parse_start = 0, parse_end = 0;
		uint32 p = 0;
		while (p < p_end) {
			dict.lookahead = chunk_mem + p;
			dict.lookahead_len = Min(chunk_feed_size, chunk_read - p);

			ASSERT(p <= parse_end);
			if (p == parse_end) {
				uint32 parse_len = parse_table(parse_nodes, model, dict, ht2, ht3, bt4, rk, mt_carry, p_end - p);
				ASSERT(parse_len > 0);
				ASSERT(p + parse_len <= p_end);
				parse_start = p;
				parse_end = p + parse_len;
			}

			ParseNode pn = parse_nodes[p - parse_start];
			if (pn.cmd == Model::CMD_Literal) {
				dict.hist[dict.hist_pos++ & dict.hist_mask] = *dict.lookahead;
				model_encode_literal(frame, model, *dict.lookahead);
				++p;
			}
			else if (pn.cmd == Model::CMD_Dict) {
				model_encode_match(frame, model, pn.delta, pn.len);
				ASSERT(dict.hist_pos >= pn.delta);

				model.rep4.Add(pn.delta);
				p += pn.len;

				while (pn.len-- > 0) {
					int y = dict.hist[(dict.hist_pos - pn.delta) & dict.hist_mask];
					ASSERT(y == *dict.lookahead);

					dict.hist[dict.hist_pos++ & dict.hist_mask] = *dict.lookahead++;
				}
			}
			else if (pn.cmd == Model::CMD_RepDelta) {
				model_encode_rep(frame, model, pn.delta, pn.len);
				uint32 delta = model.rep4.table[pn.delta];
				ASSERT(dict.hist_pos >= delta);

				model.rep4.Add(delta);
				p += pn.len;

				while (pn.len-- > 0) {
					int y = dict.hist[(dict.hist_pos - delta) & dict.hist_mask];
					ASSERT(y == *dict.lookahead);

					dict.hist[dict.hist_pos++ & dict.hist_mask] = *dict.lookahead++;
				}
			}
			else {
				ASSERT(false);
			}
		}

		ASSERT(!frame.NeedsFlush());
		uint32 frame_written = frame.Flush();
		write_pos += frame_written;

		fwrite(frame_mem, 1, write_pos, fout);
		file_pos_writer += write_pos;
		write_pos = 0;

		if ((file_pos_writer - print_pos) >> 17) {
			clock_t t1 = (clock() - t0) / CLOCKS_PER_SEC;
			if (t1) {
				uint32 est_time_left = Max(2u, uint32(((flen - file_pos_reader) * t1) / file_pos_reader));
				print_fill(printf("Working... %" PRIu64 " -> %" PRIu64 " ~%d seconds left", file_pos_reader, file_pos_writer, est_time_left));
			}
			else {
				print_fill(printf("Working... %" PRIu64 " -> %" PRIu64, file_pos_reader, file_pos_writer));
			}

			print_pos = file_pos_writer;
		}

		if (p_end < chunk_read) {
			memmove(chunk_mem, chunk_mem + p_end, chunk_read - p_end);
			chunk_read -= p_end;

			uint32 read_add = fread(chunk_mem + chunk_feed_size - p_end, 1, p_end, fin);
			crc32 = crc32_calc(chunk_mem + chunk_feed_size - p_end, read_add, crc32);

			file_pos_reader += read_add;
			chunk_read += read_add;
		}
		else {
			chunk_read = fread(chunk_mem, 1, chunk_feed_size, fin);
			crc32 = crc32_calc(chunk_mem, chunk_read, crc32);

			file_pos_reader += chunk_read;
		}
	}

	fwrite(frame_mem, 1, write_pos, fout);
	file_pos_writer += write_pos;

	frame_mem[0] = 0;
	frame_mem[1] = 0;
	frame_mem[2] = 0;
	frame_mem[3] = 0;
	fwrite(frame_mem, 1, 4, fout);
	file_pos_writer += 4;

	print_fill(printf("Working... %" PRIu64 " -> %" PRIu64, file_pos_reader, file_pos_writer));
	printf("\nDone (input CRC32 %X, %.2f sec)\n", crc32, (clock() - t0) / double(CLOCKS_PER_SEC));

	ht2.Release();
	ht3.Release();
	bt4.Release();
	rk.Release();

	delete[] chunk_mem;
	delete[] frame_mem;
	delete[] dict.hist;
	delete[] freqs_buf;
}

void decode_file(FILE *fin, FILE *fout) {
	byte header[4];
	ASSERT(fread(header, 1, 4, fin) == 4);
	uint32 hist_bits = (header[0] << 8) + header[1];
	uint32 frame_bits = (header[2] << 8) + header[3];

	ASSERT(hist_bits >= 12);
	ASSERT(hist_bits <= 28);
	ASSERT(frame_bits >= 12);
	ASSERT(frame_bits <= 20);

	RingDictionary dict;
	dict.hist_bits = hist_bits;
	dict.hist_mask = (1u << dict.hist_bits) - 1;
	dict.hist = new byte[1u << dict.hist_bits];
	dict.hist_pos = 0;

	const uint32 window_size = 1 << hist_bits;

	DecodeFrame frame;
	Model model;
	model_init(model);

	printf("Model: %d KB\n", (sizeof(model) + 1023) >> 10);
	printf("Dictionary: %d KB\n", (dict.hist_mask + 1 + 1023) >> 10);
	printf("Frame: %d KB\n", ((1u << frame_bits) + 1023) >> 10);

	clock_t t0 = clock();

	byte write_buf[0x1200];
	uint32 write_pos = 0;

	uint64 file_pos_reader = 0, file_pos_writer = 0, print_pos = 0;

	const uint32 frame_size = 1u << frame_bits;
	byte *frame_mem = new byte[frame_size];
	uint32 frame_read = fread(frame_mem, 1, frame_size, fin);
	file_pos_reader += frame_read;

	print_fill(printf("Working... %" PRIu64 " -> %" PRIu64, file_pos_reader, file_pos_writer));

	uint32 crc32 = 0;

	while (true) {
		uint32 consumed_frame_size = frame.Init(frame_mem, frame_mem + frame_read);
		if (consumed_frame_size == uint32(-1)) {
			break;
		}

		ASSERT(consumed_frame_size <= frame_read);

		if (dict.hist_pos >= 2 * window_size) {
			dict.Shift(window_size);
		}

		while (frame.num_ops > 0) {
			if (write_pos & ~0xFFF) {
				file_pos_writer += write_pos;
				if (fout) {
					fwrite(write_buf, 1, write_pos, fout);
				}
				crc32 = crc32_calc(write_buf, write_pos, crc32);
				write_pos = 0;
			}

			int cmd = model_decode_cmd(frame, model);
			if (cmd == Model::CMD_Literal) {
				int y = model_decode_literal(frame, model);

				dict.hist[dict.hist_pos++ & dict.hist_mask] = y;
				write_buf[write_pos++] = y;
			}
			else if (cmd == Model::CMD_Dict) {
				uint32 lv = model_decode_lv(frame, model);
				uint32 dv = model_decode_dv(frame, model, lv);
				++dv;
				lv += get_match_min(dv);

				model.rep4.Add(dv);
				ASSERT(dict.hist_pos >= dv);
				while (lv-- > 0) {
					int y = dict.hist[(dict.hist_pos - dv) & dict.hist_mask];
					dict.hist[dict.hist_pos++ & dict.hist_mask] = y;
					write_buf[write_pos++] = y;
				}
			}
			else if (cmd == Model::CMD_RepDelta) {
				byte rep_idx = frame.ReadBits(2);
				uint32 lv = model_decode_lv(frame, model);
				uint32 dv = model.rep4.table[rep_idx];
				lv += get_match_min(dv);

				model.rep4.Add(dv);
				ASSERT(dict.hist_pos >= dv);
				while (lv-- > 0) {
					int y = dict.hist[(dict.hist_pos - dv) & dict.hist_mask];
					dict.hist[dict.hist_pos++ & dict.hist_mask] = y;
					write_buf[write_pos++] = y;
				}
			}
		}

		memmove(frame_mem, frame_mem + consumed_frame_size, frame_read - consumed_frame_size);
		frame_read -= consumed_frame_size;
		uint32 frame_add = fread(frame_mem + frame_read, 1, frame_size - frame_read, fin);
		file_pos_reader += frame_add;
		frame_read += frame_add;

		if ((file_pos_writer - print_pos) >> 23) {
			print_fill(printf("Working... %" PRIu64 " -> %" PRIu64, file_pos_reader, file_pos_writer));
			print_pos = file_pos_writer;
		}
	}

	if (write_pos) {
		file_pos_writer += write_pos;
		if (fout) {
			fwrite(write_buf, 1, write_pos, fout);
		}
		crc32 = crc32_calc(write_buf, write_pos, crc32);
	}

	print_fill(printf("Working... %" PRIu64 " -> %" PRIu64, file_pos_reader, file_pos_writer));
	printf("\nDone (output CRC32 %X, %.2f sec)\n", crc32, (clock() - t0) / double(CLOCKS_PER_SEC));

	delete[] frame_mem;
	delete[] dict.hist;
}

int _tolower(int c) { return c | 0x20; }

void _tolower(char *v) {
	while (*v) {
		*v = _tolower(*v);
		v++;
	}
}

int main(int argc, char **argv) {
	printf("NLZM 1.03 - Written by Nauful\n");

#ifdef MSVC
#	ifdef DEBUG
	int dbg_flag = _CRTDBG_REPORT_FLAG | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF;
#	else
	int dbg_flag = _CRTDBG_REPORT_FLAG | _CRTDBG_LEAK_CHECK_DF;
#	endif
	_CrtSetDbgFlag(dbg_flag);
#endif

	crc32_init();
	log2_init();
	cdf_init_mixins();

	for (int i = 0; i < 80; i++) {
		print_fill_table[i] = ' ';
}
	print_fill_table[80] = '\0';

	uint32 hist_bits = 22;

	FILE *fin = nullptr, *fout = nullptr;
	while (argc >= 2 && *argv[1] == '-') {
		char *arg = argv[1];
		argv++;
		argc--;

		_tolower(arg);
		while (*arg == '-') {
			++arg;
		}

		if (!strncmp(arg, "window:", 7)) {
			hist_bits = Clamp<uint32>(atoi(arg + 7), 15u, 28u);
			printf("Window bits: %d\n", hist_bits);
		}
		else {
			printf("Unrecognized flag %s\n", arg);
			return -1;
		}
	}

	if (argc == 4 && _tolower(argv[1][0]) == 'c') {
		if ((fout = fopen(argv[3], "rb"))) {
			printf("Error: %s already exists\n", argv[3]);
			fclose(fout);
			return -1;
		}

		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		fout = fopen(argv[3], "wb");
		if (!fout) {
			printf("Error: %s file does not exist\n", argv[3]);
			fclose(fin);
			return -1;
		}

		encode_file(fin, fout, hist_bits);

		fclose(fin);
		fclose(fout);
	}
	else if (argc == 4 && _tolower(argv[1][0]) == 'd') {
		if ((fout = fopen(argv[3], "rb"))) {
			printf("Error: %s already exists\n", argv[3]);
			fclose(fout);
			return -1;
		}

		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		fout = fopen(argv[3], "wb");
		if (!fout) {
			printf("Error: %s file does not exist\n", argv[3]);
			fclose(fin);
			return -1;
		}

		decode_file(fin, fout);
		fclose(fin);
		fclose(fout);
	}
	else if (argc == 3 && _tolower(argv[1][0]) == 't') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		decode_file(fin, nullptr);
		fclose(fin);
	}
	else if (argc == 3 && _tolower(argv[1][0]) == 'h') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		uint32 crc32 = crc32_file(fin);
		fclose(fin);
		printf("%X\n", crc32);
	}
	else {
		printf("Commands:\n"
			"\t[flags] c [input] [output] - Compress input file to output file (best parser)\n"
			"\td [input] [output] - Decompress input file to output file\n"
			"\tt [input] - Decompress input file in memory\n"
			"\th [input] - Calculate CRC32 for input file\n"
			"Flags:\n"
			"\t-window:bits = Maximum window size in bits, default 22 (4 MB), min 15, max 28 (32 KB to 256 MB)\n");
	}

#ifdef MSVC
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}
