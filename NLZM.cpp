// NLZM 1.02 64-bit - Written by Nauful
// Released into the public domain.
// Please drop a comment if you find this useful.

#ifndef _CRT_DISABLE_PERFCRIT_LOCKS
#	define _CRT_DISABLE_PERFCRIT_LOCKS
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
#	define _CRT_SECURE_NO_WARNINGS
#endif

#include <ctime>
#include <intrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef HAS_MAIN
#	ifdef _MSC_VER
#		define MSVC
#	endif
#endif

#ifdef MSVC
# define ERROR(x) { printf("\nError: %s\n", x); __debugbreak(); exit(-1); }
#else
# define ERROR(x) { printf("\nError: %s\n", x); exit(-1); }
#endif

#ifdef MSVC
# define ASSERT(x) { if (!(x)) { printf("\nAssert failed!\nCondition: %s\nLine %d\n", #x, __LINE__); __debugbreak(); exit(-1); } }
#else
# define ASSERT(x) { if (!(x)) { printf("\nAssert failed!\nCondition: %s\nLine %d\n", #x, __LINE__); exit(-1); } }
#endif

#ifdef MSVC
# define snprintf _snprintf
#	ifndef ALIGN
#		define ALIGN(x) __declspec(align(x))
#	endif
#	ifndef RESTRICT
#		define RESTRICT __restrict
#	endif
#	ifndef INLINE
#		define INLINE __forceinline
#	endif
#else
#	ifndef ALIGN
#		define ALIGN(x) __attribute__((aligned(x)))
#	endif
#	ifndef RESTRICT
#		define RESTRICT __restrict
#	endif
#	ifndef INLINE
#		define INLINE __attribute__((always_inline))
#	endif
#endif

#ifdef DEBUG
# define ASSERT_DEBUG(x) ASSERT(x)
# define ASSERT_MATCH(x) ASSERT(x)
#else
# define ASSERT_DEBUG(x)
# define ASSERT_MATCH(x)
#endif

typedef unsigned int uint;
typedef unsigned char byte;

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

typedef char int8;
typedef short int16;
typedef int int32;
typedef long long int64;

template<typename T, size_t N> constexpr size_t CountOf(T const (&)[N]) { return N; }
template <typename T> constexpr T Min(const T a, const T b) { return (a < b) ? a : b; }
template <typename T> constexpr T Max(const T a, const T b) { return (a < b) ? b : a; }

template <typename T> constexpr T Clamp(const T v, const T low, const T high) {
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

INLINE uint32 clz32(uint32 x) {
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

INLINE uint32 popcnt32(uint32 x) {
	x -= ((x >> 1) & 0x55555555);
	x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
	x = (((x >> 4) + x) & 0x0F0F0F0F);
	x += (x >> 8);
	x += (x >> 16);

	return x & 0x0000003F;
}

INLINE uint32 ctz32(uint32 x) {
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
		uint32 acc = 0;

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

const static int COUNTER_SCALE_BITS = 12;
const static int COUNTER_SCALE = 1 << COUNTER_SCALE_BITS;
const static int COUNTER_MASK = COUNTER_SCALE - 1;
const static int COUNTER_HSCALE_BITS = COUNTER_SCALE_BITS - 1;
const static int COUNTER_HSCALE = 1 << COUNTER_HSCALE_BITS;

const static int COUNTER_ADAPT1_SPEED = 5;
const static int COUNTER_ADAPT3_SPEED = 6;
const static int COUNTER_ADAPT4_SPEED = 7;

const static int COUNTER_MAX_ADAPT1_SPEED = 1 << COUNTER_ADAPT1_SPEED;
const static int COUNTER_MAX_ADAPT3_SPEED = 1 << COUNTER_ADAPT3_SPEED;
const static int COUNTER_MAX_ADAPT4_SPEED = 1 << COUNTER_ADAPT4_SPEED;

namespace RANS {
	typedef uint32 rans_t;

	const static rans_t MID = 1 << 16;

	void rans_enc_flush(rans_t r, byte** pptr) {
		uint32 x = r;

		pptr[0] -= 4;
		pptr[0][0] = (byte)(x);
		pptr[0][1] = (byte)(x >> 8);
		pptr[0][2] = (byte)(x >> 16);
		pptr[0][3] = (byte)(x >> 24);
	}

	rans_t rans_dec_init(byte** pptr) {
		uint32 x = pptr[0][0];
		x |= pptr[0][1] << 8;
		x |= pptr[0][2] << 16;
		x |= pptr[0][3] << 24;
		pptr[0] += 4;

		return x;
	}

	rans_t rans_enc_put(rans_t x, byte** pptr, uint32 start, uint32 freq) {
		uint32 x_max = ((MID >> COUNTER_SCALE_BITS) << 16);
		x_max *= freq;

		if (x >= x_max) {
			*--(pptr[0]) = (byte)x;
			*--(pptr[0]) = (byte)(x >> 8);
			x >>= 16;
		}

		return ((x / freq) << COUNTER_SCALE_BITS) + (x % freq) + start;
	}

	INLINE rans_t rans_dec_consume(rans_t x, uint32 start, uint32 freq) {
		return freq * (x >> COUNTER_SCALE_BITS) + (x & COUNTER_MASK) - start;
	}

	INLINE rans_t rans_dec_normalize(rans_t x, byte** pptr) {
		if (x < MID) {
			x = (x << 16) + (pptr[0][0] << 8) + pptr[0][1];
			*pptr += 2;
		}

		return x;
	}

	static int16 mixin16[16][16];
	static int16 mixin9[9][16];
	static int16 mixin4[5][8];

	template<int num_syms, int rate, int scale_max, int dim_out>
	void init_mixin_table(int16(&mixin)[num_syms][dim_out]) {
		const int mixin_bias = (1 << rate) - 1 - num_syms;
		ASSERT(mixin_bias > 0);

		for (int y = 0; y < num_syms; y++) {
			for (int x = 0; x <= y; x++) {
				mixin[y][x] = x;
			}

			for (int x = y + 1; x < num_syms; x++) {
				mixin[y][x] = scale_max + x + mixin_bias;
			}

			for (int x = num_syms; x < dim_out; x++) {
				mixin[y][x] = scale_max;
			}
		}
	}

	template<int num_syms, int scale_max>
	void init_cdf(int16(&cdf)[num_syms]) {
		for (int i = 0; i <= num_syms; i++) {
			cdf[i] = (i << scale_bits) / num_syms;
		}
	}

	void init() {
		init_mixin_table<5, COUNTER_ADAPT1_SPEED, COUNTER_SCALE, 8>(mixin4);
		init_mixin_table<9, COUNTER_ADAPT3_SPEED, COUNTER_SCALE, 16>(mixin9);
		init_mixin_table<16, COUNTER_ADAPT4_SPEED, COUNTER_SCALE, 16>(mixin16);
	}
}

struct BitDecoder {
	const static int BLOCK_MAX_SIZE = 0x4000;
	const static int HEADER_SIZE = 20;

	uint32 bit_code_mask[32];

	byte buf[BLOCK_MAX_SIZE];
	uint32 buf_n, buf_offset;

	byte* rptr_bitcode;
	byte* rptr_probcode;
	uint32 bit_stream;
	int bit_stream_bits;

	int bitcode_total_bits_encoded;
	int rans_num_buffered, rans_total_bytes_written, rans_index;

	FILE* f;

	RANS::rans_t rs[4];

	BitDecoder(FILE* f) : f(f), buf_n(0) {}

	void Clear() {
		rptr_bitcode = nullptr;
		rptr_probcode = nullptr;

		bit_stream = 0;
		bit_stream_bits = 0;
		bitcode_total_bits_encoded = 0;

		rans_num_buffered = 0;
		rans_total_bytes_written = 0;
	}

	void ReadBlock() {
		Clear();

		if (buf_offset > 0) {
			if (buf_offset < BLOCK_MAX_SIZE) {
				memmove(buf, buf + buf_offset, BLOCK_MAX_SIZE - buf_offset);
				buf_offset = BLOCK_MAX_SIZE - buf_offset;
			}
			else {
				buf_offset = 0;
			}

			buf_n = fread(buf + buf_offset, 1, BLOCK_MAX_SIZE - buf_offset, f);
			buf_n += buf_offset;
		}
		else {
			buf_n = fread(buf, 1, BLOCK_MAX_SIZE, f);
		}

		buf_offset = 0;
		ASSERT(buf_n >= 4);

		int block_bytes = (buf[0] << 24) + (buf[1] << 16) + (buf[2] << 8) + buf[3];
		if (block_bytes == 0) {
			return;
		}

		ASSERT(block_bytes >= HEADER_SIZE + 4);
		ASSERT(buf_n >= block_bytes);

		bitcode_total_bits_encoded = (buf[4] << 24) + (buf[5] << 16) + (buf[6] << 8) + buf[7];
		int bitcode_total_bytes_written = (buf[8] << 24) + (buf[9] << 16) + (buf[10] << 8) + buf[11];
		rans_num_buffered = (buf[12] << 24) + (buf[13] << 16) + (buf[14] << 8) + buf[15];
		rans_total_bytes_written = (buf[16] << 24) + (buf[17] << 16) + (buf[18] << 8) + buf[19];

		rptr_bitcode = buf + HEADER_SIZE;
		rptr_probcode = buf + HEADER_SIZE + bitcode_total_bytes_written;

		while (bit_stream_bits < 24) {
			bit_stream = (bit_stream << 8) + *rptr_bitcode++;
			bit_stream_bits += 8;
		}

		rs[0] = RANS::rans_dec_init(&rptr_probcode);
		rs[1] = RANS::rans_dec_init(&rptr_probcode);
		rs[2] = RANS::rans_dec_init(&rptr_probcode);
		rs[3] = RANS::rans_dec_init(&rptr_probcode);

		rans_index = 0;
		buf_offset = block_bytes;
	}

	void Init() {
		for (int i = 0; i < 32; i++) {
			bit_code_mask[i] = (1 << i) - 1;
		}

		buf_offset = 0;

		Clear();
		ReadBlock();
	}

	INLINE uint16 DecodeSlot() {
		if (rans_num_buffered <= rans_index) {
			ASSERT(rans_num_buffered == rans_index);
			ReadBlock();
		}

		return rs[rans_index & 3] & COUNTER_MASK;
	}

	INLINE RANS::rans_t* GetSlot() {
		if (rans_num_buffered <= rans_index) {
			ASSERT(rans_num_buffered == rans_index);
			ReadBlock();
		}

		return rs + (rans_index & 3);
	}

	INLINE void DecodeRange(RANS::rans_t* r, uint16 start, uint16 freq) {
		RANS::rans_t x = RANS::rans_dec_consume(*r, start, freq);
		*r = RANS::rans_dec_normalize(x, &rptr_probcode);

		++rans_index;
	}

	uint32 DecodeBits(int num_bits) {
		ASSERT_DEBUG(num_bits <= 24);

		if (bitcode_total_bits_encoded <= 0) {
			ASSERT(bitcode_total_bits_encoded == 0);
			ReadBlock();
		}

		ASSERT_DEBUG(num_bits <= bit_stream_bits);
		uint32 v = bit_stream >> (bit_stream_bits - num_bits);

		bit_stream_bits -= num_bits;
		bit_stream &= bit_code_mask[bit_stream_bits];

		int num_refill = (31 - bit_stream_bits) >> 3;
		int num_refill_bits = num_refill << 3;

		uint32 c24 = *reinterpret_cast<uint32*>(rptr_bitcode);
		c24 = ((c24 & 0xFF0000) >> 16) | (c24 & 0xFF00) | ((c24 & 0xFF) << 16);
		bit_stream = bit_stream << num_refill_bits;
		bit_stream |= c24 >> (24 - num_refill_bits);
		bit_stream_bits |= num_refill_bits;
		rptr_bitcode += num_refill;

		bitcode_total_bits_encoded -= num_bits;

		return v;
	}
};

struct BitCoder {
	const static int WRITE_BLOCK_BYTES = ((15 * BitDecoder::BLOCK_MAX_SIZE) / 16) - 256;
	const static int WRITE_BLOCK_BITS = WRITE_BLOCK_BYTES << 3;

	const static int RANS_BUFFER_SIZE = BitDecoder::BLOCK_MAX_SIZE * 8;

	byte buf[BitDecoder::BLOCK_MAX_SIZE];
	uint32* rans_buffer;

	byte* wptr_bitcode;
	uint32 bit_stream;
	int bit_stream_bits;
	int bitcode_total_bits_encoded, bitcode_total_bytes_written;
	int rans_est_bits_encoded, rans_num_buffered;

	FILE* f;

	BitCoder(FILE* f) : f(f) {
		rans_buffer = new uint32[RANS_BUFFER_SIZE];
	}

	void Flush() {
		WriteBlock();

		int buf_terminal_size[4] = { 0 };
		fwrite(buf_terminal_size, 1, 32, f);

		delete[] rans_buffer;
	}

	void Reset() {
		wptr_bitcode = buf + BitDecoder::HEADER_SIZE;

		bit_stream = 0;
		bit_stream_bits = 0;
		bitcode_total_bits_encoded = 0;
		bitcode_total_bytes_written = 0;

		rans_est_bits_encoded = 0;
		rans_num_buffered = 0;
	}

	void Init() {
		Reset();
	}

	void WriteBlock() {
		while (bit_stream_bits > 0) {
			*wptr_bitcode++ = byte(bit_stream >> 24);
			bit_stream <<= 8;
			bit_stream_bits -= 8;
			++bitcode_total_bytes_written;
		}

		byte* wptr_probcode = buf + BitDecoder::BLOCK_MAX_SIZE - 1;
		RANS::rans_t rs[4] = { RANS::MID, RANS::MID, RANS::MID, RANS::MID };
		for (int i = rans_num_buffered - 1; i >= 0; i--) {
			rs[i & 3] = RANS::rans_enc_put(rs[i & 3], &wptr_probcode, rans_buffer[i] & 0xFFFF, rans_buffer[i] >> 16);
			ASSERT(wptr_probcode > wptr_bitcode);
		}

		RANS::rans_enc_flush(rs[3], &wptr_probcode);
		RANS::rans_enc_flush(rs[2], &wptr_probcode);
		RANS::rans_enc_flush(rs[1], &wptr_probcode);
		RANS::rans_enc_flush(rs[0], &wptr_probcode);
		ASSERT(wptr_probcode > wptr_bitcode);

		int rans_total_bytes_written = (buf + BitDecoder::BLOCK_MAX_SIZE) - wptr_probcode;
		memmove(buf + BitDecoder::HEADER_SIZE + bitcode_total_bytes_written,
			buf + BitDecoder::BLOCK_MAX_SIZE - rans_total_bytes_written,
			rans_total_bytes_written);

		int block_bytes = BitDecoder::HEADER_SIZE + bitcode_total_bytes_written + rans_total_bytes_written;
		ASSERT(block_bytes <= BitDecoder::BLOCK_MAX_SIZE);

		buf[0] = block_bytes >> 24;
		buf[1] = block_bytes >> 16;
		buf[2] = block_bytes >> 8;
		buf[3] = block_bytes;

		buf[4] = bitcode_total_bits_encoded >> 24;
		buf[5] = bitcode_total_bits_encoded >> 16;
		buf[6] = bitcode_total_bits_encoded >> 8;
		buf[7] = bitcode_total_bits_encoded;

		buf[8] = bitcode_total_bytes_written >> 24;
		buf[9] = bitcode_total_bytes_written >> 16;
		buf[10] = bitcode_total_bytes_written >> 8;
		buf[11] = bitcode_total_bytes_written;

		buf[12] = rans_num_buffered >> 24;
		buf[13] = rans_num_buffered >> 16;
		buf[14] = rans_num_buffered >> 8;
		buf[15] = rans_num_buffered;

		buf[16] = rans_total_bytes_written >> 24;
		buf[17] = rans_total_bytes_written >> 16;
		buf[18] = rans_total_bytes_written >> 8;
		buf[19] = rans_total_bytes_written;

		fwrite(buf, 1, block_bytes, f);
		Reset();
	}

	bool NeedWrite() const {
		return (bitcode_total_bytes_written << 3) + (rans_est_bits_encoded >> LOG2_LUT_SCALE_BITS) >= WRITE_BLOCK_BITS ||
			rans_num_buffered == RANS_BUFFER_SIZE;
	}

	void EncodeRange(uint16 start, uint16 freq) {
		ASSERT(freq > 0);
		ASSERT(start + freq <= COUNTER_SCALE);

		rans_est_bits_encoded += log2_lut[freq >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
		rans_buffer[rans_num_buffered++] = (freq << 16) + start;

		if (NeedWrite()) {
			WriteBlock();
		}
	}

	void EncodeBits(uint32 v, int num_bits) {
		ASSERT(num_bits <= 24);
		ASSERT(bit_stream_bits + num_bits <= 32);

		int lsh = 32 - bit_stream_bits - num_bits;
		ASSERT_DEBUG(lsh >= 0);

		bit_stream += v << lsh;
		bit_stream_bits += num_bits;
		bitcode_total_bits_encoded += num_bits;

		while (bit_stream_bits & ~7) {
			*wptr_bitcode++ = byte(bit_stream >> 24);
			bit_stream <<= 8;
			bit_stream_bits -= 8;
			++bitcode_total_bytes_written;
		}

		if (NeedWrite()) {
			WriteBlock();
		}
	}
};

const static int MATCH_MIN_BASE = 2;
const static int MATCH_LIMIT = MATCH_MIN_BASE + (1 << 12);

INLINE uint32 get_match_min(uint32 dist) {
	uint32 match_min = 2;
	if (dist & ~0xFF) { ++match_min; }
	if (dist & ~0xFFFF) { ++match_min; }
	if (dist & ~0xFFFFFF) { ++match_min; }

	return match_min;
}

template<int rate>
static void cdf_update8(int16(&cdf)[9], int16* updates) {
#ifdef USE_SSE
	__m128i _cdf0 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf));
	__m128i _upd0 = _mm_load_si128(reinterpret_cast<const __m128i*>(updates));
	_upd0 = _mm_srai_epi16(_mm_sub_epi16(_upd0, _cdf0), rate);
	_cdf0 = _mm_add_epi16(_cdf0, _upd0);
	_mm_store_si128(reinterpret_cast<__m128i*>(cdf), _cdf0);
#else
	for (int i = 0; i < 8; i++) {
		cdf[i] += (updates[i] - cdf[i]) >> rate;
	}
#endif
}

static int cdf_lookup8(const int16(&cdf)[9], uint16 f) {
#ifdef USE_SSE
	__m128i _cdf0 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf));

	__m128i _f = _mm_cvtsi32_si128(f);
	_f = _mm_shuffle_epi32(_mm_unpacklo_epi16(_f, _f), 0);
	__m128i _cmp0 = _mm_cmpgt_epi16(_cdf0, _f);
	int idx_mask = _mm_movemask_epi8(_cmp0) | 0x10000;

	int idx = ctz32(idx_mask) >> 1;
	return idx - 1;
#else
	int r = 0;

	if (f >= cdf[4]) { r += 4; }
	if (f >= cdf[2 + r]) { r += 2; }
	if (f >= cdf[1 + r]) { r += 1; }

	return r;
#endif
}

template<int rate>
static void cdf_update16(int16(&cdf)[17], int16* updates) {
#ifdef USE_SSE
	__m128i _cdf0 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf));
	__m128i _cdf1 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf + 8));

	__m128i _upd0 = _mm_load_si128(reinterpret_cast<const __m128i*>(updates));
	__m128i _upd1 = _mm_load_si128(reinterpret_cast<const __m128i*>(updates + 8));

	_upd0 = _mm_srai_epi16(_mm_sub_epi16(_upd0, _cdf0), rate);
	_upd1 = _mm_srai_epi16(_mm_sub_epi16(_upd1, _cdf1), rate);

	_cdf0 = _mm_add_epi16(_cdf0, _upd0);
	_cdf1 = _mm_add_epi16(_cdf1, _upd1);

	_mm_store_si128(reinterpret_cast<__m128i*>(cdf), _cdf0);
	_mm_store_si128(reinterpret_cast<__m128i*>(cdf + 8), _cdf1);
#else
	for (int i = 0; i < 16; i++) {
		cdf[i] += (updates[i] - cdf[i]) >> rate;
	}
#endif
}

static int cdf_lookup16(const int16(&cdf)[17], uint16 f) {
#ifdef USE_SSE
	__m128i _cdf0 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf));
	__m128i _cdf1 = _mm_load_si128(reinterpret_cast<const __m128i*>(cdf + 8));

	__m128i _f = _mm_cvtsi32_si128(f);
	_f = _mm_shuffle_epi32(_mm_unpacklo_epi16(_f, _f), 0);

	__m128i _cmp0 = _mm_cmpgt_epi16(_cdf0, _f);
	__m128i _cmp1 = _mm_cmpgt_epi16(_cdf1, _f);
	int idx_mask = _mm_movemask_epi8(_mm_packs_epi16(_cmp0, _cmp1)) | 0x10000;

	int idx = ctz32(idx_mask);
	return idx - 1;
#else
	int r = 0;

	if (f >= cdf[8]) { r += 8; }
	if (f >= cdf[4 + r]) { r += 4; }
	if (f >= cdf[2 + r]) { r += 2; }
	if (f >= cdf[1 + r]) { r += 1; }

	return r;
#endif
}

const static int16 mix1[2] = { (1 << COUNTER_SCALE_BITS) - (1 << COUNTER_ADAPT1_SPEED) - 1, (1 << COUNTER_ADAPT1_SPEED) + 1 };
const static int16 def_cdf1[3] = { 0, 1 << (COUNTER_SCALE_BITS - 1), COUNTER_SCALE };

const static int16 def_cdf2[9] = {
	0 << (COUNTER_SCALE_BITS - 2), 1 << (COUNTER_SCALE_BITS - 2), 2 << (COUNTER_SCALE_BITS - 2), 3 << (COUNTER_SCALE_BITS - 2),
	4 << (COUNTER_SCALE_BITS - 2), COUNTER_SCALE, COUNTER_SCALE, COUNTER_SCALE, COUNTER_SCALE
};

const static int16 def_cdf8[9] = {
	0 << (COUNTER_SCALE_BITS - 3), 1 << (COUNTER_SCALE_BITS - 3), 2 << (COUNTER_SCALE_BITS - 3), 3 << (COUNTER_SCALE_BITS - 3),
	4 << (COUNTER_SCALE_BITS - 3), 5 << (COUNTER_SCALE_BITS - 3), 6 << (COUNTER_SCALE_BITS - 3), 7 << (COUNTER_SCALE_BITS - 3),
	8 << (COUNTER_SCALE_BITS - 3)
};

const static int16 def_cdf16[17] = {
	0 << (COUNTER_SCALE_BITS - 4), 1 << (COUNTER_SCALE_BITS - 4), 2 << (COUNTER_SCALE_BITS - 4), 3 << (COUNTER_SCALE_BITS - 4),
	4 << (COUNTER_SCALE_BITS - 4), 5 << (COUNTER_SCALE_BITS - 4), 6 << (COUNTER_SCALE_BITS - 4), 7 << (COUNTER_SCALE_BITS - 4),
	8 << (COUNTER_SCALE_BITS - 4), 9 << (COUNTER_SCALE_BITS - 4), 10 << (COUNTER_SCALE_BITS - 4), 11 << (COUNTER_SCALE_BITS - 4),
	12 << (COUNTER_SCALE_BITS - 4), 13 << (COUNTER_SCALE_BITS - 4), 14 << (COUNTER_SCALE_BITS - 4), 15 << (COUNTER_SCALE_BITS - 4),
	16 << (COUNTER_SCALE_BITS - 4)
};

struct Bit1 {
	int16 cdf[3];

	void Init() {
		for (int i = 0; i < 3; i++) {
			cdf[i] = def_cdf1[i];
		}
	}

	INLINE uint16 Price(int y) {
		return log2_lut[(cdf[y + 1] - cdf[y]) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
	}

	void Encode(int y, BitCoder& bc) {
		bc.EncodeRange(cdf[y], cdf[y + 1] - cdf[y]);
		cdf[1] += (mix1[y] - cdf[1]) >> COUNTER_ADAPT1_SPEED;
	}

	int Decode(BitDecoder& bc) {
		RANS::rans_t* rs = bc.GetSlot();

		int y = (*rs & COUNTER_MASK) >= cdf[1];
		bc.DecodeRange(rs, cdf[y], cdf[y + 1] - cdf[y]);

		cdf[1] += (mix1[y] - cdf[1]) >> COUNTER_ADAPT1_SPEED;
		return y;
	}
};

struct BitNibble2 {
	ALIGN(16) int16 cdf[9];
	byte weight;

	void Init() {
		weight = 1;
		for (int i = 0; i < 9; i++) {
			cdf[i] = def_cdf2[i];
		}
	}

	INLINE byte Weight() const { return weight; }
	INLINE uint16 Price(int y) { return log2_lut[(cdf[y + 1] - cdf[y]) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)]; }

	void Encode(int y, BitCoder& bc) {
		bc.EncodeRange(cdf[y], cdf[y + 1] - cdf[y]);
		cdf_update8<COUNTER_ADAPT1_SPEED>(cdf, RANS::mixin4[y]);
		weight += weight < COUNTER_MAX_ADAPT1_SPEED;
	}

	int Decode(BitDecoder& bc) {
		RANS::rans_t* rs = bc.GetSlot();

		int y = cdf_lookup8(cdf, *rs & COUNTER_MASK);
		bc.DecodeRange(rs, cdf[y], cdf[y + 1] - cdf[y]);

		cdf_update8<COUNTER_ADAPT1_SPEED>(cdf, RANS::mixin4[y]);
		weight += weight < COUNTER_MAX_ADAPT1_SPEED;
		return y;
	}
};

struct BitNibble3 {
	ALIGN(16) int16 cdf[9];

	void Init() {
		for (int i = 0; i < 9; i++) {
			cdf[i] = def_cdf8[i];
		}
	}

	INLINE uint16 Price(int y) {
		return log2_lut[(cdf[y + 1] - cdf[y]) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
	}

	void Encode(int y, BitCoder& bc) {
		bc.EncodeRange(cdf[y], cdf[y + 1] - cdf[y]);
		cdf_update8<COUNTER_ADAPT3_SPEED>(cdf, RANS::mixin9[y]);
	}

	int Decode(BitDecoder& bc) {
		RANS::rans_t* rs = bc.GetSlot();

		int y = cdf_lookup8(cdf, *rs & COUNTER_MASK);
		bc.DecodeRange(rs, cdf[y], cdf[y + 1] - cdf[y]);

		cdf_update8<COUNTER_ADAPT3_SPEED>(cdf, RANS::mixin9[y]);
		return y;
	}

	void EncodeReverse(int y, BitCoder& bc) { Encode(y, bc); }
	int DecodeReverse(BitDecoder& bc) { return Decode(bc); }
};

struct BitNibble4 {
	ALIGN(16) int16 cdf[17];

	void Init() {
		for (int i = 0; i < 17; i++) {
			cdf[i] = def_cdf16[i];
		}
	}

	INLINE uint32 Price(int y) {
		return log2_lut[(cdf[y + 1] - cdf[y]) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
	}

	void Encode(int y, BitCoder& bc) {
		bc.EncodeRange(cdf[y], cdf[y + 1] - cdf[y]);
		cdf_update16<COUNTER_ADAPT4_SPEED>(cdf, RANS::mixin16[y]);
	}

	int Decode(BitDecoder& bc) {
		RANS::rans_t* rs = bc.GetSlot();

		int y = cdf_lookup16(cdf, *rs & COUNTER_MASK);
		bc.DecodeRange(rs, cdf[y], cdf[y + 1] - cdf[y]);

		cdf_update16<COUNTER_ADAPT4_SPEED>(cdf, RANS::mixin16[y]);
		return y;
	}

	void EncodeReverse(int y, BitCoder& bc) { Encode(y, bc); }
	int DecodeReverse(BitDecoder& bc) { return Decode(bc); }
};

struct BitNibble8 {
	BitNibble4 th, tl[16];

	void Init() {
		th.Init();
		for (int i = 0; i < 16; i++) {
			tl[i].Init();
		}
	}

	uint32 Price(int y) {
		int nh = y >> 4;
		int nl = y & 15;

		return th.Price(nh) + tl[nh].Price(nl);
	}

	void Encode(int y, BitCoder& bc) {
		int nh = y >> 4;
		int nl = y & 15;

		th.Encode(nh, bc);
		tl[nh].Encode(nl, bc);
	}

	int Decode(BitDecoder& bc) {
		int nh = th.Decode(bc);
		int nl = tl[nh].Decode(bc);

		return (nh << 4) + nl;
	}
};

struct BitNibble6 {
	BitNibble3 th, tl[8];

	void Init() {
		th.Init();
		for (int i = 0; i < 8; i++) {
			tl[i].Init();
		}
	}

	uint32 Price(int y) {
		int nh = y >> 3;
		int nl = y & 7;

		return th.Price(nh) + tl[nh].Price(nl);
	}

	void Encode(int y, BitCoder& bc) {
		int nh = y >> 3;
		int nl = y & 7;

		th.Encode(nh, bc);
		tl[nh].Encode(nl, bc);
	}

	int Decode(BitDecoder& bc) {
		int nh = th.Decode(bc);
		int nl = tl[nh].Decode(bc);

		return (nh << 3) + nl;
	}
};

const static int LITERAL_MAX_POSITION_BITS = 4;
const static int LITERAL_MAX_POSITION_SIZE = 1 << LITERAL_MAX_POSITION_BITS;

const static int CMD_LITERAL = 0;
const static int CMD_DICT_MATCH = 1;
const static int CMD_REP4 = 2;
//const static int CMD_ROLZ_MATCH = 3;

struct RepState {
	int entry[4];

	void Init() {
		entry[0] = 1;
		entry[1] = 2;
		entry[2] = 3;
		entry[3] = 4;
	}

	void Add(int dist) {
		int i = 0;
		while (i < 4) {
			if (entry[i] == dist) {
				break;
			}

			++i;
		}

		if (i < 4) {
			for (int j = i - 1; j >= 0; j--) {
				entry[j + 1] = entry[j];
			}

			entry[0] = dist;
		}
		else {
			entry[3] = entry[2];
			entry[2] = entry[1];
			entry[1] = entry[0];
			entry[0] = dist;
		}
	}
};

//const static int ROLZ_BITS = 9;
//const static int ROLZ_SIZE = 1 << ROLZ_BITS;
//const static int ROLZ_MASK = ROLZ_SIZE - 1;

struct State {
	BitNibble8* lit;
	BitNibble2 cmd[4][LITERAL_MAX_POSITION_SIZE];

	BitNibble2 len_marker;
	BitNibble3 len_low[LITERAL_MAX_POSITION_SIZE];
	BitNibble3 len_mid[LITERAL_MAX_POSITION_SIZE];
	BitNibble8 len_high;
	BitNibble6 dist_slots[4];
	BitNibble4 dist_low;

	RepState rep;
	BitNibble2 rep_dist;

	//RolzState rolz;

	byte match_hist;

	uint32 cached_price_cmd[LITERAL_MAX_POSITION_SIZE][4];
	uint32* cached_price_lit;
	uint32 cached_price_length[LITERAL_MAX_POSITION_SIZE][273];
	uint32 cached_price_dist_slot[4][64];
	uint32 cached_price_dist_low[16];
	uint32 cached_price_rep_dist[4];

	int pos_align_bits, context_bits;

	void Init(int pos_align_bits, int context_bits) {
		this->pos_align_bits = pos_align_bits;
		this->context_bits = context_bits;

		ASSERT(pos_align_bits <= 4);
		ASSERT(context_bits <= 8);
		ASSERT(pos_align_bits + context_bits <= 8);

		lit = new BitNibble8[1 << (pos_align_bits + context_bits)];
		cached_price_lit = new uint32[(1 << pos_align_bits) << (8 + context_bits)];

		rep.Init();
		match_hist = 0;

		int pos_align_size = 1 << pos_align_bits;
		int context_size = 1 << context_bits;

		for (int i = 0; i < pos_align_size; i++) {
			for (int j = 0; j < context_size; j++) {
				lit[(i << context_bits) + j].Init();
			}
		}

		for (int s = 0; s < 4; s++) {
			for (int i = 0; i < pos_align_size; i++) {
				cmd[s][i].Init();
			}
		}

		len_marker.Init();
		for (int i = 0; i < pos_align_size; i++) { len_low[i].Init(); }
		for (int i = 0; i < pos_align_size; i++) { len_mid[i].Init(); }
		len_high.Init();

		for (int i = 0; i < 4; i++) { dist_slots[i].Init(); }
		dist_low.Init();

		rep_dist.Init();

		//rolz.Init();

		printf("State: %d KB\n", uint32(sizeof(State) +
			((sizeof(BitNibble8) << pos_align_bits) << context_bits) +
			((sizeof(uint32) << pos_align_bits) << (8 + context_bits)) + 1023) >> 10);
	}

	void Release() {
		//rolz.Release();
		delete[] lit;
		delete[] cached_price_lit;
	}

	void ShiftLeft(int shift) {
		//rolz.ShiftLeft(shift);
	}

	void NextStateMatch() { match_hist = ((match_hist << 1) + 1) & 3; }
	void NextStateLiteral() { match_hist = (match_hist << 1) & 3; }

	void UpdatePriceCache() {
		int pos_align_size = 1 << pos_align_bits;

		for (int pb = 0; pb < pos_align_size; ++pb) {
			for (int c = 0; c < 4; c++) {
				cached_price_cmd[pb][c] = 0;
				uint32 total_weight = 0;

				for (int _match_hist = 0; _match_hist < 4; _match_hist++) {
					uint32 w = cmd[_match_hist][pb].Weight();

					total_weight += w;
					cached_price_cmd[pb][c] += cmd[_match_hist][pb].Price(c) * w;
				}

				cached_price_cmd[pb][c] /= total_weight;
			}
		}

		memset(cached_price_lit, -1, sizeof(uint32) * (1 << pos_align_bits) << (8 + context_bits));

		uint32 lp0 = len_marker.Price(0);
		uint32 lp1 = len_marker.Price(1);
		uint32 lp2 = len_marker.Price(2);
		uint32 lp3 = len_marker.Price(3);
		for (int lv = 0; lv < 273; lv++) {
			for (int pb = 0; pb < pos_align_size; ++pb) {
				if (lv < 8) {
					cached_price_length[pb][lv] = lp0 + len_low[pb].Price(lv);
				}
				else if (lv < 16) {
					cached_price_length[pb][lv] = lp1 + len_mid[pb].Price(lv - 8);
				}
				else if (lv < 272) {
					cached_price_length[pb][lv] = lp2 + len_high.Price(lv - 16);
				}
				else {
					cached_price_length[pb][lv] = lp3 + (12 << LOG2_LUT_SCALE_BITS);
				}
			}
		}

		for (int slot = 0; slot < 4; slot++) {
			for (int y = 0; y < 64; y++) {
				cached_price_dist_slot[slot][y] = dist_slots[slot].Price(y);
			}
		}

		for (int i = 0; i < 16; i++) {
			cached_price_dist_low[i] = dist_low.Price(i);
		}

		for (int i = 0; i < 4; i++) {
			cached_price_rep_dist[i] = rep_dist.Price(i);
		}
	}

	INLINE int ToLitIndex(int pb, int ctx) {
		return (pb << context_bits) + (ctx >> (8 - context_bits));
	}

	uint32 PriceLiteralCached(int y, int pb, int ctx) {
		int lit_idx = ToLitIndex(pb, ctx);
		uint32& lit_price = cached_price_lit[(lit_idx << 8) + y];
		if (lit_price == uint32(-1)) {
			lit_price = lit[lit_idx].Price(y);
		}

		return cached_price_cmd[pb][CMD_LITERAL] + lit_price;
	}


	uint32 PriceRepCached(uint32 lv, int rep_idx, int pb, int ctx) {
		uint32 sum = cached_price_cmd[pb][CMD_REP4];
		sum += cached_price_length[pb][Min(lv, 272U)];
		sum += cached_price_rep_dist[rep_idx];

		return sum;
	}

	uint32 PriceMatchCached(uint32 lv, uint32 dv, int pb, int ctx) {
		ASSERT(dv < 0x7FFFFFFF);
		uint32 sum = cached_price_cmd[pb][CMD_DICT_MATCH];
		sum += cached_price_length[pb][Min(lv, 272U)];

		int dist_idx = lv < 4 ? lv : 3;
		if (dv < 4) {
			sum += cached_price_dist_slot[dist_idx][dv];
		}
		else {
			int nb = clz32(dv) + 1;

			int add_bits = nb - 2;
			int top = dv >> add_bits;
			int add_mask = (1 << add_bits) - 1;
			int slot = ((nb - 1) << 1) + (top & 1);

			sum += cached_price_dist_slot[dist_idx][slot];
			if (add_bits < 4) {
				sum += add_bits << LOG2_LUT_SCALE_BITS;
			}
			else {
				sum += (add_bits - 4) << LOG2_LUT_SCALE_BITS;
				sum += cached_price_dist_low[dv & 15];
			}
		}

		return sum;
	}
};

const static int MATCH_MAX_RESULTS = 256;
const static int MATCH_MAX_RESULT_ROWS = 2 * MATCH_MAX_RESULTS;

struct MatchFinder_Base {
	uint32 match_results[2 * MATCH_MAX_RESULT_ROWS + 1];

	virtual uint32 MinLength() const = 0;
	virtual void ShiftLeft(int shift) = 0;
	virtual uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) = 0;
	virtual void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) = 0;
	virtual uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) = 0;
};

static INLINE uint8 value8(const byte* b) { return *b; }
static INLINE uint16 value16(const byte* b) { return *reinterpret_cast<const uint16*>(b); }
static INLINE uint32 value24(const byte* b) { return *reinterpret_cast<const uint32*>(b) & 0xFFFFFFu; }
static INLINE uint32 value32(const byte* b) { return *reinterpret_cast<const uint32*>(b); }
static INLINE uint64 value64(const byte* b) { return *reinterpret_cast<const uint64*>(b); }

static INLINE uint32 hash4(uint32 x) { return x * 123456791u; }

static int match_finder_padding = 16;

INLINE static uint try_match(const byte* buf, uint32 mp, uint32 p, uint32 p_end) {
	uint len = 0;

#ifdef MSVC
	unsigned long r = 0;
	const uint64* p0 = reinterpret_cast<const uint64*>(buf + p);
	const uint64* p1 = reinterpret_cast<const uint64*>(buf + mp);
	while (p + len < p_end && len < MATCH_LIMIT) {
		uint64 v = *p0 ^ *p1;
		if (v) {
			_BitScanForward64(&r, v);
			len += r >> 3;
			if (p + len > p_end) {
				len = p_end - p;
			}

			if (len >= MATCH_LIMIT) {
				return MATCH_LIMIT - 1;
			}

			return len;
		}
		++p0;
		++p1;
		len += 8;
	}

	if (p + len > p_end) {
		len = p_end - p;
	}

#else
	while (p + len + 4 < p_end && len < MATCH_LIMIT && *reinterpret_cast<const uint32*>(buf + p + len) == *reinterpret_cast<const uint32*>(buf + mp + len)) {
		len += 4;
	}

	while (p + len < p_end && len < MATCH_LIMIT && buf[p + len] == buf[mp + len]) {
		++len;
	}
#endif

	if (len >= MATCH_LIMIT) {
		return MATCH_LIMIT - 1;
	}

	return len;
}

struct HR2 : public MatchFinder_Base {
	int hash_bits;
	uint32* table;
	uint32* cache;

	HR2(int hash_bits)
		: hash_bits(hash_bits) {
		table = nullptr;
		cache = nullptr;
	}

	HR2() : HR2(0) {}

	uint32 MinLength() const override { return 2; }

	void Init() {
		table = new uint32[2 << hash_bits];
		cache = new uint32[1 << hash_bits];

		const int table_mem = sizeof(uint32) * 2 << hash_bits;
		const int cache_mem = sizeof(uint32) << hash_bits;
		memset(table, -1, table_mem);
		memset(cache, -1, cache_mem);
		printf("HR2: %d KB\n", (table_mem + cache_mem + 1023) >> 10);
	}

	void Release() {
		delete[] table;
		delete[] cache;
	}

	void ShiftLeft(int shift) override {
		uint32* row = table;
		uint32* rows_end = table + (2 << hash_bits);

		while (row < rows_end) {
			*row = (*row < shift || *row == uint32(-1)) ? uint32(-1) : (*row - shift);
			++row;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return;
		}

		uint16 v2 = value16(buf + p);
		uint32 row_base_idx = hash4(v2) >> (32 - hash_bits);
		uint32* row = table + (2 * row_base_idx);
		uint32* cache_line = cache + row_base_idx;

		*cache_line = (*cache_line << 16) | v2;
		row[1] = row[0];
		row[0] = p;
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		int best_len = MinLength() - 1;
		int num_results = 0;

		uint16 v2 = value16(buf + p);
		uint32 row_base_idx = hash4(v2) >> (32 - hash_bits);
		uint32* row = table + (2 * row_base_idx);
		uint32* cache_line = cache + row_base_idx;

		if ((*cache_line & 0xFFFF) == v2 && row[0] != uint32(-1)) {
			if (p - row[0] < max_dist) {
				match_results[num_results++] = 2;
				match_results[num_results++] = row[0];
			}
		}
		else if (((*cache_line >> 16) & 0xFFFF) == v2 && row[1] != uint32(-1)) {
			if (p - row[1] < max_dist) {
				match_results[num_results++] = 2;
				match_results[num_results++] = row[1];
			}
		}

		*cache_line = (*cache_line << 16) | v2;
		row[1] = row[0];
		row[0] = p;

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		int best_len = MinLength() - 1;
		int num_results = 0;

		uint16 v2 = value16(buf + p);
		uint32 row_base_idx = hash4(v2) >> (32 - hash_bits);
		const uint32* row = table + (2 * row_base_idx);
		const uint32* cache_line = cache + row_base_idx;

		if ((*cache_line & 0xFFFF) == v2 && row[0] != uint32(-1)) {
			if (p - row[0] <= max_dist) {
				match_results[num_results++] = 2;
				match_results[num_results++] = row[0];
			}
		}
		else if (((*cache_line >> 16) & 0xFFFF) == v2 && row[1] != uint32(-1)) {
			if (p - row[1] <= max_dist) {
				match_results[num_results++] = 2;
				match_results[num_results++] = row[1];
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}
};

struct HR3 : public MatchFinder_Base {
	int hash_bits, row_width;
	uint32* table;
	uint32* cache;
	byte* heads;

	HR3(int hash_bits, int row_width)
		: hash_bits(hash_bits), row_width(row_width) {
		table = nullptr;
		cache = nullptr;
		heads = nullptr;
	}

	HR3() : HR3(0, 0) {}

	uint32 MinLength() const override { return 3; }

	void Init() {
		table = new uint32[row_width << hash_bits];
		cache = new uint32[row_width << hash_bits];
		heads = new byte[1 << hash_bits];

		const int table_mem = sizeof(uint32) * row_width << hash_bits;
		memset(table, -1, table_mem);
		memset(cache, -1, table_mem);
		memset(heads, 0, 1 << hash_bits);
		printf("HT3: %d KB\n", (2 * table_mem + (1 << hash_bits) + 1023) >> 10);
	}

	void Release() {
		delete[] table;
		delete[] cache;
		delete[] heads;
	}

	void ShiftLeft(int shift) override {
		uint32* row = table;
		uint32* rows_end = table + (row_width << hash_bits);

		while (row < rows_end) {
			*row = (*row < shift || *row == uint32(-1)) ? uint32(-1) : (*row - shift);
			++row;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		if (p_end - p < MinLength()) {
			return;
		}

		uint32 v3 = value24(buf + p);
		uint32 row_base_idx = hash4(v3) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;
		uint32* cache_row = cache + row_base_idx * row_width;

		byte& head = heads[row_base_idx];
		if (++head >= row_width) {
			head -= row_width;
		}

		row[head] = p;
		cache_row[head] = v3;
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		uint32 v3 = value24(buf + p);
		uint32 row_base_idx = hash4(v3) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;
		uint32* cache_row = cache + row_base_idx * row_width;

		byte& head = heads[row_base_idx];
		for (int ridx = 0; ridx < row_width; ridx++) {
			int ri = head - ridx;
			if (ri < 0) {
				ri += row_width;
			}

			uint32 mp = row[ri];
			uint32 v = cache_row[ri];

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (v == v3 && buf[p + best_len] == buf[mp + best_len]) {
				uint len = 3 + try_match(buf, mp + 3, p + 3, p_end);

				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		if (++head >= row_width) {
			head -= row_width;
		}

		row[head] = p;
		cache_row[head] = v3;

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		uint32 v3 = value24(buf + p);
		uint32 row_base_idx = hash4(v3) >> (32 - hash_bits);
		const uint32* row = table + row_base_idx * row_width;
		const uint32* cache_row = cache + row_base_idx * row_width;

		const byte head = heads[row_base_idx];
		for (int ridx = 0; ridx < row_width; ridx++) {
			int ri = head - ridx;
			if (ri < 0) {
				ri += row_width;
			}

			uint32 mp = row[ri];
			uint32 v = cache_row[ri];
			uint32 eh = hash4(v) >> (32 - hash_bits);

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (v == v3 && buf[p + best_len] == buf[mp + best_len]) {
				uint len = 3 + try_match(buf, mp + 3, p + 3, p_end);

				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}
};

struct HR4 : public MatchFinder_Base {
	int hash_bits, row_width;
	uint32* table;
	uint32* cache;
	byte* heads;

	HR4(int hash_bits, int row_width)
		: hash_bits(hash_bits), row_width(row_width) {
		table = nullptr;
		cache = nullptr;
		heads = nullptr;
	}

	HR4() : HR4(0, 0) {}

	uint32 MinLength() const override { return 4; }

	void Init() {
		table = new uint32[row_width << hash_bits];
		cache = new uint32[row_width << hash_bits];
		heads = new byte[1 << hash_bits];

		const int table_mem = sizeof(uint32) * row_width << hash_bits;
		memset(table, -1, table_mem);
		memset(cache, -1, table_mem);
		memset(heads, 0, 1 << hash_bits);
		printf("HR4: %d KB\n", (2 * table_mem + (1 << hash_bits) + 1023) >> 10);
	}

	void Release() {
		delete[] table;
		delete[] cache;
		delete[] heads;
	}

	void ShiftLeft(int shift) override {
		uint32* row = table;
		uint32* rows_end = table + (row_width << hash_bits);

		while (row < rows_end) {
			*row = (*row < shift || *row == uint32(-1)) ? uint32(-1) : (*row - shift);
			++row;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		if (p_end - p < MinLength()) {
			return;
		}

		uint32 v4 = value32(buf + p);
		uint32 row_base_idx = hash4(v4) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;
		uint32* cache_row = cache + row_base_idx * row_width;

		byte& head = heads[row_base_idx];
		if (++head >= row_width) {
			head -= row_width;
		}

		row[head] = p;
		cache_row[head] = v4;
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		uint32 v4 = value32(buf + p);
		uint32 row_base_idx = hash4(v4) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;
		uint32* cache_row = cache + row_base_idx * row_width;

		byte& head = heads[row_base_idx];
		for (int ridx = 0; ridx < row_width; ridx++) {
			int ri = head - ridx;
			if (ri < 0) {
				ri += row_width;
			}

			uint32 mp = row[ri];
			uint32 v = cache_row[ri];

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (v == v4 && buf[p + best_len] == buf[mp + best_len]) {
				uint len = 4 + try_match(buf, mp + 4, p + 4, p_end);
				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		if (++head >= row_width) {
			head -= row_width;
		}

		row[head] = p;
		cache_row[head] = v4;

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		uint32 v4 = value32(buf + p);
		uint32 row_base_idx = hash4(v4) >> (32 - hash_bits);
		const uint32* row = table + row_base_idx * row_width;
		const uint32* cache_row = cache + row_base_idx * row_width;

		byte& head = heads[row_base_idx];
		for (int ridx = 0; ridx < row_width; ridx++) {
			int ri = head - ridx;
			if (ri < 0) {
				ri += row_width;
			}

			uint32 mp = row[ri];
			uint32 v = cache_row[ri];

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (v == v4 && buf[p + best_len] == buf[mp + best_len]) {
				uint len = 4 + try_match(buf, mp + 4, p + 4, p_end);
				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}
};

struct HT4 : public MatchFinder_Base {
	int hash_bits, row_width;
	uint32* table;

	HT4(int hash_bits, int row_width)
		: hash_bits(hash_bits), row_width(row_width) {
		table = nullptr;
	}

	HT4() : HT4(0, 0) {}

	uint32 MinLength() const override { return 4; }

	void Init() {
		table = new uint32[row_width << hash_bits];

		const int table_mem = sizeof(uint32) * row_width << hash_bits;
		memset(table, -1, table_mem);
		printf("HT4: %d KB\n", (table_mem + 1023) >> 10);
	}

	void Release() {
		delete[] table;
	}

	void ShiftLeft(int shift) override {
		uint32* row = table;
		uint32* rows_end = table + (row_width << hash_bits);

		while (row < rows_end) {
			*row = (*row < shift || *row == uint32(-1)) ? uint32(-1) : (*row - shift);
			++row;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		if (p_end - p < MinLength()) {
			return;
		}

		uint32 v4 = value32(buf + p);
		uint32 row_base_idx = hash4(v4) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;

		uint32 fill_in = p;
		for (int ri = 0; ri < row_width; ri++) {
			uint32 mp = *row;
			*row++ = fill_in;
			fill_in = mp;
		}
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		uint32 v4 = value32(buf + p);
		uint32 row_base_idx = hash4(v4) >> (32 - hash_bits);
		uint32* row = table + row_base_idx * row_width;

		uint32 fill_in = p;
		for (int ri = 0; ri < row_width; ri++) {
			uint32 mp = *row;
			*row++ = fill_in;
			fill_in = mp;

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (buf[p + best_len] == buf[mp + best_len]) {
				uint len = try_match(buf, mp, p, p_end);
				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		int best_len = MinLength() - 1;
		int num_results = 0;

		uint32 h4 = hash4(value32(buf + p));
		uint32 row_base_idx = row_width * (h4 >> (32 - hash_bits));
		const uint32* row = table + row_base_idx * row_width;

		for (int ri = 0; ri < row_width; ri++) {
			uint32 mp = row[ri];

			if (mp >= p || p - mp > max_dist) {
				break;
			}

			if (buf[p + best_len] == buf[mp + best_len]) {
				int len = try_match(buf, mp, p, p_end);
				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = mp;
					if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
						break;
					}
				}
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}
};

struct BT4 : public MatchFinder_Base {
	const static int MAX_TESTS = 128;

	uint32 hash_bits, window_size;
	uint32* head;
	uint32* tree;

	BT4(uint32 window_size, uint32 hash_bits)
		: window_size(window_size), hash_bits(hash_bits) {
		head = nullptr;
		tree = nullptr;
	}

	BT4() : BT4(0, 0) {}

	uint32 MinLength() const override { return 4; }

	void Init() {
		head = new uint32[1 << hash_bits];
		tree = new uint32[2 * window_size];

		const uint64 table_mem = sizeof(uint32) * ((1ULL << hash_bits) + (2ULL * window_size));
		memset(head, -1, sizeof(uint32) << hash_bits);
		memset(tree, -1, sizeof(uint32) * 2 * window_size);
		printf("BT4: %lld KB\n", (table_mem + 1023) >> 10);
	}

	void Release() {
		delete[] head;
		delete[] tree;
	}

	void ShiftLeft(int shift) override {
		uint32* row_head = head;
		const uint32* head_end = head + (1 << hash_bits);

		while (row_head < head_end) {
			*row_head = (*row_head < shift || *row_head == uint32(-1)) ? uint32(-1) : (*row_head - shift);
			++row_head;
		}

		uint32 child_offset = shift << 1;
		uint32* tree_head = tree + child_offset;
		const uint32* tree_end = tree + (2U * window_size);
		while (tree_head < tree_end) {
			uint32* child = tree_head - child_offset;
			*child = (*tree_head < shift || *tree_head == uint32(-1)) ? uint32(-1) : (*tree_head - shift);
			++tree_head;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		if (p_end - p < MinLength()) {
			return;
		}

		ASSERT_DEBUG((p << 1) < (2 * window_size));
		uint32* right = tree + (p << 1);
		uint32* left = right + 1;

		uint left_len = 0, right_len = 0;

		uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
		uint32 sp = head[h4];
		head[h4] = p;

		int tests_left = MAX_TESTS;
		while (sp != uint32(-1) && p > sp
			&& p - sp <= max_dist
			&& tests_left--) {
			uint32* pair = tree + (sp << 1);

			uint len = Min(left_len, right_len);
			if (buf[sp + len] == buf[p + len]) {
				++len;
				len += try_match(buf, sp + len, p + len, p_end);

				if (p_end - p == len) {
					break;
				}
			}

			if (buf[sp + len] < buf[p + len]) {
				*right = sp;
				right = pair + 1;
				sp = *right;
				right_len = len;
			}
			else {
				*left = sp;
				left = pair;
				sp = *left;
				left_len = len;
			}
		}

		*left = -1;
		*right = -1;
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		uint best_len = MinLength() - 1;
		int num_results = 0;

		ASSERT_DEBUG((p << 1) < (2 * window_size));
		uint32* right = tree + (p << 1);
		uint32* left = right + 1;

		uint left_len = 0, right_len = 0;

		uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
		uint32 sp = head[h4];
		head[h4] = p;

		int tests_left = MAX_TESTS;
		while (sp != uint32(-1) && p > sp
			&& p - sp <= max_dist
			&& tests_left--) {
			uint32* pair = tree + (sp << 1);

			uint len = Min(left_len, right_len);
			if (buf[sp + len] == buf[p + len]) {
				++len;
				len += try_match(buf, sp + len, p + len, p_end);

				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = sp;
				}

				if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
					break;
				}
			}

			if (buf[sp + len] < buf[p + len]) {
				*right = sp;
				right = pair + 1;
				sp = *right;
				right_len = len;
			}
			else {
				*left = sp;
				left = pair;
				sp = *left;
				left_len = len;
			}
		}

		*left = -1;
		*right = -1;

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		if (p_end - p < MinLength()) {
			return match_results;
		}

		int best_len = MinLength() - 1;
		int num_results = 0;

		ASSERT_DEBUG((p << 1) < (2 * window_size));
		const uint32* right = tree + (p << 1);
		const uint32* left = right + 1;

		uint left_len = 0, right_len = 0;

		uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
		uint32 sp = head[h4];

		int tests_left = MAX_TESTS;
		while (sp != uint32(-1) && p > sp
			&& p - sp <= max_dist
			&& tests_left--) {
			const uint32* pair = tree + (sp << 1);

			uint len = Min(left_len, right_len);
			if (buf[sp + len] == buf[p + len]) {
				++len;
				len += try_match(buf, sp + len, p + len, p_end);

				if (len > best_len) {
					best_len = len;
					match_results[num_results++] = len;
					match_results[num_results++] = sp;
				}

				if (p_end - p == len || num_results == MATCH_MAX_RESULT_ROWS) {
					break;
				}
			}

			if (buf[sp + len] < buf[p + len]) {
				right = pair + 1;
				sp = *right;
				right_len = len;
			}
			else {
				left = pair;
				sp = *left;
				left_len = len;
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}
};

struct RK256 : public MatchFinder_Base {
	const static int RK_MIN_LENGTH = 64;
	const static int BLOCK_BITS = 8;
	const static int BLOCK_SIZE = 1 << BLOCK_BITS;
	const static int BLOCK_MASK = BLOCK_SIZE - 1;

	const static uint32 ADDH = 0x2F0FD693u;

	//uint32 remh = 1; for (int i = 0; i < BLOCK_SIZE; i++) { remh *= addh; }
	const static uint32 REMH = 0x0E4EA401u;

	INLINE uint32 rolling_hash_add(uint32 p, int y) { return y + p * ADDH; }
	INLINE uint32 rolling_hash_add_remove(uint32 p, int yr, int yl) { return yr + p * ADDH - yl * REMH; }

	uint32 hash_bits, window_size;
	uint16* cache;
	uint32* table;

	uint32 rh, rh_end;

	RK256(uint32 window_size, uint32 hash_bits)
		: window_size(window_size), hash_bits(hash_bits) {
		cache = nullptr;
		table = nullptr;
	}

	RK256() : RK256(0, 0) {}

	uint32 MinLength() const override { return RK_MIN_LENGTH; }

	void Init() {
		cache = new uint16[1 << hash_bits];
		table = new uint32[1 << hash_bits];

		memset(cache, -1, sizeof(uint16) << hash_bits);
		memset(table, -1, sizeof(uint32) << hash_bits);
		printf("RK256: %lld KB\n", (((sizeof(uint32) + sizeof(uint16)) << hash_bits) + 1023) >> 10);

		rh = 0;
		rh_end = 0;
	}

	void Release() {
		delete[] cache;
		delete[] table;
	}

	void ShiftLeft(int shift) override {
		uint32* row_head = table;
		const uint32* head_end = table + (1 << hash_bits);

		while (row_head < head_end) {
			*row_head = (*row_head < shift || *row_head == uint32(-1)) ? uint32(-1) : (*row_head - shift);
			++row_head;
		}

		if (rh_end >= shift + BLOCK_SIZE) {
			rh_end -= shift;
		}
		else {
			rh = 0;
			rh_end = 0;
		}
	}

	void Update(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		while (rh_end < BLOCK_SIZE && rh_end < p_end) {
			rh = rolling_hash_add(rh, buf[rh_end]);
			++rh_end;

			if (!(rh_end & BLOCK_MASK)) {
				cache[rh >> (32 - hash_bits)] = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
				table[rh >> (32 - hash_bits)] = rh_end;
			}
		}

		while (p >= rh_end && rh_end < p_end) {
			rh = rolling_hash_add_remove(rh, buf[rh_end], buf[rh_end - BLOCK_SIZE]);
			++rh_end;

			if (!(rh_end & BLOCK_MASK)) {
				cache[rh >> (32 - hash_bits)] = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
				table[rh >> (32 - hash_bits)] = rh_end;
			}
		}

		while (rh_end >= BLOCK_SIZE && rh_end - p < BLOCK_SIZE && rh_end < p_end) {
			rh = rolling_hash_add_remove(rh, buf[rh_end], buf[rh_end - BLOCK_SIZE]);
			++rh_end;

			if (!(rh_end & BLOCK_MASK)) {
				cache[rh >> (32 - hash_bits)] = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
				table[rh >> (32 - hash_bits)] = rh_end;
			}
		}
	}

	uint32* FindAndUpdate(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;

		uint best_len = MinLength() - 1;
		int num_results = 0;

		while (rh_end < BLOCK_SIZE && rh_end < p_end) {
			rh = rolling_hash_add(rh, buf[rh_end]);
			++rh_end;

			if (!(rh_end & BLOCK_MASK)) {
				cache[rh >> (32 - hash_bits)] = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
				table[rh >> (32 - hash_bits)] = rh_end;
			}
		}

		while (p >= rh_end && rh_end < p_end) {
			rh = rolling_hash_add_remove(rh, buf[rh_end], buf[rh_end - BLOCK_SIZE]);
			++rh_end;

			uint16& cache_end = cache[rh >> (32 - hash_bits)];
			uint32& hist_end = table[rh >> (32 - hash_bits)];
			uint16 hash_cur = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
			if (hist_end < rh_end && hist_end >= BLOCK_SIZE && cache_end == hash_cur) {
				uint32 sp = hist_end - BLOCK_SIZE;
				uint32 mp = rh_end - BLOCK_SIZE;

				ASSERT(p >= mp);
				ASSERT(p > sp);

				uint32 pos_delta = p - mp;
				sp += pos_delta;
				mp += pos_delta;

				if (p > sp
					&& p - sp < max_dist && !num_results) {
					int len = try_match(buf, sp, mp, p_end);
					if (len > best_len
						&& num_results < MATCH_MAX_RESULT_ROWS) {
						best_len = len;
						match_results[num_results++] = len;
						match_results[num_results++] = sp;
					}
				}
			}

			if (!(rh_end & BLOCK_MASK)) {
				cache_end = hash_cur;
				hist_end = rh_end;
			}
		}

		while (rh_end >= BLOCK_SIZE && rh_end - p < BLOCK_SIZE && rh_end < p_end) {
			rh = rolling_hash_add_remove(rh, buf[rh_end], buf[rh_end - BLOCK_SIZE]);
			++rh_end;

			uint16& cache_end = cache[rh >> (32 - hash_bits)];
			uint32& hist_end = table[rh >> (32 - hash_bits)];
			uint16 hash_cur = hash4(value32(buf + rh_end - BLOCK_SIZE)) >> 16;
			if (hist_end < rh_end && hist_end >= BLOCK_SIZE && cache_end == hash_cur) {
				uint32 sp = hist_end - BLOCK_SIZE;
				uint32 mp = rh_end - BLOCK_SIZE;

				ASSERT(p >= mp);
				ASSERT(p > sp);

				uint32 pos_delta = p - mp;
				sp += pos_delta;
				mp += pos_delta;

				if (p > sp
					&& p - sp < max_dist && !num_results) {
					int len = try_match(buf, sp, mp, p_end);
					if (len > best_len
						&& num_results < MATCH_MAX_RESULT_ROWS) {
						best_len = len;
						match_results[num_results++] = len;
						match_results[num_results++] = sp;
					}
				}
			}

			if (!(rh_end & BLOCK_MASK)) {
				cache_end = hash_cur;
				hist_end = rh_end;
			}
		}

		match_results[num_results] = -1;
		return match_results;
	}

	uint32* Find(const byte* buf, uint32 p, uint32 p_end, uint32 max_dist) override {
		match_results[0] = -1;
		return match_results;
	}
};

const static uint32 TABLE_MATCH_SLOTS = 18;
struct MatchNode {
	uint32 match_pos[TABLE_MATCH_SLOTS];
	uint32 max_match_pos;
	uint32 max_len;
};

struct MatcherOpt {
	const static int INCOMPRESSIBLE_STEP = 16;
	const static int INCOMPRESSIBLE_MAX_STEP_BITS = 6;
	const static int NICE_MATCH = 128;
	const static int NICE_MATCH_UNTIL = 16;

	HR2 mf2;
	HR3 mf3;
	BT4 mf4;
	RK256 mf_rk;
	int dry_count, dry_step;

	MatcherOpt(uint32 window_size) {
		mf2 = HR2(14);
		mf3 = HR3(16, 4);
		mf4 = BT4(window_size, 18);
		mf_rk = RK256(window_size, 19);
		dry_count = dry_step = 0;
	}

	void Init() {
		mf2.Init();
		mf3.Init();
		mf4.Init();
		mf_rk.Init();
		dry_count = dry_step = 0;
	}

	void Release() {
		mf2.Release();
		mf3.Release();
		mf4.Release();
		mf_rk.Release();
	}

	void ShiftLeft(int shift) {
		mf2.ShiftLeft(shift);
		mf3.ShiftLeft(shift);
		mf4.ShiftLeft(shift);
		mf_rk.ShiftLeft(shift);
	}

	void UpdateBlock(MatchNode* opt_table, const byte* window, int start, int end, int peek_end, int history_size) {
		for (int p = start; p < peek_end; ++p) {
			MatchNode* node = opt_table + p - start;
			node->max_len = MATCH_MIN_BASE - 1;
			node->max_match_pos = -1;
			memset(node->match_pos, -1, sizeof(uint32) * TABLE_MATCH_SLOTS);
		}

		int run_step = 1;
		int p = start;
		while (p < end) {
			MatchNode* node = opt_table + p - start;

			const uint32* result_set[4] = {
				mf2.FindAndUpdate(window, p, peek_end, history_size),
				mf3.FindAndUpdate(window, p, peek_end, history_size),
				mf4.FindAndUpdate(window, p, peek_end, history_size),
				mf_rk.FindAndUpdate(window, p, peek_end, history_size)
			};

			for (int i = 0; i < 4; i++) {
				const uint32* mresults = result_set[i];
				while (*mresults != uint32(-1)) {
					uint32 len = *mresults++;
					uint32 mp = *mresults++;

					ASSERT_MATCH(p - mp <= history_size);
					ASSERT_MATCH(!memcmp(window + p, window + mp, len));
					ASSERT_MATCH(p + len <= peek_end);

					if (len >= MATCH_LIMIT) { len = MATCH_LIMIT - 1; }
					if (len <= node->max_len) { continue; }

					for (int c = MATCH_MIN_BASE; c <= len && c < TABLE_MATCH_SLOTS + MATCH_MIN_BASE; c++) {
						uint32& hist_p = node->match_pos[c - MATCH_MIN_BASE];
						if (hist_p == uint32(-1) || hist_p < mp) {
							node->match_pos[c - MATCH_MIN_BASE] = mp;
						}
					}

					node->max_len = len;
					node->max_match_pos = mp;
				}
			}

			if (node->max_len < MATCH_MIN_BASE) {
				if (++dry_count >= INCOMPRESSIBLE_STEP) {
					dry_step += dry_step < INCOMPRESSIBLE_MAX_STEP_BITS;
					dry_count = 0;
					run_step = 1 << dry_step;
				}
			}
			else {
				dry_count = 0;
				dry_step >>= 1;
				run_step = 1 << dry_step;
			}

			if (node->max_len >= NICE_MATCH) {
				uint32 len = node->max_len;
				uint32 mp = node->max_match_pos;

				++p, ++mp;
				--len;
				bool can_expand = true;
				while (p < end && len >= NICE_MATCH_UNTIL) {
					uint32* _mresults = mf_rk.FindAndUpdate(window, p, peek_end, history_size);
					while (*_mresults != uint32(-1)) {
						uint32 _len = *_mresults++;
						uint32 _mp = *_mresults++;
						ASSERT_MATCH(p - _mp <= history_size);
						ASSERT_MATCH(!memcmp(window + p, window + _mp, _len));
						ASSERT_MATCH(p + _len <= peek_end);

						if (_len >= MATCH_LIMIT) { _len = MATCH_LIMIT - 1; }

						if (_len > len) {
							len = _len;
							mp = _mp;
						}
					}

					if (can_expand && p + len < peek_end && (can_expand = window[mp + len] == window[p + len])) {
						++len;
					}

					ASSERT_MATCH(!memcmp(window + p, window + mp, len));

					MatchNode* n2 = opt_table + p - start;
					n2->max_len = len;
					n2->max_match_pos = mp;
					for (int c = MATCH_MIN_BASE; c <= len && c < TABLE_MATCH_SLOTS + MATCH_MIN_BASE; c++) {
						n2->match_pos[c - MATCH_MIN_BASE] = mp;
					}

					++p, ++mp;
					--len;
				}
			}
			else {
				p += run_step;
			}
		}

		while (p < peek_end) {
			MatchNode* node = opt_table + p - start;
			const uint32* result_set[4] = {
				mf2.Find(window, p, peek_end, history_size),
				mf3.Find(window, p, peek_end, history_size),
				mf4.Find(window, p, peek_end, history_size),
				mf_rk.Find(window, p, peek_end, history_size)
			};

			for (int i = 0; i < 4; i++) {
				const uint32* mresults = result_set[i];
				while (*mresults != uint32(-1)) {
					uint32 len = *mresults++;
					uint32 mp = *mresults++;

					ASSERT_MATCH(p - mp <= history_size);
					ASSERT_MATCH(!memcmp(window + p, window + mp, len));
					ASSERT_MATCH(p + len <= peek_end);

					if (len >= MATCH_LIMIT) { len = MATCH_LIMIT - 1; }
					if (len <= node->max_len) { continue; }

					for (int c = MATCH_MIN_BASE; c <= len && c < TABLE_MATCH_SLOTS + MATCH_MIN_BASE; c++) {
						uint32& hist_p = node->match_pos[c - MATCH_MIN_BASE];
						if (hist_p == uint32(-1) || hist_p < mp) {
							node->match_pos[c - MATCH_MIN_BASE] = mp;
						}
					}

					node->max_len = len;
					node->max_match_pos = mp;
				}
			}

			if (node->max_len < MATCH_MIN_BASE) {
				if (++dry_count >= INCOMPRESSIBLE_STEP) {
					dry_step += dry_step < INCOMPRESSIBLE_MAX_STEP_BITS;
					dry_count = 0;
					run_step = 1 << dry_step;
				}
			}
			else {
				dry_count = 0;
				dry_step >>= 1;
				run_step = 1 << dry_step;
			}

			if (node->max_len >= NICE_MATCH) {
				uint32 len = node->max_len;
				uint32 mp = node->max_match_pos;

				++p;
				--len;
				++mp;
				bool can_expand = true;
				while (p < peek_end && len >= NICE_MATCH_UNTIL) {
					MatchNode* n2 = opt_table + p - start;

					if (can_expand && p + len < peek_end && (can_expand = window[mp + len] == window[p + len])) {
						++len;
					}

					ASSERT_MATCH(!memcmp(window + p, window + mp, len));

					n2->max_len = len;
					n2->max_match_pos = mp;
					for (int c = MATCH_MIN_BASE; c <= len && c < TABLE_MATCH_SLOTS + MATCH_MIN_BASE; c++) {
						n2->match_pos[c - MATCH_MIN_BASE] = mp;
					}

					++p;
					--len;
					++mp;
				}
			}
			else {
				p += run_step;
			}
		}
	}
};

struct MatcherFast {
	HR2 mf2;
	HR3 mf3;
	HT4 mf4;
	RK256 mf_rk;

	MatcherFast(uint32 window_size) {
		mf2 = HR2(14);
		mf3 = HR3(16, 8);
		mf4 = HT4(18, 24);
		mf_rk = RK256(window_size, 19);
		//mf4 = HC4(window_size, 16, 24, 4);
		//mf4 = BT4(window_size, 16, 8);
	}

	void Init() {
		mf2.Init();
		mf3.Init();
		mf4.Init();
		mf_rk.Init();
	}

	void Release() {
		mf2.Release();
		mf3.Release();
		mf4.Release();
		mf_rk.Release();
	}

	void ShiftLeft(int shift) {
		mf2.ShiftLeft(shift);
		mf3.ShiftLeft(shift);
		mf4.ShiftLeft(shift);
		mf_rk.ShiftLeft(shift);
	}

	uint32 PeekLazyLength(const byte* window, int p, int end, int history_size) {
		const uint32* result_set[2] = {
			mf2.Find(window, p, end, history_size),
			mf3.Find(window, p, end, history_size)
		};

		uint32 max_len = 0;
		for (int i = 0; i < 2; i++) {
			const uint32* mresults = result_set[i];
			while (*mresults != uint32(-1)) {
				uint32 len = *mresults++;
				uint32 mp = *mresults++;

				if (len > max_len) {
					max_len = len;
				}
			}
		}
		return max_len;
	}

	void Update(const byte* window, int p, int end, int history_size, uint32& best_pos, uint32& best_len) {
		mf2.Update(window, p, end, history_size);
		mf3.Update(window, p, end, history_size);
		mf4.Update(window, p, end, history_size);
		mf_rk.Update(window, p, end, history_size);
	}

	void FindAndUpdate(const byte* window, int p, int end, int history_size, uint32& best_pos, uint32& best_len) {
		const uint32* result_set[4] = {
			mf2.FindAndUpdate(window, p, end, history_size),
			mf3.FindAndUpdate(window, p, end, history_size),
			mf4.FindAndUpdate(window, p, end, history_size),
			mf_rk.FindAndUpdate(window, p, end, history_size)
		};

		best_len = 0;
		uint32 cur_lv = 0;

		for (int i = 0; i < 4; i++) {
			const uint32* mresults = result_set[i];
			while (*mresults != uint32(-1)) {
				uint32 len = *mresults++;
				uint32 mp = *mresults++;

				ASSERT_MATCH(p - mp <= history_size);
				ASSERT_MATCH(!memcmp(window + p, window + mp, len));
				ASSERT_MATCH(p + len <= end);

				if (len >= MATCH_LIMIT) { len = MATCH_LIMIT - 1; }

				uint32 lv = get_match_min(p - mp + 1);
				if (len < lv) {
					continue;
				}

				lv = len - lv;
				if (len > best_len) {
					best_len = len;
					best_pos = mp;
				}
			}
		}

		if (best_len >= 2) {
			if (PeekLazyLength(window, p + 1, end, history_size) > best_len) {
				best_len = 0;
			}
			else if (PeekLazyLength(window, p + 2, end, history_size) > best_len + 1) {
				best_len = 0;
			}
		}
	}
};

struct ParsePath {
	byte cmd;
	uint32 local_cost, forward_cost;

	union {
		byte rep_index;
		//uint16 rolz_index;
	};

	uint32 dict_match_pos;

	uint32 match_len;
	//uint32 rolz_len_diff;
};

void opt_parse_path_matches(const byte* window, ParsePath* path, const MatchNode* match_table, State& state, int opt_start_at, int block_end) {
	for (int op = block_end - 1; op >= opt_start_at; --op) {
		int pb = state.pos_align_bits ? (op & state.pos_align_bits) : 0;
		int ctx = op > 0 ? window[op - 1] : 0;
		ParsePath* cn = path + op - opt_start_at;

		uint32 cost_next = path[op - opt_start_at + 1].forward_cost;
		cn->forward_cost = cost_next + cn->local_cost;

		const MatchNode* opt = match_table + op - opt_start_at;
		if (opt->max_len < MATCH_MIN_BASE) {
			continue;
		}

		uint32 max_match_pos = opt->max_match_pos;
		int max_len = opt->max_len;

		if (max_len >= MATCH_LIMIT) {
			max_len = MATCH_LIMIT - 1;
		}

		if (op + max_len > block_end) {
			max_len = block_end - op;
		}

		int test_step = (max_len - MATCH_MIN_BASE) >> 3;
		test_step += test_step == 0;

		int test_len = max_len;
		while (test_len >= MATCH_MIN_BASE) {
			int node_slot = test_len - MATCH_MIN_BASE;

			uint32 match_pos = node_slot < TABLE_MATCH_SLOTS ? opt->match_pos[node_slot] : max_match_pos;
			uint32 dv = op - match_pos - 1;

			ASSERT_MATCH(!memcmp(window + op, window + match_pos, test_len));
			ASSERT_MATCH(match_pos + test_len <= block_end);

			if (test_len < get_match_min(dv)) {
				--test_len;
				continue;
			}

			uint32 lv = test_len - get_match_min(dv);
			ASSERT(lv < 0x7FFFFF);

			uint32 local_cost = state.PriceMatchCached(lv, dv, pb, ctx);
			uint32 forward_cost = local_cost + path[op - opt_start_at + test_len].forward_cost;
			if (forward_cost < cn->forward_cost) {
				cn->cmd = CMD_DICT_MATCH;
				cn->local_cost = local_cost;
				cn->forward_cost = forward_cost;
				cn->match_len = test_len;
				cn->dict_match_pos = match_pos;
			}

			if (test_len >= MatcherOpt::NICE_MATCH) {
				break;
			}

			test_len -= test_step;
		}
	}
}

void opt_parse_path_matches_into_rep4(const byte* window, ParsePath* path, const MatchNode* match_table, State& state, int opt_current_p, int opt_start_at, int opt_next_at, int block_end) {
	RepState forward_rep = state.rep;

	int op = opt_current_p;
	while (op < opt_next_at) {
		int pb = state.pos_align_bits ? (op & state.pos_align_bits) : 0;
		int ctx = op > 0 ? window[op - 1] : 0;

		ParsePath* cn = path + op - opt_start_at;
		ASSERT(cn->cmd == CMD_LITERAL || cn->cmd == CMD_DICT_MATCH);
		for (int i = 0; i < 4; i++) {
			uint32 rp = op - forward_rep.entry[i];
			if (rp >= op) {
				continue;
			}

			int rlen = try_match(window, rp, op, block_end);
			if (op - opt_start_at + rlen > block_end) {
				rlen = block_end - op + opt_start_at;
			}

			if (rlen >= MATCH_LIMIT) {
				rlen = MATCH_LIMIT - 1;
			}

			if (rlen < MATCH_MIN_BASE) {
				continue;
			}

			ASSERT_MATCH(op - opt_start_at + rlen <= block_end);

			for (int test_len = rlen; test_len >= MATCH_MIN_BASE; test_len--) {
				uint32 lv = test_len - MATCH_MIN_BASE;
				uint32 local_cost = state.PriceRepCached(lv, i, pb, ctx);
				uint32 forward_cost = local_cost + path[op - opt_start_at + test_len].forward_cost;

				if (forward_cost < cn->forward_cost) {
					cn->cmd = CMD_REP4;
					cn->local_cost = local_cost;
					cn->forward_cost = forward_cost;
					cn->dict_match_pos = rp;
					cn->match_len = test_len;
					cn->rep_index = i;
					break;
				}
			}
		}

		if (cn->cmd == CMD_REP4) {
			uint32 dist = forward_rep.entry[cn->rep_index];
			ASSERT(op - cn->dict_match_pos == dist);
			ASSERT(dist < op);

			forward_rep.Add(dist);
			op += cn->match_len;
		}
		else if (cn->cmd == CMD_DICT_MATCH) {
			forward_rep.Add(op - cn->dict_match_pos);
			op += cn->match_len;
		}
		else {
			++op;
		}
	}
}

void opt_parse_path_fill_literals(const byte* window, ParsePath* path, const MatchNode* match_table, State& state, int opt_start_at, int block_end) {
	for (int op = block_end - 1; op >= opt_start_at; --op) {
		int pb = state.pos_align_bits ? (op & state.pos_align_bits) : 0;
		int ctx = op > 0 ? window[op - 1] : 0;

		ParsePath* cn = path + op - opt_start_at;

		uint32 local_cost = state.PriceLiteralCached(window[op], pb, ctx);
		uint32 cost_next = path[op - opt_start_at + 1].forward_cost;

		cn->cmd = CMD_LITERAL;
		cn->local_cost = local_cost;
		cn->forward_cost = cost_next + cn->local_cost;
	}
}

void opt_parse_path(const byte* window, ParsePath* path, const MatchNode* match_table, State& state, int opt_current_p, int opt_start_at, int opt_next_at, int block_end) {
	path[block_end - opt_start_at].cmd = CMD_LITERAL;
	path[block_end - opt_start_at].forward_cost = 0;
	path[block_end - opt_start_at].local_cost = 0;

	opt_parse_path_fill_literals(window, path, match_table, state, opt_start_at, block_end);
	opt_parse_path_matches(window, path, match_table, state, opt_start_at, block_end);
	opt_parse_path_matches_into_rep4(window, path, match_table, state, opt_current_p, opt_start_at, opt_next_at, block_end);
}

// Forward optimal parsing, better at handling state but requires more memory (prev links + visited nodes)
//void opt_parse_path_alt(const byte* window, ParsePath* _path, const MatchNode* match_table, State& state, int opt_current_p, int opt_start_at, int opt_next_at, int block_end) {
//	const int MAX_STATE_SIZE = 0x100;
//	const int MAX_STATE_MASK = MAX_STATE_SIZE - 1;
//
//	ParsePath table[OPT_BLOCK_SIZE + OPT_BLOCK_PEEK + 1];
//	uint32 prev[OPT_BLOCK_SIZE + OPT_BLOCK_PEEK + 1];
//	RepState rep_state[MAX_STATE_SIZE];
//
//	const int stack_mem = (sizeof(table) + sizeof(prev) + sizeof(rep_state) + 1023) >> 10;
//
//	memset(table, 0, sizeof(table));
//	memset(prev, 0, sizeof(prev));
//
//	int path_len = block_end - opt_start_at;
//	int path_offs = opt_current_p - opt_start_at;
//
//	ASSERT(path_len <= OPT_BLOCK_SIZE + OPT_BLOCK_PEEK + 1);
//
//	table[path_len].cmd = CMD_LITERAL;
//	table[path_len].carried_cost = 0;
//	table[path_len].local_cost = 0;
//
//	for (int op = 0; op < path_offs; op++) {
//		table[op].carried_cost = -1;
//		prev[op] = 0;
//	}
//
//	rep_state[path_offs & MAX_STATE_MASK] = state.rep;
//	table[path_offs].carried_cost = 0;
//	prev[path_offs] = -1;
//
//	for (int op = path_offs + 1; op <= path_len; op++) {
//		table[op].carried_cost = -1;
//	}
//
//	for (int op = path_offs; op < path_len; op++) {
//		int pb = state.pos_align_bits ? (op & state.pos_align_bits) : 0;
//		int ctx = op > 0 ? window[op - 1] : 0;
//
//		uint32 literal_cost = state.PriceLiteralCached(window[op], pb, ctx);
//		if (literal_cost + table[op].carried_cost < table[1 + op].carried_cost) {
//			prev[1 + op] = op;
//			table[1 + op].cmd = CMD_LITERAL;
//			table[1 + op].carried_cost = literal_cost + table[op].carried_cost;
//
//			rep_state[(1 + op) & MAX_STATE_MASK] = rep_state[op & MAX_STATE_MASK];
//		}
//
//		for (int rep_index = 0; rep_index < 4; rep_index++) {
//			uint32 rp = op + opt_start_at - rep_state[op & MAX_STATE_MASK].entry[rep_index];
//			if (rp >= op + opt_start_at) {
//				continue;
//			}
//
//			int max_len = try_match(window, rp, op + opt_start_at, block_end);
//			if (max_len >= MAX_STATE_SIZE) { max_len = MAX_STATE_MASK; }
//			if (op + max_len > path_len) { max_len = path_len - op; }
//
//			if (max_len < MATCH_MIN_BASE) {
//				continue;
//			}
//
//			int test_len = max_len;
//			while (test_len >= MATCH_MIN_BASE) {
//				ASSERT_MATCH(!memcmp(window + op + opt_start_at, window + rp, test_len));
//				ASSERT_MATCH(rp + test_len <= opt_start_at + path_len);
//				uint32 lv = test_len - MATCH_MIN_BASE;
//
//				uint32 rep_cost = state.PriceRepCached(lv, rep_index, pb, ctx);
//				if (rep_cost + table[op].carried_cost < table[op + test_len].carried_cost) {
//					prev[test_len + op] = op;
//					table[op + test_len].cmd = CMD_REP4;
//					table[op + test_len].rep_index = rep_index;
//					table[op + test_len].dict_match_pos = rp;
//					table[op + test_len].carried_cost = rep_cost + table[op].carried_cost;
//
//					rep_state[(test_len + op) & MAX_STATE_MASK] = rep_state[op & MAX_STATE_MASK];
//					rep_state[(test_len + op) & MAX_STATE_MASK].Add(op + opt_start_at - rp);
//				}
//				--test_len;
//			}
//		}
//
//		const MatchNode* opt = match_table + op;
//		if (opt->max_len >= MATCH_MIN_BASE) {
//			uint32 max_match_pos = opt->max_match_pos;
//			int max_len = opt->max_len;
//
//			if (max_len >= MAX_STATE_SIZE) { max_len = MAX_STATE_MASK; }
//			if (op + max_len > path_len) { max_len = path_len - op; }
//
//			int test_len = max_len;
//			while (test_len >= MATCH_MIN_BASE) {
//				int node_slot = test_len - MATCH_MIN_BASE;
//
//				uint32 match_pos = node_slot < TABLE_MATCH_SLOTS ? opt->match_pos[node_slot] : max_match_pos;
//				uint32 dv = op + opt_start_at - match_pos - 1;
//
//				ASSERT_MATCH(!memcmp(window + op + opt_start_at, window + match_pos, test_len));
//				ASSERT_MATCH(match_pos + test_len <= opt_start_at + path_len);
//
//				if (test_len < get_match_min(dv)) {
//					--test_len;
//					continue;
//				}
//
//				uint32 lv = test_len - get_match_min(dv);
//				ASSERT(lv < 0x7FFFFF);
//
//				uint32 dict_cost = state.PriceMatchCached(lv, dv, pb, ctx);
//				if (dict_cost + table[op].carried_cost < table[op + test_len].carried_cost) {
//					prev[test_len + op] = op;
//					table[op + test_len].cmd = CMD_DICT_MATCH;
//					table[op + test_len].dict_match_pos = match_pos;
//					table[op + test_len].carried_cost = dict_cost + table[op].carried_cost;
//
//					rep_state[(test_len + op) & MAX_STATE_MASK] = rep_state[op & MAX_STATE_MASK];
//					rep_state[(test_len + op) & MAX_STATE_MASK].Add(op + opt_start_at - match_pos);
//				}
//				--test_len;
//			}
//		}
//	}
//
//	int op = path_len;
//	while (op > path_offs) {
//		ParsePath* cn = table + op;
//		ParsePath* cn_from = table + prev[op];
//		ASSERT(prev[op] >= 0);
//		ASSERT(prev[op] < op);
//		ASSERT(cn->carried_cost > cn_from->carried_cost);
//
//		_path[prev[op]] = *cn;
//		_path[prev[op]].match_len = op - prev[op];
//		_path[prev[op]].local_cost = cn->carried_cost - cn_from->carried_cost;
//
//		op = prev[op];
//	}
//}

void encode_dist_value(State& state, BitCoder& code, int pb, uint32 dv, uint32 lv) {
	ASSERT(dv < 0x7FFFFFFF);

	int dist_idx = lv < 4 ? lv : 3;
	if (dv < 4) {
		state.dist_slots[dist_idx].Encode(dv, code);
	}
	else {
		int nb = clz32(dv) + 1;

		int add_bits = nb - 2;
		int top = dv >> add_bits;
		int add_mask = (1 << add_bits) - 1;
		int slot = ((nb - 1) << 1) + (top & 1);

		state.dist_slots[dist_idx].Encode(slot, code);
		if (add_bits < 4) {
			code.EncodeBits(dv & add_mask, add_bits);
		}
		else {
			if (add_bits > 4) {
				code.EncodeBits((dv & add_mask) >> 4, add_bits - 4);
			}
			state.dist_low.EncodeReverse(dv & 15, code);
		}
	}
}

void encode_length_value(State& state, BitCoder& code, int pb, uint32 lv) {
	ASSERT(lv < 0x7FFFFF);

	if (lv < 8) {
		state.len_marker.Encode(0, code);
		state.len_low[pb].Encode(lv, code);
	}
	else if (lv < 16) {
		state.len_marker.Encode(1, code);
		state.len_mid[pb].Encode(lv - 8, code);
	}
	else if (lv < 272) {
		state.len_marker.Encode(2, code);
		state.len_high.Encode(lv - 16, code);
	}
	else {
		ASSERT(lv < 4368);
		state.len_marker.Encode(3, code);
		code.EncodeBits(lv - 272, 12);
	}
}

void compress_opt(FILE* fin, FILE* fout, byte history_bits, int pos_align_bits, int context_bits, uint32 crc32) {
	fseek(fin, 0, SEEK_END);
	int64 total_file_size = _ftelli64(fin);
	fseek(fin, 0, SEEK_SET);

	if (history_bits < 16) {
		history_bits = 16;
	}

	while (history_bits > 16 && total_file_size < (1 << (history_bits - 1))) {
		--history_bits;
	}

	const uint32 history_size = 1 << history_bits;
	const uint32 buffer_io_size = Clamp(history_size / 4, 1U << 16, 1U << 22);
	const uint32 buffer_forward_size = MATCH_LIMIT;
	const uint32 window_size = history_size + buffer_io_size + buffer_forward_size;
	const uint32 slide_size = history_size + buffer_io_size;

	const int opt_block_size = 0x2000;
	const int opt_block_after = opt_block_size + 0x100;

	byte* window = new byte[window_size + match_finder_padding];
	memset(window, 0, window_size + match_finder_padding);

	State state;
	state.Init(pos_align_bits, context_bits);

	MatcherOpt matcher(window_size);
	matcher.Init();

	MatchNode* match_table = new MatchNode[opt_block_after];
	memset(match_table, 0, sizeof(MatchNode) * opt_block_after);

	ParsePath* path = new ParsePath[opt_block_after + 1];
	memset(path, 0, sizeof(ParsePath) * (opt_block_after + 1));

	printf("Window size: %d KB\n", (window_size + match_finder_padding + 1023) >> 10);
	printf("Parser: %d KB\n", uint32((sizeof(ParsePath) + sizeof(MatchNode)) * opt_block_after + 1023) >> 10);
	printf("Bit coder: %llu KB\n", (sizeof(BitCoder) + sizeof(uint16) * BitCoder::RANS_BUFFER_SIZE + 1023) >> 10);

	BitCoder code(fout);
	code.Init();

	code.EncodeBits(history_bits, 8);

	code.EncodeBits(total_file_size & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 16) & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 32) & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 48) & 0xFFFFU, 16);

	code.EncodeBits(crc32 & 0xFFFFU, 16);
	code.EncodeBits((crc32 >> 16) & 0xFFFFU, 16);

	code.EncodeBits(pos_align_bits, 4);
	code.EncodeBits(context_bits, 4);

	int window_length = fread(window, 1, window_size, fin);
	int p = 0;

	char status_buf[70];
	clock_t last_message_time = 0, start_time = clock();
	snprintf(status_buf, sizeof(status_buf), "Working...");
	printf("%-79s\r", status_buf);

	while (true) {
		if (!window_length || p >= slide_size) {
			if (window_length >= buffer_io_size) {
				memmove(window, window + buffer_io_size, window_length - buffer_io_size);

				matcher.ShiftLeft(buffer_io_size);
				state.ShiftLeft(buffer_io_size);

				window_length -= buffer_io_size;
				p -= buffer_io_size;
			}

			int add = fread(window + window_length, 1, window_size - window_length, fin);
			window_length += add;
		}

		if (p >= window_length) {
			break;
		}

		int p_end = slide_size;
		int p_available = window_length;
		if (p_end > p_available) {
			p_end = p_available;
		}

		int opt_start_at = 0, opt_next_at = 0;
		while (p < p_end) {
			if (p >= opt_next_at) {
				clock_t cur_time = clock();
				if (cur_time - last_message_time > CLOCKS_PER_SEC / 4) {
					snprintf(status_buf, sizeof(status_buf), "Working... %lld KB -> %lld KB, %lld KB/sec",
						(_ftelli64(fin) + p - p_end + 1023) >> 10,
						(_ftelli64(fout) + 1023) >> 10,
						int64(CLOCKS_PER_SEC)* ((_ftelli64(fin) + p - p_end + 1023) >> 10) / (1 + cur_time - start_time));
					printf("%-79s\r", status_buf);

					last_message_time = cur_time;
				}

				opt_start_at = p & ~(opt_block_size - 1);
				opt_next_at = opt_start_at + opt_block_size;
				if (opt_next_at > p_end) {
					opt_next_at = p_end;
				}

				int block_end = opt_start_at + opt_block_after;
				if (block_end > p_end) {
					block_end = p_end;
				}

				state.UpdatePriceCache();
				matcher.UpdateBlock(match_table, window, opt_start_at, opt_next_at, block_end, history_size);
				opt_parse_path(window, path, match_table, state, p, opt_start_at, opt_next_at, block_end);

				ASSERT(!path[block_end - opt_start_at].forward_cost);
				ASSERT(!path[block_end - opt_start_at].cmd);
			}

			int pb = state.pos_align_bits ? (p & state.pos_align_bits) : 0;
			int ctx = p > 0 ? window[p - 1] : 0;
			int y = window[p];

			const ParsePath* cur = path + p - opt_start_at;

			int adv = 1;
			if (cur->cmd == CMD_REP4) {
				adv = cur->match_len;

				uint32 match_pos = p - state.rep.entry[cur->rep_index];
				ASSERT_MATCH(match_pos == cur->dict_match_pos);
				ASSERT_MATCH(p + cur->match_len <= window_size);
				ASSERT_MATCH(!memcmp(window + p, window + match_pos, cur->match_len));
				state.rep.Add(state.rep.entry[cur->rep_index]);

				ASSERT(cur->match_len >= MATCH_MIN_BASE);
				adv = cur->match_len;
				uint32 lv = cur->match_len - MATCH_MIN_BASE;

				state.cmd[state.match_hist][pb].Encode(CMD_REP4, code);
				state.NextStateMatch();

				encode_length_value(state, code, pb, lv);
				state.rep_dist.Encode(cur->rep_index, code);
			}
			else if (cur->cmd == CMD_DICT_MATCH) {
				adv = cur->match_len;

				ASSERT_MATCH(p + cur->match_len <= window_size);
				ASSERT_MATCH(!memcmp(window + p, window + cur->dict_match_pos, cur->match_len));

				state.rep.Add(p - cur->dict_match_pos);

				uint32 dv = p - cur->dict_match_pos - 1;
				uint32 min_match_len = get_match_min(dv);
				ASSERT(cur->match_len >= min_match_len);

				uint32 lv = cur->match_len - min_match_len;

				state.cmd[state.match_hist][pb].Encode(CMD_DICT_MATCH, code);
				state.NextStateMatch();
				encode_length_value(state, code, pb, lv);
				encode_dist_value(state, code, pb, dv, lv);
			}
			else {
				auto& lit_coder = state.lit[state.ToLitIndex(pb, ctx)];
				state.cmd[state.match_hist][pb].Encode(CMD_LITERAL, code);
				state.NextStateLiteral();
				lit_coder.Encode(y, code);
			}

			p += adv;
		}
	}

	state.Release();
	matcher.Release();
	code.Flush();

	delete[] path;
	delete[] match_table;
	delete[] window;

	clock_t cur_time = clock();

	snprintf(status_buf, sizeof(status_buf), "Finished... %lld KB -> %lld KB, %lld KB/sec, %.3f sec",
		(_ftelli64(fin) + 1023) >> 10,
		(_ftelli64(fout) + 1023) >> 10,
		int64(CLOCKS_PER_SEC)* ((_ftelli64(fin) + 1023) >> 10) / (1 + cur_time - start_time),
		(1 + cur_time - start_time) / double(CLOCKS_PER_SEC));
	printf("%-79s\r\nInput: %lld\nOutput: %lld\n", status_buf, _ftelli64(fin), _ftelli64(fout));
}

void compress_fast(FILE* fin, FILE* fout, byte history_bits, int pos_align_bits, int context_bits, uint32 crc32) {
	fseek(fin, 0, SEEK_END);
	int64 total_file_size = _ftelli64(fin);
	fseek(fin, 0, SEEK_SET);

	if (history_bits < 16) {
		history_bits = 16;
	}

	while (history_bits > 16 && total_file_size < (1 << (history_bits - 1))) {
		--history_bits;
	}

	const uint32 history_size = 1 << history_bits;
	const uint32 buffer_io_size = Clamp(history_size / 4, 1U << 16, 1U << 22);
	const uint32 buffer_forward_size = MATCH_LIMIT;
	const uint32 window_size = history_size + buffer_io_size + buffer_forward_size;
	const uint32 slide_size = history_size + buffer_io_size;

	const int opt_block_size = 0x2000;
	const int opt_block_after = opt_block_size + 0x100;

	byte* window = new byte[window_size + match_finder_padding];
	memset(window, 0, window_size + match_finder_padding);

	State state;
	state.Init(pos_align_bits, context_bits);

	MatcherFast matcher(window_size);
	matcher.Init();

	printf("Window size: %d KB\n", (window_size + match_finder_padding + 1023) >> 10);
	printf("State: %d KB\n", uint32(sizeof(state) + 1023) >> 10);
	printf("Bit coder: %llu KB\n", (sizeof(BitCoder) + sizeof(uint16) * BitCoder::RANS_BUFFER_SIZE + 1023) >> 10);

	BitCoder code(fout);
	code.Init();

	code.EncodeBits(history_bits, 8);

	code.EncodeBits(total_file_size & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 16) & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 32) & 0xFFFFU, 16);
	code.EncodeBits((total_file_size >> 48) & 0xFFFFU, 16);

	code.EncodeBits(crc32 & 0xFFFFU, 16);
	code.EncodeBits((crc32 >> 16) & 0xFFFFU, 16);

	code.EncodeBits(pos_align_bits, 4);
	code.EncodeBits(context_bits, 4);

	int window_length = fread(window, 1, window_size, fin);
	int p = 0;

	char status_buf[70];
	clock_t last_message_time = 0, start_time = clock();
	snprintf(status_buf, sizeof(status_buf), "Working...");
	printf("%-79s\r", status_buf);

	while (true) {
		if (!window_length || p >= slide_size) {
			if (window_length >= buffer_io_size) {
				memmove(window, window + buffer_io_size, window_length - buffer_io_size);

				matcher.ShiftLeft(buffer_io_size);
				state.ShiftLeft(buffer_io_size);

				window_length -= buffer_io_size;
				p -= buffer_io_size;
			}

			int add = fread(window + window_length, 1, window_size - window_length, fin);
			window_length += add;
		}

		if (p >= window_length) {
			break;
		}

		int p_end = slide_size;
		int p_available = window_length;
		if (p_end > p_available) {
			p_end = p_available;
		}

		int opt_start_at = 0, opt_next_at = 0;
		while (p < p_end) {
			if (p >= opt_next_at) {
				clock_t cur_time = clock();
				if (cur_time - last_message_time > CLOCKS_PER_SEC / 4) {
					snprintf(status_buf, sizeof(status_buf), "Working... %lld KB -> %lld KB, %lld KB/sec",
						(_ftelli64(fin) + p - p_end + 1023) >> 10,
						(_ftelli64(fout) + 1023) >> 10,
						int64(CLOCKS_PER_SEC)* ((_ftelli64(fin) + p - p_end + 1023) >> 10) / (1 + cur_time - start_time));
					printf("%-79s\r", status_buf);

					last_message_time = cur_time;
				}

				opt_start_at = p & ~(opt_block_size - 1);
				opt_next_at = opt_start_at + opt_block_size;
				if (opt_next_at > p_end) {
					opt_next_at = p_end;
				}

				int block_end = opt_start_at + opt_block_after;
				if (block_end > p_end) {
					block_end = p_end;
				}
			}

			int pb = state.pos_align_bits ? (p & state.pos_align_bits) : 0;
			int ctx = p > 0 ? window[p - 1] : 0;
			int y = window[p];

			uint32 match_pos = 0, match_len = 0;
			matcher.FindAndUpdate(window, p, p_available, history_size, match_pos, match_len);

			//uint32 match_pos = match_table[p - opt_start_at].max_match_pos;
			//uint32 match_len = match_table[p - opt_start_at].max_len;
			//if (match_table[p - opt_start_at + 1].max_len > match_len ||
			//	match_table[p - opt_start_at + 2].max_len > match_len) {
			//	match_len = 0;
			//}

			uint32 dv = p - match_pos - 1;
			uint32 min_match_len = get_match_min(dv);

			int rep_idx = -1;
			for (int i = 0; i < 4; i++) {
				uint32 rp = p - state.rep.entry[i];
				if (rp >= p) {
					continue;
				}

				int rlen = try_match(window, rp, p, p_end);
				if (p - opt_start_at + rlen > opt_block_after) { rlen = opt_block_after - p + opt_start_at; }
				if (rlen >= MATCH_LIMIT) { rlen = MATCH_LIMIT - 1; }

				if (rlen >= MATCH_MIN_BASE) {
					ASSERT_MATCH(p - opt_start_at + rlen <= opt_block_after);

					if (rlen > match_len) {
						min_match_len = MATCH_MIN_BASE;
						match_len = rlen;
						rep_idx = i;
					}
				}
			}

			int adv = 1;
			if (rep_idx > -1) {
				match_pos = p - state.rep.entry[rep_idx];
				state.rep.Add(state.rep.entry[rep_idx]);
				ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

				adv = match_len;

				uint32 lv = match_len - min_match_len;

				state.cmd[state.match_hist][pb].Encode(CMD_REP4, code);
				state.NextStateMatch();

				encode_length_value(state, code, pb, lv);
				state.rep_dist.Encode(rep_idx, code);
			}
			else if (match_pos < p && min_match_len <= match_len) {
				ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

				adv = match_len;
				state.rep.Add(p - match_pos);

				uint32 lv = match_len - min_match_len;

				state.cmd[state.match_hist][pb].Encode(CMD_DICT_MATCH, code);
				state.NextStateMatch();

				encode_length_value(state, code, pb, lv);
				encode_dist_value(state, code, pb, dv, lv);
			}
			else {
				auto& lit_coder = state.lit[state.ToLitIndex(pb, ctx)];

				state.cmd[state.match_hist][pb].Encode(CMD_LITERAL, code);
				state.NextStateLiteral();
				lit_coder.Encode(y, code);
			}

			for (int i = 1; i < adv - 1; i++) {
				matcher.Update(window, p + i, p_available, history_size, match_pos, match_len);
			}

			p += adv;
		}
	}

	state.Release();
	matcher.Release();
	code.Flush();

	delete[] window;

	clock_t cur_time = clock();

	snprintf(status_buf, sizeof(status_buf), "Finished... %lld KB -> %lld KB, %lld KB/sec, %.3f sec",
		(_ftelli64(fin) + 1023) >> 10,
		(_ftelli64(fout) + 1023) >> 10,
		int64(CLOCKS_PER_SEC)* ((_ftelli64(fin) + 1023) >> 10) / (1 + cur_time - start_time),
		(1 + cur_time - start_time) / double(CLOCKS_PER_SEC));
	printf("%-79s\r\nInput: %lld\nOutput: %lld\n", status_buf, _ftelli64(fin), _ftelli64(fout));
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

uint32 crc32_calc(const byte* window, int64 n, uint64 crci) {
	uint64 crc = crci ^ 0xFFFFFFFF;
	const byte* next = window;

	while (n && (byte(next) & 15)) {
		crc = crc32_table[0][(crc ^ (*next++)) & 0xFF] ^ (crc >> 8);
		--n;
	}

	while (n >= 16) {
		crc ^= *reinterpret_cast<const uint64*>(next);
		uint64 high = *reinterpret_cast<const uint64*>(next + 8);

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

uint32 crc_file(FILE* f) {
	fseek(f, 0, SEEK_SET);

	byte buf[0x10000];
	int n;
	uint32 crc = 0;
	while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
		crc = crc32_calc(buf, n, crc);
	}

	return crc;
}

void decompress(FILE* fin, FILE* fout, uint32& unpacked_crc32, uint32& running_crc32) {
	BitDecoder code(fin);
	code.Init();

	byte history_bits = code.DecodeBits(8);

	int64 total_file_size = code.DecodeBits(16);
	total_file_size += code.DecodeBits(16) << 16;
	total_file_size += uint64(code.DecodeBits(16)) << 32ULL;
	total_file_size += uint64(code.DecodeBits(16)) << 48ULL;

	unpacked_crc32 = code.DecodeBits(16);
	unpacked_crc32 += code.DecodeBits(16) << 16;

	uint32 pos_align_bits = code.DecodeBits(4);
	uint32 context_bits = code.DecodeBits(4);

	ASSERT(history_bits >= 16);
	ASSERT(history_bits <= 30);

	const uint32 history_size = 1 << history_bits;
	const uint32 buffer_io_size = Clamp(history_size / 8, 1U << 15, 1U << 22);
	const uint32 window_size = history_size + buffer_io_size + MATCH_LIMIT;
	const uint32 slide_size = history_size + buffer_io_size;

	byte* window = new byte[window_size + match_finder_padding];
	memset(window, 0, window_size + match_finder_padding);

	State state;
	state.Init(pos_align_bits, context_bits);

	printf("Window size: %d KB\n", (window_size + match_finder_padding + 1023) >> 10);
	printf("Bit decoder: %llu KB\n", (sizeof(BitDecoder) + 1023) >> 10);

	int wp = 0;
	int64 rem = total_file_size;

	char status_buf[70];
	clock_t last_message_time = 0, start_time = clock();
	snprintf(status_buf, sizeof(status_buf), "Working...");
	printf("%-79s\r", status_buf);

	while (rem > 0) {
		while (wp >= slide_size) {
			clock_t cur_time = clock();
			if (cur_time - last_message_time > CLOCKS_PER_SEC / 4) {
				snprintf(status_buf, sizeof(status_buf), "Working... %lld KB -> %lld KB, %lld KB/sec",
					(_ftelli64(fin) + 1023) >> 10,
					fout ? ((_ftelli64(fout) + wp + 1023) >> 10) : 0,
					fout ? (((_ftelli64(fout) + wp + 1023) >> 10)* int64(CLOCKS_PER_SEC) / (1 + cur_time - start_time)) : 0);
				printf("%-79s\r", status_buf);

				last_message_time = cur_time;
			}

			running_crc32 = crc32_calc(window, buffer_io_size, running_crc32);

			if (fout) {
				fwrite(window, 1, buffer_io_size, fout);
			}

			memmove(window, window + buffer_io_size, wp - buffer_io_size);
			wp -= buffer_io_size;
		}

		int pb = state.pos_align_bits ? (wp & state.pos_align_bits) : 0;
		int ctx = wp > 0 ? window[wp - 1] : 0;

		int cmd = state.cmd[state.match_hist][pb].Decode(code);

		int adv = 1;
		if (cmd == CMD_DICT_MATCH || cmd == CMD_REP4) {
			state.NextStateMatch();

			uint32 lv = 0;
			int len_marker = state.len_marker.Decode(code);
			if (len_marker == 0) {
				lv = state.len_low[pb].Decode(code);
			}
			else if (len_marker == 1) {
				lv = 8 + state.len_mid[pb].Decode(code);
			}
			else if (len_marker == 2) {
				lv = 16 + state.len_high.Decode(code);
			}
			else {
				lv = 272 + code.DecodeBits(12);
			}

			int dist, length;
			if (cmd == CMD_REP4) {
				int rep_idx = state.rep_dist.Decode(code);
				dist = state.rep.entry[rep_idx];
				length = lv + MATCH_MIN_BASE;

				state.rep.Add(dist);
			}
			else {
				int dist_idx = lv < 4 ? lv : 3;
				uint32 dv = state.dist_slots[dist_idx].Decode(code);
				if (dv >= 4) {
					uint32 add_bits = (dv >> 1) - 1;
					dv = (2 + (dv & 1)) << add_bits;

					if (add_bits < 4) {
						dv += code.DecodeBits(add_bits);
					}
					else {
						if (add_bits > 4) {
							dv += code.DecodeBits(add_bits - 4) << 4;
						}
						dv += state.dist_low.DecodeReverse(code);
					}
				}

				dist = dv + 1;
				length = lv + get_match_min(dv);

				state.rep.Add(dist);
			}
			ASSERT(wp - dist >= 0);
			ASSERT(wp + length <= window_size);

			adv = length;
			while (length-- > 0) {
				int y = window[wp - dist];
				//state.rolz.Update(ctx, wp, 1);
				window[wp++] = y;
				//ctx = y;
			}
		}
		else {
			state.NextStateLiteral();
			//state.rolz.Update(ctx, wp, 1);

			int y = state.lit[state.ToLitIndex(pb, ctx)].Decode(code);
			window[wp++] = y;
		}

		rem -= adv;
	}

	if (wp) {
		running_crc32 = crc32_calc(window, wp, running_crc32);

		if (fout) {
			fwrite(window, 1, wp, fout);
		}
	}

	state.Release();
	delete[] window;

	clock_t cur_time = clock();
	snprintf(status_buf, sizeof(status_buf), "Finished... %lld KB -> %lld KB, %lld KB/sec, %.3f sec",
		(_ftelli64(fin) + 1023) >> 10,
		fout ? ((_ftelli64(fout) + 1023) >> 10) : 0,
		fout ? (((_ftelli64(fout) + 1023) >> 10)* int64(CLOCKS_PER_SEC) / (1 + cur_time - start_time)) : 0,
		(1 + cur_time - start_time) / double(CLOCKS_PER_SEC));

	printf("%-79s\nInput: %lld\nOutput: %lld\n", status_buf, _ftelli64(fin), fout ? _ftelli64(fout) : 0);
}

int _tolower(int c) { return c | 0x20; }

void _tolower(char* v) {
	while (*v) {
		*v = _tolower(*v);
		v++;
	}
}

int main(int argc, char** argv) {
	log2_init();
	crc32_init();
	RANS::init();

	FILE* fin = nullptr;
	FILE* fout = nullptr;

	int window_bits = 20;
	int pos_bits = 0;
	int ctx1_bits = 2;

	while (argc >= 2 && *argv[1] == '-') {
		char* arg = argv[1];
		argv++;
		argc--;

		_tolower(arg);
		if (!strncmp(arg, "-window:", 8)) {
			window_bits = Clamp(atoi(arg + 8), 16, 30);
			printf("Window bits: %d\n", window_bits);
		}
		else if (!strncmp(arg, "-posbits:", 9)) {
			pos_bits = Clamp(atoi(arg + 9), 0, 4);
			printf("Position bits: %d\n", pos_bits);
		}
		else if (!strncmp(arg, "-ctx1:", 6)) {
			ctx1_bits = Clamp(atoi(arg + 6), 0, 8);
			printf("ctx1 bits: %d\n", ctx1_bits);
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

		printf("Calculating crc... ");
		uint32 crc32 = crc_file(fin);
		printf("%X\n", crc32);
		fseek(fin, 0, SEEK_SET);

		if (_tolower(argv[1][1]) == 'f') {
			compress_fast(fin, fout, window_bits, pos_bits, ctx1_bits, crc32);
		}
		else {
			compress_opt(fin, fout, window_bits, pos_bits, ctx1_bits, crc32);
		}

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

		uint32 unpacked_crc32 = 0, running_crc32 = 0;
		decompress(fin, fout, unpacked_crc32, running_crc32);
		fclose(fin);
		fclose(fout);

		printf("Expected CRC32 %X, was %X\n", unpacked_crc32, running_crc32);
		return unpacked_crc32 == running_crc32;
	}
	else if (argc == 3 && _tolower(argv[1][0]) == 't') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		uint32 unpacked_crc32 = 0, running_crc32 = 0;
		decompress(fin, nullptr, unpacked_crc32, running_crc32);
		fclose(fin);

		printf("Expected CRC32 %X, was %X\n", unpacked_crc32, running_crc32);
		return unpacked_crc32 == running_crc32;
	}
	else if (argc == 3 && _tolower(argv[1][0]) == 'b') {
		while (true) {
			fin = fopen(argv[2], "rb");
			if (!fin) {
				printf("Error: %s file does not exist\n", argv[2]);
				return -1;
			}

			uint32 unpacked_crc32 = 0, running_crc32 = 0;
			decompress(fin, nullptr, unpacked_crc32, running_crc32);
			fclose(fin);

			printf("Expected CRC32 %X, was %X\n", unpacked_crc32, running_crc32);
		}
	}
	else if (argc == 3 && _tolower(argv[1][0]) == 'h') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		printf("Calculating crc... ");
		uint32 crc32 = crc_file(fin);
		fclose(fin);
		printf("%X\n", crc32);
	}
	else {
		printf("NLZM 1.02 64-bit - Written by Nauful\n"
			"Options:\n"
			"\t[flags] c [input] [output] - Compress input file to output file (best)\n"
			"\t[flags] cf [input] [output] - Compress input file to output file (fast)\n"
			"\td [input] [output] - Decompress input file to output file\n"
			"\tt [input] - Decompress input file in memory and test crc\n"
			"\th [input] - CRC32 for input file\n"
			"Flags:\n"
			"\t-window:bits = Window size in bits, default 20 (1 MB), min 16, max 30 (64 KB to 1 GB)\n"
			"\t-posbits:bits = Position alignment bits, default 0, max 4\n"
			"\t-ctx1:bits = Bits used for o1 prediction, default 2, max 8\n"
			"\t\tposbits + ctx1 limited to 8"
		);
	}

	return 0;
}
