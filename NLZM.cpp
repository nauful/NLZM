// NLZM 64-bit - Written by Nauful
// Released into the public domain.
// Please drop a comment if you find this useful.

#ifndef _CRT_DISABLE_PERFCRIT_LOCKS
#	define _CRT_DISABLE_PERFCRIT_LOCKS
#endif

#ifndef _CRT_SECURE_NO_WARNINGS
#	define _CRT_SECURE_NO_WARNINGS
#endif

#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef HAS_MAIN
#	ifdef _MSC_VER
#		define MSVC
#	endif
#endif

#ifdef MSVC
#	include <intrin.h>
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

namespace Release_NLZM {
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
		unsigned long v0 = 0;
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
		return popcnt32((x & uint32(-int32(x))) - 1);
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

	static uint32 CRC32_TABLE[16][256];
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
			CRC32_TABLE[0][n] = crc;
		}

		for (n = 0; n < 256; n++) {
			crc = CRC32_TABLE[0][n];

			for (k = 1; k < 16; k++) {
				crc = CRC32_TABLE[0][crc & 0xff] ^ (crc >> 8);
				CRC32_TABLE[k][n] = crc;
			}
		}
	}

	uint32 crc32_calc(const byte* window, int64 n, uint64 crci) {
		uint64 crc = crci ^ 0xFFFFFFFF;
		const byte* next = window;

		while (n && (byte(next) & 15)) {
			crc = CRC32_TABLE[0][(crc ^ (*next++)) & 0xFF] ^ (crc >> 8);
			--n;
		}

		while (n >= 16) {
			crc ^= *reinterpret_cast<const uint64*>(next);
			uint64 high = *reinterpret_cast<const uint64*>(next + 8);

			crc =
				CRC32_TABLE[15][crc & 0xFF] ^
				CRC32_TABLE[14][(crc >> 8) & 0xFF] ^
				CRC32_TABLE[13][(crc >> 16) & 0xFF] ^
				CRC32_TABLE[12][(crc >> 24) & 0xFF] ^
				CRC32_TABLE[11][(crc >> 32) & 0xFF] ^
				CRC32_TABLE[10][(crc >> 40) & 0xFF] ^
				CRC32_TABLE[9][(crc >> 48) & 0xFF] ^
				CRC32_TABLE[8][crc >> 56] ^
				CRC32_TABLE[7][high & 0xFF] ^
				CRC32_TABLE[6][(high >> 8) & 0xFF] ^
				CRC32_TABLE[5][(high >> 16) & 0xFF] ^
				CRC32_TABLE[4][(high >> 24) & 0xFF] ^
				CRC32_TABLE[3][(high >> 32) & 0xFF] ^
				CRC32_TABLE[2][(high >> 40) & 0xFF] ^
				CRC32_TABLE[1][(high >> 48) & 0xFF] ^
				CRC32_TABLE[0][high >> 56];

			next += 16;
			n -= 16;
		}

		while (n) {
			crc = CRC32_TABLE[0][(crc ^ (*next++)) & 0xFF] ^ (crc >> 8);
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

	struct BufferedFile {
		const static int BUF_SIZE = 0x1000;
		byte buf[BUF_SIZE];
		int p, n;
		FILE* f;

		BufferedFile(FILE* f) : p(0), n(-1), f(f) {}

		INLINE void Put(int y) {
			buf[p++] = y;
			if (p == BUF_SIZE) {
				fwrite(buf, 1, BUF_SIZE, f);
				p = 0;
			}
		}

		void Flush() {
			if (p) {
				fwrite(buf, 1, p, f);
				p = 0;
			}
		}

		INLINE int Get() {
			if (p >= n) {
				if (n == 0) {
					return -1;
				}

				p = 0;
				n = fread(buf, 1, BUF_SIZE, f);
				if (n == 0) {
					return -1;
				}
			}

			return buf[p++];
		}
	};

	struct BitDecoder {
		const static uint32 TOP = 1 << 24;

		uint32 range, code;
		BufferedFile* f;

		BitDecoder(BufferedFile* f) : f(f) {}

		void Init() {
			code = 0;
			range = 0xFFFFFFFFU;
			for (int i = 0; i < 5; i++) {
				code <<= 8;
				code += f->Get();
			}
		}

		inline void Normalize() {
			while (range < TOP) {
				code <<= 8;
				range <<= 8;

				code += f->Get();
			}
		}

		inline int Decode12P(uint16 p) {
			ASSERT_DEBUG(p < 0x1000);

			uint32 bound = p * (range >> 12);
			int y;
			if (code < bound) {
				range = bound;
				y = 1;
			}
			else {
				code -= bound;
				range -= bound;
				y = 0;
			}

			Normalize();
			return y;
		}

		uint32 DecodeBits(int num_bits) {
			uint32 v = 0;
			for (int i = num_bits; i; i--) {
				range >>= 1;
				uint32 t = (code - range) >> 31;
				code -= range & (t - 1);
				v = (v << 1) | (1 - t);

				if (range < TOP) {
					code <<= 8;
					range <<= 8;

					code += f->Get();
				}
			}

			return v;
		}
	};

	struct BitCoder {
		const static uint32 LOW_THRESHOLD = 0xFF000000U;

		uint32 _ff;
		byte _overflow_cache;

		uint64 low;
		uint32 range;

		BufferedFile* f;

		BitCoder(BufferedFile* f) : f(f) {}

		void Init() {
			low = 0;
			range = 0xFFFFFFFF;
			_ff = 1;
			_overflow_cache = 0;
		}

		void Normalize() {
			if (uint32(low) < LOW_THRESHOLD || uint32(low >> 32)) {
				f->Put(byte(_overflow_cache + byte(low >> 32)));
				while (--_ff != 0) {
					f->Put(byte(0xFF + byte(low >> 32)));
				}

				_overflow_cache = byte(uint32(low) >> 24);
			}

			_ff++;
			low = (uint32)low << 8;
		}

		void Flush() {
			for (int i = 0; i < 5; i++) {
				Normalize();
			}

			f->Flush();
		}

		void Encode12P(uint16 p, int y) {
			ASSERT_DEBUG(p > 0 && p < 0xFFF);

			uint32 bound = p * (range >> 12);
			if (y) {
				range = bound;
			}
			else {
				low += bound;
				range -= bound;
			}

			while (range < BitDecoder::TOP) {
				range <<= 8;
				Normalize();
			}
		}

		void EncodeBits(uint32 v, int num_bits) {
			ASSERT_DEBUG(num_bits > 0);
			ASSERT_DEBUG(v < (1 << num_bits));

			for (num_bits--; num_bits >= 0; num_bits--) {
				range >>= 1;
				low += range & (0 - ((v >> num_bits) & 1));
				if (range < BitDecoder::TOP) {
					range <<= 8;
					Normalize();
				}
			}
		}
	};

	INLINE uint32 get_match_min(uint32 dist) {
		uint32 match_min = 2;
		if (dist & ~0xFF) { ++match_min; }
		if (dist & ~0xFFFF) { ++match_min; }
		if (dist & ~0xFFFFFF) { ++match_min; }

		return match_min;
	}

	const static int MATCH_MIN_BASE = 2;
	const static int MATCH_LIMIT = MATCH_MIN_BASE + (1 << 12);

	const static int COUNTER_SCALE_BITS = 12;
	const static int COUNTER_SCALE = 1 << COUNTER_SCALE_BITS;
	const static int COUNTER_MASK = COUNTER_SCALE - 1;
	const static int COUNTER_SCALE_HALF = 1 << (COUNTER_SCALE_BITS - 1);

	const static int COUNTER_ADAPT_SPEED = 5;

	struct BitCounter {
		uint16 c;

		void Init() {
			c = COUNTER_SCALE_HALF;
		}

		inline uint32 Price(int b) const {
			if (b) {
				return log2_lut[c >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
			}
			else {
				return log2_lut[(COUNTER_MASK ^ c) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
			}
		}

		inline void Encode(int b, BitCoder& code) {
			code.Encode12P(c, b);
			if (b) {
				c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
			}
			else {
				c -= c >> COUNTER_ADAPT_SPEED;
			}
		}

		inline int Decode(BitDecoder& code) {
			int b = code.Decode12P(c);
			if (b) {
				c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
			}
			else {
				c -= c >> COUNTER_ADAPT_SPEED;
			}

			return b;
		}
	};

	template<int NUM_BITS>
	struct BitTree {
		uint16 tree[1 << NUM_BITS];

		void Init() {
			for (int i = 0; i < (1 << NUM_BITS); i++) { tree[i] = COUNTER_SCALE_HALF; }
		}

		uint32 Price(int y) const {
			uint32 sum = 0;
			int p = 1;
			for (int j = NUM_BITS - 1; j >= 0; j--) {
				int b = (y >> j) & 1;

				const uint16& c = tree[p - 1];
				if (b) {
					sum += log2_lut[c >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
				}
				else {
					sum += log2_lut[(COUNTER_MASK ^ c) >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
				}

				p = (p << 1) + b;
			}

			return sum;
		}

		uint32 PriceReverse(int y) const {
			uint32 sum = 0;
			int p = 1;
			for (int j = 0; j < NUM_BITS; j++) {
				int b = (y >> j) & 1;

				const uint16& c = tree[p - 1];
				sum += log2_lut[c >> (COUNTER_SCALE_BITS - LOG2_LUT_SIZE_BITS)];
				p = (p << 1) + b;
			}

			return sum;
		}

		void Encode(int y, BitCoder& code) {
			int p = 1;
			for (int j = NUM_BITS - 1; j >= 0; j--) {
				int b = (y >> j) & 1;

				uint16& c = tree[p - 1];
				code.Encode12P(c, b);
				if (b) {
					c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
				}
				else {
					c -= c >> COUNTER_ADAPT_SPEED;
				}

				p = (p << 1) + b;
			}
		}

		void EncodeReverse(int y, BitCoder& code) {
			int p = 1;
			for (int j = 0; j < NUM_BITS; j++) {
				int b = (y >> j) & 1;

				uint16& c = tree[p - 1];
				code.Encode12P(c, b);
				if (b) {
					c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
				}
				else {
					c -= c >> COUNTER_ADAPT_SPEED;
				}

				p = (p << 1) + b;
			}
		}

		int Decode(BitDecoder& code) {
			int p = 1;
			for (int j = NUM_BITS - 1; j >= 0; j--) {
				uint16& c = tree[p - 1];
				int b = code.Decode12P(c);
				if (b) {
					c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
				}
				else {
					c -= c >> COUNTER_ADAPT_SPEED;
				}

				p = (p << 1) + b;
			}

			return p - (1 << NUM_BITS);
		}

		int DecodeReverse(BitDecoder& code) {
			int p = 1;
			int sym = 0;
			for (int j = 0; j < NUM_BITS; j++) {
				uint16& c = tree[p - 1];
				int b = code.Decode12P(c);
				if (b) {
					c += (COUNTER_MASK ^ c) >> COUNTER_ADAPT_SPEED;
				}
				else {
					c -= c >> COUNTER_ADAPT_SPEED;
				}

				p = (p << 1) + b;
				sym += b << j;
			}

			return sym;
		}
	};

	const static int LITERAL_MAX_POSITION_BITS = 4;
	const static int LITERAL_MAX_POSITION_SIZE = 1 << LITERAL_MAX_POSITION_BITS;

	struct State {
		BitTree<8>* lit;
		BitCounter is_match[4][LITERAL_MAX_POSITION_SIZE];
		BitCounter is_rep[4][LITERAL_MAX_POSITION_SIZE];

		BitCounter len_marker[3];
		BitTree<3> len_low[LITERAL_MAX_POSITION_SIZE];
		BitTree<3> len_mid[LITERAL_MAX_POSITION_SIZE];
		BitTree<8> len_high;
		BitTree<6> dist_slots[4];
		BitTree<4> dist_low;
		int rep0;
		byte match_hist;

		uint32 cached_price_match[LITERAL_MAX_POSITION_SIZE][2];
		uint32 cached_price_rep[LITERAL_MAX_POSITION_SIZE][2];
		uint32* cached_price_lit;
		uint32 cached_price_length[273];
		uint32 cached_price_dist_slot[4][64];
		uint32 cached_price_dist_low[16];

		int pos_align_bits, context_bits;

		void Init(int pos_align_bits, int context_bits) {
			this->pos_align_bits = pos_align_bits;
			this->context_bits = context_bits;

			ASSERT(pos_align_bits <= 4);
			ASSERT(context_bits <= 8);
			ASSERT(pos_align_bits + context_bits <= 8);

			lit = new BitTree<8>[(1 << pos_align_bits) << context_bits];
			cached_price_lit = new uint32[(1 << pos_align_bits) << (8 + context_bits)];

			match_hist = 0;
			rep0 = 0;

			int pos_align_size = 1 << pos_align_bits;
			int context_size = 1 << context_bits;

			for (int i = 0; i < pos_align_size; i++) {
				for (int j = 0; j < context_size; j++) {
					lit[(i << context_bits) + j].Init();
				}
			}

			for (int s = 0; s < 4; s++) {
				for (int i = 0; i < pos_align_size; i++) {
					is_match[s][i].Init();
				}

				for (int i = 0; i < pos_align_size; i++) {
					is_rep[s][i].Init();
				}
			}

			for (int i = 0; i < 3; i++) { len_marker[i].Init(); }
			for (int i = 0; i < pos_align_size; i++) { len_low[i].Init(); }
			for (int i = 0; i < pos_align_size; i++) { len_mid[i].Init(); }
			len_high.Init();

			for (int i = 0; i < 4; i++) { dist_slots[i].Init(); }
			dist_low.Init();

			printf("State: %d KB\n", uint32(sizeof(State) +
				((sizeof(BitTree<8>) << pos_align_bits) << context_bits) +
				((sizeof(uint32) << pos_align_bits) << (8 + context_bits)) + 1023) >> 10);
		}

		void Release() {
			delete[] lit;
			delete[] cached_price_lit;
		}

		void NextStateMatch() { match_hist = ((match_hist << 1) + 1) & 3; }
		void NextStateLiteral() { match_hist = (match_hist << 1) & 3; }

		void UpdatePriceCache() {
			int pos_align_size = 1 << pos_align_bits;

			for (int pb = 0; pb < pos_align_size; ++pb) {
				for (int b = 0; b < 2; b++) {
					cached_price_match[pb][b] = is_match[match_hist][pb].Price(b);
					cached_price_rep[pb][b] = is_rep[match_hist][pb].Price(b);
				}
			}

			memset(cached_price_lit, -1, sizeof(uint32) * (1 << pos_align_bits) << (8 + context_bits));

			uint32 frag_lp0 = len_marker[0].Price(1);
			uint32 frag_lp1 = len_marker[1].Price(1);
			uint32 lp0 = len_marker[0].Price(0);
			uint32 lp1 = frag_lp0 + len_marker[1].Price(0);
			uint32 lp2 = frag_lp0 + frag_lp1 + len_marker[2].Price(0);
			uint32 lp3 = frag_lp0 + frag_lp1 + len_marker[2].Price(1);
			for (int pb = 0; pb < pos_align_size; ++pb) {
				for (int lv = 0; lv < 273; lv++) {
					if (lv < 8) {
						cached_price_length[lv] = lp0 + len_low[pb].Price(lv);
					}
					else if (lv < 16) {
						cached_price_length[lv] = lp1 + len_mid[pb].Price(lv - 8);
					}
					else if (lv < 272) {
						cached_price_length[lv] = lp2 + len_high.Price(lv - 16);
					}
					else {
						cached_price_length[lv] = lp3 + (12 << LOG2_LUT_SCALE_BITS);
					}
				}
			}

			for (int slot = 0; slot < 4; slot++) {
				for (int y = 0; y < 64; y++) {
					cached_price_dist_slot[slot][y] = dist_slots[slot].Price(y);
				}
			}

			for (int y = 0; y < 16; y++) {
				cached_price_dist_low[y] = dist_low.Price(y);
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

			return cached_price_match[pb][0] + lit_price;
		}

		uint32 PriceRep0Cached(uint32 lv, int pb, int ctx) {
			uint32 sum = cached_price_match[pb][1] + cached_price_rep[pb][1];
			sum += cached_price_length[Min(lv, 272U)];

			return sum;
		}

		uint32 PriceMatchCached(uint32 lv, uint32 dv, int pb, int ctx) {
			ASSERT(dv < 0x7FFFFFFF);
			uint32 sum = cached_price_match[pb][1] + cached_price_rep[pb][0];
			sum += cached_price_length[Min(lv, 272U)];

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

	static inline uint8 value8(const byte* b) { return *b; }
	static inline uint16 value16(const byte* b) { return *reinterpret_cast<const uint16*>(b); }
	static inline uint32 value24(const byte* b) { return *reinterpret_cast<const uint32*>(b) & 0xFFFFFFU; }
	static inline uint32 value32(const byte* b) { return *reinterpret_cast<const uint32*>(b); }
	static inline uint64 value64(const byte* b) { return *reinterpret_cast<const uint64*>(b); }

	static inline uint32 hash4(uint32 x) { return x * 123456791U; }

	static int match_finder_padding = 16;

	INLINE static uint try_match(const byte* buf, uint32 mp, uint32 p, uint32 p_end) {
		uint len = 0;

#ifdef MSVC
		unsigned long r = 0;
		const uint64* p0 = reinterpret_cast<const uint64*>(buf + p);
		const uint64* p1 = reinterpret_cast<const uint64*>(buf + mp);
		while (p + len < p_end) {
			uint64 v = *p0 ^ *p1;
			if (v) {
				_BitScanForward64(&r, v);
				len += r >> 3;
				if (p + len > p_end) {
					len = p_end - p;
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
		while (p + len + 4 < p_end && *reinterpret_cast<const uint32*>(buf + p + len) == *reinterpret_cast<const uint32*>(buf + mp + len)) {
			len += 4;
		}

		while (p + len < p_end && buf[p + len] == buf[mp + len]) {
			++len;
		}
#endif

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
			uint32* head_end = head + (1 << hash_bits);

			while (row_head < head_end) {
				*row_head = (*row_head < shift || *row_head == uint32(-1)) ? uint32(-1) : (*row_head - shift);
				++row_head;
			}

			uint32 child_offset = shift << 1;
			uint32* tree_head = tree + child_offset;
			uint32* tree_end = tree + (window_size << 1U);
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

			ASSERT((p << 1) < (2 * window_size));
			uint32* right = tree + (p << 1);
			uint32* left = right + 1;

			uint left_len = 0, right_len = 0;

			uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
			uint32 sp = head[h4];
			head[h4] = p;

			while (sp != uint32(-1) && p > sp
				&& p - sp <= max_dist) {
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

			ASSERT((p << 1) < (2 * window_size));
			uint32* right = tree + (p << 1);
			uint32* left = right + 1;

			uint left_len = 0, right_len = 0;

			uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
			uint32 sp = head[h4];
			head[h4] = p;

			while (sp != uint32(-1) && p > sp
				&& p - sp <= max_dist) {
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

			ASSERT((p << 1) < (2 * window_size));
			const uint32* right = tree + (p << 1);
			const uint32* left = right + 1;

			uint left_len = 0, right_len = 0;

			uint32 h4 = hash4(value32(buf + p)) >> (32 - hash_bits);
			uint32 sp = head[h4];

			while (sp != uint32(-1) && p > sp
				&& p - sp <= max_dist) {
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

	const static uint32 TABLE_MATCH_SLOTS = 18;
	struct ParseNode {
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
		int dry_count, dry_step;

		MatcherOpt(uint32 window_size) {
			mf2 = HR2(10);
			mf3 = HR3(14, 2);
			mf4 = BT4(window_size, 16);
			dry_count = dry_step = 0;
		}

		void Init() {
			mf2.Init();
			mf3.Init();
			mf4.Init();
			dry_count = dry_step = 0;
		}

		void Release() {
			mf2.Release();
			mf3.Release();
			mf4.Release();
		}

		void ShiftLeft(int shift) {
			mf2.ShiftLeft(shift);
			mf3.ShiftLeft(shift);
			mf4.ShiftLeft(shift);
		}

		void UpdateBlock(ParseNode* opt_table, const byte* window, int start, int end, int peek_end, int history_size) {
			for (int p = start; p < peek_end; ++p) {
				ParseNode* node = opt_table + p - start;
				node->max_len = MATCH_MIN_BASE - 1;
				node->max_match_pos = -1;
				memset(node->match_pos, -1, sizeof(uint32) * TABLE_MATCH_SLOTS);
			}

			int run_step = 1;
			int p = start;
			while (p < end) {
				ParseNode* node = opt_table + p - start;

				const uint32* result_set[3] = {
					mf2.FindAndUpdate(window, p, peek_end, history_size),
					mf3.FindAndUpdate(window, p, peek_end, history_size),
					mf4.FindAndUpdate(window, p, peek_end, history_size)
				};

				for (int i = 0; i < 3; i++) {
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
						ParseNode* n2 = opt_table + p - start;

						if (can_expand && p + len < peek_end && (can_expand = window[mp + len] == window[p + len])) {
							++len;
						}

						ASSERT_MATCH(!memcmp(window + p, window + mp, len));

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
				ParseNode* node = opt_table + p - start;
				const uint32* result_set[3] = {
					mf2.Find(window, p, peek_end, history_size),
					mf3.Find(window, p, peek_end, history_size),
					mf4.Find(window, p, peek_end, history_size)
				};

				for (int i = 0; i < 3; i++) {
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
						ParseNode* n2 = opt_table + p - start;

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

		MatcherFast(uint32 window_size) {
			mf2 = HR2(10);
			mf3 = HR3(14, 4);
			mf4 = HT4(16, 32);
			//mf4 = HC4(window_size, 16, 24, 4);
			//mf4 = BT4(window_size, 16, 8);
		}

		void Init() {
			mf2.Init();
			mf3.Init();
			mf4.Init();
		}

		void Release() {
			mf2.Release();
			mf3.Release();
			mf4.Release();
		}

		void ShiftLeft(int shift) {
			mf2.ShiftLeft(shift);
			mf3.ShiftLeft(shift);
			mf4.ShiftLeft(shift);
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
		}

		void FindAndUpdate(const byte* window, int p, int end, int history_size, uint32& best_pos, uint32& best_len) {
			const uint32* result_set[3] = {
				mf2.FindAndUpdate(window, p, end, history_size),
				mf3.FindAndUpdate(window, p, end, history_size),
				mf4.FindAndUpdate(window, p, end, history_size)
			};

			best_len = 0;
			uint32 cur_lv = 0;

			for (int i = 0; i < 3; i++) {
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
		uint32 cost, match_pos, match_len;
	};

	void opt_parse_path(const byte* window, ParsePath* path, const ParseNode* opt_table, State& state, int opt_start_at, int block_end) {
		path[block_end - opt_start_at].cost = 0;
		path[block_end - opt_start_at].match_len = 0;
		path[block_end - opt_start_at].match_pos = uint32(-1);

		int op = block_end - 1;
		while (op >= opt_start_at) {
			int pb = state.pos_align_bits ? (op & state.pos_align_bits) : 0;
			int ctx = op > 0 ? window[op - 1] : 0;

			const ParseNode* opt = opt_table + op - opt_start_at;
			ParsePath* cn = path + op - opt_start_at;

			uint32 cost_literal = state.PriceLiteralCached(window[op], pb, ctx);
			uint32 cost_next = path[op - opt_start_at + 1].cost + cost_literal;

			uint32 max_match_pos = opt->max_match_pos;
			int max_len = opt->max_len;

			if (max_len >= MATCH_LIMIT) {
				max_len = MATCH_LIMIT - 1;
			}

			if (op + max_len > block_end) {
				max_len = block_end - op;
			}

			cn->cost = cost_next;
			cn->match_len = 0;
			cn->match_pos = uint32(-1);

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

				uint32 cost_match = state.PriceMatchCached(lv, dv, pb, ctx);
				cost_match += path[op - opt_start_at + test_len].cost;
				if (cost_match < cn->cost) {
					cn->cost = cost_match;
					cn->match_len = test_len;
					cn->match_pos = match_pos;
				}

				if (test_len >= MatcherOpt::NICE_MATCH) {
					break;
				}

				test_len -= test_step;
			}
			--op;
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

		ParseNode* match_table = new ParseNode[opt_block_after];
		memset(match_table, 0, sizeof(ParseNode) * opt_block_after);

		ParsePath* path = new ParsePath[opt_block_after + 1];
		memset(path, 0, sizeof(ParsePath) * (opt_block_after + 1));

		printf("Window size: %d KB\n", (window_size + match_finder_padding + 1023) >> 10);
		printf("Parser: %d KB\n", uint32((sizeof(ParsePath) + sizeof(ParseNode)) * opt_block_after + 1023) >> 10);

		BufferedFile buf_file(fout);
		BitCoder code(&buf_file);
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
		//uint64 total_sum = 0;

		char status_buf[70];
		clock_t last_message_time = 0, start_time = clock();
		snprintf(status_buf, sizeof(status_buf), "Working...");
		printf("%-79s\r", status_buf);

		while (true) {
			if (!window_length || p >= slide_size) {
				if (window_length >= buffer_io_size) {
					memmove(window, window + buffer_io_size, window_length - buffer_io_size);
					matcher.ShiftLeft(buffer_io_size);
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
					opt_parse_path(window, path, match_table, state, opt_start_at, block_end);
					ASSERT(!path[block_end - opt_start_at].cost);
					ASSERT(!path[block_end - opt_start_at].match_len);
				}

				int pb = state.pos_align_bits ? (p & state.pos_align_bits) : 0;
				int ctx = p > 0 ? window[p - 1] : 0;
				int y = window[p];

				uint32 match_pos = path[p - opt_start_at].match_pos;
				uint32 match_len = path[p - opt_start_at].match_len;

				// Comparison against lazy parser decisions
				//uint32 match_pos = match_table[p - opt_start_at].max_match_pos;
				//uint32 match_len = match_table[p - opt_start_at].max_len;
				//if (match_table[p - opt_start_at + 1].max_len > match_len ||
				//	match_table[p - opt_start_at + 2].max_len > match_len) {
				//	match_len = 0;
				//}

				uint32 dv = p - match_pos - 1;
				uint32 min_match_len = get_match_min(dv);

				int rep = 0;
				if (state.rep0 > 0) {
					int rp = p - state.rep0;
					int rlen = try_match(window, rp, p, p_end);
					if (p - opt_start_at + rlen > opt_block_after) { rlen = opt_block_after - p + opt_start_at; }
					if (rlen >= MATCH_LIMIT) { rlen = MATCH_LIMIT - 1; }

					if (rlen >= MATCH_MIN_BASE) {
						ASSERT_MATCH(p - opt_start_at + rlen <= opt_block_after);

						int c0 = path[p - opt_start_at].cost;
						int c1 = state.PriceRep0Cached(rlen - MATCH_MIN_BASE, pb, ctx) + path[p - opt_start_at + rlen].cost;

						if (c1 < c0) {
							min_match_len = MATCH_MIN_BASE;
							match_len = rlen;
							rep = 1;
						}
					}
				}

				int adv = 1;
				if (rep) {
					ASSERT(state.rep0 > 0);
					match_pos = p - state.rep0;

					ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

					adv = match_len;

					uint32 lv = match_len - min_match_len;
					//total_sum += state.PriceRep0(lv, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(1, code);
					state.is_rep[state.match_hist][pb].Encode(1, code);
					state.NextStateMatch();
					if (lv < 8) {
						state.len_marker[0].Encode(0, code);
						state.len_low[pb].Encode(lv, code);
					}
					else if (lv < 16) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(0, code);
						state.len_mid[pb].Encode(lv - 8, code);
					}
					else if (lv < 272) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(0, code);
						state.len_high.Encode(lv - 16, code);
					}
					else {
						ASSERT(lv < 4368);
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(1, code);
						code.EncodeBits(lv - 272, 12);
					}
				}
				else if (match_pos < p && min_match_len <= match_len) {
					ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

					adv = match_len;
					state.rep0 = p - match_pos;

					uint32 lv = match_len - min_match_len;
					//total_sum += state.PriceMatch(lv, dv, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(1, code);
					state.is_rep[state.match_hist][pb].Encode(0, code);
					state.NextStateMatch();
					if (lv < 8) {
						state.len_marker[0].Encode(0, code);
						state.len_low[pb].Encode(lv, code);
					}
					else if (lv < 16) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(0, code);
						state.len_mid[pb].Encode(lv - 8, code);
					}
					else if (lv < 272) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(0, code);
						state.len_high.Encode(lv - 16, code);
					}
					else {
						ASSERT(lv < 4368);
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(1, code);
						code.EncodeBits(lv - 272, 12);
					}

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
				else {
					auto& lit_coder = state.lit[state.ToLitIndex(pb, ctx)];
					//total_sum += state.PriceLiteral(y, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(0, code);
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

		const int opt_block_size = 0x1000;
		const int opt_block_after = opt_block_size + 0x100;

		byte* window = new byte[window_size + match_finder_padding];
		memset(window, 0, window_size + match_finder_padding);

		State state;
		state.Init(pos_align_bits, context_bits);

		MatcherFast matcher(window_size);
		matcher.Init();

		printf("Window size: %d KB\n", (window_size + match_finder_padding + 1023) >> 10);
		printf("State: %d KB\n", uint32(sizeof(state) + 1023) >> 10);

		BufferedFile buf_file(fout);
		BitCoder code(&buf_file);
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
		//uint64 total_sum = 0;

		char status_buf[70];
		clock_t last_message_time = 0, start_time = clock();
		snprintf(status_buf, sizeof(status_buf), "Working...");
		printf("%-79s\r", status_buf);

		while (true) {
			if (!window_length || p >= slide_size) {
				if (window_length >= buffer_io_size) {
					memmove(window, window + buffer_io_size, window_length - buffer_io_size);
					matcher.ShiftLeft(buffer_io_size);
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

				uint32 dv = p - match_pos - 1;
				uint32 min_match_len = get_match_min(dv);

				int rep = 0;
				if (state.rep0 > 0) {
					int rp = p - state.rep0;
					int rlen = try_match(window, rp, p, p_end);
					if (p - opt_start_at + rlen > opt_block_after) { rlen = opt_block_after - p + opt_start_at; }
					if (rlen >= MATCH_LIMIT) { rlen = MATCH_LIMIT - 1; }

					if (rlen >= MATCH_MIN_BASE) {
						ASSERT_MATCH(p - opt_start_at + rlen <= opt_block_after);

						if (rlen > match_len + 1) {
							min_match_len = MATCH_MIN_BASE;
							match_len = rlen;
							rep = 1;
						}
					}
				}

				int adv = 1;
				if (rep) {
					ASSERT(state.rep0 > 0);
					match_pos = p - state.rep0;

					ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

					adv = match_len;

					uint32 lv = match_len - min_match_len;
					//total_sum += state.PriceRep0(lv, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(1, code);
					state.is_rep[state.match_hist][pb].Encode(1, code);
					state.NextStateMatch();
					if (lv < 8) {
						state.len_marker[0].Encode(0, code);
						state.len_low[pb].Encode(lv, code);
					}
					else if (lv < 16) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(0, code);
						state.len_mid[pb].Encode(lv - 8, code);
					}
					else if (lv < 272) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(0, code);
						state.len_high.Encode(lv - 16, code);
					}
					else {
						ASSERT(lv < 4368);
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(1, code);
						code.EncodeBits(lv - 272, 12);
					}
				}
				else if (match_pos < p && min_match_len <= match_len) {
					ASSERT_MATCH(!memcmp(window + p, window + match_pos, match_len));

					adv = match_len;
					state.rep0 = p - match_pos;

					uint32 lv = match_len - min_match_len;
					//total_sum += state.PriceMatch(lv, dv, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(1, code);
					state.is_rep[state.match_hist][pb].Encode(0, code);
					state.NextStateMatch();
					if (lv < 8) {
						state.len_marker[0].Encode(0, code);
						state.len_low[pb].Encode(lv, code);
					}
					else if (lv < 16) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(0, code);
						state.len_mid[pb].Encode(lv - 8, code);
					}
					else if (lv < 272) {
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(0, code);
						state.len_high.Encode(lv - 16, code);
					}
					else {
						ASSERT(lv < 4368);
						state.len_marker[0].Encode(1, code);
						state.len_marker[1].Encode(1, code);
						state.len_marker[2].Encode(1, code);
						code.EncodeBits(lv - 272, 12);
					}

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
				else {
					auto& lit_coder = state.lit[state.ToLitIndex(pb, ctx)];
					//total_sum += state.PriceLiteral(y, pb, ctx);

					state.is_match[state.match_hist][pb].Encode(0, code);
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

	void decompress(FILE* fin, FILE* fout, uint32& unpacked_crc32, uint32& running_crc32) {
		BufferedFile buf_file(fin);
		BitDecoder code(&buf_file);
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

			int is_match = state.is_match[state.match_hist][pb].Decode(code);

			int adv = 1;
			if (is_match) {
				int is_rep = state.is_rep[state.match_hist][pb].Decode(code);
				state.NextStateMatch();

				uint32 lv = 0;
				if (state.len_marker[0].Decode(code) == 0) {
					lv = state.len_low[pb].Decode(code);
				}
				else if (state.len_marker[1].Decode(code) == 0) {
					lv = 8 + state.len_mid[pb].Decode(code);
				}
				else if (state.len_marker[2].Decode(code) == 0) {
					lv = 16 + state.len_high.Decode(code);
				}
				else {
					lv = 272 + code.DecodeBits(12);
				}

				int dist, length;
				if (is_rep) {
					dist = state.rep0;
					length = lv + MATCH_MIN_BASE;
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

					state.rep0 = dist;
				}
				ASSERT(wp - dist >= 0 && wp + length <= window_size);

				adv = length;
				while (length-- > 0) {
					int y = window[wp - dist];
					window[wp++] = y;
				}
			}
			else {
				state.NextStateLiteral();

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
}

int main(int argc, char** argv) {
	Release_NLZM::log2_init();
	Release_NLZM::crc32_init();

	FILE* fin = nullptr;
	FILE* fout = nullptr;

	int window_bits = 20;
	int pos_bits = 0;
	int ctx1_bits = 2;

	while (argc >= 2 && *argv[1] == '-') {
		char* arg = argv[1];
		argv++;
		argc--;

		Release_NLZM::_tolower(arg);
		if (!strncmp(arg, "-window:", 8)) {
			window_bits = Release_NLZM::Clamp(atoi(arg + 8), 16, 30);
			printf("Window bits: %d\n", window_bits);
		}
		else if (!strncmp(arg, "-posbits:", 9)) {
			pos_bits = Release_NLZM::Clamp(atoi(arg + 9), 0, 4);
			printf("Position bits: %d\n", pos_bits);
		}
		else if (!strncmp(arg, "-ctx1:", 6)) {
			ctx1_bits = Release_NLZM::Clamp(atoi(arg + 6), 2, 8);
			printf("ctx1 bits: %d\n", ctx1_bits);
		}
		else {
			printf("Unrecognized flag %s\n", arg);
			return -1;
		}
	}

	if (argc == 4 && Release_NLZM::_tolower(argv[1][0]) == 'c') {
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
			printf("Error: %s file does not exist\n", argv[2]);
			fclose(fin);
			return -1;
		}

		printf("Calculating crc... ");
		uint32 crc32 = Release_NLZM::crc_file(fin);
		printf("%X\n", crc32);
		fseek(fin, 0, SEEK_SET);

		if (Release_NLZM::_tolower(argv[1][1]) == 'f') {
			Release_NLZM::compress_fast(fin, fout, window_bits, pos_bits, ctx1_bits, crc32);
		}
		else {
			Release_NLZM::compress_opt(fin, fout, window_bits, pos_bits, ctx1_bits, crc32);
		}

		fclose(fin);
		fclose(fout);
	}
	else if (argc == 4 && Release_NLZM::_tolower(argv[1][0]) == 'd') {
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
			printf("Error: %s file does not exist\n", argv[2]);
			fclose(fin);
			return -1;
		}

		uint32 unpacked_crc32 = 0, running_crc32 = 0;
		Release_NLZM::decompress(fin, fout, unpacked_crc32, running_crc32);
		fclose(fin);
		fclose(fout);

		printf("Expected CRC32 %X, was %X\n", unpacked_crc32, running_crc32);
		return unpacked_crc32 == running_crc32;
	}
	else if (argc == 3 && Release_NLZM::_tolower(argv[1][0]) == 't') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		uint32 unpacked_crc32 = 0, running_crc32 = 0;
		Release_NLZM::decompress(fin, nullptr, unpacked_crc32, running_crc32);
		fclose(fin);

		printf("Expected CRC32 %X, was %X\n", unpacked_crc32, running_crc32);
		return unpacked_crc32 == running_crc32;
	}
	else if (argc == 3 && Release_NLZM::_tolower(argv[1][0]) == 'h') {
		fin = fopen(argv[2], "rb");
		if (!fin) {
			printf("Error: %s file does not exist\n", argv[2]);
			return -1;
		}

		printf("Calculating crc... ");
		uint32 crc32 = Release_NLZM::crc_file(fin);
		fclose(fin);
		printf("%X\n", crc32);
	}
	else {
		printf("NLZM 64-bit - Written by Nauful\n"
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
