NLZM is a nibbled ANS general-purpose file compressor that implements mostly-optimal parsing LZ with a sliding window, exhaustive dictionaries up to 256 MB, long-range match finding and a context-based statistical model for decision coding.

#### Results
- enwik8.txt to 25,056,894 bytes and decompresses in 0.81 seconds.
- enwik9.txt to 206,548,020 bytes and decompresses in 7.16 seconds.

##### 1.03 update:
- New forward graph parser with early exits.
- Replaced sliding dictionary with a ring dictionary + lookahead to reduce memory usage and simplify matchers.
- Improved model for long-range match finder, carrying matches and repeated offsets.
- General performance improvements (~12% decompression speed improvement).
- Greatly simplified source, reduced lines of code by ~40%.

##### 1.02 update:
- Added long-range matchfinder (RK256, Rabin-Karp with 256 byte indexing and rolling matcher).
- Simplified probability-coded decision sequences with CDFs.
- Replaced literal/match decision with commands (unlimited operation types vs just two).
- Improved cost estimation, rep model.
- Started on forward-pass optimal parse, to eventually replace reverse optimal parse, for better handling stateful models (e.g. better rep, LZP, ROLZ).

#### 1.01 update:
- Replaced range coding stream with chunks of two-state rANS + direct bit stream to improve decompression speed (~50-70% decompression speed improvement).
- Chunks are 16 KB, end after roughly 15 KB when direct bits counted + log2[rANS freq] approaches limit.
- Replaced bit tree with nibbles (2x4 for byte literals, 2x3 for 64-bit dist slots, etc) + CDF.
- SIMD for 3-bit, 4-bit CDF updates and look-up (compile with -DUSE_SSE).
- Branch free counter (~5% decompression speed improvement)
- Limit BT4 match attempts to 128 per position.
