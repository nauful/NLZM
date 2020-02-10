NLZM is a general-purpose file compressor that implements mostly-optimal parsing LZ with a sliding window, exhaustive dictionaries up to 1 GB, long-range match finding, cyclical + o0 + o1 context modelling, and a context-based statistical model for decision coding.

Results

enwik8.txt to 25,226,391 bytes and decompresses in 0.95 seconds.

enwik9.txt to 203,649,613 bytes and decompresses in 8.64 seconds.

1.02 update:
- Added long-range matchfinder (RK256, Rabin-Karp with 256 byte indexing and rolling matcher).
- Simplified probability-coded decision sequences with CDFs.
- Replaced literal/match decision with commands (unlimited operation types vs just two).
- Improved cost estimation, rep model.
- Started on forward-pass optimal parse, to eventually replace reverse optimal parse, for better handling stateful models (e.g. better rep, LZP, ROLZ).

1.01 update:
- Replaced range coding stream with chunks of two-state rANS + direct bit stream to improve decompression speed (~50-70% decompression speed improvement).
- Chunks are 16 KB, end after roughly 15 KB when direct bits counted + log2[rANS freq] approaches limit.
- Replaced bit tree with nibbles (2x4 for byte literals, 2x3 for 64-bit dist slots, etc) + CDF.
- SIMD for 3-bit, 4-bit CDF updates and look-up (compile with -DUSE_SSE).
- Branch free counter (~5% decompression speed improvement)
- Limit BT4 match attempts to 128 per position.

Semi-optimal parsing is done by finding all matches (exhaustively) for the block size, then finding (but not updating search) matches for block size + nice length, then traced by finding the minimum cost from block size + nice length to block start. Not terminating matches at block size allows some (better) parsing decisions when matches could cross block boundaries. Block start is always aligned to block size so that the match finder is updated with every position. When a decision table is calculated for a block, the current costs of encoding (literal cost, match cost, etc) are estimated (log2[prob] bps, ~1-3% error compared to final compressed file size with 8k blocks) and used for path finding.

There are two optimizations for semi-optimal parsing:
- In the match finder (forward pass), if a sufficiently long match is found, the next position will try to continue that match or reduce length by one, and ends when length is less than an escape threshold to not spend extra time in the binary tree match finder.
- In the path finder (reverse pass), a maximum of 8 different lengths forward (step size of (max length - min length) / 8) are tested to not waste time exhausting all the options when a few represent most of the optimization.

To speed up passing incompressible data, if no matches are found for N bytes, then every 2nd byte is tested, 4th, 8th until another match is found and the step size is decreased again. Speed drops off rapidly with larger dictionary sizes and semi-optimal parsing due to the exhaustive match finder.

Most of my target files are plain-text natural language texts, structured binary data (sensor data, rows and hierarchies of floats, doubles, ints, repetitive strings and labels) between 300 KB and 200 MB, compressed once and decompressed often. Default window size is 1 MB and default memory usage during decompression is not significantly more because my target environment only has ~4-10 MB RAM available for decompression and ~20-32 MB for compression. Memory for compression is 12N bytes (optimal) + IO overhead or 9 MB (fast, reduced from 17 MB) + N + IO overhead.
