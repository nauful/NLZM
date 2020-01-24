NLZM is a general-purpose file compressor that implements semi-optimal parsing LZ with a sliding window, exhaustive dictionaries up to 1 GB, cyclical + o0 + o1 context modelling, and a context-based statistical model for branch decisions.

Results (c to compress, cf to compress fast, d to decompress, t to test in-memory decompression and CRC32):
```
nlzm -window:27 c enwik8.txt out8
100000000 -> 25116875 in 102.080 sec

nlzm -window:30 c enwik9 out9
1000000000 -> 202561110 in 10467.483 sec

nlzm -window:27 c enwik9 out9_27
1000000000 -> 211524600 in 1427.124 sec

nlzm t out8
25116875 in 1.637 sec

nlzm t out9_27
211524600 in 15.529 sec

nlzm -window:27 cf enwik8.txt out8f
100000000 -> 33460723 in 9.730 sec

nlzm -window:30 cf enwik9 out9f
1000000000 -> 293022504 in 89.565 sec
```

Semi-optimal parsing is done by finding all matches (exhaustively) for the block size, then finding (but not updating search) matches for block size + nice length, then traced by finding the minimum cost from block size + nice length to block start. Not terminating matches at block size allows some (better) parsing decisions when matches could cross block boundaries. Block start is always aligned to block size so that the match finder is updated with every position. When a decision table is calculated for a block, the current costs of encoding (literal cost, match cost, etc) are estimated (log2[prob] bps, ~1-3% error compared to final compressed file size with 8k blocks) and used for path finding.

There are two optimizations for semi-optimal parsing:
- In the match finder (forward pass), if a sufficiently long match is found, the next position will try to continue that match or reduce length by one, and ends when length is less than an escape threshold to not spend extra time in the binary tree match finder.
- In the path finder (reverse pass), a maximum of 8 different lengths forward (step size of (max length - min length) / 8) are tested to not waste time exhausting all the options when a few represent most of the optimization.

To speed up passing incompressible data, if no matches are found for N bytes, then every 2nd byte is tested, 4th, 8th until another match is found and the step size is decreased again. Speed drops off rapidly with larger dictionary sizes and semi-optimal parsing due to the exhaustive match finder.

Most of my target files are plain-text natural language texts, structured binary data (sensor data, rows and hierarchies of floats, doubles, ints, repetitive strings and labels) between 300 KB and 200 MB, compressed once and decompressed often. Default window size is 1 MB and default memory usage during decompression is not significantly more because my target environment only has ~4-10 MB RAM available for decompression and ~20-32 MB for compression. Memory for compression is 12N bytes (optimal) + IO overhead or 9 MB (fast, reduced from 17 MB) + N + IO overhead.
