Based on:
2018-07-04 : Igor Pavlov : Public domain
PPMd var.I (2002): Dmitry Shkarin

Allocator:
- Use doubly linked list to manage memory blocks.
- Use a free bit instead of stamp to mark blocks as free or in use.
- Block list root doesn't count numbers of units.
- Block list allocation is included in the units space.
- Freeing units defragments, replacing GlueFreeBlocks.
- num_units_available tracks the total number of available units.
- Text area shrinking moves to unit allocation time.
- Simplified text area expansion.

- Cleaned up and simplified model handling.
- Moved encoding and decoding methods to PPMD::EncodeFile and PPMD::DecodeFile.
- Moved global states into classes.