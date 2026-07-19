# SIMD block mask with proof-scoped unchecked indexing

This branch builds on the first exact SIMD-mask finalist. It keeps the
32-interval block index, start-order end storage, block minima, and
`Lapper<u32, T>::find_block_mask`.

The loop uses `vld1q_u32_x2` to load eight adjacent starts or stops into two
registers per iteration. It narrows the two four-lane comparison vectors into
eight `u16` lanes, ANDs them with weights `[1, 2, 4, 8, 16, 32, 64, 128]`, and
uses one `vaddvq_u16` reduction to produce eight mask bits. On the measured
Apple M3 binary this compiles to paired `ldp` loads, one `uzp1.8h`, and one
`addv.8h` per eight lanes.

Every candidate block takes one exact route:

1. `max_end <= query.start`: jump over the all-miss block.
2. `min_end > query.start`: all active lanes pass the end test, so return the
   sorted-start prefix directly.
3. Otherwise: calculate one NEON overlap mask, save it in the iterator, and
   consume its low set bits in forward order across `next()` calls.

There is no score, sampling heuristic, data classifier, or query mode switch.
The scalar fallback is correct but was not performance-tested. The measured
SIMD path is specialized to `u32` on AArch64.

Fifteen paired trials against the paired-load, two-reduction version, alternating
execution order and using two unreported warmups per sample:

| Case | Two reductions | Eight-lane reduction | Paired median change | Faster pairs |
|---|---:|---:|---:|---:|
| `1-2` | 3.380 ms | 3.221 ms | -4.4% | 15/15 |
| `7-3` | 48.212 ms | 47.225 ms | -2.6% | 13/15 |
| `8-7` | 573.589 ms | 570.430 ms | -1.1% | 13/15 |

Both binaries were built with `-C target-cpu=native`. The eight-lane reduction
improved the paired median in all three cases.

This branch additionally removes redundant bounds checks for block metadata and
the two input slices after checking `block_start < intervals.len()`. Construction
keeps all sidecars at the same length, each block number is derived from a valid
32-entry boundary, and `block_end` is clamped to the interval length. Result
yields remain safely indexed.

The inlined AArch64 query function falls from 442 to 382 static instructions and
from nine bounds-panic edges to two. Fifteen alternating-order pairs against the
safe eight-lane branch measured:

| Case | Safe query | Unchecked-index query | Paired median change | Faster pairs |
|---|---:|---:|---:|---:|
| `1-2` | 3.251 ms | 3.008 ms | -7.7% | 15/15 |
| `7-3` | 47.459 ms | 43.504 ms | -8.1% | 15/15 |
| `8-7` | 573.038 ms | 566.117 ms | -1.0% | 10/15 |

In a fresh five-library run this version beat SuperIntervals in total time on
`1-2` (6.026 vs 6.129 ms) and `7-3` (66.783 vs 67.951 ms), and trailed on dense
`8-7` (579.622 vs 548.992 ms).

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
