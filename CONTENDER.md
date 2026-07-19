# SIMD block mask with paired loads and eight-lane reduction

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

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
