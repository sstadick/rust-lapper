# SIMD block mask with weighted lane reduction and paired loads

This branch builds on the first exact SIMD-mask finalist. It keeps the
32-interval block index, start-order end storage, block minima, and
`Lapper<u32, T>::find_block_mask`.

The weighted reduction remains unchanged: each four-lane NEON result is ANDed
with bit weights `[1, 2, 4, 8]` and reduced with `vaddvq_u32`. The follow-up
uses `vld1q_u32_x2` to load eight adjacent starts or stops into two registers
per iteration. This lets AArch64 select a multi-register `LD1` form while the
same two four-lane masks are assembled in forward order.

Every candidate block takes one exact route:

1. `max_end <= query.start`: jump over the all-miss block.
2. `min_end > query.start`: all active lanes pass the end test, so return the
   sorted-start prefix directly.
3. Otherwise: calculate one NEON overlap mask, save it in the iterator, and
   consume its low set bits in forward order across `next()` calls.

There is no score, sampling heuristic, data classifier, or query mode switch.
The scalar fallback is correct but was not performance-tested. The measured
SIMD path is specialized to `u32` on AArch64.

Fifteen paired trials against the weighted single-load version, alternating
execution order and using two unreported warmups per sample:

| Case | Single-load query | Paired-load query | Paired median change |
|---|---:|---:|---:|
| `1-2` | 3.595 ms | 3.443 ms | -4.6% |
| `7-3` | 51.008 ms | 49.410 ms | -3.1% |

A separate nine-repeat dense run measured 579.273 ms for the single-load
version and 579.313 ms for paired loads, which is neutral.

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
