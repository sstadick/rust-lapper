# SIMD block mask with weighted lane reduction

This branch builds on the first exact SIMD-mask finalist. It keeps the
32-interval block index, start-order end storage, block minima, and
`Lapper<u32, T>::find_block_mask`.

The new change is deliberately small: each four-lane NEON result is ANDed with
bit weights `[1, 2, 4, 8]` and reduced with `vaddvq_u32`. This replaces four
lane extractions plus their scalar shifts and ORs with one horizontal sum and
one scalar shift/OR.

Every candidate block takes one exact route:

1. `max_end <= query.start`: jump over the all-miss block.
2. `min_end > query.start`: all active lanes pass the end test, so return the
   sorted-start prefix directly.
3. Otherwise: calculate one NEON overlap mask, save it in the iterator, and
   consume its low set bits in forward order across `next()` calls.

There is no score, sampling heuristic, data classifier, or query mode switch.
The scalar fallback is correct but was not performance-tested. The measured
SIMD path is specialized to `u32` on AArch64.

Seven-repeat same-session means against the original mask prototype:

| Case | Original mask query | Weighted query | Change |
|---|---:|---:|---:|
| `1-2` | 4.126 ms | 3.551 ms | -13.9% |
| `7-3` | 55.261 ms | 49.575 ms | -10.3% |
| `8-7` | 604.449 ms | 567.982 ms | -6.0% |

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
