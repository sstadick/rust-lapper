# Always-on SIMD block-mask contender

This branch archives the first exact SIMD-mask finalist. It keeps the
32-interval block index and adds start-order end storage, minimum ends per
block, and `Lapper<u32, T>::find_block_mask`.

Every candidate block takes one exact route:

1. `max_end <= query.start`: jump over the all-miss block.
2. `min_end > query.start`: all active lanes pass the end test, so return the
   sorted-start prefix directly.
3. Otherwise: calculate one NEON overlap mask, save it in the iterator, and
   consume its low set bits in forward order across `next()` calls.

There is no score, sampling heuristic, data classifier, or query mode switch.
The scalar fallback is correct but was not performance-tested. The measured
SIMD path is specialized to `u32` on AArch64.

Recorded nine-repeat medians from one direct comparison session:

| Case | Build | Query | Total |
|---|---:|---:|---:|
| `1-2` | 3.231 ms | 4.100 ms | 7.295 ms |
| `7-3` | 26.764 ms | 55.286 ms | 82.012 ms |
| `8-7` | 36.926 ms | 603.311 ms | 640.183 ms |

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
