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
| `1-2` | 3.614 ms | 3.482 ms | -3.2% | 12/15 |
| `7-3` | 49.464 ms | 48.376 ms | -2.1% | 14/15 |
| `8-7` | 586.479 ms | 586.697 ms | -0.02% | 8/15 |

The dense result is effectively neutral rather than evidence of a speedup.

For 1,956,864 intervals, the extra `u32` start-order ends cost about 7.47 MiB
and block minima cost about 0.23 MiB. Metadata must be rebuilt after mutation
before production integration.
