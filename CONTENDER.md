# Worked portable SIMD block index

This branch integrates the best forward-order contender into rust-lapper's
normal API. `find()` and `seek()` retain their existing signatures and lazy,
ascending-start output. There is no workload score, data classifier, mode flag,
or alternate public query method.

The implementation keeps one canonical start-sorted interval vector and adds
fixed 32-entry block facts:

- a maximum end and a next-greater-block link for exact forward skips;
- a minimum end for an exact all-ends-pass route;
- a prefix maximum for finding the first possible block;
- ends in start order for mixed-block masks.

Mixed blocks use NEON on AArch64, AVX2 when available on x86-64, and an exact
scalar fallback elsewhere. The SIMD dispatch covers all primitive integer
coordinate types from 8 through 64 bits, including `usize` and `isize`.
Custom `PrimInt` implementations use the scalar mask.

`insert()`, `merge_overlaps()`, construction, and serde deserialization all run
the same derived-index rebuild. Signed coordinates, minimum-value queries,
`count()` at a type maximum, and depth spans crossing zero have dedicated tests.
Serde keeps the original six-field representation and reconstructs the new
private sidecars when reading it.

The hot loop uses `get_unchecked` only after a length invariant has proved the
private sidecars cover the candidate block. Returned intervals remain safely
indexed. The assembly audit found no reason to add handwritten assembly: the
AArch64 compiler output already contains paired loads, narrowing, and horizontal
reduction, while forced-Haswell output contains the expected AVX2 compares and
movemask instructions.

Fresh native five-library total medians:

| Case | SuperIntervals | Worked Lapper | COITrees | rust-bio IITree | rust-bio AVL |
|---|---:|---:|---:|---:|---:|
| `1-2` | 6.063 ms | **5.729 ms** | 10.415 ms | 13.311 ms | 42.062 ms |
| `7-3` | 68.881 ms | **66.385 ms** | 89.782 ms | 149.975 ms | 344.913 ms |
| `8-7` | **550.372 ms** | 588.766 ms | 834.219 ms | 1259.076 ms | 2140.363 ms |

Alternating direct trials confirmed the shape: Lapper won total time on 15/15
`1-2` pairs, was 1.76% faster on `7-3`, and trailed by about 9.5% total on
dense `8-7`. These are Apple M3 AArch64 measurements; AVX2 was compiled and
instruction-audited but not timed on native Intel or AMD hardware.

See [`PORTABLE_SIMD_INDEX.md`](PORTABLE_SIMD_INDEX.md) for the full design,
safety argument, compatibility notes, tests, measurements, and identity audit.
Use [`AARCH64_U32_HAND_TYPE.md`](AARCH64_U32_HAND_TYPE.md) on this branch as the
single hand-typing course from original scalar Lapper through the full portable
implementation. At its ARM checkpoints, the intentionally narrow answer code is
the diff from `contender/block-32` to `tutorial/aarch64-u32-hand-typed`.
