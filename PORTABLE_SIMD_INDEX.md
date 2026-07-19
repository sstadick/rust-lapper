# Portable forward-order SIMD index

## Goal

This branch takes the fastest retained rust-lapper contender and makes it a
complete implementation rather than a benchmark-only `u32` experiment.

The constraints are deliberate:

- keep `Lapper::find()` and `Lapper::seek()` as the query API;
- yield the same interval references in ascending start order;
- use the same exact algorithm for every workload;
- rebuild every derived array after supported mutation;
- support signed and unsigned primitive coordinates from 8 through 64 bits;
- accelerate AArch64 and x86-64 without making other targets incorrect; and
- use `unsafe` only where a documented invariant removes measured hot-loop work.

CPU feature dispatch is not a query mode. A query always follows the same block
algorithm. Only the instruction sequence used to calculate an exact mixed-block
mask changes with the CPU and coordinate type.

## Public behavior

The production path is the ordinary one:

```rust
for interval in lapper.find(query_start, query_stop) {
    // Same borrowed Interval, same ascending-start order.
}

for interval in lapper.seek(query_start, query_stop, &mut cursor) {
    // Same cursor-facing API.
}
```

The experimental `find_block_mask()` method is gone. No mode setting or score is
exposed. `Lapper`, `Interval`, `IterFind`, mutation methods, iterators, and serde
field names remain available as before.

The generic implementation now requires `I: 'static` so it can use `TypeId` to
select only exact built-in integer representations before any SIMD pointer cast.
This is invisible for all primitive coordinate types. It is a theoretical source
compatibility constraint for an exotic lifetime-carrying custom `PrimInt`.

## Derived layout

The canonical `intervals` vector remains sorted by `(start, stop)`. The index does
not reorder that vector. Construction derives:

| Array | Length | Purpose |
|---|---:|---|
| `starts` | `n` | Start binary searches and SIMD input |
| `stops` | `n` | Globally sorted ends for `count()` |
| `stops_by_start` | `n` | End paired with each start for a mask lane |
| `block_max_ends` | `ceil(n / 32)` | Prove a whole block misses |
| `block_min_ends` | `ceil(n / 32)` | Prove every active lane's end passes |
| `block_index` | `ceil(n / 32)` | Next block to the right with a greater maximum end |
| `block_prefix_max_ends` | `ceil(n / 32)` | Binary-search the first possible block |

The block size remains 32. It maps one block result to one `u32` mask, amortizes
metadata, and was the best compromise across the measured sparse and dense cases.
Making SIMD lanes wider changes how quickly a mixed block is classified; it does
not remove the metadata, output, and partial-block tradeoffs that made 64-entry
blocks regress in an earlier trial.

`block_index` is built with a reverse monotonic stack. If block `b` has
`max_end <= query.start`, its link points to the first later block with a strictly
greater maximum. Every skipped intermediate block has a maximum no greater than
block `b`, so all of them are also exact misses.

## One always-on traversal

`find()` binary-searches `block_prefix_max_ends` to skip the leading prefix whose
intervals all end at or before `query.start`. `seek()` retains its cursor behavior,
then rounds the candidate position down to its containing 32-entry block. The
normal block loop has three exact routes:

```text
block.max_end <= query.start
    -> all 32 miss; follow the next-greater block link

block.min_end > query.start
    -> every end passes; return starts < query.stop as a dense prefix

otherwise
    -> starts and ends are mixed; calculate and save an exact overlap mask
```

Before any route, the first block start is compared with `query.stop`. Starts are
globally sorted, so failure ends the iterator. In the dense route a partition
point finds the active start prefix. In the mixed route the mask checks both
half-open overlap predicates:

```text
interval.stop > query.start
interval.start < query.stop
```

This removes the old active-prefix special case from the mixed path. SIMD checks
the candidate slice that exists, including start failures; those lanes simply
become zero.

## Why order is unchanged

The iterator stores one mixed-block mask. Each `next()` call removes its least
significant set bit:

```rust
let lane = self.mask.trailing_zeros() as usize;
self.mask &= self.mask - 1;
return Some(&self.inner.intervals[self.mask_block_start + lane]);
```

The lowest bit is the lowest position in the start-sorted block. Blocks are only
visited to the right, dense prefixes are walked from their first item, and skip
links only point right. Therefore results remain references into the canonical
vector in ascending start order. SIMD changes how a block answer is calculated,
not the observable traversal identity.

## SIMD backends and integer widths

`src/simd.rs` owns one private operation:

```rust
overlap_mask(starts, stops, query_start, query_stop) -> u32
```

It always returns the same lane-to-bit mapping. Exact `TypeId` checks allow SIMD
only for these built-ins:

```text
u8 i8 u16 i16 u32 i32 u64 i64 usize isize
```

Any other `PrimInt` takes the scalar implementation. Partial final vectors also
use a small scalar tail.

### AArch64 NEON

NEON is part of baseline AArch64. Its registers are 128 bits, so one register has
16 byte lanes, 8 halfword lanes, 4 word lanes, or 2 doubleword lanes. Paired
loads process two adjacent registers when possible.

| Coordinate | Values per paired load | Truth-to-bit packing |
|---|---:|---|
| `u8` / `i8` | 32 | two 16-lane weighted byte reductions |
| `u16` / `i16` | 16 | narrow to bytes, then weighted reduction |
| `u32` / `i32` | 8 | narrow to halfwords, then one `addv.8h` |
| `u64` / `i64` | 4 | narrow to words, then one weighted reduction |

The central `u32` trick is:

```rust
let eight_truths = vcombine_u16(
    vmovn_u32(overlaps_low),
    vmovn_u32(overlaps_high),
);
let weights = vld1q_u16([1, 2, 4, 8, 16, 32, 64, 128].as_ptr());
let bits = vaddvq_u16(vandq_u16(eight_truths, weights));
```

A comparison produces all ones for true and zero for false. Narrowing preserves
that truth. AND keeps a lane's power-of-two weight only when it is true; the
horizontal sum is therefore the bit mask, not a match count.

### x86-64 AVX2

Every iterator performs the standard runtime AVX2 feature check once. AVX2 uses
256-bit registers and native movemask instructions:

| Coordinate | Values per vector | Mask extraction |
|---|---:|---|
| `u8` / `i8` | 32 | `vpmovmskb` |
| `u16` / `i16` | 16 | `vpmovmskb`, then compact duplicate byte bits |
| `u32` / `i32` | 8 | `vmovmskps` on comparison bits |
| `u64` / `i64` | 4 | `vmovmskpd` on comparison bits |

AVX2 integer greater-than comparisons are signed. Unsigned inputs are XORed with
their type's sign bit before comparing. That order-preserving transform maps the
unsigned domain onto the signed domain without changing the overlap predicate.

### Other CPUs and custom integers

The scalar backend calculates the same `u32` mask lane by lane. This keeps the
algorithm and output identical on targets such as `wasm32-wasip1` and on x86-64
machines without AVX2. It is a correctness fallback, not a workload-dependent
choice.

SVE/SVE2 and AVX-512 are possible later backends, but are not present here. Their
predicate/mask facilities could classify a 32-entry block in fewer vector groups.
They need native hardware measurements before becoming a recommendation.

## Mutation, construction, and serde

All sidecars are derived state. `rebuild_derived()` is the single construction
path and is called by:

- `Lapper::new()`;
- `insert()` after its sorted insertion;
- `merge_overlaps()` after it replaces the interval vector; and
- serde deserialization through `Lapper::new()`.

The rebuild uses one bulk `unzip()` for starts and start-order stops, clones the
stops once for global sorting, then builds block summaries. Restoring this bulk
path recovered build time that an earlier incremental prototype had lost.

With `with_serde`, serialization deliberately writes the original six fields:
`intervals`, `starts`, `stops`, `max_len`, `cov`, and `overlaps_merged`. New
private sidecars are not serialized. Deserialization treats `intervals` as the
source of truth and rebuilds every sidecar, then restores the cached coverage and
merge flag. Existing serialized shape therefore does not acquire architecture-
or implementation-specific metadata.

## Signed-coordinate corrections

Generalizing beyond unsigned coordinates exposed three assumptions in the old
code:

- `seek()` now uses `checked_sub(max_len)` and falls back to `I::min_value()`;
- `count()` uses `partition_point(stop <= start)` instead of forming `start + 1`,
  which can overflow at the coordinate maximum; and
- `depth()` has an explicit initialization flag instead of using zero as a
  sentinel, so an interval spanning negative through positive coordinates begins
  at its actual start.

These are scalar correctness repairs needed by the wider type support; they are
not query heuristics.

## Bounds-check and unsafe audit

Release assembly showed that the SIMD helpers already use pointer vector loads;
there were no hidden per-lane Rust slice checks in the vector core. The remaining
hot-loop checks were on private sidecar indexing.

The accepted branch uses `get_unchecked` for those sidecars after one public
invariant boundary:

1. `block_start < starts.len()` is checked.
2. Constructors, supported mutation, and deserialization rebuild `starts` and
   `stops_by_start` to equal lengths.
3. They also build one entry in every block metadata vector for each 32-entry
   start block.
4. `block_end` is clamped to `starts.len()`.
5. A link is either a valid later block or the block-count sentinel.

The block link is guarded by a debug assertion. Returned `intervals` remain
safely indexed, so unsupported direct edits to the public vector can at worst
produce stale behavior or a panic, not turn a result yield into unchecked memory
access. Directly appending to `intervals` is ignored by the private derived bound;
supported mutation must continue to use `insert()` or `merge_overlaps()`.

For the specialized AArch64 query, this change reduced the inlined function from
442 to 382 static instructions and bounds-panic edges from nine to two. Making
the final result access unchecked reduced it further to 363/zero but improved
queries by only about 0.3-0.6%, so that broader unsafe surface was rejected.

Handwritten assembly was also rejected. LLVM already emits the intended AArch64
`ldp`, `uzp1.8h`, and `addv.8h` sequence. A forced Haswell build contains
`vpmovmskb`, `vmovmskps`, `vmovmskpd`, `vpcmpgt*`, `vpxor`, and `vpand`. Inline
assembly would add review and register-allocation risk without removing a known
instruction in these cores.

The detailed audit is in the bakeoff's `ASSEMBLY_BOUNDS_AUDIT.md`.

## Correctness and portability checks

The branch passes:

```text
native cargo test --all-features
native cargo clippy --all-targets --all-features -- -D warnings
x86_64-apple-darwin cargo test --all-features --test portable_index
x86_64-apple-darwin cargo check --all-features
wasm32-wasip1 cargo clippy --lib --all-features -- -D warnings
```

Coverage includes:

- randomized `find()`, `seek()`, and `count()` against brute force;
- exact reference order for all ten primitive integer types;
- negative intervals and a query at the signed minimum;
- insert and merge rebuilding across more than four blocks;
- serde round trips with derived metadata reconstruction;
- signed depth across zero;
- safety when the public interval vector is directly extended;
- exhaustive AVX2 `u16` movemask compaction; and
- forced calls to each AVX2 primitive function when AVX2 is available.

The x86-64 test executable ran under Rosetta on the Apple host. Rosetta did not
advertise AVX2, so the backend's instructions were verified in forced-Haswell
assembly rather than performance-tested or executed there. Native Intel/AMD CI
remains required before merging an x86 performance claim.

## Performance record

All times below are medians in milliseconds. The five-library worked run used
the normal `find()` API and native Apple M3 compilation.

| Case | SuperIntervals | Worked Lapper | COITrees | rust-bio IITree | rust-bio AVL |
|---|---:|---:|---:|---:|---:|
| `1-2` | 6.063 | **5.729** | 10.415 | 13.311 | 42.062 |
| `7-3` | 68.881 | **66.385** | 89.782 | 149.975 | 344.913 |
| `8-7` | **550.372** | 588.766 | 834.219 | 1259.076 | 2140.363 |

Those totals combine independently reported build and query medians. Direct
alternating Lapper/SuperIntervals process pairs give the more reliable close-call
interpretation:

| Case | Worked Lapper versus SuperIntervals total | Pair wins |
|---|---:|---:|
| `1-2` | -7.71% paired median | 15 / 15 |
| `7-3` | -1.76% paired median | 12 / 15 |
| `8-7` | +9.50% paired median | 0 / 10 |

On `7-3`, Lapper's query remained about 13% slower but its cheaper build erased
that difference for one build plus one query batch. On dense `8-7`, output work
dominates and SuperIntervals retains a clear advantage.

The worked generic integration stayed within about 1.3% of the specialized
unchecked `u32` branch in paired query trials:

| Case | Worked query change versus specialized | Worked wins |
|---|---:|---:|
| `1-2` | -2.14% | 11 / 15 |
| `7-3` | +1.53% | 1 / 15 |
| `8-7` | +1.31% | 2 / 10 |

For orientation, the original v1.3.0 query medians were 5.648, 5049.867, and
884.546 ms on `1-2`, `7-3`, and `8-7`. Those original-versus-worked values came
from separate controlled runs and should not be read as paired microbenchmark
percentages. Raw files and methods live under the bakeoff's
`results/forward-simd-2026-07-18/` directory.

## Why this is not AIList or another renamed index

This conclusion is structural, not based on names such as "maximum end."

| Implementation | Canonical layout | Index shape | Query direction / state |
|---|---|---|---|
| This Lapper | One vector sorted by start | Fixed positional block summaries and forward next-greater links | Forward blocks; query-local low-bit-first mask |
| AIList | Intervals decomposed into multiple containment-derived lists | Per-list running maximum ends | Searches/scans each component backward |
| NCList | Intervals organized by containment | Nested containment hierarchy | Hierarchical sublist traversal |
| SuperIntervals | Separate starts, ends, payloads, and branch arrays | Previous-greater/equal end links | Begins near the query boundary and follows links backward |
| COITrees | Augmented tree nodes in cache-oriented order | Centered interval tree | Tree traversal over augmented nodes |
| rust-bio IITree / AVL | Tree-owned interval nodes | Interval-tree augmentation | Tree traversal |

This branch does not perform AIList's containment decomposition, create component
lists, or scan those lists backward. It does not build NCList's containment
hierarchy. It does not use SuperIntervals' previous-greater link direction or
separate canonical storage. It does not reorder the intervals into a tree.

What is shared is general interval-index vocabulary: sorting endpoints, storing
end extrema, skipping regions proved unable to overlap, and using SIMD to evaluate
independent predicates. Those shared techniques do not make the data structures
the same. This comparison is not a patent-clearance or exhaustive prior-art
opinion.

Primary references used for the structural check:

- [AIList paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6901075/)
- [AIList C construction and query](https://github.com/databio/AIList/blob/366fe0a6f5c77cdd165c1452f49adf0ecbbef4e6/src/AIList.c#L130-L317)
- [gtars Rust AIList](https://github.com/databio/gtars/blob/cd75994966f039aa52d9359d64a3a201eee25de7/gtars-overlaprs/src/ailist.rs#L92-L375)
- [SuperIntervals 0.3.6 source](https://github.com/biodatageeks/sequila-native/blob/745d40f77da7ced5d540f9285eb5123ba12682ff/sequila/sequila-core/superintervals/src/superintervals.rs)
- [COITrees 0.4.0 non-SIMD source](https://github.com/dcjones/coitrees/blob/4afeeebe6e5d3229c2c191d5676ed71310d953e0/src/nosimd.rs)

## Branch map

| Branch | Commit | Purpose |
|---|---|---|
| `baseline/v1.3.0` | `c545d40` | Original rust-lapper release state |
| `contender/per-interval` | `0af2d23` | Reconstructed next-greater index per interval |
| `contender/block-32` | `559c183` | Preserved always-on 32-entry block index |
| `contender/simd-mask-narrow8-ld1x2` | `84872e7` | Safe specialized AArch64 mask |
| `contender/simd-mask-narrow8-unchecked-index` | `b6ba527` | Specialized best with audited unchecked sidecars |
| `tutorial/aarch64-u32-hand-typed` | `dfc5f01` | Narrow answer key for typing the ARM path once |
| `reference/aarch64-u32-hand-typed` | `b6ba527` | Exact full-reference comparison point |
| `worked/portable-simd-index` | this branch | Mutation-safe, multi-type, multi-architecture integration |

The tutorial starts from `contender/block-32` and intentionally excludes generic
dispatch, AVX2, mutation rebuilding, and unchecked indexing. That keeps the first
hand-typed ARM pass readable. The worked branch is the complete candidate, not the
tutorial branch.

## Remaining decisions

- Measure AVX2 on native Intel and AMD hardware.
- Add native big-endian and non-Apple AArch64 CI if those are support targets.
- Measure memory/cache behavior across coordinate and payload widths.
- Decide whether the theoretical new `I: 'static` constraint is acceptable.
- Treat SVE2 and AVX-512 as separate CPU backends only after native evidence.
- Preserve safe result indexing unless a new measurement justifies expanding the
  unsafe proof.
