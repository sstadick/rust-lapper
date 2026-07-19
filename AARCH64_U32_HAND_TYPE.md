# Hand-type the AArch64 `u32` mask path

This is a deliberately narrow exercise: add the forward-order AArch64 NEON mask
to the preserved 32-block rust-lapper prototype. The answer branch avoids generic
type dispatch, x86, mutation rebuilding, and unchecked indexing so each machine-
level idea is visible once.

## Set up the exercise

Create your own branch from the exact 32-block source:

```bash
git switch contender/block-32
git switch -c practice/aarch64-u32-mask
```

Use the tutorial branch only as the answer key:

```bash
git diff tutorial/aarch64-u32-hand-typed -- src/lib.rs tests/block_index.rs
```

Run this checkpoint after every stage:

```bash
cargo test --test block_index
```

The test checks the count, the exact interval references, and forward order
against brute force over randomized inputs.

## 1. Store ends in start order

The canonical interval vector and `starts` already share start order. Add one
sidecar that preserves each interval's stop in that same order:

```rust
stops_by_start: Vec<I>,
```

Build it before sorting the separate `stops` vector used by `count()`:

```rust
let (starts, mut stops): (Vec<_>, Vec<_>) =
    intervals.iter().map(|x| (x.start, x.stop)).unzip();
let stops_by_start = stops.clone();
```

Why: a SIMD lane must load the start and stop for the same interval. The sorted
`stops` array cannot provide that pairing.

## 2. Add exact block facts

For every 32-entry block, retain both its minimum and maximum stop. Also build a
prefix maximum over block maxima.

The three facts have separate jobs:

```text
max_end <= query.start  -> every lane misses; follow the skip link
min_end > query.start   -> every active lane's end passes
otherwise               -> the block is mixed; calculate a mask
prefix_max <= start     -> binary-search directly to the first possible block
```

These are proofs, not a workload score.

## 3. Give the iterator memory

Add state for one query-local answer:

```rust
mask_block_start: usize,
mask: u32,
dense_next: usize,
dense_end: usize,
```

At the top of `next()`, drain a dense prefix first, then a saved mask:

```rust
if self.mask != 0 {
    let lane = self.mask.trailing_zeros() as usize;
    self.mask &= self.mask - 1;
    return Some(&self.inner.intervals[self.mask_block_start + lane]);
}
```

`trailing_zeros()` selects the lowest lane. `mask & (mask - 1)` erases that lane.
Those two operations are the forward-order contract.

## 4. Make a scalar mask before using NEON

Write and test this version first:

```rust
fn overlap_mask(
    starts: &[u32],
    stops: &[u32],
    query_start: u32,
    query_stop: u32,
) -> u32 {
    let mut mask = 0;
    for lane in 0..starts.len() {
        if stops[lane] > query_start && starts[lane] < query_stop {
            mask |= 1 << lane;
        }
    }
    mask
}
```

Do not move on until the randomized test passes. SIMD must be only another way to
compute this exact integer.

## 5. Type one four-lane NEON comparison

On AArch64, one 128-bit register holds four `u32` lanes:

```rust
let lane_starts = vld1q_u32(starts.as_ptr().add(lane));
let lane_stops = vld1q_u32(stops.as_ptr().add(lane));

let overlapping = vandq_u32(
    vcgtq_u32(lane_stops, query_start_v),
    vcgtq_u32(query_stop_v, lane_starts),
);
```

Read the comparisons literally:

```text
interval.stop > query.start
query.stop > interval.start
```

Each true lane is `0xffff_ffff`; each false lane is zero.

## 6. Turn four truth lanes into four bits

Type the weighted reduction once:

```rust
let weights = vld1q_u32([1_u32, 2, 4, 8].as_ptr());
let bits = vaddvq_u32(vandq_u32(overlapping, weights));
mask |= bits << lane;
```

A true lane keeps its weight and a false lane keeps zero. Adding
`[1, 0, 4, 8]` gives `13`, whose binary representation is `1101`. The sum is the
mask; it is not a count.

## 7. Load and pack eight lanes

Replace two four-lane iterations with one paired load:

```rust
let lane_starts = vld1q_u32_x2(starts.as_ptr().add(lane));
let lane_stops = vld1q_u32_x2(stops.as_ptr().add(lane));
```

Compare `.0` and `.1` independently. Then narrow and join their all-ones/all-zero
answers:

```rust
let overlapping = vcombine_u16(
    vmovn_u32(overlapping0),
    vmovn_u32(overlapping1),
);
let weights = vld1q_u16(
    [1_u16, 2, 4, 8, 16, 32, 64, 128].as_ptr(),
);
let bits = vaddvq_u16(vandq_u16(overlapping, weights));
mask |= u32::from(bits) << lane;
```

Narrowing does not lose truth: the low 16 bits of all ones are still all ones.
This lets one `addv.8h` produce eight bits instead of two `addv.4s` reductions.

Keep the four-lane loop for a partial final block, and keep a scalar tail for a
length not divisible by four.

## 8. Inspect what the compiler made

Ask rustc to retain assembly for the release library:

```bash
RUSTFLAGS='-C target-cpu=native' \
  cargo rustc --release --lib -- --emit=asm
rg -n 'ldp|uzp1|addv' target/release/deps/rust_lapper-*.s
```

The measured Apple M3 path contains paired `ldp` vector loads, `uzp1.8h`, and
`addv.8h`. Handwritten assembly is unnecessary for this core.

## Stop here once

At this point you have typed every important ARM operation and can compare the
whole answer against `tutorial/aarch64-u32-hand-typed`. The production-oriented
`worked/portable-simd-index` branch deliberately moves type dispatch, mutation
rebuilding, x86 AVX2, scalar fallback, and audited `get_unchecked` calls into
separate helpers. Those abstractions are useful after the mechanics are familiar;
they are noise during this first pass.
