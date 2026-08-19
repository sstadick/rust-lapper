# Hand-type the whole forward index: original Lapper to portable SIMD

This is the single end-to-end exercise for understanding the worked
rust-lapper contender. It starts with the original start-sorted linear iterator,
adds the first per-interval skip index, replaces it with the 32-entry block
index, adds exact block routes and a saved scalar mask, types the AArch64 NEON
path one operation at a time, and finishes with the portable, mutation-safe
implementation.

Do not begin with the SIMD intrinsics. SIMD is only a faster implementation of
one integer-valued function. The data structure, skip proofs, iterator state,
and output order must already make sense before that function is replaced.

Every checkpoint has the same four gates:

1. **Type it:** the smallest useful source change.
2. **State the invariant:** the fact that makes the change correct.
3. **Work an example:** calculate the result without running Rust.
4. **Pass the gate:** compare with a preserved branch and run an exact-order test.

## Set up a practice worktree

Keep this worked branch open as the guide and create a second worktree from the
original release state:

```bash
cd /Users/sethstadick/dev/super_intervals_more/bakeoff/sandboxes/rust-lapper-contenders

git worktree add \
  -b practice/hand-typed-forward-index \
  ../rust-lapper-hand-typed \
  baseline/v1.3.0

cd ../rust-lapper-hand-typed
```

The preserved answer checkpoints are:

| Checkpoint | Branch | What it isolates |
|---|---|---|
| Original | `baseline/v1.3.0` | Start-sorted scalar Lapper |
| First index | `contender/per-interval` | One next-greater-end link per interval |
| Block index | `contender/block-32` | One summary and link per 32 intervals |
| First mask | `contender/simd-mask` | Exact routes, saved mask, first NEON implementation |
| Weighted mask | `contender/simd-mask-addv` | Four truth lanes reduced to four bits |
| Paired loads | `contender/simd-mask-addv-ld1x2` | Eight values loaded per loop |
| Eight-lane pack | `contender/simd-mask-narrow8-ld1x2` | Eight truths reduced once |
| Bounds audit | `contender/simd-mask-narrow8-unchecked-index` | Proof-scoped unchecked sidecars |
| Narrow ARM answer | `tutorial/aarch64-u32-hand-typed` | The ARM-only portion without portable integration |
| Full answer | `worked/portable-simd-index` | Normal API, all primitive types, mutation, serde, and fallbacks |

Use this as the progress sheet:

- [ ] 0. Install the forward brute-force oracle.
- [ ] 1. Explain original `max_len` search and its pathological case.
- [ ] 2. Build and prove per-interval next-greater links.
- [ ] 3. Replace them with 32-entry block summaries and links.
- [ ] 4. Binary-search block prefix maxima.
- [ ] 5. Derive the exact skip, dense, and mixed routes.
- [ ] 6. Save and drain a scalar block mask.
- [ ] 7. Reproduce that mask with four-lane NEON comparisons.
- [ ] 8. Pack four truth lanes with weighted horizontal addition.
- [ ] 9. Pair loads, narrow truth, and pack eight lanes once.
- [ ] 10. Find the intended instructions in release assembly.
- [ ] 11. Move the mask state behind normal `find()` and `seek()`.
- [ ] 12. Add exact type and CPU dispatch with scalar fallback.
- [ ] 13. Rebuild all derived state after mutation and deserialize.
- [ ] 14. Repair the scalar algorithms for signed coordinates.
- [ ] 15. Remove only bounds checks covered by private invariants.

At any point, inspect the historical delta rather than copying the whole file:

```bash
git diff baseline/v1.3.0..contender/per-interval -- src/lib.rs
git diff contender/per-interval..contender/block-32 -- src/lib.rs
git diff contender/block-32..contender/simd-mask -- src/lib.rs tests/block_index.rs
```

The historical index and SIMD branches through checkpoint 12 were immutable
query experiments: their older mutation methods do not rebuild every new
sidecar. Use `new()` plus queries while learning those stages. Checkpoint 13 is
where the exercise becomes mutation-safe.

## The contract that never changes

Intervals and queries are half-open:

```text
interval = [interval.start, interval.stop)
query    = [query.start, query.stop)
```

They overlap only when both strict comparisons pass:

```rust
interval.start < query_stop && interval.stop > query_start
```

Equality means no overlap. An interval ending at `query_start` is to the left;
an interval starting at `query_stop` is to the right.

The observable identity of this design is also fixed:

- `intervals` is the one canonical vector;
- it is sorted by `(start, stop)` during construction;
- queries yield borrowed entries from that vector;
- results appear in increasing vector position, hence start order; and
- no index may omit a true overlap or invent a false one.

The index may skip, classify, or compare many entries together. It may not
reorder the canonical vector or return results backward.

## Checkpoint 0: make brute force the executable specification

Before adding an index, add an integration test that calculates the answer in the
most obvious way and compares values in order. This is more important than a
count-only test: two iterators can return the same count while returning the
wrong entries or order.

Create `tests/block_index.rs`:

```rust
use rust_lapper::{Interval, Lapper};

#[test]
fn block_index_matches_forward_brute_force() {
    let mut state = 0x1234_5678_u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (state >> 32) as u32
    };

    let mut intervals = vec![Interval {
        start: 0,
        stop: 1_000_000,
        val: 0usize,
    }];
    for value in 1..5000 {
        let start = next() % 1_000_000;
        let len = 1 + next() % 1000;
        intervals.push(Interval {
            start,
            stop: start + len,
            val: value,
        });
    }
    let lapper = Lapper::new(intervals);

    for _ in 0..20_000 {
        let start = next() % 1_000_000;
        let stop = start + 1 + next() % 2000;
        let got: Vec<_> = lapper.find(start, stop).map(|iv| iv.val).collect();
        let expected: Vec<_> = lapper
            .intervals
            .iter()
            .filter(|iv| iv.start < stop && iv.stop > start)
            .map(|iv| iv.val)
            .collect();
        assert_eq!(got, expected, "query {start}..{stop}");
    }
}
```

### Invariant

The expected vector is produced by filtering the canonical start-sorted vector
without skipping or reordering. It is the oracle for every later iterator.

### Work it by hand

For intervals `[0, 3)`, `[2, 8)`, `[5, 7)`, and `[9, 12)`, query `[3, 6)`:

```text
[0, 3):  0 < 6, but 3 > 3 is false
[2, 8):  2 < 6 and 8 > 3 -> hit
[5, 7):  5 < 6 and 7 > 3 -> hit
[9, 12): 9 < 6 is false; all later starts also fail

answer: the second interval, then the third interval
```

### Pass the gate

```bash
cargo test --test block_index
```

Do not proceed until this passes on `baseline/v1.3.0`.

## Checkpoint 1: understand the original search

Original Lapper stores:

```text
intervals  canonical intervals sorted by (start, stop)
starts     starts sorted globally
stops      stops sorted globally, no longer paired with starts
max_len    length of the longest interval
```

`find(query_start, query_stop)` computes a conservative earliest possible start:

```rust
let earliest_start = query_start
    .checked_sub(&self.max_len)
    .unwrap_or_else(zero::<I>);
let off = Self::lower_bound(earliest_start, &self.intervals);
```

Why it is safe: no interval starting before `query_start - max_len` can be long
enough to reach `query_start`. From `off`, the iterator checks intervals one by
one and stops forever when `interval.start >= query_stop`.

### Invariant

Sorted starts prove the right-hand stopping condition. `max_len` proves the
left-hand starting condition.

### Work it by hand

If `query_start = 1_000` and `max_len = 40`, no interval starting before `960`
can reach the query. But one extremely long interval makes `max_len` enormous,
pushes `earliest_start` back toward zero, and forces nearly every query to scan a
large prefix. That is the pathological case the new index must repair.

### Pass the gate

Be able to answer these before continuing:

```text
Why can the iterator stop on interval.start >= query_stop?
Why can it not stop merely because one interval.stop <= query_start?
Why does one global long interval make max_len a weak lower bound?
```

## Checkpoint 2: add one next-greater link per interval

The first index asks: if interval `i` ends too early, which later interval is the
first one whose end is strictly greater?

Add one sidecar:

```rust
jump_index: Vec<usize>,
```

Build it in linear time with a monotonic stack, scanning right to left:

```rust
let interval_count = intervals.len();
let mut jump_index = vec![interval_count; interval_count];
let mut stack = Vec::<usize>::new();

for i in (0..interval_count).rev() {
    while stack
        .last()
        .is_some_and(|j| intervals[*j].stop <= intervals[i].stop)
    {
        stack.pop();
    }
    jump_index[i] = stack.last().copied().unwrap_or(interval_count);
    stack.push(i);
}
```

Why the stack gives the *nearest* greater end:

- before processing `i`, the stack contains undominated positions to its right;
- viewed from the top outward, their ends are strictly increasing;
- anything popped has an end no greater than `end[i]` and is farther right, so
  `i` is both closer and at least as capable of reaching a future query;
- after those dominated entries are removed, the top is the closest surviving
  position with an end strictly greater than `end[i]`.

That domination argument is why construction is linear: every position is
pushed once and popped at most once.

Then replace the scalar miss step with the link:

```rust
while self.off < self.inner.intervals.len() {
    let interval = &self.inner.intervals[self.off];

    if interval.start >= self.stop {
        break;
    }
    if interval.stop > self.start {
        self.off += 1;
        return Some(interval);
    }

    self.off = self.inner.jump_index[self.off];
}
```

### Invariant

Suppose `jump_index[i] = j`. Every index strictly between `i` and `j` has
`stop <= intervals[i].stop`. The link is used only after proving
`intervals[i].stop <= query_start`, so every skipped interval also ends at or
before the query and cannot overlap.

The link points right, so forward order is preserved.

### Work it by hand

For start-ordered end values:

```text
index: 0  1  2  3  4
end:   7  3  5  2  9
link:  4  2  4  4  sentinel
```

With `query_start = 6`, a miss at index 1 follows `1 -> 2 -> 4`:

```text
end[1] = 3 <= 6
end[2] = 5 <= 6
end[4] = 9 > 6
```

Index 3 is skipped by `2 -> 4` because its end `2` is no greater than end `5`.

### Pass the gate

```bash
cargo test --test block_index
git diff contender/per-interval -- src/lib.rs
```

This historical branch is an immutable-query experiment. Its `insert()` and
`merge_overlaps()` paths do not rebuild the new link array. Do not treat that as
production-correct; mutation is repaired in checkpoint 13.

## Checkpoint 3: amortize the index over fixed blocks

The per-interval index costs one `usize` per interval and performs irregular link
work for individual misses. Replace it with one summary for each 32 consecutive
intervals:

```rust
const INDEX_BLOCK_SIZE: usize = 32;

block_index: Vec<usize>,
block_max_ends: Vec<I>,
block_prefix_max_ends: Vec<I>,
```

First calculate each block's maximum end:

```rust
let mut block_max_ends = Vec::new();
for block in intervals.chunks(INDEX_BLOCK_SIZE) {
    let mut max_end = block[0].stop;
    for interval in &block[1..] {
        max_end = std::cmp::max(max_end, interval.stop);
    }
    block_max_ends.push(max_end);
}
```

Run the same next-greater construction over block maxima:

```rust
let block_count = block_max_ends.len();
let mut block_index = vec![block_count; block_count];
let mut stack = Vec::<usize>::new();

for block in (0..block_count).rev() {
    while stack
        .last()
        .is_some_and(|next| block_max_ends[*next] <= block_max_ends[block])
    {
        stack.pop();
    }
    block_index[block] = stack.last().copied().unwrap_or(block_count);
    stack.push(block);
}
```

At a block boundary, skip only when the entire block is an exact miss:

```rust
let block = block_start / INDEX_BLOCK_SIZE;
if self.inner.block_max_ends[block] <= self.start {
    self.next_block_start =
        self.inner.block_index[block] * INDEX_BLOCK_SIZE;
    continue;
}
```

### Invariant

`block_max_end <= query_start` proves every end in the block fails the strict
`end > query_start` overlap test. If the link jumps from block `b` to block `j`,
all intermediate block maxima are no greater than `max_end[b]`, so they also
fail.

### Work it by hand

Use a paper-only block size of four. The real implementation remains 32.

| Block | Start-ordered intervals | Minimum end | Maximum end |
|---|---|---:|---:|
| 0 | `[0,4) [2,30) [5,7) [8,9)` | 4 | 30 |
| 1 | `[10,12) [13,14) [15,18) [19,21)` | 12 | 21 |
| 2 | `[22,60) [24,25) [27,35) [40,44)` | 25 | 60 |

The block maxima are `[30, 21, 60]`, so the links are `[2, 2, sentinel]`.
For query `[23, 28)`, block 1 has `max_end = 21 <= 23` and can jump directly to
block 2.

### Pass the gate

```bash
cargo test --test block_index
git diff contender/block-32 -- src/lib.rs
```

Explain why a block maximum can prove an all-miss block but cannot prove which
individual intervals hit inside a candidate block.

## Checkpoint 4: add the prefix maximum entry search

`find()` should not enter block zero for every query. Build a monotonically
nondecreasing prefix over block maxima:

```rust
let mut block_prefix_max_ends = block_max_ends.clone();
for block in 1..block_prefix_max_ends.len() {
    block_prefix_max_ends[block] = std::cmp::max(
        block_prefix_max_ends[block - 1],
        block_prefix_max_ends[block],
    );
}
```

Binary-search it at query creation:

```rust
let first_block = self
    .block_prefix_max_ends
    .partition_point(|max_end| *max_end <= query_start);
let next_block_start = first_block * INDEX_BLOCK_SIZE;
```

### Invariant

If `prefix_max[b] <= query_start`, every interval in every block through `b`
ends too early. The first prefix value greater than `query_start` is therefore
the first block that could contain an overlap.

### Work it by hand

For block maxima `[30, 21, 60]`, the prefix maxima are `[30, 30, 60]`.

```text
query_start = 23 -> partition point 0 -> block 0 may contain [2,30)
query_start = 31 -> partition point 2 -> blocks 0 and 1 are impossible
query_start = 61 -> partition point 3 -> no block can overlap
```

### Pass the gate

Add assertions for those three partition points or write them on paper before
running:

```bash
cargo test --test block_index
```

Notice that final `find()` no longer depends on global `max_len`; `seek()` still
uses it to advance a cursor for sorted query streams.

## Checkpoint 5: add exact candidate-block routes

Maximum ends prove an all-miss route. Add the complementary minimum end:

```rust
block_min_ends: Vec<I>,
```

Calculate minimum and maximum together. Also preserve ends in start order:

```rust
let (starts, mut stops): (Vec<_>, Vec<_>) =
    intervals.iter().map(|iv| (iv.start, iv.stop)).unzip();
let stops_by_start = stops.clone();
stops.sort();
```

These two end arrays are intentionally different:

```text
stops           globally sorted; used by count()
stops_by_start  lane i belongs to starts[i]; used by block masks
```

Every candidate block now has three exact routes:

```text
max_end <= query_start
    -> all miss; follow the forward block link

min_end > query_start
    -> every end passes; return only starts < query_stop

otherwise
    -> end results are mixed; calculate every lane exactly
```

For the all-ends-pass route, sorted starts identify one dense prefix:

```rust
let starts = &self.inner.starts[block_start..block_end];
let active_len = starts.partition_point(|start| *start < self.stop);
self.dense_next = block_start;
self.dense_end = block_start + active_len;
```

### Invariant

`min_end > query_start` proves the end half of overlap for every lane. It does
not prove the start half. The active prefix is still needed because starts at or
after `query_stop` do not overlap.

This answers the active-prefix question precisely:

- dense route: yes, find the active start prefix;
- mixed route: no separate prefix is needed once the mask checks both predicates;
- whole query: stop forever when the first start of a block is at or beyond
  `query_stop`.

### Work it by hand

Return to query `[23, 28)` in the three toy blocks:

```text
block 0: min 4 <= 23 < max 30 -> mixed
block 1: max 21 <= 23         -> jump
block 2: min 25 > 23          -> every end passes
```

Block 2 starts are `[22, 24, 27, 40]`. Its active prefix is the first three
lanes because `22`, `24`, and `27` are less than `28`; `40` is not.

### Pass the gate

Before adding SIMD, write down why both inequalities are strict:

```text
max_end == query_start -> all miss
min_end == query_start -> not all ends pass
start == query_stop    -> lane misses
```

## Checkpoint 6: save one scalar block mask

Add a temporary `find_block_mask()` iterator while retaining normal `find()` as
the control. This is scaffolding for the experiment, not a permanent API or a
mode switch.

Give the iterator memory:

```rust
next_block_start: usize,
mask_block_start: usize,
mask: u32,
dense_next: usize,
dense_end: usize,
start: u32,
stop: u32,
```

Start with a scalar mask that is obviously equivalent to brute force:

```rust
fn overlap_mask(
    starts: &[u32],
    stops: &[u32],
    query_start: u32,
    query_stop: u32,
) -> u32 {
    let mut mask = 0_u32;
    for lane in 0..starts.len() {
        if stops[lane] > query_start && starts[lane] < query_stop {
            mask |= 1 << lane;
        }
    }
    mask
}
```

Drain the lowest saved bit across later `next()` calls:

```rust
if self.mask != 0 {
    let lane = self.mask.trailing_zeros() as usize;
    self.mask &= self.mask - 1;
    return Some(&self.inner.intervals[self.mask_block_start + lane]);
}
```

The full safe state machine is:

```rust
loop {
    if self.dense_next < self.dense_end {
        let index = self.dense_next;
        self.dense_next += 1;
        return Some(&self.inner.intervals[index]);
    }

    if self.mask != 0 {
        let lane = self.mask.trailing_zeros() as usize;
        self.mask &= self.mask - 1;
        return Some(&self.inner.intervals[self.mask_block_start + lane]);
    }

    let block_start = self.next_block_start;
    if block_start >= self.inner.starts.len()
        || self.inner.starts[block_start] >= self.stop
    {
        return None;
    }

    let block = block_start / INDEX_BLOCK_SIZE;
    if self.inner.block_max_ends[block] <= self.start {
        self.next_block_start =
            self.inner.block_index[block] * INDEX_BLOCK_SIZE;
        continue;
    }

    let block_end =
        (block_start + INDEX_BLOCK_SIZE).min(self.inner.starts.len());

    if self.inner.block_min_ends[block] > self.start {
        let starts = &self.inner.starts[block_start..block_end];
        let active_len =
            starts.partition_point(|lane_start| *lane_start < self.stop);
        self.dense_next = block_start;
        self.dense_end = block_start + active_len;
        self.next_block_start = if active_len == starts.len() {
            block_end
        } else {
            self.inner.starts.len()
        };
        continue;
    }

    self.mask_block_start = block_start;
    self.next_block_start = block_end;
    self.mask = overlap_mask(
        &self.inner.starts[block_start..block_end],
        &self.inner.stops_by_start[block_start..block_end],
        self.start,
        self.stop,
    );
}
```

### Invariant

Bit `k` describes lane `k` in one block. `trailing_zeros()` selects the lowest
remaining lane, and `mask & (mask - 1)` erases exactly that bit. Blocks only move
right. Therefore output remains in canonical vector order.

### Work it by hand

For toy block 0 and query `[23, 28)`:

```text
lane 0 [0,4)  -> 0
lane 1 [2,30) -> 1
lane 2 [5,7)  -> 0
lane 3 [8,9)  -> 0

mask = 0b0010
trailing_zeros(mask) = 1
next() returns canonical index block_start + 1
```

Across all blocks the result indices are `[1, 8, 9, 10]`, still in start order.

### Pass the gate

Temporarily extend the oracle test:

```rust
let got_mask: Vec<_> = lapper
    .find_block_mask(start, stop)
    .map(|iv| iv.val)
    .collect();
assert_eq!(got_mask, expected);
```

Then run:

```bash
cargo test --test block_index
```

Do not type an intrinsic until the scalar mask passes.

## Checkpoint 7: compute four `u32` lanes with NEON

On AArch64, one 128-bit NEON register holds four `u32` values. Inside an
`unsafe` block, broadcast the query bounds once, then load and compare four
paired intervals:

```rust
use std::arch::aarch64::*;

let query_start_v = vdupq_n_u32(query_start);
let query_stop_v = vdupq_n_u32(query_stop);

let lane_starts = vld1q_u32(starts.as_ptr().add(lane));
let lane_stops = vld1q_u32(stops.as_ptr().add(lane));

let overlaps = vandq_u32(
    vcgtq_u32(lane_stops, query_start_v),
    vcgtq_u32(query_stop_v, lane_starts),
);
```

Read that literally:

```text
lane_stops > query_start
query_stop > lane_starts
```

NEON comparison results are not scalar booleans. A true `u32` lane contains
`0xffff_ffff`; a false lane contains zero.

For the first version, shift truth down to zero/one and extract each lane:

```rust
let bits = vshrq_n_u32::<31>(overlaps);
mask |= vgetq_lane_u32::<0>(bits) << lane;
mask |= vgetq_lane_u32::<1>(bits) << (lane + 1);
mask |= vgetq_lane_u32::<2>(bits) << (lane + 2);
mask |= vgetq_lane_u32::<3>(bits) << (lane + 3);
```

Keep a scalar tail for a final block whose length is not divisible by four.

### Invariant

The SIMD function must return the exact same `u32` as checkpoint 6. It does not
get a new definition of overlap, block membership, or order.

### Work it by hand

Truth lanes `[true, false, true, true]` become:

```text
[0xffff_ffff, 0, 0xffff_ffff, 0xffff_ffff]
shift by 31 -> [1, 0, 1, 1]
packed mask -> 0b1101
```

### Pass the gate

```bash
cargo test --test block_index
git diff contender/simd-mask -- src/lib.rs tests/block_index.rs
```

On non-AArch64, retain the scalar mask under `cfg` so the branch still compiles.

## Checkpoint 8: replace lane extraction with weighted bits

Each true lane is already all ones. Give each lane the numeric value of its final
bit:

```rust
let weights = vld1q_u32([1_u32, 2, 4, 8].as_ptr());
let bits = vaddvq_u32(vandq_u32(overlaps, weights));
mask |= bits << lane;
```

This is the complete explanation of the apparent magic:

```text
truths:       [all ones, zero, all ones, all ones]
weights:      [1,        2,    4,        8       ]
after AND:    [1,        0,    4,        8       ]
horizontal +: 13
binary 13:    0b1101
```

The sum is a mask, not a count. Powers of two cannot carry into each other
because each weight appears at most once.

### Invariant

Lane `k` contributes exactly `1 << k` if and only if it overlaps.

### Work it by hand

Calculate these without Rust:

```text
[false, true, false, true] -> 2 + 8 = 10 -> 0b1010
[true, true, true, true]   -> 1 + 2 + 4 + 8 = 15 -> 0b1111
[false, false, false, false] -> 0
```

### Pass the gate

```bash
cargo test --test block_index
git diff contender/simd-mask-addv -- src/lib.rs
```

## Checkpoint 9: load and pack eight lanes

First load two adjacent NEON registers at once:

```rust
let lane_starts = vld1q_u32_x2(starts.as_ptr().add(lane));
let lane_stops = vld1q_u32_x2(stops.as_ptr().add(lane));
```

Compare `.0` and `.1` independently. The intermediate paired-load checkpoint
uses two four-lane reductions:

```rust
let bits0 = vaddvq_u32(vandq_u32(overlaps0, weights4));
let bits1 = vaddvq_u32(vandq_u32(overlaps1, weights4));
mask |= (bits0 | (bits1 << 4)) << lane;
```

Then narrow the two all-ones/all-zero vectors into eight `u16` truth lanes and
reduce only once:

```rust
let eight_truths = vcombine_u16(
    vmovn_u32(overlaps0),
    vmovn_u32(overlaps1),
);
let weights8 = vld1q_u16(
    [1_u16, 2, 4, 8, 16, 32, 64, 128].as_ptr(),
);
let bits = vaddvq_u16(vandq_u16(eight_truths, weights8));
mask |= u32::from(bits) << lane;
```

Narrowing preserves truth: the low 16 bits of all ones are still all ones, and
zero stays zero.

Retain this loop structure:

```text
while at least 8 lanes remain -> paired eight-lane path
while at least 4 lanes remain -> single four-lane path
while any lanes remain        -> scalar tail
```

### Invariant

The eight-lane result owns eight consecutive mask bits beginning at `lane`.
Neither paired loads nor narrowing changes which interval belongs to a lane.

### Work it by hand

For truth lanes:

```text
[true, false, true, true, false, false, true, false]
```

the kept weights are `[1, 0, 4, 8, 0, 0, 64, 0]`, whose sum is `77`, or
`0b0100_1101`.

### Pass the gate

```bash
cargo test --test block_index
git diff contender/simd-mask-addv-ld1x2 -- src/lib.rs
git diff contender/simd-mask-narrow8-ld1x2 -- src/lib.rs
```

## Checkpoint 10: inspect the generated AArch64 instructions

Do not add handwritten assembly because the intrinsic source looks verbose.
Inspect what LLVM actually emitted:

```bash
RUSTFLAGS='-C target-cpu=native' \
  cargo rustc --release --lib -- --emit=asm

rg -n 'ldp|uzp1|addv' target/release/deps/rust_lapper-*.s
```

On the measured Apple M3 path:

```text
ldp q..., q...  paired vector loads
uzp1.8h         narrow/join the relevant halfwords
addv.8h         horizontal weighted reduction
```

### Invariant

Intrinsics define semantics; assembly inspection tells you the final cost. Add
inline assembly only when a specific unwanted instruction remains and a measured
replacement wins. That condition was not met here.

### Work it by hand

For one eight-lane group, account for the conceptual work before reading the
assembly:

```text
two vector registers of starts
two vector registers of stops
two end comparisons
two start comparisons
two ANDs
one narrow/join
one weighted AND
one horizontal reduction
```

Then identify which instructions LLVM combined or folded. In particular, verify
that you do not see eight scalar loads or eight scalar branches.

### Pass the gate

Find the paired loads and one eight-halfword reduction in your release assembly.
The randomized correctness test is still mandatory; assembly shape cannot prove
the lane-to-bit mapping.

## Checkpoint 11: make the mask iterator the normal API

The temporary `find_block_mask()` method was useful for A/B testing. Remove it
and put its state into the existing generic `IterFind`. `find()` initializes the
first block from prefix maxima:

```rust
let off = self
    .block_prefix_max_ends
    .partition_point(|max_end| *max_end <= start)
    * INDEX_BLOCK_SIZE;

IterFind {
    inner: self,
    next_block_start: off,
    mask_block_start: 0,
    mask: 0,
    dense_next: 0,
    dense_end: 0,
    backend: detect_backend(),
    start,
    stop,
}
```

`seek()` still uses its sorted-query cursor and `max_len`, but it must round down
to the containing block:

```rust
next_block_start: (*cursor / INDEX_BLOCK_SIZE) * INDEX_BLOCK_SIZE,
```

Rounding up would miss earlier lanes in the cursor's block that can still overlap
the query. The mask will reject lanes that are actually too early.

### Invariant

The API has one query algorithm. Backend selection changes how a mixed mask is
calculated, not whether the block index is used.

### Work it by hand

If `cursor = 45` and block size is 32, the containing block begins at 32, not 64.
A long interval at lane 35 may overlap even though the cursor has advanced to 45.

### Pass the gate

Remove the temporary `got_mask` test path. The ordinary call must now pass the
same oracle:

```bash
cargo test --test block_index
rg -n 'find_block_mask' src tests
```

The final `rg` should return nothing.

## Checkpoint 12: separate the mask backend and support integer widths

Move mask generation to private `src/simd.rs`. Keep one semantic entry point:

```rust
pub(crate) fn overlap_mask<I>(
    backend: MaskBackend,
    starts: &[I],
    stops: &[I],
    query_start: I,
    query_stop: I,
) -> u32
where
    I: PrimInt + 'static,
```

Select only CPU capability:

```rust
pub(crate) enum MaskBackend {
    Scalar,
    #[cfg(target_arch = "aarch64")]
    Neon,
    #[cfg(target_arch = "x86_64")]
    Avx2,
}

pub(crate) fn detect_backend() -> MaskBackend {
    #[cfg(target_arch = "aarch64")]
    {
        return MaskBackend::Neon;
    }

    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            return MaskBackend::Avx2;
        }
    }

    #[allow(unreachable_code)]
    MaskBackend::Scalar
}
```

Use exact `TypeId` checks before reinterpreting a generic slice as a primitive
slice. This is why the generic query implementation acquires `I: 'static`:

```rust
if TypeId::of::<I>() == TypeId::of::<u32>() {
    let starts = unsafe {
        std::slice::from_raw_parts(
            starts.as_ptr().cast::<u32>(),
            starts.len(),
        )
    };
    // Convert stops and query bounds under the same exact-type proof.
    return unsafe { neon::mask_u32(starts, stops, query_start, query_stop) };
}
```

Do not dispatch from `size_of::<I>()` alone. A custom integer wrapper can have
the same size without having a primitive integer's layout. Unknown `PrimInt`
implementations must use the scalar mask.

The AArch64 packing map is:

| Type | Values in paired 128-bit loads | Packing idea |
|---|---:|---|
| `u8` / `i8` | 32 | weighted byte halves |
| `u16` / `i16` | 16 | narrow truth to bytes |
| `u32` / `i32` | 8 | narrow truth to halfwords |
| `u64` / `i64` | 4 | narrow truth to words |

On x86-64 AVX2:

| Type | Values per 256-bit vector | Bit extraction |
|---|---:|---|
| `u8` / `i8` | 32 | `_mm256_movemask_epi8` |
| `u16` / `i16` | 16 | byte movemask, then compact duplicate bits |
| `u32` / `i32` | 8 | `_mm256_movemask_ps` |
| `u64` / `i64` | 4 | `_mm256_movemask_pd` |

AVX2 integer greater-than operations are signed. For unsigned types, XOR both
operands with the sign bit before the signed comparison:

```text
unsigned order: 0 ............ MAX
XOR sign bit:   MIN_SIGNED ... MAX_SIGNED
```

The transform preserves order while moving it into the signed domain.

### Invariant

Every backend returns the same `u32` lane mask. CPU dispatch is not workload
selection. `usize`/`isize` are mapped only to their exact 64-bit representation
on these 64-bit targets.

### Work it by hand

Explain why an AVX2 `u16` comparison produces two identical sign bits in the
byte movemask for each true lane, and why those bits must be compacted to one bit
per interval.

### Pass the gate

```bash
cargo test --all-features
cargo test --target x86_64-apple-darwin --all-features --lib
cargo clippy --target wasm32-wasip1 --lib --all-features -- -D warnings
```

The Rosetta x86 test may take the scalar backend if AVX2 is not advertised. A
forced-Haswell assembly build is a separate instruction check, not native x86
performance evidence.

## Checkpoint 13: make every derived array rebuildable

At this point the query is fast but historical mutation paths leave new metadata
stale. Treat `intervals` as the source of truth and centralize all sidecar
construction:

```rust
fn rebuild_derived(&mut self) {
    let (starts, stops_by_start): (Vec<_>, Vec<_>) = self
        .intervals
        .iter()
        .map(|interval| (interval.start, interval.stop))
        .unzip();

    self.starts = starts;
    self.stops = stops_by_start.clone();
    self.stops_by_start = stops_by_start;

    self.max_len = self
        .intervals
        .iter()
        .map(|interval| {
            interval
                .stop
                .checked_sub(&interval.start)
                .unwrap_or_else(zero::<I>)
        })
        .max()
        .unwrap_or_else(zero::<I>);

    self.stops.sort();

    self.block_max_ends.clear();
    self.block_min_ends.clear();
    let block_count = self.intervals.len().div_ceil(INDEX_BLOCK_SIZE);

    for block in self.intervals.chunks(INDEX_BLOCK_SIZE) {
        let mut max_end = block[0].stop;
        let mut min_end = max_end;
        for interval in &block[1..] {
            max_end = std::cmp::max(max_end, interval.stop);
            min_end = std::cmp::min(min_end, interval.stop);
        }
        self.block_max_ends.push(max_end);
        self.block_min_ends.push(min_end);
    }

    self.block_index.clear();
    self.block_index.resize(block_count, block_count);
    let mut stack = Vec::<usize>::new();
    for block in (0..block_count).rev() {
        while stack.last().is_some_and(|next| {
            self.block_max_ends[*next] <= self.block_max_ends[block]
        }) {
            stack.pop();
        }
        self.block_index[block] =
            stack.last().copied().unwrap_or(block_count);
        stack.push(block);
    }

    self.block_prefix_max_ends
        .clone_from(&self.block_max_ends);
    for block in 1..block_count {
        self.block_prefix_max_ends[block] = std::cmp::max(
            self.block_prefix_max_ends[block - 1],
            self.block_prefix_max_ends[block],
        );
    }
}
```

Preserve the existing `sort_unstable` feature branches around the two sorts in
the actual source.

Call `rebuild_derived()` from:

```text
new()
insert(), after inserting into canonical order
merge_overlaps(), after replacing intervals
deserialize(), by constructing through new()
```

Do not serialize CPU- or implementation-specific sidecars. Custom serde writes
the original six fields and rebuilds the new private fields from `intervals` on
read.

The serialization side deliberately names only the original fields:

```rust
let mut state = serializer.serialize_struct("Lapper", 6)?;
state.serialize_field("intervals", &self.intervals)?;
state.serialize_field("starts", &self.starts)?;
state.serialize_field("stops", &self.stops)?;
state.serialize_field("max_len", &self.max_len)?;
state.serialize_field("cov", &self.cov)?;
state.serialize_field("overlaps_merged", &self.overlaps_merged)?;
state.end()
```

Deserialize those six fields into a helper, but trust only the canonical
intervals as structural input:

```rust
let serialized = SerializedLapper::<I, T>::deserialize(deserializer)?;
let mut lapper = Self::new(serialized.intervals);
lapper.cov = serialized.cov;
lapper.overlaps_merged = serialized.overlaps_merged;
Ok(lapper)
```

The serialized `starts`, `stops`, and `max_len` fields remain in the helper for
wire compatibility, but `new()` recalculates them and every new sidecar. This
also prevents stale or inconsistent serialized metadata from entering the unsafe
query proof.

### Invariant

After every supported mutation:

```text
starts.len() == stops_by_start.len()
metadata.len() == ceil(starts.len() / 32)
starts[i] and stops_by_start[i] describe canonical interval i
every block link is a valid later block or the block-count sentinel
```

### Work it by hand

Insert a new interval at canonical index 17. Every later lane position changes,
and potentially every block minimum, maximum, link, and prefix maximum changes.
Updating only `starts` and `stops` cannot be correct. A full rebuild is simple and
appropriate because `insert()` is already documented as inefficient.

### Pass the gate

Add tests that cross several block boundaries:

```bash
cargo test --all-features insert_rebuilds_every_query_index
cargo test --all-features merge_rebuilds_every_query_index
cargo test --all-features serde_keeps_the_v1_six_field_representation
```

## Checkpoint 14: remove unsigned assumptions

Supporting signed primitive coordinates requires more than writing signed SIMD
comparisons.

Remove the `Unsigned` bound, then repair the scalar assumptions below.

First, a valid signed interval can have a mathematical length larger than
`I::max_value()`. Record that case while rebuilding instead of treating the
failed subtraction as a zero-length interval:

```rust
self.max_len = zero::<I>();
self.max_len_overflowed = false;
for interval in &self.intervals {
    match interval.stop.checked_sub(&interval.start) {
        Some(length) => self.max_len = std::cmp::max(self.max_len, length),
        None if interval.stop >= interval.start => self.max_len_overflowed = true,
        None => {}
    }
}
```

For ordinary `seek()` calls, saturation belongs at the type minimum, not zero.
When `max_len_overflowed` is set, use the exact prefix index because no value of
`I` can provide a conservative length bound:

```rust
if self.max_len_overflowed {
    *cursor = self
        .block_prefix_max_ends
        .partition_point(|max_end| *max_end <= start)
        * INDEX_BLOCK_SIZE;
} else {
    let earliest_start = start
        .checked_sub(&self.max_len)
        .unwrap_or_else(I::min_value);
    // Retain the normal monotonic-cursor search.
}
```

For `count()`, avoid `start + 1`, which can overflow at the coordinate maximum:

```rust
let ends_before_or_at_start = self
    .stops
    .partition_point(|interval_stop| *interval_stop <= start);
let starts_before_stop = Self::bsearch_seq(stop, &self.starts);

let count = starts_before_stop - ends_before_or_at_start;
```

The source expresses the same arithmetic through total length and excluded
suffix count. The identity is:

```text
overlaps = intervals with start < query_stop
         - intervals with stop <= query_start
```

For `depth()`, zero cannot be an uninitialized sentinel. Add an explicit
`initialized: bool` so a coverage interval crossing zero starts at its real
negative coordinate. After advancing to the merged interval's stop, break
before constructing another one-unit query; that stop may be `I::max_value()`.

### Invariant

Zero is an ordinary coordinate. Overflow behavior and iterator initialization
must not give it special semantic meaning.

### Work it by hand

For ends `[-10, -2, 4, 9]` and query start `-2`, exactly the first two ends are
`<= -2` and cannot overlap. No `+1` conversion is required.

### Pass the gate

```bash
cargo test --all-features every_primitive_integer_type_and_block_tail_matches_forward_brute_force
cargo test --all-features signed_seek_saturates_at_the_coordinate_minimum
cargo test --all-features signed_seek_keeps_intervals_whose_length_exceeds_the_coordinate_type
cargo test --all-features signed_depth_crosses_zero_once
cargo test --all-features signed_depth_stops_at_the_coordinate_maximum
```

## Checkpoint 15: remove only proven bounds checks

Finish the safe implementation first and inspect its optimized assembly. The
SIMD helpers already use pointer loads; they do not retain a Rust slice check per
vector lane. The remaining useful reduction was private sidecar indexing in the
block loop.

The proof boundary is:

```rust
let block_start = self.next_block_start;
if block_start >= self.inner.starts.len() {
    return None;
}
```

After construction or supported mutation, private invariants prove:

```text
stops_by_start.len() == starts.len()
each metadata vector has ceil(starts.len() / 32) entries
block = block_start / 32 names the metadata for this boundary
block_end is clamped to starts.len()
```

Only then replace repeated private accesses:

```rust
let first_start = unsafe {
    *self.inner.starts.get_unchecked(block_start)
};
let max_end = unsafe {
    *self.inner.block_max_ends.get_unchecked(block)
};
let starts = unsafe {
    self.inner.starts.get_unchecked(block_start..block_end)
};
```

Keep returned interval indexing safe:

```rust
return Some(&self.inner.intervals[index]);
```

The public `intervals` vector could be edited directly, which has always made
derived metadata stale. Bounding unsafe traversal by the private `starts.len()`
prevents direct growth from extending unchecked sidecars; safe result indexing
turns direct shrinkage into a possible panic rather than unchecked result access.

### Invariant

Every unsafe access must cite a private length relation established by
`rebuild_derived()`. "The benchmark did not crash" is not a safety proof.

### Work it by hand

If `starts.len() = 65`, metadata length is `ceil(65 / 32) = 3`. Valid block
starts are 0, 32, and 64; their block numbers are 0, 1, and 2. The final
`block_end` is `min(96, 65) = 65`, so the final SIMD slice has one lane.

### Pass the gate

```bash
cargo test --all-features
cargo clippy --all-targets --all-features -- -D warnings
git diff contender/simd-mask-narrow8-unchecked-index -- src/lib.rs
```

The measured specialized AArch64 function fell from 442 to 382 static
instructions and from nine bounds-panic edges to two. Making result yields
unchecked reached 363/zero but improved only about 0.3-0.6%, so that broader
unsafe surface was rejected.

## Final trace: run the complete search in your head

Use the toy blocks from checkpoint 3 and query `[23, 28)`.

1. Prefix maxima `[30, 30, 60]` return block 0 as the first possible block.
2. Block 0 has `min=4`, `max=30`, so it is mixed.
3. Its mask is `0b0010`; the iterator saves it.
4. `next()` drains lane 1 and returns canonical interval `[2,30)`.
5. The mask is empty, so the iterator enters block 1.
6. Block 1 has `max=21 <= 23`; its next-greater link jumps to block 2.
7. Block 2 has `min=25 > 23`; every end passes.
8. The start prefix below query stop 28 contains lanes 0, 1, and 2.
9. Three later `next()` calls return `[22,60)`, `[24,25)`, and `[27,35)`.
10. Lane 3 starts at 40, so the globally sorted starts prove the query is done.

The returned canonical indices are `[1, 8, 9, 10]`. That single trace exercises
the mixed mask, saved iterator state, forward block skip, dense prefix, and global
right-hand stopping proof.

## The query algorithms side by side

The sidecars support three distinct public query shapes:

| API | Entry strategy | Result strategy |
|---|---|---|
| `find(start, stop)` | Binary-search `block_prefix_max_ends` | Lazily enumerate borrowed intervals through block routes |
| `seek(start, stop, cursor)` | Reuse a monotonic-query cursor bounded by `max_len`, or use the exact block-prefix search if a signed length overflowed, then round down to a block | Use the same lazy block iterator as `find()` |
| `count(start, stop)` | Binary-search global `starts` and global `stops` | Subtract the two excluded endpoint populations without enumeration |

`find()` answers arbitrary queries independently. `seek()` is the same overlap
problem with an extra promise from the caller: query starts arrive in sorted
order, so a cursor can avoid repeating the entry binary search. It still returns
the same references in the same order.

`count()` uses a different identity because it does not need interval values:

```text
starts_before_query_stop       = number of interval starts < query_stop
ends_at_or_before_query_start  = number of interval stops <= query_start

overlap_count = starts_before_query_stop - ends_at_or_before_query_start
```

For a valid half-open query, every interval in the second population is also in
the first, so the subtraction is nonnegative. `union_and_intersect()` and
`depth()` compose the enumerating APIs; they do not introduce another overlap
index.

## The final data structure at a glance

| State | Order | Purpose |
|---|---|---|
| `intervals` | `(start, stop)` | Canonical storage and returned references |
| `starts` | Same as `intervals` | Prefix searches, stopping, SIMD starts |
| `stops` | Globally sorted | `count()` binary search |
| `stops_by_start` | Same as `intervals` | SIMD ends paired with lanes |
| `block_max_ends` | One per 32 | Prove an all-miss block |
| `block_min_ends` | One per 32 | Prove every end passes |
| `block_index` | One per 32 | First later block with greater max end |
| `block_prefix_max_ends` | Running max | Binary-search first possible block |
| `max_len` | Scalar | Conservative cursor movement in `seek()` |
| `max_len_overflowed` | Scalar flag | Select exact prefix entry when no `I`-sized length bound is conservative |

Complexity is:

```text
build:  O(n log n) sorting + O(n) sidecar construction
space:  O(n) coordinate sidecars + O(n / 32) block metadata
count:  O(log n)
find:   O(log blocks) entry search + candidate block work + output count
seek:   cursor reuse for sorted queries + the same block iterator
```

## Why this remains Lapper-shaped

The final index does not decompose intervals into containment-derived lists like
AIList, build an NCList containment hierarchy, follow SuperIntervals'
previous-greater links backward, or reorder nodes into a tree. It retains one
canonical start-sorted vector, fixed positional block summaries, forward links,
and a query-local mask drained low bit first.

Maximum-end summaries, monotonic stacks, binary search, and SIMD comparison are
general techniques. The structure's identity comes from how they are composed
and from the preserved forward iterator contract.

## Final comprehension checkpoint

Do not look at the answers above while explaining these aloud:

1. Why is `end == query_start` a miss?
2. What does a next-greater link prove about every skipped entry?
3. Why did per-interval links become 32-entry block links?
4. Why is a prefix maximum needed in addition to each block maximum?
5. Why does the dense route still need an active start prefix?
6. Why does the mixed route not need a separate active prefix?
7. Why must `stops` and `stops_by_start` both exist?
8. Why does clearing the lowest set bit preserve forward order?
9. Why is the weighted horizontal sum a mask rather than a count?
10. Why is narrowing all-ones truth from `u32` to `u16` lossless?
11. Why must `seek()` round its cursor down to a block boundary?
12. What makes CPU dispatch different from a workload mode switch?
13. Which invariant makes each `get_unchecked` access valid?
14. Why are returned intervals still safely indexed?
15. Which arrays must change after one insertion, and why is rebuilding simpler?

Once you can answer all fifteen comprehension questions and reproduce the final
toy trace, read
[`PORTABLE_SIMD_INDEX.md`](PORTABLE_SIMD_INDEX.md) for the engineering summary
and compare your complete practice worktree with
`worked/portable-simd-index`.
