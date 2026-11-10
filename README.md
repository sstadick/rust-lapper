# rust-lapper

<p align="center">
  <a href="https://github.com/sstadick/rust-lapper/actions?query=workflow%3Aci"><img src="https://github.com/sstadick/rust-lapper/workflows/ci/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/rust-lapper.svg" alt="license">
  <a href="https://crates.io/crates/rust-lapper"><img src="https://img.shields.io/crates/v/rust-lapper.svg?colorB=319e8c" alt="Version info"></a><br>
</p>

[Documentation](https://docs.rs/rust-lapper)
[Crates.io](https://crates.io/crates/rust-lapper)

rust-lapper 2 is stable. Add `rust-lapper = "2"` to your dependencies.

This was originally a Rust port of Brent Pedersen's
[nim-lapper](https://github.com/brentp/nim-lapper). `find()` and `seek()` return
lazy borrowed iterators in ascending start order, so normal iterator adaptors
work without collecting results first.

All stored intervals and query ranges use half-open `[start, stop)` semantics.
`Lapper` keeps its intervals sorted by start and indexes them in blocks of 32.
A max tree over the blocks' end positions lets a query jump over any run of
blocks that cannot overlap, in time logarithmic in the length of the run.
Blocks that might overlap are checked 32 intervals at a time with NEON on
AArch64, runtime-detected AVX2 on x86-64, and an exact scalar fallback
elsewhere. The same algorithm handles both ordinary datasets and datasets with
long intervals that engulf many shorter intervals.

The `count()` method uses the
[BITS algorithm](https://academic.oup.com/bioinformatics/article/29/1/1/273289)
to count overlaps with two binary searches.

## Query backends

The block index and SIMD backends are private implementation details: there is
no mode flag, alternate query method, or architecture-specific API. `find()`
and `seek()` use the index directly, and `depth()`, `union_and_intersect()`,
`union()`, and `intersect()` use it through them. `count()` is an independent
BITS implementation.

| Target | Mixed-block backend | Selection |
|---|---|---|
| AArch64 | 128-bit NEON | Baseline for the architecture |
| x86-64 with AVX2 | 256-bit AVX2 | Runtime detected once per iterator |
| x86-64 without AVX2 | Scalar | Automatic fallback |
| Other architectures | Scalar | Automatic fallback |

NEON and AVX2 cover `u8`, `i8`, `u16`, `i16`, `u32`, `i32`, `u64`, `i64`,
`usize`, and `isize`. The same block algorithm uses exact scalar masks for
`u128`, `i128`, custom `PrimInt` types, and partial vector tails.

## Minimum Supported Rust Version

rust-lapper 2 supports Rust 1.59 and newer. Rust 1.59 is the first stable
release that provides the AArch64 intrinsics used by the NEON query backend.

Query coordinates must be `'static` so private dispatch code can use `TypeId`
before reinterpreting primitive integer slices for SIMD. This includes every
primitive integer and ordinary owned custom numeric type; it does not require a
`Lapper` value to live for the entire program. Non-primitive `PrimInt` types use
the scalar mask implementation.

## Mutation

Use `insert()` and `merge_overlaps()` for coordinate or structural changes so
the private query index is rebuilt. Directly changing `Lapper::intervals`
coordinates or length leaves derived metadata stale; changing payload values is
safe.

## Serde Support

`rust-lapper` supports serialization with serde for `Lapper` and `Interval` objects:

```toml
[dependencies]
rust-lapper = { version = "2", features = ["with_serde"] }
```

See `examples/serde.rs` for a brief example.

## Benchmarks

The table shows one index build and one query batch for the three overlap cases
of the
[September 2025 polars-bio benchmark](https://biodatageeks.org/polars-bio/blog/2025/09/05/interval-operations-benchmark--update-september-2025/),
timed after the input is loaded. Each case builds one `Lapper<u32, ()>` per
chromosome and sums `find(start, stop).count()` over every query interval.

| Case | Indexed intervals | Query intervals | Overlaps | Build | Query | Total |
|---|---:|---:|---:|---:|---:|---:|
| `1-2` | 438,694 | 198,621 | 54,246 | 5.8 ms | 4.6 ms | 10.4 ms |
| `7-3` | 1,956,864 | 1,194,285 | 4,408,383 | 54.1 ms | 55.1 ms | 109.2 ms |
| `8-7` | 1,194,285 | 9,944,559 | 307,184,634 | 57.2 ms | 690.5 ms | 747.7 ms |

Times are medians on one core of an AMD Ryzen 9 3950X with AVX2, built with
Rust 1.95, thin LTO, and `-C target-cpu=native`. Benchmark results are workload-
and hardware-specific. `cargo bench` runs the crate's synthetic benchmarks.

## Example

```rust
use rust_lapper::{Interval, Lapper};

type Iv = Interval<usize, u32>;
fn main() {
    // create some fake data
    let data: Vec<Iv> = vec![
        Iv {
            start: 70,
            stop: 120,
            val: 0,
        }, // a long interval
        Iv {
            start: 10,
            stop: 15,
            val: 0,
        },
        Iv {
            start: 10,
            stop: 15,
            val: 0,
        }, // exact overlap
        Iv {
            start: 12,
            stop: 15,
            val: 0,
        }, // inner overlap
        Iv {
            start: 14,
            stop: 16,
            val: 0,
        }, // overlap end
        Iv {
            start: 40,
            stop: 45,
            val: 0,
        },
        Iv {
            start: 50,
            stop: 55,
            val: 0,
        },
        Iv {
            start: 60,
            stop: 65,
            val: 0,
        },
        Iv {
            start: 68,
            stop: 71,
            val: 0,
        }, // overlap start
        Iv {
            start: 70,
            stop: 75,
            val: 0,
        },
    ];

    // make lapper structure
    let mut lapper = Lapper::new(data);

    // Find every interval that overlaps [11, 15).
    // For queries in nondecreasing start order, seek() can reuse a caller-owned cursor.
    assert_eq!(
        lapper.find(11, 15).collect::<Vec<&Iv>>(),
        vec![
            &Iv {
                start: 10,
                stop: 15,
                val: 0
            },
            &Iv {
                start: 10,
                stop: 15,
                val: 0
            }, // exact overlap
            &Iv {
                start: 12,
                stop: 15,
                val: 0
            }, // inner overlap
            &Iv {
                start: 14,
                stop: 16,
                val: 0
            }, // overlap end
        ]
    );
    assert_eq!(lapper.count(11, 15), 4);

    // Merge overlapping regions to simplify queries that only depend on whether
    // any interval overlaps.
    lapper.merge_overlaps();
    assert_eq!(
        lapper.find(11, 15).collect::<Vec<&Iv>>(),
        vec![&Iv {
            start: 10,
            stop: 16,
            val: 0
        },]
    );

    // Get the number of positions covered by the interval collection.
    assert_eq!(lapper.cov(), 73);

    // Get the union and intersection lengths of two interval collections.
    let data = vec![
        Iv {
            start: 5,
            stop: 15,
            val: 0,
        },
        Iv {
            start: 48,
            stop: 80,
            val: 0,
        },
    ];
    let (union, intersect) = lapper.union_and_intersect(&Lapper::new(data));
    assert_eq!(union, 88);
    assert_eq!(intersect, 27);

    // Get the depth at each position covered by the lapper
    for interval in lapper.depth().filter(|x| x.val > 2) {
        println!(
            "Depth at {} - {}: {}",
            interval.start, interval.stop, interval.val
        );
    }

}
```

## Release Notes

- Unreleased: Skip runs of blocks that cannot overlap with a max tree. A skip
  is logarithmic in its length for every data shape, including many short
  intervals nested under a long one.
- `2.0.0`: Promote `2.0.0-beta.1` to stable with no implementation changes.
  Includes the portable SIMD block index and signed-coordinate support.
  Requires Rust 1.59 or newer and the `I: 'static` coordinate bound.
- `2.0.0-beta.1`: Begin the rust-lapper 2 beta with a portable SIMD block
  index, signed coordinates, Rust 1.59 as the MSRV, and the `I: 'static`
  coordinate bound.
- `1.3.0`: Add the `sort_unstable` feature flag for allocation-sensitive sorting thanks to @jameslkingsley.
- `1.1.0`: Added insert functionality thanks to @zaporter
- `1.0.0`: Add serde support via the `with_serde` feature flag
- `0.5.0`: Make Interval start/stop generic
- `0.4.3`: Remove leftover print statement
- `0.4.2`: Bugfix in to update starts/stops vectors when overlaps merged
- `0.4.0`: Addition of the BITS count algorithm.
