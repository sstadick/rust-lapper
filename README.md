# rust-lapper

<p align="center">
  <a href="https://github.com/sstadick/rust-lapper/actions?query=workflow%3Aci"><img src="https://github.com/sstadick/rust-lapper/workflows/ci/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/rust-lapper.svg" alt="license">
  <a href="https://crates.io/crates/rust-lapper"><img src="https://img.shields.io/crates/v/rust-lapper.svg?colorB=319e8c" alt="Version info"></a><br>
</p>

[Documentation](https://docs.rs/rust-lapper)
[Crates.io](https://crates.io/crates/rust-lapper)

This was originally a Rust port of Brent Pedersen's
[nim-lapper](https://github.com/brentp/nim-lapper). `find()` and `seek()` return
lazy borrowed iterators in ascending start order, so normal iterator adaptors
work without collecting results first.

All stored intervals and query ranges use half-open `[start, stop)` semantics.
`Lapper` keeps its intervals sorted by start and builds a fixed 32-interval
block index that can skip regions proven not to overlap. Mixed blocks use NEON
on AArch64, runtime-detected AVX2 on x86-64, and an exact scalar fallback
elsewhere. The same algorithm handles both ordinary and pathological datasets
with long intervals that engulf many shorter intervals.

The `count()` method uses the
[BITS algorithm](https://academic.oup.com/bioinformatics/article/29/1/1/273289)
to count overlaps with two binary searches.

## API and algorithm compatibility

The block index and SIMD backends are private implementation details: there is
no mode flag, alternate query method, or architecture-specific API. Existing
call patterns for `find()`, `seek()`, `count()`, `cov()`, `set_cov()`,
`merge_overlaps()`, `depth()`, `union_and_intersect()`, `union()`, and
`intersect()` retain their range semantics and return types. `count()` remains
the independent BITS implementation; methods that use `find()` or `seek()`
internally automatically share the exact indexed query path.

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
the scalar mask implementation. The bound is the API-breaking change that makes
this a major release.

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

The retained v2 release measurements, raw samples, compiler flags, and pinned
competitor revisions live in
[lapper_bakeoff](https://github.com/sstadick/lapper_bakeoff/tree/main/results/avx2-2026-07-28).
On an AMD Ryzen 9 3950X with AVX2, the new implementation improved total time
over rust-lapper 1.3.0 by 37.30%, 98.89%, and 34.64% on the three retained
article cases. All implementations returned identical overlap counts.

Benchmark results are workload- and hardware-specific; use the linked harness
and raw data when making comparisons.

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

- `2.0.0`: Replace the longest-interval scan with a portable SIMD block index,
  add signed coordinates, declare Rust 1.59 as the MSRV, and accept the
  `I: 'static` coordinate bound.
- `1.3.0`: Add the `sort_unstable` feature flag for allocation-sensitive sorting thanks to @jameslkingsley.
- `1.1.0`: Added insert functionality thanks to @zaporter
- `1.0.0`: Add serde support via the `with_serde` feature flag
- `0.5.0`: Make Interval start/stop generic
- `0.4.3`: Remove leftover print statement
- `0.4.2`: Bugfix in to update starts/stops vectors when overlaps merged
- `0.4.0`: Addition of the BITS count algorithm.
