# rust-lapper

<p align="center">
  <a href="https://github.com/sstadick/rust-lapper/actions?query=workflow%3Aci"><img src="https://github.com/sstadick/rust-lapper/workflows/ci/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/rust-lapper.svg" alt="license">
  <a href="https://crates.io/crates/rust-lapper"><img src="https://img.shields.io/crates/v/rust-lapper.svg?colorB=319e8c" alt="Version info"></a><br>
</p>

[Documentation](https://docs.rs/rust-lapper)
[Crates.io](https://crates.io/crates/rust-lapper)

This is a rust port of Brent Pendersen's
[nim-lapper](https://github.com/brentp/nim-lapper). It has a few notable
differences, mostly that the find and seek methods both return
iterators, so all adaptor methods may be used normally.

This crate works well for most interval data that does not include very long
intervals that engulf a majority of other intervals. It is still fairly
comparable to other methods. If you absolutely need time guarantees in the
worst case, see [COItres](https://github.com/dcjones/coitrees) and [IITree](https://docs.rs/bio/0.32.0/bio/data_structures/interval_tree/struct.ArrayBackedIntervalTree.html).

However, on more typical datasets, this crate is between 4-10x faster
than other interval overlap methods.

It should also be noted that the `count` method is agnostic to data
type, and should be about as fast as it is possible to be on any
dataset. It is an implementation of the [BITS
algorithm](https://academic.oup.com/bioinformatics/article/29/1/1/273289)

## Serde Support

`rust-lapper` supports serialization with serde for `Lapper` and `Interval` objects:

```toml
[dependencies]
rust-lapper = { version = "*", features = ["with_serde"] }
```

See `examples/serde.rs` for a brief example.

## Benchmarks

Benchmarking interval tree-ish datastructures is hard
Please see the
[interval_bakeoff](https://github.com/sstadick/interval_bakeoff) project
for details on how the benchmarks were run... It's not fully baked yet
though, and is finiky to run.

Command to run:

```
./target/release/interval_bakeoff fake -a -l RustLapper -l
RustBio -l NestedInterval -n50000 -u100000

# This equates to the following params:
# num_intervals	50000
# universe_size	100000
# min_interval_size	500
# max_interval_size	80000
# add_large_span	true (universe spanning)
```

Set A / b Creation Times

| crate/method     | A time   | B time   |
| ---------------- | -------- | -------- |
| rust_lapper      | 15.625ms | 31.25ms  |
| nested_intervals | 15.625ms | 15.625ms |
| bio              | 15.625ms | 31.25ms  |

100% hit rate (A vs A)

| crate/method                       | mean time  | intersection |
| ---------------------------------- | ---------- | ------------ |
| rust_lapper/find                   | 4.78125s   | 1469068763   |
| rust_lapper/count                  | 15.625ms   | 1469068763   |
| nested_intervals/query_overlapping | 157.4375s  | 1469068763   |
| bio/find                           | 33.296875s | 1469068763   |


Sub 100% hit rate (A vs B)

| crate/method                       | mean time  | intersection |
| ---------------------------------- | ---------- | ------------ |
| rust_lapper/find                   | 531.25ms   | 176488436    |
| rust_lapper/count                  | 15.625ms   | 176488436    |
| nested_intervals/query_overlapping | 11.109375s | 196090092    |
| bio/find                           | 4.3125s    | 176488436    |

[nested_intervals](https://docs.rs/nested_intervals/0.2.0/nested_intervals/)
[rust-bio](https://docs.rs/bio/0.28.2/bio/)
*Note that rust-bio has a new interval tree structure which should be faster than what is shown here*

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
        }, // max_len = 50
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

    // Iterator based find to extract all intervals that overlap 6..7
    // If your queries are coming in start sorted order, use the seek method to retain a cursor for
    // a big speedup.
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

    // Merge overlaping regions within the lapper to simplifiy and speed up quries that only depend
    // on 'any
    lapper.merge_overlaps();
    assert_eq!(
        lapper.find(11, 15).collect::<Vec<&Iv>>(),
        vec![&Iv {
            start: 10,
            stop: 16,
            val: 0
        },]
    );

    // Get the number of positions covered by the lapper tree:
    assert_eq!(lapper.cov(), 73);

    // Get the union and intersect of two different lapper trees
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

- `1.1.0`: Added insert functionality thanks to @zaporter
- `0.4.0`: Addition of the BITS count algorithm.
- `0.4.2`: Bugfix in to update starts/stops vectors when overlaps merged
- `0.4.3`: Remove leftover print statement
- `0.5.0`: Make Interval start/stop generic
- `1.0.0`: Add serde support via the `with_serde` feature flag
