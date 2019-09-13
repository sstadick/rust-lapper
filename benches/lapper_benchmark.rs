#[macro_use]
extern crate criterion;
extern crate rand;
extern crate rust_lapper;

use bio::data_structures::interval_tree::IntervalTree;
use coitree::{COITree, IntervalNode};
use cpu_time::ProcessTime;
use criterion::black_box;
use criterion::Criterion;
use nested_intervals::IntervalSet;
use rand::prelude::*;
use rand::Rng;
use rust_lapper::{Interval, Lapper};
use std::ops::Range;
use std::time::Duration;

type Iv = Interval<bool>;

fn randomi(imin: usize, imax: usize) -> usize {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: usize, range_max: usize, size_min: usize, size_max: usize) -> Vec<Iv> {
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let s = randomi(0, range_max);
        let e = s + randomi(size_min, size_max);
        result.push(Interval {
            start: s,
            stop: e,
            val: false,
        });
    }
    result
}

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>, Vec<Range<u32>>, Vec<Range<u32>>) {
    let n = 3_000_000;
    let chrom_size = 100_000_000;
    let min_interval_size = 500;
    let max_interval_size = 80000;
    let intervals = make_random(n, chrom_size, min_interval_size, max_interval_size);
    let other_intervals = make_random(n, 10 * chrom_size, 1, 2);
    let nested_intervals = intervals
        .iter()
        .map(|x| x.start as u32..x.stop as u32)
        .collect();
    let nested_other_intervals = other_intervals
        .iter()
        .map(|x| x.start as u32..x.stop as u32)
        .collect();
    (
        intervals,
        other_intervals,
        nested_intervals,
        nested_other_intervals,
    )
}

pub fn query(c: &mut Criterion) {
    let (intervals, other_intervals, nested_intervals, nested_other_intervals) =
        make_interval_set();
    let mut bad_nested_intervals = nested_intervals.clone();
    bad_nested_intervals.push(0..90_000_000);
    let mut bad_intervals = intervals.clone();

    // Make Lapper intervals
    let lapper = Lapper::new(intervals);
    let other_lapper = Lapper::new(other_intervals);
    bad_intervals.push(Iv {
        start: 0,
        stop: 90_000_000,
        val: false,
    });
    let bad_lapper = Lapper::new(bad_intervals);
    // Make COITree intervals
    let coi_intervals: Vec<IntervalNode<bool>> = nested_intervals
        .clone()
        .into_iter()
        .map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
        .collect();
    let coi_tree_intervals_preserved = coi_intervals.clone();
    let coi_other_intervals: Vec<IntervalNode<bool>> = nested_other_intervals
        .clone()
        .into_iter()
        .map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
        .collect();
    let coi_bad_itervals: Vec<IntervalNode<bool>> = bad_nested_intervals
        .clone()
        .into_iter()
        .map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
        .collect();
    let coi_tree = COITree::new(coi_intervals);
    let coi_bad_tree = COITree::new(coi_bad_itervals);

    // Make NesteInterval / Bio Intervals
    let mut nested_interval_set = IntervalSet::new(&nested_intervals).unwrap();
    let mut nested_other_interval_set = IntervalSet::new(&nested_other_intervals).unwrap();
    let mut nested_bad_interval_set = IntervalSet::new(&bad_nested_intervals).unwrap();

    let mut bio_interval_tree = IntervalTree::new();
    nested_intervals
        .iter()
        .for_each(|x| bio_interval_tree.insert(x, x.start));
    let mut bio_other_interval_tree = IntervalTree::new();
    nested_other_intervals
        .iter()
        .for_each(|x| bio_other_interval_tree.insert(x, x.start));
    let mut bio_bad_interval_tree = IntervalTree::new();
    bad_nested_intervals
        .iter()
        .for_each(|x| bio_bad_interval_tree.insert(x, x.start));

    let start = ProcessTime::now();
    println!("Starting timer");
    let mut count = 0;
    for x in lapper.iter() {
        count += lapper.find(x.start, x.stop).count();
    }
    let cpu_time: Duration = start.elapsed();
    println!("rust-lapper: 100% hit rate: {:#?}", cpu_time);
    println!("Found {}", count);

    let start = ProcessTime::now();
    println!("Starting timer");
    let mut count = 0;
    for x in coi_tree_intervals_preserved.iter() {
        count += coi_tree.find(x.start, x.stop).count();
    }
    let cpu_time: Duration = start.elapsed();
    println!("COITree: 100% hit rate: {:#?}", cpu_time);
    println!("Found {}", count);

    let start = ProcessTime::now();
    println!("Starting timer");
    let mut count = 0;
    for x in nested_intervals.iter() {
        count += nested_interval_set.query_overlapping(x).iter().count();
    }
    let cpu_time: Duration = start.elapsed();
    println!("Nested Intervals: 100% hit rate: {:#?}", cpu_time);
    println!("Found {}", count);

    let start = ProcessTime::now();
    println!("Starting timer");
    let mut count = 0;
    for x in nested_intervals.iter() {
        count += bio_interval_tree.find(x).count();
    }
    let cpu_time: Duration = start.elapsed();
    println!("Bio Intervals: 100% hit rate: {:#?}", cpu_time);
    println!("Found {}", count);

    println!(">>>>>>>>>>>>>> Done");

    let mut comparison_group = c.benchmark_group("Bakeoff");
    comparison_group.bench_function("rust-lapper: find with 100% hit rate", |b| {
        b.iter(|| {
            for x in lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });

    comparison_group.bench_function("rust-lapper: find with below 100% hit rate", |b| {
        b.iter(|| {
            for x in other_lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });
    comparison_group.bench_function("rust-lapper: seek with 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            for x in lapper.iter() {
                lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });
    comparison_group.bench_function("rust-lapper: seek with below 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            for x in other_lapper.iter() {
                lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });
    comparison_group.bench_function(
        "rust-lapper: find with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bad_lapper.find(x.start, x.stop).count();
                }
            });
        },
    );
    comparison_group.bench_function(
        "rust-lapper: seek with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                let mut cursor = 0;
                for x in other_lapper.iter() {
                    bad_lapper.seek(x.start, x.stop, &mut cursor).count();
                }
            });
        },
    );
    comparison_group.bench_function("COITree: query_overlapping with 100% hit rate", |b| {
        b.iter(|| {
            for x in coi_tree_intervals_preserved.iter() {
                coi_tree.find(x.start, x.stop).count();
            }
        });
    });
    comparison_group.bench_function("COITree: query_overlapping with below 100% hit rate", |b| {
        b.iter(|| {
            for x in coi_other_intervals.iter() {
                coi_tree.find(x.start, x.stop).count();
            }
        });
    });
    comparison_group.bench_function(
        "COITree: query_overlapping with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in coi_other_intervals.iter() {
                    coi_bad_tree.find(x.start, x.stop).count();
                }
            });
        },
    );
    comparison_group.bench_function(
        "nested_intervals: query_overlapping with 100% hit rate",
        |b| {
            b.iter(|| {
                for x in nested_intervals.iter() {
                    nested_interval_set.query_overlapping(x).iter().count();
                }
            });
        },
    );
    comparison_group.bench_function(
        "nested_intervals: query_overlapping with below 100% hit rate",
        |b| {
            b.iter(|| {
                for x in nested_other_intervals.iter() {
                    nested_interval_set.query_overlapping(x).iter().count();
                }
            });
        },
    );
    comparison_group.bench_function("nested_intervals: query_overlapping with below 100% hit rate - chromosome spanning interval", |b| {
        b.iter(|| {
            for x in nested_other_intervals.iter() {
                nested_bad_interval_set.query_overlapping(x).iter().count();
            }
        });
    });
    comparison_group.bench_function("rust-bio: find with 100% hit rate", |b| {
        b.iter(|| {
            for x in nested_intervals.iter() {
                bio_interval_tree.find(x).count();
            }
        });
    });
    comparison_group.bench_function("rust-bio: find with below 100% hit rate", |b| {
        b.iter(|| {
            for x in nested_other_intervals.iter() {
                bio_interval_tree.find(x).count();
            }
        });
    });
    comparison_group.bench_function(
        "rust-bio: find with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in nested_other_intervals.iter() {
                    bio_bad_interval_tree.find(x).count();
                }
            });
        },
    );
    comparison_group.finish();

    let mut group = c.benchmark_group("Lapper find overlaps");
    group.bench_function("find with 100% hit rate", |b| {
        b.iter(|| {
            for x in lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });

    group.bench_function("find with below 100% hit rate", |b| {
        b.iter(|| {
            for x in other_lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });

    group.bench_function(
        "find with 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in lapper.iter() {
                    bad_lapper.find(x.start, x.stop).count();
                }
            });
        },
    );

    group.bench_function("seek with 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            for x in lapper.iter() {
                lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });

    group.bench_function("seek with below 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            for x in other_lapper.iter() {
                lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });
    group.finish();
}

criterion_group!(benches, query);
criterion_main!(benches);
