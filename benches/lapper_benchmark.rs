#[macro_use]
extern crate criterion;
extern crate rand;
extern crate rust_lapper;

use cpu_time::ProcessTime;
use criterion::black_box;
use criterion::Criterion;
use rand::prelude::*;
use rand::Rng;
use rust_lapper::{Interval, Lapper};
use std::ops::Range;
use std::time::Duration;

type Iv = Interval<bool>;

fn randomi(imin: u32, imax: u32) -> u32 {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: usize, range_max: u32, size_min: u32, size_max: u32) -> Vec<Iv> {
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

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>) {
    //let n = 3_000_000;
    let n = 1_000;
    let chrom_size = 100_000_000;
    let min_interval_size = 500;
    let max_interval_size = 80000;
    let intervals = make_random(n, chrom_size, min_interval_size, max_interval_size);
    let other_intervals = make_random(n, 10 * chrom_size, 1, 2);
    (intervals, other_intervals)
}

pub fn query(c: &mut Criterion) {
    let (intervals, other_intervals) = make_interval_set();
    // Make Lapper intervals
    let mut bad_intervals = intervals.clone();
    let lapper = Lapper::new(intervals);
    let other_lapper = Lapper::new(other_intervals);
    bad_intervals.push(Iv {
        start: 0,
        stop: 90_000_000,
        val: false,
    });
    let bad_lapper = Lapper::new(bad_intervals);
    //let start = ProcessTime::now();
    //println!("Starting timer");
    //let mut count = 0;
    //for x in lapper.iter() {
    //count += lapper.find(x.start, x.stop).count();
    //}
    //let cpu_time: Duration = start.elapsed();
    //println!("100% hit rate: {:#?}", cpu_time);
    //println!("Found {}", count);

    let mut comparison_group = c.benchmark_group("Bakeoff");
    comparison_group.bench_function("rust-lapper: find with 100% hit rate", |b| {
        b.iter(|| {
            for x in lapper.iter() {
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
    comparison_group.bench_function("rust-lapper: count with 100% hit rate", |b| {
        b.iter(|| {
            for x in lapper.iter() {
                lapper.count(x.start, x.stop);
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
    comparison_group.bench_function("rust-lapper: seek with below 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            for x in other_lapper.iter() {
                lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });
    comparison_group.bench_function("rust-lapper: count with below 100% hit rate", |b| {
        b.iter(|| {
            for x in other_lapper.iter() {
                lapper.count(x.start, x.stop);
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
    comparison_group.bench_function(
        "rust-lapper: count with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bad_lapper.count(x.start, x.stop);
                }
            });
        },
    );
    comparison_group.finish();
}

criterion_group!(benches, query);
criterion_main!(benches);
