#[macro_use]
extern crate criterion;
extern crate rust_lapper;
extern crate rand;

use criterion::Criterion;
use criterion::black_box;
use rust_lapper::{Interval, Lapper};
use rand::prelude::*;
use rand::Rng;
use std::time::Duration;
use cpu_time::ProcessTime;

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
        result.push(Interval{start: s, stop: e, val: false});
    }
    result
}

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>){
    let n = 3_000_000;
    //let n = 1_000;
    let chrom_size = 100_000_000;
    let min_interval_size = 500;
    let max_interval_size = 80000;
    let intervals = make_random(n, chrom_size, min_interval_size, max_interval_size);
    let other_intervals = make_random(n, 10 * chrom_size, 1, 2);
    (intervals, other_intervals)
}

pub fn query(c: &mut Criterion) {
    let (intervals, other_intervals) = make_interval_set();
    let lapper = Lapper::new(intervals);
    let other_lapper = Lapper::new(other_intervals);

    let start = ProcessTime::now();
    println!("Starting timer");
    let mut count = 0;
    for x in lapper.iter() {
        count += lapper.find(x.start, x.stop).count();
    }
    let cpu_time: Duration = start.elapsed();
    println!("100% hit rate: {:#?}", cpu_time);
    println!("Found {}", count);

    c.bench_function("find with 100% hit rate", |b| {
        b.iter(|| {
            let mut count = 0;
            for x in lapper.iter() {
                count += lapper.find(x.start, x.stop).count();
            }
        });
    });

    c.bench_function("find with below 100% hit rate", |b| {
        b.iter(|| {
            let mut count = 0;
            for x in other_lapper.iter() {
                count += lapper.find(x.start, x.stop).count();
            }
        });
    });

    c.bench_function("seek with 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            let mut count = 0;
            for x in lapper.iter() {
                count += lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });

    c.bench_function("seek with below 100% hit rate", |b| {
        b.iter(|| {
            let mut cursor = 0;
            let mut count = 0;
            for x in other_lapper.iter() {
                count += lapper.seek(x.start, x.stop, &mut cursor).count();
            }
        });
    });
}


criterion_group!(benches, query);
criterion_main!(benches);

