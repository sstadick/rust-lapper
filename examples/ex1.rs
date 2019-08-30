use rust_lapper::{Interval, Lapper};
use std::cmp;

type Iv = Interval<u32>;
fn main() {
    // create some fake data
    let data: Vec<Iv> = (0..20)
        .step_by(5)
        .map(|x| Iv {
            start: x,
            stop: x + 2,
            val: 0,
        })
        .collect();
    println!("{:#?}", data);

    // make lapper structure
    let laps = Lapper::new(data);

    assert_eq!(
        laps.find(6, 11).next(),
        Some(&Iv {
            start: 5,
            stop: 7,
            val: 0
        })
    );

    let mut sim: i32 = 0;
    let mut cursor = 0;
    // Calculate the overlap between the query and the found intervals, sum total overlap
    for i in (0..10).step_by(3) {
        sim += laps
            .seek(i, i + 2, &mut cursor)
            .map(|iv| cmp::min(i + 2, iv.stop) - cmp::max(i, iv.start))
            .sum::<i32>();
    }
    assert_eq!(sim, 10);
}
