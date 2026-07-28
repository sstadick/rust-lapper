use num_traits::PrimInt;
use rust_lapper::{Interval, Lapper};

fn assert_queries_match<I>(lapper: &Lapper<I, usize>, query_min: i64, query_max: i64)
where
    I: PrimInt + Ord + Clone + Send + Sync + 'static,
{
    let mut cursor = 0;
    for raw_start in query_min..query_max {
        let start = I::from(raw_start).unwrap();
        let stop = I::from(raw_start + 7).unwrap();
        let expected: Vec<_> = lapper
            .intervals
            .iter()
            .filter(|interval| interval.start < stop && interval.stop > start)
            .map(|interval| interval.val)
            .collect();
        let found: Vec<_> = lapper
            .find(start, stop)
            .map(|interval| interval.val)
            .collect();
        let sought: Vec<_> = lapper
            .seek(start, stop, &mut cursor)
            .map(|interval| interval.val)
            .collect();
        assert_eq!(found, expected, "find query {raw_start}..{}", raw_start + 7);
        assert_eq!(
            sought,
            expected,
            "seek query {raw_start}..{}",
            raw_start + 7
        );
        assert_eq!(
            lapper.count(start, stop),
            expected.len(),
            "count query {raw_start}..{}",
            raw_start + 7
        );
    }
}

fn exercise_unsigned<I>()
where
    I: PrimInt + Ord + Clone + Send + Sync + 'static,
{
    let intervals = (0..96)
        .map(|value| {
            let raw_start = (value * 2) % 181;
            let raw_stop = raw_start + 1 + (value * 7) % 19;
            Interval {
                start: I::from(raw_start).unwrap(),
                stop: I::from(raw_stop).unwrap(),
                val: value,
            }
        })
        .collect();
    assert_queries_match(&Lapper::new(intervals), 0, 190);
}

fn exercise_signed<I>()
where
    I: PrimInt + Ord + Clone + Send + Sync + 'static,
{
    let intervals = (0..96)
        .map(|value| {
            let raw_start = (value * 2) as i64 % 121 - 60;
            let raw_stop = raw_start + 1 + (value * 7) as i64 % 19;
            Interval {
                start: I::from(raw_start).unwrap(),
                stop: I::from(raw_stop).unwrap(),
                val: value,
            }
        })
        .collect();
    assert_queries_match(&Lapper::new(intervals), -70, 70);
}

#[test]
fn every_primitive_integer_type_matches_forward_brute_force() {
    exercise_unsigned::<u8>();
    exercise_unsigned::<u16>();
    exercise_unsigned::<u32>();
    exercise_unsigned::<u64>();
    exercise_unsigned::<u128>();
    exercise_unsigned::<usize>();
    exercise_signed::<i8>();
    exercise_signed::<i16>();
    exercise_signed::<i32>();
    exercise_signed::<i64>();
    exercise_signed::<i128>();
    exercise_signed::<isize>();
}

#[test]
fn insert_rebuilds_every_query_index() {
    let mut lapper = Lapper::<u32, usize>::new(Vec::new());
    for value in (0..129).rev() {
        let start = (value * 13 % 257) as u32;
        lapper.insert(Interval {
            start,
            stop: start + 1 + (value * 11 % 31) as u32,
            val: value,
        });
    }
    assert_queries_match(&lapper, 0, 280);
}

#[test]
fn merge_rebuilds_every_query_index() {
    let intervals = (0..160)
        .map(|value| {
            let start = (value * 3) as u32;
            Interval {
                start,
                stop: start + 5 + (value % 9) as u32,
                val: value,
            }
        })
        .collect();
    let mut lapper = Lapper::new(intervals);
    lapper.merge_overlaps();
    assert_queries_match(&lapper, 0, 500);
}

#[test]
fn signed_seek_saturates_at_the_coordinate_minimum() {
    let lapper = Lapper::new(vec![
        Interval {
            start: i8::MIN,
            stop: -100,
            val: 0_usize,
        },
        Interval {
            start: -110,
            stop: -90,
            val: 1,
        },
    ]);
    let mut cursor = 0;
    let found: Vec<_> = lapper
        .seek(-125, -105, &mut cursor)
        .map(|interval| interval.val)
        .collect();
    assert_eq!(found, vec![0, 1]);
}

#[test]
fn signed_depth_crosses_zero_once() {
    let lapper = Lapper::new(vec![Interval {
        start: -2_i16,
        stop: 3,
        val: (),
    }]);
    let depth: Vec<_> = lapper.depth().collect();
    assert_eq!(
        depth,
        vec![Interval {
            start: -2,
            stop: 3,
            val: 1
        }]
    );
}

#[test]
fn direct_structural_growth_cannot_extend_unchecked_indexing() {
    let mut lapper = Lapper::new(vec![Interval {
        start: 1_u32,
        stop: 3,
        val: 0_usize,
    }]);
    lapper.intervals.push(Interval {
        start: 4,
        stop: 6,
        val: 1,
    });
    let found: Vec<_> = lapper.find(0, 10).map(|interval| interval.val).collect();
    assert_eq!(found, vec![0]);
}
