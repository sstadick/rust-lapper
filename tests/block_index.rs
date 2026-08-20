use rust_lapper::{Interval, Lapper};

#[test]
fn block_index_matches_forward_brute_force() {
    let mut state = 0x1234_5678_u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (state >> 32) as u32
    };

    let mut intervals = Vec::new();
    intervals.push(Interval {
        start: 0,
        stop: 1_000_000,
        val: 0usize,
    });
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
