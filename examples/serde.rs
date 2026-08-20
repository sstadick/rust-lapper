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

    let encoded = bincode::serialize(&lapper).unwrap();
    let decoded: Lapper<usize, u32> = bincode::deserialize(&encoded[..]).unwrap();
    assert_eq!(decoded.intervals, lapper.intervals);
    assert_eq!(decoded.cov(), lapper.cov());
    assert_eq!(decoded.count(11, 15), lapper.count(11, 15));
}
