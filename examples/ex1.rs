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
