//! This module provides a simple data-structure for fast interval searches.
//! ## Features
//! - Extremely fast on most genomic datasets. (3-4x faster than other methods)
//! - Extremely fast on in order queries. (10x faster than other methods)
//! - Parallel friendly. Queries are on an immutable structure, even for seek
//! - Consumer / Adapter paradigm, Iterators are returned and serve as the main API for interacting
//! with the lapper
//!
//! ## Details:
//!
//! ```text
//!       0  1  2  3  4  5  6  7  8  9  10 11
//! (0,10]X  X  X  X  X  X  X  X  X  X
//! (2,5]       X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (5,9]                X  X  X  X
//! (8,11]                        X  X  X
//!
//! Query: (8, 11]
//! Answer: ((0,10], (5,9], (8,11])
//! ```
//!
//! The main methods are `find` and `seek` where the latter uses a cursor and is very fast for
//! cases when the queries are sorted. This is another innovation in this library that allows an
//! additional ~50% speed improvement when consecutive queries are known to be in sort order.
//!
//! The overlap function for this assumes a zero based genomic coordinate system. So [start, stop)
//! is not inclusive of the stop position for neither the queries, nor the Intervals.
//!
//! Lapper does not use an interval tree, instead, it operates on the assumtion that most intervals are
//! of similar length; or, more exactly, that the longest interval in the set is not long compred to
//! the average distance between intervals.
//!
//! For cases where this holds true (as it often does with genomic data), we can sort by start and
//! use binary search on the starts, accounting for the length of the longest interval. The advantage
//! of this approach is simplicity of implementation and speed. In realistic tests queries returning
//! the overlapping intervals are 1000 times faster than brute force and queries that merely check
//! for the overlaps are > 5000 times faster.
//!
//! When this is not the case, if possible in your scenario, use merge_overlaps first, and then use
//! find or seek.
//!
//! # Examples
//!
//! ```rust
//!    use rust_lapper::{Interval, Lapper};
//!    use std::cmp;
//!    type Iv = Interval<u32>;
//!
//!    // create some fake data
//!    let data: Vec<Iv> = (0..20).step_by(5).map(|x| Iv{start: x, stop: x + 2, val: 0}).collect();
//!    println!("{:#?}", data);
//!
//!    // make lapper structure
//!    let laps = Lapper::new(data);
//!
//!    assert_eq!(laps.find(6, 11).next(), Some(&Iv{start: 5, stop: 7, val: 0}));
//!    
//!    // Demonstration of seek function. By passing in the &mut cursor, seek can have thread local
//!    // cursors going
//!    let mut sim: usize = 0;
//!    let mut cursor = 0;
//!    // Calculate the overlap between the query and the found intervals, sum total overlap
//!    for i in (0..10).step_by(3) {
//!        sim += laps
//!            .seek(i, i + 2, &mut cursor)
//!            .map(|iv| cmp::min(i + 2, iv.stop) - cmp::max(i, iv.start))
//!            .sum::<usize>();
//!    }
//!    assert_eq!(sim, 4);
//! ```
//#![feature(asm)]
use std::cmp::Ordering::{self};
use std::collections::VecDeque;

/// Represent a range from [start, stop)
/// Inclusive start, exclusive of stop
#[derive(Eq, Debug, Clone)]
pub struct Interval<T: Eq + Clone> {
    pub start: usize,
    pub stop: usize,
    pub val: T,
}

/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[derive(Debug)]
pub struct Lapper<T: Eq + Clone> {
    /// List of intervasl
    pub intervals: Vec<Interval<T>>,
    /// The length of the longest interval
    max_len: usize,
    /// A cursor to hold the position in the list in between searches with `seek` method
    cursor: usize,
    /// The calculated number of positions covered by the intervals
    cov: Option<usize>,
    /// Whether or not overlaps have been merged
    pub overlaps_merged: bool,
}

impl<T: Eq + Clone> Interval<T> {
    /// Compute the intsect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Interval<T>) -> usize {
        std::cmp::min(self.stop, other.stop)
            .checked_sub(std::cmp::max(self.start, other.start))
            .unwrap_or(0)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: usize, stop: usize) -> bool {
        self.start < stop && self.stop > start
    }
}

impl<T: Eq + Clone> Ord for Interval<T> {
    #[inline]
    fn cmp(&self, other: &Interval<T>) -> Ordering {
        if self.start < other.start {
            Ordering::Less
        } else if other.start < self.start {
            Ordering::Greater
        } else {
            self.stop.cmp(&other.stop)
        }
    }
}

impl<T: Eq + Clone> PartialOrd for Interval<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl<T: Eq + Clone> PartialEq for Interval<T> {
    #[inline]
    fn eq(&self, other: &Interval<T>) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}

impl<T: Eq + Clone> Lapper<T> {
    /// Create a new instance of Lapper by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    pub fn new(mut intervals: Vec<Interval<T>>) -> Self {
        intervals.sort();
        let mut max_len = 0;
        for interval in intervals.iter() {
            let i_len = interval.stop.checked_sub(interval.start).unwrap_or(0);
            if i_len > max_len {
                max_len = i_len;
            }
        }
        Lapper {
            intervals,
            max_len,
            cursor: 0,
            cov: None,
            overlaps_merged: false,
        }
    }

    /// Get the number over intervals in Lapper
    #[inline]
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if lapper is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Get the number of positions covered by the intervals in Lapper. This provides immutable
    /// access if it has already been set, or on the fly calculation.
    #[inline]
    pub fn cov(&self) -> usize {
        match self.cov {
            None => self.calculate_coverage(),
            Some(cov) => cov,
        }
    }

    /// Get the number fo positions covered by the intervals in Lapper and store it. If you are
    /// going to be using the coverage, you should set it to avoid calculating it over and over.
    pub fn set_cov(&mut self) -> usize {
        let cov = self.calculate_coverage();
        self.cov = Some(cov);
        cov
    }

    /// Calculate teh actual coverage behind the scenes.
    fn calculate_coverage(&self) -> usize {
        let mut moving_interval = Interval {
            start: 0,
            stop: 0,
            val: 0,
        };
        let mut cov = 0;

        for interval in self.intervals.iter() {
            // If it overlaps, embrace, extend, extinguish
            if moving_interval.overlap(interval.start, interval.stop) {
                moving_interval.start = std::cmp::min(moving_interval.start, interval.start);
                moving_interval.stop = std::cmp::max(moving_interval.stop, interval.stop);
            } else {
                // add the set and move on
                cov += moving_interval.stop - moving_interval.start;
                moving_interval.start = interval.start;
                moving_interval.stop = interval.stop;
            }
        }
        // add in the last bit
        cov += moving_interval.stop - moving_interval.start;
        cov
    }

    /// Return an iterator over the intervals in Lapper
    #[inline]
    pub fn iter(&self) -> IterLapper<T> {
        IterLapper {
            inner: self,
            pos: 0,
        }
    }

    /// Merge any intervals that overlap with eachother within the Lapper. This is an easy way to
    /// speed up queries.
    pub fn merge_overlaps(&mut self) {
        let mut stack: VecDeque<&mut Interval<T>> = VecDeque::new();
        let mut ivs = self.intervals.iter_mut();
        if let Some(first) = ivs.next() {
            stack.push_back(first);
            for interval in ivs {
                let mut top = stack.pop_back().unwrap();
                if top.stop < interval.start {
                    stack.push_back(top);
                    stack.push_back(interval);
                } else if top.stop < interval.stop {
                    top.stop = interval.stop;
                    //stack.pop_back();
                    stack.push_back(top);
                } else {
                    // they were equal
                    stack.push_back(top);
                }
            }
            self.overlaps_merged = true;
            self.intervals = stack
                .into_iter()
                .map(|x| Interval {
                    start: x.start,
                    stop: x.stop,
                    val: x.val.clone(),
                })
                .collect();
        }
    }

    /// Determine the first index that we should start checking for overlaps for via a binary
    /// search.
    #[inline]
    fn lower_bound(&self, start: usize) -> usize {
        let mut size = self.intervals.len();
        let mut low = 0;

        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = low + half;
            let other_low = low + other_half;
            let v = &self.intervals[probe];
            size = half;
            low = if v.start < start { other_low } else { low }
        }
        low
    }

    /// Find the union and the intersect of two lapper objects.
    /// Union: The set of positions found in both lappers
    /// Intersect: The number of positions where both lappers intersect. Note that a position only
    /// counts one time, multiple Intervals covering the same position don't add up.
    #[inline]
    pub fn union_and_intersect(&self, other: &Self) -> (usize, usize) {
        let mut cursor: usize = 0;

        if !self.overlaps_merged || !other.overlaps_merged {
            let mut intersections: Vec<Interval<bool>> = vec![];
            for self_iv in self.iter() {
                for other_iv in other.seek(self_iv.start, self_iv.stop, &mut cursor) {
                    let start = std::cmp::max(self_iv.start, other_iv.start);
                    let stop = std::cmp::min(self_iv.stop, other_iv.stop);
                    intersections.push(Interval {
                        start,
                        stop,
                        val: true,
                    });
                }
            }
            let mut temp_lapper = Lapper::new(intersections);
            temp_lapper.merge_overlaps();
            temp_lapper.set_cov();
            let union = self.cov() + other.cov() - temp_lapper.cov();
            (union, temp_lapper.cov())
        } else {
            let mut intersect = 0;
            for c1_iv in self.iter() {
                for c2_iv in other.seek(c1_iv.start, c1_iv.stop, &mut cursor) {
                    let local_intersect = c1_iv.intersect(c2_iv);
                    intersect += local_intersect;
                }
            }
            let union = self.cov() + other.cov() - intersect;
            (union, intersect)
        }
    }

    /// Find the intersect of two lapper objects.
    /// Intersect: The number of positions where both lappers intersect. Note that a position only
    /// counts one time, multiple Intervals covering the same position don't add up
    #[inline]
    pub fn intersect(&self, other: &Self) -> usize {
        self.union_and_intersect(other).1
    }

    /// Find the union of two lapper objects.
    #[inline]
    pub fn union(&self, other: &Self) -> usize {
        self.union_and_intersect(other).0
    }

    #[inline]
    pub fn depth(&self) -> IterDepth<T> {
        let mut merged_lapper = Lapper::new(
            self.intervals
                .iter()
                .map(|i| Interval {
                    start: i.start,
                    stop: i.stop,
                    val: true,
                })
                .collect::<Vec<Interval<bool>>>(),
        );
        merged_lapper.merge_overlaps();
        let merged_len = merged_lapper.intervals.len();
        IterDepth {
            inner: self,
            merged: merged_lapper,
            curr_merged_pos: 0,
            curr_pos: 0,
            cursor: 0,
            end: merged_len,
        }
    }

    /// Find all intervals that overlap start .. stop
    #[inline]
    pub fn find(&self, start: usize, stop: usize) -> IterFind<T> {
        IterFind {
            inner: self,
            off: self.lower_bound(start.checked_sub(self.max_len).unwrap_or(0)),
            end: self.intervals.len(),
            start,
            stop,
        }
    }

    /// Find all intevals that overlap start .. stop. This method will work when queries
    /// to this lapper are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the lapper object mutable, and thus use the same lapper accross threads.
    #[inline]
    pub fn seek<'a>(&'a self, start: usize, stop: usize, cursor: &mut usize) -> IterFind<'a, T> {
        if *cursor == 0 || (*cursor < self.intervals.len() && self.intervals[*cursor].start > start)
        {
            *cursor = self.lower_bound(start.checked_sub(self.max_len).unwrap_or(0));
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start < start.checked_sub(self.max_len).unwrap_or(0)
        {
            *cursor += 1;
        }

        IterFind {
            inner: self,
            off: *cursor,
            end: self.intervals.len(),
            start,
            stop,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, T>
where
    T: Eq + Clone + 'a,
{
    inner: &'a Lapper<T>,
    off: usize,
    end: usize,
    start: usize,
    stop: usize,
}

impl<'a, T: Eq + Clone> Iterator for IterFind<'a, T> {
    type Item = &'a Interval<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.off < self.inner.intervals.len() {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.stop) {
                return Some(interval);
            } else if interval.start >= self.stop {
                break;
            }
        }
        None
    }
}

/// Depth Iterator
#[derive(Debug)]
pub struct IterDepth<'a, T>
where
    T: Eq + Clone + 'a,
{
    inner: &'a Lapper<T>,
    merged: Lapper<bool>,   // A lapper that is the merged_lapper of inner
    curr_merged_pos: usize, // Current start position in current interval
    curr_pos: usize,        // In merged list of non-overlapping intervals
    cursor: usize,          // cursor for seek over inner lapper
    end: usize,             // len of merged
}

impl<'a, T: Eq + Clone> Iterator for IterDepth<'a, T> {
    type Item = Interval<u32>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let mut interval: &Interval<bool> = &self.merged.intervals[self.curr_pos];
        if self.curr_merged_pos == 0 {
            self.curr_merged_pos = interval.start;
        }
        if interval.stop == self.curr_merged_pos {
            if self.curr_pos + 1 != self.end {
                self.curr_pos += 1;
                self.curr_merged_pos = interval.start;
                interval = &self.merged.intervals[self.curr_pos];
            } else {
                return None;
            }
        }
        let start = self.curr_merged_pos;
        let depth_at_point = self
            .inner
            .seek(
                self.curr_merged_pos,
                self.curr_merged_pos + 1,
                &mut self.cursor,
            )
            .count();
        let mut new_depth_at_point = depth_at_point;
        while new_depth_at_point == depth_at_point && self.curr_merged_pos != interval.stop {
            self.curr_merged_pos += 1;
            new_depth_at_point = self
                .inner
                .seek(
                    self.curr_merged_pos,
                    self.curr_merged_pos + 1,
                    &mut self.cursor,
                )
                .count();
        }
        Some(Interval {
            start,
            stop: self.curr_merged_pos,
            val: depth_at_point as u32,
        })
    }
}
/// Lapper Iterator
pub struct IterLapper<'a, T>
where
    T: Eq + Clone + 'a,
{
    inner: &'a Lapper<T>,
    pos: usize,
}

impl<'a, T: Eq + Clone> Iterator for IterLapper<'a, T> {
    type Item = &'a Interval<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            self.inner.intervals.get(self.pos - 1)
        }
    }
}

impl<T: Eq + Clone> IntoIterator for Lapper<T> {
    type Item = Interval<T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, T: Eq + Clone> IntoIterator for &'a Lapper<T> {
    type Item = &'a Interval<T>;
    type IntoIter = std::slice::Iter<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<T>> {
        self.intervals.iter()
    }
}

impl<'a, T: Eq + Clone> IntoIterator for &'a mut Lapper<T> {
    type Item = &'a mut Interval<T>;
    type IntoIter = std::slice::IterMut<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval<u32>;
    fn setup_nonoverlapping() -> Lapper<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                stop: x + 10,
                val: 0,
            })
            .collect();
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_overlapping() -> Lapper<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                stop: x + 15,
                val: 0,
            })
            .collect();
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_badlapper() -> Lapper<u32> {
        let data: Vec<Iv> = vec![
            Iv{start: 70, stop: 120, val: 0}, // max_len = 50
            Iv{start: 10, stop: 15, val: 0},
            Iv{start: 10, stop: 15, val: 0}, // exact overlap
            Iv{start: 12, stop: 15, val: 0}, // inner overlap
            Iv{start: 14, stop: 16, val: 0}, // overlap end
            Iv{start: 40, stop: 45, val: 0},
            Iv{start: 50, stop: 55, val: 0},
            Iv{start: 60, stop: 65, val: 0},
            Iv{start: 68, stop: 71, val: 0}, // overlap start
            Iv{start: 70, stop: 75, val: 0},
        ];
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_single() -> Lapper<u32> {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            stop: 35,
            val: 0,
        }];
        let lapper = Lapper::new(data);
        lapper
    }

    // Test that a query stop that hits an interval start returns no interval
    #[test]
    fn test_query_stop_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        assert_eq!(None, lapper.seek(15, 20, &mut cursor).next())
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        assert_eq!(None, lapper.seek(30, 35, &mut cursor).next())
    }

    // Test that a query that overlaps the start of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 25).next());
        assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next())
    }

    // Test that a query that overlaps the stop of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(25, 35).next());
        assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next())
    }

    // Test that a query that is enveloped by interval returns interval
    #[test]
    fn test_interval_envelops_query() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(22, 27).next());
        assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next())
    }

    // Test that a query that envolops an interval returns that interval
    #[test]
    fn test_query_envolops_interval() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 35).next());
        assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next())
    }

    #[test]
    fn test_overlapping_intervals() {
        let lapper = setup_overlapping();
        let mut cursor = 0;
        let e1 = Iv {
            start: 0,
            stop: 15,
            val: 0,
        };
        let e2 = Iv {
            start: 10,
            stop: 25,
            val: 0,
        };
        assert_eq!(vec![&e1, &e2], lapper.find(8, 20).collect::<Vec<&Iv>>());
        assert_eq!(
            vec![&e1, &e2],
            lapper.seek(8, 20, &mut cursor).collect::<Vec<&Iv>>()
        );
    }

    #[test]
    fn test_merge_overlaps() {
        let mut lapper = setup_badlapper();
        let expected: Vec<&Iv> = vec![
            &Iv{start: 10, stop: 16, val: 0},
            &Iv{start: 40, stop: 45, val: 0},
            &Iv{start: 50, stop: 55, val: 0},
            &Iv{start: 60, stop: 65, val: 0},
            &Iv{start: 68, stop: 120, val: 0}, // max_len = 50
        ];
        lapper.merge_overlaps();
        assert_eq!(expected, lapper.iter().collect::<Vec<&Iv>>())
        
    }

    #[test]
    fn test_lapper_cov() {
        let mut lapper = setup_badlapper();
        let before = lapper.cov();
        lapper.merge_overlaps();
        let after = lapper.cov();
        assert_eq!(before, after);

        let mut lapper = setup_nonoverlapping();
        lapper.set_cov();
        assert_eq!(lapper.cov(), 50);
    }

    #[test]
    fn test_interval_intersects() {
        let i1 = Iv{start: 70, stop: 120, val: 0}; // max_len = 50
        let i2 = Iv{start: 10, stop: 15, val: 0};
        let i3 = Iv{start: 10, stop: 15, val: 0}; // exact overlap
        let i4 = Iv{start: 12, stop: 15, val: 0}; // inner overlap
        let i5 = Iv{start: 14, stop: 16, val: 0}; // overlap end
        let i6 = Iv{start: 40, stop: 50, val: 0};
        let i7 = Iv{start: 50, stop: 55, val: 0};
        let i_8 = Iv{start: 60, stop: 65, val: 0};
        let i9 = Iv{start: 68, stop: 71, val: 0}; // overlap start
        let i10 = Iv{start: 70, stop: 75, val: 0};

        assert_eq!(i2.intersect(&i3), 5); // exact match
        assert_eq!(i2.intersect(&i4), 3); // inner intersect
        assert_eq!(i2.intersect(&i5), 1); // end intersect
        assert_eq!(i9.intersect(&i10), 1); // start intersect
        assert_eq!(i7.intersect(&i_8), 0); // no intersect
        assert_eq!(i6.intersect(&i7), 0); // no intersect stop = start
        assert_eq!(i1.intersect(&i10), 5); // inner intersect at start
    }

    #[test]
    fn test_union_and_intersect() {
        let data1: Vec<Iv> = vec![
            Iv{start: 70, stop: 120, val: 0}, // max_len = 50
            Iv{start: 10, stop: 15, val: 0}, // exact overlap
            Iv{start: 12, stop: 15, val: 0}, // inner overlap
            Iv{start: 14, stop: 16, val: 0}, // overlap end
            Iv{start: 68, stop: 71, val: 0}, // overlap start
        ];
        let data2: Vec<Iv> = vec![

            Iv{start: 10, stop: 15, val: 0},
            Iv{start: 40, stop: 45, val: 0},
            Iv{start: 50, stop: 55, val: 0},
            Iv{start: 60, stop: 65, val: 0},
            Iv{start: 70, stop: 75, val: 0},
        ];
        
        let (mut lapper1, mut lapper2) = (Lapper::new(data1), Lapper::new(data2)) ;
        // Should be the same either way it's calculated
        let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        lapper1.merge_overlaps();
        lapper1.set_cov();
        lapper2.merge_overlaps();
        lapper2.set_cov();

        // Should be the same either way it's calculated
        let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
    }

    #[test]
    fn test_find_overlaps_in_large_intervals() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, stop: 8, val: 0},
            Iv{start: 1, stop: 10, val: 0}, 
            Iv{start: 2, stop: 5, val: 0}, 
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 4, stop: 7, val: 0},
            Iv{start: 5, stop: 8, val: 0},
            Iv{start: 8, stop: 8, val: 0},
            Iv{start: 9, stop: 11, val: 0},
            Iv{start: 10, stop: 13, val: 0},
            Iv{start: 100, stop: 200, val: 0},
            Iv{start: 110, stop: 120, val: 0},
            Iv{start: 110, stop: 124, val: 0},
            Iv{start: 111, stop: 160, val: 0},
            Iv{start: 150, stop: 200, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.find(8, 11).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 1, stop: 10, val: 0}, 
            &Iv{start: 9, stop: 11, val: 0},
            &Iv{start: 10, stop: 13, val: 0},
        ]);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 100, stop: 200, val: 0},
            &Iv{start: 111, stop: 160, val: 0},
            &Iv{start: 150, stop: 200, val: 0},
        ]);

    }

    #[test]
    fn test_depth_sanity() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, stop: 10, val: 0},
            Iv{start: 5, stop: 10, val: 0}
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<u32>>>();
        assert_eq!(found, vec![
                   Interval{start: 0, stop: 5, val: 1},
                   Interval{start: 5, stop: 10, val: 2}
        ]);
    }

    #[test]
    fn test_depth_hard() {
        let data1: Vec<Iv> = vec![
            Iv{start: 1, stop: 10, val: 0},
            Iv{start: 2, stop: 5, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 5, stop: 8, val: 0},
            Iv{start: 9, stop: 11, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<u32>>>();
        assert_eq!(found, vec![
                   Interval{start: 1, stop: 2, val: 1},
                   Interval{start: 2, stop: 3, val: 2},
                   Interval{start: 3, stop: 8, val: 5},
                   Interval{start: 8, stop: 9, val: 1},
                   Interval{start: 9, stop: 10, val: 2},
                   Interval{start: 10, stop: 11, val: 1},
        ]);
    }
    // BUG TESTS - these are tests that came from real life

    // Test that it's not possible to induce index out of bounds by pushing the cursor past the end
    // of the lapper.
    #[test]
    fn test_seek_over_len() {
        let lapper = setup_nonoverlapping();
        let single = setup_single();
        let mut cursor: usize = 0;

        for interval in lapper.iter() {
            for o_interval in single.seek(interval.start, interval.stop, &mut cursor) {
                println!("{:#?}", o_interval);
            }
        }
    }

    // Test that if lower_bound puts us before the first match, we still return a match
    #[test]
    fn test_find_over_behind_first_match() {
        let lapper = setup_badlapper();
        let e1 = Iv {start: 50, stop: 55, val: 0};
        let found = lapper.find(50, 55).next();
        assert_eq!(found, Some(&e1));
    }

    // When there is a very long interval that spans many little intervals, test that the little
    // intevals still get returne properly
    #[test]
    fn test_bad_skips() {
        let data = vec![
            Iv{start:25264912, stop: 25264986, val: 0},	
            Iv{start:27273024, stop: 27273065	, val: 0},
            Iv{start:27440273, stop: 27440318	, val: 0},
            Iv{start:27488033, stop: 27488125	, val: 0},
            Iv{start:27938410, stop: 27938470	, val: 0},
            Iv{start:27959118, stop: 27959171	, val: 0},
            Iv{start:28866309, stop: 33141404	, val: 0},
        ];
        let lapper = Lapper::new(data);

        let found = lapper.find(28974798, 33141355).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start:28866309, stop: 33141404	, val: 0},
        ])

    }
}
