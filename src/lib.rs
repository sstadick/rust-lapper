//! (Docs from nim-lapper, by Brent Pendersen)
//! This module provides a simple data-structure for fast interval searches. It does not use an
//! interval tree, instead, it operates on the assumtion that most intervals are of similar length;
//! or, more exactly, that the longest interval in the set is not long compred to the average
//! distance between intervals. On any dataset where that is not the case, this method will not
//! perform well. For cases where this holds true (as it often does with genomic data), we can sort
//! by start and use binary search on the starts, accounting for the length of the longest
//! interval. The advantage of this approach is simplicity of implementation and speed. In
//! realistic tests queries returning the overlapping intervals are 1000 times faster than brute
//! force and queries that merely check for the overlaps are > 5000 times faster.
//!
//! The main methods are `find` and `seek` where the latter uses a cursor and is very fast for
//! cases when the queries are sorted. This is another innovation in this library that allows an
//! additional ~50% speed improvement when consecutive queries are known to be in sort order.
//!
//! The overlap function for this assumes a zero based genomic coordinate system. So [start, stop)
//! is not inclusive of the stop position for neither the queries, nor the Intervals.
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
//!    let mut sim: i32 = 0;
//!    let mut cursor = 0;
//!    // Calculate the overlap between the query and the found intervals, sum total overlap
//!    for i in (0..10).step_by(3) {
//!        sim += laps
//!            .seek(i, i + 2, &mut cursor)
//!            .map(|iv| cmp::min(i + 2, iv.stop) - cmp::max(i, iv.start))
//!            .sum::<i32>();
//!    }
//!    assert_eq!(sim, 3);
//! ```
// TODO: Add benchmarks
use std::cmp::Ordering;

#[derive(Eq, Debug)]
pub struct Interval<T: Eq> {
    pub start: i32,
    pub stop: i32,
    pub val: T,
}

/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[derive(Debug)]
pub struct Lapper<T: Eq> {
    pub intervals: Vec<Interval<T>>,
    max_len: i32,
    cursor: usize,
}

impl<T: Eq> Ord for Interval<T> {
    #[inline]
    fn cmp(&self, other: &Interval<T>) -> Ordering {
        if self.start < other.start {
            return Ordering::Less;
        } else if other.start < self.start {
            return Ordering::Greater;
        } else {
            return self.stop.cmp(&other.stop);
        }
    }
}

impl<T: Eq> PartialOrd for Interval<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl<T: Eq> PartialEq for Interval<T> {
    #[inline]
    fn eq(&self, other: &Interval<T>) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}

impl<T: Eq> Interval<T> {
    fn overlap(&self, start: i32, stop: i32) -> bool {
        self.start < stop && self.stop > start
    }
}

impl<T: Eq> Lapper<T> {
    pub fn new(mut intervals: Vec<Interval<T>>) -> Self {
        intervals.sort();
        let mut max_len = 0;
        for interval in intervals.iter() {
            if interval.stop - interval.start > max_len {
                max_len = interval.stop - interval.start;
            }
        }
        let lapper = Lapper {
            intervals,
            max_len,
            cursor: 0,
        };
        lapper
    }

    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    pub fn iter<'a>(&'a self) -> IterLapper<'a, T> {
        IterLapper {
            inner: self,
            pos: 0,
        }
    }

    fn lower_bound(&self, start: i32) -> usize {
        let mut result = 0;
        let mut count = self.intervals.len();
        let mut step: usize;
        let mut pos: usize;

        while count != 0 {
            step = count / 2;
            pos = result + step;
            if self.intervals[pos].start < start {
                result = pos + 1;
                count -= step + 1;
            } else {
                count = step;
            }
        }
        result
    }

    /// Find all intervals that overlap start .. stop
    pub fn find<'a>(&'a self, start: i32, stop: i32) -> IterFind<'a, T> {
        IterFind {
            inner: self,
            off: self.lower_bound(start - self.max_len),
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
    pub fn seek<'a>(&'a self, start: i32, stop: i32, cursor: &mut usize) -> IterFind<'a, T> {
        if *cursor == 0 || (*cursor < self.intervals.len() && self.intervals[*cursor].start > start)
        //if *cursor == 0 || self.intervals[*cursor].start > start {
        {
            *cursor = self.lower_bound(start - self.max_len);
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start < (start - self.max_len)
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
pub struct IterFind<'a, T>
where
    T: Eq + 'a,
{
    inner: &'a Lapper<T>,
    off: usize,
    end: usize,
    start: i32,
    stop: i32,
}

impl<'a, T: Eq> Iterator for IterFind<'a, T> {
    type Item = &'a Interval<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.off >= self.end {
            None
        } else {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.stop) {
                Some(interval)
            } else {
                None
            }
        }
    }
}

/// Lapper Iterator
pub struct IterLapper<'a, T>
where
    T: Eq + 'a,
{
    inner: &'a Lapper<T>,
    pos: usize,
}

impl<'a, T: Eq> Iterator for IterLapper<'a, T> {
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

impl<T: Eq> IntoIterator for Lapper<T> {
    type Item = Interval<T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, T: Eq> IntoIterator for &'a Lapper<T> {
    type Item = &'a Interval<T>;
    type IntoIter = std::slice::Iter<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<T>> {
        self.intervals.iter()
    }
}

impl<'a, T: Eq> IntoIterator for &'a mut Lapper<T> {
    type Item = &'a mut Interval<T>;
    type IntoIter = std::slice::IterMut<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
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
}
