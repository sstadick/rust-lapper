//! This module provides a simple data-structure for fast interval searches.
//! ## Features
//! - Extremely fast on most genomic datasets. (3-4x faster than other methods)
//! - Extremely fast on in order queries. (10x faster than other methods)
//! - Extremely fast intersections count method based on the
//! [BITS](https://arxiv.org/pdf/1208.3407.pdf) algorithm
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
//! Most interaction with this crate will be through the [`Lapper`](struct.Lapper.html) struct
//! The main methods are [`find`](struct.Lapper.html#method.find),
//! [`seek`](struct.Lapper.html#method.seek), and [`count`](struct.Lapper.html#method.count)
//! where both `seek` and `count` are special cases allowing for very fast queries in certain scenarios.
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
//! `find` or `seek`. The `count` method will be fast in all scenarios.
//!
//! # Examples
//!
//! ```rust
//!    use rust_lapper::{Interval, Lapper};
//!    use std::cmp;
//!    type Iv = Interval<usize, u32>;
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
use num_traits::{
    identities::{one, zero},
    PrimInt, Unsigned,
};
use std::cmp::Ordering;

#[cfg(feature = "with_serde")]
use serde::{Deserialize, Serialize};

/// Represent a range from [start, stop)
/// Inclusive start, exclusive of stop (start <= x < end)
#[cfg_attr(feature = "with_serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone)]
pub struct Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    pub start: I,
    pub stop: I,
    pub val: T,
}

/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[cfg_attr(feature = "with_serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone)]
pub struct Lapper<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    /// List of intervals
    pub intervals: Vec<Interval<I, T>>,
    /// Sorted list of start positions,
    starts: Vec<I>,
    /// Sorted list of end positions,
    stops: Vec<I>,
    /// The length of the longest interval
    max_len: I,
    /// The calculated number of positions covered by the intervals
    cov: Option<I>,
    /// Whether or not overlaps have been merged
    pub overlaps_merged: bool,
}

impl<I, T> Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    /// Compute the intsect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Interval<I, T>) -> I {
        std::cmp::min(self.stop, other.stop)
            .checked_sub(std::cmp::max(&self.start, &other.start))
            .unwrap_or_else(zero::<I>)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: I, stop: I) -> bool {
        self.start < stop && self.stop > start
    }
}

impl<I, T> PartialEq for Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    #[inline]
    fn eq(&self, other: &Interval<I, T>) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}

impl<I, T> PartialOrd for Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(match self.start.cmp(&other.start) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.stop.cmp(&other.stop),
        })
    }
}

impl<I, T> Lapper<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Clone + Send + Sync,
{
    /// Create a new instance of Lapper by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// ```
    pub fn new(mut intervals: Vec<Interval<I, T>>) -> Self {
        intervals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let (mut starts, mut stops): (Vec<_>, Vec<_>) =
            intervals.iter().map(|x| (x.start, x.stop)).unzip();
        starts.sort();
        stops.sort();
        let mut max_len = zero::<I>();
        for interval in intervals.iter() {
            let i_len = interval
                .stop
                .checked_sub(&interval.start)
                .unwrap_or_else(zero::<I>);
            if i_len > max_len {
                max_len = i_len;
            }
        }
        Lapper {
            intervals,
            starts,
            stops,
            max_len,
            cov: None,
            overlaps_merged: false,
        }
    }

    /// Insert a new interval after the Lapper has been created. This is very
    /// inefficient and should be avoided if possible.
    ///
    /// SIDE EFFECTS: This clears cov() and overlaps_merged
    /// meaning that those will have to be recomputed after a insert
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data : Vec<Interval<usize, usize>>= vec!{
    ///     Interval{start:0,  stop:5,  val:1},
    ///     Interval{start:6,  stop:10, val:2},
    /// };
    /// let mut lapper = Lapper::new(data);
    /// lapper.insert(Interval{start:0, stop:20, val:5});
    /// assert_eq!(lapper.len(), 3);
    /// assert_eq!(lapper.find(1,3).collect::<Vec<&Interval<usize,usize>>>(),
    ///     vec![
    ///         &Interval{start:0, stop:5, val:1},
    ///         &Interval{start:0, stop:20, val:5},
    ///     ]
    /// );
    ///
    /// ```
    pub fn insert(&mut self, elem: Interval<I, T>) {
        let starts_insert_index = Self::bsearch_seq(elem.start, &self.starts);
        let stops_insert_index = Self::bsearch_seq(elem.stop, &self.stops);
        let intervals_insert_index = Self::bsearch_seq_ref(&elem, &self.intervals);
        let i_len = elem.stop.checked_sub(&elem.start).unwrap_or_else(zero::<I>);
        if i_len > self.max_len {
            self.max_len = i_len;
        }
        self.starts.insert(starts_insert_index, elem.start);
        self.stops.insert(stops_insert_index, elem.stop);
        self.intervals.insert(intervals_insert_index, elem);
        self.cov = None;
        self.overlaps_merged = false;
    }

    /// Get the number over intervals in Lapper
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.len(), 4);
    /// ```
    #[inline]
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if lapper is empty
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data: Vec<Interval<usize, bool>> = vec![];
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.is_empty(), true);
    /// ```
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Get the number of positions covered by the intervals in Lapper. This provides immutable
    /// access if it has already been set, or on the fly calculation.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.cov(), 25);
    #[inline]
    pub fn cov(&self) -> I {
        match self.cov {
            None => self.calculate_coverage(),
            Some(cov) => cov,
        }
    }

    /// Get the number fo positions covered by the intervals in Lapper and store it. If you are
    /// going to be using the coverage, you should set it to avoid calculating it over and over.
    pub fn set_cov(&mut self) -> I {
        let cov = self.calculate_coverage();
        self.cov = Some(cov);
        cov
    }

    /// Calculate the actual coverage behind the scenes.
    fn calculate_coverage(&self) -> I {
        let mut moving_interval = Interval {
            start: zero::<I>(),
            stop: zero::<I>(),
            val: zero::<I>(),
        };
        let mut cov = zero::<I>();

        for interval in self.intervals.iter() {
            // If it overlaps, embrace, extend, extinguish
            if moving_interval.overlap(interval.start, interval.stop) {
                moving_interval.start = std::cmp::min(moving_interval.start, interval.start);
                moving_interval.stop = std::cmp::max(moving_interval.stop, interval.stop);
            } else {
                // add the set and move on
                cov = cov + (moving_interval.stop - moving_interval.start);
                moving_interval.start = interval.start;
                moving_interval.stop = interval.stop;
            }
        }
        // add in the last bit
        cov = cov + (moving_interval.stop - moving_interval.start);
        cov
    }

    /// Return an iterator over the intervals in Lapper
    #[inline]
    pub fn iter(&self) -> IterLapper<I, T> {
        IterLapper {
            inner: self,
            pos: 0,
        }
    }

    /// Merge any intervals that overlap with eachother within the Lapper. This is an easy way to
    /// speed up queries.
    pub fn merge_overlaps(&mut self) {
        let mut merged: Vec<Interval<I, T>> = Vec::new();

        for interval in &self.intervals {
            match merged.last_mut() {
                // If there is an overlap; extend the last interval to cover the new interval
                Some(last) if last.stop > interval.start => {
                    last.stop = std::cmp::max(last.stop, interval.stop);
                }
                // No overlap; add the new interval as is
                _ => merged.push(interval.clone()),
            }
        }

        self.intervals = merged;
        self.update_auxiliary_structures();
        self.overlaps_merged = true;
    }

    /// Processes overlapping intervals by splitting and merging them based on a custom merge function.
    pub fn merge_overlaps_with<F>(&mut self, merge_fn: F)
    where
        F: Fn(&T, &T) -> T,
    {
        let mut merged: Vec<Interval<I, T>> = Vec::new();

        for interval in &self.intervals {
            match merged.last_mut() {
                // If there is an overlap; extend the last interval to cover the new interval
                Some(last) if last.stop > interval.start => {
                    last.stop = std::cmp::max(last.stop, interval.stop);
                    last.val = merge_fn(&last.val, &interval.val);
                }
                // No overlap; add the new interval as is
                _ => merged.push(interval.clone()),
            }
        }

        self.intervals = merged;
        self.update_auxiliary_structures();
        self.overlaps_merged = true;
    }

    /// Divides a set of overlapping intervals into non-overlapping intervals,
    /// aggregating associated data for each resulting interval using a custom merge function.
    // Based on: https://stackoverflow.com/questions/628837/how-to-divide-a-set-of-overlapping-ranges-into-non-overlapping-ranges
    pub fn divide_overlaps_with<F>(&mut self, merge_fn: F)
    where
        F: Fn(&[&T], std::ops::Range<I>) -> T,
    {
        // Create start and end events for each interval
        let mut events: Vec<(I, bool, I, usize)> = Vec::new();
        for (index, interval) in self.intervals.iter().enumerate() {
            events.push((interval.start, true, interval.stop, index)); // Start event
            events.push((interval.stop, false, interval.start, index)); // End event
        }

        events.sort_by(|a, b| {
            a.0.cmp(&b.0) // First, sort by endpoint
                .then_with(|| (!a.1 as u8).cmp(&(!b.1 as u8))) // Then, start events before end events
                .then_with(|| a.2.cmp(&b.2)) // Finally, sort by the other endpoint if needed
        });

        let mut active_indices: Vec<usize> = Vec::new();
        let mut ranges: Vec<Interval<I, T>> = Vec::new();
        let mut current_start: Option<I> = None;

        for (endpoint, is_start, _, index) in events {
            // Handle the start of an interval
            if is_start {
                if let Some(start) = current_start {
                    // Merge and push the interval if it doesn't overlap directly with its predecessor
                    if endpoint > start && active_indices.len() > 0 {
                        let values = active_indices
                            .iter()
                            .map(|&i| &self.intervals[i].val)
                            .collect::<Vec<_>>();
                        ranges.push(Interval {
                            start,
                            stop: endpoint,
                            val: merge_fn(&values, start..endpoint),
                        });
                    }
                }

                // Add index to active intervals
                active_indices.push(index);
            }
            // Handle the end of an interval
            else {
                // Create an interval up to the current endpoint
                if let Some(start) = current_start {
                    if endpoint > start && active_indices.len() > 0 {
                        let values = active_indices
                            .iter()
                            .map(|&i| &self.intervals[i].val)
                            .collect::<Vec<_>>();
                        ranges.push(Interval {
                            start,
                            stop: endpoint,
                            val: merge_fn(&values, start..endpoint),
                        });
                    }
                }

                // Remove ended interval
                active_indices.retain(|&i| i != index);
            }

            // Update the start for a new or continued interval
            current_start = Some(endpoint);
        }

        self.intervals = ranges;
        self.update_auxiliary_structures();
        self.overlaps_merged = true;
    }

    /// Helper method to update starts, stops, and max_len based on the current state of intervals.
    fn update_auxiliary_structures(&mut self) {
        let (mut starts, mut stops): (Vec<_>, Vec<_>) =
            self.intervals.iter().map(|iv| (iv.start, iv.stop)).unzip();
        starts.sort();
        stops.sort();
        self.starts = starts;
        self.stops = stops;
        self.max_len = self
            .intervals
            .iter()
            .map(|iv| iv.stop.checked_sub(&iv.start).unwrap_or_else(zero::<I>))
            .max()
            .unwrap_or_else(zero::<I>);
    }

    /// Determine the first index that we should start checking for overlaps for via a binary
    /// search.
    /// Assumes that the maximum interval length in `intervals` has been subtracted from
    /// `start`, otherwise the result is undefined
    #[inline]
    pub fn lower_bound(start: I, intervals: &[Interval<I, T>]) -> usize {
        let mut size = intervals.len();
        let mut low = 0;

        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = low + half;
            let other_low = low + other_half;
            let v = &intervals[probe];
            size = half;
            low = if v.start < start { other_low } else { low }
        }
        low
    }

    #[inline]
    pub fn bsearch_seq<K>(key: K, elems: &[K]) -> usize
    where
        K: PartialEq + PartialOrd,
    {
        Self::bsearch_seq_ref(&key, elems)
    }

    #[inline]
    pub fn bsearch_seq_ref<K>(key: &K, elems: &[K]) -> usize
    where
        K: PartialEq + PartialOrd,
    {
        if elems.is_empty() {
            return 0;
        }
        if elems[0] > *key {
            return 0;
        }
        let mut high = elems.len();
        let mut low = 0;

        while high - low > 1 {
            let mid = (high + low) / 2;
            if elems[mid] < *key {
                low = mid;
            } else {
                high = mid;
            }
        }
        high
    }

    /// Find the union and the intersect of two lapper objects.
    /// Union: The set of positions found in both lappers
    /// Intersect: The number of positions where both lappers intersect. Note that a position only
    /// counts one time, multiple Intervals covering the same position don't add up.
    /// ``` rust
    /// use rust_lapper::{Lapper, Interval};
    /// type Iv = Interval<u32, u32>;
    /// let data1: Vec<Iv> = vec![
    ///     Iv{start: 70, stop: 120, val: 0}, // max_len = 50
    ///     Iv{start: 10, stop: 15, val: 0}, // exact overlap
    ///     Iv{start: 12, stop: 15, val: 0}, // inner overlap
    ///     Iv{start: 14, stop: 16, val: 0}, // overlap end
    ///     Iv{start: 68, stop: 71, val: 0}, // overlap start
    /// ];
    /// let data2: Vec<Iv> = vec![
    ///
    ///     Iv{start: 10, stop: 15, val: 0},
    ///     Iv{start: 40, stop: 45, val: 0},
    ///     Iv{start: 50, stop: 55, val: 0},
    ///     Iv{start: 60, stop: 65, val: 0},
    ///     Iv{start: 70, stop: 75, val: 0},
    /// ];
    ///
    /// let (mut lapper1, mut lapper2) = (Lapper::new(data1), Lapper::new(data2)) ;
    /// // Should be the same either way it's calculated
    /// let (union, intersect) = lapper1.union_and_intersect(&lapper2);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// let (union, intersect) = lapper2.union_and_intersect(&lapper1);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// lapper1.merge_overlaps();
    /// lapper1.set_cov();
    /// lapper2.merge_overlaps();
    /// lapper2.set_cov();
    ///
    /// // Should be the same either way it's calculated
    /// let (union, intersect) = lapper1.union_and_intersect(&lapper2);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// let (union, intersect) = lapper2.union_and_intersect(&lapper1);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// ```
    #[inline]
    pub fn union_and_intersect(&self, other: &Self) -> (I, I) {
        let mut cursor: usize = 0;

        if !self.overlaps_merged || !other.overlaps_merged {
            let mut intersections: Vec<Interval<I, bool>> = vec![];
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
            let mut intersect = zero::<I>();
            for c1_iv in self.iter() {
                for c2_iv in other.seek(c1_iv.start, c1_iv.stop, &mut cursor) {
                    let local_intersect = c1_iv.intersect(c2_iv);
                    intersect = intersect + local_intersect;
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
    pub fn intersect(&self, other: &Self) -> I {
        self.union_and_intersect(other).1
    }

    /// Find the union of two lapper objects.
    #[inline]
    pub fn union(&self, other: &Self) -> I {
        self.union_and_intersect(other).0
    }

    /// Return the contiguous intervals of coverage, `val` represents the number of intervals
    /// covering the returned interval.
    ///
    /// # Examples
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.depth().collect::<Vec<Interval<usize, usize>>>(), vec![
    ///             Interval { start: 0, stop: 5, val: 1 },
    ///             Interval { start: 5, stop: 20, val: 2 },
    ///             Interval { start: 20, stop: 25, val: 1 }]);
    /// ```
    #[inline]
    pub fn depth(&self) -> IterDepth<I, T> {
        let mut merged_lapper = Lapper::new(
            self.intervals
                .iter()
                .map(|i| Interval {
                    start: i.start,
                    stop: i.stop,
                    val: true,
                })
                .collect::<Vec<Interval<I, bool>>>(),
        );
        merged_lapper.merge_overlaps();
        let merged_len = merged_lapper.intervals.len();
        IterDepth {
            inner: self,
            merged: merged_lapper,
            curr_merged_pos: zero::<I>(),
            curr_pos: 0,
            cursor: 0,
            end: merged_len,
        }
    }

    /// Count all intervals that overlap start .. stop. This performs two binary search in order to
    /// find all the excluded elements, and then deduces the intersection from there. See
    /// [BITS](https://arxiv.org/pdf/1208.3407.pdf) for more details.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(lapper.count(5, 11), 2);
    /// ```
    #[inline]
    pub fn count(&self, start: I, stop: I) -> usize {
        let len = self.intervals.len();
        let mut first = Self::bsearch_seq(start, &self.stops);
        let last = Self::bsearch_seq(stop, &self.starts);
        //println!("{}/{}", start, stop);
        //println!("pre start found in stops: {}: {}", first, self.stops[first]);
        //println!("pre stop found in starts: {}", last);
        //while last < len && self.starts[last] == stop {
        //last += 1;
        //}
        while first < len && self.stops[first] == start {
            first += 1;
        }
        let num_cant_after = len - last;
        len - first - num_cant_after
        //println!("{:#?}", self.starts);
        //println!("{:#?}", self.stops);
        //println!("start found in stops: {}", first);
        //println!("stop found in starts: {}", last);
    }

    /// Find all intervals that overlap start .. stop
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(lapper.find(5, 11).count(), 2);
    /// ```
    #[inline]
    pub fn find(&self, start: I, stop: I) -> IterFind<I, T> {
        IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            ),
            start,
            stop,
        }
    }

    /// Find all intevals that overlap start .. stop. This method will work when queries
    /// to this lapper are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the lapper object mutable, and thus use the same lapper accross threads.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// let mut cursor = 0;
    /// for i in lapper.iter() {
    ///    assert_eq!(lapper.seek(i.start, i.stop, &mut cursor).count(), 1);
    /// }
    /// ```
    #[inline]
    pub fn seek<'a>(&'a self, start: I, stop: I, cursor: &mut usize) -> IterFind<'a, I, T> {
        if *cursor == 0 || (*cursor < self.intervals.len() && self.intervals[*cursor].start > start)
        {
            *cursor = Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            );
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start
                < start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>)
        {
            *cursor += 1;
        }

        IterFind {
            inner: self,
            off: *cursor,
            start,
            stop,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    off: usize,
    start: I,
    stop: I,
}

impl<'a, I, T> Iterator for IterFind<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    #[inline]
    // interval.start < stop && interval.stop > start
    fn next(&mut self) -> Option<Self::Item> {
        while self.off < self.inner.intervals.len() {
            //let mut generator = self.inner.intervals[self.off..].iter();
            //while let Some(interval) = generator.next() {
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
pub struct IterDepth<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    merged: Lapper<I, bool>, // A lapper that is the merged_lapper of inner
    curr_merged_pos: I,      // Current start position in current interval
    curr_pos: usize,         // In merged list of non-overlapping intervals
    cursor: usize,           // cursor for seek over inner lapper
    end: usize,              // len of merged
}

impl<'a, I, T> Iterator for IterDepth<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = Interval<I, I>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let mut interval: &Interval<I, bool> = &self.merged.intervals[self.curr_pos];
        if self.curr_merged_pos == zero::<I>() {
            self.curr_merged_pos = interval.start;
        }
        if interval.stop == self.curr_merged_pos {
            if self.curr_pos + 1 != self.end {
                self.curr_pos += 1;
                interval = &self.merged.intervals[self.curr_pos];
                self.curr_merged_pos = interval.start;
            } else {
                return None;
            }
        }
        let start = self.curr_merged_pos;
        let depth_at_point = self
            .inner
            .seek(
                self.curr_merged_pos,
                self.curr_merged_pos + one::<I>(),
                &mut self.cursor,
            )
            .count();
        let mut new_depth_at_point = depth_at_point;
        while new_depth_at_point == depth_at_point && self.curr_merged_pos < interval.stop {
            self.curr_merged_pos = self.curr_merged_pos + one::<I>();
            new_depth_at_point = self
                .inner
                .seek(
                    self.curr_merged_pos,
                    self.curr_merged_pos + one::<I>(),
                    &mut self.cursor,
                )
                .count();
        }
        Some(Interval {
            start,
            stop: self.curr_merged_pos,
            val: I::from(depth_at_point).unwrap(), // from usize should always work
        })
    }
}
/// Lapper Iterator
pub struct IterLapper<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    pos: usize,
}

impl<'a, I, T> Iterator for IterLapper<'a, I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            self.inner.intervals.get(self.pos - 1)
        }
    }
}

impl<I, T> IntoIterator for Lapper<I, T>
where
    T: Clone + Send + Sync,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = Interval<I, T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, I, T> IntoIterator for &'a Lapper<I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;
    type IntoIter = std::slice::Iter<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<I, T>> {
        self.intervals.iter()
    }
}

impl<'a, I, T> IntoIterator for &'a mut Lapper<I, T>
where
    T: Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a mut Interval<I, T>;
    type IntoIter = std::slice::IterMut<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<I, T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval<usize, u32>;
    fn setup_nonoverlapping() -> Lapper<usize, u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                stop: x + 10,
                val: 0,
            })
            .collect();
        Lapper::new(data)
    }

    fn setup_overlapping() -> Lapper<usize, u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                stop: x + 15,
                val: 0,
            })
            .collect();
        Lapper::new(data)
    }

    fn setup_badlapper() -> Lapper<usize, u32> {
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
        Lapper::new(data)
    }

    fn setup_single() -> Lapper<usize, u32> {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            stop: 35,
            val: 0,
        }];
        Lapper::new(data)
    }

    // Test that inserting data ends up with the same lapper (nonoverlapping)
    #[test]
    fn insert_equality_nonoverlapping() {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                stop: x + 10,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.stops, insert_lapper.stops);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (overlapping)
    #[test]
    fn insert_equality_overlapping() {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                stop: x + 15,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.stops, insert_lapper.stops);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data half with new and half with insert
    // ends up with the same lapper
    #[test]
    fn insert_equality_half_and_half() {
        let data: Vec<Iv> = (0..100)
            .step_by(1)
            .map(|x| Iv {
                start: x,
                stop: x + 15,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let (new_data, insert_data) = data.split_at(50);
        let mut insert_lapper = Lapper::new(new_data.to_vec());
        let mut insert_data = insert_data.to_vec();
        insert_data.reverse();
        for elem in insert_data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.stops, insert_lapper.stops);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (badlapper)
    #[test]
    fn insert_equality_badlapper() {
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
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.stops, insert_lapper.stops);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (single)
    #[test]
    fn insert_equality_single() {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            stop: 35,
            val: 0,
        }];
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.stops, insert_lapper.stops);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that a query stop that hits an interval start returns no interval
    #[test]
    fn test_query_stop_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        assert_eq!(None, lapper.seek(15, 20, &mut cursor).next());
        assert_eq!(lapper.find(15, 20).count(), lapper.count(15, 20));
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        assert_eq!(None, lapper.seek(30, 35, &mut cursor).next());
        assert_eq!(lapper.find(30, 35).count(), lapper.count(30, 35));
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
        assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next());
        assert_eq!(lapper.find(15, 25).count(), lapper.count(15, 25));
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
        assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next());
        assert_eq!(lapper.find(25, 35).count(), lapper.count(25, 35));
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
        assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next());
        assert_eq!(lapper.find(22, 27).count(), lapper.count(22, 27));
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
        assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next());
        assert_eq!(lapper.find(15, 35).count(), lapper.count(15, 35));
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
        assert_eq!(lapper.count(8, 20), 2);
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
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.merge_overlaps();
        assert_eq!(expected, lapper.iter().collect::<Vec<&Iv>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len())
    }

    #[test]
    fn test_merge_overlaps_with() {
        let mut lapper = setup_badlapper();
        let expected: Vec<&Iv> = vec![
            &Iv{ start: 10, stop: 16, val: 3 }, // 3 overlaps, initial val = 0, +3 overlaps
            &Iv{ start: 40, stop: 45, val: 0 }, // No overlap, val remains 0
            &Iv{ start: 50, stop: 55, val: 0 }, // No overlap, val remains 0
            &Iv{ start: 60, stop: 65, val: 0 }, // No overlap, val remains 0
            &Iv{ start: 68, stop: 120, val: 2 }, // 2 overlaps, initial val = 0, +2 overlaps
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
         lapper.merge_overlaps_with( |a: &u32, _b: &u32| -> u32 { *a + 1 });
        assert_eq!(expected, lapper.iter().collect::<Vec<&Iv>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
    }

    #[test]
    fn test_divide_overlaps_with() {
        let mut lapper : Lapper<u32, String> = Lapper::new(vec![
            Interval { start: 1, stop: 5, val: String::from("a") },
            Interval { start: 3, stop: 7, val: String::from("b") },
            Interval { start: 6, stop: 9, val: String::from("c") },
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 1, stop: 3, val: String::from("a") },
            Interval { start: 3, stop: 5, val: String::from("a, b") },
            Interval { start: 5, stop: 6, val: String::from("b") },
            Interval { start: 6, stop: 7, val: String::from("b, c") },
            Interval { start: 7, stop: 9, val: String::from("c") },
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len()); 
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len()); 
    }

    // Testcase from:
    // https://stackoverflow.com/questions/628837/how-to-divide-a-set-of-overlapping-ranges-into-non-overlapping-ranges
    #[test]
    fn test_divide_overlaps_with_stackoverflow() {
        let mut lapper: Lapper<u32, String> = Lapper::new(vec![
            Interval { start: 0, stop: 100, val: String::from("a") },
            Interval { start: 0, stop: 75, val: String::from("b") },
            Interval { start: 75, stop: 80, val: String::from("d") },
            Interval {start: 95, stop: 150, val: String::from("c")},
            Interval {start: 120, stop: 130, val: String::from("d")},
            Interval {start: 160, stop: 175, val: String::from("e")},
            Interval {start: 165, stop: 180, val: String::from("a")},
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 0, stop: 75, val: String::from("b, a") },
            Interval { start: 75, stop: 80, val: String::from("a, d") },
            Interval {start: 80, stop: 95, val: String::from("a")},
            Interval {start: 95, stop: 100, val: String::from("a, c")},
            Interval {start: 100, stop: 120, val: String::from("c")},
            Interval {start: 120, stop: 130, val: String::from("c, d")},
            Interval {start: 130, stop: 150, val: String::from("c")},
            Interval {start: 160, stop: 165, val: String::from("e")},
            Interval {start: 165, stop: 175, val: String::from("e, a")},
            Interval {start: 175, stop: 180, val: String::from("a")},
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len()); 
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len()); 
    }

    #[test]
    fn test_divide_overlaps_with_contained_interval() {
        let mut lapper = Lapper::new(vec![
            Interval { start: 1, stop: 10, val: String::from("a") },
            Interval { start: 3, stop: 7, val: String::from("b") },
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 1, stop: 3, val: String::from("a") },
            Interval { start: 3, stop: 7, val: String::from("a, b") },
            Interval { start: 7, stop: 10, val: String::from("a") },
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len()); 
    }

    #[test]
    fn test_divide_overlaps_with_non_overlapping_intervals() {
        let mut lapper = Lapper::new(vec![
            Interval { start: 1, stop: 2, val: String::from("a") },
            Interval { start: 3, stop: 4, val: String::from("b") },
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 1, stop: 2, val: String::from("a") },
            Interval { start: 3, stop: 4, val: String::from("b") },
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
    }

    #[test]
    fn test_divide_overlaps_with_partial_overlap() {
        let mut lapper = Lapper::new(vec![
            Interval { start: 1, stop: 4, val: String::from("a") },
            Interval { start: 3, stop: 6, val: String::from("b") },
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 1, stop: 3, val: String::from("a") },
            Interval { start: 3, stop: 4, val: String::from("a, b") },
            Interval { start: 4, stop: 6, val: String::from("b") },
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
    }

    #[test]
    fn test_divide_overlaps_with_exact_overlap() {
        let mut lapper = Lapper::new(vec![
            Interval { start: 1, stop: 4, val: String::from("a") },
            Interval { start: 1, stop: 4, val: String::from("b") },
            Interval { start: 3, stop: 6, val: String::from("c") },
        ]);
        let expected: Vec<Interval<u32, String>> = vec![
            Interval { start: 1, stop: 3, val: String::from("a, b") },
            Interval { start: 3, stop: 4, val: String::from("a, b, c") },
            Interval { start: 4, stop: 6, val: String::from("c") },
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.divide_overlaps_with(|overlap, _| overlap
            .iter()
            .fold(String::new(), |acc, x| acc + x + ", ")
            .trim_end_matches(", ").to_string());
        assert_eq!(expected, lapper.iter().cloned().collect::<Vec<_>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
    }

    // This test was added because this breakage was found in a library user's code, where after
    // calling merge_overlaps(), the find() call returned an empty iterator.
    #[test]
    fn test_merge_overlaps_find() {
        let data = vec![
                Iv{start: 2, stop: 3, val: 0},
                Iv{start: 3, stop: 4, val: 0},
                Iv{start: 4, stop: 6, val: 0},
                Iv{start: 6, stop: 7, val: 0},
                Iv{start: 7, stop: 8, val: 0},
        ];
        let mut lapper = Lapper::new(data);

        let found = lapper.find(7, 9).collect::<Vec<&Interval<_,_>>>();
        assert_eq!(found, vec![
            &Iv{start:7, stop: 8, val: 0},
        ]);

        // merge_overlaps should merge all intervals to one, which should be returned in the find call.
        lapper.merge_overlaps();

        let found = lapper.find(7, 9).collect::<Vec<&Interval<_,_>>>();
        assert_eq!(found, vec![
            &Iv{start:2, stop: 8, val: 0},
        ]);
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
        assert_eq!(lapper.count(8, 11), 3);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 100, stop: 200, val: 0},
            &Iv{start: 111, stop: 160, val: 0},
            &Iv{start: 150, stop: 200, val: 0},
        ]);

        assert_eq!(lapper.count(145, 151), 3);
    }

    #[test]
    fn test_depth_sanity() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, stop: 10, val: 0},
            Iv{start: 5, stop: 10, val: 0}
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
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
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
        assert_eq!(found, vec![
                   Interval{start: 1, stop: 2, val: 1},
                   Interval{start: 2, stop: 3, val: 2},
                   Interval{start: 3, stop: 8, val: 5},
                   Interval{start: 8, stop: 9, val: 1},
                   Interval{start: 9, stop: 10, val: 2},
                   Interval{start: 10, stop: 11, val: 1},
        ]);
    }
    #[test]
    fn test_depth_harder() {
        let data1: Vec<Iv> = vec![
            Iv{start: 1, stop: 10, val: 0},
            Iv{start: 2, stop: 5, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 5, stop: 8, val: 0},
            Iv{start: 9, stop: 11, val: 0},
            Iv{start: 15, stop: 20, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
        assert_eq!(found, vec![
                   Interval{start: 1, stop: 2, val: 1},
                   Interval{start: 2, stop: 3, val: 2},
                   Interval{start: 3, stop: 8, val: 5},
                   Interval{start: 8, stop: 9, val: 1},
                   Interval{start: 9, stop: 10, val: 2},
                   Interval{start: 10, stop: 11, val: 1},
                   Interval{start: 15, stop: 20, val: 1},
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
        assert_eq!(lapper.find(50, 55).count(), lapper.count(50,55));
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
        ]);
        assert_eq!(lapper.count(28974798, 33141355), 1);
    }

    // #[test]
    // fn serde_test() {
    //     let data = vec![
    //         Iv{start:25264912, stop: 25264986, val: 0},
    //         Iv{start:27273024, stop: 27273065	, val: 0},
    //         Iv{start:27440273, stop: 27440318	, val: 0},
    //         Iv{start:27488033, stop: 27488125	, val: 0},
    //         Iv{start:27938410, stop: 27938470	, val: 0},
    //         Iv{start:27959118, stop: 27959171	, val: 0},
    //         Iv{start:28866309, stop: 33141404	, val: 0},
    //     ];
    //     let lapper = Lapper::new(data);

    //     let serialized = bincode::serialize(&lapper).unwrap();
    //     let deserialzed: Lapper<usize, u32> = bincode::deserialize(&serialized).unwrap();

    //     let found = deserialzed.find(28974798, 33141355).collect::<Vec<&Iv>>();
    //     assert_eq!(found, vec![
    //         &Iv{start:28866309, stop: 33141404	, val: 0},
    //     ]);
    //     assert_eq!(deserialzed.count(28974798, 33141355), 1);
    // }

}
