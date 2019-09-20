## 9/20/19
- Maybe add max_ends list to the Lapper. Then binary search for upper
  bounds like AIList does, and iterate backwards checking against the
  max_end to deteremine overlap... It could make the hopper actually
  work as well since binary searches in a local area would be useful

## 9/3/19
- Added a worst case mode to allow for setting max misses and then
  recalculating the binary search from there. I think that for anything
  except maybe whole chromosome intervals, it's faster to not use it.
- Updated docs
## 9/2/19
- Corrected the range on the IterFind to properly stop looking once stop
  has been hit.
- Added benchmarks for find and seek
- Sped up the intesection method when we know that overlaps have been
  merged.
- Figure out the possible speedup with miss counting!!
- Fix up docs to include features list / make less block of text
## 9/1/19
- Corrected the intersection_and_union function to work regardless of
  the comparison direction. It's a little slower now and in need of
  benchmarking. I think there is a gain to be had if we know that both
  lappers have had overlaps merged, I'm just now sure how to make it
  work yet.
- Copy the version of intersect_and_union from chromcomp, it works and
  is faster
- Next big todo: add benchmarking!

## 8/31/19
- Moved a bunch of functionality in, like merge overlaps and detect
  intersects and coverage calculation.
- I would like to add some whole Lapper functions like a union an
  intersect. -- did this... but needs fine tuning
- I would like to add benchmarks and identify the slow stuff
- Something is now really slow I think. Like maybe the coverage
  calculation. Bench it and see. I also think the checked_sub might be
  slowing things down. -- figured it out. HashSet was a dumb idea for
  calculation of coverage and has been changed to be faster. 
