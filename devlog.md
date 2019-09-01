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
