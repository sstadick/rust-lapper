# Always-on 32-interval block contender

This branch is an exact archive of the live `feat/index` working-tree source at
the start of the SIMD experiments. The source hash was:

```text
ece99b9751f9aecbcdd6758415690072901ca4934d16dacfd703fcf7b7eef0ee  src/lib.rs
```

It stores maximum ends, next-greater links, and monotonic prefix maxima for
fixed 32-interval blocks. The prefix maxima choose the first candidate block;
the block links skip later all-miss blocks. Both optimizations are always on,
and `find` still yields references in ascending start order.

Recorded nine-repeat medians from the same direct harness used for the SIMD
comparison:

| Case | Build | Query | Total |
|---|---:|---:|---:|
| `1-2` | 3.263 ms | 5.875 ms | 9.124 ms |
| `7-3` | 26.742 ms | 87.010 ms | 113.716 ms |
| `8-7` | 36.949 ms | 608.013 ms | 644.990 ms |

The archived mutation caveat still applies: block metadata must be rebuilt
after `insert` and `merge_overlaps` before production integration.
