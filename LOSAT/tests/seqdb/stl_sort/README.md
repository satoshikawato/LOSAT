# libstdc++ `std::sort` permutations

`cases.txt`: the order in which libstdc++'s `std::sort` (and `std::partial_sort` over the whole
range, its fallback) leaves ranges `from:to_open` under NCBI's `CRange_Less` (`CASE 0`) and
`CRange_ReverseLess` (`CASE 1`), for the test `writedb::mask_data::range_set::tests` (`stl_sort`).
Many of the ranges are empty and differ, which those comparators treat as equivalent, so the
permutation is the library's own. Made with g++ 13.3.0 (the algorithm is the one of GCC 8.5, which
built makeblastdb 2.17.0):

```
g++ -O2 -o sort_cases sort_cases.cpp
python3 antiqsort.py | ./sort_cases > cases.txt
```

`sort_cases.cpp` prints each case before (`CASE <reverse> <n> ...`) and after the sort (`SORTED ...`,
`HEAP ...`): first six whose ranges it draws with `std::mt19937(7)` (2026-11-13; n = 40-305, too
shallow for the fallback), then the two inputs that `antiqsort.py` writes (2026-11-14), which reach
introsort's depth limit: McIlroy's adversary ("A Killer Adversary for Quicksort", 1999), mirrored,
against a model of `__introsort_loop`, stopped at the limit. The elements still undecided there are
only known to be below the decided ones, so they become empty ranges (every third) and small
non-empty ones: the heap sort of the fallback orders 97 ranges of the first (n = 120, 32 of them
empty) and 123 of the second (n = 150, 41 empty), and the order it leaves the empty ones in differs
from the one without the fallback.
