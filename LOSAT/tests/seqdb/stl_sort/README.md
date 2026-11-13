# libstdc++ `std::sort` permutations

`cases.txt`: the order in which libstdc++'s `std::sort` (and `std::partial_sort` over the whole
range, its fallback) leaves ranges `from:to_open` under NCBI's `CRange_Less` (`CASE 0`) and
`CRange_ReverseLess` (`CASE 1`), for the test `writedb::mask_data::range_set::tests` (`stl_sort`).
Many of the ranges are empty and differ, which those comparators treat as equivalent, so the
permutation is the library's own. Made on 2026-11-13 with g++ 13.2 (the algorithm is the one of
GCC 8.5, which built makeblastdb 2.17.0) from a 40-line C++ program that draws the ranges with
`std::mt19937(7)` and prints each case before (`CASE <reverse> <n> ...`) and after the sort
(`SORTED ...`, `HEAP ...`).
