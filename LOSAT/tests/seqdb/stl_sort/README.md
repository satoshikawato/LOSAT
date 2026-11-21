# Observed orders of NCBI's range sort

`makeblastdb -mask_data` combines the Seq-locs of one Seq-id with `CSeq_loc::Add`, which sorts the
ranges of each Seq-id and strand with C++ `sort` under `CRange_Less` or `CRange_ReverseLess`
(`c++/src/objects/seqloc/Seq_loc.cpp:4283-4322`, `4565-4572`). `cases.txt` holds orders that a C++
standard library leaves such ranges in, written by `sort_cases.cpp`:

- `CASE <reverse> <n> <from:to_open> ...`: the input, `reverse` 1 for `CRange_ReverseLess`;
- `SORTED <from:to_open> ...`: the order `std::sort` left them in;
- `HEAP <from:to_open> ...`: the order `std::partial_sort` over all of them left them in.

`sort_cases.cpp` prints six random cases (`std::mt19937` seed 7), then one case for each `CASE`
line of its standard input; the last two cases of `cases.txt` came from such lines. Built with
g++ 13.3.0 (`g++-13 -O2 sort_cases.cpp`) and given the last two `CASE` lines of `cases.txt`, it
prints `cases.txt` byte for byte (2026-11-20).

The ranges that neither comparator orders are two whole ranges, two equal ranges, or two empty
ranges (`to_open <= from`); only the empty ones can differ, and the cases hold many different
ones. The test `the_ranges_take_the_observed_orders`
(`LOSAT/src/writedb/mask_data/range_set_tests.rs`) checks that LOSAT's sort, `sort_as_observed`,
leaves every case in the `SORTED` order and that `heap_order` leaves it in the `HEAP` order. How
that sort was written from observed comparator calls: `LOSAT/tests/unit/sort_observed/README.md`.
