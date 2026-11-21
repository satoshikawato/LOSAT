# The calls that `std::list::sort` makes to its comparator (observed)

`-sorthits` and `-sorthsps` follow NCBI's `CAlignFormatUtil::SortSeqalignForSortableFormat`
(`c++/src/objtools/align_format/align_format_util.cpp`). It sorts `std::list`s of hits and HSPs
with NCBI's comparators, and some of those comparators have two properties that make the call
order visible in the output:

- Some are not strict weak orders. `SortHitByTotalScoreDescending` uses `>=`.
  `SortHspByPercentIdentityDescending` and `SortHitByPercentIdentityDescendingEx` order some lists
  in a cycle when an HSP has no identities.
- Some sort both hits' HSPs each time they are called.

The order NCBI prints therefore depends on which elements its list sort compares, and in which
order. `list_sort` (`LOSAT/src/report/hit_sort.rs`) hands the comparator the same elements, in the
same argument order and sequence, as `std::list::sort` did here.

The calls were observed, not taken from any library's source. `observe_list_sort.cpp` sorts a
`std::list<int>` of the input positions 0..n-1 with a comparator that records every call it
receives:

- the two elements it was handed (first and second argument);
- the answer it gave;
- the order of the list after the sort.

Nothing else was looked at: no library header or source, no preprocessed output, and no
disassembly or debugger trace of the library.

## Files

- `observe_list_sort.cpp`: the observation program. Its modes are:
  - `tree <n>`: every answer sequence for a list of n elements. It replays each prefix and
    branches each new call on false and true.
  - `rand <seed> <cases> <maxn>`: random answers.
  - `list <seed> <permille> <n>...`: one list of each length, with random answers.
- `list_sort_calls.txt`: the fixture of the unit test `list_sort_calls_as_observed`
  (`LOSAT/src/report/hit_sort_tests.rs`). It has one line per sort:
  - `T <n> <answers> <calls> <order>`, from `tree`;
  - `R <i> <n> <permille> <answers> <calls> <order>`, from `list`.

  The fields are:
  - answers: 0/1 in call order;
  - calls: `a:b` for each call, where `a` and `b` are the input positions of the first and second
    argument;
  - order: the input positions after the sort;
  - `-` stands for an empty field.

  The fixture holds every answer sequence for lists of 0 to 5 elements (154 sorts) and random
  answers for 28 lists of 6 to 257 elements: 182 sorts, 7,563 calls.

## How it was produced (2026-11-21)

Compiler: `g++-13 (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0`. The `.comment` section of NCBI's BLAST+
2.17.0 Linux `blastn` names GCC 13.2.0.

```sh
g++-13 -O2 -std=c++17 -o observe_list_sort observe_list_sort.cpp
( for n in 0 1 2 3 4 5; do ./observe_list_sort tree $n; done
  ./observe_list_sort list 20261121 500 6 7 8 9 10 11 12 13 15 16 17 31 32 33 63 64 65 100 127 128 129 257
  ./observe_list_sort list 20261122 100 16 33 100
  ./observe_list_sort list 20261123 900 16 33 100 ) > list_sort_calls.txt
```

## What was observed

The rule, in words:

- **Short lists.** A list of fewer than two elements is left as it is, without a call.
- **Building runs.** The sort reads the elements from the front, one at a time, and keeps those it
  has read as runs. A run is a stretch of consecutive elements, merged in order. The runs'
  lengths are distinct powers of two, and the longer runs are nearer the front. Each element read
  becomes a new run of one at the back. Then, while the run before the last one has the same
  length as the last one, the two are merged.
- **Finishing.** When every element is read, the runs are merged from the back: the last run with
  the one before it, then that run with the one before it, up to the first run.
- **Merging.** Two neighbouring runs are merged from their fronts. While both have elements left,
  the comparator is called with the first element left of the run behind as its first argument
  and the first element left of the run in front as its second. The answer true takes the
  element of the run behind; false takes the element of the run in front. When one run is used
  up, the rest of the other follows without a call.

The larger observation behind the rule, kept outside the repository, matched it with no exception
in the elements handed over, their argument order, the sequence of the calls, or the final order:

- every answer sequence for lists of 0 to 9 elements (409,114 sorts, 8,190,114 calls);
- random answers (p(true) from 0.02 to 0.98) for:
  - 50,000 lists of up to 300 elements (25,904,346 calls);
  - 2,000 lists of up to 5,000 elements (27,061,211 calls);
  - 100 lists of up to 20,000 elements (7,144,069 calls).

In all: 461,214 sorts and 68,299,740 calls.

The NCBI-shaped comparators were also run through `sort_query_hits`, on lists of 0 to 300 hits
or HSPs with many ties (24,080 lists). Every one ended in the order `std::list::sort` left,
including:

- the order of each hit's HSPs after the comparators' own sorts;
- 2,947 lists that `SortHspByPercentIdentityDescending` orders in a cycle.
