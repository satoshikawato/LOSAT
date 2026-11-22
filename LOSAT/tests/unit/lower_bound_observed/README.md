# The order in which `std::lower_bound` calls its comparator (observed)

BLASTN's `qcovus` column (`uniq_seq_coverage` in `LOSAT/src/algorithm/blastn/tabular.rs`)
follows NCBI's `GetUniqSeqCoverage`, which builds `CRangeCollection<TSeqPos>` objects with
`x_CombineWith` and `x_Subtract` (`c++/include/util/range_coll.hpp:321-376`). Both call
`std::lower_bound` with `PRangeLessPos` (a range is less than a position when
`R.GetToOpen() <= Pos`). A query gap at the first or last subject letter of a part gives a query
range whose `TSeqPos` arithmetic wrapped, and the collection is then not in the order of
`GetToOpen()`. For such a collection the position `lower_bound` returns depends on which
elements it compares and in which order, so `lower_bound_by` (`LOSAT/src/utils/lower_bound.rs`)
asks about the same elements in the same order as `std::lower_bound` did here. The same search
stands for NCBI's other calls of `lower_bound` whose elements nothing checks for order: the units
of a binary windowmasker file (`CSeqMaskerUsetArray::get_info`, `LOSAT/src/utils/winmask/istat.rs`)
and the `indexes` of a Seq-table's sparse index (`CSeqTable_sparse_index::GetIndexAt`,
`LOSAT/src/algorithm/pssm_search/annot_index/table.rs`).

The order was observed, not taken from any library's source: `observe_lower_bound.cpp` passes
`std::lower_bound` a comparator that records every call it receives (the index of the element and
the answer it gave) and the position returned. Nothing else was looked at: no library header or
source, no preprocessed output, no disassembly or debugger trace of the library.

## Files

- `observe_lower_bound.cpp`: the observation program (modes `tree`, `data`, `one`; see its header).
- `probe_trees.txt`: the fixture of the unit test `lower_bound_asks_in_the_observed_order`
  (`LOSAT/src/utils/lower_bound_tests.rs`). One
  line per window: `tree <n> <start> <total> <preorder tree>`, the window being the `n` elements
  from index `start` of a vector of `total` elements. The tree lists every sequence of answers
  the comparator can give: `p<index>` is a call with the element at `index` (followed by the tree
  after the answer false, then the tree after the answer true), `r<position>` the position
  returned (an index of the vector). Windows: 0 to 64 elements and 127, 128, 129, 255, 256, 257,
  511, 512, 513 elements from index 0; 0 to 20 elements from indices 1, 7 and 100 (137 trees, 5,535
  answer sequences).

## How it was produced (2026-11-20)

Compiler: `g++-13 (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0` (the `.comment` section of NCBI's
BLAST+ 2.17.0 Linux `blastn` names GCC 13.2.0).

```sh
g++-13 -O2 -std=c++17 -o observe_lower_bound observe_lower_bound.cpp
( for n in $(seq 0 64) 127 128 129 255 256 257 511 512 513; do
    ./observe_lower_bound tree $n 0 $n | grep '^tree'; done
  for s in 1 7 100; do for n in $(seq 0 20); do
    ./observe_lower_bound tree $n $s $((s + n + 3)) | grep '^tree'; done; done ) > probe_trees.txt
```

The unordered collections of `lower_bound_of_unordered_ranges_is_the_observed_one` are
`./observe_lower_bound one <start> <key> <to_open...>`:

```text
one 0 45 40 50 0 60          -> calls 2:1 3:0, returns 3
one 0 25 10 30 0 40 50       -> calls 2:1 4:0 3:0, returns 3
one 1 25 10 30 0 40 50       -> calls 3:0 2:1, returns 3
one 0 50 100 0 0 0 0 0 0 0   -> calls 4:1 6:1 7:1, returns 8
one 1 100 50 0 120 4 90 200  -> calls 3:1 5:0 4:1, returns 5
```

## What was observed

The rule, in words: the search keeps a window of `n` elements from `lo`, at first the whole
range. While the window is not empty, it calls the comparator with the element `n / 2` places into
the window (rounded down). When the comparator answers true (the element is less), the window
becomes the elements after that one; otherwise it becomes the elements before it. The search
returns `lo` of the empty window. Every window of `n` elements had exactly `n` distinct calls in
its tree and `n + 1` returned positions; no call repeated an element within one search or fell
outside the window.

The larger observation behind the rule, kept outside the repository, matched it with no
exception: every answer sequence for 0 to 1,100 elements from index 0 and for 0 to 300 elements
from indices 1, 2, 3, 7, 100 and 1,000 (2,907 trees, 876,450 calls, 879,357 returned positions),
and 2,200,001 calls of `std::lower_bound` with `PRangeLessPos`'s answer on random collections
(in order, with ties, out of order, and built by `x_CombineWith` and `x_Subtract` with wrapped
ranges; up to 5,000 elements; 1,013,350 from a start above 0; 675,941 out of order for their key).
An independent review (2026-11-21, the same compiler and options) observed NCBI's two other
forms of the call: no comparator on a `vector<unsigned>` with a `size_t` value (`GetIndexAt`), and
`operator<` on `const entry *` pointers (`get_info`). They followed the same rule: 80,088 answer
sequences (0 to 140 elements from indices 0, 1, 7 and 100) and 21,570,705 searches (every vector
of 0 to 7 values from 0 to 7, every value 0 to 8; 12,881,916 out of order for their value), 0
mismatches.
