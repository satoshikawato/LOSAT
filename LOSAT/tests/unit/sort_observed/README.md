# The calls that C++ `std::sort` makes to its comparator (observed)

`makeblastdb -mask_data` combines the Seq-locs of one Seq-id with `CSeq_loc::Add`, which sorts the
ranges of each Seq-id and strand with C++ `sort` under `CRange_Less` or `CRange_ReverseLess`
(`c++/src/objects/seqloc/Seq_loc.cpp:4283-4322`, `4565-4572`). Those comparators leave whole
ranges, equal ranges and empty ranges (`to_open <= from`) unordered, and the C++ standard does not
say where `sort` leaves elements it does not order. The order of different empty ranges changes
the merge that follows, so the masks written depend on what NCBI's `sort` does with them.

`sort_as_observed` (`LOSAT/src/writedb/mask_data/range_set.rs`) hands the comparator the same
elements, in the same argument order and sequence, as the `sort` of g++ 13 did here, and leaves
the same order. `heap_order` does the same for `std::partial_sort` over a whole range, which `sort`
was seen to do on a stretch after a long chain of splits.

The calls were observed, not taken from any library's source. `observe_sort.cpp` sorts the input
positions 0..n-1 with a comparator that records every call it receives:

- the two positions it was handed (first and second argument);
- the answer it gave;
- the order of the positions after the sort.

The positions sit in the middle of a larger buffer whose other cells hold numbers outside 0..n-1,
so a read outside the range would reach the comparator and be counted. Nothing else was looked at:
no library header or source, no preprocessed output, and no disassembly or debugger trace of the
library.

## Files

- `observe_sort.cpp`: the observation program. Its modes are:
  - `tree <alg> <n>`: every answer sequence that some strict weak order gives for n elements. It
    replays each prefix and, at each new call, branches on false and true where both are
    consistent with the answers so far.
  - `keys <alg> <seed> <cases> <maxn>`: random strict weak orders with many ties (small domains,
    sorted and reversed runs, organ pipes, saws, all equal, ...).
  - `ranges <alg> <seed> <cases> <maxn>`: NCBI's two comparators on ranges shaped like masks:
    exactly two different empty ranges, one empty range repeated, none, or a mixture with whole
    ranges; ranges with shared starts and duplicates.
  - `adv <alg> <seed> <cases> <maxn>`: a comparator that decides lazily: every element starts
    undecided and above every decided one, and when two undecided elements meet one of them is
    decided as the next value. This drives `sort` into long chains of splits and so into its heap
    phase. Every other case turns the undecided elements into ties at a random call.
  - `file <alg>`: one sort per input line of keys.

  `alg` is `s` for `std::sort` and `p` for `std::partial_sort(first, last, last)`; `S` and `P` do
  the same on an element type that carries a string and a shared pointer, and give the same
  records.
- `sort_calls.txt`: the fixture of the unit test `the_sort_makes_the_observed_calls`
  (`LOSAT/src/writedb/mask_data/range_set_tests.rs`). It has one line per sort:
  `<tag> <alg> <n> <answers> <calls> <order> <oob> [<extra>]`.
  - tag: `T` tree, `K` keys, `G` ranges, `A` adversary;
  - answers: 0/1 in call order;
  - calls: `a:b` for each call, where `a` and `b` are the input positions of the first and second
    argument;
  - order: the input positions after the sort;
  - oob: calls handed a cell outside the range (0 on every line);
  - extra: the shape and keys (`K`), `<reverse>/<from:to_open,...>` (`G`: the test also sorts these
    ranges with `range_less` and checks the answers), the call at which the adversary made ties
    (`A`, -1 for none);
  - `-` stands for an empty field.

  The fixture holds every answer sequence for 0 to 5 elements for both `sort` and `partial_sort`
  (333 sorts), and 58 random sorts: 50 `sort` runs of 21 to 92 elements (keys, ranges, adversary;
  10 of them reach the heap phase) and 8 `partial_sort` runs of 6 to 62 elements (keys): 391 sorts,
  24,798 calls.

## What was seen (the rule `sort_as_observed` follows)

- Fewer than two elements: no call.
- While a stretch has more than 16 elements, it is split: first the middle one of the elements at
  its positions 1, len/2 and len-1 (three or two calls) is swapped to its front; then a left mark
  from position 1 moves right while `(element, front)` is true, a right mark from the end moves
  left one step and then while `(front, element)` is true; if the left mark is before the right
  one, the two elements are swapped, the left mark steps on and the scan goes on, else the left
  mark is the cut. The part right of the cut is handled first, then the left part.
- A chain of splits has room for `2 * floor(log2(n))` splits, n being the whole length; a stretch
  of more than 16 elements that finds none left is ordered through a heap (`heap_order`, the calls
  that `std::partial_sort` makes over that stretch).
- Last, an insertion pass over the whole range: among the first 16 positions an element is first
  compared with the front `(element, front)`, and goes to the front if that is true; otherwise, and
  always from position 16 on, it moves left while `(element, left neighbour)` is true.

Each detail is pinned by the records: changing the split threshold to 15 or 17, the room for
splits by one either way, the middle position, the order of the two parts, the choice among three
tied candidates or the child chosen in the heap makes the model disagree with them.

## How it was produced (2026-11-21)

Compiler: `g++-13 (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0`. The `.comment` section of NCBI's
BLAST+ 2.17.0 Linux `blastn` names GCC 13.2.0.

```sh
g++-13 -O2 -std=c++17 -o observe_sort observe_sort.cpp
( for n in 0 1 2 3 4 5; do ./observe_sort tree s $n; done
  for n in 0 1 2 3 4 5; do ./observe_sort tree p $n; done
  ./observe_sort keys s 20261121 400 90 | awk '$3 >= 17' | head -16
  ./observe_sort ranges s 20261122 400 90 | awk '$3 >= 17' | head -20
  ./observe_sort adv s 20261123 400 100 | awk '$3 >= 17 && !seen[$0]++' | head -14
  ./observe_sort keys p 20261124 200 60 | awk '$3 >= 6' | head -8 ) > sort_calls.txt
```

Larger record sets (not in the repository) were replayed through the same Rust functions with no
mismatch in calls, argument order, number of calls or final order:

- every answer sequence of `sort` for 0 to 10 elements: 4,037,914 sorts, 124,682,208 calls;
- every answer sequence of `partial_sort` for 0 to 9 elements: 606,357 sorts, 16,760,726 calls;
- 15 random sets of keys, ranges and adversary answers of up to 20,000 elements, for `sort` and
  `partial_sort`, two of them on the element type with a string and a shared pointer: 176,760
  sorts, 310,231,613 calls; 17,434 of the `sort` runs reached the heap phase.

That is 4,821,031 sorts and 451,674,547 calls in all. No call was handed a cell outside the range.
