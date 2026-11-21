- makeblastdb `-mask_data`: the ranges that a combination of masks sorts (`CSeq_loc::Add`) are put
  in order by LOSAT's own code, written from the comparator calls and the orders that the C++
  `sort` of g++ 13 was observed to make (`LOSAT/tests/unit/sort_observed/`), in place of a port of
  the C++ library's sort (AGENTS.md rule 4: no GPL-mixed code). The masks written are unchanged and
  stay NCBI's, also where one Seq-id's ranges hold different empty intervals.
