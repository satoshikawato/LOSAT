- `blastn -use_index true` follows NCBI BLAST+ 2.17.0's loader of the MegaBLAST index: `-index_name`
  (old-style indices) is accepted, and every case where NCBI uses no index gives NCBI's bytes (its
  error lines about index volumes it skips, its warning, then the search without the index: no
  index found, a first index volume whose header does not read, a `-word_size` below the index's
  minimum, a task other than megablast, `-subject`).
- `blastn -use_index true` also searches with the index where NCBI does, and gives NCBI's bytes: the
  one-hit and two-hit (`-window_size`) searches, and an index of only some database volumes, whose
  other volumes are scanned (a later query batch finds no seeds in some index volumes there, as in
  NCBI). It refuses `-db_hard_mask` with an index, and a later query batch of an index of only some
  volumes where NCBI's search writes to freed memory.
