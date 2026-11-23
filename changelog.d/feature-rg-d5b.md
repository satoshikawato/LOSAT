- `blastn -use_index true` follows NCBI BLAST+ 2.17.0's loader of the MegaBLAST index: `-index_name`
  (old-style indices) is accepted, and every case where NCBI uses no index gives NCBI's bytes (its
  error lines about index volumes it skips, its warning, then the search without the index: no
  index found, a first index volume whose header does not read, a `-word_size` below the index's
  minimum, a task other than megablast, `-subject`).
- `blastn -use_index true` also searches with the index where NCBI does, and gives NCBI's bytes: the
  one-hit and two-hit (`-window_size`) searches, and an index of only some database volumes, whose
  other volumes are scanned (a later query batch finds no seeds in some index volumes there, as in
  NCBI's `-num_threads 1`, at every thread count). With `-db_hard_mask` the search with the index
  reports the alignments in the unmasked parts of each sequence, the same at every
  `-num_threads`, where NCBI BLAST+ takes some of their seeds from another part of the sequence or
  from the next sequence (see the manual's Differences page). It refuses an index that does not
  fit the database (where NCBI reads past the index or the sequence), and a later query batch of
  an index of only some volumes where NCBI's search writes to freed memory.
