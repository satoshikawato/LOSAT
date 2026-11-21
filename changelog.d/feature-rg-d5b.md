- `blastn -use_index true` follows NCBI BLAST+ 2.17.0's loader of the MegaBLAST index: `-index_name`
  (old-style indices) is accepted, and every case where NCBI uses no index gives NCBI's bytes (its
  error lines about index volumes it skips, its warning, then the search without the index: no
  index found, a first index volume whose header does not read, a `-word_size` below the index's
  minimum, a task other than megablast, `-subject`).
  A search that NCBI runs with the index is still refused.
