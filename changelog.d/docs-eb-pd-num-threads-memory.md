- Documentation: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 1.27 records in exception 2 that with more than 64 threads a
  search in which NCBI BLAST+ 2.17.0's memory does not grow with `-num_threads` (`-subject`, or a `tblastn` database of
  one sequence) can use more memory in LOSAT: up to 1.21 times NCBI's peak memory at 96 threads and 1.44 times at 128,
  while LOSAT finished 8 to 69 times sooner (the Owner's choice of 2026-11-22). The manual's differences and
  performance pages say so.
