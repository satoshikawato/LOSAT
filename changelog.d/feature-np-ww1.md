- LOSAT Web: BLASTN's `-window_masker_db` (WindowMasker). "Filters and Masking" gains NCBI's
  "Species-specific repeats" mark with a file field for the WindowMasker counts (NCBI's
  `wmasker.obinary` or a `windowmasker -mk_counts` obinary / binary file); the search masks the
  queries as the command line does, with the same output and the same errors. The counts are read
  once and kept for the next BLASTN searches of the session; the run details show the file's name
  and SHA-256. ABI v2 gains four exports for named input files (`file_begin`, `file_chunk`,
  `file_end`, `file_release`; still version 2).
