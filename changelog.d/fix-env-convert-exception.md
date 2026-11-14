- Environment values and location parts that NCBI BLAST+ cannot read as a number now end the run
  as NCBI BLAST+ 2.17.0 ends it, instead of LOSAT's own refusal (exit status 1): a `CHUNK_SIZE`,
  `OVERLAP_CHUNK_SIZE`, `BATCH_SIZE` or `PRE_FETCH_SEQS_LIMIT` such as `abc`, `1.5`, ` 5` or
  `2147483648`, and a `-query_loc` / `-subject_loc` part such as `a-10`, give NCBI's
  `Error: NCBI C++ Exception:` report with exit status 255, at the point where NCBI reads the value
  (after the `-outfmt 0` prolog for the query split sizes and the prefetch limit, before anything for
  `BATCH_SIZE` and the location; `blastn` reads `BATCH_SIZE` after the prolog), with NCBI's error
  archive for `-outfmt 11`.
