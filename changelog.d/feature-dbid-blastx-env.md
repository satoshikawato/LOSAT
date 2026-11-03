- BLASTX now refuses the environment variable `PRE_FETCH_SEQS_LIMIT` with a value that NCBI BLAST+ cannot convert
  to an int (an empty value too), with `-db` and `-subject` and in every format but outfmt 11, and `BL2SEQ_LEGACY`
  with `-subject`, as BLASTP, TBLASTN and TBLASTX refuse them. NCBI BLAST+ 2.17.0 BLASTX stops before its first
  report in the first case (exit status 255; with `-mt_mode 1` and `-db` it writes nothing and exits 0) and
  searches each subject on its own, with another report, in the second; LOSAT wrote its usual report with exit
  status 0. Outfmt 11, whose archive NCBI writes without reading the value, and `BL2SEQ_LEGACY` with `-db`, which
  NCBI ignores, are written as before. LOSAT refuses before the search, so an error that NCBI raises in the first
  batch's search (a query record without letters, `-query_loc` past every record: exit status 3) is this refusal
  with these variables set.
