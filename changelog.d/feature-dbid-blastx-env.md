- BLASTX now refuses the environment variable `PRE_FETCH_SEQS_LIMIT` with a value that NCBI BLAST+ cannot convert
  to an int (an empty value too), with `-db` and `-subject` and in every format but outfmt 11, and `BL2SEQ_LEGACY`
  with `-subject`, as BLASTP, TBLASTN and TBLASTX refuse them. NCBI BLAST+ 2.17.0 BLASTX stops before its first
  report with exit status 255 in the first case and searches each subject on its own, with another report, in the
  second; LOSAT wrote its usual report with exit status 0. Outfmt 11, whose archive NCBI writes without reading the
  value, and `BL2SEQ_LEGACY` with `-db`, which NCBI ignores, are written as before.
