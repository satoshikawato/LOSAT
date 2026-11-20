- BLASTP, TBLASTN, TBLASTX and BLASTN now write outfmt 11 with the environment variable `PRE_FETCH_SEQS_LIMIT` set to
  a value that NCBI BLAST+ cannot convert to an int (an empty value too), and BLASTP, PSIBLAST, TBLASTN and TBLASTX
  now search `-db` with `BL2SEQ_LEGACY` set, as NCBI BLAST+ 2.17.0 does: its archive is written without reading the
  value, and a database search ignores `BL2SEQ_LEGACY`. LOSAT refused both before. The value is still refused in
  the other formats and by PSIBLAST, which NCBI stops before its first report (exit status 255), and
  `BL2SEQ_LEGACY` is still refused with `-subject`, where NCBI searches each subject on its own with another
  report.
