- Database searches: a sequence with a definition line without a Seq-id beside lines with Seq-ids (a
  hand-edited database; makeblastdb writes none) is reported as NCBI BLAST+ 2.17.0 reports it instead of
  the run being refused: `Unknown` in the description line where NCBI's walk of the lines meets that line
  (up to 4096 bytes), the skipped id in place of each alignment's heading, outfmt 5's `Unknown` hit where
  that line has a title after the first titled line, the XML2 and JSON descriptions up to that line, and
  `Unknown`, 0 and `N/A` in the tabular fields. LOSAT refuses `sallacc` of such a sequence (NCBI stops with
  a null pointer, exit 255), and such a sequence in the organism report (outfmt 18) and in `rpsblast` and
  `rpstblastn`.
- Database searches: where LOSAT refuses a tabular row at which NCBI BLAST+ stops or prints an uninitialised
  number (`sallacc` above; `sacc`, `saccver` or `staxid` of a sequence whose lookup fails), the reports of
  the query batches before it are written first, as NCBI writes them, and then the warnings that NCBI
  writes of the failing batch before that row (its reader messages and the warnings of its queries up to
  the row's); `blastp`, `psiblast`, `tblastn` and `tblastx` wrote no report before the refusal, and no
  program wrote all of those warnings.
