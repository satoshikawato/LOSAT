- `blastn`, `tblastx` and `tblastn` (the latter with `-searchsp` or `OLD_FSC`) write NCBI BLAST+ 2.17.0's report (apart from the version line, which names LOSAT)
  of a run whose queries have no letters to search over a database whose statistics count 0
  sequences (an alias file's `NSEQ 0`, or its `LENGTH` or `STATS_TOTLEN` over volumes without
  sequences); LOSAT refused the database before reading the queries. A run with a query to search
  is still refused, where NCBI BLAST+ divides by the count of 0 and crashes, now in the setup of
  the first query batch NCBI searches, as `blastp`, `blastx`, `psiblast` and the RPS programs do.
