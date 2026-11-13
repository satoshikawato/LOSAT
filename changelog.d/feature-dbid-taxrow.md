- `-db` with `PRE_FETCH_SEQS_LIMIT` turning a batch's prefetch off, in every search program: the taxonomy
  database's messages of a report that reads its subjects one by one (a damaged `taxdb.bti` or `taxdb.btd`)
  now come where NCBI BLAST+ reads each subject: before the tabular row that reads it (outfmt 6, 7, 10), its
  description line or, for a subject the table does not show, its alignment's heading (outfmt 0), and its
  `Hit` (outfmt 5), as NCBI BLAST+ 2.17.0 writes them where standard output and standard error go to one
  file. LOSAT wrote them before the query's report (standard output is unchanged). A row that shows a plural
  taxonomy name opens the taxonomy database where its subject's read did not, as NCBI BLAST+'s does.
