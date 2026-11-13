- `-db` with the single-file XML2 and JSON documents (outfmt 16 and 15), in every search program that writes
  them: the taxonomy database's messages of a query's `Hit`s (a damaged `taxdb.bti` or `taxdb.btd`) now come
  after the document's start and, in outfmt 15, after the `,` before the query's report, where NCBI BLAST+
  2.17.0 writes them when standard output and standard error go to one file. LOSAT wrote them before the
  document's start and the earlier queries' reports (standard output and standard error are each unchanged).
