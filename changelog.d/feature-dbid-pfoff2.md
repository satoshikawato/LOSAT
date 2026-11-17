- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: the tabular formats (`-outfmt 6`, `7` and `10`) are now written as NCBI BLAST+ writes them, in BLASTP,
  psiblast, BLASTN, BLASTX, TBLASTN and TBLASTX. A field list with an id, title or taxonomy field (the default fields
  too) writes that row with the scores, zero counts and coordinates and the ids, titles and taxonomy of the query's
  last row that read its sequence (`Unknown` and `N/A` before one), adds the Seq-align id to `sallseqid`, `sallgi` and
  `sallacc`, and goes on; a field list that reads the sequence otherwise (`slen`, the aligned letters, or a count
  computed from them) writes the rows before it and stops with `BLAST Database error: DB contains no accession info.`
  (exit status 2). No taxonomy text comes for such a row. LOSAT refused such a run; it still refuses the first row of
  a query's report with `sacc`, `saccver` or `staxid`, where NCBI BLAST+ crashes or prints an uninitialised number.
