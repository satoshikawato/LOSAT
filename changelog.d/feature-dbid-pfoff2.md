- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: the tabular formats (`-outfmt 6`, `7`, `10` and `20`) and the organism report (`-outfmt 18`) now give the
  results of NCBI BLAST+ instead of LOSAT's refusal, in BLASTP, psiblast (tabular formats), BLASTN, BLASTX, TBLASTN
  and TBLASTX. A field list with an id, title, strand or taxonomy field (the default fields too) writes that row with
  the scores, zero counts and coordinates and the ids, titles and taxonomy of the query's last row that read its
  sequence (`Unknown` and `N/A` before one), adds the Seq-align id to `sallseqid`, `sallgi` and `sallacc`, and goes
  on; a field list that reads the sequence otherwise (`slen`, the aligned letters, a count computed from them, or in
  BLASTX, TBLASTN and TBLASTX the subject's frame) writes the rows before it and stops with `LOSAT database error: DB
  contains no accession info.` (exit status 2), and `nident` alone stops only at a row with identities. The organism
  report leaves such a sequence out and goes on. No taxonomy text comes for such a row, and a damaged or empty
  taxonomy database is read only for the sequences that were read. LOSAT still refuses the first row of a query's
  report with `sacc`, `saccver` or `staxid`, where NCBI BLAST+ crashes or prints an uninitialised number.
