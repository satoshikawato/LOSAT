- BLASTN, TBLASTN, BLASTX, TBLASTX and BLASTP `-db` with `PRE_FETCH_SEQS_LIMIT=2147483647` and
  outfmt 8, 9 or 12 now prefetch the reported subjects as NCBI BLAST+ 2.17.0 does at that value
  (its prefetch returns before it tests the format; the archive, outfmt 11, is written without
  it): with `taxdb.bti` on the search path they write NCBI's messages about the taxonomy database
  after the reader's messages of each batch, as for outfmt 0, 5, 6 and 7, and refuse the taxonomy
  databases and tax ids those formats refuse. LOSAT wrote no message there, with exit status 0.
  A subject with parsed Seq-ids is read as NCBI's prefetch reads it by its Seq-align id: a lookup
  that fails (NCBI writes `Error pre-fetching sequence data`), another sequence whose definition
  lines do not hold the Seq-align id (NCBI stops, `oid headers do not contain target gi/seq_id`),
  and in BLASTN and TBLASTX a sequence that NCBI gets no Bioseq of with an alignment on the minus
  strand of the query and the plus strand of the subject (NCBI's prefetch fails to read its
  length) are refused, with or without a taxonomy database, and so are copies and other sequences
  whose taxonomy database texts can differ from what NCBI reads, as in outfmt 6. The same check of
  another sequence without the Seq-align id now also applies to psiblast's list beyond
  `-max_target_seqs`.
