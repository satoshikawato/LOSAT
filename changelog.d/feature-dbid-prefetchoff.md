- BLASTN, TBLASTN, BLASTX and TBLASTX `-db` with `taxdb.bti` on the search path and
  `PRE_FETCH_SEQS_LIMIT` turning NCBI's prefetch of a query batch off (0, or below the batch's
  number of alignments: one for each HSP of a gapped search, one for each subject of an ungapped
  one, as TBLASTX's) now write NCBI BLAST+ 2.17.0's messages about the taxonomy database after each
  query's own warnings, for the subjects its report shows, and none for a tabular report whose
  fields do not read the subjects (such as `6 qseqid evalue`; for BLASTN also `nident`, `mismatch`
  and `pident`; in TBLASTN, BLASTX and TBLASTX `sframe` and `frames` read them with a field of the
  alignment's coordinates), as NCBI's formatters fetch the subjects one query at a time. LOSAT wrote
  them before the warnings of the batch's queries, and also for such tabular reports, with exit
  status 0. Where every batch's prefetch is off (0 or a negative value) and the report reads no
  subject, the refusal of a database whose taxids NCBI's taxonomy cache would mix up no longer
  applies, as NCBI looks no taxid up.
