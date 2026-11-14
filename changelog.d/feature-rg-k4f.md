- BLASTP, BLASTX, TBLASTN and TBLASTX with a `-matrix` name that NCBI BLAST+ 2.17.0's BLAST XML
  report does not spell (such as `-matrix blosum62`) and outfmt 13 to 16 against `-subject` (and
  TBLASTN `-in_pssm`) stop where NCBI stops: at the batch of the second query, or at the end of a
  run of one query, instead of the first batch, with the same search strategy requests
  (`-export_search_strategy`), reader messages and error. BLASTX no longer refuses a query title
  of such a run that its XML could not write: NCBI writes none. TBLASTX reports the errors of the
  batch's search and prefetch (`The average subject length is too short`, `PRE_FETCH_SEQS_LIMIT`)
  before that error, as NCBI does.
