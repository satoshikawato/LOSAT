- BLASTN: an alignment that NCBI BLAST+ keeps outside its query context (from an ungapped
  extension past an end of the context, such as with `-xdrop_ungap 1e10`) whose query coordinates
  pass the query's end is no longer refused: every report shows it as NCBI BLAST+ 2.17.0 does,
  with `-subject` and from a database, on both strands. Its query letters past the query's end show
  the gap character `-` (outfmt 0 to 5, 15 and 16, the tabular `qseq` and `btop`), the identities
  that the reports count from the letters count them as mismatches, and SAM's `SEQ` holds the
  letters before the end (none on the reverse strand). Where its query coordinates leave the query
  on the left, NCBI BLAST+ stops with an exception with `-subject`; LOSAT now refuses it there
  after what NCBI BLAST+ writes before: the start of the XML2 and JSON documents, outfmt 20's line
  of field names, and the reports of the query batches before (README, What LOSAT refuses).
