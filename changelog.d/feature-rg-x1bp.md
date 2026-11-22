- BLASTN: an alignment that NCBI BLAST+ keeps outside its query context (from an ungapped
  extension past an end of the context, such as with `-xdrop_ungap 1e10`) whose query coordinates
  pass the query's end is no longer refused: every report shows it as NCBI BLAST+ 2.17.0 does,
  with `-subject` and from a database, on both strands. Its query letters are those of the query's
  record, also outside a `-query_loc` interval, and past the record's end the gap character `-`
  (outfmt 0 to 5, 13 to 16, the tabular `qseq` and `btop`), which the identities that the reports
  count from the letters count as mismatches; SAM's `SEQ` leaves out the letters past the end. With
  `-query_loc`, an ungapped alignment that leaves the interval on the left but not the query's
  record is shown too, with the record's letters before the interval. Where its query coordinates
  leave the query's record on the left, NCBI BLAST+ stops with an exception with `-subject`; LOSAT
  now refuses it there after what NCBI BLAST+ writes before: the start of the XML2 and JSON
  documents, outfmt 20's line of field names, and the reports of the query batches before (manual,
  What LOSAT refuses).
- BLASTN: with `-ungapped` and a `-query_loc` range whose end passes the query's end, `qcovhsp` of
  a subject's only alignment counts its letters against the range as given, as NCBI BLAST+ 2.17.0
  does (LOSAT counted them against the searched letters); a subject's several alignments and
  gapped ones keep the searched letters.
