- When a `tblastn` search of several query batches (`BATCH_SIZE`, or a first query long enough to fill a batch) stops
  with an error in a later batch, LOSAT now writes the results of the batches before it first, as NCBI BLAST+
  2.17.0 prints each batch's results before it searches the next; before, those results were
  not written. This covers each `tblastn` error at a step where NCBI BLAST+ crashes, does not finish or stops: an
  alignment on a subject shorter than 12 bases, two `-window_size 0` alignments on which NCBI BLAST+ never finishes,
  the overflow banks of the compressed lookup table, a query split that NCBI splits again, and with `-in_pssm` a
  window past the PSSM or a gapped start past the subject's translation (these write the report's header first).
  Searches without such an error are unchanged.
