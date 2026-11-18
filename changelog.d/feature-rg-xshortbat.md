- When a `tblastn` search of several query batches (`BATCH_SIZE`, or queries whose letters fill a batch) stops with an
  error in a later batch, LOSAT now writes the results of the batches before it first, as NCBI BLAST+ 2.17.0 prints
  each batch's results before it searches the next; before, those results were not written. This covers each
  `tblastn` error at a step where NCBI BLAST+ crashes, does not finish or stops: an alignment on a subject shorter
  than 12 bases, two `-window_size 0` alignments on which NCBI BLAST+ never finishes, the overflow banks of the
  compressed lookup table, a query chunk that NCBI splits again or that holds no query (also with a negative
  `CHUNK_SIZE`), and with `-in_pssm` a window past the PSSM, a gapped start past the subject's translation, a traceback
  that starts past the end of its frame or a PSSM query that `-seg` masks completely (these write the prologs
  first). Searches without such an error are unchanged.
- `tblastn` stops at such a query chunk also when a seqid list leaves the database without a sequence to search, as
  NCBI BLAST+ 2.17.0 does; before, it wrote that batch's empty results. With `MIN_SEQ_LENGTH 1` or a database
  without sequences, the error of the database's first read and the query chunk stops now come in NCBI BLAST+'s
  order.
