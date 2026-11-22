- When a `tblastn` search of several query batches (`BATCH_SIZE`, or queries whose letters fill a batch) stops with an
  error in a later batch, LOSAT now writes the results of the batches before it first, as NCBI BLAST+ 2.17.0 prints
  each batch's results before it searches the next; before, those results were not written. This covers each
  `tblastn` error at a step where NCBI BLAST+ crashes, does not finish or stops: an alignment on a subject shorter
  than 12 bases, two `-window_size 0` alignments on which NCBI BLAST+ never finishes, the overflow banks of the
  compressed lookup table, a query chunk that NCBI splits again or that holds no query (also with a negative
  `CHUNK_SIZE`), and with `-in_pssm` a window past the PSSM, a gapped start past the subject's translation, a traceback
  that starts past the end of its frame or a PSSM query that `-seg` masks completely (these write the prologs
  first). Searches without such an error are unchanged.
- `tblastn` stops at such a query chunk also when the database is left without a subject (a seqid list that keeps
  none, an alias over volumes without sequences), as NCBI BLAST+ 2.17.0 does; before, it wrote that batch's empty
  results. With a database whose shortest sequence is under 3 bases (`MIN_SEQ_LENGTH` 1 or 2), one without sequences
  or one whose sequence file is missing, the database's error and the query chunk stops now come in NCBI BLAST+'s
  order.
- A negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` whose pieces leave gaps between them, which NCBI
  BLAST+ searches with the hits there lost and LOSAT refuses, now stops `tblastn` where the search of that batch's
  first piece with a valid query would read the subjects, after the results of the query batches before and after
  NCBI BLAST+'s own errors of that piece; with no subject to search nothing is lost, and the batch is reported. Before,
  it stopped before any output.
- A `tblastn` batch that LOSAT refuses for two queries with one ID (`-parse_deflines`) and that NCBI BLAST+ splits
  into pieces gives NCBI BLAST+'s error of the batch's setup (subjects too short for `-searchsp`, for example) only
  where NCBI BLAST+'s pieces reach that setup, and the refusal where they stop first.
