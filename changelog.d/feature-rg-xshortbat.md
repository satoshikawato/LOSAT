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
- With `-searchsp` or `OLD_FSC`, a `tblastn` search of a database of volumes without sequences (an alias over them)
  now fails with NCBI BLAST+ 2.17.0's `Database memory map file error`; before, LOSAT gave `The average subject length
  is too short`, an error that NCBI BLAST+ checks only after that read.
- `blastp` and `psiblast` now take the same steps as `tblastn` with split queries: a negative `CHUNK_SIZE` above a
  negative `OVERLAP_CHUNK_SIZE` whose pieces leave gaps stops the search where the batch's first piece with searchable
  query letters would read the subjects, after the results of the query batches before; with no subject to search (a
  seqid list that keeps none) or no such piece, the batch is reported as NCBI BLAST+ 2.17.0 reports it. Before, these
  searches stopped before the batch was searched.
- With a split query, `blastp` and `psiblast` now stop at a piece that NCBI BLAST+ would split again, or one left
  without a query, before the database's errors (no sequences, a missing sequence file) and `-searchsp`'s average
  subject length error, as NCBI BLAST+ 2.17.0 does; a batch whose queries are all unsearchable is no longer refused
  for its pieces (NCBI BLAST+ does not search it).
- A split `blastp` or `tblastn` batch none of whose pieces has searchable query letters now fails with NCBI BLAST+'s
  `Traceback failed` where `-searchsp` leaves the subjects too short, and is reported where the database's sequence
  file is missing (NCBI BLAST+ reads nothing for it); a later batch that reads the database then fails after its
  report. Before, LOSAT gave the average subject length error or the database's error at that batch.
- A `blastp` batch that LOSAT refuses for two queries with one ID and that NCBI BLAST+ splits now gives NCBI BLAST+'s
  error of the batch's setup where its pieces reach it, as `tblastn` does.
