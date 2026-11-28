- `-db` searches of `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` with `-outfmt` 0 to 7, 10, 13 to 16, 18 or
  20 (`-outfmt` 0 to 4 with and without `-html`), `blastn` with `-outfmt 17`, and `-outfmt` 8, 9 and 12 with
  `PRE_FETCH_SEQS_LIMIT=2147483647` now run on a database whose stored titles hold bytes outside the printable ASCII
  range: the output and NCBI BLAST+'s error messages for each such byte are NCBI's, where LOSAT refused them before.
  NCBI writes them where each query batch's prefetch loads the sequence (in a batch of more than 200 sequences, in
  NCBI's pieces of 100), again at each row whose fields read the title, ids or taxonomy of the definition line, at
  each description line (twice) and alignment heading (once) of the pairwise reports, at each hit of the XML and
  JSON reports and at each alignment (twice) of the organism report. These runs are still refused where a query
  batch asks for a sequence that NCBI may still keep from an earlier batch or not, depending on the order of objects
  in its memory, or one whose header shares NCBI's header cache slot with such a sequence; and, until LOSAT
  reproduces them, where the prefetch can be off (`PRE_FETCH_SEQS_LIMIT` set to another value), on a database with
  parsed ids, where NCBI's object manager settings (`OBJMGR_*`) are set, and in a query batch whose prefetch makes a
  damaged taxonomy database write a message (its first lookup, or one that NCBI's cache does not hold or may not hold).
