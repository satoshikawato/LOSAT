- `-db` lists of version 5 databases made with and without `makeblastdb -parse_seqids` (`-db "a b"`, or an alias file
  that lists them): where a query batch reports a parsed Seq-id (an accession, a name, a giim) that a database of the
  list cannot look up, every search program now stops as NCBI BLAST+ does, with `LOSAT database error: DB contains no
  accession info.` (exit status 2), after the header and the reports of the batches before it, in every format that
  reads the reported sequences (the pairwise report, the tabular formats, XML, XML2 and JSON documents and side files,
  SAM, and the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`); LOSAT refused such a run before. With
  `PRE_FETCH_SEQS_LIMIT` turning the prefetch off, psiblast's PSSM drops such a sequence's letters with `Failed to
  retrieve sequence` and goes on, as NCBI's does, where the reports do not read the sequence. Still refused: such a
  Seq-id that the reports read with the prefetch off, and a stopping batch whose other sequences make the taxonomy
  database write a message or fail on their ids.
- blastp and psiblast: where a later query batch stops the run (also at a header file that NCBI cannot read), the
  taxonomy database's messages of the batches before it are written again before the error, as NCBI writes them.
