- LOSAT Web: a query file whose reader fails on a later record (a line that NCBI's reader
  refuses, "there's a line that doesn't look like plausible data") now runs as on the command line,
  where it was refused before the run. The input card lists the records before that record and
  warns with the command line's message, and the run searches the query batches before it and ends
  with that message. Under `-parse_deflines`, an error of an earlier query batch (two records with
  one ID) comes first, as on the command line. A query of comment lines only runs and ends with
  NCBI's `Empty CBlastQueryVector`. In the library, `read_input_records` keeps the query records
  read before a reader's error and the error (`QueryEnd`), and `run` and `run_batches` of them end
  with it in its batch, after the reports of the batches before; `read_inputs` is unchanged.
  Under `-parse_deflines` a query ID repeated in the batch that the reader's error ends no longer
  takes the place of that error (BLASTX's `run` and every program's `run_batches`), and BLASTX's
  `run` and BLASTN's and BLASTX's `run_batches` now refuse a repeated ID after the batches before
  it in outfmt 6, as the command line does. LOSAT Web: such an input under `-parse_deflines` is no
  longer refused when its reader stops before the end of the index scan; in a Combined query an
  input whose first line is not a defline and that the reader stops at (a byte order mark) is
  refused after another input, as before; the warning says what to correct and names the run's
  line, and the inputs after it say that they are not searched.
