- LOSAT Web: a query file whose reader fails on a later record (a line that NCBI's reader
  refuses, "there's a line that doesn't look like plausible data") now runs as on the command line,
  where it was refused before the run. The input card lists the records before that record and
  warns with the command line's message, and the run searches the query batches before it and ends
  with that message. Under `-parse_deflines`, an error of an earlier query batch (two records with
  one ID) comes first, as on the command line. A query of comment lines only runs and ends with
  NCBI's `Empty CBlastQueryVector`. In the library, `read_input_records` keeps the query records
  read before a reader's error and the error (`QueryEnd`), and `run` and `run_batches` of them end
  with it in its batch, after the reports of the batches before; `read_inputs` is unchanged.
