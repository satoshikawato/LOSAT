- LOSAT Web: a query file whose reader fails on a later record (a line that NCBI's reader
  refuses, "there's a line that doesn't look like plausible data") now runs as on the command line,
  where it was refused before the run. The input card lists the records before that record and
  warns that a run of the input ends with the command line's message and gives no results; the run
  ends with that message. Under `-parse_deflines`, an error of an earlier query batch (two records
  with one ID) comes first, as on the command line. A query of comment lines only runs and ends with
  NCBI's `Empty CBlastQueryVector`. In a Combined query the warning names the run's line, which
  counts the lines of the inputs before it (those without records too), and the input's own line
  where the two differ; the inputs after it say that they are not searched. An input whose first
  line is not a defline cannot follow another input in a Combined query, as before: its card now
  says so, also when the reader stops at that line (a byte order mark), and the inputs after it no
  longer say that the search stops there. Such an input under `-parse_deflines` is no longer
  refused when its reader stops before the end of the index scan.
- Library: `read_input_records` keeps the query records read before a reader's error and the error
  (`QueryEnd`), and `run` and `run_batches` of them end with it in its batch, after the reports of
  the batches before; `read_inputs` is unchanged. Under `-parse_deflines` with `-subject`, a query
  ID repeated in a batch that NCBI fails (the reader's error, queries without letters, psiblast's
  batch of two queries) no longer takes the place of that batch's error, in BLASTX's `run` and in
  every program's `run_batches`. BLASTX's `run` and BLASTN's and BLASTX's `run_batches` refuse a
  repeated ID after the batches before it in outfmt 6, and TBLASTX's `run_batches` once it has
  searched the batch of the second ID, an error of a batch up to it coming first, as each command
  line does.
