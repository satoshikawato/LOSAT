- LOSAT Web: a query file whose reader fails on a later record (a line that NCBI's reader
  refuses, "there's a line that doesn't look like plausible data") now runs as on the command line,
  where it was refused before the run. The input card lists the records before that record and
  warns that a run that reaches it fails there with the command line's message, and that only the
  records searched before it keep their results (none in a short query); the run ends with that
  message, after the hits of the query batches before the error's. Under `-parse_deflines`, an
  error of an earlier query batch (two records
  with one ID) comes first, as on the command line. A query of comment lines only runs and ends with
  NCBI's `Empty CBlastQueryVector`. In a Combined query the warning names the run's line, which
  counts the lines of the inputs before it (those without records too), and the input's own line
  where the two differ; the inputs after it say that they are not searched. An input whose first
  line is not a defline cannot follow another input in a Combined query, as before: its card now
  says so, also when the reader stops at that line (a byte order mark), and the inputs after it no
  longer say that the search stops there. Such an input under `-parse_deflines` is no longer
  refused when its reader stops before the end of the index scan. ABI v2 is version 11 (*register*,
  *inputs* and *scan* give `reader_error`).
- Library: `read_input_records` keeps the query records read before a reader's error and the error
  (`QueryEnd`), and `run` and `run_batches` of them end with it in its batch, after the reports of
  the batches before; `read_inputs` is unchanged. Under `-parse_deflines` with `-subject`, a query
  ID repeated in a batch that NCBI fails (the reader's error, queries without letters, psiblast's
  batch of two queries) no longer takes the place of that batch's error, in BLASTX's `run` and in
  every program's `run_batches`. Every program's `run_batches` stops at the batch of a repeated ID
  in the command line's form, after the batches before it.
- Under `-parse_deflines`, a query ID repeated in a later query batch is now refused at that batch,
  after the reports of the batches before (in the first batch, after the prologs), in every program
  and output format, as the command line writes them when NCBI BLAST+ stops at that batch. An error
  that NCBI BLAST+ 2.17.0 gives in the search of an earlier batch now comes first, as NCBI's: for
  example psiblast with a query of only `X` before the repeated ID, or subjects too short for
  `-searchsp`, which LOSAT refused before with the repeated ID. Before, LOSAT wrote nothing before
  the refusal in BLASTP, TBLASTN, PSIBLAST and TBLASTX (TBLASTX's outfmt 0: its prolog only), and in
  BLASTN's and BLASTX's XML and JSON reports. The library's `run` and `run_batches` do the same.
- The batch of a repeated query ID now stops at one point in every program. Where the batch of
  NCBI BLAST+ 2.17.0 holds the same queries as LOSAT's (NCBI counts the repeated query into its
  batch by the first query's letters), an error of NCBI that does not depend on the repeated
  query's letters now comes instead of the refusal, as NCBI gives it: for example subjects too
  short for `-searchsp`, where `blastp` and `tblastn` refused. With `-export_search_strategy`,
  `blastn`, `blastx` and `tblastx` no longer write a search strategy for the refused batch; with
  one batch the file stays empty, as when NCBI stops at that batch.
- Under `-parse_deflines`, the query batches are now those of NCBI BLAST+ 2.17.0 when a query ID is
  repeated: like NCBI, LOSAT counts the repeated query into its batch by the first query's letters.
  A reader error after a short repeat now comes in the batch where NCBI meets it, where LOSAT
  refused the repeat before; a repeat that fills NCBI's batch is refused in that batch, where LOSAT
  gave the reader error of a later record. NCBI's errors of the repeat's batch now come instead of
  the refusal also where queries follow the repeat in the batch (subjects too short for
  `-searchsp`), and with a database whose sequence file cannot be read (`blastp`, `blastx` and
  `blastn` now give NCBI's "Database memory map file error", exit status 3, where they refused).
  In `blastx`, a query with the ID of an earlier query that `-query_loc` skips is now refused, where
  LOSAT searched it and NCBI stops with an error.
