- pylosat: a query file whose reader fails on a later record now searches as the command line
  does. The records before it are searched one query batch at a time, and the search ends with the
  reader's message in the batch of that record. `batches` gives the batches before it, and
  `run_to_file` writes their reports, as the command line's `-out` does; then `SearchError` is
  raised. Before, the binding read every query first and raised `InputError` before any search. A
  query record that the reader refuses, the first one too, now raises `SearchError`; a file that
  does not open, or a subject record that the reader refuses, still raises `InputError`. A query
  of comment lines only ends with NCBI's `Empty CBlastQueryVector`, as on the command line.
