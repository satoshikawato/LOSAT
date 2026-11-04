- LOSAT Web: records left out and a subrange (`-query_loc` / `-subject_loc`) are no longer lost when
  you pass through a mode that reads that input with the other FASTA reader kind (for example
  BLASTN, BLASTP, BLASTN); each source keeps its record table of both kinds, so the way back reads
  nothing again. A run with an empty query completes as the command line does, and its results say
  "The query has no sequences." with the engine's warning instead of "No queries match the query
  filters." After the browser removes the tab's storage (site data cleared), the next run makes it
  again, and a run that cannot says the storage was removed and to reload, not the browser's own
  words. A failed run's message ends with one full stop.
