- LOSAT Web: with a database folder picked, `-logfile taxdb.bti` (and `taxdb.btd`,
  `taxonomy4blast.sqlite3`) failed in the page where the shown command runs on the command line. The
  folder's taxonomy files were also placed at the top of the run's directory, so the run saw files
  that the user's directory does not have. They now stay where the user has them, and the engine finds
  them through the shown `BLASTDB=<folder>` as the command line does, so every `-logfile` and `-out`
  name behaves as on the command line. Taxonomy lookups (`-taxids`, `sscinames`) give the same results.
