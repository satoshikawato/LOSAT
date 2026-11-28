- LOSAT Web: with a database folder picked, `-logfile taxdb.btd` and `-logfile taxonomy4blast.sqlite3`
  failed in the page where the shown command runs on the command line. The folder's taxonomy files
  were also placed at the top of the run's directory, so the run saw files that the user's directory
  does not have. They now stay where the user has them, and the engine finds them through the shown
  `BLASTDB=<folder>` as the command line does. Taxonomy lookups (`-taxids`, `sscinames`) give the
  same results.
- `-logfile` named like a taxonomy file that the run looks for is now refused on the command line as
  in LOSAT Web, before the file is made: `taxdb.bti` (and the `taxdb.btd` beside the one found) for
  every search of a database, for a `-subject` search with `sscinames`-like fields or `-outfmt 18`,
  and for `dustmasker`, `segmasker`, `windowmasker`, `makembindex`, `blastdbcmd` and `blastdbcheck`
  with a database; `taxonomy4blast.sqlite3` with `-taxids` or `-taxidlist`. The run looks for these names in its
  working directory first, so it used to find its own empty log there instead of your taxonomy
  files and lost the names. LOSAT Web now refuses the same names, also where it let them run before
  (`-subject -outfmt 18`, `-db` with `-outfmt 5`, `8`, `11` or `"6 qseqid pident"`); `-subject` with
  other formats still runs with them, as on the command line.
