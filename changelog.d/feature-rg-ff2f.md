- All programs (command line, library API and LOSAT Web): a `-logfile` naming a file that the strategy
  of `-import_search_strategy` makes the run read (its database, or the command line's `-db`, opened as
  the type of the strategy program's subjects; its repeat filter database; the counts of its window
  masker taxid; on the command line and in the library API also its forced MegaBLAST index) is now
  refused as the run's other inputs are, with exit status 1 before the file is opened or created;
  NCBI BLAST+ 2.17.0 opens the log before it reads the strategy and appends its messages to that file.
  A strategy read from the standard input or a pipe is not checked: it is read once, by the import. In
  the library API the database and the PSSM or alignment file given as `Inputs` count too.
- A run that imports a strategy no longer refuses a `-logfile` naming a file of an option that the
  import ignores: the list files (`-seqidlist`, `-gilist`, `-taxidlist`, ...) without `-db` or
  `-subject`, BLASTN's `-filtering_db`, `-window_masker_db` and `-window_masker_taxid`, `-index_name`
  without `-use_index true`, and `-in_pssm`. NCBI BLAST+ 2.17.0 warns that they cannot be overridden
  and does not read them. LOSAT Web still refuses a log over a window masker counts file that the
  page holds, as it refuses any file of the run's directory.
- blastn and blastx running a TBLASTX strategy, and tblastx running a BLASTN or BLASTX one, name the
  running program in the refusal of a subject field whose lookup fails, in the taxonomy database
  warning written before that refusal (`Warning: [blastn] ...`, as NCBI BLAST+ 2.17.0 writes it), and
  in LOSAT's other refusals of the run; they named the strategy's program.
