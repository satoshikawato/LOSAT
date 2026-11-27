- All programs (command line, library API and LOSAT Web): a `-logfile` naming a file that the strategy
  of `-import_search_strategy` makes the run read (its database, its repeat filter database, the counts
  of its window masker taxid or its forced MegaBLAST index) is now refused as the run's other inputs are,
  with exit status 1 before the file is opened or created; NCBI BLAST+ 2.17.0 opens the log before it
  reads the strategy and appends its messages to that file. In the library API the database and the
  PSSM or alignment file given as `Inputs` count too.
- blastn and blastx running a TBLASTX strategy, and tblastx running a BLASTN or BLASTX one, name the
  running program in the taxonomy database warning (`Warning: [blastn] ...`) and in the refusal of a
  subject field whose lookup fails, as NCBI BLAST+ 2.17.0 does; they named the strategy's program.
