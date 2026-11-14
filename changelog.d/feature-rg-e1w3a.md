- LOSAT Web's engine runs the NCBI C++ Toolkit options `-logfile`, `-conffile` and `-dryrun` in blastn, blastp, blastx,
  tblastn and tblastx as the command line does. They were refused. The file of `-conffile` is a file that the page gives
  the engine, the `-logfile` file comes back with the run's outputs, and a dry run checks its arguments and writes no
  output. The page's controls for them come in a later change.
