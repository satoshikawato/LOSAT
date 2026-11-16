- Engine: the drivers of a BLASTP call's preliminary stage (the choice between threads and one
  thread, and the RPS-BLAST, PSSM, query block, database and subject set branches) move out of
  `run_call` into a module of their own, which reads one set of the call's values. No output
  changes.
