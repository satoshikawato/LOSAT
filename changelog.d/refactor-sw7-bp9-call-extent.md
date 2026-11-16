- Engine: the setup of a BLASTP call's statistics and scan (the database's length and number
  of subjects, the search spaces of the setup and of the preliminary stage, the preliminary
  cutoffs, the statistics of each subject, the hit lists that the search fills and the scan's
  X-drops and sizes) moves out of `run_call` into the call's setup module. No output changes.
