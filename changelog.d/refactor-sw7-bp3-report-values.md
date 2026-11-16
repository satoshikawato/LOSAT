- Engine: the steps of a BLASTP call between its final hit lists and the rows of its reports
  (the residues that the reports read of a database's subjects, the check of the prefetch
  environment, the subjects that the reports show, psiblast's cut and PSSMs with the batches'
  prefetch, and the labels of the reports) move out of `run_call` into the functions of its
  report module. No output changes.
