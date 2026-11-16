- Engine: the refusals of BLASTN's query-anchored reports (outfmt 1 to 4) and its SAM writer
  (outfmt 17) move out of `post_process_hits_and_write` into their own module, and the subjects'
  names that three of the writers take are made in one place. No output changes.
