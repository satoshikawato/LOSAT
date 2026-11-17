- Engine: the gapped half of a BLASTP subject's preliminary search (the order and the chaining of
  its ungapped hits, the gapped extensions with their exact retry, the purge of the HSPs with
  common endpoints and the subject's HSP lists with their E-values) moves out of `run_call`'s
  kernel into a module of its own. No output changes.
