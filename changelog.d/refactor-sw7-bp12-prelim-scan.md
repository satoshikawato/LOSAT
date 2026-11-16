- Engine: the scan half of a BLASTP subject's preliminary search (the scan of its words with the
  two-hit or one-hit test and the ungapped extensions, for every matrix, lookup table and
  `-lcase_masking`) moves out of `run_call` into the module of its gapped half, so the
  per-subject kernel is in one place. No output changes.
