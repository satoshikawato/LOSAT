- The library API and pylosat search a PSSM query, psiblast's and tblastn's `-in_pssm`, with the
  command line's bytes: `Inputs::pssm(path, &subject)` and `Inputs::pssm_database(path, db)` in
  Rust, `Pssm.from_file(path)` in Python, against records or a database, with psiblast's rounds and
  PSSM files too. The one query is the PSSM's sequence, query 0 of the results, in one batch of
  `run_batches` (line RG, Q6b). `-in_pssm` was refused, and is now an input rather than an argv
  option. LOSAT Web keeps refusing tblastn's `-in_pssm` (line PS, PSW1).
