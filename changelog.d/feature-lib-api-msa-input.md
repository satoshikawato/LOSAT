- The library API and pylosat search psiblast's alignment query, `-in_msa`, with the command
  line's bytes: `Inputs::msa(path, &subject)` and `Inputs::msa_database(path, db)` in Rust,
  `Alignment.from_file(path)` in Python, against records or a database, with psiblast's rounds and
  PSSM files too. The one query is the alignment's master, query 0 of the results, in one batch of
  `run_batches`; the alignment is read as `-in_msa` reads it, with the same errors. `-msa_master_idx`
  and `-ignore_msa_master` are options of the argv (`msa_master_idx` and `ignore_msa_master` in
  Python, listed by `describe`); `-in_msa` is an input rather than an argv option, as `-in_pssm`
  (line LIB, item 1.6). The command line's search code does not change.
