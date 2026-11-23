- Each option's help in the library API's `describe` (`Parameter::help`), in `pylosat` (the docstrings
  of the `Search` keywords) and in LOSAT Web's search form (the "?" help) is now the text that
  `losat <program> -help` prints for it (NCBI BLAST+ 2.17.0's text, in LOSAT's words: `Expectation
  value (E) threshold for saving hits. Default = 10`), so the command line, the library and the Web
  describe an option in the same words, and the options that had no help (`-lcase_masking`,
  `-strand`, `-ungapped`, `-use_sw_tback`, ...) have one. The options come in the order of the
  usage line of `-help` (blastn: `-task`, `-db`, `-dbsize`, ...). The options listed, their
  defaults and their choices are unchanged. psiblast's `-out_pssm` and `-out_ascii_pssm` are
  `ValueKind::Path` (pylosat: `path`), as for NCBI's other output files, where they were `Text`.
