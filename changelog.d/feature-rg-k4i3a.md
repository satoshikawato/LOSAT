- The library API's `describe`, `pylosat.describe` and LOSAT Web's option list read every part of an
  option from NCBI's argument tables, which the command line's reader reads: the help of an option
  (`Parameter::help`, the pylosat docstrings of the `Search` keywords, the "?" help of LOSAT Web's search
  form) is NCBI BLAST+ 2.17.0's text of the option, the one that `-help` prints, in LOSAT's words
  (`Expectation value (E) threshold for saving hits. Default = 10`), where it was LOSAT's own text, and
  the options come in the order of `-help`. The options listed, their defaults and their choices are
  unchanged. psiblast's `-out_pssm` and `-out_ascii_pssm` are `ValueKind::Path` (pylosat: `path`), as
  NCBI's output files, where they were `Text`.
