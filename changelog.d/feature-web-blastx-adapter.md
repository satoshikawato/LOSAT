- LOSAT Web ABI v2 runs BLASTX: `describe`, `validate`, `register` and `run` take `blastx` (outfmt 0, 6 and 7, the
  HSP records, a nucleotide query and protein subjects read by NCBI's reader). `validate` checks the options of each
  search's formats with the run's messages (`blastx::check_options`).
