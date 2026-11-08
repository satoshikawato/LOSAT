- Python binding: `pylosat/_options.pyi` again matches the engine's table of options: BLASTP's `-task` accepts
  `blastp-short`, and the help texts of `-evalue`, `-threshold`, `-word_size`, `-window_size`, `-matrix` and `-gapopen` name
  its defaults (#937). The `Python` workflow now runs for every change under `LOSAT/src`, since the stubs come from the whole
  command line, not only `api/`.
