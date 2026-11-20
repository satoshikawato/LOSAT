- `-import_search_strategy` in `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`: runs the text strategy
  files that NCBI BLAST+ and LOSAT write (`-export_search_strategy`) with NCBI's output, including the
  strategy's queries, subjects or database, its options and the command line's overrides; `blastn` and
  `blastx` run a `tblastx` strategy and `tblastx` a `blastn` or `blastx` one as NCBI does. A strategy's
  subjects are searched once whatever `-num_threads` says (NCBI prints every alignment once per thread, or
  stops). The library API takes a strategy as bytes (`Search::read_strategy`, `Inputs::strategy`), and
  LOSAT Web's ABI takes it as file kind 3. Binary and XML strategy files, and those of `psiblast`,
  `rpsblast` and `rpstblastn`, are refused for now.
