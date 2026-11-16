- The library API and pylosat write psiblast's PSSM files (`-out_pssm`, `-out_ascii_pssm`, with
  `-save_pssm_after_last_round` and `-save_each_pssm`) where the command line writes them, byte for
  byte; `-save_each_pssm` numbers them over the whole search, also over `run_batches`' batches of one
  query each, as NCBI numbers them over its run (line RG, Q6a). They were refused. A call of `run` or
  `run_collect` whose formats take two searches (`-num_descriptions` or `-num_alignments` with a
  report and a table) and that writes PSSM files is refused: one call per format writes the command
  line's files.
