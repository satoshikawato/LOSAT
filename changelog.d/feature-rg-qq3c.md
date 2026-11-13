- BLASTX's ungapped search (`-ungapped` with `-max_intron_length` at most 4) links its HSPs with
  the same even-gap linking as PSIBLAST's `-sum_stats true` (`core::link_hsps`), and BLASTX and
  TBLASTX compute the linking cutoffs (`CalculateLinkHSPCutoffs`) with one function. Outputs do not
  change; new frozen NCBI rows `bx.q3c.*` cover linked sets across frames, small search spaces,
  same-range HSPs in two frames and `-max_intron_length` with and without `-ungapped`.
