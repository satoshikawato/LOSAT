- BLASTX's ungapped search (`-ungapped` with `-max_intron_length` at most 4) links its HSPs with
  the same even-gap linking as PSIBLAST's `-sum_stats true` (`core::link_hsps`), and BLASTX and
  TBLASTX compute the linking cutoffs (`CalculateLinkHSPCutoffs`) with one function. Outputs do not
  change, except for a linked set of more than 32,767 HSPs, whose size now wraps as NCBI's `Int2`
  count does (E-value 0 for a count below 1) in PSIBLAST too. New frozen NCBI rows `bx.q3c.*` cover
  linked sets across frames, small search spaces, same-range HSPs in two frames and
  `-max_intron_length` with and without `-ungapped`.
