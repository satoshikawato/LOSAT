- The library API and pylosat run psiblast's rounds (`-num_iterations` other than 1): the reports
  are the command line's bytes, each HSP record has its round (`HspRecord::round`, pylosat's `round`
  column), and each query its rounds' statistics and whether they converged (`Outcome::rounds`,
  `Batch::rounds`, `QueryRounds`); `QueryHits::rounds` gives each round's HSPs and statistics, and
  `QueryHits::round` and `SubjectHits::round` tell a view's round (line RG, Q5). They were refused.
- The library's statistics of a PSSM's search (psiblast's later rounds) are the PSSM's gapped
  Karlin-Altschul values, as NCBI's (`kbp_gap` is `kbp_gap_psi` there).
