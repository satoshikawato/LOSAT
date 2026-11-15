- The library API and pylosat run psiblast's rounds (`-num_iterations` other than 1): the reports
  are the command line's bytes, each HSP record has its round (`HspRecord::round`, pylosat's `round`
  column), and each query its rounds' statistics and whether they converged (`Outcome::rounds`,
  `Batch::rounds`, `QueryRounds`, `QueryHits::rounds`) (line RG, Q5). They were refused.
