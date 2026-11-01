- The variables `LOSAT_X_BXDEFER` and `LOSAT_X_BXDEFERMIN` are gone: no default blastx path reached the rule they
  chose (the masked AVX2 batch of the fused seed loop, which still runs with a value of `LOSAT_X_BXDIAG` where
  `LOSAT_X_BXCROSS`'s cell test does not), and `LOSAT_X_BXDIAGMIN` has one built-in value, 670, for every cell test
  (the `[TIMING] blastx_diag` line lost `defer_min_hits_per_cell=` and `tables_avx2=`). Default results and paths do not change.
