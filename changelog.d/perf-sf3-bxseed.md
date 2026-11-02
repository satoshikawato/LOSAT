- The variables `LOSAT_X_BXDEFER` and `LOSAT_X_BXDEFERMIN` are gone: no default blastx path reached the rule they
  chose, and `LOSAT_X_BXDIAGMIN` has one built-in value, 670, for every cell test (the `[TIMING] blastx_diag` line
  lost `defer_min_hits_per_cell=` and `tables_avx2=`). Results do not change anywhere, and paths do not change on
  x86-64 with AVX2 and POPCNT. What changes: `LOSAT_X_BXCROSS=0`, or an AVX2 CPU without POPCNT, now takes the
  direct loop on every table when `LOSAT_X_BXDIAG` is not set, and `LOSAT_X_BXDIAG=<value>` or
  `LOSAT_X_BXDIAGSHADOW` with `LOSAT_X_BXCROSS=0` takes the masked AVX2 batch (`mode=avx2-defer-masked`) with the
  cutoff 670 instead of 900.
