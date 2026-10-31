- The variable `LOSAT_X_SUMPMEMO` is gone: the per-thread reuse of the sum P-value of a set of linked HSPs
  (BLASTX, TBLASTN and TBLASTX: the last two calls of `blast_sum_p_calc` with the same `r` and bits of `s`)
  was already the default and now runs unconditionally, so `=0` no longer integrates every value again.
  Outputs do not change.
