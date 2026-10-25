- The SHADOW modes of blastn's engine switches `LOSAT_X_BNDBCUT`, `LOSAT_X_BNDBSTREAM`, `LOSAT_X_BNTBRANGE`,
  `LOSAT_X_CTXO1`, `LOSAT_X_QBATCHPAR`, `LOSAT_X_QGROUPPAR`, `LOSAT_X_QPREPPAR`, `LOSAT_X_QSTREAM`,
  `LOSAT_X_RPTSTREAM` and `LOSAT_X_QDUSTPAR` are removed (`LOSAT_X_<NAME>SHADOW` is ignored, and the `[X_SHADOW]`
  lines of these switches are no longer printed). The switches stay in the strict set with their `=0` paths; the
  `[X_STATS]` lines are unchanged. `LOSAT_X_QSTREAM` applies its rule (`off=megablast`, `batch`, `slots`) in every
  run. Outputs do not change. The checks that ran only under a shadow are unit tests of the live functions against the
  `=0` ones, or were dropped with the shadow.
