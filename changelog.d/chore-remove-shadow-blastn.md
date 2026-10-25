- The SHADOW modes of blastn's engine switches `LOSAT_X_BNDBCUT`, `LOSAT_X_BNDBSTREAM`, `LOSAT_X_BNTBRANGE`,
  `LOSAT_X_CTXO1`, `LOSAT_X_QBATCHPAR`, `LOSAT_X_QGROUPPAR`, `LOSAT_X_QPREPPAR`, `LOSAT_X_QSTREAM`,
  `LOSAT_X_RPTSTREAM` and `LOSAT_X_QDUSTPAR` are removed (`LOSAT_X_<NAME>SHADOW` is ignored, and the `[X_SHADOW]`
  lines of these switches are no longer printed). The switches stay in the strict set with their `=0` paths; the
  `[X_STATS]` lines are unchanged. `LOSAT_X_QSTREAM` applies its rule (`off=megablast`, `batch`, `slots`) in every
  run. Outputs do not change. The checks that ran only under a shadow are unit tests of the live functions against the
  `=0` ones, or were dropped with the shadow.
- The SHADOW modes of blastn's lookup and seed switches `LOSAT_X_BNUNGAPCOUNT`, `LOSAT_X_DCCSR`, `LOSAT_X_DIAGLV`,
  `LOSAT_X_DISCLUTPAR`, `LOSAT_X_MBFILL`, `LOSAT_X_MBLUTKEEP`, `LOSAT_X_PAIRPARDISC`, `LOSAT_X_PAIRPARDISCARR` and
  `LOSAT_X_PAIRPARNA` are removed (`LOSAT_X_<NAME>SHADOW` is ignored, and the `[X_SHADOW]` and `[X_PAIRPAR*SHADOW]`
  lines of these switches are no longer printed). The switches stay in the strict set with their `=0` paths; the
  `[TIMING]` lines are unchanged. A one-thread search without a pool no longer runs the strips stages of
  `LOSAT_X_PAIRPARDISC`, `LOSAT_X_PAIRPARDISCARR` and `LOSAT_X_PAIRPARNA` on copies (only `LOSAT_X_BNDCONE` runs the
  dc-megablast stages there). Outputs do not change. The checks that ran only under a shadow are unit tests of the
  live functions against the `=0` ones, or were dropped with the shadow.
