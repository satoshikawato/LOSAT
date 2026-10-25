- The SHADOW modes of `LOSAT_X_BPCLUT`, `LOSAT_X_BPLKPFLAT`, `LOSAT_X_BPLUT`, `LOSAT_X_BPLKPAHEAD` and
  `LOSAT_X_TBXLUT` are removed (`LOSAT_X_BPCLUTSHADOW`, `LOSAT_X_BPLKPFLATSHADOW`, `LOSAT_X_BPLUTSHADOW`,
  `LOSAT_X_BPLKPAHEADSHADOW` and `LOSAT_X_TBXLUTSHADOW` are ignored, and the `[X_SHADOW]` lines of these switches
  are no longer printed). The switches stay in the strict set with their `=0` paths. With `LOSAT_X_STATS`, the
  live counters that the shadow lines showed are printed as `[X_STATS] BPLKPFLAT tables= recycled=
  recycled_hits=` and `[X_STATS] BPLKPAHEAD blocks= prebuilt= dropped=` (only when a counter is not zero).
  Outputs do not change. The check that a block's values built ahead (query frames, table, contexts) are the
  values of its own build is now a unit test in `x_lookup_ahead.rs`; the other checks that ran only under a
  shadow are the unit tests that compare the live builders with the reference builders.
- The SHADOW modes of `LOSAT_X_BPDBCUT`, `LOSAT_X_BPLISTEDLEAN`, `LOSAT_X_BPQBLOCKFMT`, `LOSAT_X_BPREDODBSUBJ`,
  `LOSAT_X_BPSUBJMETA` and `LOSAT_X_BPSUBJRUN` are removed (`LOSAT_X_BPDBCUTSHADOW`, `LOSAT_X_BPLISTEDLEANSHADOW`,
  `LOSAT_X_BPQBLOCKFMTSHADOW`, `LOSAT_X_BPREDODBSUBJSHADOW`, `LOSAT_X_BPSUBJMETASHADOW` and
  `LOSAT_X_BPSUBJRUNSHADOW` are ignored, and the `[X_SHADOW]` lines of these switches are no longer printed). The
  switches stay in the strict set with their `=0` paths; `LOSAT_X_BPQBLOCKFMTSHADOW` no longer turns the other
  block-search switches on. Outputs do not change. The checks that ran only under a shadow are now unit tests
  (the lean titles and letters against the records, the diagonal offset and table state of a run's later subjects)
  or are dropped where the reports are compared with the `=0` path and the reach rows show that the path ran.
- The SHADOW modes of `LOSAT_X_BPGAP`, `LOSAT_X_BPREDOPAR`, `LOSAT_X_BPREDOWINEXPAND`, `LOSAT_X_BPSEED`,
  `LOSAT_X_BPUNGAPSENT`, `LOSAT_X_BPUNGAP`, `LOSAT_X_PAIRWISEPAR` and `LOSAT_X_BPCOMPOREUSE` are removed
  (`LOSAT_X_BPGAPSHADOW`, `LOSAT_X_BPREDOPARSHADOW`, `LOSAT_X_BPREDOWINEXPANDSHADOW`, `LOSAT_X_BPSEEDSHADOW`,
  `LOSAT_X_BPUNGAPSENTSHADOW`, `LOSAT_X_BPUNGAPSHADOW`, `LOSAT_X_PAIRWISEPARSHADOW` and
  `LOSAT_X_BPCOMPOREUSESHADOW` are ignored, and the `[X_SHADOW]` lines of these switches are no longer
  printed). The switches stay in the strict set with their `=0` paths; a shadow variable no longer turns its
  switch on, and `LOSAT_X_BPSEEDSHADOW` no longer turns the fused loop on for `LOSAT_X_BPDIAG` and the
  switches that act in it. With `LOSAT_X_STATS`, `[X_STATS] BPHITBF subjects=` counts the subjects the
  branch-free loop scanned (with the `xstats` feature, only when not zero); `[X_STATS] BPGAP` stays. Outputs do
  not change. The comparisons that ran only under a shadow are now unit tests (the redo windows computed
  largest first against the `par_iter` windows, the windows that expand compact lists ahead against the full
  expansion) or are dropped where the reports are compared with the `=0` path and the reach rows or
  `[TIMING]` lines show that the path ran.
