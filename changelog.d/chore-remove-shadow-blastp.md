- The SHADOW modes of `LOSAT_X_BPCLUT`, `LOSAT_X_BPLKPFLAT`, `LOSAT_X_BPLUT`, `LOSAT_X_BPLKPAHEAD` and
  `LOSAT_X_TBXLUT` are removed (`LOSAT_X_BPCLUTSHADOW`, `LOSAT_X_BPLKPFLATSHADOW`, `LOSAT_X_BPLUTSHADOW`,
  `LOSAT_X_BPLKPAHEADSHADOW` and `LOSAT_X_TBXLUTSHADOW` are ignored, and the `[X_SHADOW]` lines of these switches
  are no longer printed). The switches stay in the strict set with their `=0` paths. With `LOSAT_X_STATS`, the
  live counters that the shadow lines showed are printed as `[X_STATS] BPLKPFLAT tables= recycled=
  recycled_hits=` and `[X_STATS] BPLKPAHEAD blocks= prebuilt= dropped=` (only when a counter is not zero).
  Outputs do not change. The check that a block's values built ahead (query frames, table, contexts) are the
  values of its own build is now a unit test in `x_lookup_ahead.rs`; the other checks that ran only under a
  shadow are the unit tests that compare the live builders with the reference builders.
