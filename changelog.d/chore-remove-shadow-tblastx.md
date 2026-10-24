- The SHADOW modes of `LOSAT_X_TBXEXTBATCH`, `LOSAT_X_TBXHITBF`, `LOSAT_X_TBXKMERGE`, `LOSAT_X_TBXPIPE`,
  `LOSAT_X_TBXQFRAMES`, `LOSAT_X_TBXQSLICE`, `LOSAT_X_TBXRPTPAR` and `LOSAT_X_TBXUNGAP` are removed
  (`LOSAT_X_<NAME>SHADOW` is ignored, and the `[X_SHADOW]` lines of these switches are no longer printed). The
  switches stay in the strict set with their `=0` paths. The `[X_STATS] TBXEXTBATCH` and `[X_STATS] TBXQSLICE`
  lines, which only the shadows printed, are gone with their counters. Outputs do not change. The pooled
  `engine_counts` path of outfmt 8 and 12 (`LOSAT_X_TBXRPTPAR`) is now compared with `=0` in
  `tblastx_report_par`; the other checks that ran only under a shadow are the unit tests of the live functions.
