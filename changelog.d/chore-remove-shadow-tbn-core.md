- The SHADOW modes of `LOSAT_X_TBNHITBF`, `LOSAT_X_TBNKEEPFEW`, `LOSAT_X_TBNLINKWIN`, `LOSAT_X_TBNSCANBYTE`,
  `LOSAT_X_TBNUNITXLATE`, `LOSAT_X_ERFABSORB`, `LOSAT_X_SEGTRIMMEMO` and `LOSAT_X_DBIDXREAD` are removed
  (`LOSAT_X_<NAME>SHADOW` is ignored, and the `[X_SHADOW]` lines of these switches are no longer printed). The
  switches stay in the strict set with their `=0` paths; the `[X_STATS] TBNLINKWIN` line loses `shadow_heads=`.
  Outputs do not change. The checks that ran only under a shadow are unit tests of the live functions against the
  `=0` ones.
