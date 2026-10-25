- The variables `LOSAT_X_ERFABSORB` and `LOSAT_X_SEGTRIMMEMO` are gone: Spouge's ErfC returning 2.0 for an
  argument below -6 without evaluating `ErfImpl` (all five programs; the value is the reference's, bound argued
  in `stats/spouge.rs`) and the per-thread memo of the SEG trim of a segment (BLASTX's query SEG and TBLASTN's
  redo windows) were already the default and now run unconditionally, so `=0` no longer returns to the
  previous path. The unit test of the memo asserts that the overlapping windows hit it. Outputs do not change.
