- The variables `LOSAT_X_TBNHITBF`, `LOSAT_X_TBNLINKWIN` and `LOSAT_X_TBNSCANBF` are gone: the TBLASTN
  two-hit loop that takes the too-far and overlap outcomes without a branch, the uneven-gap linking that
  tests only the candidates of a head set that can be admissible, and the lookup scan without a branch on
  the presence vector were already the default and now run unconditionally. The call of NCBI's scan loop
  from the unit path and the `[X_STATS] TBNSCANBF` counts of it (`units_ref`, `calls_ref`; the line is now
  `units= calls=`) were deleted, and so was `LOSAT/tests/tblastn_hitbf.rs`, which compared the report with
  `LOSAT_X_TBNHITBF=0`; the other three tests of these switches lost their `=0` arms. `XTwoHit::pairs` stays
  as the body that the branch-free loop hands the other hits to, and NCBI's whole-window loop of the
  linking call stays as the reference of its unit test. Outputs do not change.
