- The variables `LOSAT_X_TBXFRAMES` and `LOSAT_X_TBXUNGAP` are gone: the TBLASTX subject frames made in
  windows and freed after the seed stage (the frames of a subject's letters only when HSPs reach the
  re-evaluation) and the byte score table of the two-hit ungapped extension were already the default and now
  run unconditionally. The whole-sequence frame making that `LOSAT_X_TBXFRAMES=0` ran was deleted (its
  helpers stay as the reference of the unit tests); the calls the byte table does not cover (an X-drop outside
  1..=127, a word that does not fit) and the `debug_ext` output still run the reference two-hit extension, now
  a rule of the table walks. `LOSAT_X_TBXNEXT`'s test no longer has a `LOSAT_X_TBXFRAMES=0` arm. Outputs do not
  change.
