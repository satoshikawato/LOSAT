- The variables `LOSAT_X_TBNPAR`, `LOSAT_X_TBNQSIDE`, `LOSAT_X_TBNSTATE`, `LOSAT_X_TBNXLATE` and
  `LOSAT_X_LUTSPLIT` are gone, with `LOSAT_X_TBNPARSHADOW` and `LOSAT_X_TBNXLATESHADOW`: the
  unit-wise preliminary stage and parallel redo of TBLASTN, the query frames and lookup table built
  once per query set, the scan state kept by a thread across subjects, the six subject frames read
  from codon tables and the table-only / contexts-only lookup builds were already the default and
  now run unconditionally, so neither `=0` nor `LOSAT_X_DEFAULT=0` reaches the previous code. Gone:
  the serial preliminary stage in the search (it stays as the oracle of unit tests), the
  per-subject query frames and lookup (`XQuerySide` no longer optional), the fresh diagonal table
  per subject, the letter-wise subject translation of the search, both shadow comparisons with
  their `[X_SHADOW]` / `[X_STATS] TBNXLATE` lines, the `resolved_made` counter of
  `[X_STATS] TBNSTORE`, the test file `tblastn_subject_translate.rs` and the tests that compared
  with the previous code. Outputs do not change.
