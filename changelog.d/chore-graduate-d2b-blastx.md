- The variables `LOSAT_X_BXLAZYCTX`, `LOSAT_X_BXLEAN`, `LOSAT_X_BXLUT`, `LOSAT_X_BXSCAN`,
  `LOSAT_X_BXSEGMEMO`, `LOSAT_X_BXPAR`, `LOSAT_X_BXPOOL`, `LOSAT_X_BXBATCH` and `LOSAT_X_BXCHUNK` are gone, with
  `LOSAT_X_BXLUTSHADOW`, `LOSAT_X_BXSCANSHADOW`, `LOSAT_X_BXSEGMEMOSHADOW` and the tuning variable
  `LOSAT_X_BXWAVE` (it acted only without `BXCHUNK`): the BLASTX paths they chose were already the default
  and now run unconditionally, so neither `=0` nor `LOSAT_X_DEFAULT=0` reaches the previous code. The context
  lookup after the diagonal test, the empty-initial-HSP return and the passive-observer skip, the in-place
  lookup table, the direct seed loop, the memoised redo SEG, the three shadow comparisons with their
  `[X_SHADOW]` lines (BXLUT, BXSCAN, BXHITBF's, BXSEGMEMO), the locked-stderr and per-batch fallbacks of
  `BXBATCH` / `BXPOOL`, and the one-subject-per-worker wave went; the reference `word_finder` (diagnostic
  stages and tests), the chain builder `build_lookup_from_prepared` (diagnostics mode, wide backbone
  indices, the tblastx and blastp tables) and the serial stages that 1 thread runs stay. `BXPAR`, `BXPOOL`,
  `BXBATCH` and `BXCHUNK` keep their rules (2+ threads, several query batches, waves of 64 subjects per
  worker) and move from the strict set to `GRADUATED`; the other five leave the set. Outputs do not change.
