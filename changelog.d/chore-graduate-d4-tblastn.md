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
- The variables `LOSAT_X_TBNSSIDE`, `LOSAT_X_TBNBATCH`, `LOSAT_X_TBNOIDCHUNK`,
  `LOSAT_X_TBNOIDNEXT`, `LOSAT_X_TBNOIDUNITS`, `LOSAT_X_SEGSHARE` (with `LOSAT_X_SEGSHARESHADOW`)
  and `LOSAT_X_DEKKER` are gone; these conditional and platform switches keep their rules (four or
  more query batches for the shared subject plans and the side-by-side batches, chunks of subjects
  taken by free threads and units of a long subject on the pool at 2+ threads, the shared SEG store
  of the parallel blastp redo, the Dekker product on wasm32) and their `x_reach!` marks. Gone: the
  entries of chunks and of one subject per thread in the tblastn subject loop
  (`x_chunks_per_entry`, the `pool_entries` counter; the line now prints `pool_entries=1`), the
  gate that sent a one-thread search to the batch loop only with `TBNSSIDE`, the SEGSHARE shadow
  comparison and its counter, and the tests that set the removed variables to `0`. Outputs do not
  change.
