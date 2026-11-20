- The variables `LOSAT_X_HSPINPLACE`, `LOSAT_X_LINKSORT`, `LOSAT_X_LINKSLICES`, `LOSAT_X_LUTARENA` and
  `LOSAT_X_SEEDRELEASE` are gone, with `LOSAT_X_HSPINPLACESHADOW`, `LOSAT_X_LINKSORTSHADOW` and
  `LOSAT_X_LINKSLICESSHADOW`: the TBLASTX HSP list sorted and linked without a second copy, the one final
  sort and the parallel index sorts of the linking call, the linking groups linked as slices of the one sorted
  list, the protein lookup-table chains carved from an arena and the early release of the seed-stage tables
  were already the default and now run unconditionally. The reference index-replay sort copy, the cut of the
  sorted list into one Vec per group with its assembly and serial loop, the `Vec` lookup chains, the off state
  of the release holders, the shadow comparisons, their `[X_SHADOW]` lines and counters and the unit tests that
  compared against them were deleted. The reference linking body that `LOSAT_DEBUG_CUTOFFS` and
  `LOSAT_DIAGNOSTICS` select stays, and so do the `[X_STATS]` lines of the paths that remain (printed whenever
  `LOSAT_X_STATS` is set). Outputs do not change.
- The variables `LOSAT_X_PAIRPAR`, `LOSAT_X_PAIRLUT` and `LOSAT_X_SEEDBUCKET` are gone, with
  `LOSAT_X_PAIRPARSHADOW`, `LOSAT_X_PAIRLUTSHADOW` and `LOSAT_X_SEEDBUCKETSHADOW`: the TBLASTX unit path
  and its direct group tables, the bucketed two-hit order for large diagonal tables, and the blastn and
  megablast one-hit strips and pool lookup fill all run by their rules without a variable. The rules stay
  (`LOSAT_X_PAIRPARGATE`'s cost rule, the gates and guard of the unit path, the 2^21-cell minimum of the
  bucketed order with `LOSAT_X_SEEDBUCKET_MIN_CELLS`, the pool and query-block rules of the strips) and so do
  their fallbacks: for TBLASTX the reference frame loop, `scan_chunk_plain` and `build_groups`, for blastn
  the queued fill `x_fill_mb_queued` and the reference seed scan. The serial table fill without a pool stays
  as the harness of the unit tests that compare live fills. The
  shadow comparisons, the `[X_SHADOW] PAIRPAR`, `PAIRLUT`, `SEEDBUCKET`, `TBXUNITMASK` and `MBLUTONCE` lines
  and the `[X_PAIRPARSHADOW] blastn:` line, the one-hit stash path of the blastn strips, and the tests that
  existed only for the variables were deleted. `LOSAT_X_DEFAULT=0` no longer turns these paths off: the
  `LOSAT_X_PAIRPARGATE` cost rule goes off with it, so every batch the gates and the guard let through takes the
  unit path. Outputs do not change.
