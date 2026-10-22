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
