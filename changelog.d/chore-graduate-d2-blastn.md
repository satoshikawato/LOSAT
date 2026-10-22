- The variables `LOSAT_X_GREEDYFAST`, `LOSAT_X_GREEDYFUSE`, `LOSAT_X_GREEDYPRELIM`, `LOSAT_X_ITREEFAST`,
  `LOSAT_X_PURGEFAST` and `LOSAT_X_CTXFAST` are gone, with `LOSAT_X_GREEDYSHADOW`, `LOSAT_X_GREEDYFUSESHADOW`,
  `LOSAT_X_GREEDYPRELIMSHADOW`, `LOSAT_X_ITREESHADOW` and `LOSAT_X_CTXSHADOW`: the megablast greedy rows, the
  fused score-only greedy loop, the one-byte subject copy of the preliminary stage, the interval-tree side
  indexes, the purge without per-removal shifts and the query-offset-to-context runs were already the default
  and now run unconditionally. The NCBI-literal shifting loops of the endpoint purge, the shadow comparisons,
  their `[X_SHADOW]` / `[X_STATS]` lines and counters, the four `#[cfg(test)]` mode hooks and the unit tests
  that compared against them were deleted; the ported loops the default still takes stay (the distance loop
  for packed subjects and traceback calls, the tree walk while tracing is on, the per-base context table for
  overlapping contexts and the `LOSAT_X_CTXO1SHADOW` comparison). `tests/blastn_greedy_fuse.rs`, which only
  tested the fused loop's variables, was deleted. Outputs do not change.
- The variables `LOSAT_X_GREEDYSPEC`, `LOSAT_X_AHEAD`, `LOSAT_X_MBDELAY` and `LOSAT_X_MBBATCH` are gone, with
  `LOSAT_X_AHEADSHADOW`: the look-ahead of the ordered gapped loops (window 16, now the constant
  `xahead::WINDOW`), the speculative megablast traceback and the delayed (x86_64) or batched (other targets,
  64 words, now `X_MB_BATCH`) megablast lookup fill keep their rules (a parallel search without
  `LOSAT_DEBUG_COORDS`, the owner loop of `LOSAT_X_QGROUPPAR`, the CPU) and move from the strict set to
  `GRADUATED`. The unqueued word loop of `build_mb_lookup`, which nothing reached any more, and the AHEAD
  shadow comparison were deleted. Outputs do not change.
