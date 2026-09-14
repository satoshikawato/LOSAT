# Independent C-X2 exploration audit

Read-only ncbi_parity_auditor recomputed all33 actual raw files against their three official NCBI2.17.0+ oracles. All match. All per-result/usage fields, six summaries, medians/ranges and maximum RSS reproduce exactly. One warmup and three alternating pairs produce18 timed samples, zero exclusions and no recorded overlaps.

All30 LOSAT invocations use the same B1-artifacts runner path and Node --no-liftoff --no-wasm-tier-up. Candidate Wasm SHA prefix8d59fdc matches its build output. Six diagnostic logs contain matching eight unique worker IDs through successful exit. Main medians change+0.84%,+0.02%,−0.07%; all declared time/RSS nonregression bounds pass. This supports further qualification, not a whole-search speedup claim or final adoption.
