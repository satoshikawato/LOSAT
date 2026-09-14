# Independent C-NA64 TBLASTX n2/n4 boundary review

Read-only review by `ncbi_parity_auditor`. All 15 actual outputs, totalling 26,072,800 bytes, were checked: three fresh official oracles and exactly three inputs × baseline/candidate × n2/n4 diagnostics. Every output is byte-identical to its oracle; all records are untimed with zero exclusions.

The cases are MelaMJNV/PemoMJNVA, MjeNMV/MelaMJNV and AP027280 self, TBLASTX outfmt 6, query/db gencode 1, NCBI BLAST+ 2.17.0+. No genetic-code exception applies. Actual artifacts, fixtures, oracle, Node, runners and all command/result/usage files match their hashes: B1 `20d604bb…`, C-NA64 `5b4f4283…`, common B1 runner paths, Node 24.21.0 with `--no-liftoff --no-wasm-tier-up`.

All 12 diagnostics have the requested pool/effective size of 2 or 4 and no caller participation. All 36 workers and 144 lifecycle events have matching tid sets, ordered spawn_attempt → spawned → ready → exited per tid, and exit code zero. Source manifests bind 589 files × three trees; 154 production hashes and the fixed driver/plan/completed log match.

Processes do not overlap on CLOCK_MONOTONIC. Minimum internal gap: 0.056763 s; gap from n1 group: 10.137373 s; gap to reuse group: 9.767626 s. Seven untimed records have realtime diagnostic disagreement; maximum +1.224480 s for AP candidate n4. Maximum boottime difference: 9.041 microseconds. All records are retained without adjustment/exclusion.

An empty engineering-checks list is expected because no timed samples were requested. This verifies correctness and thread boundaries, not time/RSS adoption or the separate reuse group.
