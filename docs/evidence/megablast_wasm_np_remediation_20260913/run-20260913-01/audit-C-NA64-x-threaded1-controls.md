# Independent C-NA64 TBLASTX threaded n1 review

Read-only review by `ncbi_parity_auditor`. All 45 records (3 official oracles, 6 diagnostic, 6 warmup, 30 timed) and 78,218,400 raw output bytes match the corresponding oracle exactly; zero exclusions. All three five-pair conditions pass the original time and RSS allowances.

| Input | Baseline median (s) | Candidate median (s) | Reduction (%) | Baseline / candidate max RSS (bytes) |
|---|---:|---:|---:|---:|
| MelaMJNV/PemoMJNVA | 4.900043522 | 4.943301525 | -0.8828085 | 180326400 / 181665792 |
| MjeNMV/MelaMJNV | 19.301373712 | 19.384823411 | -0.4323511 | 185716736 / 186605568 |
| AP027280 self | 38.063448156 | 37.286719928 | 2.0406145 | 197435392 / 199954432 |

All ranges overlap. All six diagnostics request/effectively use one thread, with pool size zero and no host workers. Actual B1 `20d604bb…` and candidate `5b4f4283…` Wasm artifacts, common B1 runner paths, Node 24.21.0 and the two TBLASTX Node flags match the records. Source manifests bind 589 files in each of three trees. Actual argv, usage, repeated order, raw bytes and non-overlap were checked.

Minimum process gap: 0.045276 s; next boundary group starts 10.137373 s later. Realtime-clock diagnostic disagreement occurs in 14 records (8 timed), retained without exclusion or correction. Maximum realtime delta: diagnostic AP candidate +2.298264 s; timed AP baseline repetition 5 +2.211111 s. Maximum boottime difference: 32.273 microseconds. CLOCK_MONOTONIC remains the elapsed-time authority.

This supports a 2.04% median reduction for AP027280 at threaded n1, not a broad or statistically established TBLASTX speedup.
