# C-X2 threaded n1 and n2/n4 independent review

Recorded from the read-only `ncbi_parity_auditor` review by `/root/audit_megablast` on 2026-09-13 UTC. The reviewer inspected actual saved raw files, not only the declared PASS fields.

Threaded n1 has 45 records, including 30 measured runs. All raw files equal the fresh oracle; medians, ranges and maximum RSS recompute exactly from one warmup and five alternating repetitions. Pool 0, effective one thread, nonparallel stages and no host workers are verified.

| Main case | Baseline median (s) | Candidate median (s) | Reduction | Time/RSS bounds |
|---|---:|---:|---:|---|
| MelaMJNV.PemoMJNVA.tlosatx | 5.590681 | 5.712184 | -2.1733% | PASS/PASS |
| MjeNMV.MelaMJNV.tlosatx | 21.057822 | 21.311891 | -1.2065% | PASS/PASS |
| AP027280.AP027280.tlosatx | 41.328267 | 41.040117 | 0.6972% | PASS/PASS |

The n2/n4 boundary group has 15 records, including 12 diagnostics and no timed samples. All raw files match; requested pools and matching unique worker IDs are verified through spawn, readiness and exit code zero. These records make no performance claim.

Both groups match the declared 154 Rust sources, Cargo inputs, 49 tooling snapshots, actual artifacts, Node flags and common runner paths. The reviewer found no recorded process overlaps. C-X2 remains a final-sort data-movement reduction; no whole-search speedup is claimed. Long gencode4 and integrated gates were still pending at this review.
