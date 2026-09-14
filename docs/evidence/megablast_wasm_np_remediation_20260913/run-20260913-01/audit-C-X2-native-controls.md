# Independent C-X2 native controls audit

Read-only ncbi_parity_auditor verified all 87 actual outputs against NCBI 2.17.0+ and the earlier threaded outputs. All are byte-identical under outfmt 6, gencode 1/1, with no exception. The run has three oracles, twelve diagnostics and 72 cold calls: six conditions, one warmup plus five alternating pairs, 60 measured calls, zero exclusions and no recorded overlapping execution.

| Fixture | Native n1 B1 → C-X2 median (s) | Native n8 B1 → C-X2 median (s) |
|---|---:|---:|
| MelaMJNV/PemoMJNVA | 3.753476 → 3.736233 (+0.46%) | 2.246867 → 2.231921 (+0.67%) |
| MjeNMV/MelaMJNV | 17.227884 → 17.076840 (+0.88%) | 10.459179 → 10.457215 (+0.02%) |
| AP027280 self | 31.250408 → 31.422908 (−0.55%) | 19.006661 → 18.686575 (+1.68%) |

All sample/result/usage fields, medians, ranges, maximum RSS and engineering tolerances reproduce exactly. Every time/RSS condition passes. Native candidate 994837f… and B1 ae494b43… match actual artifacts; source 154 files, Cargo inputs and saved tooling 49 files match their hashes. The argv calls native executables directly, differing only in executable and output destination. Diagnostics show pool 0 and no parallel selection for n1, pool 8 for n8.

The auditor found a descriptive metadata defect: an old fixed boundary sentence mentions Node/Wasm costs even for native runs. Actual argv and measured values are native process launch-through-exit. Historical metadata is retained and corrected in interpretation by native-timing-boundary-amendment.json. The current harness now records explicit boundary_by_kind text; its timing code is unchanged.

These findings support native nonregression only. Threaded n1, remaining boundaries/genetic-code cases, and the final integrated combination remain separate gates. No whole-search TBLASTX speedup is claimed.
