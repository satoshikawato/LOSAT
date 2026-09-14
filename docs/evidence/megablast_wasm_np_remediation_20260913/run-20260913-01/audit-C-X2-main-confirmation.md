# Independent C-X2 main confirmation audit

Read-only ncbi_parity_auditor verified all 45 actual outputs against their official NCBI 2.17.0+ oracles. All match byte for byte under outfmt 6 and gencode 1/1; no exception is used. One warmup plus five alternating pairs gives 30 measured samples, zero exclusions and no recorded overlapping executions.

| Fixture | B1 median (s) | C-X2 median (s) | Change | B1 maximum RSS (B) | C-X2 maximum RSS (B) |
|---|---:|---:|---:|---:|---:|
| MelaMJNV/PemoMJNVA | 4.144775 | 4.182248 | −0.90% | 331300864 | 344023040 |
| MjeNMV/MelaMJNV | 13.546580 | 13.367785 | +1.32% | 345059328 | 340176896 |
| AP027280 self | 24.980909 | 25.083464 | −0.41% | 359641088 | 353034240 |

All sample/result/usage records, ranges, medians, maximum RSS and engineering tolerances reproduce exactly. All three time/RSS nonregression conditions pass. Source 154 files, Cargo inputs and runner snapshot 49 files match their hashes. Both versions use identical runner paths and Node flags. Candidate Wasm 8d59fdc… matches the built artifact. All six diagnostic runs contain eight identical worker IDs through successful exit.

This supports only the declared main-condition parity/nonregression result. It does not show whole-search speedup, substitute for native/threaded n1 controls, or certify final integrated adoption.

The auditor also confirmed the capacity diagnostic draft corrections: exact NCBI 666–669 citation, AtomicU64 cumulative counters, and explicit null/N/A interpretation for unmeasured baseline fields. BLASTN task-specific outfmt 7 comparisons pass separate task arguments and output paths; cmdline_flags.cpp:69 was verified. These are static design findings, not executed capacity/final-format results.
