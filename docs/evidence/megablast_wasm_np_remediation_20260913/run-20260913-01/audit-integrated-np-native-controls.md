# Independent integrated native BLASTN/BLASTP controls

The required read-only `ncbi_parity_auditor` supports all eight declared median-time and RSS gates. No concrete blocker was found. Native BLASTP speedup is not supported.

All 116 records (four official oracle, 16 diagnostic, 16 warmup, 80 measured) were independently checked: actual raw bytes/SHA-256, result/usage/ordered argv, five-pair order, medians/ranges and maximum RSS. There are zero exclusions and every raw output matches exactly.

| Input | Threads | B1 median (s) | Integrated median (s) | Reduction |
|---|---:|---:|---:|---:|
| LC738874/LC738870 | 1 | 0.491106 | 0.497652 | −1.33% |
| LC738874/LC738870 | 8 | 0.350307 | 0.369669 | −5.53% |
| AP027202/LC738875 | 1 | 11.473342 | 9.075481 | 20.90% |
| AP027202/LC738875 | 8 | 6.859844 | 4.463187 | 34.94% |
| AP027078/AP027131 | 1 | 28.694368 | 30.000106 | −4.55% |
| AP027078/AP027131 | 8 | 5.403868 | 5.571826 | −3.11% |
| AP027132/NZ_CP006932 | 1 | 43.468891 | 45.382233 | −4.40% |
| AP027132/NZ_CP006932 | 8 | 8.358210 | 8.516648 | −1.90% |

The short LC n8 increase is 19.36 ms, within the unchanged 50 ms bound. Long native n1 BLASTP is 1.306 s / 1.913 s slower, leaving only 129 ms / 260 ms to the 5% bounds. Both native n1 BLASTP baseline/candidate ranges are disjoint: these data show a slowdown and must not be described as no difference.

The auditor verified blastn/task-blastn and blastp outfmt 6 argv, all eight fixture hashes, actual official NCBI 2.17.0+ binary hashes, the 154 frozen production-source hashes, Cargo inputs and 49 runner hashes. No db_gencode exception is used in this group. Native B1 (`ae494b43…`) and integrated (`f6e27376…`) actual artifact hashes match their records; full hashes remain in metadata.

Each repeat reverses the whole placement order as recorded; AB/BA and exact monotonic endpoints establish non-overlap. Diagnostic n1/n8 pool/stage records match actual stderr. All four oracle bytes also match the already audited integrated threaded-n8 main group. Native cold time includes startup, search, output and termination. Maximum RSS passes every declared bound.

This result is limited to these eight native conditions. Remaining controls, measured reuse, scratch-capacity and formal frozen certification are separate scopes.
