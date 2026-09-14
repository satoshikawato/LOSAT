# Independent integrated task controls

The required read-only `ncbi_parity_auditor` supports all nine raw-output, median nonregression and maximum-RSS gates. No additional 5% speed floor is imposed on controls.

All 135 records (nine official oracle, 18 diagnostic, 18 warmup, 90 measured) were independently verified, with zero exclusions and exact raw bytes. All five pairs per condition contribute to the saved medians, ranges and maximum RSS. The table is retained in `integrated-performance-progress.md`; original checks and all samples remain under `integrated-np-task-controls/`.

The auditor verified actual task/outfmt 6 argv, 14 fixture hashes, official NCBI 2.17.0+ binary hashes, 154 production-source hashes, Cargo inputs and 49 runner hashes. Actual B1 (`20d604bb…`) and integrated (`5b4f4283…`) Wasm hashes match, both using Node 24.21.0 without extra flags and the same actual B1-artifacts runners. No db_gencode exception applies. All nine newly produced oracle bytes match the already audited component-control oracles.

Every one of the 18 diagnostics has the requested n8 pool. Actual stderr confirms worker IDs 1–8 with spawn_attempt→spawned→ready→exited(code=0). All argv/result/usage records, AB/BA ordering and exact monotonic endpoints establish non-overlap.

All raw samples remain present, including the measured megablast single-run maxima: EDL933/Sakai baseline 4.232277 s and Sakai/MG1655 candidate 3.538064 s. Passing a median gate does not imply every sample became faster.

## Clock diagnostic observations

Seven records, including four measured records, have realtime_clock_agreement=false. All 135 records pass monotonic/boottime agreement; their maximum absolute difference is 9.457 microseconds. Under the existing declared clock policy in wasm_performance.py, realtime and GNU elapsed remain adjustable-clock diagnostics. No sample was excluded, corrected or retried. The cause of the realtime difference was not established by the audit.

Paths below are relative to `integrated-np-task-controls/`; each directory contains its actual result/raw files.

| Directory | Phase | Realtime minus monotonic (ms) |
|---|---|---:|
| diagnostic/PeseMJNV.PemoMJNVB.losatn.blastn/baseline-threaded-n8 | diagnostic | 568.995350 |
| cold/repeat-0/NZ_CP006932.NZ_CP006932.losatn.blastn-candidate-threaded-n8 | warmup | 478.255213 |
| cold/repeat-0/AP027131.NZ_CP006932.losatp-candidate-threaded-n8 | warmup | 497.086535 |
| cold/repeat-3/AP027131.NZ_CP006932.losatp-candidate-threaded-n8 | measured | 486.145143 |
| cold/repeat-3/NZ_CP006932.NZ_CP006932.losatn.blastn-baseline-threaded-n8 | measured | 486.083001 |
| cold/repeat-4/Sakai.MG1655.losatn.megablast-candidate-threaded-n8 | measured | 486.348228 |
| cold/repeat-5/EDL933.Sakai.losatn.megablast-baseline-threaded-n8 | measured | 543.899062 |

This is the completed nine-control scope only. Thread boundaries, measured reuse, capacity and formal certification remain separate evidence.
