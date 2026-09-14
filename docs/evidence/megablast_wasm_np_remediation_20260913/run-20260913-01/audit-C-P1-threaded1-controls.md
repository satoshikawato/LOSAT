# Independent C-P1 threaded n1 control audit

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-14. No inconsistencies found. This supports output parity and speed improvement on the two declared inputs and this runtime condition.

| Input | Baseline median s | Candidate median s | Improvement | Baseline max RSS bytes | Candidate max RSS bytes |
|---|---:|---:|---:|---:|---:|
| AP027078/AP027131 |36.894490|31.730372|14.00%|165392384|165150720|
| AP027132/NZ_CP006932 |55.970711|47.786443|14.62%|177381376|176766976|

All30 actual output files (2 oracle,4 diagnostic,24 cold) equal the official oracle and preceding main-confirmation raw bytes. Each condition has1 warmup and5 measured pairs,20 measured samples total, with no exclusions. Alternating order, median/range/maxRSS and GNU time entries agree; there is no recorded process overlap. Four diagnostic logs show requested1, pool0, effective1, zero worker spawns and every stage's parallel flag false.

The154-file candidate snapshot, build inputs, fixtures,49 saved harness files and actual artifact/Node/shared runner/oracle hashes agree. The benchmark harness itself matches main-confirmation. Candidate source is the unchanged `ae3d1cc…` traceback source and the threaded artifact is `f7f598fc…`, with full identities retained in metadata. Oracle: BLASTP2.17.0+, local query/subject, outfmt6, num_threads1. No biological exception applies.

Task controls were still running at review time. Reuse, linear memory, final integration and frozen certification are outside this audit.

The reviewer also verified the C-X1 test-only follow-up at source SHA `b6203cb9…`: the group execution body remains text-identical to B1; the added direct-group assertion covers fixed addresses, original-owner updates and subsequent-group scratch reuse. That candidate remains unbuilt and unaccepted. No tests, builds, benchmarks or edits were performed in this audit.
