# Independent C-P1 task-control and component selection audit

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-14. Selection for integration is supported; full G-P adoption remains pending final gates.

All45 actual task outputs equal the current oracle bytes:3 oracle,6 diagnostic and36 cold records, including30 measured samples. Each of3 conditions retains1 warmup and5 measured pairs. Exclusions are empty. All medians, ranges, maximum RSS, GNU time CPU/elapsed values, result records, argv and output sinks agree. Six n8 diagnostics show8 unique worker IDs and complete spawn_attempt→spawned→ready→exited sequences with pool/effective8. No recorded overlap occurs within or across the four C-P1 groups (163 total records,110 measured samples).

| Input | Baseline median s | Candidate median s | Improvement | Baseline max RSS bytes | Candidate max RSS bytes |
|---|---:|---:|---:|---:|---:|
| AP027131/NZ_CP006932 |11.237340664|9.865743037|12.205714%|324472832|340131840|
| WSSV/PajaWSV |1.274843040|1.199806882|5.885913%|277520384|282271744|
| SicyWSV/CoBV |0.771622380|0.716136002|7.190872%|275644416|275632128|

Time and RSS gates pass. The long input's maximum RSS rises4.83%; WSSV rises1.71%, both below the allowed bounds. The154-file source snapshot,49 saved harness files, build inputs,6 fixtures and actual Wasm/Node/shared-runner/oracle hashes match; shared identities equal the prior main group.

The component selection rests on both n8 main improvements11.76%/12.45% and all controls passing. Native n1 remains4.56%/3.43% slower. Reuse/linear memory, formats, n2/n4 and final combined source/artifacts are pending; frozen certification is separate.

Post-control logs show routing15 tests OK, evidence22 tests OK and bash syntax exit0. The candidate's focused release tests show20 passed/1 existing diagnostic ignored; final native/threaded builds succeeded. This local tooling scope does not certify the full normal or serial-compatibility workflow. The root-integrated P1 source was subsequently verified byte-identical to isolated source `ae3d1cc…`.
