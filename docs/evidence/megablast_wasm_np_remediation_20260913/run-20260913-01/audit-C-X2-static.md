# Independent C-X2 source audit

Read-only `ncbi_parity_auditor` reviewed the actual isolated source, saved patch/identity and focused test log. Linking SHA-256: `879a36c15c0bb4747562ba5829466e5f304b7c9b9d0445894e53845af3bf976a`. Of154 Rust files, only linking.rs differs from B1; Cargo inputs match B1. Twelve focused tests passed.

The initial owned sort, frame grouping, complete owned group kernel and indexed parallel reduction remain unchanged. Final singleton handling, stable reverse/forward reference sorting and replay preserve B1 semantics, including duplicate visits. NCBI link_hsps.c:990–994,1080–1085 supports pointer sorting and replay order. No concrete ownership/ordering defect was found. Stable ties mean preservation of B1 behavior, not a claim that C qsort is stable.

The auditor noted a cosmetic test reference comment attached to the preceding closing brace at line2867. It will be separated at integration and explicitly recorded as a comment-only difference; final tests/builds will use the integrated source. This source audit does not certify performance, real fixtures, or the integrated combination.
