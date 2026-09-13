# Independent read-only P2 baseline audit

The repository-required `ncbi_parity_auditor` reviewed the completed evidence without running builds or benchmarks. It re-read all 41 outputs, recomputed hashes, compared actual bytes with each fixture oracle, and independently recomputed counters/capacity/timing sums from raw diagnostic lines. All passed. The 21 threaded/instrumented records have the requested lifecycle counts, ordered lifecycle events, and zero worker exit codes. The parser yields nine PASS diagnostic records and two expected n1 INACTIVE records.

Verified: AP027202 has 4,214 batches, 59,119 computed, 59,087 consumed, 32 discarded; its nine zero-job and six one-job batches total 2.042171 ms at n8. LC738874 has 186 batches, 2,761 computed, 2,751 consumed, ten discarded and no zero/one-job batches. NZ self retains 572,262,953 bytes (545.752 MiB) of speculative scratch capacity, which is not RSS. LC738873 has 48 discarded precomputations and one ordered-replay fallback traceback.

The isolated candidate changes only threaded-Wasm batch16 to batch32; original-index storage and sequential containment/materialization/insertion remain intact. Batch counts do not prove synchronization dominates. A larger window changes extra speculative work and scratch-slot assignment; scratch must remain independent of prior alignment lengths/bytes. NCBI common-endpoint replacement (`blast_itree.c:273–291`) can make a previously contained HSP eligible, so the fallback must remain.

Before any acceptance, the auditor requires candidate raw parity/lifecycle across thread counts, full/partial windows, equal/contained HSPs, negative strands, multiple query/subject shapes, short/long alignments and repeated/error paths. LC738873 fallback and NZ self memory/trap coverage are explicit controls, alongside megablast. Adoption additionally needs ordinary release A/B repetitions, regression and memory gates, and separate candidate discarded-work/capacity diagnostics. Diagnostic process times include formatting and stderr overhead. Reserve/reallocation/physical-copy counts remain unmeasured.

This audit approves the baseline evidence, not the candidate or a release certification.
