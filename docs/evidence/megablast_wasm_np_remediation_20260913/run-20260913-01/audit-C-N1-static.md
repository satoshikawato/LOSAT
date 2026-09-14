# Independent C-N1 source review

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), 2026-09-13. Read-only; no edits/builds/measurements.

Reviewed candidate `purge_endpoints.rs` SHA256 `0ada014bf3e45568f81fe413a12fcbef1fd370652fb76ce9c7b2d75ac30f232b`, NCBI `blast_hits.c:2496–2499,2522–2525`, and caller sequence `run.rs:9704 → 9752 → 9938`.

Both rotations preserve every slot of the previous Rust take-loop and NCBI pointer shift. After removal the slice is `[None,A,…,Z]`; rotation produces `[A,…,Z,None]`, then the same tail HSP is installed. The inner guard and decrement guarantee a nonempty slice, including a last-active-slot removal. Existing tail entries lie outside that range. Repeated removal preserves active order and reverse-removal tail order. This holds for trimmed Some and deleted None with purge=false and for purge=true. The second sort covers only the active prefix; final flatten preserves live order and the same extra_start. No comparison, coordinate, precision, HSP content or parallel reduction changed.

The new regression exercises trim and delete in both endpoint passes, active count three, tail order after interspersed NULLs, and edit scripts. NCBI `unit_tests/api/blasthits_unit_test.cpp:1163` invokes BLASTP (effectively purge=true); the new test adds nucleotide trimmed-tail coverage. Explicit common-end last-slot deletion and purge=false last-slot trimming are not separate fixtures, but the same static proof covers them; this is not a blocker. Existing score-tie coverage reaches common-start/purge=true last-slot removal.

The focused release log contains ten passing tests: six BLASTN and four TBLASTX. The reviewer rehashed all sixteen fair-exploration cold outputs; both inputs match across baseline/candidate. This review supports source equivalence, not performance adoption. Five-pair controls and final target gates were still running/pending.

## Completed control audit

The same independent reviewer subsequently recalculated all fourteen control configurations from samples: main-controls 86 process records / six configurations, native8-controls 30 / two, task-controls 90 / six, each with one excluded warmup and five measured pairs. All medians, full sample values, ranges, maximum measured RSS, deltas and plan tolerances agree. All time/RSS nonregression checks pass.

The reviewer compared actual bytes of all206 control outputs and22 fair-exploration outputs with their saved official oracles, with no mismatch or exclusion; verified alternating order and non-overlap from monotonic timestamps; and checked normalized argv, Node flags, thread counts, runner paths and actual artifact/runner/fixture hashes. All four runs have identical candidate source hashes. The reviewer supports only the component-limited disposition in C-N1-disposition.json, not G-N completion. Stale running-status text in REPORT/STATUS/performance-investigation was corrected following the review.
