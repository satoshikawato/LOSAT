# Independent audit scope

The custom read-only `ncbi_parity_auditor` (`/root/audit_measurement`) reviewed this task. It did not build, run tests or collect performance samples.

- P1 adoption/control/local exploration: independently recomputed the stored summaries and read all150 successful actual outputs (71 adoption,62 control,17 valid local exploration). Oracle bytes and recorded hashes agree. General TurboFan adoption fails memory and other-task controls; the local candidate is below the10% threshold.
- P2 baseline: see `audit-p2-baseline.md`. All41 actual outputs, lifecycle and diagnostic counters were checked. Candidate invariant review confirms the same speculative DP set/replay semantics within the proposed batch change; the rejected exploration is not release adoption.
- P3 source: independently verified the contiguous-row invariant, first-cell fence behavior and NCBI ownership. The within-search SEG-cache idea is conditionally valid only after the query-dependent gate, but it was not implemented. Source analysis does not establish a measured10% improvement.
- P1 long follow-up: all6 processes, actual bytes/hashes, four n8 worker lifecycles and10 fresh code4 database files verified. All are untimed and have no adoption summary.
- P1 reuse follow-up: all13 processes and40 invocations (10 warmup/30 measured), job bindings, artifact identity, raw bytes, worker event order and successful exits, close and timing boundaries verified. All10 summary medians and memory/preparation values recompute. Thread order is[1,8,8,1,1,8,8,1]. One AB session with three measured invocations per condition is exploratory; no cold/adoption/leak/browser claim.
- P4: see `audit-p4.md`; all33 megablast outputs and the final158 source/four artifact identities independently checked. The real simple comparison's five output bytes also agree. Existing megablast failures remain PARTIAL.

The final measurement-file review found missing monotonic/boottime admission checks in two legacy consumers. The current source now rejects inconsistent cold/profile timing while preserving raw audit semantics, and marks warm invocation samples CLOCK_INCONSISTENT so no affected median is emitted. The added integration-style synthetic test exercises five consumer cases. Audited P1 observations are unaffected. The final read-only follow-up accepted both consumer gates and the five-case regression wiring with no further findings, and verified36 passing Python tests in `p0-clock-consumer-tests.log`.

This scope is not a complete formal release or platform GateA/GateB certification. Raw archive packaging is verified by local hash checks, not delegated performance tests.
