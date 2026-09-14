# Independent TBLASTX fixed-payload design review

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-14. Design support only; no build, execution, timing or adoption claim.

Keeping the incoming `Vec<UngappedHit>` fixed and moving exclusive references from `iter_mut()` through the three linking sorts and frame groups corresponds to NCBI `link_hsps.c:476–486,553–558,990–994,1012–1085`. The group function's input/output lifetime must be independent of context and thread-local pool borrows. Pools retain only values and indices. Comparators, stable ties, group-index parallel reduction, inner group indices, link-ID assignment and E-value update order must stay unchanged.

Replay must preserve its existing empty/singleton shortcut before non-head filtering, missing-ID handling, duplicate-ID last-write mapping, hop limit and repeated-index traversal. Keeping final payload clones makes repeated replay indices safe without duplicating mutable references. A singleton now needs one final clone where the old shortcut moved it; for N>1 the three sorts remove6N clones and frame transfers move references. Final replay clones remain.

Required focused coverage: all3 comparators and successive sorts with ties, fixed payload addresses, complete fields and E-value bits, serial/parallel frame boundaries, singleton and repeated replay. Current LC threshold and long gencode4 regressions remain; the two NCBI-required linking invocations are unchanged. Actual speed and memory require measurement.
