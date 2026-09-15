# Adoption decision

Status: **ADOPTED by explicit user instruction on 2026-09-15.**

User instruction: 「じゃあいいじゃん。採用」

The adopted implementation is the measured integrated version: C0's M0/N0
correctness corrections plus M1, N1, N2, N3, P1, P2 and X2. Its source and
artifact identities are recorded in `final-build-bindings.json`. The corrected
reuse-harness child-count assertion is also retained. Conditional candidates
previously skipped remain skipped.

The user accepted the implementation after reviewing the three-pair results,
the six time-guard failures, the improvements on long-running searches, and
the AP027132/NZ_CP006932 BLASTP native/Wasm comparison. This decision accepts
the observed performance tradeoffs. It does not change any measured value or
turn a failed or unrun check into a PASS.

- Final measured cold results: time guards 39/45, RSS guards 45/45; all 270
  included timed outputs match. See `integrated-three-pair-results.md`.
- Existing validation: 625 Rust tests passed / 3 ignored; clippy and format
  checks passed. The independent source/correctness/performance review is
  retained in `INDEPENDENT_REVIEW.md`.
- Measurements remain stopped at exactly the first three complete A/B pairs.
  Final control timings, body/reuse timing and memory windows, and secondary
  TurboFan timings remain unrun. No further benchmarking is required for this
  adoption decision.
- Existing parity and reactor/reuse issues retain their documented scope.
  This decision does not certify the v0.1.0 release or amend PR5/PR6 contracts.

This decision supersedes the pending-adoption assessment in the historical
`measurement-stop.json`, its copy in `measurement-plan.json`, and the frozen
execution archive. Those records and the auditor's nonregression findings
remain unchanged. No production code changed to record this decision.
