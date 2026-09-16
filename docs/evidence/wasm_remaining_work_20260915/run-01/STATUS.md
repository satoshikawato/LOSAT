# Completion status

**Completed: H1a is applied and accepted. X1 is rejected because it is slower.**
X2, M1 and H1b did not meet their implementation eligibility conditions.
The Wasm/Native 1.20 ratio is a goal, not an adoption prerequisite, as confirmed
by the user during execution.

See [ADOPTION.md](ADOPTION.md) for measurements and decisions,
[DIAGNOSIS.md](DIAGNOSIS.md) for the n8 investigation, and
[REUSE.md](REUSE.md) for bounded module/reactor reuse results.

| Stage | Outcome |
|---|---|
| B0 | Frozen dirty-tree inputs; fresh-oracle raw parity 75/75 |
| D0 | Same group work at n1/n2/n4/n8; actual V8 code inspected; TurboFan diagnostic removes reversal on both primary inputs |
| X1 | Isolated implementation; state/trace/stage/long-code4 comparisons pass; all 4 default body medians worsen; Rust change rejected |
| X2 | Capacity observed; dominant physical-memory/allocation cost not established; deferred |
| M1 | Only 0.070–0.134% of uncompressed calls reach 16 matches; no SIMD implementation adopted |
| H1a | Five host/test files applied; serial preparation 7.541→4.700ms; primary complete-process 51.205→50.399ms |
| H1b | Raw inspection/guard/guarded compile measured; retain guard and defer parser changes |
| V0 | Final serial 15/15 raw; all four artifact kinds; reactor recovery; reuse 192/192 and all 24 session/case guards pass |

## Validation

- Standard `cargo test`: 625 passed, 0 failed, 3 existing ignored tests.
- `cargo clippy --release`:PASS; `cargo fmt --check`:PASS.
- Node host harness: 19 passed outside the sandbox after both old/new were blocked
  by sandbox `spawnSync EPERM`. No code workaround was added for that restriction.
- Final five host-file hashes match the measured/tested H1 candidate. All 158
  Rust/configuration/lock/build inputs remain identical to the initial baseline.
- Native/serial/threaded baseline artifact hashes remain unchanged after checks.
- Independent `ncbi_parity_auditor` review reconciled final files, raw oracle
  outputs, all 192 reuse outputs and all 24 session/case guards.

The initial release-test attempt omitted one fixture from the isolated snapshot;
its unchanged original was restored. The next release-test attempt exposed the
baseline rlib/cdylib output collision and abort/unwind panic-strategy conflict.
Both failures are retained. The documented standard `cargo test` then passed in
a separate target directory. No Rust profile or production code was changed.

## Measurement scope

Main performance runs used ordinary Node, the same filesystem and identical
artifacts/runners where applicable, with one warmup per version and three
alternating A/B pairs. Timing runs had no concurrent builds, correctness tests,
other benchmarks or compression. Matched coarse entry/exit clocks **were used**
for Native/Wasm body medians; inner-loop counters and stage diagnostics were
excluded from adoption samples. TurboFan-fixed single searches are supplementary
diagnostics and cannot establish a default-runtime speedup.

Retained threaded times are baseline references for unchanged Rust, not a fresh
measurement of the final shared host inspector. Module reuse uses fresh instances;
same-instance final close is recorded outside invocation samples. No unlimited
repeat or reactor-memory repair claim is made.

## Remaining scope

Normal-Node TBLASTX remains above the 1.20 reference goal. The existing invalid
TBLASTX-query status difference and reactor-memory issues remain unresolved.
PR5 Gate A / official-platform Gate B certification, browser timing/QA, Wasmtime
and wasm32-unknown-unknown timing were not performed in this host-only adoption.
Browser release artifacts were not changed or published. The prior adoption,
expected files and registered platform fingerprints remain unchanged.

## Handoff

Suggested commit title: **Reuse validated Wasm modules in serial WASI hosts**

Reuse the inspected serial Module while preserving identity JSON, ABI rejection,
instance ownership and exit status. Keep the independently evaluated regression
and performance evidence with the change. No commit, push, tag or publication
was performed.
