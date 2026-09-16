# I4: TBLASTX n8 improvement

**Adopted and implemented in production. Final focused checks and independent evidence review passed.**

The ordinary-Node cold n8 slowdown is resolved on both primary fixtures. Mje n8 body time falls 51.254%; AP self falls 53.388%. Exact output is preserved. A 1.20 Native ratio is neither an adoption nor a completion condition; the decision uses improvement over the current baseline and the required correctness/non-regression checks.

## Change and cause

`LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs` is the only production change relative to the frozen baseline. Threaded WASI calls the large-gap predecessor scan once per scan-eligible HSP, providing later entries at which V8 can use completed optimized code. A borrowed backward prefix removes redundant helper addressing in the inspected optimized code. Native/plain WASI retain the original indexed loop in the original caller; one macro owns selection and tracing. Candidates, strict ties, arithmetic order, HSP updates, scheduling and result order are preserved.

NCBI authority: `c++/src/algo/blast/core/link_hsps.c:812–895`, particularly the predecessor loop at 827–863 and backward next_larger construction at 675–684/876–885. I1 PC/profile evidence in `DIAGNOSIS.md` supports the tier mechanism as a material contributor, not the sole cause of all elapsed time. I4 selected threaded group/helper WAT is byte-identical to I2b; whole artifacts differ. No profiling/compiler-forcing flags enter acceptance timings.

## Ordinary Node: fresh process

Node 26.8.2, V8 14.6.202.34-node.28, Rust 1.92.0, Intel i9-14900HX/WSL2, affinity 0–31. One warmup/version and fixed three AB/BA/AB pairs. Body covers identical Rust dispatch preparation/search/output boundaries; process time additionally includes startup and shutdown.

| Input | Threads | Baseline body s | I4 body s | Change | Baseline process s | I4 process s |
|---|---:|---:|---:|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | 8 | 22.223515 | 10.833122 | -51.254% | 22.400608 | 11.035249 |
| MjeNMV.MelaMJNV.tlosatx | 1 | 14.793121 | 14.100169 | -4.684% | 14.950732 | 14.226376 |
| AP027280.AP027280.tlosatx | 8 | 37.739457 | 17.590966 | -53.388% | 37.883066 | 17.736080 |
| AP027280.AP027280.tlosatx | 1 | 28.858794 | 26.316143 | -8.811% | 29.012081 | 26.461091 |

I4 n8 is 23.17% shorter than I4 n1 on Mje and 33.15% shorter on AP self. This cold-Node result is distinct from warmed browser/reuse results.

## Native and plain-WASI controls

| Input | Target | Threads | Baseline body s | I4 body s | Change |
|---|---|---:|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | native | 8 | 7.128258 | 7.277330 | +2.091% |
| MjeNMV.MelaMJNV.tlosatx | native | 1 | 12.183987 | 12.115818 | -0.559% |
| AP027280.AP027280.tlosatx | native | 8 | 13.746009 | 13.664597 | -0.592% |
| AP027280.AP027280.tlosatx | native | 1 | 22.966468 | 23.070087 | +0.451% |
| MjeNMV.MelaMJNV.tlosatx | serial | 1 | 13.890113 | 13.993189 | +0.742% |
| AP027280.AP027280.tlosatx | serial | 1 | 26.537369 | 26.745370 | +0.784% |

All process/body time and peak-RSS guards pass. Small positive control deltas are retained, not described as improvements.

## Compiled-module reuse

Two independent AB/BA sessions; each version/session has one warmup and three measured calls. All raw/time/memory guards pass.

| Input | Threads | Session | Baseline s | I4 s | Change |
|---|---:|---|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | 8 | AB | 9.337476 | 8.828089 | -5.455% |
| MjeNMV.MelaMJNV.tlosatx | 8 | BA | 8.697376 | 8.457998 | -2.752% |
| MjeNMV.MelaMJNV.tlosatx | 1 | AB | 13.804216 | 13.634049 | -1.233% |
| MjeNMV.MelaMJNV.tlosatx | 1 | BA | 13.788888 | 13.195333 | -4.305% |
| AP027280.AP027280.tlosatx | 8 | AB | 16.785167 | 15.762611 | -6.092% |
| AP027280.AP027280.tlosatx | 8 | BA | 16.657609 | 15.494707 | -6.981% |
| AP027280.AP027280.tlosatx | 1 | AB | 27.530119 | 25.134989 | -8.700% |
| AP027280.AP027280.tlosatx | 1 | BA | 27.257685 | 25.249807 | -7.366% |

## Same-instance reactor

Two independent AB/BA sessions; each version/session has one warmup and three measured calls. All raw/time/memory guards pass.

| Input | Threads | Session | Baseline s | I4 s | Change |
|---|---:|---|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | 8 | AB | 8.688774 | 8.548302 | -1.617% |
| MjeNMV.MelaMJNV.tlosatx | 8 | BA | 8.439761 | 8.575962 | +1.614% |
| MjeNMV.MelaMJNV.tlosatx | 1 | AB | 13.852748 | 14.038701 | +1.342% |
| MjeNMV.MelaMJNV.tlosatx | 1 | BA | 14.097917 | 13.726973 | -2.631% |
| AP027280.AP027280.tlosatx | 8 | AB | 17.027592 | 15.620948 | -8.261% |
| AP027280.AP027280.tlosatx | 8 | BA | 16.808534 | 15.677629 | -6.728% |
| AP027280.AP027280.tlosatx | 1 | AB | 28.412301 | 26.382654 | -7.144% |
| AP027280.AP027280.tlosatx | 1 | BA | 26.773823 | 26.714014 | -0.223% |

## Actual browser consumer, warmed repeats

Two independent AB/BA sessions; each version/session has one warmup and three measured calls. All raw/time/memory guards pass.

| Input | Threads | Session | Baseline s | I4 s | Change |
|---|---:|---|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | 8 | AB | 8.352405 | 8.177525 | -2.094% |
| MjeNMV.MelaMJNV.tlosatx | 8 | BA | 8.310595 | 8.291675 | -0.228% |
| MjeNMV.MelaMJNV.tlosatx | 1 | AB | 13.761270 | 14.034240 | +1.984% |
| MjeNMV.MelaMJNV.tlosatx | 1 | BA | 13.739005 | 14.097250 | +2.608% |
| AP027280.AP027280.tlosatx | 8 | AB | 16.773155 | 15.912155 | -5.133% |
| AP027280.AP027280.tlosatx | 8 | BA | 16.755940 | 15.938045 | -4.881% |
| AP027280.AP027280.tlosatx | 1 | AB | 28.386135 | 25.125195 | -11.488% |
| AP027280.AP027280.tlosatx | 1 | BA | 28.027685 | 25.440380 | -9.231% |

Node reuse includes full search/output and worker completion; reactor includes input/result ownership operations. RSS and linear-memory budgets are checked separately. Browser uses Chromium 149.0.7827.55 and the frozen real gbdraw consumer, measuring the runLosatPairsParallel promise including transfer/dispatch/search/result completion. Each version/session starts a fresh browser, then retains normal same-page/module/worker ownership. This is warmed component timing, not cold-browser startup or whole-SPA timing. Browser Mje n1 is approximately 2.0–2.6% slower; it passes the original non-regression guard and is not claimed as an improvement.

## Correctness and integration

- 79 command conditions match raw NCBI output: 60 Native/threaded, 15 serial, 4 n2/n4. Source fixtures/options and output hashes are indexed in `I4-functional-index.json`.
- Long AP027131/AP027133 code4: all four Native n1 / serial n1 / threaded n1,n8 results match the saved **NCBI database oracle**. Its argv uses **-db, -query_gencode 4, -db_gencode 4**; makeblastdb uses AP027133. LOSAT uses local -subject with the approved code4 semantics. The archive/inputs/oracle executable hashes were verified; the oracle was reused from round 01, not rerun here. See `recovered-run01/long-code4/manifest.json` and `I4-long-code4/`.
- 802 extracted scan cases × Native/serial/threaded × trace off/on agree. Complete group transcripts agree for 352 cases, 1,927 rounds and 11,494 visits per target/trace mode, including f64 bits, links and HSP state.
- Standard Rust tests: 627 passed, 3 ignored, 0 failed. Candidate clippy/format pass. Four ABI identities and serial API error/recovery pass. Source-identical root integration also passes the two focused tests, clippy and formatting in `I4-integrated-checks.json`.
- Browser raw primary/control checks and short-fixture threaded/nonisolated-serial lifecycle checks pass repeated results, invalid-input recovery, cancellation recovery and pagehide. External requests were blocked; none attempted. Observed 8/1 worker counts are JS construction/termination observations, not proof of physical teardown or hot-loop cancellation. This is component QA, not whole-application or release-bundle certification.
- `I4-integration.json` binds the single production edit; all 158 production build inputs match the tested candidate. Production linking SHA256: `7d23846ca169d5dadf355e0f2de8eefbc58d7f348aef64b258ae2eb5f07cf9ce`.

## Evidence and limits

Every fixed study preserves warmups, all samples, argv, source/artifact hashes, raw outputs and failures. Time regression guard remains max(5%,50ms); RSS/linear-memory guard remains max(10%,16MiB). All seven study policies are COMPLETE and their observed-background guards pass. Normal editor/OS background was observed each second; this is not an exclusive-machine study. New/exited PID CPU can be undercounted, workload-name detection is heuristic, and I/O/cache/frequency effects or subsecond processes cannot be excluded. The background bounds are observed values, not strict external-CPU upper bounds.

I1 reuse failure, I2b contaminated studies and I3 serial failure remain preserved and are not pooled. `I4-serial-git-diagnostic` stays ineligible. The first I4 browser lifecycle invocation had an unregistered fixture name and failed in the harness before LOSAT; its separate record remains `HARNESS_INPUT_ERROR`.

No unlimited reactor-memory stability, existing statistically invalid-query exit-status resolution, Wasmtime performance, full browser UI or deployed-bundle claim is made. Existing deferred issues remain separate. No publishing/deployment occurred.

## Review separation

- Production: one Rust file, with the two new unit tests in its test module.
- Host/runtime: no new host changes in I4; prior accepted H1 baseline preserved.
- Evidence tooling: isolated differential/build/measurement/browser scripts and one-second environment observer.
- Documentation: this result record, candidate policies, plan/status updates.
- Generated evidence: frozen sources, built artifacts, raw outputs, emitted code and logs under run-03.

Independent review scope and remaining audit limits: [I4-INDEPENDENT-REVIEW.md](I4-INDEPENDENT-REVIEW.md).
