# Independent C-NA64 TBLASTX reuse review — time FAIL

Read-only ncbi_parity_auditor review confirms 17 parent records, eight processes and 144 jobs (120 timed). All 250,298,880 actual job-output bytes match the official oracles. Exact argv/artifacts/Node/source, 1,152 worker lifecycles and 4,608 events match. Zero samples are excluded.

| Input | Mode | Session 0 time increase | Session 1 time increase |
|---|---|---:|---:|
| Mela/Pemo | command | 28.80% FAIL | 0.08% PASS |
| Mje/Mela | command | 20.85% FAIL | 2.48% PASS |
| AP self | command | 15.00% FAIL | 0.04% PASS |
| Mela/Pemo | reactor | 27.87% FAIL | 6.14% FAIL |
| Mje/Mela | reactor | 54.42% FAIL | 10.32% FAIL |
| AP self | reactor | 32.44% FAIL | 15.63% FAIL |

All 12 actual-array time verdicts agree with the saved evaluator. Session-0 B/C ranges are disjoint for all six conditions. All four actual-process RSS comparisons pass max(10%,16MiB), within the 4GiB budget. Time allowance remains max(5%,50ms); the memory scope decision does not waive it.

The same B1 slows between sessions: command median wall +18.84–20.28% and CPU +14.82–17.78%; reactor wall +13.11–19.70% and CPU +12.10–14.70%. Constant baseline speed is not established. This does not cancel the failed gate.

Candidate session-0 reactor AP repeat1 actually takes 88.0949448 s: instantiation 0.000056 s, runPair plus I/O/API copy 88.0941703 s, post-return worker wait 0.000708 s; Node-process CPU is 144.576477 s. Neither ending waits nor adjustable clocks explain it. The 88 s observation is not the selected median, so removing it would not resolve the AP failure.

All commands use NCBI BLAST+ 2.17.0+, TBLASTX/outfmt6/gencode1/n8, without an exception. Actual command/reactor bytes, shared runner paths and Node flags match; 589 files × three source trees and 87 launch bindings match. Native alignment configuration cannot explain a Wasm change because the actual Wasm bytes are identical. The eight AB/BA processes do not overlap internally.

Read/hash sessions53632 and98235 were reported to overlap portions of the run. No dedicated monotonic start/end record exists, so exact overlapping jobs and effect sizes cannot be established. This is a possible external influence, not a proved cause or basis for exclusion/correction. Cold and reuse boundaries differ; subtracting cold startup costs does not isolate a cause.

Actual reactor final-three linear-memory stability is11/12. Candidate session1 Mje grows from96,600,064 to102,039,552 bytes at the final call. Command plateau flags are N/A sentinels. Preserve this observation under the approved common-runtime follow-up, separately from the time failure.

The outer ledger stopped with exit1 at x-reuse. Runtime scope disposition was not generated, and the six subsequent groups have not executed. Bulk audit stopped when the separate12-job stage diagnostic began. That diagnostic is not adoption evidence or a replacement for this FAIL.
