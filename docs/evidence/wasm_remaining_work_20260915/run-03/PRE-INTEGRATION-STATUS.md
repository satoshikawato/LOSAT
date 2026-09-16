# Round 03: normal-runtime TBLASTX n8 reversal

## Current status (2026-09-16)

Production Rust is unchanged. I1 failed its fixed reuse guard despite large
normal-runtime cold n8 improvements. I4 is the current isolated candidate;
see `I4-POLICY.md`, `I3-POLICY.md` and the preserved I2b evidence in `DIAGNOSIS.md`,
`I2b-POLICY.md`, `I2b-CODE.md` and `I2b-independent-static-review.md`.

- I2b: 802 scan cases × 3 targets × trace off/on pass exactly.
- I2b: six fresh-oracle fixture comparisons at n1/n8 pass raw bytes and thread
  contracts (12 runs). Actual generated code removes the identified extra
  address calculation; no all-bounds-eliminated claim is made.
- I1 and I2b: 352 complete group cases × 3 targets × trace off/on pass
  exactly, including helper/link/HSP state and floating-point bits.
- I2b first reuse study: AB 9.472270 → 9.294128 s; BA 9.518283 →
  10.033301 s, which numerically fails the 5% guard. A separate gbdraw
  benchmark started at boot time 1123.08 s, overlapping both BA processes
  (candidate 1119.38–1161.53 s, baseline 1162.09–1214.74 s). All original
  data and the numerical failure remain saved under `I2b-module-mje-n8/`.
  The overlapping session cannot establish adoption or a candidate regression.
  A new, separate exclusive fixed set subsequently passed: AB 9.408776 →
  9.025415 s, BA 9.197112 → 9.072569 s; all time/RSS/linear-memory guards.
  Its 16 raw outputs and seven-child lifecycles were independently audited.
- I2b full-group state comparisons now pass on all three targets, trace off/on.
  Raw/thread verification totals 79 successful command conditions (12 initial,
  45 control matrix, 18 primary Native/serial, 4 n2/n4).
- Long AP027131/AP027133 code4 passes Native n1, serial n1 and threaded n1/n8.
  The round-01 DB oracle archive, query/subject and oracle executable hashes
  were verified before reuse; this is not a newly run oracle.
- Chromium 149.0.7827.55: primary n8 raw outputs and the BLASTP n8 control pass.
  Actual consumer threaded/nonisolated-serial lifecycle checks pass repeated
  results, invalid-input recovery, in-flight cancellation recovery and pagehide.
  Lifecycle uses one short TBLASTX fixture. Cancellation occurs after sending
  run; it does not prove interruption inside the hot loop. Worker maxima 8/1
  count JS construction/termination calls, not physical teardown or CPU activity.
  External requests were blocked;
  none were attempted. These are LOSAT consumer-component checks, not whole-app UI
  or release-bundle certification.
- Focused tests pass; release library suite 433 passed/1 ignored; clippy and
  cargo fmt checks pass. Standard cargo test passes 627 tests, with 3 ignored
  and no failures. All four command/reactor ABI shapes match baseline; serial
  reactor unsupported-thread errors and subsequent recovery also match.
  Production Rust is still unchanged pending performance acceptance.
- Independent functional audit closed at 2026-09-16 02:44:26 UTC: all 360
  candidate source files matched their manifest; the 79 raw outputs and thread
  contracts, saved code4 oracle binding, state hashes, ABI and test logs were
  rechecked. The four long-code4 runs prove raw parity only; their child-worker
  lifecycle was not separately instrumented. No new parity defect was found.
  Performance acceptance remains pending.
- I2b serial Mje screen numerically fails (+9.53% body), but the last
  candidate sample overlaps a newly started external benchmark for 10.254939 s.
  It is invalid for acceptance; all samples and the numerical failure remain.
  I3 limits `inline(never)` to threaded WASI and requests serial inlining.
  Its 802-case extracted scan state/trace differential passes on all three
  targets. The harness does not include the inlining attributes; actual serial
  emission and 15 raw-oracle conditions separately validate that boundary.
  All five artifacts are built; four ABI shapes and serial API recovery pass.
  Chromium primary n8 inputs and the BLASTP n8 control pass on I3. The threaded
  group/helper WAT is byte-identical to I2b, but whole artifacts differ.
  Independent I3 source/build/functional review closed 03:15:50 UTC. Clippy
  and format pass; standard Rust suite is 627 passed/3 ignored/0 failed.
  Long code4 raw output matches the saved DB oracle on all four conditions.
  I2b timings are not transferred to I3.
- I3 is rejected: clean fixed Mje serial body median 13.540885 → 16.847026 s
  (+24.416%), all raw equal, RSS guard passes. All samples are preserved.
  I4 restores the indexed caller at source level and shares the selection body
  through one macro. Selection/trace tokens are preserved after explicit macro
  argument normalization. All three target scan differentials and the six-case
  serial raw gate pass. Actual timed serial group WAT has the baseline
  instruction structure, with 201 constant values and 7 memory offsets still
  different. Independent source/emission audit closed 03:43:29 UTC.
  I4 diagnostic-only serial screening with recorded external editor Git indexing
  passes the numerical guards: body 13.848079 → 14.115543 s (+1.931%),
  process +2.090%, RSS −0.259%; eight raw outputs match. This background-load
  study is not eligible for acceptance. The separate clean first screen is pending.
- I4 full functional checks now include 352 complete group cases × 3 targets ×
  trace off/on (1,927 rounds and 11,494 visits per run), all exactly matching
  each target baseline. Standard Rust tests: 627 passed, 3 ignored, no failures;
  clippy and format pass. All four ABI shapes and serial error recovery pass.
  Browser primary/raw and threaded/nonisolated-serial lifecycle checks pass;
  the first lifecycle invocation used an unregistered fixture name and failed
  in harness lookup before LOSAT, then was corrected in a separately saved run.
  Long code4 raw output matches all four conditions. Actual threaded group/helper
  WAT matches I2b byte for byte; no whole-artifact or speed transfer is claimed.
- Remaining: current-candidate cold n1/n8, Native/serial performance controls,
  remaining command/reactor reuse, browser performance, final source integration
  and audit. All timing waits for other active
  repository builds/tests/browser jobs to finish. The last I4 command controls
  and independent read-only audit are still running.

## Prospective decision policy

The user identified the unresolved n8/n1 reversal as the priority. I1 extracts
the unchanged large-gap predecessor loop into an actual Wasm function call per
HSP, while requesting Native inlining. This revisits C-X3 with the currently
adopted b0 fast path, preserved caller initial values, and Native inlining;
it is not a new algorithm or an accepted optimization.

First complete the boundary unit tests and fresh six-fixture raw/thread gate.
Then measure the two primary TBLASTX inputs (MjeNMV/MelaMJNV and AP027280 self)
under normal Node, with n8 followed by n1. Each condition uses one warmup per
version and three fixed AB/BA/AB pairs. Both versions receive the same coarse
Rust dispatch clocks and identical host, environment and affinity. No other
build, test, profiling, compression or browser job runs during timing.

These first timings screen the hypothesis; they do not authorize adoption.
All samples, raw outputs, clocks, RSS, source bindings, executable hashes and
failures are retained under this directory as they are produced. A failed
condition is not extended until it passes. The time regression limit is
max(5%, 50 ms), and the RSS limit is max(10%, 16 MiB). The 1.20 Wasm/Native
ratio is a goal, not an adoption requirement. Claim reversal resolved only
where normal-runtime candidate n8 is faster than candidate n1 for the same
input. Report each input separately.

Only a promising candidate proceeds to complete internal state/FP/ordering
comparisons, Native and serial controls, compiled-module and reactor reuse
(two AB/BA sessions, each warmup plus three samples), the relevant full parity
sweep, and actual browser validation. Inspect the emitted Wasm function/call
and separate diagnostic execution evidence before claiming a tier mechanism.
Compiler flags used in diagnostics never enter adoption measurements.

Round 02 source and raw measurements were lost at environment reset; its
progress notes cannot serve as adoption evidence. See ../run-02/RECOVERY.md.
Round 01 archives remain the preserved source of earlier diagnostic evidence.

## Source ownership

- `work/baseline-source`: frozen current baseline, including accepted H1 hosts.
- `work/I1-source`: isolated candidate; the actual production tree is unchanged.
- `work/*body-source`: matched measurement-only dispatch clocks.
- `work/*/artifacts`: actual built files, copied and hashed immediately.
- `/tmp/losat-run03-*-cache`: disposable build caches only.

NCBI reference: `c++/src/algo/blast/core/link_hsps.c:827-863`, original
predecessor visit/jump/selection loop. Search candidates, ties, coordinates,
floating-point evaluation, scheduling and output order remain unchanged.
