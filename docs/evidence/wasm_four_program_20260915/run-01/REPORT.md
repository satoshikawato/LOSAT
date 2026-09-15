# Four-program Wasm optimization — run-01

Status: **ADOPTED by explicit user instruction on 2026-09-15; measurement ended after three complete pairs.**

[Adoption decision](ADOPTION.md): retain the integrated implementation with the
observed performance tradeoffs. The six time-guard failures and unrun scope
remain recorded below.

## Authority and frozen input

This implements `docs/wasm_four_program_optimization_plan_20260915.md`. The dirty
working tree at HEAD `1de348ad677af0290ddb08c8d95f7ea52975f5d4` is the baseline;
HEAD alone does not identify it. `baseline-inputs.tar.gz` and
`baseline-inputs.json` freeze the required source, configuration, fixtures and
runners. `initial-status.txt` / `initial-user.diff` identify pre-existing changes.
Those changes are preserved.

NCBI C/C++ source authority is the local ncbi-blast checkout at
`598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison executables are NCBI
BLAST+ 2.17.0+. NCBI is used only as a validation oracle. This run does not update
PR5 frozen output, PR6 platform fingerprints, or release certification.

The execution workspace is `/tmp/losat-four-program-20260915`.
`execution-evidence.tar.gz` and its SHA256 manifest retain source snapshots,
artifacts, raw outputs, commands, failures and interrupted records. Node is **26.8.2**, so older
Node 24 measurements are not this run's baseline.

## Candidate ledger

| ID | Change | NCBI owner | Current status |
|---|---|---|---|
| M0 correction | Preserve already-written nonaffine scratch on fence return | `core/greedy_align.c:523–526,571–576` | Failing regression reproduced and fixed; separate from optimization |
| M1 | Nonaffine rows refer directly to base/pool storage | `core/greedy_align.c:71–77,100–125,666–678` | ADOPTED in integrated version by user; raw/state gates pass |
| N1 | Carry prepared coordinates with speculative traceback | `core/blast_traceback.c:436–475`; `core/blast_gapalign.c:3323,4163` | ADOPTED in integrated version by user; raw/unit gates pass |
| N2 | Reuse batch jobs, results and per-slot vectors; batch16 unchanged | `core/blast_traceback.c:403–405,509–513,583–612` | ADOPTED in integrated version by user; raw gates pass |
| N3 | Precompute the exact BLASTNA four-base expression per position | `core/na_ungapped.c:262–349` | ADOPTED in integrated version by user; raw/exhaustive gates pass |
| P1 | Lazily construct one immutable initial matrix per search; clone completed values | `core/blast_kappa.c:2228–2236`; `composition_adjustment.c:Blast_Int4MatrixFromFreq` | ADOPTED in integrated version by user; raw/matrix bits pass |
| P2 | Share single-query preparation and reuse independently owned DP/composition scratch | `core/blast_kappa.c:2309–2336,3329–3334,3493–3504` | ADOPTED in integrated version by user; raw gates pass |
| X1 | Stop unchanged maximum-tree propagation | `core/link_hsps.c:605–650` | SKIPPED: tree updates about 0.1% of native summed group time |
| X2 | Defer large-gap coordinate reads until sum qualifies | `core/link_hsps.c:827–861` | ADOPTED in integrated version by user; source/state/raw gates and all ten TBLASTX cold time/RSS guards pass |
| P3 | Diagonal offset reuse | `core/blast_extend.c:162–184` | SKIPPED: no proof for scheduler-dependent subject order/rollover |
| P4 | Bounded speculative redo | `core/blast_kappa.c:3525–3539` | SKIPPED: all 32 prepared matches are needed in the declared single-query fixture |
| M2 / P5 | Further mismatch/DP kernel changes | Respective greedy/ALIGN_EX owners | Conditional investigation; no new kernel change justified yet |

The P2 range cache remains per match: retaining it across subjects would retain
subject copies. Mutable composition and DP state remain private to each Rayon
work unit. Initializer count is measured, not equated with thread count.

## Completed correctness evidence

- M1: 4,800 fresh/reused scratch cases match the corrected M0 baseline, including
  retained cells, scores, coordinates, edit operations and seed state. Forward /
  reverse, score-only / traceback, ambiguity, fence and increasing retry budgets
  are covered. See `scratch-comparison-fresh.json`.
- N1: the ordered speculative regression crosses the 90,000-base subject shift
  boundary and has at least 20 unequal gapped matches (batch16 boundary).
- N3: 65,536 BLASTNA four-letter combinations and all four phases; 785,200 actual
  approximate-extension comparisons include short/context boundaries and exact
  fallback thresholds. See `N3-extension-test.log`.
- P1: 530 query matrix/frequency/lambda records have equal bit hashes; initial
  matrix construction falls from 530 to one. The diagnostic `callback_new=true`
  is a literal, **not a measurement of callback independence**; fresh callback
  ownership is source-reviewed. See `P1-bits-comparison.json`.
- X2: 128 actual linking groups yield 13,014 equal transcript lines, including
  helper visits, predecessor state, `xsum.to_bits()`, final chains/order and
  E-value bits. Trace-target execution also has equal state and 214 equal full diagnostic
  lines (Cargo build headers excluded). See `X2-state-comparison.json` and `X2-fulltrace-comparison.json`.
- Final integrated Rust tests: 625 passed / 3 ignored; clippy all-targets /
  all-features with `-D warnings` and `cargo fmt --check` pass. The final
  integrated main/control matrix has 130 passing records. Each version also
  has 317 runtime records, 150 supported edge records, 12 expected unsupported
  records, nine long code-4 comparisons and 25 browser lifecycle records.
  See `validation-review.json`; this explicitly preserves the failed probes below.

An early scratch comparison accidentally reused a Cargo test artifact across
copied crates. That transcript is **invalid and excluded**. Fresh separate
compilations and different binary hashes establish the replacement evidence;
see `verification-correction.md`. No conclusion relies on the invalid run.

## Diagnostics and generated Wasm

On EDL933/Sakai, baseline nonaffine greedy performs 23,921 calls and 813,486 row
allocations; estimated copied bytes are an upper bound of 1,235,110,504. NZ self
has 4,178 calls / 136,531 rows / 183,812,488 estimated bytes. These are diagnostic
counts, not elapsed-time improvement claims.

For large BLASTN, 59,119 speculative and 59,087 ordered preparations occur; the
ordered preparations all have already-prepared results on this input. The
fallback is retained, but this fixture does not demonstrate its activation.

For MjeNMV/MelaMJNV, both native and threaded X0 diagnostics visit 4,047,645,435
large-gap helpers, of which 1,534,074,932 (37.9004%) fail the sum predicate.
INDEX1 takes 84.83% / 82.62% of **summed group time**. This is not a wall-time
breakdown. The existing tree is inexpensive in the native diagnostic, so X1
was skipped.

`baseline-linking.s` loads coordinates before the sum branch. `X2-linking.s`
loads only sum / next_larger in the rejection loop and defers coordinate loads
to contenders. Total static load count does not decrease; the improvement
hypothesis concerns loads executed on rejected visits. Independent read-only
source/assembly review found no change to jump order, comparisons, predecessor
selection, or floating-point evaluation order.

## Long code-4 boundary

The fresh AP027131/AP027133 full-length code-4 **NCBI database oracle** exceeded
300 seconds. That is recorded as TIMEOUT, not a LOSAT discrepancy or a PASS.
Dependent comparisons were stopped. The single fixed 1,800-second correctness-only
retry on CPU31 completed successfully (617.97 monotonic seconds; the adjustable
GNU/realtime clock reported 671.84 seconds). Both versions match the oracle at native n1/n2/n4/n8, explicit serial
Wasm n1 and threaded Wasm n1/n2/n4/n8. The approved local-subject genetic-code exception
is applied narrowly; no other difference is waived.

## Browser baseline

The current gbdraw application was served locally without modifying its source,
with this run's baseline command Wasm files overlaid at its existing URLs.
Chromium initially failed under the tool sandbox; the identical local check
passed after sandbox escalation. See `browser-baseline.json`.

Normal browser settings, desktop/mobile contexts, nonisolated serial and
isolated threaded execution, two invocations each of blastn/megablast/tblastx/
code4 tblastx/blastp, threaded n8, and cancellation followed by recovery all
match the fresh native outputs. All external requests were blocked and no
external attempt or page error occurred. This proves the consumer boundary;
it is not a separate rebuilt gbdraw distribution certification or a five-pair
browser speed claim. Both final C0 and integrated browser lifecycle/budget checks also pass; see
`browser-C0-lifecycle.json` and `browser-integrated-lifecycle.json`.

## Initial individual measurements

[Initial results](individual-performance-initial.md) cover N1/N2/N3/M1/P1/P2/X2.
All seven individual candidates were measured. These initial results do not override the six final integrated time-guard failures. No individual effects are summed.

Initial command build records bind argv and artifact SHA256 but did not record
source hashes at build time. Later candidate snapshots include test-only edits;
these are not contemporaneous source bindings. Historical rebuilds were not
performed. Final C0/integrated bindings in `final-build-bindings.json` are clean.

## Measurement contract

The original `measurement-plan.json` declared one warmup and five alternating A/B pairs,
monotonic whole-process time and peak RSS, CPU affinity0–7, identical runner
paths, diagnostics disabled. Native n8 and threaded command Wasm n8 are the
individual conditions; default Node flags are primary. Only the code4 oracle
on CPU31 may overlap. Compilers, browser and other validation jobs finish
before adoption samples. One additional five-pair cohort is permitted only
for inconclusive results, and all ten pairs must then be evaluated together.

C0 corrects two baseline defects: M0 scratch persistence and the BLASTN multi-subject cutoff described below. Final integration measurements use C0; original individual measurements retain their original artifact identities.

The user's final instruction ends all measurement and changes the final sample
count to three. `measurement-stop.json` records that amendment; original raw
metadata is preserved, including its last RUNNING checkpoint. No jobs remain
running. The first three complete pairs are used at all 45 conditions; 54
completed repeat-4 records and the interrupted invocation remain retained and
excluded. No condition or sample was selected by speed.

## Final three-pair result and decision

[Full cold results](integrated-three-pair-results.md) and
[included samples](integrated-three-pair-summary.json): 270 timed processes,
all raw-equal; time guard **39/45**, peak-RSS guard **45/45**. The independent
NCBI auditor recomputed the same result. The five-percent target is met in
18/45 conditions. This is an optimization evaluation, not v0.1.0 release
certification.

TBLASTX passes both guards in all ten conditions. Threaded n8 median time
improves 16.32% on MjeNMV/MelaMJNV (42.592→35.640 s) and 14.40% on AP027280 self
(74.282→63.585 s). This supports those measured cold gains, not a Native-equivalent
engine/body or reuse claim.

Six time guards fail: AP027078/AP027131 BLASTP threaded n8 (+7.36%),
Mela/Pemo BLASTN serial n1 (+7.70%) and threaded n8 (+6.21%), Mj/Ml BLASTN
native n8 (+11.00%), NZ self megablast native n1 (+64.95%, +163 ms), and
single-query/32-match BLASTP threaded n8 (+13.73%). Their cause is unresolved;
the small sample count does not justify dismissing them as noise.
The user subsequently instructed 「じゃあいいじゃん。採用」 after reviewing the
results and tradeoffs. The integrated implementation is **ADOPTED**; see
[the decision](ADOPTION.md). All-condition nonregression remains unproven,
with these six time-guard failures retained. No further measurements or
production code changes were made after the stop instruction.

Final integrated control timings, body timings, compiled-module reuse timings,
same-reactor timing/memory windows, and secondary TurboFan timings are **NOT RUN,
canceled by user**. Existing completed control parity and individual timings
remain recorded; they do not substitute for these unrun conditions. Existing
reactor-memory and TBLASTX reuse-time issues remain open. A separate open
statistically-invalid-query TBLASTX exit-status defect is detailed below.

## Additional baseline correction: subject-set cutoffs

The word7/unequal-subject boundary fixture found 18 extra HSPs in the frozen
baseline (two per query/subject pair). This is a baseline discrepancy, not a
performance candidate regression. The ordinary NCBI CLI `-subject` path uses
database mode (`blast_app_util.cpp:204–210`), sums all subject lengths
(`seqsrc_multiseq.cpp:175–180`), and skips `BLAST_OneSubjectUpdateParameters`
when that total is nonzero (`blast_engine.c:1434–1445`). Hit and initial-word
cutoffs use the existing query context's effective search space
(`blast_parameters.c:925–946,370–374`).

LOSAT recomputed both cutoffs for each subject individually. On this fixture it
admitted a score17 ungapped seed at cutoff17, then extended it to an extra
score22 gapped HSP. The proper set-level cutoff is19. C0 uses the already-owned
context search space for both cutoffs. Fresh C0 native output exactly matches
the oracle SHA256 `7c1541325718a0328964b02d93419f64742c7b925bddbff659f50ee24fe5452c`;
the frozen-output regression passes at n1/n4. C0's 60 affected matrix records
and the integrated 130-record matrix pass.
An independent read-only NCBI source audit confirms the boundary. This creates
no new exception and adds no BL2SEQ_LEGACY mode.

## Reuse harness correction

The old reuse assertion expected N child workers. An untimed n2 reproduction
returned correct raw bytes and exactly one started/ready/exited child, then
failed `1 !== 2`. The current runtime contract includes the caller, so the
expected child count is N−1. The one-line test correction passes against the
same baseline artifact; the runtime is unchanged. Original/fixed transcripts
are retained under `reuse-contract-old` and `reuse-contract-fixed`.

## Final verification

Root `cargo test`: 625 passed / 3 ignored; `cargo clippy --all-targets
--all-features -- -D warnings` and `cargo fmt --check` pass. Both C0 and integrated
versions pass the existing runtime harness, including raw format routes,
unsupported-option rejection, worker faults, API and reactor lifecycle checks.
Both actual-consumer browser checks pass 25 records in Chromium149.0.7827.55;
no external request or page error occurs. Successful n8 jobs spawn seven child
threads, concurrent jobs stay within their budget, and pagehide terminates all
tracked LOSAT parent workers. Cancellation and invalid-FASTA errors recover.

The supplemental edge driver initially requested unsupported CLI `-strand`.
Its rejection is retained in `*-edge`; reverse-coordinate coverage uses the
existing default-both-strands search in `*-edge-remainder`. The original
all-C TBLASTX query is valid FASTA but has statistically invalid contexts.
NCBI emits a warning and an empty successful result; baseline/C0/integrated
LOSAT returns `no valid query contexts`. This **open baseline exit-status
parity defect** is retained in `*-edge-remainder/tblastx-nohit`. It is not an
approved exception or a candidate regression. The unchanged guard executes
before the X2 linker. NCBI's missing wrapper behavior is in
`api/local_blast.cpp:177–224`; a separate correctness port is required.
Normal valid-query/no-hit behavior is separately tested in the qualified
`*-edge-tail/tblastx-nohit` case. No all-input or release certification claim
is made. The independent NCBI auditor accepts this scoped boundary.

N2 diagnostic batches preserve job order and scatter indices in186/4,214
batches, including9 empty batches and145 batches with fewer than eight jobs
on the large input. Jobs/results capacities stay at16 and each slot at4.
Actual native allocator instrumentation shows total alloc+realloc calls
42,539→40,159 (short) and500,734→446,667 (large). Peak requested bytes differ
by90 and26 bytes. These are untimed diagnostics; requested bytes are not RSS
or allocator usable bytes. Raw outputs remain equal. See
`N2-batch-comparison.json` and `N2-allocator-counts.json`.

## Secondary TurboFan cold condition

The plan also requires a distinct TBLASTX host-flag condition.
`final-turbofan-declaration.json` declares equal
`--no-liftoff --no-wasm-tier-up` flags for baseline and candidate. Individual X2
uses threaded n8 on two main inputs and one control; integrated main inputs
use native/threaded n1/n8 and serial n1, with the control at native/threaded n8.
These secondary cold measurements were queued but never started; the queue was canceled by the user.
Default Node remains the primary adoption runtime. TurboFan body/reuse and
browser timings are outside this secondary condition; no such claim is made.
