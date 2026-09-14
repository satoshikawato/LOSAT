# Remediation report — user-directed provisional adoption

The current candidate fixes the megablast initial X-drop discrepancy, keeps
BLASTN HSP payloads stationary during endpoint pruning, reduces BLASTP
traceback-row bookkeeping, and replaces TBLASTX final payload sorts with
reference sorts. Standard Wasm distribution and verification use threaded
artifacts. The isolated native loop-alignment configuration passes all 14
fresh native performance conditions. The user explicitly requested stopping measurements and provisional adoption.
The qualified native configuration is now in root, whose 589 guarded files
match the C-NA64 snapshot exactly. This does not turn the TBLASTX reuse gate
into PASS: 9/12 time conditions fail, while all 144 job outputs/lifecycles and
4/4 RSS/budget comparisons pass. Six subsequent groups and the final native
megablast recheck remain unexecuted. See the user decision and root-adoption
records. No further measurement was run after the stop instruction.

## New reuse timing failure

The completed C-NA64 X reuse run retains all 17 parent records and 144 jobs.
Session-0 command median increases are 15.00% (AP), 28.80% (Mela) and 20.85%
(Mje); reactor increases are 32.44%, 27.87% and 54.42%. All three session-1
command comparisons pass; all three reactor comparisons fail. The same B1
artifact also slows between sessions, so environment and reused-state effects
need investigation. One candidate AP reactor job takes 88.094945 s, including
88.094170 s inside runPair plus I/O/API copies, 0.000056 s instantiation and
0.000708 s post-return worker wait, with 144.576477 s Node-process CPU. A wait
only or adjustable-clock explanation is insufficient.

Read-only bulk leaf/hash checks overlapped portions of the run; their exact
monotonic boundaries were not separately recorded. This is disclosed as a
possible external influence, not a proved cause or a reason to discard samples.
No record is excluded, corrected or reclassified as a timing PASS. The fixed
AP-only 12-job stage diagnostic completes with both versions about25–27s per
call; it differs from the original three-case interleaving and cannot replace
the failed cohort. A subsequent fixed36-job interleaved diagnostic was stopped
by user instruction (exit130) with partial evidence retained. The predeclared
additional-five confirmation was never run. The original 37 cold conditions below remain valid
observations; they do not waive the separate reuse-time gate.

## Measured performance

Positive reduction means less elapsed time. Each cold condition has one warmup
and five measured samples per version, with alternating order, exact actual raw
output checks, and zero excluded samples. Baseline is parity-correct B1. These
are observed medians on the recorded machine/runtime, not universal guarantees.

| Threaded Wasm n8 input | Baseline (s) | Candidate (s) | Reduction |
|---|---:|---:|---:|
| BLASTN LC738874/LC738870 | 0.802802 | 0.803456 | −0.08% (+0.65 ms) |
| BLASTN AP027202/LC738875 | 11.574859 | 6.207426 | 46.37% |
| BLASTP AP027078/AP027131 | 7.449706 | 6.618647 | 11.16% |
| BLASTP AP027132/NZ_CP006932 | 10.786954 | 9.556909 | 11.40% |
| TBLASTX MelaMJNV/PemoMJNVA | 4.120223 | 4.130872 | −0.26% |
| TBLASTX MjeNMV/MelaMJNV | 13.576849 | 13.923771 | −2.56% |
| TBLASTX AP027280 self | 25.556511 | 25.562424 | −0.02% |

The user set a 5% main-input floor and explicitly waived that floor for short
BLASTN. Its earlier standalone 3.84% reduction is not substituted for the final
near-zero result. Large BLASTN and both BLASTP inputs meet 5%. The time
nonregression allowance remains `max(5% of baseline median, 50 ms)` for every
condition; RSS allowance remains `max(10% of baseline maximum, 16 MiB)`.

| Mode / input | Native n1 reduction | Native n8 reduction | Wasm n1 reduction |
|---|---:|---:|---:|
| BLASTN short LC | −0.33% | 6.13% | 0.03% |
| BLASTN large AP/LC | 19.85% | 33.65% | 27.90% |
| BLASTP AP027078/AP027131 | −3.80% | −3.28% | 14.95% |
| BLASTP AP027132/NZ_CP006932 | −3.35% | −0.82% | 14.64% |
| TBLASTX MelaMJNV/PemoMJNVA | 1.58% | 2.55% | −0.88% |
| TBLASTX MjeNMV/MelaMJNV | 0.03% | 1.17% | −0.43% |
| TBLASTX AP027280 self | 1.12% | 0.48% | 2.04% |

All these time/RSS comparisons pass their predeclared allowances. Native
BLASTP is slower; its Wasm gains do not establish native gains. Native short-LC
n8 saves 21.84 ms with overlapping ranges. All native TBLASTX and Wasm n1
TBLASTX ranges overlap. The AP self Wasm n1 median falls from 38.063448 to
37.286720 s. Some TBLASTX conditions are faster; there is no general TBLASTX
speedup claim. Additional task controls include three megablast inputs whose
Wasm n8 changes range from −0.89% to +0.21%, all within nonregression limits.

The native values above use freshly qualified C-NA64. Existing integrated Wasm
observations apply because rebuilt C-NA64 command/reactor binaries and all four
packaged JS files are byte-identical. Baseline and candidate use the same actual
runner paths and Node 24.21.0. Only TBLASTX receives the identical pair of Node
flags `--no-liftoff --no-wasm-tier-up`; N/P receive neither. Cold process,
compiled-module reuse and same-instance reactor reuse are separate boundaries.

## Implementation and correctness

- Megablast uses NCBI's zero-initialized initial ungapped X-drop option; blastn
  keeps its separate default. The three current megablast cases and all 14
  existing BLASTN task cases match official NCBI raw output. Frozen Sakai
  release expectations are preserved, including their separate contract failure.
- BLASTN endpoint pruning stores `Option<BlastnHsp>` payloads in one arena and
  sorts/copies pointer-sized handles. Deletions still drop payloads at the same
  pruning point; survivors transfer once. Stable comparators, active prefixes,
  scripts and output order remain intact. The focused endpoint filter runs
  eight BLASTN and four existing TBLASTX tests.
- BLASTP fills initialized traceback cells in reserved rows and sets row length
  once. Fence/X-drop boundaries, forward/reverse paths, scores, scripts and
  fresh/reused scratch are covered. Main profiles have zero preliminary retry
  counts; they are not evidence for restricted-to-exact retry coverage.
- TBLASTX retains the initial owned sort, frame groups, indexed linking kernel,
  singleton handling and final replay clones. Reference sorting removes 4N
  final-sort payload clones for N>1. Twelve focused tests cover payloads,
  floating-point bits, ties, addresses and replay indices. Temporary reference
  vectors are search-local and do not enter TLS, globals or retained results.
  Original payloads, references and the replay output can coexist; fewer clones
  do not by themselves prove lower peak memory.
- Native x86_64 configuration aligns loops to 64 bytes with maximum padding 63.
  All Rust source bytes are unchanged from the integrated candidate. Ordinary
  Cargo build produces the exact previously diagnosed native artifact. The
  normalized dominant linking-function comparison preserves all 5,572 non-NOP
  instructions and control flow; this is not a whole-binary equivalence proof
  or proof of one unique microarchitectural cause. Fresh paired native output,
  time and RSS checks determine acceptance.
- Standard build/comparison/CI routes use threaded command/reactor artifacts.
  Explicit serial compatibility remains because read-only inspection found
  current gbdraw AUTO/fallback consumers. gbdraw was neither changed nor
  certified. Default packages reject stale serial artifacts.

TBLASTX local-subject non-default `db_gencode` retains the approved project
behavior. The protected long code-4 case uses the corresponding NCBI database
oracle, not the differing NCBI local-subject semantics. Other behavior has no
new parity exception. The LC threshold cases retain exact raw output at E-values
10, 100 and 10000. These threshold and long code-4 successes apply to their
previously recorded artifacts; the final C-NA64 rechecks were not executed.
NCBI remains an external validation oracle only.

Integrated Rust verification passes 621 tests, with three pre-existing ignored
cases, plus format and Clippy. Default/compatibility contract checkers record
256/373 parent runs, including oracle and expected-failure invocations. Each
main threaded reactor checker has 144 API observations: 108 successes and 36
expected failures, including recovery. These counts are not all successful
searches. Final n2/n4 diagnostics check exact raw outputs and actual worker
lifecycles; they are untimed, not performance medians.

## Repeated memory and the approved scope decision

The original N/P reuse cohort has 192 jobs. It passes all raw outputs/lifecycles,
15/16 invocation-time comparisons and 4/4 whole-process RSS comparisons. The
short-LC reactor session-0 increase is 52.390055 ms against a 50 ms allowance.
One predeclared additional-five cohort passes all eight additional and eight
pooled comparisons and both RSS checks. Ten pooled samples come from two
separate five-sample instance cohorts, not ten consecutive calls. The original
failed record remains unchanged.

The separate fixed 16-call observation has 256 untimed jobs. Its final-eight
linear-memory condition fails in 6/16 series across baseline and candidate;
`FAIL_FIXED_OBSERVATION` is retained. No observation window was extended until
it passed. The user explicitly authorized treating this shared runtime issue
separately while preserving memory nonregression for the performance changes.
The decision is recorded in `runtime-memory-scope-decision.json`; the unresolved
runtime work is in `docs/wasm_reactor_memory_followup_20260914.md`.

A separate finite allocator diagnostic contains 96 untimed jobs with validated
malloc-family ABI coverage, zero accounting errors, exact raw outputs and all
worker exits. Post-call retained usable-byte maxima/final values and block
maxima are equal in each of eight baseline/candidate case-session pairs:

| Input | Retained usable bytes | Live blocks | Retained result bytes |
|---|---:|---:|---:|
| short LC BLASTN | 1,419,328 | 40 | 200,989 |
| large BLASTN | 6,334,528 | 40 | 4,238,551 |
| BLASTP P078 | 1,616,940 | 41 | 385,994 |
| BLASTP P132 | 2,141,132 | 41 | 530,562 |

These values include the current result and exclude allocator headers/free
chunks. They are neither RSS nor linear memory. Later linear growth appears
in 15/16 diagnostic series without a larger observed post-call live maximum.
One uninstrumented fixed-window P132 session-0 candidate linear maximum is
655,360 bytes larger; it is disclosed, not treated as zero growth. Actual RSS
and configured budgets remain within the agreed checks. Instrumentation can
change scheduling/layout; these finite observations do not prove an unlimited
bound. TBLASTX allocator live bytes remain unmeasured/N/A; its new owners are
search-local. Final capacity checks were not run; X reuse time failures remain as described above.

Source inspection identifies delayed crossbeam epoch reclamation as a possible
contributor, not a numerical attribution of the observed allocations. No epoch
flush, permanent pool or hidden production cleanup workaround was added.

## Preserved failures and release boundary

The original integrated native TBLASTX candidate fails four of six timing
conditions. C-X3 loop extraction and C-X4 dead-field removal were rejected.
Their source, failed outputs and measurements are preserved; none entered root.
The new C-NA64 native measurements qualify a distinct artifact and do not
rewrite earlier failures. Its copied source-manifest cwd field is corrected in
a separate metadata record. Current mutable Cargo fingerprints are not used to
infer historical B1 library build flags.

All measurements retain realtime/GNU clock diagnostics without adjusting or
excluding samples. Elapsed acceptance uses CLOCK_MONOTONIC with outer BOOTTIME
agreement. Individual reactor jobs do not have a BOOTTIME observation. Logical
copy bytes and aggregate parallel stage times are not memory-bus measurements
or wall-time fractions.

Frozen PR5 bytes, registered native platform fingerprints and the 46 previously
deferred cases remain unchanged. Corrected Sakai output differs from the frozen
LOSAT contract; formal release recertification requires a separately reviewed
authority version. This work does not claim full release certification.

The fixed C-NA64 ledger remains stopped at its recorded failure. The user
explicitly adopted the current source provisionally and stopped additional
measurements. Separate production/test/tooling/docs diff review and the local
source/artifact archive preserve the handoff. The user subsequently requested a
commit and push including every tested-build prerequisite. The source-to-index
binding is recorded in `commit-build-input-binding.json`. This Git handoff does
not constitute release certification or a gbdraw migration.
