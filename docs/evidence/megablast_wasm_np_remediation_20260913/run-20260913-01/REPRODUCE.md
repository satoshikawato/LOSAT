# Reproducing the remediation evidence

The user stopped measurements and provisionally adopted the current source.
Do not automatically resume unfinished or prepared measurement drivers.
Use the exact saved source and artifacts, not the original dirty HEAD. The
final `handoff-archive.json`, `handoff-file-map.json`, archive verification and
`SHA256SUMS` identify the saved bytes. Until those records exist, archive handoff
is still pending. This directory's `STATUS.md` distinguishes executed checks
from prepared drivers.

## Source and artifact identity

`manifest.json`, `initial-status.z`, `initial-tracked.diff.gz` and
`B0-source.tar.gz` retain the starting state. B0/B1 snapshots contain the
recorded source/config/tooling set, not the later complete integrated test tree.
`integrated-source.json` freezes 589 files, including 154 Rust files and the
regular fixture tree. `C-NA64-source.json` freezes the final candidate: only
`LOSAT/.cargo/config.toml` differs from the original integrated source. Its
copied cwd field is corrected by `C-NA64-source-metadata-correction.json`.
When present, `C-NA64-root-adoption.json` verifies every one of those 589 root
files against the qualified snapshot. It does not replace the earlier manifest.

`change-review/inventory.json` and separate patches compare this task's changes
with captured initial bytes. Production, inline tests, verification/tooling,
documentation and the raw oracle fixture are distinguished. Source diff views
may normalize CRLF for readability; actual source and search output hashes do
not. Extensive pre-existing user changes must not be reverted.

The local archive contains `work/` source snapshots, measured LOSAT artifacts,
runners, drivers and logs; shared `repo/LOSAT/tests/fasta/`; and `evidence/`.
Regular fixture trees remain regular. Existing fixture symlinks become portable
relative links. The file map records every original path, size and hash; archive
verification checks every member without extraction. Cargo caches and NCBI
executables/libraries are external prerequisites, not implementation dependencies.

## Rebuild and run sequence

Use Rust 1.92, locked dependencies, the recorded WASI targets, Node 24.21.0 and
NCBI BLAST+ 2.17.0+ comparison tools. Exact versions, hashes, CPU/runtime data and
argv are saved in each group's metadata. Release settings are opt-level 3, LTO,
one codegen unit and panic abort; Wasm has `+simd128`. Native x86_64 loop
alignment belongs to the candidate Cargo configuration. Current mutable Cargo
fingerprints do not reconstruct historical B1 library build flags; see
`C-NA64-profile-binding-note.json`.

From the qualified source root, build normal artifacts into fresh directories:

```bash
python LOSAT/tests/build_wasi_artifacts.py --target-dir /tmp/losat-replay-target --output-dir /tmp/losat-replay-artifacts --node /path/to/node
```

Add `--include-serial` only for explicit compatibility builds, using a separate
output directory. Standard output rejects stale serial artifacts. Native,
command and reactor target directories stay separate as recorded in the ledgers.

The recorded planned sequence is below. Only native qualification, build
validation and the first three remaining groups executed; X reuse stops the
ledger with time failure. Later steps, including the automatic full-gate
summary/adoption drivers, were not executed. Root adoption instead follows the
explicit user provisional-adoption record and complete589-file hash comparison.
The source guards in old drivers intentionally describe pre-adoption root;
they must not be blindly rerun after the configuration adoption.



1. `qualify_na64_native.py`: six TBLASTX and eight N/P cold native conditions,
   using the frozen matrix in `C-NA64-native-matrix-planned.json`.
2. `validate_na64_builds.py`: format, Rust tests, Clippy and actual threaded
   builds. `C-NA64-wasm-equivalence.json` proves command/reactor and four JS
   files are byte-identical to the earlier integrated bundle. Only exact-byte
   Wasm observations carry forward; native is measured freshly.
3. `run_na64_remaining.py`: the nine fixed, sequential groups in
   `C-NA64-remaining-guarded-planned.json`, with bound source/artifact/driver
   identities and exclusive logs. The steps include X n1/n2/n4, reuse,
   thresholds, long code-4, daily routes and capacity diagnostics.
4. `na64_final_megablast.py`: three fresh official oracles and twelve untimed
   native n1/n8 comparisons for the fixed megablast cases, after those nine
   groups. Its own plan and validation record are separate.
5. `summarize_na64.py` with `na64_summary_proof.py`: verify actual leaf files,
   full argv, exact condition/repeat sets, ordering, elapsed/RSS arithmetic,
   source identities and the approved runtime boundary. Independent review
   precedes `adopt_na64_root.py` and final archive creation/verification.

Drivers retain original absolute machine paths. To replay, extract a snapshot,
map paths consistently and choose new output directories. They are not claimed
to relocate automatically. Never rerun an orchestration driver over the saved
evidence. A prepared plan is not proof that its command executed.

## Timing and raw output

Each benchmark group saves metadata, full ordered command records, actual raw
outputs, diagnostics, every warmup/measurement, exclusions and run status.
Cold adoption conditions use one warmup and five paired measurements in
alternating order. Reconstruct argv from the group's saved `work/*-argv.json`
or planned ledger. Use the same actual runner path and filesystem placement
for baseline/candidate. Heavy builds/tests/benchmarks run sequentially.

Only TBLASTX gets both `--no-liftoff` and `--no-wasm-tier-up`, identically for
both versions; BLASTN/BLASTP get neither. CLOCK_MONOTONIC is the elapsed
acceptance clock, with outer BOOTTIME agreement. Adjustable realtime/GNU
measurements remain diagnostics; no samples are adjusted or selectively
removed. Individual reactor jobs have no BOOTTIME observation.

The user authorized a 5% main speed floor and waived it only for short BLASTN.
Time/RSS nonregression remains `max(5% of baseline median, 50 ms)` and
`max(10% of baseline maximum RSS, 16 MiB)`. Raw output comparisons never sort,
normalize or overwrite the official oracle bytes.

## Reuse and memory scope

Command reuse is one compiled module with a fresh instance per job. Reactor
reuse retains one instance per case. Keep two AB/BA sessions, interleaved case
order, one warmup, five timed calls and actual worker lifecycle records.
`evaluate_reuse.py` is the original strict evaluator. Replaying original N/P
records must retain its one time failure and two final-three plateau failures.

The permitted additional-five timing cohort in `reactor-supplemental-plan.json`
passes all additional and pooled comparisons. Ten pooled values come from two
separate five-sample instance cohorts, not ten consecutive calls. The fixed
memory amendment binds four cases × 16 calls × two versions × two sessions,
256 untimed jobs, with a constant final-eight requirement. Six of sixteen
series fail. `FAIL_FIXED_OBSERVATION` remains unchanged; no window is extended
until it passes.

The user explicitly approved a common-runtime follow-up while retaining memory
nonregression: `runtime-memory-scope-decision.json` and
[the open runtime follow-up](../../../wasm_reactor_memory_followup_20260914.md).
`evaluate_reuse_runtime_scope.py` preserves the original time/RSS/raw/budget and
lifecycle checks and records absolute plateau separately. It does not rewrite
old failures. Observed N/P retained-byte/block equality and the disclosed
655,360-byte P132 linear maximum increase are in
`runtime-memory-np-nonregression-evidence.json`.

The allocator diagnostic is separate from production and adoption timing.
Its plan fixes four cases, six calls, two versions and two sessions: 96 jobs.
`allocator_accounting_v2.rs`, r2 build/smoke records, link maps and
`work/allocator-reference/` bind actual malloc-family ABI coverage. The first
selftest failure and its corrected opaque-call selftest remain separately
recorded; the measured wrapper itself was unchanged. Search observations occur
after raw-result copying, freeing the five input buffers and actual worker
exit, with LAST_RESULT still live. No extra clear/selftest occurs in search
instances. These finite usable-byte counts are not RSS, linear memory or an
unlimited memory bound. TBLASTX allocator live bytes remain unmeasured/N/A.

Capacity evidence uses isolated instrumented source snapshots, not production
artifacts. `C-NA64-capacity-build-and-cycles.json`, diagnostic patches/source
identities and `C-NA64-capacity-summary.json` retain actual capacities and four
repeated scratch cycles for each version. Diagnostic durations are ineligible
for adoption. Owner maxima are neither a simultaneous heap peak nor allocation
counts. The long code-4 check reuses unchanged fresh C-X2 database-oracle bytes
with bound provenance; it makes no timing claim. Threshold checks use fresh
explicit oracles.

## Preserved history and certification boundary

The original `final_remaining.py` stops at N/P reuse. Its original strict
`resume_remaining_after_reuse.py` was not run. The separate
`continue_unaffected_integrated_checks.py` stops at original native X time
regressions. Neither is the active final continuation. C-X3 and C-X4 are rejected
isolated experiments; their unexecuted downstream plans are not acceptance
records. Native CPU profiles and six alignment diagnostic comparisons explain
why C-NA64 was investigated; they do not replace paired adoption measurements.

Frozen PR5 expected bytes and registered NCBI platform fingerprints are
unchanged. Corrected Sakai output differs from the old frozen LOSAT contract,
requiring separately reviewed authority/versioning for formal recertification.
The 46 deferred native cases remain deferred. No new biological search behavior
or parity exception is added. gbdraw was inspected, not modified or certified.
