# LOSATP Acceleration Plan

## Scope

This plan targets LOSATP / pure-Rust `blastp` performance while preserving
bit-perfect NCBI BLAST+ output parity. NCBI BLAST source remains the only
behavioral reference. No runtime path may call NCBI BLAST binaries, libraries,
FFI, or subprocesses.

Current LOSAT references:

- `LOSAT/src/algorithm/blastp/blast_engine.rs:2953-4298`: preliminary
  search setup, subject processing closure, parallel/serial subject execution,
  deterministic hitlist merge, timing instrumentation, and output formatting.
- `LOSAT/src/algorithm/blastp/blast_engine.rs:1184-1224`: scan wrapper around
  the amino-acid lookup scanner.
- `LOSAT/src/algorithm/blastp/extension.rs:1-357`: one-hit and two-hit
  ungapped extension helpers.
- `LOSAT/src/algorithm/blastp/gapalign.rs:77-128,1737-1852`: BLOSUM62 score
  dispatch, reusable gapped alignment scratch, and score-only preliminary
  gapped alignment.
- `LOSAT/src/algorithm/blastp/kappa.rs`: composition-based redo alignment.
- `LOSAT/tests/run_blastp_comparison.sh`: current tabular parity harness.
- `LOSAT/tests/run_blastp_threads_comparison.sh`: serial-vs-threaded BLASTP
  parity harness.
- `LOSAT/src/algorithm/tblastx/lookup/compressed.rs` and
  `LOSAT/tests/compressed_lookup.rs`: compressed amino-acid lookup construction
  scaffolding for future `blastp-fast` runtime integration.

Primary NCBI references:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554`:
  subject iteration, `s_BlastSearchEngineCore`, preliminary e-value reaping,
  and HSP stream write order.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1705`:
  setup of per-run scoring, extension, hit-saving, diagnostics, and
  `BlastGapAlignStruct`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-123`:
  `s_BlastAaScanSubject` emission order and offset-pair buffer behavior.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:440-588`:
  two-hit word finder and diagonal-state updates.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:713-800`:
  one-hit word finder and diagonal-state updates.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:831-957`:
  left/right ungapped extension loops.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1020-1177`:
  `s_BlastAaExtendOneHit` and `s_BlastAaExtendTwoHit`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3700-4007`:
  preliminary protein gapped extension loop, restricted alignment retry, and
  redo-index handling.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319`:
  `s_BlastProtGappedAlignment`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2430-2479`:
  BLASTP composition-adjusted preliminary redo parameters.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2939-3128`:
  `Blast_RedoAlignmentCore` and redone-match heap initialization.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:753-1355`:
  compressed amino-acid lookup table implementation for the future
  `blastp-fast` path.

## Baseline Observations

- LOSATP already has subject-level parallel execution under the default
  `parallel` feature. It is only used when `-num_threads` resolves to more
  than one and there are multiple subjects.
- The CLI default is still `-num_threads 1`, so existing comparison scripts
  exercise the serial path unless explicitly changed.
- The parallel path computes per-subject `BlastpSubjectResult`s, sorts them by
  `subject_index`, and then updates global preliminary hitlists. This is the
  right parity shape because it avoids Rayon scheduling order in final merge.
- `docs/multithreading_opportunities.md` is stale for BLASTP; it still marks
  subject-level BLASTP parallelism as not implemented.
- Pure-Rust BLASTP currently rejects `word_size != 3`, non-BLOSUM62/11/1
  scoring, unsupported composition modes, `--ungapped`, and `-use_sw_tback`.
  Those guardrails should remain until each NCBI behavior is ported.

## Progress Snapshot (2026-05-10)

This snapshot records the current state after implementing the safe portions of
Phases 0-4 under the parity gate requested for this acceleration work.

- Worktree state before this pass already included Phase 0/1 scaffolding,
  BLOSUM62 direct score helpers, and compressed lookup construction files.
  This pass extended Phase 2 scratch reuse, Phase 3 gapped-score dispatch,
  and Phase 3 ungapped-extension BLOSUM62 dispatch.
- Pre-change baseline from the current dirty worktree:
  - `cargo build`: passed.
  - `blastp` outfmt 6 vs NCBI references: 9 total diff lines across 4 cases:
    `AP027078.AP027131`, `AP027131.AP027133`, `AP027132.AP027133`,
    `AP027133.NZ_CP006932`.
  - `-num_threads 1` vs `-num_threads 8`: 0 mismatches.
  - threaded LOSATP vs NCBI: same 4 mismatch cases as the serial baseline.
- After Phase 2 scratch reuse:
  - `cargo build`: passed.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 9 total diff lines.
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases.
- After Phase 3 gapped BLOSUM62 dispatch:
  - `cargo build`: passed.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 9 total diff lines.
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases.
- After Phase 3 ungapped-extension BLOSUM62 dispatch:
  - Routed BLASTP one-hit/two-hit ungapped extension through BLOSUM62-specialized
    helper instantiations while preserving NCBI `aa_ungapped.c` loop order,
    x-drop checks, and scoring update timing.
  - `rustfmt src/algorithm/blastp/extension.rs`: passed.
  - `cargo test --lib test_extend`: 6 passed.
  - `cargo build`: passed.
  - `cargo build --release`: passed.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 9 total diff lines.
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases.
- After additional Phase 2 scratch reuse:
  - `cargo build --release`: passed.
  - Reused per-subject preliminary interval-tree storage with NCBI-equivalent
    query/subject bounds reset before preliminary gapped extension.
  - Reused BLASTP chaining node/keep buffers in subject scratch, matching
    NCBI's `gap_align->chaining->nodes` reuse shape.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 9 total diff lines.
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases.
- Phase 0 timing smoke run:
  - `LOSAT_TIMING=1` WSSV/PajaWSV serial and threaded outputs were
    byte-identical to the NCBI outfmt 6 reference.
  - Timing logs were written under `.codex_tmp/blastp_timing/`.
- Phase 4 construction-only tests:
  - `cargo test --test compressed_lookup`: 4 passed.
  - Runtime `blastp-fast` remains unsupported until compressed scan integration
    and full NCBI comparison cases are complete.
- Additional Phase 3 BLOSUM62 dispatch pass:
  - Routed Kappa redo identity/positive recounting through BLOSUM62-specialized
    helper instantiations while preserving NCBI `blast_hits.c` comparison order.
  - Routed preliminary gapped-start window scoring through a BLOSUM62-specialized
    helper while preserving NCBI `blast_gapalign.c` rolling-window order.
  - `cargo test --lib kappa`: 20 passed.
  - `cargo test --lib gapalign`: 21 passed, 1 ignored diagnostic.
  - `cargo build`: passed after each step.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 4 mismatch cases
    (`git diff --no-index` aggregate stayed 67 lines in this pass's metric).
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases after each step.
- Additional Phase 3 final-alignment recount pass:
  - Routed BLASTP alignment-view identity/positive recounting through
    BLOSUM62-specialized helper instantiations while preserving NCBI
    `blast_hits.c` match/positive comparison order.
  - Corrected the local alignment-view gapped unit test to use NCBISTDAA residue
    constants and NCBI `eGapAlignDel` gap-in-query semantics.
  - `rustfmt src/algorithm/blastp/alignment.rs`: passed.
  - `cargo test --lib alignment`: 26 passed.
  - `cargo build`: passed.
  - `blastp` outfmt 6 vs NCBI references: unchanged at 4 mismatch cases
    (`git diff --no-index` aggregate stayed 67 lines).
  - `run_blastp_threads_comparison.sh`: serial/threaded mismatches 0;
    threaded/NCBI mismatches unchanged at 4 cases.

## Phase 0: Benchmark and Timing Harness

Goal: identify where time is actually spent before changing hot loops.

Status: implemented and smoke-tested.

Implementation:

- Add `LOSAT_TIMING=1` guarded timing to `blastp` only, following the style
  already used in BLASTN/TBLASTX.
- Record at least these stages:
  - query encoding and SEG setup
  - subject encoding
  - lookup construction
  - subject scan plus ungapped extension
  - preliminary gapped extension
  - preliminary hitlist merge and sort
  - Kappa/composition redo
  - final traceback/stat refresh needed for output
  - output formatting
- Keep all timing writes on stderr and fully disabled when `LOSAT_TIMING` is
  unset.

Acceptance:

- `LOSAT_TIMING` disabled output is byte-identical to current output.
- One representative serial run and one multi-thread run have timing breakdowns
  attached in `tests/losat_out/` or an untracked temp directory.

Current result:

- Implemented in `blast_engine.rs` with `LOSAT_TIMING=1` guarded stderr output.
- Representative WSSV/PajaWSV serial and threaded timing logs were saved in
  `.codex_tmp/blastp_timing/`.
- Timing-disabled parity was checked indirectly by the full outfmt 6 comparison
  sweeps; diff lines did not increase from the pre-change baseline.

## Phase 1: Multi-Thread Parity Harness

Goal: make the existing subject-level parallel path easy to validate.

Status: implemented for outfmt 6.

Implementation:

- Add a comparison script, `LOSAT/tests/run_blastp_threads_comparison.sh`.
- For each existing BLASTP case, run LOSATP with `-num_threads 1` and
  `-num_threads 0` or a fixed count such as `-num_threads 8`.
- Diff serial LOSATP vs threaded LOSATP first, then diff threaded LOSATP vs the
  existing NCBI reference output.
- Add outfmt 0 and outfmt 7 variants if Phase 0 shows formatting or final
  traceback dominates.

Acceptance:

- `-num_threads 1` and `-num_threads N` produce byte-identical output for
  outfmt 6 on all current BLASTP comparison cases.
- Threaded LOSATP output remains byte-identical to NCBI references.
- Any failure is classified before optimization work continues.

Current result:

- Serial LOSATP and threaded LOSATP are byte-identical for all current outfmt 6
  BLASTP comparison cases.
- Threaded-vs-NCBI still has the same 4 mismatch cases already present in the
  serial baseline; no new threaded-only divergence was introduced.
- Outfmt 0 and outfmt 7 threaded variants are still future work.

## Phase 2: Thread-Local Scratch Reuse

Goal: remove allocation churn without changing algorithm timing or order.

Status: partially implemented and parity-gated.

Implementation:

- Extend `BlastpSubjectScratch` to own reusable per-subject buffers that are
  currently created inside `process_subject`.
- First candidates:
  - `Vec<InitHSP>`
  - `Vec<BlastpPreliminaryHsp>`
  - `Vec<Vec<BlastpPreliminaryHsp>>` or another deterministic split buffer for
    per-query preliminary hits
  - `GapAlignScratch` for preliminary score-only gapped alignment
  - `BlastIntervalTree` reset/rebuild support, if it can be reset without
    changing containment behavior
- Preserve NCBI order:
  - scan order from `s_BlastAaScanSubject`
  - init-HSP score sort before chaining/gapped extension
  - redo-index loop behavior from `blast_gapalign.c:3700-4007`
  - subject-result merge by increasing subject index

Acceptance:

- Serial and threaded BLASTP parity scripts pass.
- Timing shows reduced allocation-sensitive stages without changing hit counts,
  scores, e-values, coordinates, or output order.

Current result:

- Implemented reusable `BlastpSubjectScratch` storage for:
  - `Vec<InitHSP>`
  - BLASTP chaining nodes and keep/drop state
  - `Vec<BlastpPreliminaryHsp>`
  - `Vec<Vec<BlastpPreliminaryHsp>>`
  - restricted-alignment `found` / `restricted` buffers
  - preliminary `GapAlignScratch`
- `BlastIntervalTree` storage is reused per subject by resetting the same
  query/subject bounds that NCBI passes to `Blast_IntervalTreeInit`.
- Parity gate held: outfmt 6 NCBI differences did not increase from the
  recorded baseline, and serial/threaded mismatches remained 0.
- Further allocation profiling is still needed before declaring Phase 2 fully
  complete.

## Phase 3: BLOSUM62 Hot-Loop Specialization

Goal: reduce score lookup overhead in the currently supported default BLASTP
path.

Status: partially implemented and parity-gated.

Implementation:

- Keep generic scoring functions for future matrix support, but route the
  validated BLOSUM62/11/1 path through direct BLOSUM62 functions.
- Candidate sites:
  - `blastp_standard_score` in `blast_engine.rs`
  - `blastp_standard_score` and `BlastpScoreMatrix::score` in `gapalign.rs`
  - `residue_score` in `extension.rs`
  - `blastp_standard_score` in `kappa.rs`
- Do not change loop order, branch timing, early termination, x-drop logic, or
  integer/floating-point formulas.
- Add or update tests only around refactored helper behavior; broad comparison
  tests run after the sweep.

Acceptance:

- All BLASTP comparison scripts pass.
- Phase 0 timing shows a measurable reduction in scan/ungapped and/or gapped
  alignment stages.

Current result:

- Implemented direct BLOSUM62 score lookup in already-supported helper paths.
- Added a dedicated `BlastpScoreMatrix::Blosum62` dispatch path for gapped
  alignment so the default BLOSUM62/11/1 path avoids repeated generic matrix
  dispatch in hot DP loops.
- Added BLOSUM62-specialized instantiations for BLASTP one-hit/two-hit
  ungapped extension helpers in `extension.rs`, matching NCBI
  `aa_ungapped.c` matrix-access timing while removing the Rust-side matrix
  branch from the inner residue loop.
- Added BLOSUM62-specialized instantiations for Kappa redo identity/positive
  recounting and preliminary gapped-start window scoring, matching NCBI
  `blast_hits.c` and `blast_gapalign.c` matrix-access timing while removing
  Rust-side matrix branches from those inner residue loops.
- Added BLOSUM62-specialized instantiations for BLASTP alignment-view
  identity/positive recounting, matching NCBI `blast_hits.c` comparison order
  while removing the Rust-side matrix branch from that recount loop.
- Parity gate held: outfmt 6 NCBI diff count remained 9 total lines, and
  serial/threaded mismatches remained 0.
- Remaining Phase 3 candidates include scanner-side score paths and release
  timing measurements.
- Release-build timing comparisons are still needed before declaring a measured
  speedup.

## Phase 4: Compressed Amino-Acid Lookup for `blastp-fast`

Goal: unlock NCBI's `blastp-fast` path instead of trying to emulate it with
ad hoc filtering.

Status: lookup construction scaffolding implemented; runtime integration not
enabled.

Implementation:

- Port only from NCBI compressed lookup references:
  - `blast_aalookup.c:753-1355`
  - `blast_aalookup.h:186-322`
- Add Rust structures for compressed alphabet, backbone cells, overflow cells,
  and compressed index computation.
- Implement lookup construction before scan integration:
  - exact/neighbor word insertion
  - sorted compressed matrix traversal
  - overflow chain behavior
  - finalization layout
- Add a compressed scan path only after lookup table unit tests match NCBI.
- Remove or narrow the current `word_size != 3` rejection only for NCBI-backed
  supported cases.

Acceptance:

- Unit tests cover compressed index computation, overflow insertion order, and
  representative neighbor counts against NCBI-derived cases.
- `blastp-fast` comparison cases match NCBI before advertising support.
- Unsupported compressed cases still fail fast.

Current result:

- Construction primitives exist for compressed alphabet, backbone cells,
  overflow cells, compressed index computation, neighbor traversal, and PV
  finalization.
- `cargo test --test compressed_lookup` passes 4 construction-focused tests.
- Compressed subject scan integration is not complete.
- `blastp-fast` and `word_size > 4` runtime paths remain intentionally
  unsupported; do not narrow the guard until compressed scan and full comparison
  cases are complete.

## Phase 5: Kappa / Composition Redo Profiling

Goal: optimize or parallelize only if Phase 0 proves Kappa redo is a dominant
cost.

Status: not started.

Implementation:

- Keep this after Phase 1-3 because CBS redo is more parity-sensitive than
  subject-level preliminary search.
- If needed, evaluate batching by query or HSPList only against NCBI
  `Blast_RedoAlignmentCore` / `Blast_RedoAlignmentCore_MT` behavior.
- Do not share mutable composition workspace across workers.
- Do not change heap admission order, e-value comparison order, or local
  scaling factor calculations.

Acceptance:

- Any CBS optimization has dedicated serial-vs-threaded parity checks.
- WSSV/PajaWSV pairwise and tabular outputs remain byte-identical.

## Verification Commands

Run only after the relevant parity sweep is complete:

```bash
cd LOSAT && cargo build --release
cd LOSAT/tests && bash run_blastp_comparison.sh
cd LOSAT/tests && bash run_blastp_pairwise_comparison.sh
cd LOSAT/tests && bash run_blastp_outfmt7_comparison.sh
cd LOSAT/tests && bash run_blastp_custom_fields_comparison.sh
cd LOSAT/tests && bash run_blastp_threads_comparison.sh
```

For this pass, outfmt 6 and threaded comparison sweeps were run in debug builds
because the acceptance check was parity preservation, not final release timing.

## Risk Controls

- Do not introduce approximate algorithms, fast-math, FMA-dependent behavior,
  relaxed floating-point ordering, or threshold shortcuts.
- Do not parallelize inside one subject until NCBI order-sensitive diagonal and
  redo behavior is fully understood and tested.
- Do not remove unsupported-feature errors until the matching NCBI source path
  is ported and covered.
- Every code change must include an immediate NCBI reference comment with file
  path and line numbers.
- If a planned optimization cannot be tied to NCBI source, drop it.

## Recommended Order

1. Phase 0 timing instrumentation. Done for BLASTP.
2. Phase 1 multi-thread parity harness. Done for outfmt 6.
3. Phase 2 scratch reuse. Partially done; continue allocation profiling.
4. Phase 3 BLOSUM62 hot-loop specialization. Partially done; gapped alignment
   and ungapped extension dispatch are parity-gated. Continue scanner-side hot
   score paths, Kappa redo scoring paths, and release timing.
5. Phase 4 compressed lookup / `blastp-fast` port. Construction tests pass;
   compressed scan and runtime support remain.
6. Phase 5 CBS redo work only if profiling justifies it. Not started.
