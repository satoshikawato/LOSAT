# Wasm Multithreading Improvement Plan

## Summary

LOSAT's Wasm multithreading currently has two separate problems:

- Plain `wasm32-wasip1` builds are intentionally serial. Real threading requires
  the `wasm32-wasip1-threads` target, `--features wasm-threads`, and the
  `tests/run_losat_wasi_threads.js` runner.
- Some threaded runs still do little or no parallel work. BLASTP can clamp the
  requested thread count back to one under WASI, while BLASTN/TBLASTX often have
  too few subject or chunk work items for the worker pool to help.

The first implementation goal is to make `-num_threads N` actually create and
use N workers in threaded Wasm when the algorithm has enough independent work.
The second goal is to widen the independent work units without changing NCBI
BLAST output ordering, scoring, filtering, pruning, or formatting.

## NCBI References

Every code change from this plan must include the local NCBI reference comment
immediately above the modified Rust or JavaScript logic.

- `ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222`
  caps `num_threads` to `CSystemInfo::GetCpuCount()` and stores the effective
  thread count.
- `ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3224-3228`
  warns that `num_threads` is ignored when local `-subject` is used.
- `ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88`
  calls `SetNumberOfThreads(num_threads)` when `num_threads > 1`.
- `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1475`
  iterates subject sequences, runs the search core per subject, and writes HSP
  lists through the NCBI stream model.
- `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584`
  splits long subject sequences into `MAX_DBSEQ_LEN` chunks and merges chunk
  HSP lists in order.
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-614`
  performs TBLASTX amino-acid subject scanning and advances the two-hit
  diagonal table state after the scan unit.
- `ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-3035`
  merges HSP lists for split subject chunks.

## Current LOSAT Findings

- `LOSAT/Cargo.toml:37-40` documents that real Wasm threading is opt-in only.
- `LOSAT/.cargo/config.toml:26` defines `build-web-threaded` for
  `wasm32-wasip1-threads --features wasm-threads`.
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:565-585` and
  `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4361-4381` force
  `num_threads = 1` for non-threaded Wasm builds.
- `LOSAT/src/algorithm/blastp/blast_engine.rs:950-956` clamps the requested
  BLASTP thread count to `num_cpus::get()`. Under WASI this can report one CPU,
  so `-num_threads 8` can become one thread before Rayon is reached.
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:2742-2761` parallelizes
  TBLASTX primarily over `subjects_raw.par_iter()`. A local `-subject` FASTA
  with one record gives the pool only one subject-level work item.
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1374-1395` currently
  walks scan interiors serially. The `LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS` path
  changes chunk sizing but does not itself parallelize those interiors.
- `tests/run_losat_wasi_threads.js:166-212` implements the WASI
  `thread-spawn` import using Node workers. This is the correct runner for
  threaded command Wasm.

## Phase 0: Measurement Guardrails

Add a repeatable measurement path before changing algorithms.

- Add a diagnostic mode, gated behind an environment variable, that prints:
  effective thread count, target family, threaded-Wasm feature state, Rayon pool
  size, number of subject records, number of subject chunks, and number of
  worker jobs submitted.
- Keep diagnostics off by default to preserve output parity.
- Extend `tests/run_wasm_comparison.sh` or add a focused script that verifies:
  the threaded artifact exports `_start` and `wasi_thread_start`, imports
  `wasi/thread-spawn`, and actually spawns workers for `-num_threads N`.
- Record `real/user/sys` plus output hash for native `n1`, native `nN`, Wasm
  serial, and Wasm threaded runs.

Acceptance for Phase 0:

- A BLASTP threaded-Wasm `-num_threads 8` run reports whether it was capped to
  one before Rayon.
- BLASTN and TBLASTX threaded-Wasm runs report how many useful jobs were
  submitted, not just how many workers were spawned.
- All diagnostics remain stderr-only and disabled by default.

## Phase 1: Fix BLASTP Effective Thread Count Under WASI

Problem: `blastp_effective_num_threads()` uses `num_cpus::get()` as a hard cap.
In threaded WASI this value can be one even though the runner can create worker
threads.

Plan:

1. Update `tests/run_losat_wasi_threads.js` to provide a runner thread-cap env
   value, for example `LOSAT_WASI_THREAD_CAP`, defaulting to
   `os.availableParallelism()` when available.
2. Update BLASTP's effective-thread calculation so that, only under
   `all(target_arch = "wasm32", feature = "wasm-threads")`, it uses the runner
   cap instead of `num_cpus::get()`.
3. Keep native behavior unchanged: native BLASTP should still cap requested
   threads to `num_cpus::get()` to match the NCBI rule.
4. Preserve deterministic reduction order. Worker results must continue to be
   stored in subject index order and replayed in NCBI subject iteration order.

Acceptance for Phase 1:

- `LOSAT_WASI_THREADS_DEBUG=1 ... blastp -num_threads 8` shows worker spawn.
- BLASTP `n1` and `n8` Wasm output files are byte-identical for the focused
  fixtures.
- Multi-subject protein fixtures show `user > real` for threaded Wasm and a
  measurable wall-time improvement over Wasm `n1`.

## Phase 2: Reduce Threading Overhead For Small Wasm Jobs

Problem: Node worker startup and Wasm instance initialization can dominate short
runs. For small BLASTN/BLASTP/TBLASTX cases, threaded Wasm can be equal or
slower even when workers are created correctly.

Plan:

1. Add per-program work-size thresholds for threaded Wasm only:
   subject count, total residues/bases, estimated chunk count, and query count.
2. Route tiny jobs to the serial path even if `-num_threads > 1`.
3. For BLASTN, avoid the writer thread in threaded Wasm when the job count is
   below the threshold; reduce in-process using the same final ordering.
4. Keep thresholds conservative and diagnostic-visible so they can be tuned from
   benchmark data rather than guesses.

Acceptance for Phase 2:

- Small Wasm `n8` fixtures are not slower than Wasm `n1` by more than a narrow
  fixed-overhead tolerance.
- Larger fixtures still enter the parallel path.
- Output remains byte-identical across `n1` and `nN`.

## Phase 3: Increase BLASTN Useful Parallel Work

Problem: BLASTN currently parallelizes subject records and long subject chunks.
A single local subject under `MAX_DBSEQ_LEN` often provides only one heavy work
item.

Plan:

1. First validate the existing NCBI subject-split path on sequences above
   `MAX_DBSEQ_LEN`, including merge ordering and common-endpoint pruning.
2. If the long-subject path is correct, make its chunk batching more explicit in
   diagnostics and avoid batching choices that leave workers idle.
3. Investigate query-context or range partitioning only after mapping the
   corresponding NCBI traceback and filtering order. Do not introduce a new
   partitioning strategy until the NCBI-equivalent merge, purge, and sorting
   points are documented with source line references.

Acceptance for Phase 3:

- Long BLASTN subjects show actual chunk-level parallelism in threaded Wasm.
- `n1` and `nN` outputs are byte-identical.
- Any new partitioning has a documented NCBI reference and parity test before it
  becomes default.

## Phase 4: Increase TBLASTX Useful Parallel Work

Problem: TBLASTX speedup is currently limited for local single-subject inputs.
Subject-level parallelism helps only when there are multiple subject records.
The current scan-chunk experiment is serial inside the scan loop and disables
subject-level parallelism while active.

Plan:

1. Treat `LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS` as diagnostic only until it has a
   real parallel reducer and parity proof.
2. Validate `LOSAT_TBLASTX_PARALLEL_CHUNKS` on sequences that actually split
   into multiple `MAX_DBSEQ_LEN` chunks. If parity holds, consider enabling it
   by default for threaded Wasm and native threaded builds.
3. Investigate frame-level parallelism only after checking NCBI context/frame
   ordering in `blast_engine.c` and `aa_ungapped.c`. The two-hit diagonal state
   and `Blast_ExtendWordExit` timing are the main risks.
4. Keep final merge ordered by NCBI subject/frame/chunk iteration order before
   chaining, e-value calculation, culling, and output formatting.

Acceptance for Phase 4:

- Existing TBLASTX parity fixtures stay byte-identical for `n1` and `nN`.
- Long split-subject fixtures show worker utilization in threaded Wasm.
- For short single-subject TBLASTX, the expected no-speedup case is explicitly
  reported by diagnostics until a parity-proven frame or scan partition exists.

## Phase 5: Test Matrix

Run tests only after the corresponding implementation phase is complete.

- Build checks:
  - `cargo build --release --target wasm32-wasip1 --no-default-features`
  - `cargo build --release --target wasm32-wasip1-threads --features wasm-threads`
- Artifact checks:
  - serial Wasm must not import `wasi/thread-spawn`
  - threaded Wasm must export `wasi_thread_start`
  - threaded Wasm must import `wasi/thread-spawn`
- Parity checks:
  - BLASTP: protein fixtures from `LOSATP_CASES` in `tests/run_wasm_comparison.sh`
  - BLASTN: megablast and blastn fixtures from `tests/run_wasm_comparison.sh`
  - TBLASTX: short viral fixtures and at least one long split-subject fixture
- Timing checks:
  - compare native `n1`/`nN`, Wasm serial, and Wasm threaded
  - record effective thread count and submitted work items
  - compare output hashes before using timing numbers

## Non-Goals

- Do not call NCBI BLAST binaries, libraries, FFI, or subprocesses from LOSAT
  runtime code.
- Do not add a fallback that delegates unsupported behavior to NCBI.
- Do not change output order to gain speed.
- Do not make scan or frame partitioning default until the NCBI-equivalent
  merge and filtering order is documented and tested.

## Immediate Next Patch

The first code patch should be small and focused:

1. Add threaded-WASI effective-thread diagnostics.
2. Have `tests/run_losat_wasi_threads.js` expose a runner CPU/thread cap.
3. Change only BLASTP's threaded-Wasm effective-thread cap to use that runner
   cap.
4. Verify worker spawn and byte-identical BLASTP output for `n1` and `n8`.

This patch addresses the clearest current failure mode before touching BLASTN
or TBLASTX algorithm partitioning.

## 2026-05-24 Implementation Update

Implemented the Phase 1 BLASTP threaded-WASI fix as a real execution-path
performance change, not as an external harness-only measurement change.

Changes made:

- `LOSAT/tests/run_losat_wasi_threads.js` now publishes
  `LOSAT_WASI_THREAD_CAP` into the WASI environment. The default is Node's
  `os.availableParallelism()` when available, falling back to `os.cpus().length`.
- `LOSAT/src/algorithm/blastp/blast_engine.rs` now uses that runner-provided
  cap only for `all(target_arch = "wasm32", feature = "wasm-threads")`.
  Native BLASTP still caps with `num_cpus::get()` to preserve the NCBI
  `CSystemInfo::GetCpuCount()` rule.
- BLASTP threaded-WASI diagnostics now report the requested thread count,
  runner cap, effective thread count, Rayon pool size, subject record count,
  submitted subject jobs, and whether the prelim search entered the parallel
  path. Diagnostics remain stderr-only and disabled unless
  `LOSAT_WASI_THREADS_DEBUG=1` or `true`.

Validation performed:

- `cargo check`
- `rustfmt src/algorithm/blastp/blast_engine.rs --edition 2021`
- `cargo check --target wasm32-wasip1-threads --features wasm-threads`
- `cargo build --release --target wasm32-wasip1-threads --features wasm-threads`
- `node --check LOSAT/tests/run_losat_wasi_threads.js`
- Threaded Wasm BLASTP fixture:
  `WSSV.faa` vs `PajaWSV.faa`, `-outfmt 6`, `-num_threads 1` and
  `-num_threads 8`.

Observed threaded-WASI behavior for the fixture:

- Debug output reported `runner_thread_cap=32`, `requested_threads=8`,
  `effective_threads=8`, `rayon_pool_threads=8`, `subject_records=86`,
  `worker_jobs=86`, and `parallel=true`.
- The threaded runner spawned workers instead of clamping BLASTP back to one
  effective thread.
- `n1` and `n8` output hashes were identical:
  `ade827538e1f218b8c4b5dce88744d8495e58335144f7259bdb529801d220d10`.
- One timing sample with debug disabled:
  - `n1`: `real 4.17`, `user 4.48`, `sys 0.06`
  - `n8`: `real 2.00`, `user 7.40`, `sys 0.20`

Notes:

- Whole-crate `cargo fmt` was not used to completion because rustfmt spent
  several minutes formatting generated/large crate-root sources unrelated to
  this patch. The edited Rust file was formatted directly with rustfmt.
- The release build still reports Cargo's existing lib/bin Wasm output filename
  collision warning; this patch did not change target naming.

Next policy:

- Treat this as the baseline for Phase 2. Do not add more scripts before
  improving execution behavior.
- Add conservative threaded-WASI small-job thresholds for BLASTP first, using
  subject count and total subject residues. The goal is to avoid worker
  startup overhead on tiny fixtures while keeping larger multi-subject searches
  on the parallel path.
- Extend the same diagnostic shape to BLASTN and TBLASTX only when it supports
  a concrete scheduling decision or performance change.
- For BLASTN, continue with existing NCBI subject-chunk semantics before
  considering any query/range partitioning.
- For TBLASTX, keep scan/frame partitioning diagnostic-only until the
  NCBI-equivalent merge/filter timing is documented and parity-proven.

## 2026-05-24 Phase 2 Implementation Update

Implemented the first Phase 2 BLASTP threaded-WASI scheduling gate as an
execution-path performance change.

Changes made:

- `LOSAT/src/algorithm/blastp/blast_engine.rs` now applies a threaded-WASI-only
  parallelization decision before creating Rayon pools for BLASTP subject
  preparation, preliminary subject search, Kappa query redo, and Kappa
  single-query match redo.
- Native behavior is unchanged: native parallel builds still use the existing
  NCBI-style effective thread count and enter the parallel path whenever there
  is more than one worker job.
- Non-threaded command Wasm behavior is unchanged: `wasm32-wasip1` without
  `wasm-threads` remains serial.
- The new threaded-WASI defaults are conservative:
  - `LOSAT_BLASTP_WASI_MIN_WORK_ITEMS_PER_THREAD=2`
  - `LOSAT_BLASTP_WASI_MIN_SUBJECT_RESIDUES_PER_THREAD=4096`
- The thresholds are environment-overridable for benchmark tuning. Setting both
  values to `0` approximates the previous threaded-WASI scheduling behavior
  except for inherently single-work-item stages.
- BLASTP threaded-WASI diagnostics now include
  `total_subject_residues`, `min_worker_jobs`, `min_subject_residues`, and
  `serial_reason` for each gated stage.

Validation performed:

- `rustfmt src/algorithm/blastp/blast_engine.rs --edition 2021`
- `cargo check`
- `cargo check --target wasm32-wasip1-threads --features wasm-threads`
- `cargo build --release --target wasm32-wasip1-threads --features wasm-threads`
- Node artifact inspection confirmed the threaded artifact exports `_start` and
  `wasi_thread_start`, and imports `wasi/thread-spawn`.
- `cargo test test_merge_blastp_subject_results_replays_subject_index_order`
- Threaded Wasm BLASTP fixtures:
  - small: `WSSV.faa` vs the first `PajaWSV.faa` subject record,
    `-outfmt 6`, `-num_threads 1` and `-num_threads 8`
  - larger: `WSSV.faa` vs full `PajaWSV.faa`, `-outfmt 6`,
    `-num_threads 1` and `-num_threads 8`

Observed threaded-WASI behavior:

- Small fixture with normal thresholds:
  - `subject-prep`: `parallel=false`, `serial_reason=worker_jobs<=1`
  - `prelim`: `parallel=false`, `serial_reason=worker_jobs<=1`
  - `kappa-query-redo`: `parallel=false`,
    `serial_reason=subject_residues<threshold`
  - no worker spawn occurred for the `-num_threads 8` run
  - `n1`: `real 0.32`, `user 0.49`, `sys 0.02`
  - `n8`: `real 0.33`, `user 0.44`, `sys 0.06`
  - `n8` with thresholds disabled: `real 0.59`, `user 0.84`, `sys 0.04`
  - `n1`, `n8`, and threshold-disabled `n8` output hashes were identical:
    `0bd7c066d9bf5fe302be681797ac4bb73f865c7a0e57d0e4293747c633f652f6`
- Larger fixture:
  - `subject-prep`, `prelim`, and `kappa-query-redo` all reported
    `parallel=true`
  - debug output showed worker spawn for the gated parallel stages
  - `n1`: `real 4.31`, `user 4.67`, `sys 0.03`
  - `n8` with debug disabled: `real 2.11`, `user 7.70`, `sys 0.13`
  - `n1` and `n8` output hashes were identical:
    `ade827538e1f218b8c4b5dce88744d8495e58335144f7259bdb529801d220d10`

Notes:

- The checked commands still emit the repository's existing warning set,
  including the Cargo lib/bin Wasm output filename collision warning. This
  patch did not change those unrelated warnings.
- The new thresholds are not a parity mechanism; they only decide whether to
  pay the threaded-WASI worker startup cost. Final reduction order and output
  are unchanged.

Next policy:

- Tune the BLASTP threaded-WASI thresholds against a broader fixture matrix,
  keeping output-hash checks paired with every timing sample.
- Consider reusing or sharing a threaded-WASI Rayon pool across BLASTP stages.
  The larger fixture still creates separate worker sets for subject preparation,
  preliminary search, and Kappa redo, so pool startup overhead remains visible.
- Add BLASTN and TBLASTX scheduling gates only when tied to concrete existing
  subject/chunk work units and NCBI-equivalent merge ordering.
- Do not introduce BLASTN query/range partitioning or TBLASTX frame/scan
  partitioning until the corresponding NCBI merge, purge, filtering, and output
  ordering points are fully documented and parity-proven.
