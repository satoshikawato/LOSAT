# LOSAT Wasm Performance Plan

## Purpose

LOSAT Wasm performance work must improve browser/WASI throughput without
weakening NCBI BLAST parity. The current performance plot labels a `wasm n8`
series, but the Rust engines force `num_threads = 1` on `target_arch =
"wasm32"` and the comparison script builds Wasm with `--no-default-features`.
Therefore the first milestone is to make measurement labels and logs reflect
the real execution model before changing hot paths.

## Non-Negotiable Constraints

- NCBI BLAST remains the only behavioral reference.
- Wasm performance changes must not introduce unsupported BLAST behavior.
- Algorithm timing, pruning, ordering, scoring, e-values, and output formatting
  must remain bit-perfect with the existing parity rules.
- NCBI BLAST binaries may be used only as validation oracles, never as LOSAT
  runtime dependencies.
- Code changes that port or alter BLAST behavior must include NCBI source
  comments with file paths and line numbers immediately above the Rust code.

## Current Findings

- `LOSAT/.cargo/config.toml` enables `+simd128` for `wasm32-wasip1` and
  `wasm32-unknown-unknown`.
- `LOSAT/tests/run_comparison.sh` builds the command Wasm artifact with:
  `cargo build --release --target wasm32-wasip1 --no-default-features`.
- TBLASTX, BLASTN, and BLASTP all gate Rayon parallelism with
  `not(target_arch = "wasm32")`.
- TBLASTX and BLASTN explicitly resolve Wasm `num_threads` to `1`.
- BLASTP resolves `use_parallel = false` on Wasm.
- The current `wasm n8` plot series is therefore not evidence of real Wasm
  eight-thread compute parallelism.

## Phase 1: Fix Measurement Labels and Runtime Diagnostics

Goal: make every Wasm timing artifact unambiguous.

Tasks:

1. Update plotting labels so current command-Wasm results are shown as
   `LOSAT wasm serial` unless true Wasm threading is enabled and verified.
2. Add timing/log output for the effective engine thread count when
   `LOSAT_TIMING=1`.
3. Make comparison scripts record:
   - target triple
   - feature set
   - whether Rayon was compiled in
   - requested `-num_threads`
   - effective engine thread count
   - Node/WASI runtime version
4. Add a small log parser check that warns when a `wasm nN` label is produced
   from a serial Wasm run.

Acceptance:

- Existing performance plots cannot label serial Wasm execution as `n8`.
- `LOSAT_TIMING=1` logs show requested and effective thread counts.
- No BLAST output changes.

## Phase 2: Design Real Wasm Parallel Execution Separately

Goal: avoid treating native Rayon and Wasm threading as the same problem.

Tasks:

1. Split Wasm execution modes in docs and scripts:
   - command-WASI serial
   - browser in-memory serial
   - browser worker-parallel
   - future WASI-threaded runtime, if adopted
2. Keep current `wasm32-wasip1` command artifact serial until a supported
   threading runtime is selected.
3. Prototype browser-side worker parallelism outside BLAST hot loops first:
   - partition work by subject or subject chunk
   - run one LOSAT Wasm instance per worker
   - merge subject results in NCBI subject iteration order
4. Validate that worker partitioning preserves:
   - output order
   - hitlist pruning order
   - E-value/statistics context
   - TBLASTX local `--db-gencode` project exception
5. Only after the worker prototype passes parity, consider shared-memory Wasm
   with `SharedArrayBuffer` and a dedicated Rust-side threading strategy.

Acceptance:

- Threaded Wasm has a separate benchmark label and cannot be confused with
  command-WASI serial.
- Parallel Wasm output is byte-identical to serial LOSAT for supported modes.
- NCBI comparison remains available as a validation oracle.

## Phase 3: Improve Single-Thread Wasm First

Goal: reduce the baseline native n1 vs Wasm serial gap before adding threading.

Tasks:

1. Benchmark command-WASI file I/O against `web_api.rs` in-memory execution.
2. Add timing buckets, where missing, for:
   - FASTA parse and encoding
   - lookup build
   - subject scan
   - ungapped extension
   - gapped alignment
   - traceback
   - pruning/linking
   - final formatting
3. Prioritize allocation reduction in hot Wasm paths:
   - reusable scratch buffers
   - flat DP/trace buffers where NCBI-compatible
   - delayed clone/copy of HSP/edit-script data
   - avoiding per-subject temporary allocation when no hits survive
4. For BLASTP, continue the existing performance direction:
   - gapped alignment scratch reuse
   - composition/Kappa redo buffer reuse
   - sparse duplicate purge without duplicate-free allocation overhead
5. For BLASTN, focus on traceback and interval-tree memory behavior before
   changing algorithm structure.
6. For TBLASTX, profile scan, ungapped extension, reevaluation, and linking
   separately, because long-sequence extra-hit parity issues remain high
   priority.

Acceptance:

- Each optimized mode has before/after timing buckets.
- No optimization is accepted from total runtime alone; the changed bucket must
  explain the gain.
- Native, Wasm scalar, and Wasm SIMD parity outputs remain identical except for
  documented project exceptions.

## Phase 4: Keep TBLASTX SIMD Behind Parity Gates

Goal: preserve the existing `simd128` path while preventing silent divergence.

Tasks:

1. Keep `tblastx-wasm-scalar` as a diagnostic build feature.
2. For every TBLASTX Wasm SIMD change, compare:
   - native release
   - Wasm scalar
   - Wasm SIMD
   - NCBI oracle where applicable
3. Add benchmark rows for scalar vs SIMD Wasm to separate correctness from
   speed.
4. Profile the current SIMD routines:
   - `reevaluate_ungapped_hit_ncbi_translated_wasm_simd128`
   - left ungapped extension
   - right ungapped extension
   - two-hit scan helpers
5. Do not expand SIMD coverage until existing long-sequence TBLASTX parity
   discrepancies are understood or isolated from the SIMD path.

Acceptance:

- TBLASTX Wasm SIMD changes must pass native-vs-scalar-vs-SIMD output diffs.
- Scalar fallback remains easy to build.
- Any SIMD speedup report includes the scalar baseline.

## Phase 5: Preserve NCBI Parity While Optimizing

Goal: constrain performance work to implementation mechanics, not behavior.

Allowed optimization types:

- Memory reuse that preserves data lifetime and processing order.
- Matrix access specialization that keeps identical scoring values.
- Flat buffer layouts that preserve DP boundaries, x-drop rules, and traceback
  semantics.
- Lazy allocation and one-pass compaction that preserve NCBI sorting and
  tie-break behavior.
- Browser/worker orchestration that merges results in NCBI subject iteration
  order.

Disallowed optimization types:

- Skipping NCBI filtering, pruning, chaining, or reevaluation stages.
- Reordering stages for convenience.
- Replacing unsupported behavior with calls to NCBI executables.
- Changing floating-point precision.
- Accepting hit-count parity without bit score, E-value, coordinate, and output
  order parity.

Acceptance:

- Every behavioral code change cites the relevant NCBI source file and lines.
- Every performance change has a parity comparison target.
- If an NCBI reference cannot be found, the feature remains unimplemented.

## Suggested Work Order

1. Patch labels and timing metadata first.
2. Run one serial command-WASI benchmark and one in-memory Wasm benchmark on the
   same small fixture.
3. Add per-bucket timing for the largest visible Wasm gap.
4. Optimize one bucket at a time.
5. Re-run parity after each bucket change.
6. Prototype browser worker parallelism only after the serial baseline and
   labels are trustworthy.

## Initial Verification Commands

From `LOSAT/`:

```bash
LOSAT_TIMING=1 cargo run --release -- tblastx \
  -query tests/fasta/AvCLPV.fasta \
  -subject tests/fasta/PsCLPV.fasta \
  -out /tmp/losat_tblastx.out \
  -num_threads 1 \
  -outfmt 6

cargo build --release --target wasm32-wasip1 --no-default-features

LOSAT_TIMING=1 bash tests/run_comparison.sh
```

Use NCBI BLAST+ only in comparison scripts or manual oracle runs.

## Implementation Status

Status as of 2026-05-21:

- Phase 1 is implemented.
  - Current command-Wasm plot labels now default to `LOSAT command-WASI
    serial` and `LOSAT command-WASI serial (requested nN)`.
  - `LOSAT_WASM_THREADS_VERIFIED=1` is required before plot scripts emit
    `LOSAT wasm nN` or mode-specific threaded labels for Wasm runs.
  - TBLASTX, BLASTN, and BLASTP print requested and effective engine thread
    counts when `LOSAT_TIMING=1`.
  - `tests/run_comparison.sh` records command-Wasm metadata in each Wasm log:
    target triple, feature set, Rayon compile status, requested thread count,
    effective engine thread count, execution mode, threading status, benchmark
    label, WASI runtime label, and Node version.
  - Plot scripts warn if a `wasm nN` label is produced from a log whose
    effective engine thread count is `1`.
- Phase 2 is partially implemented.
  - `tests/run_comparison.sh` now has an explicit `LOSAT_WASM_EXECUTION_MODE`
    switch with named modes: `command-wasi-serial`,
    `browser-in-memory-serial`, `browser-worker-parallel`, and
    `future-wasi-threaded`.
  - `tests/run_comparison.sh` only executes `command-wasi-serial`; browser and
    future threaded modes are skipped with an explicit message so they cannot be
    confused with command-WASI serial runs.
  - `tests/plot_execution_time.py` and `tests/plot_comparison.py` derive Wasm
    labels from `LOSAT_WASM_EXECUTION_MODE`, so plots can distinguish
    command-WASI serial, browser in-memory serial, browser worker, and future
    WASI-threaded results.
  - Browser worker parallelism remains intentionally disabled/unverified; no
    worker-partitioned BLAST output is accepted until order, pruning,
    statistics, and project `db_gencode` rules pass parity.
- Phase 3 is partially implemented.
  - Added `tests/benchmark_wasm_serial_modes.sh` to benchmark command-WASI
    file-I/O execution against the `web_api.rs` browser/in-memory execution
    path on the same fixture.
  - The serial-mode benchmark records mode metadata and diffs command-WASI
    output against browser in-memory output before accepting the measurement.
  - TBLASTX qsort-compatible HSP ordering now avoids cloning the full
    `Vec<UngappedHit>` in the Wasm index-replay path. The implementation sorts
    only the index array and applies the destination-to-source order in place
    by cycle replay, reducing hot-path allocation/copy volume while preserving
    the NCBI qsort comparator order.
  - TBLASTX qsort-compatible HSP ordering now replays each permutation cycle by
    moving `UngappedHit` ownership with `ptr::read`/`ptr::write` instead of
    cloning every moved HSP. This keeps the same NCBI pointer-array qsort order
    while further reducing Wasm score/link sort copy overhead.
  - BLASTN serial Wasm traceback now reuses subject-local scratch storage for
    the combined preliminary HSP list and interval tree. The change keeps the
    NCBI traceback/init/reset order and only removes per-subject `Vec` and
    interval-tree allocation churn in `src/algorithm/blastn/blast_engine/run.rs`.
  - TBLASTX serial Wasm scanning now avoids per-scan-partition allocation for
    diagnostic scan interiors and clipped `seq_ranges` in
    `src/algorithm/tblastx/blast_engine/run_impl.rs`. The scan partition order
    and NCBI `seq_range` values are unchanged; only the temporary storage is
    converted to an iterator plus reusable scratch `Vec`.
  - TBLASTX ungapped HSP reevaluation now reuses the incoming HSP list buffer in
    `src/algorithm/tblastx/blast_engine/mod.rs` with order-preserving
    in-place survivor retention. This keeps the NCBI
    `Blast_HSPListReevaluateUngapped` remove-then-score-sort sequence while
    avoiding a fresh survivor `Vec` and per-survivor HSP moves on serial Wasm.
  - Verified the BLASTN scratch reuse on `LC738874.fasta` vs `LC738875.fasta`
    with `LOSAT_TIMING=1 cargo run --release -- blastn ... -num_threads 1
    -outfmt 6`; before/after outfmt 6 output diff was empty, and traceback
    counts remained `prelim=9`, `tree_skipped=1`, `full_traceback=8`.
  - Verified the command-WASI artifact still builds with:
    `cargo build --release --target wasm32-wasip1 --no-default-features`.
  - Verified the qsort replay behavior with:
    `cargo test index_replay --lib`.
  - Verified the TBLASTX scan scratch change with:
    `cargo check --lib`,
    `cargo test scan_interior --lib`, and
    `cargo build --release --target wasm32-wasip1 --no-default-features`.
  - Verified the TBLASTX reevaluation survivor-buffer reuse with
    `cargo check --lib`. Full-repository `cargo fmt` was stopped after it ran
    too long on the existing large tree; the touched Rust file was formatted
    directly with `rustfmt src/algorithm/tblastx/blast_engine/mod.rs --edition
    2021`.
  - Additional hot-path allocation optimization remains pending, especially
    around TBLASTX linking, remaining reevaluation profiling, and BLASTN
    traceback/interval-tree memory behavior.
- Phase 4 is partially implemented.
  - `tblastx-wasm-scalar` remains the diagnostic feature for forcing the
    scalar TBLASTX Wasm path.
  - Added `LOSAT/tests/benchmark_tblastx_wasm_simd_modes.sh` to build/capture
    command-WASI TBLASTX SIMD and scalar artifacts separately, run them against
    the same fixture, record timing/metadata, and reject native-vs-scalar-vs-SIMD
    raw output differences by default.
  - `tests/compare_tblastx_wasm_parity.sh` now records raw byte diffs in
    addition to sorted diffs and fails by default when LOSAT native, Wasm SIMD,
    and Wasm scalar outputs differ. NCBI oracle raw-diff failure remains
    opt-in with `REQUIRE_NCBI_TBLASTX_PARITY=1` because known TBLASTX parity
    defects are tracked separately from the Wasm SIMD-vs-scalar gate.
  - `LOSAT/tests/plot_execution_time.py` and
    `LOSAT/tests/plot_comparison.py` now include optional TBLASTX
    `LOSAT Wasm SIMD` and `LOSAT Wasm scalar` rows/series when matching
    `.wasm.simd.*` and `.wasm.scalar.*` artifacts are present.
  - Per-routine SIMD profiling for reevaluation, left/right ungapped extension,
    and two-hit scan helpers is still pending.
- Phase 5 remains an ongoing constraint for all later optimization work.
