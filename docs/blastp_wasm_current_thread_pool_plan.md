# BLASTP Wasm Current-Thread Pool Plan

## Summary

LOSAT threaded command-Wasm currently pays one extra worker for BLASTP
parallel stages. A `-num_threads N` run starts the command Wasm instance on one
JavaScript worker, then Rayon creates `N` additional WASI threads. The command
worker mostly waits while the Rayon workers do the BLASTP work, so the browser
or Node runtime sees `N + 1` workers.

The implementation goal is to make BLASTP `-num_threads N` complete with exactly
`N` JavaScript workers in threaded Wasm:

- one command/leader worker that also participates as Rayon worker index 0
- `N - 1` additional WASI thread workers spawned through `wasi/thread-spawn`

This plan is BLASTP-first. BLASTN and TBLASTX can use the same helper later,
but their code paths should not be changed until BLASTP is verified.

## Constraints

- Preserve NCBI BLASTP output parity. Parallel execution may move independent
  computation earlier, but hit-list admission, reduction, pruning, scoring, and
  output ordering must remain the existing deterministic LOSAT order that
  replays NCBI subject/query order.
- Do not call NCBI BLAST binaries, libraries, FFI, or subprocesses from runtime
  code. NCBI is reference-only.
- Every Rust or JavaScript implementation change must keep the required NCBI
  reference comments immediately above the changed logic.
- Keep plain `wasm32-wasip1` command Wasm serial. The change applies only to
  `wasm32-wasip1-threads --features wasm-threads`.
- Keep native Rayon behavior unchanged unless a native-only compile fix is
  required.

## NCBI References To Keep On Code Changes

- `ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222`
  caps `num_threads` to `CSystemInfo::GetCpuCount()` and stores the effective
  thread count.
- `ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88`
  calls `SetNumberOfThreads(num_threads)` when `num_threads > 1`.
- `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1554`
  iterates subject OIDs, runs the BLAST search core per subject, and writes HSP
  lists to the stream.
- `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1688`
  sets up per-thread BLAST preliminary-search state before MT execution.
- `ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3429-3459`
  runs composition-based redo work in an OpenMP static loop using the effective
  thread count.
- `ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3817-3850`
  consolidates thread-local composition-adjustment results before continuing.

## Current BLASTP Parallel Pool Sites

The first implementation sweep should touch these BLASTP pool creation sites:

- `LOSAT/src/algorithm/blastp/blast_engine.rs:671`
  `prepare_blastp_subjects_preserving_order()`, subject encoding.
- `LOSAT/src/algorithm/blastp/blast_engine.rs:5348`
  preliminary subject search.
- `LOSAT/src/algorithm/blastp/blast_engine.rs:5659`
  Kappa single-query match redo.
- `LOSAT/src/algorithm/blastp/blast_engine.rs:5778`
  Kappa query redo.

Do not add `.use_current_thread()` independently at each site. Rayon
`ThreadPoolBuilder::use_current_thread()` marks the current thread as worker
index 0 for that pool, and local current-thread pools are not suitable for
repeated creation on the same command worker. BLASTP needs one reusable per-run
pool, created only if a parallel BLASTP stage actually runs.

## Design

### 1. Add A BLASTP Per-Run Pool Helper

Add a small BLASTP-local helper near the existing BLASTP threaded-WASI helpers:

- `build_blastp_thread_pool(num_threads) -> Result<rayon::ThreadPool>`
- `blastp_thread_pool_spawned_workers(num_threads) -> usize` for diagnostics
  if useful

The helper should build:

```rust
let builder = rayon::ThreadPoolBuilder::new().num_threads(num_threads);

#[cfg(all(target_arch = "wasm32", feature = "wasm-threads"))]
let builder = builder.use_current_thread();

builder.build()
```

Expected behavior:

- native: `num_threads N` spawns/uses the same Rayon-managed worker pool as now
- threaded Wasm: `num_threads N` creates a Rayon pool of N workers, with the
  command worker as index 0 and only `N - 1` `wasi/thread-spawn` calls

### 2. Add Lazy Per-Run Pool State

Inside the BLASTP run path, keep:

```rust
let mut blastp_parallel_pool: Option<rayon::ThreadPool> = None;
```

Then use a small local function or closure:

```rust
let pool = get_or_build_blastp_pool(&mut blastp_parallel_pool, num_threads)?;
pool.install(|| {
    // existing par_iter work
});
```

This must be lazy. Existing threaded-WASI gates deliberately keep small jobs
serial; those serial jobs should not create a worker pool at all.

### 3. Reuse The Same Pool Across BLASTP Stages

Replace the four direct `ThreadPoolBuilder::new()` sites with the shared pool:

- subject preparation uses the shared pool only when
  `blastp_should_parallelize_subject_prep()` returns true
- preliminary search uses the shared pool when `use_parallel` is true
- Kappa match redo uses the shared pool when `use_parallel_kappa_match_redo` is
  true
- Kappa query redo uses the shared pool when `use_parallel_kappa_redo` is true

The same effective thread count should be used for all BLASTP stages in one
run. If a future stage needs a different effective count, do not build a second
current-thread pool in the same Wasm command worker; finish the first design
first and document the new requirement separately.

### 4. Preserve Reduction Order

Keep all existing ordered reducers:

- preliminary search still writes one `BlastpSubjectResult` per subject index
  and replays those slots in subject order
- Kappa match redo still precomputes expensive results in parallel, then
  replays heap admission in original local-match order
- Kappa query redo still writes `(q_idx, heap)` results and assigns them by
  query index
- subject preparation still collects prepared subjects in input order

No parallel stage should write directly to the final output or shared hit-list
heap in nondeterministic worker completion order.

### 5. Update Diagnostics

Extend the existing `LOSAT_WASI_THREADS_DEBUG=1` BLASTP diagnostics so threaded
Wasm can distinguish:

- `rayon_pool_threads=N`
- `expected_spawned_workers=N-1` when `.use_current_thread()` is active
- stage-level `parallel=true/false`
- `serial_reason` unchanged

The Node runner already logs actual `spawn tid=...` events. Acceptance should
compare the Rust expectation against runner logs.

## Implementation Steps

1. Add `build_blastp_thread_pool()` behind the same cfg used by the existing
   BLASTP parallel code.
2. Add lazy per-run pool state in the BLASTP run function after the effective
   `num_threads` calculation.
3. Change `prepare_blastp_subjects_preserving_order()` to optionally accept a
   pool reference, or split its parallel body so subject preparation can be run
   from the caller's shared pool.
4. Replace the preliminary-search pool creation with the shared pool.
5. Replace both Kappa redo pool creations with the shared pool.
6. Adjust threaded-WASI diagnostic text to report expected spawned workers.
7. Keep all non-threaded Wasm cfg branches returning serial behavior.
8. Run formatting and focused builds.

## Validation

Minimum validation after implementation:

```bash
cd LOSAT
cargo check
cargo check --target wasm32-wasip1-threads --features wasm-threads
cargo build --release --target wasm32-wasip1-threads --features wasm-threads
node --check tests/run_losat_wasi_threads.js
```

Focused behavioral validation:

```bash
cd LOSAT
LOSAT_WASI_THREADS_DEBUG=1 node tests/run_losat_wasi_threads.js \
  target/wasm32-wasip1-threads/release/LOSAT.wasm \
  blastp -query tests/fasta/WSSV.faa -subject tests/fasta/PajaWSV.faa \
  -out /tmp/losatp.wasm.n8.out -outfmt 6 -num_threads 8
```

Expected diagnostics for `-num_threads 8`:

- Rust reports `effective_threads=8`
- Rust reports `rayon_pool_threads=8`
- Rust reports `expected_spawned_workers=7`
- runner logs show seven `spawn tid=...` lines for the first parallel stage
- output remains byte-identical to the previous `n8` threaded-Wasm output and
  to `n1` where the existing BLASTP parity checks require it

Also run at least one small BLASTP fixture below the threaded-WASI thresholds.
It should remain serial and show no worker spawn.

## Risks And Mitigations

- Current-thread local pools can leak Rayon registry state. Mitigation: create
  at most one shared BLASTP pool per Wasm command run, lazily, and reuse it
  across all BLASTP parallel stages.
- Calling `pool.install()` from a thread that is already worker index 0 should
  execute directly inside the same registry. Mitigation: test a BLASTP run that
  uses both preliminary search and a Kappa redo path in the same process.
- Subject preparation currently owns its own pool inside a helper. Mitigation:
  pass the shared pool from the caller or keep subject preparation serial in the
  first patch if the lifetime change becomes too invasive.
- The change reduces worker count but may not improve small jobs. Mitigation:
  keep the existing threaded-WASI stage gates and validate no-spawn small jobs.

## Rollout Order

1. Preliminary search only, with subject preparation forced through the existing
   serial path if needed. Verify `N-1` spawns and output identity.
2. Move subject preparation onto the shared pool.
3. Move Kappa match redo onto the shared pool.
4. Move Kappa query redo onto the shared pool.
5. Remove any temporary serial fallback added only to simplify step 1.

The final state should have no BLASTP direct `ThreadPoolBuilder::new()` calls
outside `build_blastp_thread_pool()`.
