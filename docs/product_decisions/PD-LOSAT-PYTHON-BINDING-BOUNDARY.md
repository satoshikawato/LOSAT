# Product Decision: LOSAT Python binding boundary

- Decision ID: `PD-LOSAT-PYTHON-BINDING-BOUNDARY`
- Version: 1.0
- Date: 2026-11-05
- Status: The maintainer decided the placement, the rules and the name on 2026-11-05 (decisions 2,
  3 and 6 of [`docs/losat_library_surface_plan.md`](../losat_library_surface_plan.md) §10.1: outside
  `LOSAT/`, the rules of `web/`, PyPI name `pylosat`). This text takes effect with the maintainer's
  merge.

## Scope

`pylosat` is the Python binding of `LOSAT::api`. This decision governs all code under `python/`:
the Rust crate that PyO3 builds into a Python extension, the Python package, its tests and tools.

It does not change the rules in the root [`AGENTS.md`](../../AGENTS.md) for code under `LOSAT/`. A
change under `LOSAT/` made for the binding still requires NCBI source authority, NCBI reference
comments and the parity gates.

## Authorities

- NCBI BLAST+ source remains the only authority for search behavior and for the compatibility
  outputs (`AGENTS.md`, [`PD-PURE-RUST-RUNTIME-AUTHORITY`](PD-PURE-RUST-RUNTIME-AUTHORITY.md)).
- `LOSAT::api` is the only way into the engine (`docs/losat_library_surface_plan.md`). A value that
  the binding lacks is added there, under the rules of `LOSAT/`, and not computed here.

## Decision

1. **Placement.** `python/` is its own Cargo build root and Python project. Nothing under `LOSAT/`
   depends on it, so Arrow, PyO3 and pickle do not enter the engine or its dependency tree.
2. **What has no NCBI equivalent.** Python's shapes of the API (`Search` keywords, `Sequences`,
   iterators, `HspTable`), the Arrow table, pickle and type stubs have no NCBI equivalent. They are
   permitted, but only under `python/` and only while the invariants below hold.
3. **Invariants.**
   1. **Engine input.** The engine receives exactly what `LOSAT::api` takes from any caller: an argv
      that the command line's parser accepted, and FASTA bytes read by the command line's reader.
      `Search` keywords become those words and nothing else. `Sequences.from_pairs` writes whole
      FASTA records. The binding never changes sequences, options or statistics, and never splits
      one search into independent searches whose results are joined. `Search.batches` is NCBI's
      own loop over query batches (`LOSAT::api::run_batches`), not such a split.
   2. **Compatibility outputs.** Reports are produced only by the engine's formatters. The binding
      returns them byte for byte and never rewrites, filters or regenerates them.
   3. **BLAST-defined values.** Code under `python/` does not compute or format scores, E-values,
      identities, coverages, NCBI number formats or alignment text. Every value of `Hsp`, `HspTable`
      and `QueryStats` is the `LOSAT::api` value, not rounded.
   4. **Labels.** A column named as an `-outfmt` column holds that column's value for the row, with
      the difference stated in the API's documentation: `qframe` and `sframe` are null where the
      program does not translate that side. The table and pickle are binding formats and never
      claim identity with NCBI output.
   5. **Pure Rust at run time.** The extension calls no NCBI executable or library and starts no
      process. Apart from Python's own C API its runtime dependencies stay free of C libraries:
      Arrow's optional `zstd` stays off.
   6. **Process-wide state.** The binding does not call Rayon's `build_global` or change the
      process's environment, current directory or panic hook. Each search keeps its own thread
      pool, as `LOSAT::api` documents.
   7. **Research data stays local.** Sequences, headers, file names and results are never sent over
      the network.
4. **Stopping.** Python cannot stop a thread inside its process, so a search stops only where the
   API lets a caller stop: between two query batches. The GIL is released while a search runs.
   Stopping inside a batch needs NCBI's `SetInterruptCallback` in the engine, which is outside this
   decision.
5. **Panics.** The binding's release profile is a copy of `LOSAT/Cargo.toml`'s with `panic = "unwind"`
   in place of `abort`, so that a panic reaches Python as an exception instead of ending the
   interpreter. An allocation failure still aborts, as on the command line.

## Consequences

- `python/tests` compares reports, diagnostics and table rows with the `losat` executable of the
  same commit. The table's numbers that a report prints cut to a few digits stay covered by
  `LOSAT/tests/library_api.rs`.
- Publishing (wheels, PyPI) is stage 4 of the plan and needs its own review.
