# AGENTS.md — python/

Instructions for coding agents working under `python/`. The root [`AGENTS.md`](../AGENTS.md) still
governs everything under `LOSAT/`, including engine changes made for the binding.

## Authority

- [`PD-LOSAT-PYTHON-BINDING-BOUNDARY`](../docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md)
  defines what code here may and may not do.
- [`docs/losat_library_surface_plan.md`](../docs/losat_library_surface_plan.md) holds the plan of the
  library and of this binding (§4.3, §12).

## Rules

1. Reach the engine through `LOSAT::api` only. No other `LOSAT::` path, even one that compiles.
2. Do not compute or format BLAST-defined values (scores, E-values, identities, coverages, NCBI
   number formats, alignment text). Add a value that is missing to `LOSAT::api`, under the root
   rules.
3. Reports stay byte for byte what the engine wrote.
4. Give the engine only an argv and FASTA bytes. Never split one search into independent searches
   and join the results.
5. No Rayon `build_global`, no change of the environment, the current directory or the panic hook.
6. Keep what NCBI does not have here. Arrow, pickle and Python's shapes never enter `LOSAT/`.
7. Names: `Hsp` has the field names of `HspRecord`; `HspTable`'s columns have `-outfmt`'s names.
8. A new `Hsp` field goes into `of_records`, the getters, `_pylosat.pyi` and `tests/test_results.py`
   together.
9. `pylosat/_options.pyi` is generated: run `tools/generate_option_stubs.py`.
10. Keep `[profile.release]` in `Cargo.toml` a copy of `LOSAT/Cargo.toml`'s, `panic` aside, and
    `.cargo/config.toml` a copy of the x86-64 rustflags of `LOSAT/.cargo/config.toml`.

## Commands

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r tests/requirements.txt
maturin develop                                # the extension, into the venv
cargo build --manifest-path ../LOSAT/Cargo.toml --bin losat
LOSAT_BIN=../LOSAT/target/debug/losat python -m pytest tests
python tools/generate_option_stubs.py
python -m mypy --strict pylosat
cargo fmt --check && cargo clippy --locked --all-targets -- -D warnings
```

`maturin` needs `RUSTUP_TOOLCHAIN` only where the repository's pinned toolchain is missing.
