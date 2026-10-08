# LOSAT comparison entry points

The repository contains several overlapping comparison scripts. Read each selected script before running it and confirm its inputs, outputs, NCBI binary, and cleanup behavior.

From the repository root:

```bash
cd LOSAT && cargo build --release
cd LOSAT && cargo test
cd LOSAT/tests && bash run_comparison.sh
cd LOSAT/tests && bash run_all_tests.sh
bash compare_tblastx_results.sh
bash tests/compare_self_tblastx.sh
bash tests/compare_long_sequences_debug.sh
```

Native and command-Wasm release builds:

```bash
cd LOSAT
cargo build --release
cargo build --release --target wasm32-wasip1 --no-default-features
cargo build --release --target wasm32-wasip1-threads --features wasm-threads
```

Use `rg` in the pinned NCBI tree (`satoshikawato/ncbi-blast` at commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`;
the local path is machine-specific, `$NCBI_SRC` in Claude Code) before relying on a remembered function or line number:

```bash
rg -n "<function-or-symbol>" "$NCBI_SRC/c++/src/algo/blast"
```

Prefer temporary output directories. Existing `tests/*_out` files may be old evidence or release-hygiene targets; do not overwrite them blindly.
