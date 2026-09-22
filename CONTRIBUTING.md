# Contributing to LOSAT

Thank you for your interest in contributing to LOSAT!

LOSAT is a standalone, pure-Rust reimplementation of NCBI BLAST+ local sequence alignment designed for direct pairwise sequence comparisons with **bit-perfect numerical and coordinate parity**. To maintain scientific integrity and reproducibility, all contributions must adhere to the engineering principles and verification workflow outlined below.

---

## Core Engineering Principles

### 1. NCBI BLAST+ Source as the Ground Truth
- The authoritative reference for all algorithmic behavior is the official NCBI BLAST C/C++ codebase (`ncbi-blast`).
- We do not approximate algorithm behavior or introduce ad-hoc heuristics. If a behavior, parameter, or pruning rule exists in NCBI BLAST+, its logic must be ported faithfully.
- Every code modification affecting algorithm behavior must include NCBI C/C++ reference comments with the file path and line numbers immediately above the Rust code:
  ```rust
  // NCBI reference: c++/src/algo/blast/core/blast_hits.c:993-1001
  ```

### 2. Standalone Pure-Rust Implementation
- LOSAT is a standalone Rust executable and library. It must never embed, link, or invoke external NCBI BLAST+ binaries or C/C++ libraries at runtime or during the build.
- NCBI BLAST+ is used strictly as an **external validation oracle** for unit testing, integration tests, and parity verification.
- If a feature or parameter is not yet implemented in Rust, LOSAT must fail fast with an explicit "unsupported option" error rather than delegating or falling back to an external tool.

### 3. Bit-Perfect Output Parity
- Alignment boundaries, coordinates, raw scores, bit scores, E-values, and tabular fields (`-outfmt 6` and `-outfmt 7`) must match NCBI BLAST+ exactly.
- *Approved Project Exception*: In TBLASTX local pairwise searches, LOSAT explicitly respects `--db-gencode` for translating the subject sequence across all non-standard genetic codes (NCBI BLAST+ defaults to standard code in local `-subject` mode). All other scoring, Karlin-Altschul statistics, and traceback mechanics strictly adhere to NCBI.

---

## Development Checks

Before submitting changes, ensure the codebase builds cleanly, passes all lints, and formats properly:

```bash
cd LOSAT
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
```

For packaging changes:

```bash
cd LOSAT
cargo package --list
cargo publish --dry-run --locked \
  --config 'build.target-dir="/tmp/losat-cargo-publish-target"'
```

---

## Parity Evidence & Testing

Any behavioral changes to search, scoring, or filtering algorithms require comparative verification against official NCBI BLAST+ (v2.17.0+):

1. **Run Comparison Suite**:
   ```bash
   cd LOSAT/tests
   ./run_comparison.sh
   ```
2. **Record Evidence**:
   - NCBI BLAST+ version used as oracle.
   - Exact LOSAT commit SHA.
   - Command lines, query/subject FASTA inputs, and parameter flags.
   - Verified output diffs covering coordinates, bit scores, and E-values.

---

## Repository Hygiene

- Do not commit generated comparison output, trace files, temporary test outputs, or Python cache files unless deliberately adding documented canonical fixtures.
- Keep the git working tree clean of build targets, database index files, or local scratch scripts.

---

## Pull Request Checklist

Before opening a PR, please verify:
- [ ] Changes align with the documented scope or explicitly document any experimental features.
- [ ] Relevant NCBI C/C++ source citations are included in code comments for all ported algorithms.
- [ ] No runtime, build, or fallback code invokes external NCBI binaries or libraries.
- [ ] Code builds without warnings under `cargo fmt` and `cargo clippy`.
- [ ] Parity comparison tests confirm bit-identical output with NCBI BLAST+ on affected profiles.
