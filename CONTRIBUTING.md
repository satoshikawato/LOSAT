# Contributing to LOSAT

Thank you for your interest in LOSAT. LOSAT is a maintainer-led open-source project providing a lightweight, pure-Rust reimplementation of NCBI BLAST+ local sequence alignment with bit-perfect parity.

Use, inspection, reproducible bug reports, minimal edge-case test sequences, research workflows, performance benchmarks, and focused documentation corrections are welcome.

The maintainer retains final authority over product behavior, architecture, algorithmic scope, compatibility policy, and releases. Bug reports with minimal reproducible FASTA inputs and small documentation corrections may be opened directly. Anything larger must begin with an issue and wait for explicit maintainer confirmation on the intended behavior and whether a pull request is the appropriate next step.

External pull requests are optional proposals, not the project's default development path or an entitlement to merge. An accepted idea may be revised, split, independently implemented, or adopted without merging the submitted patch. Contributing does not imply roadmap priority, repository access, a merge, or maintainer status.

---

## Governance and Decision Rights

User reports and scientific feedback inform the project, but LOSAT is not governed by contributor voting or consensus. The maintainer decides:

- whether a feature or use case belongs within LOSAT's scope;
- authoritative search, scoring, and output behavior;
- public CLI options, argument parsing, and output formatting contracts;
- internal architecture, data structures, and dependency changes;
- performance optimization strategies;
- release contents, milestones, and timing; and
- whether a proposed change is sufficiently maintainable, verified, and aligned with project goals to merge.

A technically functioning pull request may still be declined if it would introduce:

- any deviation from NCBI BLAST+ algorithm ground truth;
- heuristic approximations or non-NCBI candidate pruning;
- an external runtime dependency or fallback to NCBI binaries;
- unnecessary complexity or excessive maintenance overhead;
- parallel implementation paths or uncertified options; or
- scope drift (e.g., pre-formatted database indexing, unsupported programs).

---

## Core Technical Requirements

All contributions that touch alignment or search code must strictly satisfy these non-negotiable principles:

### 1. NCBI BLAST+ Source as the Ground Truth
- The authoritative reference for all algorithmic behavior is the official NCBI BLAST C/C++ source code (`ncbi-blast`).
- We do not speculate, approximate, or introduce unverified heuristics. Every ported algorithm must faithfully follow NCBI's exact logic, scoring, and statistics.
- Code modifications affecting search or scoring behavior must include the corresponding NCBI C/C++ reference comments (file path and line numbers) immediately above the Rust code:
  ```rust
  // NCBI reference: c++/src/algo/blast/core/blast_hits.c:993-1001
  ```

### 2. Standalone Pure-Rust Implementation
- LOSAT is a pure-Rust application. It must never embed, link, or invoke external NCBI BLAST+ binaries or C/C++ shared libraries at runtime or build time.
- Official NCBI BLAST+ executables (v2.17.0+) are used strictly as an **external validation oracle** during parity testing and benchmarks.
- Unsupported options must fail fast with an explicit error message rather than degrading silently or delegating to external tools.

### 3. Bit-Perfect Output Parity
- Alignment boundaries, coordinates, raw scores, bit scores, E-values, and tabular output (`-outfmt 6` and `-outfmt 7`) must match NCBI BLAST+ exactly.
- *Approved Project Exception*: In TBLASTX local pairwise searches, LOSAT explicitly respects `--db-gencode` for translating the subject sequence across all non-standard genetic codes (NCBI BLAST+ defaults to standard code in local `-subject` mode). All other scoring, statistical calculations, and alignment logic strictly adhere to NCBI.

---

## Before You Start

- **Search Existing Discussions**: Check [open and closed issues](https://github.com/satoshikawato/LOSAT/issues) and [pull requests](https://github.com/satoshikawato/LOSAT/pulls) before submitting.
- **Discuss Non-Trivial Changes First**: Open an issue to discuss proposed features, algorithmic changes, or refactoring before writing code. Do not invest time in large implementations without maintainer agreement.
- **Unsolicited Broad Changes**: Large uncoordinated pull requests, broad refactors, non-NCBI feature additions, or AI-generated bulk rewrites may be closed without line-by-line review.

---

## Reporting a Bug or Parity Discrepancy

When reporting an issue or output discrepancy against NCBI BLAST+, please provide:

1. **Exact command lines** used for both `LOSAT` and `blastn`/`blastp`/`tblastx`.
2. **Minimal reproducible FASTA inputs** (synthetic or non-sensitive sequences).
3. **Exact diff of outputs**, highlighting discrepancies in coordinates, scores, or E-values.
4. **Environment details**: OS, architecture, LOSAT version/commit, and NCBI BLAST+ version.

Do not upload confidential, restricted, or proprietary genomic data. Minimal synthetic excerpts that reproduce the failure are preferred.

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
- [ ] An issue was opened and approved by the maintainer for non-trivial changes.
- [ ] Relevant NCBI C/C++ source citations are included in code comments for all ported algorithms.
- [ ] No runtime, build, or fallback code invokes external NCBI binaries or libraries.
- [ ] Code builds without warnings under `cargo fmt` and `cargo clippy`.
- [ ] Parity comparison tests confirm bit-identical output with NCBI BLAST+ on affected profiles.
