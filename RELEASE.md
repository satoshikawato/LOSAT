# Release Procedure

This file is the short entry point for publishing LOSAT releases. The detailed
v0.1.0 gate lives in:

- [docs/v0.1.0_release_readiness_plan.md](docs/v0.1.0_release_readiness_plan.md)
- [docs/release/v0.1.0.md](docs/release/v0.1.0.md)
- [docs/v0.1.0_scope.md](docs/v0.1.0_scope.md)

## v0.1.0 Release Candidate Gate

Do not create a release candidate while any P0 readiness item is unresolved.

Before tagging:

- Confirm the release scope in README, CHANGELOG, CLI help, and release notes.
- Confirm NCBI BLAST+ is used only as a comparison oracle, never as LOSAT
  runtime, build, feature, fallback, or unsupported-feature implementation.
- Run the native quality gate from a clean release branch.
- Run current NCBI parity fixtures for every supported scope.
- Record NCBI BLAST+ version, LOSAT commit, command lines, inputs, outputs,
  checksums, and diff summaries in the release note.
- Build native/Wasm artifacts only for targets documented in the release note.
- Record artifact checksums and smoke-test extracted artifacts.
- Keep generated comparison output and debug scratch out of release source.

## Required Local Gate

```bash
cd LOSAT
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
cargo package --list
cargo publish --dry-run --locked \
  --config 'build.target-dir="/tmp/losat-cargo-publish-target"'
```

## Required Parity Gate

NCBI BLAST+ executables may be used only in this validation role.

```bash
cd LOSAT/tests
RUN_LOSAT_WASM=0 RUN_LOSAT_WASM_THREADED=0 bash run_comparison.sh
python3 compare_blastn_parity.py
bash run_blastp_comparison.sh

cd ../..
bash tests/compare_tblastx_native_ncbi_parity.sh
bash tests/compare_self_tblastx.sh LOSAT/tests/fasta/NZ_CP006932.fasta 4
bash tests/compare_long_sequences_debug.sh \
  LOSAT/tests/fasta/AP027131.fasta \
  LOSAT/tests/fasta/AP027133.fasta \
  4
```

## Artifact Gate

Use `.github/workflows/release-readiness.yml` on the release branch to assemble
release-candidate artifacts and collect checksums. Copy the resulting metadata
into [docs/release/v0.1.0.md](docs/release/v0.1.0.md) before drafting the
GitHub Release.

## Publish

Only after the release note, artifacts, checksums, and parity evidence agree:

```bash
git tag -a v0.1.0 -m "LOSAT v0.1.0"
git push origin release/v0.1.0
git push origin v0.1.0
```

If publishing to crates.io:

```bash
cd LOSAT
cargo publish
```
