# Release Procedure

This file is the short entry point for publishing LOSAT releases. The detailed
v0.1.0 gate lives in:

- [docs/v0.1.0_release_readiness_plan.md](docs/v0.1.0_release_readiness_plan.md)
- [docs/release/v0.1.0.md](docs/release/v0.1.0.md)
- [docs/v0.1.0_scope.md](docs/v0.1.0_scope.md)
- [BLASTN v0.1.0 certification](docs/release/blastn_v0.1.0_certification.md)
- [BLASTP v0.1.0 certification](docs/release/blastp_v0.1.0_certification.md)
- [TBLASTX v0.1.0 certification](docs/release/tblastx_v0.1.0_certification.md)

## v0.1.0 Release Candidate Gate

Do not create a release candidate while any P0 readiness item is unresolved.

Before tagging:

- Confirm the release scope in README, CHANGELOG, CLI help, and release notes.
- Confirm NCBI BLAST+ is used only as a comparison oracle, never as LOSAT
  runtime, build, feature, fallback, or unsupported-feature implementation.
- Run the native quality gate from a clean release branch.
- Enforce the committed certification gates for every supported program
  profile; use older broad comparison scripts only as diagnostics.
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
cd LOSAT
cargo build --release --locked
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --fresh-paired \
  --paired-output-dir /tmp/losat-blastn-v010-certification/paired-base \
  --losat-bin target/release/LOSAT \
  --ncbi-bin /home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn
python3 tests/certify_blastn_v010.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --paired-output-dir /tmp/losat-blastn-v010-certification/paired-base \
  --exceptions tests/blastn_v010_source_exceptions.tsv

cd ..
python LOSAT/tests/audit_blastp_v010.py \
  --output-dir /tmp/losat-blastp-v010-audit/final-native
python LOSAT/tests/audit_tblastx_v010.py \
  --output-dir /tmp/losat-tblastx-v010-certification/final-native
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
