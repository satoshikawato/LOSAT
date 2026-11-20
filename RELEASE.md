# Release Procedure

This file is the short entry point for publishing LOSAT releases. The detailed
v0.1.0 gate lives in:

- [docs/archive/v0.1.0_release_readiness_plan.md](docs/archive/v0.1.0_release_readiness_plan.md)
- [docs/release/v0.1.0.md](docs/release/v0.1.0.md)
- [docs/archive/v0.1.0_scope.md](docs/archive/v0.1.0_scope.md)
- [BLASTN v0.1.0 certification](docs/release/blastn_v0.1.0_certification.md)
- [BLASTP v0.1.0 certification](docs/release/blastp_v0.1.0_certification.md)
- [TBLASTX v0.1.0 certification](docs/release/tblastx_v0.1.0_certification.md)

## CHANGELOG entries

Pull requests leave their entries in [`changelog.d/`](changelog.d/README.md), one file each, not in
`CHANGELOG.md`. Before a release names its version in `CHANGELOG.md`, run
`python3 scripts/changelog_fragments.py collect` (it moves the entries under `## [Unreleased]`,
newest first, and deletes the files) and commit both.

## v0.2.0

The [v0.2.0 release record](docs/release/v0.2.0.md) holds the scope, the platform status and the items still open. The [release test plan](docs/release/v0.2.0_certification_plan.md) lists the runs on the candidate commit and their order (section 7); v0.2.0 has its own contract (`docs/release/v0.2.0_rc_contract.json`), assembler (`LOSAT/tests/prepare_release_candidate_v020.py`) and workflow (`.github/workflows/release-readiness-v020.yml`). The v0.1.0 workflow and exact-SHA contract below are historical and must not be presented as a v0.2.0 artifact handoff.

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
- Use the exact-SHA RC contract in
  [`docs/release/v0.1.0_rc_contract.json`](docs/release/v0.1.0_rc_contract.json);
  do not copy version, target, or certification authority into another release
  script.

## Required Local Gate

```bash
cd LOSAT
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
cargo package --list
cargo package --locked \
  --config 'build.target-dir="/tmp/losat-cargo-package-target"'
cargo publish --dry-run --locked \
  --config 'build.target-dir="/tmp/losat-cargo-publish-target"'
```

## Certified Parity Authority

NCBI BLAST+ executables may be used only in this validation role.

The merged PR 5/PR 6 evidence satisfies this gate for an exact-SHA candidate
whose runtime, build, fixture, classifier, exception, and Product Decision
inputs remain unchanged. Do not rerun the commands below solely because release
documentation or artifact metadata changed. They remain the bounded
recertification entry points if the post-merge gate finds an invalidating
change; in that case, stop the ordinary RC path until the new evidence is
formally reviewed.

```bash
cd LOSAT
cargo build --release --locked
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --fresh-paired \
  --paired-output-dir /tmp/losat-blastn-v010-certification/paired-base \
  --losat-bin target/release/losat \
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

Run `.github/workflows/release-readiness.yml` manually on the branch whose head is
the exact committed candidate SHA (`--ref <branch> -f candidate_sha=<SHA>`; the
workflow checks out the dispatched commit and fails early when it is not the
candidate). Leave `run_integrated_certification` false unless the post-merge
gate identifies an output-affecting runtime/build/contract change after the
certified lineage.

The workflow must finish with one `LOSAT-v0.1.0-rc-<SHA>` handoff artifact. It
contains:

- native archives for Linux x64, Windows x64, macOS arm64, and macOS x64;
- the serial `wasm32-wasip1` command artifact;
- the Cargo source package produced by the package/publish dry runs;
- `SHA256SUMS`, per-artifact metadata, `RC-HANDOFF.json`, and
  `RC-HANDOFF.md`.

The workflow verifies exact source ancestry and CERT_TOOLCHAIN identity,
records candidate-specific binary and archive hashes, inspects target
architecture, extracts every binary archive, performs `--version` and
representative BLASTP smoke checks, and installs/smokes the source package in a
clean temporary location. Historical certification binary hashes remain
provenance identities, not cross-runner byte-equality gates; the metadata
records whether a candidate build happens to match them. Signing and
notarization remain explicitly unperformed.

Record the workflow run URL and generated `RC-HANDOFF.json` with the release
review. Do not edit an artifact or checksum after the workflow creates it.

## Publish

Only after the release note, artifacts, checksums, and parity evidence agree:

```bash
git tag -a v0.1.0 -m "LOSAT v0.1.0"
git push origin release/v0.1.0
git push origin v0.1.0
```

From v0.2.0 the release is also published on PyPI (`pylosat`, the Python package, under LOSAT's
version number) and on crates.io (`losat`). Both uploads are the maintainer's act, started by hand
after the GitHub release exists, on `main` while its head is the release (or on the release's tag),
in this order. A run publishes the commit it was started on, and an upload requires that commit to
be the source of the tag `v<version>`: the same tree, `.github/` aside, so that a fix of these
workflows after the tag needs no new tag (`.github/scripts/release_tag_source.sh`).

1. PyPI: workflow `Python wheels` (`.github/workflows/python-wheels.yml`) with `publish` = `pypi`
   (`testpypi` rehearses the same upload on TestPyPI from any commit). It builds the five wheels and
   the sdist, runs the frozen search through each, and uploads them from the environment `pypi` by
   Trusted Publishing: no token is stored. PyPI never takes a file name twice, so a wrong upload
   cannot be replaced, only yanked.
2. crates.io, last: workflow `Publish crate` (`.github/workflows/publish-crate.yml`). With `publish`
   off it packages the crate, holds it to the 10 MiB limit, compares it byte for byte with the GitHub
   release's `losat-<version>.crate` and runs `cargo publish --dry-run`; with `publish` on it then
   uploads with the token of the environment `crates-io`. A version on crates.io can be yanked but
   never deleted.

Both environments have the maintainer as required reviewer; the upload job waits for that approval.


Threading remediation builds use separate `serial-command`, `threaded-command`,
`serial-reactor`, and `threaded-reactor` target directories. The release assembler
checks that its input is a serial command before running or packaging it. The
frozen v0.1.0 artifact contract remains serial-only. Threaded/reactor integration
checks are additional regression evidence, not promotion into that contract.
The Node 18 certification record is historical; current threading gates pin
Rust 1.92.0 and Node 24.21.0 and retain their own artifact metadata.
