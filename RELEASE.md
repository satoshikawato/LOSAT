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
- Run the native quality gate from a clean work branch at the exact candidate SHA.
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

Run `.github/workflows/release-readiness.yml` manually with the exact committed
candidate SHA. Leave `run_integrated_certification` false unless the post-merge
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

Only after all gates pass, present the exact candidate SHA, handoff directory,
release notes, checksums and workflow URLs for publication approval. Do not
infer publication permission from a request to prepare or implement a release.
Never rebuild, rename, recompress or overwrite a published asset.

After explicit authorization, confirm the tag and Release do not already exist.
Use the exact reviewed commit (not an implicit HEAD or nonexistent release branch):

```bash
git tag -a v0.1.0 <reviewed-candidate-SHA> -m "LOSAT v0.1.0"
git push origin refs/tags/v0.1.0
gh release create v0.1.0 --repo satoshikawato/LOSAT --verify-tag --draft \
  --title "LOSAT v0.1.0" --notes-file <reviewed-public-release-notes>
gh release upload v0.1.0 --repo satoshikawato/LOSAT <verified-handoff-files>
```

Upload all six original artifacts, metadata/checksum sidecars, `SHA256SUMS`,
`RC-HANDOFF.json` and `RC-HANDOFF.md`. Verify uploaded names/sizes/hashes and
tag/candidate identity, then publish the reviewed draft. Download each fixed
`https://github.com/satoshikawato/LOSAT/releases/download/v0.1.0/<asset>` URL
without credentials and verify bytes and native execution on all four targets.
crates.io, Bioconda, signing and notarization are outside this distribution.

Linux release preparation assumes glibc 2.34 or newer; re-inspect the final ELF
with `readelf --version-info` and execute on the declared minimum environment.
Inspect macOS load-command deployment targets and Windows DLL imports on their
runners. Runner OS labels alone do not establish minimum supported OS versions.
Do not publish broader support claims than the measured artifact evidence.


Threading remediation builds use separate `serial-command`, `threaded-command`,
`serial-reactor`, and `threaded-reactor` target directories. The release assembler
checks that its input is a serial command before running or packaging it. The
frozen v0.1.0 artifact contract remains serial-only. Threaded/reactor integration
checks are additional regression evidence, not promotion into that contract.
The Node 18 certification record is historical; current threading gates pin
Rust 1.92.0 and Node 24.21.0 and retain their own artifact metadata.
