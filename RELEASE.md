# Release Procedure

This file is the entry point for publishing LOSAT releases. A release is cut from the head of `dev`,
tested on that one commit, promoted to `main`, and published from the branch `release/<major>.<minor>`.
The procedure LOSAT v0.1.0 was published with is kept in
[docs/archive/v0.1.0_release_procedure.md](docs/archive/v0.1.0_release_procedure.md).

## CHANGELOG entries

Pull requests leave their entries in [`changelog.d/`](changelog.d/README.md), one file each, not in
`CHANGELOG.md`. Before a release names its version in `CHANGELOG.md`, run
`python3 scripts/changelog_fragments.py collect` (it moves the entries under `## [Unreleased]`,
newest first, and deletes the files) and commit both.

## v0.2.0

The [v0.2.0 release record](docs/release/v0.2.0.md) holds the scope, the platform status and the
results on the candidate. The [release test plan](docs/release/v0.2.0_certification_plan.md) lists the
runs on the candidate commit (section 2) and their order (section 7). v0.2.0 has its own contract
(`docs/release/v0.2.0_rc_contract.json`), assembler (`LOSAT/tests/prepare_release_candidate_v020.py`)
and workflow (`.github/workflows/release-readiness-v020.yml`).

1. Candidate `C`: a head of `dev` on which the full `platform-checks` run (`set=full`) and the `nightly`
   run are green.
2. Runs on `C` (plan, section 7): the GitHub workflows on `cert/<first 9 hex of C>`, sets A, D and E on
   a Linux x64 machine with NCBI BLAST+ 2.17.0, and set J (memory and speed) on the workstation.
3. Promotion: the `promote/<date>` pull request from `dev` to `main`, merged as a merge commit `M`
   whose tree is `C`'s (`git diff --quiet C M`); `release/0.2` is created at `M`.
4. Release commit `R` on `release/0.2`: the filled contract, the results in the release record,
   `CHANGELOG.md` with the entries collected, and the evidence under `docs/evidence/release_v020/`;
   no runtime input differs from `C`. Then dispatch `release-readiness-v020.yml --ref release/0.2
   -f candidate_sha=<R>` and record its run.
5. GO (the Owner): tag `v0.2.0` on `R`, the GitHub Release with the eight artifacts, `SHA256SUMS` and
   `LOSAT-0.2.0.cdx.json`, then Bioconda ([packaging/bioconda/README.md](packaging/bioconda/README.md)).

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
