# Changelog

All notable release-facing changes are documented here.

## v0.1.0 - Unreleased

Initial release candidate for LOSAT as a Rust implementation of
BLAST-compatible local sequence alignment.

### Added

- Native CLI entry points for `blastn`, `tblastx`, and experimental `blastp`.
- TBLASTX local `-query`/`-subject` release candidate scope.
- BLASTN `megablast` and `blastn` task fixtures for parity validation.
- Experimental command-Wasm build support for `wasm32-wasip1`.
- Release readiness, scope, and release-note draft documents.
- Root release procedure checklist for v0.1.0 release candidates.
- Contributor and security policy documents for release-candidate handling.

### Changed

- User-facing documentation now states that NCBI BLAST+ is a validation oracle
  only, not a runtime dependency or fallback path.
- BLASTP documentation now reflects the implemented but experimental CLI
  subcommand instead of describing BLASTP as merely planned.
- Wasm documentation now describes WASI command builds rather than
  `wasm-bindgen` browser API stability.

### Known Limitations

- BLASTN remains experimental until current comparison fixtures prove
  bit-perfect parity for the exact release scope.
- BLASTP remains experimental in v0.1.0.
- Rust library and embeddable Web/Wasm APIs are internal only.
- Large database search workflows are outside the v0.1.0 release scope.

### Approved Deviation

- TBLASTX local `-subject` searches intentionally honor non-default
  `--db-gencode` for subject translation/search/reporting. This is the only
  approved v0.1.0 deviation from NCBI BLAST+ local `-subject` behavior.

### Verification

Final release verification is pending. The release tag must record NCBI BLAST+
version, LOSAT commit, commands, fixtures, checksums, and diff summaries in
`docs/release/v0.1.0.md`.
