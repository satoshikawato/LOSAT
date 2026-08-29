# Changelog

All notable release-facing changes are documented here.

## v0.1.0 - Unreleased

Initial release candidate for LOSAT as a Rust implementation of
BLAST-compatible local sequence alignment.

### Added

- Native CLI entry points for `blastn`, `tblastx`, and experimental `blastp`.
- TBLASTX local `-query`/`-subject` release candidate scope.
- BLASTN `megablast` and `blastn` task fixtures for parity validation.
- A machine-checkable BLASTN v0.1.0 certification gate, one-row
  source-exception specification, and durable native certification record.
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

- BLASTN support is limited to the certified v0.1.0 local query/subject profile
  for `megablast` and `blastn`. `dc-megablast`, database search, and threaded
  Wasm BLASTN certification remain pending.
- BLASTP remains experimental in v0.1.0.
- Rust library and embeddable Web/Wasm APIs are internal only.
- Large database search workflows are outside the v0.1.0 release scope.

### Approved Deviation

- TBLASTX local `-subject` searches intentionally honor non-default
  `--db-gencode` for subject translation/search/reporting. This is the only
  approved v0.1.0 deviation from NCBI BLAST+ local `-subject` behavior.

### Verification

Current local release-gate evidence is recorded in `docs/release/v0.1.0.md`.
The supported TBLASTX candidate scope has current native NCBI parity evidence.
BLASTN native certification has 13 exact source-defined cases plus one narrowly
source-compatible equal-HSP case. The raw comparator continues to report the
real `VALUE_DIFF`; the separate certification gate accepts it only under the
frozen five-key invariant boundary. The final tag still requires the remaining
release-wide gates and workflow artifact checksums.
