# Changelog

All notable release-facing changes are documented here.

## v0.1.0 - Unreleased

Initial release candidate for LOSAT as a Rust implementation of
NCBI BLAST local-sequence-alignment behavior for certified profiles.

### Added

- Native CLI entry points for `blastn`, `blastp`, and `tblastx`.
- Supported BLASTN `blastn`/`megablast` local profile backed by the durable
  14-case certification record: 13 exact source-defined cases and one narrow
  Version 1.1 source-underdetermined equal-HSP contract.
- Supported BLASTP / LOSATP gbdraw P1-P3 local standard-outfmt-6 profiles
  backed by nine fresh `EXACT_TEXT`, repeatable certification cases.
- Supported TBLASTX / TLOSATX gbdraw P1-P2 local standard-outfmt-6 profiles
  backed by 14 exact source-defined cases and six passing approved-deviation
  contracts.
- Serial `wasm32-wasip1` command-build support with the program-specific
  certified or targeted evidence recorded by the three durable certifications.
- Release readiness, scope, and release-note draft documents.
- Root release procedure checklist for v0.1.0 release candidates.
- Contributor and security policy documents for release-candidate handling.

### Changed

- User-facing documentation now states that NCBI BLAST+ is a validation oracle
  only, not a runtime dependency or fallback path.
- User-facing v0.1.0 status now reflects the certified BLASTN, BLASTP, and
  TBLASTX local profiles without broadening them to generic BLAST support.
- Wasm documentation now describes WASI command builds rather than
  `wasm-bindgen` browser API stability.

### Known Limitations

- BLASTN support is limited to the certified v0.1.0 local query/subject profile
  for `megablast` and `blastn`. `dc-megablast`, database search, and threaded
  Wasm BLASTN certification remain pending.
- BLASTP support is limited to the certified gbdraw P1-P3 local standard-outfmt-6
  profiles; alternate tasks/options, other formats, database/remote search, and
  threaded Wasm are excluded.
- TBLASTX support is limited to the certified gbdraw P1-P2 local
  standard-outfmt-6 profiles with one thread per job; unexercised options,
  database/remote search, and threaded Wasm are excluded.
- Serial Wasm evidence is program-specific and does not establish general Wasm
  support. `wasm32-wasip1-threads` remains experimental.
- Rust library and embeddable Web/Wasm APIs are internal only.
- Large database search workflows are outside the v0.1.0 release scope.

### Approved Deviation

- TBLASTX local `-subject` searches intentionally honor non-default
  `--db-gencode` for subject translation/search/reporting. This is the only
  approved v0.1.0 deviation from NCBI BLAST+ local `-subject` behavior.

### Verification

Current program-profile evidence is governed by
`docs/release/blastn_v0.1.0_certification.md`,
`docs/release/blastp_v0.1.0_certification.md`, and
`docs/release/tblastx_v0.1.0_certification.md`. The final tag still requires
final v0.1.0 release-readiness and artifact certification, including provenance,
workflow artifact checksums, and smoke tests.
