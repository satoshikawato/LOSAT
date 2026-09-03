# Changelog

All notable release-facing changes are documented here.

## v0.1.0 - Unreleased

Initial release candidate for LOSAT as a Rust implementation of
NCBI BLAST local-sequence-alignment behavior for certified profiles.

### Added

- Native CLI entry points for `blastn`, `blastp`, and `tblastx`.
- Supported BLASTN `blastn`/`megablast` local profile backed by the durable
  14-case certification record: 13 exact source-defined cases and one narrow
  Version 1.2 source-underdetermined equal-HSP contract.
- Supported BLASTP / LOSATP gbdraw P1-P3 local standard-outfmt-6 profiles
  backed by nine fresh `EXACT_TEXT`, repeatable certification cases.
- Supported TBLASTX / TLOSATX gbdraw P1-P2 local standard-outfmt-6 profiles
  backed by 14 exact source-defined cases and six passing approved-deviation
  contracts.
- Serial `wasm32-wasip1` command-build support with the program-specific
  integrated evidence covering all 41 directly applicable manifest rows.
- Application-level pure-Rust runtime ownership with zero project-authored
  production algorithm delegation findings.
- Cross-platform native certification for the complete 43-contract LOSAT
  matrix on Windows x64, macOS arm64, and macOS x64 against the frozen Linux
  canonical output.
- Exact-SHA release-candidate assembly for four native archives, serial
  `wasm32-wasip1`, and the Cargo source package, including checksums,
  provenance, architecture checks, clean extraction/install, and smoke tests.
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
- Release readiness now consumes the PR 5/PR 6 certification lineage and
  reruns the expensive integrated campaign only when the post-merge gate
  identifies an invalidating change.

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
`docs/release/tblastx_v0.1.0_certification.md`. Integrated Linux/serial-Wasm
authority is recorded in
`docs/release/pure_rust_runtime_v0.1.0_certification.md`; cross-platform native
authority is PR 6 run `33625511701`. The final tag still requires a successful
exact-SHA release-readiness run and separate publication authorization.
