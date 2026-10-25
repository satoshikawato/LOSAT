- The v0.2.0 release-candidate assembler `LOSAT/tests/prepare_release_candidate_v020.py` (with its unit
  tests) reads `docs/release/v0.2.0_rc_contract.json`: it builds and verifies the eight artifacts
  (five native targets, serial and threaded Wasm, the source crate), runs one smoke search per program
  (blastn, blastp, blastx, tblastn, tblastx) against NCBI-certified output hashes, and checks that the
  release-candidate commit descends from the certified candidate with no runtime-input change. The
  v0.1.0 contract and tooling are unchanged. No LOSAT behaviour changes.
