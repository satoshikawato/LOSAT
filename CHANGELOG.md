# Changelog

All notable release-facing changes and milestones are documented here.

## [Unreleased]

- Fixed when `blastp` ends the composition-based redo of a query early. NCBI BLAST+ ends the remaining matches of a query only when the hit lists of all queries of its query batch (10,000 residues of queries) are full. LOSAT looked at the query's own hit list alone with two or more threads and two or more queries that have hits (a query whose hit list was full got the result of searching it alone: two queries against Swiss-Prot differed in 42 of 624 rows), and at the hit lists of all queries of the input, not of the batch, with one thread or one query that has hits. With `-max_target_seqs` 500 either case needs more than 500 significant subjects for one query. The result now follows NCBI's batches at every thread count.
- The exact performance paths of the 2026-10 performance work (46 `LOSAT_X_*` switches of the "strict set") are on by default in every program, natively and in WebAssembly. Each of them gives the reference result by construction, so reports do not change. `LOSAT_X_DEFAULT=0` runs the previous paths, and `LOSAT_X_<NAME>=0` turns one switch off ([`docs/agents/diagnostics.md`](docs/agents/diagnostics.md), "Performance switches"). The Web reactors receive no environment and always run the default.
- Costs of the default: megablast keeps one byte per base of a subject while it is searched, and tblastn keeps up to 256 MiB of subject translations for a search. TBLASTX builds only the lookup tables its seed stage scans and sorts and links the HSP list of a subject without a second copy; on a 0.66 Mb by 0.61 Mb pair its peak memory is 224 MiB at 8 threads and 168 MiB at 1 thread, below the previous paths (281 and 241 MiB). On targets other than x86_64 and WebAssembly no vector X-drop kernel exists, and the scalar code runs there.

## [v0.2.0] - TLOSAN release candidate (unpublished)

- Added local `tblastn` search with 27 NCBI genetic-code IDs and outfmt 0/6/7. Code 32 and local selected subject translation follow the TBLASTN-only `PD-TLOSAN-LOCAL-GENCODE-32` decision.
- Certified the [declared Stage G fixtures](docs/evidence/tlosan_stage_g/STAGE_G_CERTIFICATION_20260926.md) on Linux x64 native, serial command-WASI, and threaded command-WASI. This does not certify arbitrary inputs/options or other platforms.
- Database search, `tblastn-fast`, PSI-TBLASTN, composition modes 1/3, and other unported options remain explicitly unsupported.
- The v0.2.0 [release-readiness decision](docs/release/v0.2.0.md) governs packaging and handoff. No tag or distribution is implied by this entry.

## [v0.1.0] - Initial Release Candidate

Initial release candidate for LOSAT as a standalone, pure-Rust reimplementation of NCBI BLAST+ local sequence alignment behavior, designed for direct pairwise FASTA comparisons (`-query` vs `-subject`) with bit-identical output parity and WebAssembly portability.

### Added

- **Supported Programs & Search Tasks**:
  - **`blastn`**: Supports pairwise nucleotide alignment with `megablast` (default, word size 28, affine/linear penalties) and traditional `blastn` (word size 11, match/mismatch 2/-3, gaps 5/2).
  - **`blastp`**: Supports pairwise protein alignment using BLOSUM62 matrix and affine gap penalties (open 11, extend 1).
  - **`tblastx`**: Supports 6-frame translated pairwise nucleotide alignment with full genetic code customization.
- **Bit-Perfect NCBI BLAST+ Parity**: Produces identical alignment boundaries, bit scores, and E-values matching official NCBI BLAST+ (v2.17.0+) across certified pairwise profiles.
- **Standard Output Formats**: Supports canonical tabular output (`-outfmt 6`) and commented tabular (`-outfmt 7`) reporting the standard 12 BLAST fields (`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`), plus pairwise text alignment (`-outfmt 0`) for `blastp`.
- **Pure-Rust Standalone Runtime**: Zero external C/C++ dependencies or shared libraries. Does not link, embed, wrap, or invoke NCBI BLAST+ binaries at runtime.
- **Direct Pairwise Alignment**: Directly aligns sequences from FASTA files via `-query` and `-subject` without requiring database formatting (`makeblastdb`).
- **WebAssembly & WASI Support**: Compiles to standalone WASI command-line binaries (`wasm32-wasip1`) and multithreaded WASI (`wasm32-wasip1-threads`) for sandboxed CLI and browser-based bioinformatics applications (such as [gbdraw](https://github.com/satoshikawato/gbdraw)).
- **Cross-Platform Native Binaries**: Pre-compiled and verified binaries for Linux (x86_64, aarch64), macOS (Apple Silicon arm64 & Intel x86_64), and Windows (x86_64).
- **Bioconda Distribution**: Officially available on Bioconda as `losat` (`conda install -c bioconda losat`) across Linux and macOS architectures.

### Changed

- User-facing documentation clarifies that NCBI BLAST+ is used strictly as an independent validation oracle during testing, not as a runtime dependency or fallback path.
- Parameter validation fails fast with explicit error messages when encountering unsupported options rather than silently degrading or falling back.
- Wasm documentation clearly distinguishes standalone WASI command builds from experimental browser library APIs.

### Approved Parity Enhancement

- **TBLASTX Subject Genetic Code (`--db-gencode`)**: In local pairwise `-subject` searches, LOSAT explicitly honors non-default `--db-gencode` for translating subject sequences (supporting genetic codes 1–33). This resolves a recognized NCBI BLAST+ limitation where local `-subject` mode defaulted to genetic code 1 regardless of `--db-gencode`. All scoring, statistical calculations, and alignment logic remain strictly identical to NCBI.

### Scope & Known Limitations

- **Database Search**: Pre-formatted BLAST databases (e.g., `.nin`, `.nhr` created by `makeblastdb`) and remote NCBI queries (`-remote`) are outside the v0.1.0 scope; comparisons are strictly pairwise FASTA (`-query` vs `-subject`).
- **Unimplemented Programs**: `blastx` and `tblastn` are not implemented in this release.
- **Discontinuous MegaBLAST**: `dc-megablast` is not yet supported.
- **Library API**: Rust crate library interfaces and embeddable Web/Wasm APIs are currently internal-only; stability is guaranteed at the CLI level.

### Verification Authority

Governed by program certification records in `docs/release/` (`blastn_v0.1.0_certification.md`, `blastp_v0.1.0_certification.md`, `tblastx_v0.1.0_certification.md`) and the pure-Rust runtime authority in `docs/release/pure_rust_runtime_v0.1.0_certification.md`.
