# Changelog

All notable release-facing changes and milestones are documented here.

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
