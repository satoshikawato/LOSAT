# LOSAT: LOcal Sequence Alignment Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.92%2B-orange.svg)](https://www.rust-lang.org)
[![WebAssembly](https://img.shields.io/badge/Wasm-WASI%20%7C%20Web-purple.svg)](#webassembly-wasm-integration)

**LOSAT** is a lightweight, pure-Rust reimplementation of the NCBI BLAST sequence alignment algorithm designed specifically for pairwise sequence-to-sequence comparisons (`-query` vs `-subject`).

It delivers bit-identical alignment scores, E-values, and coordinates matching NCBI BLAST+ without requiring external C/C++ libraries, BLAST+ installations, or pre-formatted database indices. Built for high portability, LOSAT runs natively on modern operating systems and compiles directly to WebAssembly for client-side, zero-install genomic analyses in web browsers and sandboxed runtimes.

---

## Key Highlights

- **Bit-Perfect NCBI BLAST+ Parity**: Produces identical alignment coordinates, E-values, bit scores, and tabular records matching official NCBI BLAST+ (v2.17.0+) on certified profiles.
- **Pure Rust, Zero Dependencies**: Standalone single executable. Does not wrap, link, or shell out to external NCBI binaries or libraries.
- **Direct Pairwise Alignment**: Compares FASTA files directly via `-query` and `-subject` without running `makeblastdb`.
- **WebAssembly Ready**: Compiles to WASI and web reactors with multithreading support, powering in-browser bioinformatics visualization tools like [gbdraw](https://github.com/satoshikawato/gbdraw).
- **Corrected TBLASTX Genetic Codes**: Fully respects `--db-gencode` for translated subject sequences in pairwise searches (resolving NCBI BLAST+'s default fallback to standard code).
- **High Performance**: Multithreaded execution via Rayon natively and shared-memory threading in WebAssembly (`wasm32-wasip1-threads`).

---

## Supported Programs & Tasks

| Program | Supported Tasks | Output Formats (`-outfmt`) | Description |
|:---|:---|:---|:---|
| **`blastn`** | `megablast` (default), `blastn` | `6` (tabular), `7` (commented tabular) | Nucleotide vs. nucleotide alignment |
| **`blastp`** | `blastp` | `0` (pairwise), `6`, `7` | Protein vs. protein alignment |
| **`tblastx`** | `tblastx` | `6` | Translated 6-frame nucleotide vs. nucleotide alignment |

> **Note**: Standard tabular output (`-outfmt 6`) generates the 12 canonical BLAST fields:  
> `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`.

---

## Installation

### Pre-built Binaries
Download pre-compiled binaries for Linux (x86_64), macOS (Apple Silicon & Intel), and Windows from the [Releases](https://github.com/satoshikawato/LOSAT/releases) page.

### Build from Source
Requires the [Rust toolchain](https://rustup.rs/) (edition 2021, Rust 1.92+ recommended):

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/LOSAT
cargo build --release
```
The compiled binary will be located at `target/release/LOSAT`.

---

## Quick Start

LOSAT adopts the standard NCBI single-dash command-line syntax:

### 1. Nucleotide Alignment (BLASTN / MegaBLAST)
```bash
# Fast pairwise nucleotide search with MegaBLAST
LOSAT blastn \
  -query query.fna \
  -subject subject.fna \
  -task megablast \
  -outfmt 6
```

### 2. Protein Alignment (BLASTP)
```bash
# Protein alignment with custom E-value cutoff and multithreading
LOSAT blastp \
  -query query.faa \
  -subject subject.faa \
  -evalue 1e-5 \
  -num_threads 4 \
  -outfmt 6
```

### 3. Translated Alignment (TBLASTX)
```bash
# 6-frame translated search using bacterial genetic code (Table 11)
LOSAT tblastx \
  -query contigA.fna \
  -subject contigB.fna \
  -query_gencode 11 \
  -db_gencode 11 \
  -outfmt 6
```

---

## Common Command-Line Options

| Flag | Type | Default | Description |
|:---|:---|:---|:---|
| `-query <file>` | File path | *(Required)* | Input query sequence file (FASTA) |
| `-subject <file>` | File path | *(Required)* | Input subject sequence file (FASTA) |
| `-task <string>` | String | Program default | Task: `megablast` or `blastn` (for `blastn`); `blastp` (for `blastp`) |
| `-evalue <real>` | Float | `10.0` | Expectation value (E-value) threshold |
| `-outfmt <int>` | Integer | `6` | Output format (`6`=tabular, `7`=commented tabular, `0`=pairwise) |
| `-num_threads <int>` | Integer | `1` | Number of threads to use |
| `-max_target_seqs <int>`| Integer | `500` | Maximum number of aligned target sequences to keep |
| `-max_hsps <int>` | Integer | Unlimited | Maximum number of HSPs per subject sequence |
| `-dust <args>` | String | `20 64 1` | DUST low-complexity filter for BLASTN (`yes`, `no`, or `"level window linker"`) |
| `-seg <args>` | String | `no` (blastp)<br>`12 2.2 2.5` (tblastx) | SEG low-complexity filter for BLASTP/TBLASTX (`yes`, `no`, or `"window locut hicut"`) |
| `-query_gencode <int>` | Integer | `1` | Genetic code for query translation (TBLASTX) |
| `-db_gencode <int>` | Integer | `1` | Genetic code for subject translation (TBLASTX) |

---

## Accuracy and Benchmarks

### Alignment Parity
LOSAT is continuously audited against NCBI BLAST+ 2.17.0 across comprehensive biological test sets.

![Alignment Hit Distribution](benchmarks/v0.1.0/hit_distribution.png)

LOSAT achieves exact row-by-row, coordinate-for-coordinate, and score-for-score parity with NCBI BLAST+ across BLASTN, BLASTP, and TBLASTX pairwise comparisons.

### Execution Speed
Benchmarked on an Intel Core i9-14900HX comparing NCBI BLAST+ 2.17.0, native LOSAT, and WebAssembly modes (single-threaded and 8-thread pool):

![Execution Time Benchmark](benchmarks/v0.1.0/execution_time.png)

- **Native**: Matches or exceeds NCBI BLAST+ execution speeds across all three alignment modes.
- **WebAssembly**: Threaded WASM provides near-native scaling, allowing compute-intensive genomic alignments directly in sandboxed or client-side environments.

*For complete reproducible datasets, scripts, and methodology, see the [Benchmark Documentation](benchmarks/v0.1.0/README.md).*

---

## WebAssembly (Wasm) Integration

LOSAT can be built as a standalone WASI module or embedded into web browsers:

```bash
cd LOSAT

# Build single-threaded WASI command
cargo build --release --target wasm32-wasip1 --no-default-features

# Build multithreaded WASI command (requires shared-memory runtime)
cargo build --release --target wasm32-wasip1-threads --features wasm-threads
```

LOSAT powers the client-side comparative genomic alignment engine in [gbdraw](https://github.com/satoshikawato/gbdraw), enabling interactive circular and linear genome visualization without any backend server requirements.

---

## Scope & Differences from NCBI BLAST+

- **Intended Use Case**: LOSAT is engineered for local pairwise sequence comparisons (`-query` vs `-subject`). It is not designed to replace large indexed database searches (`-db` created via `makeblastdb`). For multi-gigabase database queries against NR/NT, continue using NCBI BLAST+ or DIAMOND.
- **TBLASTX Subject Genetic Code**: In NCBI BLAST+, local `-subject` searches inadvertently default to genetic code 1 (Standard) regardless of command-line flags. LOSAT correctly translates the subject according to `--db-gencode`.
- **Fail-Fast Configuration**: Unsupported flags and tasks fail fast with clear error messages rather than silently falling back to uncertified defaults.

---

## Documentation & Verification

- [Release Readiness & Scope](docs/v0.1.0_scope.md)
- [Verification & Parity Specifications](docs/release/pure_rust_runtime_v0.1.0_certification.md)
- [CLI Migration Details](docs/cli_v2_migration.md)
- [Developer & Contributor Guidelines](AGENTS.md)

---

## References

1. Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. *Journal of Molecular Biology*, 215(3), 403–410.
2. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.
3. [NCBI C++ Toolkit Repository](https://github.com/ncbi/ncbi-cxx-toolkit-public)

---

## License

This project is licensed under the [MIT License](LICENSE).
