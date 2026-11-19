# Using LOSAT as a library

In addition to the command line tool, LOSAT runs directly inside your own code with no subprocess overhead:

- **Same engine, same outputs:** In-process APIs pass your FASTA sequences, options, PSSMs, or BLAST database names straight to the LOSAT engine using the exact CLI option parser, returning byte-identical reports alongside typed HSP hits and Karlin-Altschul statistics.
- **No external binaries required:** Searches execute completely in-memory or against local BLAST database files without calling NCBI BLAST+.

## Choose an interface

| Interface | Target Package | Supported Searches & Tools | Guide |
| :--- | :--- | :--- | :--- |
| **Python** | `pylosat` (PyO3 + Apache Arrow) | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, plus `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask` | [Python (`pylosat`)](libraries/python.md) |
| **Rust** | `losat` crate (`losat::api`) | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, plus `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker`, `run_convert2blastmask` | [Rust (`losat`)](libraries/rust.md) |
| **WebAssembly / JS / TS** | `losat` npm host (`WasmEngine`) & WASI modules | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` (browser reactors) and full CLI (WASI command builds) | [WebAssembly and Browser](libraries/webassembly.md) |

## Availability in v0.2.0

LOSAT v0.2.0 launches its binaries on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) and [Bioconda](https://anaconda.org/bioconda/losat). The **`losat`** Rust crate (crates.io) and **`pylosat`** Python package (PyPI) are scheduled for publication immediately following the v0.2.0 release; until published there, both can be used directly from the [GitHub repository](https://github.com/satoshikawato/LOSAT).
