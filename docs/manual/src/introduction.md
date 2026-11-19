<img class="logo-for-light" src="images/losat-logo.svg" alt="LOSAT" width="420">
<img class="logo-for-dark" src="images/losat-logo-white.svg" alt="LOSAT" width="420">

# Introduction

LOSAT is a pure-Rust reimplementation of NCBI BLAST+ (`2.17.0` compatible CLI and output). Run searches directly as `losat blastn`, `losat blastp`, or `losat makeblastdb` using the same options you already pass to NCBI BLAST+, or run `losat --install-shims` ([Coming from NCBI BLAST+](from-ncbi-blast.md)) so pipelines calling `blastn` or `makeblastdb` on `PATH` switch to LOSAT with zero script changes.

Because LOSAT has no C or C++ runtime dependencies, pairwise searches (`-subject`), database builds (`makeblastdb`), and database searches (`-db`) work out of the box from one static binary or library. On the release test set across the supported programs and options, outputs match NCBI BLAST+ 2.17.0 byte for byte while cutting wall time and peak memory on standard single-node workloads and multi-gigabyte databases ([Threads, speed and memory](performance.md)).

## What v0.2.0 provides

| Category | Program | What it does |
| :--- | :--- | :--- |
| **Search** | [`blastn`](programs/blastn.md) | Nucleotide vs. nucleotide (`megablast` default, `blastn`, `dc-megablast`, `blastn-short`). |
| **Search** | [`blastp`](programs/blastp.md) | Protein vs. protein (`blastp` default, `blastp-fast`, `blastp-short`). |
| **Search** | [`blastx`](programs/blastx.md) | Translated nucleotide vs. protein (`blastx` default, `blastx-fast`). |
| **Search** | [`tblastn`](programs/tblastn.md) | Protein vs. translated nucleotide (`tblastn` default, `tblastn-fast`; supports `-in_pssm`). |
| **Search** | [`tblastx`](programs/tblastx.md) | Translated nucleotide vs. translated nucleotide. |
| **Profile & Iterative** | [`psiblast`](programs/psiblast.md) | Iterative PSSM search from FASTA (`-query`) or checkpoint (`-in_pssm`). |
| **Profile & Iterative** | [`rpsblast`](programs/rpsblast.md) | Protein queries against conserved-domain PSSM databases such as CDD. |
| **Profile & Iterative** | [`rpstblastn`](programs/rpstblastn.md) | Translated nucleotide queries against profile databases such as CDD. |
| **Database & Masking** | [`makeblastdb`](programs/makeblastdb.md) | Builds v4/v5 BLAST databases directly from FASTA (supports `.gz`, `.bz2`, `.zst`). |
| **Database & Masking** | [`blastdbcmd`](programs/blastdbcmd.md) | Queries database metadata, extracts sequences, and prints mask ranges. |
| **Database & Masking** | [`dustmasker`](programs/dustmasker.md) | Computes low-complexity masking intervals for nucleotide sequences. |
| **Database & Masking** | [`segmasker`](programs/segmasker.md) | Computes low-complexity masking intervals for protein sequences. |
| **Database & Masking** | [`convert2blastmask`](programs/convert2blastmask.md) | Converts lowercase FASTA masking for `makeblastdb -mask_data`. |
| **Annotation** | [`rpsbproc`](programs/rpsbproc.md) | Post-processes `-outfmt 11` archives from NCBI's `rpsblast` or `rpstblastn` into conserved domain, site, and motif annotations. |

Every search program reads BLAST databases (`-db`) and, except `rpsblast` and `rpstblastn` (matching NCBI), FASTA subjects (`-subject`) and standard input (`-query -`).

Intentional semantic improvements and rare edge cases where LOSAT differs from NCBI BLAST+ are listed in [Differences from NCBI BLAST+](differences.md). Options that LOSAT intentionally declines to run (because they are unimplemented or trigger unsafe NCBI behaviour) fail cleanly with an error instead of silently diverging; see [What LOSAT refuses](refusals.md).

## Distribution and ecosystem

LOSAT v0.2.0 launches pre-built native and WASI binaries on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) and Linux/macOS packages on [Bioconda](https://anaconda.org/bioconda/losat), with the **`losat`** Rust crate (crates.io) and **`pylosat`** Python bindings (PyPI, with zero-copy Arrow tables) scheduled for publication immediately after v0.2.0. Browser and Node.js WebAssembly modules (`WasmEngine` and WASI CLI builds) are built from the same codebase; see [Using LOSAT as a library](library.md) and [Platforms and WebAssembly](platforms.md).

## What is not in v0.2.0

Planned for later releases:

- `deltablast`, `makeprofiledb`, `blast_formatter`, and `windowmasker`
- `psiblast -in_msa` (`-msa_master_idx`, `-ignore_msa_master`) and `-phi_pattern`
- `psiblast` output formats other than `0`, `6`, `7`, `10`, and `20` (with `-in_pssm`: `0`, `6`, `7`)
- `rpsblast` output formats other than `0`, `6`, `7`, and `10`, and `rpstblastn` output formats other than `0`–`4`, `6`, `7`, and `10` (until `-outfmt 11` lands natively, feed `losat rpsbproc` an archive produced by NCBI `rpsblast`/`rpstblastn`)
- `-html` output

LOSAT does not include NCBI's network/VDB/helper utilities (`blastdb_aliastool`, `blast_formatter_vdb`, `blast_vdb_cmd`, `blastn_vdb`, `tblastn_vdb`, `blastdbcheck`, `makembindex`, `legacy_blast.pl`, `update_blastdb.pl`, `get_species_taxids.sh`, `cleanup-blastdb-volumes.py`) or remote/GI/IPG options (`-remote`, `-entrez_query`, `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`, `-import_search_strategy`).

## How this manual is organised

- **Getting started:** [Installation](installation.md), [Quick start](quick-start.md), [Common workflows](workflows.md), [Coming from NCBI BLAST+](from-ncbi-blast.md), and [Coming from LOSAT v0.1.0](from-losat-v0-1.md).
- **Programs:** Reference pages for all 14 CLI commands with options, formats, and examples.
- **Libraries & WebAssembly:** Overviews and examples for [Python (`pylosat`)](libraries/python.md), [Rust (`losat`)](libraries/rust.md), and [WebAssembly and Browser](libraries/webassembly.md).
- **Reference:** [BLAST databases](databases.md), [Output formats](output-formats.md), [Environment variables](environment.md), [Threads, speed and memory](performance.md), [Platforms and WebAssembly](platforms.md), [Differences from NCBI BLAST+](differences.md), [What LOSAT refuses](refusals.md), and [Troubleshooting](troubleshooting.md).

Run `losat <program> -help` at any time for full command-line flags.

## Reporting issues and License

If LOSAT fails or diverges from NCBI BLAST+ 2.17.0 in a case not documented in [Differences](differences.md) or [What LOSAT refuses](refusals.md), open an issue at [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues) with the exact command, `losat --version`, and a small test case (no confidential data; report security issues via the repository's Security tab).

Released under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/main/LICENSE).
