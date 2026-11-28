<img class="logo-for-light" src="images/losat-logo.svg" alt="LOSAT" width="420">
<img class="logo-for-dark" src="images/losat-logo-white.svg" alt="LOSAT" width="420">

# Introduction

LOSAT is a BLAST+ compatible toolkit in one executable, `losat`. It takes the command lines of NCBI BLAST+ 2.17.0, so the search you already run with `blastn`, `blastp` or `makeblastdb` runs as `losat blastn`, `losat blastp` or `losat makeblastdb`, with the same options and the same output (apart from the version line, which names LOSAT). Link the program names to `losat` (see [Coming from NCBI BLAST+](from-ncbi-blast.md)) and tools that look for `blastn` on your `PATH` use LOSAT without changes.

LOSAT is written in pure Rust. Searching and building databases need nothing else: no NCBI BLAST+, no C or C++ libraries, no database to format before a pairwise search. Two things come from NCBI: the CDD data files that `rpsbproc` reads and NCBI's prebuilt databases, which are downloaded from NCBI (see [BLAST databases](databases.md)). A `-subject` search compares your query with a FASTA file directly. When you do want a database, `losat makeblastdb` writes the files that NCBI's `makeblastdb` writes, and the search programs read databases made by either tool, including NCBI's distributed ones such as `swissprot`.

## What v0.2.0 provides

| Kind | Program | What it does |
|:---|:---|:---|
| Search | [`blastn`](programs/blastn.md) | Nucleotide query against nucleotide subjects. Tasks `megablast` (default), `blastn`, `dc-megablast`, `blastn-short`, `rmblastn`. |
| Search | [`blastp`](programs/blastp.md) | Protein query against protein subjects. Tasks `blastp` (default), `blastp-fast`, `blastp-short`. |
| Search | [`blastx`](programs/blastx.md) | Translated nucleotide query against protein subjects. Tasks `blastx` (default), `blastx-fast`. |
| Search | [`tblastn`](programs/tblastn.md) | Protein query against translated nucleotide subjects. Tasks `tblastn` (default), `tblastn-fast`; a PSSM query with `-in_pssm`. |
| Search | [`tblastx`](programs/tblastx.md) | Translated nucleotide query against translated nucleotide subjects. |
| Search | [`psiblast`](programs/psiblast.md) | Position-specific iterated search, from a protein query or a PSSM. |
| Search | [`rpsblast`](programs/rpsblast.md) | Protein queries against a database of profiles such as CDD. |
| Search | [`rpstblastn`](programs/rpstblastn.md) | Translated nucleotide queries against a database of profiles. |
| Databases and masks | [`makeblastdb`](programs/makeblastdb.md) | Builds a BLAST database (version 5 by default, or 4) from FASTA, also from a `.gz`, `.bz2` or `.zst` file. |
| Databases and masks | [`blastdbcmd`](programs/blastdbcmd.md) | Reads a BLAST database: information, sequences, masks. |
| Databases and masks | [`dustmasker`](programs/dustmasker.md) | Finds low-complexity regions in nucleotide sequences. |
| Databases and masks | [`windowmasker`](programs/windowmasker.md) | Makes the unit counts file of a genome, and masks the repeats of nucleotide sequences with it. |
| Databases and masks | [`segmasker`](programs/segmasker.md) | Finds low-complexity regions in protein sequences. |
| Databases and masks | [`convert2blastmask`](programs/convert2blastmask.md) | Turns lower-case masking in FASTA into masking information for `makeblastdb`. |
| Databases and masks | [`makembindex`](programs/makembindex.md) | Builds the MegaBLAST index of a nucleotide database. |
| Databases and masks | [`makeprofiledb`](programs/makeprofiledb.md) | Builds a profile database for `rpsblast`, `rpstblastn` and DELTA-BLAST (such as CDD) from PSSM files (`-dbtype rps` and `-dbtype delta`). |
| Databases and masks | [`blastdb_aliastool`](programs/blastdb_aliastool.md) | Makes alias files that combine databases or keep a list of their sequences; converts GI and sequence id lists to binary lists. |
| Databases and masks | [`blastdbcheck`](programs/blastdbcheck.md) | Checks BLAST databases for damage: files, sequences, ISAM indexes, alias files, taxonomy ids. |
| Databases and masks | [`cleanup-blastdb-volumes.py`](programs/cleanup-blastdb-volumes.md) | Removes orphaned BLAST database volume files not listed in a database's `.pal` or `.nal` alias file. |
| Archives | [`blast_formatter`](programs/blast_formatter.md) | Writes the reports of a search again from its `-outfmt 11` BLAST archive, without searching again. |
| Domain annotation | [`rpsbproc`](programs/rpsbproc.md) | Annotates conserved domains, sites and motifs from the `-outfmt 11` results of `rpsblast` and `rpstblastn` (see its page). |

The search programs read BLAST databases (`-db`) and, except `rpsblast` and `rpstblastn` (as in NCBI), FASTA subjects (`-subject`). Each can read its query from the standard input (`-query -`).

On the release test set, LOSAT's output is the same as NCBI BLAST+ 2.17.0's, byte for byte, apart from the line that names the program's version (`BLASTN (LOSAT 0.2.0)` where NCBI writes `BLASTN 2.17.0+`) and the prefix of an error message (`LOSAT engine error:` where NCBI writes `BLAST engine error:`). The few other places where it differs, on purpose, are listed in [Differences from NCBI BLAST+](differences.md). Options LOSAT does not implement are refused with an error, never ignored; see [What LOSAT refuses](refusals.md).

## Where to get it

LOSAT v0.2.0 is published as pre-built files for Linux, macOS, Windows and WASI on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases), as packages for Linux and macOS on [Bioconda](https://anaconda.org/bioconda/losat), and as the crate `losat` on [crates.io](https://crates.io/crates/losat), which `cargo install` builds. See [Installation](installation.md).

The Rust library (`losat::api`) is the same crate, and the Python binding `pylosat` is on [PyPI](https://pypi.org/project/pylosat/). The WebAssembly engine for web pages is built from the repository. See [Using LOSAT as a library](library.md).

## What is not in v0.2.0

These are on the list for a later release:

- the program `deltablast`, and `makeprofiledb -dbtype cobalt`;
- `blast_formatter` output formats other than `0` to `4`, `6` to `12` and `20`;
- `psiblast -phi_pattern` with more than one round, with `-html`, or with output formats other than `0`, `6`, `7`, `10` and `20`,
  and in the library and pylosat;
- the `psiblast` output formats `13` to `16`.

LOSAT does not provide these NCBI programs: `blast_formatter_vdb`, `blast_vdb_cmd`, `blastn_vdb`, `tblastn_vdb`, `legacy_blast.pl`, `update_blastdb.pl` and `get_species_taxids.sh`.

Options that need NCBI's servers are not available: `-remote` and `-entrez_query` (all programs). `-import_search_strategy` runs in `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`; `psiblast`, `rpsblast` and `rpstblastn` refuse it.

## How this manual is organised

- **Getting started.** [Installation](installation.md), a five-minute [Quick start](quick-start.md), recipes for everyday tasks in [Common workflows](workflows.md), and two pages for people who already use another version: [Coming from NCBI BLAST+](from-ncbi-blast.md) and [Coming from LOSAT v0.1.0](from-losat-v0-1.md).
- **Programs.** One page for each program, with its options, input, output and examples.
- **Libraries, WebAssembly and Web GUI.** [Using LOSAT as a library](library.md), with pages for [Python (`pylosat`)](libraries/python.md), [Rust (`losat`)](libraries/rust.md), [WebAssembly and Browser](libraries/webassembly.md) and [LOSAT Web (Browser GUI)](web.md).
- **Reference.** [BLAST databases](databases.md), [Output formats](output-formats.md), [Environment variables](environment.md), [Threads, speed and memory](performance.md), [Platforms and WebAssembly](platforms.md), [Differences from NCBI BLAST+](differences.md), [What LOSAT refuses](refusals.md) and [Troubleshooting](troubleshooting.md).

Every program also describes itself: `losat <program> -help` prints the options of the version you have installed.

## Reporting problems

If LOSAT fails, or its output differs from NCBI BLAST+ 2.17.0 for the same command in a case that [Differences from NCBI BLAST+](differences.md) and [What LOSAT refuses](refusals.md) do not list, open an issue at [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues). Include the exact command line, the output of `losat --version`, and a small input that shows the problem. See [Troubleshooting](troubleshooting.md) for what else helps. Do not post confidential data. Report a security vulnerability privately, through the repository's Security tab, not in a public issue.

## Citation

If you use LOSAT in published work, cite LOSAT ([CITATION.cff](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/CITATION.cff)) together with the foundational BLAST / BLAST+ papers:

- Kawato, S. LOSAT (LOcal Sequence Alignment Tool). <https://github.com/satoshikawato/LOSAT>.
- Altschul, S. F. et al. (1990). Basic local alignment search tool. *J. Mol. Biol.*, 215(3), 403–410.
- Altschul, S. F. et al. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res.*, 25(17), 3389–3402.
- Camacho, C. et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/THIRD_PARTY_NOTICES).
