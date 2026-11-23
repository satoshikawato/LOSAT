<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/losat-logo-white.svg">
    <img alt="LOSAT" src="docs/assets/losat-logo.svg" width="420">
  </picture>
</p>

<h3 align="center">LOSAT: LOcal Sequence Alignment Tool. NCBI BLAST+ 2.17.0's local searches, written in Rust, for the command line, Rust, Python and the browser</h3>

<p align="center">
  <a href="https://satoshikawato.github.io/LOSAT/">User manual</a> •
  <a href="benchmarks/v0.2.0/README.md">v0.2.0 benchmarks</a> •
  <a href="LOSAT/README.md">Rust library</a> •
  <a href="python/README.md">Python package</a> •
  <a href="web/packages/losat/README.md">WebAssembly package</a>
</p>

---

**LOSAT** runs the local search programs and database tools of NCBI BLAST+ 2.17.0 with their command lines and writes their output. It is one executable written in Rust: it needs no BLAST+ installation, and it does not call NCBI's programs or libraries. `losat --install-shims` puts it in place of BLAST+ in an existing pipeline, so scripts that call `blastn` or `makeblastdb` run unchanged. The same engine is also a Rust library, a Python package (`pylosat`) and a WebAssembly package for the browser.

On the v0.2.0 release test set, LOSAT writes the same output as NCBI BLAST+ 2.17.0, apart from the [differences listed in the manual](docs/manual/src/differences.md) (among them the line that names the program's version, which says `BLASTN (LOSAT 0.2.0)` where NCBI writes `BLASTN 2.17.0+`), and the same output at every `-num_threads`. It searches FASTA subjects (`-subject`) on all the threads you give it, where NCBI BLAST+ uses one, reads gzip, bzip2 and Zstandard FASTA directly, and keeps its memory small on large databases: a search of `core_nt` (319.4 GB) peaked at 407 MiB of resident memory, where NCBI BLAST+ peaked at 257.4 GiB (LOSAT `095ff0c16`, an earlier development build; NCBI BLAST+ maps the database files into memory, so its figure includes file pages the system can drop).

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="benchmarks/v0.2.0/figures/fig1_standard_workloads_dark.png">
    <img alt="Wall time and peak memory of LOSAT 0.2.0 and NCBI BLAST+ 2.17.0 in 12 search modes, 1 and 16 threads" src="benchmarks/v0.2.0/figures/fig1_standard_workloads.png" width="100%">
  </picture>
</p>

*Twelve search modes on a dual Xeon Gold 6338 workstation (64 cores), LOSAT at commit `93ca2f652` against NCBI BLAST+ 2.17.0. At 16 threads LOSAT was faster in all twelve: 629 times on `tblastx` (4.72 s against 2,970.75 s), 1.7 to 18.3 times on ten other modes, and 1.02 times on `psiblast`. Where the subject is a FASTA file, LOSAT ran `-subject` and NCBI BLAST+ `-db` on a database made from the same file. Inputs, commands and every number: [benchmarks/v0.2.0](benchmarks/v0.2.0/README.md).*

---

## Installation

### Bioconda (Linux and macOS)

```bash
conda install -c bioconda -c conda-forge losat
# or
mamba install -c bioconda -c conda-forge losat
```

With pixi, add the `bioconda` and `conda-forge` channels to the workspace, then `pixi add losat`.

### Release files

The [v0.2.0 release](https://github.com/satoshikawato/LOSAT/releases) has one archive per platform, the source crate, `SHA256SUMS`, a CycloneDX SBOM (`LOSAT-0.2.0.cdx.json`) and GitHub build attestations:

| Platform | Target | Archive | Program |
| --- | --- | --- | --- |
| Linux x86_64 | `x86_64-unknown-linux-gnu` | `.tar.gz` | `losat` |
| Linux aarch64 | `aarch64-unknown-linux-gnu` | `.tar.gz` | `losat` |
| macOS arm64 | `aarch64-apple-darwin` | `.tar.gz` | `losat` |
| macOS x86_64 | `x86_64-apple-darwin` | `.tar.gz` | `losat` |
| Windows x86_64 | `x86_64-pc-windows-msvc` | `.zip` | `losat.exe` |
| WASI serial | `wasm32-wasip1` | `.tar.gz` | `losat.wasm` |
| WASI threads | `wasm32-wasip1-threads` | `.tar.gz` | `losat.wasm` |

The macOS and Windows programs are not code-signed; [Installation](docs/manual/src/installation.md) says how to run them.

### crates.io

With Rust 1.92 or later:

```bash
cargo install --locked losat
```

`--locked` builds with the dependency versions of the release build.

### Build from source

Requires Rust 1.92 or later.

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/LOSAT
cargo build --release
# the program is target/release/losat
```

### Platforms

Linux (x86_64, aarch64), macOS (arm64, x86_64), Windows (x86_64) and the WASI serial module are **Supported**: they match NCBI BLAST+ 2.17.0 on the release test set (apart from the version line, which names LOSAT). The WASI threads module is **Supported** and tested in CI, outside that comparison. LOSAT Web, the search page that runs in the browser, is a **Preview**. Details, and the one macOS difference (the last digit of a full-precision E value): [Platforms](docs/manual/src/platforms.md).

---

## Quick start

Run any program as `losat <program>`, or install links with the NCBI names so that existing scripts call LOSAT:

```bash
# Build a database and search it on 16 threads (the query can be gzip, bzip2 or Zstandard)
losat makeblastdb -in 16s_all.fna -dbtype nucl -out 16s_db
losat blastn -query reads.fna.gz -db 16s_db -evalue 1e-5 -outfmt 6 -num_threads 16

# Compare two genomes directly, without a database
losat blastn -query genome_a.fna -subject genome_b.fna -outfmt 6 -num_threads 16

# PSI-BLAST and a conserved-domain search
losat psiblast -query proteins.faa -db swissprot -num_iterations 3 -outfmt 6 -num_threads 8
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 6 -num_threads 8

# Use LOSAT under the NCBI names (blastn, blastp, makeblastdb, ...)
losat --install-shims ~/.local/share/losat-shims
export PATH="$HOME/.local/share/losat-shims:$PATH"
blastp -version
```

<a id="common-command-line-options"></a>
<a id="database-searches"></a>

Every search program reads BLAST databases (`-db`, versions 4 and 5, found through `BLASTDB` as in NCBI BLAST+) and its query from standard input (`-query -`). All except `rpsblast` and `rpstblastn` also take FASTA subjects (`-subject`), as in NCBI BLAST+. Options, output formats and examples for each program are in the [user manual](docs/manual/src/introduction.md).

---

## Programs

| Program | What it does | Output formats (`-outfmt`) | Manual |
| --- | --- | --- | --- |
| `blastn` | Nucleotide against nucleotide (`megablast`, `blastn`, `dc-megablast`, `blastn-short`, `rmblastn`) | `0`–`18`, `20` | [blastn](docs/manual/src/programs/blastn.md) |
| `blastp` | Protein against protein (`blastp`, `blastp-fast`, `blastp-short`) | `0`–`16`, `18`, `20` | [blastp](docs/manual/src/programs/blastp.md) |
| `blastx` | Translated nucleotide against protein (`blastx`, `blastx-fast`) | `0`–`16`, `18`, `20` | [blastx](docs/manual/src/programs/blastx.md) |
| `tblastn` | Protein or PSSM (`-in_pssm`) against translated nucleotide (`tblastn`, `tblastn-fast`) | `0`–`16`, `18`, `20` | [tblastn](docs/manual/src/programs/tblastn.md) |
| `tblastx` | Translated nucleotide against translated nucleotide | `0`–`16`, `18`, `20` | [tblastx](docs/manual/src/programs/tblastx.md) |
| `psiblast` | Iterated protein search with a PSSM (`-num_iterations`, `-in_pssm`, `-in_msa`, `-out_pssm`) | `0`, `5`–`11`, `20` | [psiblast](docs/manual/src/programs/psiblast.md) |
| `rpsblast` | Protein against a PSSM database such as CDD | `0`–`16`, `18`, `20` | [rpsblast](docs/manual/src/programs/rpsblast.md) |
| `rpstblastn` | Translated nucleotide against a PSSM database | `0`–`16`, `18`, `20` | [rpstblastn](docs/manual/src/programs/rpstblastn.md) |
| `makeblastdb` | Builds BLAST databases (versions 5 and 4; `-parse_seqids`, `-taxid`, `-mask_data`) from FASTA, also compressed | | [makeblastdb](docs/manual/src/programs/makeblastdb.md) |
| `blastdbcmd` | Prints sequences and information from a BLAST database (`-entry`, `-entry_batch`, `-info`, `-outfmt`) | | [blastdbcmd](docs/manual/src/programs/blastdbcmd.md) |
| `dustmasker` | Masks low-complexity nucleotide regions | intervals, FASTA, mask files for `makeblastdb` | [dustmasker](docs/manual/src/programs/dustmasker.md) |
| `windowmasker` | Makes the unit counts file of a genome (`-mk_counts`), converts it (`-convert`), masks repeats of nucleotide sequences with it (`-ustat`) | counts files, intervals, FASTA, mask files for `makeblastdb` | [windowmasker](docs/manual/src/programs/windowmasker.md) |
| `segmasker` | Masks low-complexity protein regions | intervals, FASTA, mask files for `makeblastdb` | [segmasker](docs/manual/src/programs/segmasker.md) |
| `convert2blastmask` | Turns lower-case masked FASTA into mask files for `makeblastdb` | mask files, intervals | [convert2blastmask](docs/manual/src/programs/convert2blastmask.md) |
| `makembindex` | Builds the MegaBLAST index of a nucleotide database | index files | [makembindex](docs/manual/src/programs/makembindex.md) |
| `makeprofiledb` | Builds profile databases for `rpsblast` and `rpstblastn` (such as CDD) from PSSM files | | [makeprofiledb](docs/manual/src/programs/makeprofiledb.md) |
| `rpsbproc` | Annotates conserved domains, sites and motifs from RPS-BLAST results, as NCBI's `rpsbproc` 0.5.1 | domain report | [rpsbproc](docs/manual/src/programs/rpsbproc.md) |

`rpsbproc` reads `-outfmt 11` archives, which LOSAT's `rpsblast` and `rpstblastn` write as NCBI BLAST+'s do: `losat rpsblast -query proteins.faa -db Cdd -outfmt 11 | losat rpsbproc -d cdd_data` ([rpsbproc](docs/manual/src/programs/rpsbproc.md)).

---

## Benchmarks and pipelines

Published tools and workflow managers run unchanged with LOSAT under the NCBI names (`--install-shims`), and their BLAST steps wrote the same output in every case compared (where NCBI BLAST+'s own output changes with the thread count, as in BUSCO and QIIME 2, the same as its 1-thread output), apart from the line that names the program's version in the formats that have one:

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="benchmarks/v0.2.0/figures/fig3_downstream_tools_dark.png">
    <img alt="Time of the BLAST steps of published tools with LOSAT 0.2.0 and with NCBI BLAST+ 2.17.0" src="benchmarks/v0.2.0/figures/fig3_downstream_tools.png" width="100%">
  </picture>
</p>

[benchmarks/v0.2.0](benchmarks/v0.2.0/README.md) has all four figures (standard searches, scaling with threads, queries and database size, published tools, WebAssembly), the inputs, and the measured values as TSV files.

---

## Libraries and WebAssembly

- **Rust (`losat`):** `losat::api` runs a search on FASTA records in memory or on a BLAST database and returns the command line's reports byte for byte, plus one typed record per HSP. Add it with `cargo add losat`. See [LOSAT/README.md](LOSAT/README.md) and the [Rust library page](docs/manual/src/libraries/rust.md).
- **Python (`pylosat`):** searches from Python, with the hits as an Arrow table for pandas or Polars. Install it with `pip install pylosat`. See [python/README.md](python/README.md) and the [Python page](docs/manual/src/libraries/python.md).
- **WebAssembly & LOSAT Web:** the npm package `losat` (built from this repository) runs the engine in web pages and Node.js, and the browser preview LOSAT Web runs on it; `losat.wasm` is the command for WASI hosts. See [web/packages/losat/README.md](web/packages/losat/README.md), the [WebAssembly page](docs/manual/src/libraries/webassembly.md), and the [LOSAT Web guide](docs/manual/src/web.md).

---

<a id="scope--differences-from-ncbi-blast"></a>

## Documentation

- **[User manual](https://satoshikawato.github.io/LOSAT/)** (source in `docs/manual/src/`):
  - [Introduction](docs/manual/src/introduction.md) • [Quick start](docs/manual/src/quick-start.md) • [Common workflows](docs/manual/src/workflows.md) • [Installation](docs/manual/src/installation.md)
  - [Coming from NCBI BLAST+](docs/manual/src/from-ncbi-blast.md) • [Upgrading from LOSAT v0.1](docs/manual/src/from-losat-v0-1.md) • [Threads, speed and memory](docs/manual/src/performance.md)
  - [Differences from NCBI BLAST+](docs/manual/src/differences.md) • [What LOSAT refuses](docs/manual/src/refusals.md) • [Databases](docs/manual/src/databases.md) • [Output formats](docs/manual/src/output-formats.md) • [Environment variables](docs/manual/src/environment.md) • [LOSAT Web](docs/manual/src/web.md)
- **Not in v0.2.0:** `deltablast`, `makeprofiledb -binary` and `-dbtype delta` or `cobalt`, `blast_formatter`, `psiblast -phi_pattern`, and `-html` with `-outfmt 1` to `4` or in `psiblast`, `rpsblast` and `rpstblastn`; options that need NCBI's servers or its GI and IPG lists (`-remote`, `-entrez_query`, `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`) and `-import_search_strategy`. LOSAT stops with an error naming the option before it searches. See [What LOSAT refuses](docs/manual/src/refusals.md).
- [v0.2.0 release record](docs/release/v0.2.0.md) • [Contributing](CONTRIBUTING.md)

## Citation & References

If you use LOSAT in published work, cite LOSAT ([CITATION.cff](CITATION.cff)) together with the foundational BLAST / BLAST+ papers:

1. Kawato, S. LOSAT (LOcal Sequence Alignment Tool). Software repository: <https://github.com/satoshikawato/LOSAT>. See [CITATION.cff](CITATION.cff).
2. Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. *Journal of Molecular Biology*, 215(3), 403–410.
3. Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., & Lipman, D. J. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research*, 25(17), 3389–3402.
4. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.
5. [NCBI C++ Toolkit](https://github.com/ncbi/ncbi-cxx-toolkit-public), whose BLAST+ code LOSAT follows.

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES).
