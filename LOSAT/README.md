# LOSAT

LOSAT is a pure-Rust reimplementation of NCBI BLAST+ 2.17.0. It runs BLAST+ command lines and matches NCBI BLAST+ 2.17.0 on its [release test set](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/release/v0.2.0.md) (apart from the version line and error prefixes, which name LOSAT). This crate builds the `losat` command line and the library `losat::api`.

- **User manual:** [satoshikawato.github.io/LOSAT](https://satoshikawato.github.io/LOSAT/)
- **Repository:** [github.com/satoshikawato/LOSAT](https://github.com/satoshikawato/LOSAT)
- **Issues:** [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues)

## Programs

| Program | What it does |
|:---|:---|
| `blastn` | Nucleotide search; tasks `megablast` (default), `blastn`, `dc-megablast`, `blastn-short`, `rmblastn` |
| `blastp` | Protein search; tasks `blastp` (default), `blastp-fast`, `blastp-short` |
| `blastx` | Translated nucleotide query against proteins; tasks `blastx` (default), `blastx-fast` |
| `tblastn` | Protein query against translated nucleotides; tasks `tblastn` (default), `tblastn-fast`; a PSSM query with `-in_pssm` |
| `tblastx` | Translated nucleotide query against translated nucleotides |
| `psiblast` | Iterative protein search from a FASTA query, an alignment (`-in_msa`), or a PSSM (`-in_pssm`); writes PSSMs (`-out_pssm`, `-out_ascii_pssm`) |
| `deltablast` | Iterative protein search whose first round uses the conserved domains of the query (`-rpsdb`, default `cdd_delta`) |
| `rpsblast` | Protein queries against a profile database such as CDD |
| `rpstblastn` | Translated nucleotide queries against a profile database |
| `rpsbproc` | Conserved-domain annotation of RPS-BLAST results, as NCBI's `rpsbproc` 0.5.1 |
| `makeblastdb` | Builds BLAST databases from FASTA, also from `.gz`, `.bz2` and `.zst` files |
| `blastdbcmd` | Reads sequences and information from BLAST databases |
| `dustmasker` | Masks low-complexity regions of nucleotide sequences |
| `segmasker` | Masks low-complexity regions of protein sequences |
| `windowmasker` | Masks repeats of nucleotide sequences with a unit counts file (`-ustat`) |
| `convert2blastmask` | Turns lower-case masking in FASTA into masking information for `makeblastdb -mask_data` |
| `makembindex` | Builds the MegaBLAST index of a nucleotide database |
| `makeprofiledb` | Builds profile databases for `rpsblast` and `rpstblastn` (such as CDD) from PSSM files |

Options that need NCBI's servers (`-remote`, `-entrez_query`) or GI lists (`-gilist`), and the few that LOSAT does not run, stop with an error instead of changing the results. See [What LOSAT refuses](https://satoshikawato.github.io/LOSAT/refusals.html).

## Installing the command line

LOSAT v0.2.0 is published on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) (Linux x86_64 and aarch64, macOS arm64 and x86_64, Windows x86_64, WASI serial and threads, and this source crate), on [crates.io](https://crates.io/crates/losat) and on [Bioconda](https://anaconda.org/bioconda/losat) (Linux and macOS).

```bash
conda install -c bioconda -c conda-forge losat
# or, with Rust 1.92 or later:
cargo install --locked losat
```

`--locked` builds with the dependency versions of the release build.

To build it from source, with Rust 1.92 or later:

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/LOSAT
cargo build --release
```

The binary is `target/release/losat`. Run a program as `losat <program>`, or make the NCBI names with `losat --install-shims` so that existing scripts call LOSAT:

```bash
losat blastn -query query.fna -subject subject.fna -outfmt 6
losat makeblastdb -in proteins.faa -dbtype prot -out proteins
losat psiblast -query query.faa -db proteins -num_iterations 3 -outfmt 6
```

## Using `losat` as a Rust library

Add the crate with `cargo add losat`, or in `Cargo.toml`:

```toml
[dependencies]
losat = "0.2"
```

The library is the module `losat::api`. Every other module is the command line or the engine: hidden from the documentation and outside the compatibility promise. A search runs in your process with the command line's options and returns the command line's reports byte for byte, plus one typed record per HSP (`HspRecord`) and the statistics of each query (`QueryStats`):

```rust
use std::io::Write as _;
use losat::api::{read_fasta, run_collect, validate, Inputs, Program, Role};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let query = read_fasta(Program::Blastp, Role::Query, &std::fs::read("query.faa")?)?;
    let subject = read_fasta(Program::Blastp, Role::Subject, &std::fs::read("subject.faa")?)?;

    // The command line's options, checked by its parser (without the inputs, -out and -outfmt).
    let search = validate(Program::Blastp, &["-evalue", "1e-5", "-num_threads", "4"])?;
    let inputs = Inputs::new(&query, &subject).named("query.faa", "subject.faa");

    let outcome = run_collect(&search, inputs, &["6", "0"])?;
    std::io::stdout().write_all(&outcome.reports[0].bytes)?; // the -outfmt 6 report

    for hsp in &outcome.hsps {
        println!("query {} subject {}: {} bits, E = {:e}", hsp.q_idx, hsp.s_idx, hsp.bit_score, hsp.e_value);
    }
    Ok(())
}
```

- `validate(program, &argv)` checks the options and returns a reusable `Search`; `describe(program)` lists the options.
- `Inputs::database(&query, "swissprot")` searches a BLAST database, found as the command line finds `-db`. `Inputs::pssm` and `Inputs::pssm_database` take the PSSM of `-in_pssm` (`psiblast`, `tblastn`). `Inputs::msa` and `Inputs::msa_database` take a multiple sequence alignment for `psiblast` (`-in_msa`).
- `run_collect` returns the reports, the HSP records, the statistics and the diagnostics. `run_batches` gives them one of NCBI's query batches at a time. `run` writes the reports to your writers, and `run_to_file` is `-out` (and the way to `-outfmt 13` and `14`).
- `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker`, `run_windowmasker`, `run_convert2blastmask` and `run_makembindex` run the tools with your reader and writers as their standard streams.
- The library writes nothing to standard output or standard error, and its errors carry the command line's messages.

`cargo doc --no-deps --open` in `LOSAT/` shows the API. The [Rust library page](https://satoshikawato.github.io/LOSAT/libraries/rust.html) of the user manual has more.

## Cargo features

| Feature | Default | What it does |
|:---|:---|:---|
| `parallel` | on | Searches with more than one thread (`-num_threads`), through Rayon |
| `wasm-threads` | off | Turns on `parallel` for the `wasm32-wasip1-threads` target |

With `default-features = false` the build has no Rayon and searches on one thread. The other features are for LOSAT's own diagnostics.

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/THIRD_PARTY_NOTICES).
