# LOSAT

A pure-Rust reimplementation of [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) (`2.17.0` compatible CLI and output).

LOSAT replaces the standard BLAST+ executables (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`, `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`, `rpsbproc`) with a single Rust binary and embeddable library crate (`losat`), producing byte-identical outputs across standard tabular, pairwise, XML/ASN.1, and PSI/CDD profile workflows while substantially reducing runtime and peak resident memory.

- **Documentation & User Manual:** [satoshikawato.github.io/LOSAT](https://satoshikawato.github.io/LOSAT/)
- **Source Repository & Benchmarks:** [github.com/satoshikawato/LOSAT](https://github.com/satoshikawato/LOSAT)
- **API Reference:** [docs.rs/losat](https://docs.rs/losat)
- **Python Package (`pylosat`):** [pypi.org/project/pylosat](https://pypi.org/project/pylosat/)

## Supported Programs

| Category | Program | Notes |
| :--- | :--- | :--- |
| **Search** | `blastn` | Tasks `megablast` (default), `blastn`, `dc-megablast`, `blastn-short` |
| **Search** | `blastp` | Tasks `blastp` (default), `blastp-fast`, `blastp-short` |
| **Search** | `blastx` | Tasks `blastx` (default), `blastx-fast` |
| **Search** | `tblastn` | Tasks `tblastn` (default), `tblastn-fast`; supports `-in_pssm` |
| **Search** | `tblastx` | Translated nucleotide vs. translated nucleotide |
| **Profile & Iterative** | `psiblast` | Iterative search with `-query` or `-in_pssm` (`-out_pssm`, `-out_ascii_pssm`) |
| **Profile & Iterative** | `rpsblast` | Protein queries against conserved domain PSSM databases (such as CDD) |
| **Profile & Iterative** | `rpstblastn` | Translated nucleotide queries against profile databases |
| **Database & Masking** | `makeblastdb` | Builds BLAST v4/v5 databases directly from FASTA (supports `.gz`/`.bz2`/`.zst` inputs) |
| **Database & Masking** | `blastdbcmd` | Extracts sequences, metadata, and masks from BLAST databases |
| **Database & Masking** | `dustmasker` | Masks low-complexity nucleotide regions |
| **Database & Masking** | `segmasker` | Masks low-complexity protein regions |
| **Database & Masking** | `convert2blastmask` | Converts lowercase FASTA masking for `makeblastdb -mask_data` |
| **Annotation** | `rpsbproc` | Annotates domain hits, superfamilies, sites, and motifs from CDD archives |

Options that require NCBI network access (`-remote`, `-entrez_query`), GI/IPG server lists (`-gilist`, `-ipglist`), or unimplemented flags fail loudly instead of silently changing semantics. See [What LOSAT refuses](https://satoshikawato.github.io/LOSAT/refusals.html) in the user manual.

## Installing the CLI

The v0.2.0 release is distributed through [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) (Linux `x86_64`/`aarch64`, macOS Intel/Apple silicon, Windows `x86_64`, and WASI) and [Bioconda](https://anaconda.org/bioconda/losat), with crates.io publishing right after v0.2.0:

```bash
# Bioconda (Linux and macOS)
conda install -c bioconda -c conda-forge losat

# Or via Cargo (once published to crates.io after v0.2.0, requires Rust >= 1.92)
cargo install losat
```

Run any BLAST+ subcommand through `losat <program>`, or symlink `blastn`, `blastp`, `psiblast`, etc. to `losat` for drop-in pipeline compatibility:

```bash
losat blastn -query reads.fna -subject reference.fna -outfmt 6
losat makeblastdb -in proteins.faa -dbtype prot -out proteins_db
losat psiblast -query query.faa -db proteins_db -num_iterations 3 -outfmt 6
losat rpsblast -query proteins.faa -db Cdd -evalue 0.01 -outfmt 6
```

## Using `losat` as a Rust Library

Add `losat` to `Cargo.toml` (requires Rust **1.92** or later; default features enable multi-threading via `parallel`):

```toml
[dependencies]
losat = "0.2"
```

To build without Rayon multi-threading (e.g., single-thread or `wasm32` targets):

```toml
[dependencies]
losat = { version = "0.2", default-features = false }
```

The stable, supported public library surface lives in `losat::api`. Searches run directly in-process over FASTA bytes, prebuilt BLAST databases, or PSSM checkpoints, returning exact BLAST report bytes alongside typed `HspRecord` structs and Karlin-Altschul `QueryStats`.

```rust
use std::io::Write as _;
use losat::api::{read_fasta, run_collect, validate, Inputs, Program, Role};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let query = read_fasta(Program::Blastp, Role::Query, &std::fs::read("query.faa")?)?;
    let subject = read_fasta(Program::Blastp, Role::Subject, &std::fs::read("subject.faa")?)?;

    // Validated with the exact CLI option parser (excluding input/output routing flags).
    let search = validate(Program::Blastp, &["-evalue", "1e-5", "-num_threads", "4"])?;
    let inputs = Inputs::new(&query, &subject).named("query.faa", "subject.faa");

    let outcome = run_collect(&search, inputs, &["6", "0"])?;
    std::io::stdout().write_all(&outcome.reports[0].bytes)?;

    for hsp in &outcome.hsps {
        println!(
            "q[{}] vs s[{}]: {} bits (E={:.2e})",
            hsp.q_idx, hsp.s_idx, hsp.bit_score, hsp.e_value
        );
    }
    Ok(())
}
```

### Main `losat::api` Building Blocks

- **Programs (`losat::api::Program`)**: `Blastn`, `Blastp`, `Blastx`, `Tblastn`, `Tblastx`, and `Psiblast` in library search calls (plus `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker`, `run_convert2blastmask` for database and masking workflows).
- **In-memory & streaming searches**:
  - `validate(program, args)` builds a reusable validated `Search`.
  - `run_collect(&search, inputs, formats)` collects report bytes, `Vec<HspRecord>`, per-query `QueryStats`, PSI-BLAST `QueryRounds`, and CLI diagnostics.
  - `run_batches(&search, inputs)` yields one NCBI query batch (`Batch`) at a time so callers can stream results and cancel early.
  - `run(&search, inputs, writers)` and `run_to_file(&search, inputs, format, path)` stream reports directly to writers or disk (required for `-outfmt 13`/`14` side files).
  - `describe(program)` exposes the CLI option schema programmatically.
- **Database and PSSM inputs**: `Inputs::database(&query, "swissprot")` searches formatted BLAST databases (following `BLASTDB` and `.ncbirc` resolution), while `Inputs::pssm` and `Inputs::pssm_database` search saved ASN.1 PSSMs from `psiblast` / `tblastn`.

Internal engine modules outside `losat::api` are `pub(crate)` and hidden from public docs. See the [Rust library guide](https://satoshikawato.github.io/LOSAT/libraries/rust.html) for complete details.

## Cargo Features

| Feature | Default | Description |
| :--- | :--- | :--- |
| `parallel` | Yes | Enables multi-threaded search via `rayon` (`-num_threads`) |
| `wasm-threads` | No | Enables `parallel` for `wasm32-wasip1-threads` reactor/command builds |

## License

Released under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/main/LICENSE).
