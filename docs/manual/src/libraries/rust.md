# Rust (`losat`)

The `losat` crate provides both the `losat` binary and an embeddable Rust library under [`losat::api`](https://docs.rs/losat).

Every internal module outside `losat::api` is `pub(crate)` and outside the semver compatibility promise. Minimum supported Rust version (MSRV) is **1.92**.

## Adding `losat` to your project

The `losat` crate is scheduled for publication on crates.io immediately following the v0.2.0 release:

```toml
[dependencies]
losat = "0.2"
```

Until crates.io publication completes (or to track git), depend on the repository subdirectory directly:

```toml
[dependencies]
losat = { git = "https://github.com/satoshikawato/LOSAT.git" }
```

Disable the default `parallel` feature if you want a single-threaded build without `rayon`:

```toml
[dependencies]
losat = { version = "0.2", default-features = false }
```

| Feature | Default | Purpose |
| :--- | :--- | :--- |
| `parallel` | Yes | Multi-threaded search execution via `rayon` (`-num_threads`) |
| `wasm-threads` | No | Enables `parallel` for `wasm32-wasip1-threads` targets |

## Quick start

```rust
use std::io::Write as _;
use losat::api::{read_fasta, run_collect, validate, Inputs, Program, Role};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let query = read_fasta(Program::Blastp, Role::Query, &std::fs::read("query.faa")?)?;
    let subject = read_fasta(Program::Blastp, Role::Subject, &std::fs::read("subject.faa")?)?;

    let search = validate(Program::Blastp, &["-evalue", "1e-5", "-num_threads", "4"])?;
    let inputs = Inputs::new(&query, &subject).named("query.faa", "subject.faa");

    let outcome = run_collect(&search, inputs, &["6", "0"])?;
    std::io::stdout().write_all(&outcome.reports[0].bytes)?;

    println!("Matched {} HSPs", outcome.hsps.len());
    Ok(())
}
```

A complete runnable example lives at [`LOSAT/examples/library_search.rs`](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LOSAT/examples/library_search.rs):

```bash
cd LOSAT
cargo run --release --example library_search -- blastp query.faa subject.faa -evalue 1e-5
```

## Core `losat::api` types and functions

- **Supported programs (`Program`)**: `Blastn`, `Blastp`, `Blastx`, `Tblastn`, `Tblastx`, and `Psiblast`.
- **Reading and validating inputs**:
  - `read_fasta(program, role, bytes)` parses FASTA bytes with the CLI's reader (`Role::Query` or `Role::Subject`).
  - `validate(program, options)` validates CLI options (all search flags except `-query`, `-subject`, `-db`, `-in_pssm`, `-out`, `-outfmt`) into a reusable `Search`.
  - `describe(program)` inspects available options, value types, defaults, and choices.
  - `Inputs::new(&query, &subject)` pairs in-memory FASTA records; `.named(qname, sname)` sets the labels printed in reports.
  - `Inputs::database(&query, "swissprot")` searches formatted BLAST databases (`-db` resolution via cwd / `BLASTDB` / `.ncbirc`).
  - `Inputs::pssm` and `Inputs::pssm_database` load `-in_pssm` checkpoints for `psiblast` and `tblastn`.
- **Executing searches**:
  - `run_collect(&search, inputs, formats)` collects `reports` (`Vec<Report>` with byte slices and span offsets), `hsps` (`Vec<HspRecord>`), `query_stats`, `rounds`, and `diagnostics`.
  - `run_batches(&search, inputs)` yields an iterator of per-batch `Batch` items for bounded-memory processing or early stopping.
  - `run(&search, inputs, writers)` streams formatted reports directly into caller-provided `std::io::Write` destinations.
  - `run_to_file(&search, inputs, format, path)` writes `-out path` reports (and companion files for `-outfmt 13`/`14`).
- **Database and masking tools**: `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker`, and `run_convert2blastmask`.

The library writes nothing to `stdout` or `stderr` by itself: CLI diagnostics and warnings return in `Outcome::diagnostics` / `Batch::diagnostics`. Note one edge difference: `-export_search_strategy` requires the CLI binary, and when passing pre-parsed records without retained raw definition lines for `-outfmt 11`, titles normalize from `>` + title (see [Differences from NCBI BLAST+](../differences.md#the-library)).
