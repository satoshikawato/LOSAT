# Rust (`losat`)

The crate `losat` holds both the `losat` command line and a library. The library is the module `losat::api`: it searches FASTA records in memory or a BLAST database with the command line's options, and returns the command line's reports byte for byte, plus one typed record per HSP and the statistics of each query.

`losat::api` is the crate's public surface. Every other module is the command line or the engine: hidden from the documentation and outside the compatibility promise. The crate needs Rust 1.92 or later.

## Adding `losat` to your project

The crate is on [crates.io](https://crates.io/crates/losat). Add it with `cargo add losat`, or in `Cargo.toml`:

```toml
[dependencies]
losat = "0.2"
# or, with a clone of the repository next to your project:
# losat = { path = "../LOSAT/LOSAT" }
```

`cargo doc --no-deps --open` in the repository's `LOSAT/` directory shows the API.

### Features

| Feature | Default | What it does |
|:---|:---|:---|
| `parallel` | on | Searches with more than one thread (`-num_threads`), through Rayon |
| `wasm-threads` | off | Turns on `parallel` for the `wasm32-wasip1-threads` target |

Without `parallel` (`default-features = false`), the build has no Rayon and searches on one thread, as the serial WebAssembly builds do. The crate's other features are for LOSAT's own diagnostics.

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
    std::io::stdout().write_all(&outcome.reports[0].bytes)?; // the -outfmt 6 report

    for hsp in &outcome.hsps {
        println!("query {} subject {}: {} bits, E = {:e}", hsp.q_idx, hsp.s_idx, hsp.bit_score, hsp.e_value);
    }
    Ok(())
}
```

The repository has a program that does this for any search program, [`LOSAT/examples/library_search.rs`](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LOSAT/examples/library_search.rs). From the `LOSAT/` directory:

```bash
cargo run --release --example library_search -- blastp query.faa subject.faa -evalue 1e-5
```

It prints the `-outfmt 6` report.

## What `losat::api` offers

- **Programs.** `Program::Blastn`, `Blastp`, `Blastx`, `Tblastn`, `Tblastx` and `Psiblast`.
- **Options.** `validate(program, &argv)` takes the command line's options, without `-query`, `-subject`, `-db`, `-in_pssm`, `-in_msa`, `-out` and `-outfmt`, checks them as the command line does (with NCBI's argument reader), and returns a `Search` you can reuse. `describe(program)` lists the options with their kinds, defaults and choices.
- **Inputs.** `read_fasta(program, role, bytes)` reads FASTA bytes with the command line's reader (`Role::Query` or `Role::Subject`). `Inputs::new(&query, &subject)` pairs the records; `.named(query_name, subject_name)` sets the names that the reports show where the command line shows its file names. `Inputs::database(&query, "swissprot")` searches a BLAST database, found as the command line finds `-db` (`BLASTDB`, `.ncbirc`). `Inputs::pssm` and `Inputs::pssm_database` take the PSSM file of `-in_pssm` (`psiblast` and `tblastn`). `Inputs::msa` and `Inputs::msa_database` take the alignment file of `psiblast`'s `-in_msa`: the search reads it as the command line does, and its one query is the alignment's master (query 0). `-msa_master_idx` and `-ignore_msa_master` are options of the argv (`validate` accepts them; the command line takes them with `-in_msa` only, so a search of another query fails with the command line's argument error).
- **Searches.**
  - `run_collect(&search, inputs, &formats)` returns an `Outcome`: the reports (`reports[i].bytes`, and `reports[i].spans`, the byte range of each HSP's row or section), the HSP records (`hsps`), the statistics of each query (`query_stats`: effective search space, length adjustment, lambda, K and H), the `psiblast` rounds (`rounds`) and the diagnostics.
  - `run_batches(&search, inputs)` gives the same results one of NCBI's query batches at a time. Your program keeps one batch in memory, reports progress between batches, and stops the search by not asking for the next.
  - `run(&search, inputs, &mut reports, &mut diagnostics)` writes each report to a writer you give it, as the search makes it.
  - `run_to_file(&search, inputs, outfmt, path, &mut diagnostics)` is the command line's `-out`, and the way to `-outfmt 13` and `14` with their side files.
- **Tools.** `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker`, `run_windowmasker`, `run_convert2blastmask` and `run_makembindex` run the tools with your reader as their standard input and your writers as their output. A help word (`-h`, `-help`) writes the tool's help to your output and returns `Ok(())`, as the command line prints it and exits with status 0; a search has no output for it, so `validate` returns the help as `Error::Argv`.
- **Errors and output.** Errors (`losat::api::Error`) carry the command line's message for the same failure as bytes (`message()`), also bytes that are not UTF-8, by the rules of `Error`'s documentation (an argument error (`Error::Argv`) is exactly what the command line writes to its standard error: NCBI's short usage, a line of 72 `=`, an empty line and NCBI's two error lines; with `-logfile`, a call's error is the usage, NCBI's line `Error in command-line arguments. See error logs for more details.` and the line of `=`, and the two error lines go to the file of `-logfile`, as on the command line, while `validate`, which opens no file, gives the error as the command line writes it without `-logfile`; a search's message of NCBI's keeps its newlines; an input's error, a tool's and that of a start that ends the program are the command line's text without its last newlines; LOSAT's own refusals have no `Error: ` and no last newline; text that a search writes before the failure is not in it); `to_string_lossy()` and `Display` give the text with one U+FFFD per maximal ill-formed subpart. The library writes nothing to standard output or standard error; what the command line writes to its standard error comes back as diagnostics, as bytes (`Outcome::diagnostics`, `Batch::diagnostics`).

The numbers of an HSP record are the engine's, not rounded: each field that names an `-outfmt` column holds that column's value for the HSP's row. The ids (`query_id`, `subject_id`) are an `Id`: the bytes of the `qaccver` and `saccver` columns, also those that are not UTF-8. Where NCBI writes no `saccver` for a subject (its lookup fails), `subject_id` is the id by the query's rules; for a database sequence, the first word of its title. `Id` has no `Display`: `to_str()` gives the text of an id that is UTF-8, and `to_string_lossy()` text with U+FFFD. Read text to show from a report.

## Differences from the command line

- The library does not write a search strategy file: it refuses `-export_search_strategy`. The command line writes it.
- Records that you build yourself with `FastaRecord::new` carry no definition line as read. In `-outfmt 11`, their definition line is `>` and the title. Records from `read_fasta` keep it, as the command line does (see [Differences from NCBI BLAST+](../differences.md#the-library)).
