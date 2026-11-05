//! A search through the library API, from two FASTA files read into memory.
//!
//! ```text
//! cargo run --release --example library_search -- blastp query.faa subject.faa -evalue 1e-5
//! ```
//!
//! Prints the outfmt 6 report (the bytes `losat blastp -query query.faa -subject subject.faa
//! -outfmt 6 -evalue 1e-5` writes) to standard output and, with `LIBRARY_SEARCH_RECORDS=1`, one
//! line per HSP record after it. Writes the search's warnings (what the command line writes to
//! its standard error) and a failure's message to standard error, and nothing else: the library
//! itself writes nothing there (`tests/library_api.rs` checks it).

use std::io::Write as _;

use LOSAT::api::{read_fasta, run_collect, validate, Error, Inputs, Program, Role};

fn main() {
    if let Err(error) = search() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

fn search() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [program, query_path, subject_path, options @ ..] = args.as_slice() else {
        return Err("usage: library_search <program> <query.fa> <subject.fa> [options...]".into());
    };
    let program: Program = program.parse().map_err(|error: Error| error.to_string())?;
    let read = |path: &str, role: Role| -> Result<_, String> {
        let bytes = std::fs::read(path).map_err(|error| format!("{path}: {error}"))?;
        read_fasta(program, role, &bytes).map_err(|error| error.to_string())
    };
    let query = read(query_path, Role::Query)?;
    let subject = read(subject_path, Role::Subject)?;
    let options: Vec<&str> = options.iter().map(String::as_str).collect();
    let search = validate(program, &options).map_err(|error| error.to_string())?;
    // The reports name the inputs as the command line names its files.
    let inputs = Inputs::new(&query, &subject).named(query_path, subject_path);
    let outcome = run_collect(&search, inputs, &["6"]).map_err(|error| error.to_string())?;
    std::io::stderr()
        .write_all(&outcome.diagnostics)
        .map_err(|error| error.to_string())?;
    let mut stdout = std::io::stdout().lock();
    stdout
        .write_all(&outcome.reports[0].bytes)
        .map_err(|error| error.to_string())?;
    if std::env::var_os("LIBRARY_SEARCH_RECORDS").is_some() {
        for hsp in &outcome.hsps {
            writeln!(
                stdout,
                "# hsp {} q{} s{} rank {} raw {} bits {} e {} q {}-{} s {}-{} row {:?}",
                hsp.index,
                hsp.q_idx,
                hsp.s_idx,
                hsp.rank,
                hsp.raw_score,
                hsp.bit_score,
                hsp.e_value,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
                outcome.reports[0].spans[hsp.index]
            )
            .map_err(|error| error.to_string())?;
        }
    }
    Ok(())
}
