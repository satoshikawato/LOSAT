//! A search through the library API, with the arguments of the command line.
//!
//! ```text
//! cargo run --release --example library_command_line -- blastn -db nt -query query.fna -outfmt 6
//! ```
//!
//! Reads the query and subject files into memory, searches with `losat::api::run`, or
//! `run_to_file` with `-out`, and writes the search's warnings and a failure's message to
//! standard error, so that a run compares with `losat <program>` with the same arguments
//! (`tests/db_search_fixtures.py check --library`). `makeblastdb`, `blastdbcmd`, `dustmasker`,
//! `segmasker` and `convert2blastmask` go to their functions, with the standard streams of this process
//! (`tests/masker_fixtures.rs`). A failure ends with status 1, whatever the command line's
//! status for it.

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use losat::api::{
    read_input_file, read_inputs, run, run_blastdbcmd, run_convert2blastmask, run_dustmasker,
    run_makeblastdb, run_segmasker, run_to_file, validate, Error, Inputs, Program, Role,
};

fn main() -> ExitCode {
    match search() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let mut stderr = std::io::stderr();
            let _ = stderr.write_all(error.message().as_bytes());
            if !error.message().is_empty() && !error.message().ends_with('\n') {
                let _ = stderr.write_all(b"\n");
            }
            ExitCode::FAILURE
        }
    }
}

/// A tool has written the failure to standard error, so `main` adds nothing; an error of its
/// argv is the parser's.
fn reported(result: Result<(), Error>) -> Result<(), Error> {
    match result {
        Err(Error::Search(_)) => Err(Error::Search(String::new())),
        other => other,
    }
}

fn search() -> Result<(), Error> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((program, words)) = args.split_first() else {
        return Err(Error::Argv(
            "usage: library_command_line <program> -query <file> -db <db> [options...]".to_string(),
        ));
    };
    // The tools take the command line's words as they are.
    let words_of_tool: Vec<&str> = words.iter().map(String::as_str).collect();
    let tool = match program.as_str() {
        "blastdbcmd" => run_blastdbcmd,
        "makeblastdb" => run_makeblastdb,
        "dustmasker" => run_dustmasker,
        "segmasker" => run_segmasker,
        "convert2blastmask" => run_convert2blastmask,
        _ => return search_program(program, words),
    };
    let (mut stdin, mut stdout, mut stderr) =
        (std::io::stdin(), std::io::stdout(), std::io::stderr());
    let result = tool(&words_of_tool, &mut stdin, &mut stdout, &mut stderr);
    let _ = stdout.flush();
    reported(result)
}

/// A search of `program` with the command line's words.
fn search_program(program: &str, words: &[String]) -> Result<(), Error> {
    let program: Program = program.parse()?;
    let (mut db, mut subject, mut query, mut out, mut outfmt) = (None, None, None, None, "0");
    let mut options: Vec<&str> = Vec::new();
    let mut words = words.iter().map(String::as_str);
    while let Some(word) = words.next() {
        match word {
            "-db" => db = words.next(),
            "-subject" => subject = words.next(),
            "-query" => query = words.next(),
            "-out" => out = words.next(),
            "-outfmt" => outfmt = words.next().unwrap_or(outfmt),
            option => options.push(option),
        }
    }
    let query = query.ok_or_else(|| Error::Argv("-query is required".to_string()))?;
    let search = validate(program, &options)?;
    let query_bytes = read_input_file(Role::Query, Path::new(query))?;
    let subject_bytes = match subject {
        Some(subject) => read_input_file(Role::Subject, Path::new(subject))?,
        None => Vec::new(),
    };
    let (query_records, subject_records) = read_inputs(&search, &query_bytes, &subject_bytes)?;
    let inputs = match (db, subject) {
        (Some(db), _) => Inputs::database(&query_records, db).named(query, db),
        (None, Some(subject)) => {
            Inputs::new(&query_records, &subject_records).named(query, subject)
        }
        (None, None) => {
            return Err(Error::Argv("-db or -subject is required".to_string()));
        }
    };
    let mut stderr = std::io::stderr();
    match out {
        Some(path) => run_to_file(&search, inputs, outfmt, Path::new(path), &mut stderr),
        None => {
            let mut stdout = std::io::stdout();
            run(&search, inputs, &mut [(outfmt, &mut stdout)], &mut stderr)
        }
    }
}
