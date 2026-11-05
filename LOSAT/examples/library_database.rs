//! A database search through the library API, with the arguments of the command line.
//!
//! ```text
//! cargo run --release --example library_database -- blastn -db nt -query query.fna -outfmt 6
//! ```
//!
//! Reads the query file into memory, searches the database of `-db` with `LOSAT::api::run` and
//! writes the report to standard output and the search's warnings and a failure's message to
//! standard error, so that a run compares with `losat <program>` with the same arguments
//! (`tests/db_search_fixtures.py check --library`). A failure ends with status 1, whatever the
//! command line's status for it.

use std::io::Write as _;
use std::path::Path;
use std::process::ExitCode;

use LOSAT::api::{read_input_file, read_inputs, run, validate, Error, Inputs, Program, Role};

fn main() -> ExitCode {
    match search() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let mut stderr = std::io::stderr();
            let _ = stderr.write_all(error.message().as_bytes());
            if !error.message().ends_with('\n') {
                let _ = stderr.write_all(b"\n");
            }
            ExitCode::FAILURE
        }
    }
}

fn search() -> Result<(), Error> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((program, words)) = args.split_first() else {
        return Err(Error::Argv(
            "usage: library_database <program> -db <db> -query <file> [options...]".to_string(),
        ));
    };
    let program: Program = program.parse()?;
    let (mut db, mut query, mut outfmt) = (None, None, "0");
    let mut options: Vec<&str> = Vec::new();
    let mut words = words.iter().map(String::as_str);
    while let Some(word) = words.next() {
        match word {
            "-db" => db = words.next(),
            "-query" => query = words.next(),
            "-outfmt" => outfmt = words.next().unwrap_or(outfmt),
            option => options.push(option),
        }
    }
    let (Some(db), Some(query)) = (db, query) else {
        return Err(Error::Argv("-db and -query are required".to_string()));
    };
    let search = validate(program, &options)?;
    let bytes = read_input_file(Role::Query, Path::new(query))?;
    let (records, _) = read_inputs(&search, &bytes, b"")?;
    let inputs = Inputs::database(&records, db).named(query, db);
    let mut stdout = std::io::stdout();
    let mut stderr = std::io::stderr();
    run(&search, inputs, &mut [(outfmt, &mut stdout)], &mut stderr)
}
