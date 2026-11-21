//! A search through the library API, with the arguments of the command line.
//!
//! ```text
//! cargo run --release --example library_command_line -- blastn -db nt -query query.fna -outfmt 6
//! ```
//!
//! Reads the query and subject files into memory (not for `-dryrun`, as the command line; a PSSM
//! query, `-in_pssm`, is the file's name, `Inputs::pssm`; so is psiblast's alignment query,
//! `-in_msa`, `Inputs::msa`), searches with `losat::api::run`, or
//! `run_to_file` with `-out`, and writes the search's warnings and a failure's message to standard
//! error, so that a run compares with `losat <program>` with the same arguments
//! (`tests/db_search_fixtures.py check --library`). `makeblastdb`, `blastdbcmd`, `dustmasker`,
//! `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` go to their functions, with
//! the standard streams of this process
//! (`tests/masker_fixtures.rs`). A failure ends with status 1, whatever the command line's
//! status for it.

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use losat::api::{
    read_input_file, read_inputs, run, run_blastdbcmd, run_convert2blastmask, run_dustmasker,
    run_makeblastdb, run_makembindex, run_segmasker, run_to_file, run_windowmasker, validate,
    Error, Inputs, Program, Role,
};

fn main() -> ExitCode {
    match search() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let mut stderr = std::io::stderr();
            let _ = stderr.write_all(error.message());
            if !error.message().is_empty() && !error.message().ends_with(b"\n") {
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
        Err(Error::Search(_)) => Err(Error::Search(Vec::new())),
        other => other,
    }
}

fn search() -> Result<(), Error> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((program, words)) = args.split_first() else {
        return Err(Error::Argv(
            "usage: library_command_line <program> -query <file> -db <db> [options...]".into(),
        ));
    };
    // The tools take the command line's words as they are.
    let words_of_tool: Vec<&str> = words.iter().map(String::as_str).collect();
    let tool = match program.as_str() {
        "blastdbcmd" => run_blastdbcmd,
        "makeblastdb" => run_makeblastdb,
        "dustmasker" => run_dustmasker,
        "segmasker" => run_segmasker,
        "windowmasker" => run_windowmasker,
        "convert2blastmask" => run_convert2blastmask,
        "makembindex" => run_makembindex,
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
    let (mut pssm, mut msa) = (None, None);
    let mut options: Vec<&str> = Vec::new();
    let mut words = words.iter().map(String::as_str);
    while let Some(word) = words.next() {
        match word {
            "-db" => db = words.next(),
            "-subject" => subject = words.next(),
            "-query" => query = words.next(),
            "-in_pssm" => pssm = words.next(),
            "-in_msa" => msa = words.next(),
            "-out" => out = words.next(),
            "-outfmt" => outfmt = words.next().unwrap_or(outfmt),
            option => options.push(option),
        }
    }
    let search = validate(program, &options)?;
    // A dry run reads no input, as the command line's does not (`Search::dry_run`).
    let dry_run = search.dry_run();
    // The query is records, a PSSM file or an alignment file (`Inputs::new`, `Inputs::pssm`,
    // `Inputs::msa`), never two of them.
    let query_file = match (pssm, msa) {
        (Some(_), Some(_)) => {
            return Err(Error::Argv(
                "-in_pssm and -in_msa are two queries: give one".into(),
            ))
        }
        (Some(file), None) => Some((file, true)),
        (None, Some(file)) => Some((file, false)),
        (None, None) => None,
    };
    let query = match (query, query_file) {
        (Some(_), Some(_)) => {
            return Err(Error::Argv(
                "-query and -in_pssm or -in_msa are two queries: give one".into(),
            ))
        }
        (Some(query), None) => query,
        (None, Some(_)) => "",
        (None, None) if dry_run => "",
        (None, None) => return Err(Error::Argv("-query is required".into())),
    };
    let read = |role, name: &str| match dry_run {
        true => Ok(Vec::new()),
        false => read_input_file(role, Path::new(name)),
    };
    let query_bytes = match query_file {
        Some(_) => Vec::new(),
        None => read(Role::Query, query)?,
    };
    let subject_bytes = match subject {
        Some(subject) => read(Role::Subject, subject)?,
        None => Vec::new(),
    };
    let (query_records, subject_records) = read_inputs(&search, &query_bytes, &subject_bytes)?;
    let inputs = match (db, subject, query_file) {
        (Some(db), _, Some((file, true))) => Inputs::pssm_database(Path::new(file), db),
        (Some(db), _, Some((file, false))) => Inputs::msa_database(Path::new(file), db),
        (Some(db), _, None) => Inputs::database(&query_records, db).named(query, db),
        (None, Some(subject), Some((file, true))) => {
            Inputs::pssm(Path::new(file), &subject_records).named(query, subject)
        }
        (None, Some(subject), Some((file, false))) => {
            Inputs::msa(Path::new(file), &subject_records).named(query, subject)
        }
        (None, Some(subject), None) => {
            Inputs::new(&query_records, &subject_records).named(query, subject)
        }
        (None, None, Some((file, true))) if dry_run => Inputs::pssm(Path::new(file), &[]),
        (None, None, Some((file, false))) if dry_run => Inputs::msa(Path::new(file), &[]),
        (None, None, None) if dry_run => Inputs::new(&[], &[]),
        (None, None, _) => {
            return Err(Error::Argv("-db or -subject is required".into()));
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
