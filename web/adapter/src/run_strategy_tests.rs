//! A search strategy file in a run of LOSAT Web (line RG item F2-4): the file of
//! `-import_search_strategy NAME` registered as a kind 3 named input file, and a run of the argv
//! that has no `-query` and no `-subject` takes the strategy's own, against the bytes that NCBI
//! BLAST+ 2.17.0 writes for the same strategy on the command line
//! (`LOSAT/tests/fixtures/strategy/import/`: NCBI's own exports `ss/`, and `bytes/` through
//! `import_cases.tsv`).

use super::*;
use crate::emit::STREAM_DIAGNOSTICS;
use crate::files::{self, KIND_SEARCH_STRATEGY};
use std::collections::BTreeMap;

fn fixtures() -> String {
    concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/fixtures/strategy/import"
    )
    .to_string()
}

/// NCBI's bytes of the file `name` of `bytes/` (`-`: none) in LOSAT's words.
fn frozen(name: &str) -> Vec<u8> {
    if name == "-" {
        return Vec::new();
    }
    let bytes = std::fs::read(format!("{}/bytes/{name}", fixtures())).unwrap();
    losat::web_adapter::losat_bytes(&bytes).into_owned()
}

/// The rows of `import_cases.tsv` named `I.<program>.sub.alone`: the program, the strategy file,
/// and the files of NCBI's standard output and standard error.
fn rows() -> Vec<(String, String, String, String)> {
    let table = std::fs::read_to_string(format!("{}/import_cases.tsv", fixtures())).unwrap();
    table
        .lines()
        .filter_map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            let name = cells.first()?;
            name.ends_with(".sub.alone").then_some(())?;
            let file = cells[2].strip_prefix("-import_search_strategy {ss}/")?;
            Some((
                cells[1].to_string(),
                file.to_string(),
                cells[4].to_string(),
                cells[5].to_string(),
            ))
        })
        .collect()
}

/// Registers `bytes` as the strategy file `name` (chunks of 1,000 bytes) and returns its handle.
fn register_strategy(name: &str, bytes: &[u8]) -> u32 {
    let handle = files::begin(KIND_SEARCH_STRATEGY, name, bytes.len()).unwrap();
    for part in bytes.chunks(1000) {
        files::chunk(handle, part).unwrap();
    }
    let response = files::end(handle).unwrap();
    assert!(
        response.contains("\"kind\":\"search_strategy\""),
        "{response}"
    );
    handle
}

/// The streams that `run` emitted for `words` with the inputs of `register_inputs` (empty bytes:
/// the strategy's own).
fn run_imported(
    words: &[&str],
    query: &[u8],
    subject: &[u8],
) -> Result<BTreeMap<u32, Vec<u8>>, String> {
    let program = Program::parse(words[0]).unwrap();
    let (query_handle, subject_handle, _) = crate::store::register_inputs(words, query, subject)?;
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(
        program,
        query_handle,
        subject_handle,
        words,
        crate::run::data_loaders(words).ok(),
        |queries, query_end, subjects| run(words, queries, query_end, subjects),
    );
    crate::store::release(query_handle).unwrap();
    // A run of `-db` has the subject handle 0, which names no input.
    if subject_handle != 0 {
        crate::store::release(subject_handle).unwrap();
    }
    result??;
    let mut streams = BTreeMap::new();
    crate::emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    Ok(streams)
}

#[test]
fn a_registered_strategy_is_the_command_lines_import() {
    let rows = rows();
    assert_eq!(
        rows.len(),
        5,
        "one strategy of subject sequences per program"
    );
    for (program, file, out, err) in rows {
        let name = format!("ss/{file}");
        let bytes = std::fs::read(format!("{}/ss/{file}", fixtures())).unwrap();
        let handle = register_strategy(&name, &bytes);
        let words = [program.as_str(), "-import_search_strategy", name.as_str()];
        assert!(crate::run::inputs_together(&words), "{program}");
        let streams =
            run_imported(&words, b"", b"").unwrap_or_else(|error| panic!("{program}: {error}"));
        // Stream 0 is outfmt 0, the CLI's default format; stream 3 the CLI's standard error.
        assert_eq!(
            streams.get(&0).cloned().unwrap_or_default(),
            frozen(&out),
            "{program}: the report"
        );
        assert_eq!(
            streams
                .get(&STREAM_DIAGNOSTICS)
                .cloned()
                .unwrap_or_default(),
            frozen(&err),
            "{program}: the diagnostics"
        );
        files::release(handle).unwrap();
    }
}

#[test]
fn a_strategy_that_is_not_registered_is_a_file_that_is_not_accessible() {
    let words = ["blastp", "-import_search_strategy", "missing.ss"];
    let error = crate::run::parse_imported(&words).err().unwrap();
    assert_eq!(
        error,
        "Command line argument error: Argument \"import_search_strategy\". File is not accessible:  `missing.ss'"
    );
    // The standard input is no registered file either.
    let error = crate::run::parse_imported(&["blastp", "-import_search_strategy", "-"])
        .err()
        .unwrap();
    assert!(error.contains("File is not accessible"), "{error}");
}

#[test]
fn a_file_that_is_not_a_strategy_is_the_command_lines_error() {
    let handle = register_strategy("junk.ss", b"not a strategy");
    let error = crate::run::parse_imported(&["blastp", "-import_search_strategy", "junk.ss"])
        .err()
        .unwrap();
    assert!(error.contains("Failed to read search strategy"), "{error}");
    files::release(handle).unwrap();
}

/// blastn runs a TBLASTX strategy as the command line does (row `C.blastn.from.tblastx.sub`:
/// NCBI's output is tblastx's import with blastn's warnings), on the subjects the strategy holds
/// (the browser has no database files).
#[test]
fn a_strategy_of_another_program_runs_as_on_the_command_line() {
    let table = std::fs::read_to_string(format!("{}/import_cases.tsv", fixtures())).unwrap();
    let cells: Vec<&str> = table
        .lines()
        .find(|line| line.starts_with("C.blastn.from.tblastx.sub\t"))
        .expect("the row")
        .split('\t')
        .collect();
    let bytes = std::fs::read(format!("{}/ss/tblastx.sub.ss", fixtures())).unwrap();
    let handle = register_strategy("tblastx.sub.ss", &bytes);
    let words = ["blastn", "-import_search_strategy", "tblastx.sub.ss"];
    let (program, _, _) = crate::run::parse_imported(&words).unwrap();
    assert_eq!(program.name(), "tblastx", "the strategy's program runs");
    let streams = run_imported(&words, b"", b"").unwrap();
    assert_eq!(
        streams.get(&0).cloned().unwrap_or_default(),
        frozen(cells[4])
    );
    assert_eq!(
        streams
            .get(&STREAM_DIAGNOSTICS)
            .cloned()
            .unwrap_or_default(),
        frozen(cells[5])
    );
    files::release(handle).unwrap();
}

#[test]
fn the_files_name_is_a_path_of_the_working_directory() {
    let bytes = std::fs::read(format!("{}/ss/blastp.sub.ss", fixtures())).unwrap();
    let handle = register_strategy("dotted.ss", &bytes);
    for name in ["dotted.ss", "./dotted.ss"] {
        let words = ["blastp", "-import_search_strategy", name];
        let (program, _, supplied) = crate::run::parse_imported(&words).unwrap();
        assert_eq!(program.name(), "blastp");
        let supplied = supplied.expect("a strategy");
        assert!(!supplied.queries.is_empty() && !supplied.subjects.is_empty());
    }
    files::release(handle).unwrap();
}

#[test]
fn the_argvs_own_inputs_override_the_strategys() {
    let bytes = std::fs::read(format!("{}/ss/blastp.sub.ss", fixtures())).unwrap();
    let handle = register_strategy("override.ss", &bytes);
    let words = [
        "blastp",
        "-import_search_strategy",
        "override.ss",
        "-query",
        "q.faa",
    ];
    let (_, _, supplied) = crate::run::parse_imported(&words).unwrap();
    let supplied = supplied.unwrap();
    assert!(supplied.queries.is_empty(), "the argv's -query overrides");
    assert!(!supplied.subjects.is_empty());
    files::release(handle).unwrap();
}

#[test]
fn validate_reads_the_strategy_where_the_command_line_does() {
    let bytes = std::fs::read(format!("{}/ss/blastp.sub.ss", fixtures())).unwrap();
    let handle = register_strategy("validated.ss", &bytes);
    crate::run::validate(&["blastp", "-import_search_strategy", "validated.ss"]).unwrap();
    files::release(handle).unwrap();
    // The file released, the name is a file that is not accessible, as for any later run.
    let error =
        crate::run::validate(&["blastp", "-import_search_strategy", "validated.ss"]).unwrap_err();
    assert!(error.contains("File is not accessible"), "{error}");
}

/// A strategy whose subject is a database is the browser's run of `-db` (line WQA D12a): the
/// *start* field `database` is true, the host gives no subject bytes and the run no subject
/// handle, and the engine opens the database by the strategy's name in the run's working
/// directory, as the command line does (row `I.blastn.db.alone`; the name made absolute here, as
/// the test's working directory is not `LOSAT/`: the report names the database by its title).
#[test]
fn a_strategy_of_a_database_is_a_run_of_db() {
    let table = std::fs::read_to_string(format!("{}/import_cases.tsv", fixtures())).unwrap();
    let cells: Vec<&str> = table
        .lines()
        .find(|line| line.starts_with("I.blastn.db.alone\t"))
        .expect("the row")
        .split('\t')
        .collect();
    let text = std::fs::read_to_string(format!("{}/ss/blastn.db.ss", fixtures())).unwrap();
    let database = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/db_search/arch_snucl"
    );
    assert!(text.contains("\"tests/db_search/arch_snucl\""));
    let text = text.replace("\"tests/db_search/arch_snucl\"", &format!("\"{database}\""));
    let handle = register_strategy("blastn.db.ss", text.as_bytes());
    let words = ["blastn", "-import_search_strategy", "blastn.db.ss"];
    assert!(crate::run::searches_database(&words));
    assert!(crate::run::start_json(&words)
        .unwrap()
        .ends_with(",\"database\":true}"));
    let streams = run_imported(&words, b"", b"").unwrap();
    assert_eq!(
        streams.get(&0).cloned().unwrap_or_default(),
        frozen(cells[4])
    );
    assert_eq!(
        streams
            .get(&STREAM_DIAGNOSTICS)
            .cloned()
            .unwrap_or_default(),
        frozen(cells[5])
    );
    files::release(handle).unwrap();
}
