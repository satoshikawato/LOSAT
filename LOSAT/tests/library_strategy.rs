//! A search strategy file as bytes (`Search::read_strategy`, `Inputs::strategy`; line RG item
//! F2-4) against the bytes that NCBI BLAST+ 2.17.0 writes for `-import_search_strategy FILE`
//! (`tests/fixtures/strategy/import/import_cases.tsv`, the rows named `.alone`: the strategy
//! files are NCBI's own exports, `ss/`, and each output is stored once, `bytes/`): the report of
//! outfmt 0 and the diagnostics of a library call are the command line's standard output and
//! standard error.

use std::path::PathBuf;

use losat::api::{run, run_collect, validate, Error, Inputs, Program};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/strategy/import")
}

/// A row of `import_cases.tsv` that imports one strategy file without any other option.
struct Row {
    name: String,
    program: Program,
    file: String,
    out: String,
    err: String,
}

fn rows() -> Vec<Row> {
    let table = std::fs::read_to_string(fixtures().join("import_cases.tsv")).expect("the table");
    table
        .lines()
        .filter(|line| !(line.starts_with('#') || line.starts_with("name\t") || line.is_empty()))
        .filter_map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            let program = Program::parse(cells[1]).ok()?;
            let file = cells[2].strip_prefix("-import_search_strategy {ss}/")?;
            (cells[0].ends_with(".alone") && cells[3] == "0" && !file.contains(' ')).then(|| Row {
                name: cells[0].to_string(),
                program,
                file: file.to_string(),
                out: cells[4].to_string(),
                err: cells[5].to_string(),
            })
        })
        .collect()
}

/// NCBI's bytes of the file `name` of `bytes/` (`-`: none) in LOSAT's words.
fn frozen(name: &str) -> Vec<u8> {
    if name == "-" {
        return Vec::new();
    }
    let bytes = std::fs::read(fixtures().join("bytes").join(name)).expect("frozen bytes");
    losat::wording::losat_output_bytes(&bytes).into_owned()
}

fn strategy_bytes(file: &str) -> Vec<u8> {
    std::fs::read(fixtures().join("ss").join(file)).expect("a strategy file")
}

#[test]
fn a_strategy_as_bytes_is_the_command_lines_import() {
    let rows = rows();
    assert!(rows.len() >= 100, "{} rows", rows.len());
    let mut failed = Vec::new();
    for row in &rows {
        let search = validate(row.program, &[]).unwrap();
        let outcome = (|| -> Result<(Vec<u8>, Vec<u8>), Error> {
            let strategy = search.read_strategy(&strategy_bytes(&row.file))?;
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            run(
                &search,
                Inputs::strategy(&strategy),
                &mut [("0", &mut report)],
                &mut diagnostics,
            )?;
            Ok((report, diagnostics))
        })();
        match outcome {
            Ok((report, diagnostics)) => {
                if report != frozen(&row.out) {
                    failed.push(format!(
                        "{}: report differs: {}",
                        row.name,
                        crate::outfmt11_blastx::first_difference(&report, &frozen(&row.out))
                    ));
                } else if diagnostics != frozen(&row.err) {
                    failed.push(format!(
                        "{}: diagnostics {:?}",
                        row.name,
                        String::from_utf8_lossy(&diagnostics)
                    ));
                }
            }
            Err(error) => failed.push(format!("{}: {error}", row.name)),
        }
    }
    assert!(
        failed.is_empty(),
        "{} of {}: {failed:#?}",
        failed.len(),
        rows.len()
    );
}

#[test]
fn the_argv_cannot_import_a_file() {
    for program in [Program::Blastn, Program::Blastp, Program::Tblastx] {
        let error = validate(program, &["-import_search_strategy", "s.asn"]).unwrap_err();
        assert_eq!(
            error,
            Error::Argv(
                "-import_search_strategy is not accepted by the library API: a search strategy is given as bytes (Search::read_strategy, Inputs::strategy)"
                    .into()
            )
        );
    }
}

#[test]
fn a_file_that_is_not_a_strategy_is_the_command_lines_error() {
    let search = validate(Program::Blastp, &[]).unwrap();
    let error = search.read_strategy(b"not a strategy").unwrap_err();
    assert!(
        matches!(&error, Error::Input(message) if String::from_utf8_lossy(message).contains("Failed to read search strategy")),
        "{error:?}"
    );
}

#[test]
fn a_strategy_belongs_to_the_search_that_read_it() {
    let bytes = strategy_bytes("blastp.sub.ss");
    let search = validate(Program::Blastp, &[]).unwrap();
    let strategy = search.read_strategy(&bytes).unwrap();
    let other = validate(Program::Blastp, &["-line_length", "40"]).unwrap();
    let error = run_collect(&other, Inputs::strategy(&strategy), &["6"]).unwrap_err();
    assert!(
        error.to_string_lossy().contains("read for another search"),
        "{error}"
    );
}
