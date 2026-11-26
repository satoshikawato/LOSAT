//! A search strategy file as bytes (`Search::read_strategy`, `Inputs::strategy`; line RG item
//! F2-4) against the bytes that NCBI BLAST+ 2.17.0 writes for `-import_search_strategy FILE`
//! (`tests/fixtures/strategy/import/import_cases.tsv`, the rows named `.alone`: the strategy
//! files are NCBI's own exports, `ss/`, and each output is stored once, `bytes/`): the report of
//! outfmt 0 and the diagnostics of a library call are the command line's standard output and
//! standard error. The rows `C.<program>.from.<program>` run another program's strategy as the
//! command line does, and a strategy's subjects are searched once at any `-num_threads` (Owner
//! 2026-11-20 F2-Q1 "Q1=C").

use std::path::PathBuf;

use losat::api::{run, run_collect, validate, Error, Inputs, Program};

/// The rows `C.<program>.from.<program>` of `import_cases.tsv` that exit 0 (alone and with the
/// command line's options).
const C_ROWS: usize = 47;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/strategy/import")
}

/// A row of `import_cases.tsv` that imports one strategy file, with the command line's options
/// `options` (none in the rows named `.alone`).
struct Row {
    name: String,
    program: Program,
    file: String,
    options: Vec<String>,
    out: String,
    err: String,
}

/// The rows of `import_cases.tsv` named as `named` says that import one strategy file (its name
/// without spaces) and exit 0.
fn rows_named(named: impl Fn(&str) -> bool) -> Vec<Row> {
    let table = std::fs::read_to_string(fixtures().join("import_cases.tsv")).expect("the table");
    table
        .lines()
        .filter(|line| !(line.starts_with('#') || line.starts_with("name\t") || line.is_empty()))
        .filter_map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            let program = Program::parse(cells[1]).ok()?;
            let mut words = cells[2]
                .strip_prefix("-import_search_strategy {ss}/")?
                .split(' ');
            let file = words.next()?;
            let options: Vec<String> = words.map(|word| word.replace("%20", " ")).collect();
            (named(cells[0]) && cells[3] == "0" && !file.contains("%20")).then(|| Row {
                name: cells[0].to_string(),
                program,
                file: file.to_string(),
                options,
                out: cells[4].to_string(),
                err: cells[5].to_string(),
            })
        })
        .collect()
}

fn rows() -> Vec<Row> {
    rows_named(|name| name.ends_with(".alone"))
        .into_iter()
        .filter(|row| row.options.is_empty())
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

/// Runs each row through the library (a search of the row's options and `argv` that reads the
/// row's strategy file as bytes) and lists the rows whose report or diagnostics differ from
/// NCBI's.
fn differing(rows: &[Row], argv: &[&str]) -> Vec<String> {
    let mut failed = Vec::new();
    for row in rows {
        // The report's format is the output's name in the library API, not an option of the argv.
        let mut format = "0";
        let mut words: Vec<&str> = Vec::new();
        let mut given = row
            .options
            .iter()
            .map(String::as_str)
            .chain(argv.iter().copied());
        while let Some(word) = given.next() {
            if word == "-outfmt" {
                format = given.next().expect("the value of -outfmt");
            } else {
                words.push(word);
            }
        }
        let argv = words.as_slice();
        let search = validate(row.program, argv).unwrap();
        let outcome = (|| -> Result<(Vec<u8>, Vec<u8>), Error> {
            let strategy = search.read_strategy(&strategy_bytes(&row.file))?;
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            run(
                &search,
                Inputs::strategy(&strategy),
                &mut [(format, &mut report)],
                &mut diagnostics,
            )?;
            Ok((report, diagnostics))
        })();
        match outcome {
            Ok((report, diagnostics)) => {
                if report != frozen(&row.out) {
                    failed.push(format!(
                        "{} {argv:?}: report differs: {}",
                        row.name,
                        crate::outfmt11_blastx::first_difference(&report, &frozen(&row.out))
                    ));
                } else if diagnostics != frozen(&row.err) {
                    failed.push(format!(
                        "{} {argv:?}: diagnostics {:?}",
                        row.name,
                        String::from_utf8_lossy(&diagnostics)
                    ));
                }
            }
            Err(error) => failed.push(format!("{} {argv:?}: {error}", row.name)),
        }
    }
    failed
}

#[test]
fn a_strategy_as_bytes_is_the_command_lines_import() {
    let rows = rows();
    assert!(rows.len() >= 100, "{} rows", rows.len());
    let failed = differing(&rows, &[]);
    assert!(
        failed.is_empty(),
        "{} of {}: {failed:#?}",
        failed.len(),
        rows.len()
    );
}

/// blastn and blastx run a TBLASTX strategy, tblastx a BLASTX or a BLASTN one, as the command
/// line does (NCBI's output is the strategy's program's import with the running program's
/// warnings, those of the options the command line gives among them).
#[test]
fn a_strategy_of_another_program_runs_as_on_the_command_line() {
    let rows = rows_named(|name| name.starts_with("C.") && name.contains(".from."));
    assert_eq!(rows.len(), C_ROWS, "the rows C.*.from.*");
    assert!(rows.iter().any(|row| !row.options.is_empty()));
    let failed = differing(&rows, &[]);
    assert!(failed.is_empty(), "{failed:#?}");
}

/// A strategy's subjects are searched once at `-num_threads 4` (Owner 2026-11-20 F2-Q1 "Q1=C"),
/// as on the command line: NCBI's `-num_threads 1` bytes.
#[test]
fn a_strategys_subjects_are_searched_once_with_threads() {
    let rows = rows_named(|name| name.starts_with("I.") && name.ends_with(".sub.alone"));
    assert_eq!(rows.len(), 5, "the rows I.<program>.sub.alone");
    let failed = differing(&rows, &["-num_threads", "4"]);
    assert!(failed.is_empty(), "{failed:#?}");
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
