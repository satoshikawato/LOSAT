//! `losat dustmasker`, `losat segmasker` and `losat convert2blastmask` against NCBI BLAST+
//! 2.17.0 dustmasker 1.0.0, segmasker 1.0.0 and convert2blastmask: each row of
//! `tests/<program>/cases.tsv` (its arguments as shell words, run in `tests/<program>`; `< file`
//! is the standard input) gives NCBI's standard output, standard error and exit status, frozen
//! in `tests/<program>/ncbi/<name>.{out,err,rc}` (probes of line RG, items K1, K2 and K3,
//! 2026-11-12 and 11-13, on ws and the local official release). A name the table repeats is
//! compared on its last row, whose output the probe kept. NCBI is not run here. Every row runs
//! through the library API too (`losat::api::run_dustmasker`, `run_segmasker`,
//! `run_convert2blastmask`, by `examples/library_command_line`), with the same bytes.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::library_examples::example_binary;

/// The fixtures of `program`.
fn fixtures(program: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(program)
}

/// One row of the table: its name, the words after the program's name, and the file of
/// `< file`.
struct Case {
    name: String,
    words: Vec<String>,
    stdin: Option<String>,
}

/// The words of `line` as the probe's shell split them: blanks between words, double quotes
/// around a word with blanks or none, and `<` before the file of the standard input (no other
/// shell syntax is in the table).
fn shell_words(line: &str) -> (Vec<String>, Option<String>) {
    let mut words = Vec::new();
    let mut word: Option<String> = None;
    // The open quote: inside `"..."` a `'` is a letter, inside `'...'` a `"`.
    let mut quote: Option<char> = None;
    for c in line.chars() {
        match (quote, c) {
            (None, '"' | '\'') => {
                quote = Some(c);
                word.get_or_insert_with(String::new);
            }
            (Some(open), _) if c == open => quote = None,
            (None, ' ' | '\t') => {
                if let Some(done) = word.take() {
                    words.push(done);
                }
            }
            _ => word.get_or_insert_with(String::new).push(c),
        }
    }
    assert!(quote.is_none(), "unbalanced quotes: {line}");
    words.extend(word);
    let mut stdin = None;
    if let Some(at) = words.iter().position(|word| word == "<") {
        stdin = Some(words[at + 1].clone());
        words.drain(at..at + 2);
    }
    (words, stdin)
}

/// The rows of `program`'s table, the last of each name.
fn cases(program: &str) -> Vec<Case> {
    let table = fs::read_to_string(fixtures(program).join("cases.tsv")).unwrap();
    let mut cases: Vec<Case> = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (name, arguments) = line.split_once('\t').unwrap();
        let (words, stdin) = shell_words(arguments);
        cases.retain(|case| case.name != name);
        cases.push(Case {
            name: name.to_string(),
            words,
            stdin,
        });
    }
    cases
}

/// One run of `losat <program>` for `case` in `tests/<program>`: stdout, stderr, exit status.
/// No LOSAT, NCBI or BLAST variable, no registry file (`NCBI_CONFIG_PATH` names the fixture
/// directory, which has none).
fn run(program: &str, case: &Case) -> (Vec<u8>, Vec<u8>, i32) {
    run_binary(Path::new(env!("CARGO_BIN_EXE_losat")), program, case)
}

/// `run` of `binary <program>` (the command line, or the library's example).
fn run_binary(binary: &Path, program: &str, case: &Case) -> (Vec<u8>, Vec<u8>, i32) {
    let dir = fixtures(program);
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    let stdin = match &case.stdin {
        Some(file) => Stdio::from(fs::File::open(dir.join(file)).unwrap()),
        None => Stdio::null(),
    };
    let output = command
        .env("NCBI_CONFIG_PATH", &dir)
        .current_dir(&dir)
        .arg(program)
        .args(&case.words)
        .stdin(stdin)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// NCBI's stdout, stderr and exit status for the row `name` of `program`.
fn expected(program: &str, name: &str) -> (Vec<u8>, Vec<u8>, i32) {
    let file = |extension: &str| {
        fixtures(program)
            .join("ncbi")
            .join(format!("{name}.{extension}"))
    };
    let status = fs::read_to_string(file("rc")).unwrap();
    // The database rows name the probe's directory in NCBI's search path; here it is this one.
    let here = |bytes: Vec<u8>| -> Vec<u8> {
        let dir = fixtures(program);
        let dir = dir.to_string_lossy();
        String::from_utf8_lossy(&bytes)
            .replace("@CWD@", &dir)
            .into_bytes()
    };
    let err = fs::read(file("err")).unwrap();
    let err = if err.windows(5).any(|w| w == b"@CWD@") {
        here(err)
    } else {
        err
    };
    (
        fs::read(file("out")).unwrap(),
        err,
        status.trim().parse().unwrap(),
    )
}

/// The rows of `program`'s table whose run differs from NCBI's.
fn rows_that_differ(program: &str) -> Vec<String> {
    cases(program)
        .iter()
        .filter(|case| run(program, case) != expected(program, &case.name))
        .map(|case| case.name.clone())
        .collect()
}

#[test]
fn dustmasker_gives_ncbis_bytes() {
    let differ = rows_that_differ("dustmasker");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn segmasker_gives_ncbis_bytes() {
    let differ = rows_that_differ("segmasker");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn convert2blastmask_gives_ncbis_bytes() {
    let differ = rows_that_differ("convert2blastmask");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

// `run_dustmasker` and `run_segmasker` are the command line's run without its application start,
// which reads the file of `-conffile` and ends `-dryrun` (`ncbi_app::run_with`): the library
// refuses those two words, as `run_blastdbcmd` does. The example ends every failure with exit
// status 1. The rows of `program`'s table whose run through the library differs.
fn rows_that_differ_through_the_library(program: &str, refusal_name: &str) -> Vec<String> {
    let example = example_binary("library_command_line");
    let mut differ = Vec::new();
    let mut refused = 0;
    for case in cases(program) {
        let (out, err, status) = run_binary(&example, program, &case);
        let word = ["-conffile", "-dryrun"].into_iter().find(|word| {
            case.words
                .iter()
                .any(|given| given == word || given.starts_with(&format!("{word}=")))
        });
        let expected = match word {
            Some(word) => {
                refused += 1;
                let text = format!(
                    "Error: the NCBI C++ Toolkit option {word} runs on the command line of \
                     LOSAT's {refusal_name} only: the library API and the Web engine run no NCBI \
                     application start, which reads it\n"
                );
                (Vec::new(), text.into_bytes(), 1)
            }
            None => {
                let (out, err, status) = expected(program, &case.name);
                (out, err, i32::from(status != 0))
            }
        };
        if (out, err, status) != expected {
            differ.push(case.name.clone());
        }
    }
    assert!(refused > 0);
    differ
}

#[test]
fn dustmasker_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("dustmasker", "DUSTMASKER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn segmasker_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("segmasker", "SEGMASKER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn convert2blastmask_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("convert2blastmask", "CONVERT2BLASTMASK");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn the_shell_words_of_the_table_are_read_as_the_shell_reads_them() {
    let words = |line: &str| shell_words(line);
    assert_eq!(
        words("-in in/basic.fa -window \" 5\""),
        (
            vec![
                "-in".to_string(),
                "in/basic.fa".into(),
                "-window".into(),
                " 5".into()
            ],
            None
        )
    );
    assert_eq!(
        words("-outfmt \"\""),
        (vec!["-outfmt".to_string(), String::new()], None)
    );
    assert_eq!(
        words("-masking_options 'a\"b' \"c'd\" ''"),
        (
            vec![
                "-masking_options".to_string(),
                "a\"b".into(),
                "c'd".into(),
                String::new()
            ],
            None
        )
    );
    assert_eq!(
        words("-in - < in/basic.fa"),
        (
            vec!["-in".to_string(), "-".into()],
            Some("in/basic.fa".into())
        )
    );
}
