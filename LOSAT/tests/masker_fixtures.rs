//! `losat dustmasker`, `losat segmasker` and `losat convert2blastmask` against NCBI BLAST+
//! 2.17.0 dustmasker 1.0.0, segmasker 1.0.0 and convert2blastmask: each row of
//! `tests/<program>/cases.tsv` (its arguments as shell words, run in `tests/<program>`; `< file`
//! is the standard input) gives NCBI's standard output, standard error and exit status, frozen
//! in `tests/<program>/ncbi/<name>.{out,err,rc}` (probes of line RG, items K1, K2 and K3,
//! 2026-11-12 and 11-13, on ws and the local official release). A name the table repeats is
//! compared on its last row, whose output the probe kept. NCBI is not run here. Every row runs
//! through the library API too (`LOSAT::api::run_dustmasker`, `run_segmasker`,
//! `run_convert2blastmask`, by `examples/library_command_line`), with the same bytes.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::blastdbcmd_info::ncbi_paths_on_windows;
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

/// Stdout, stderr and exit status of a run.
type Streams = (Vec<u8>, Vec<u8>, i32);

/// One run of `losat <program>` for `case` in `tests/<program>`: stdout, stderr, exit status.
/// No LOSAT, NCBI or BLAST variable, no registry file (`NCBI_CONFIG_PATH` names the fixture
/// directory, which has none).
fn run(program: &str, case: &Case) -> Streams {
    run_binary(Path::new(env!("CARGO_BIN_EXE_losat")), program, case)
}

/// `run` of `binary <program>` (the command line, or the library's example).
fn run_binary(binary: &Path, program: &str, case: &Case) -> Streams {
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

/// NCBI's stdout, stderr and exit status for `case` of `program`. NCBI's bytes come from Linux;
/// on Windows they are the ones NCBI's code writes there:
/// - a directory does not open as a file (`CArg_InputFile::x_Open`: the `ifstream::open` of
///   MSVC's library fails on a directory, `is_open()` is false, the exception
///   `CArgException::eNoFile`), so `-in` naming one is "File is not accessible" as a file that
///   does not exist, the row `in_missing` with the name of `case` (segmasker and
///   convert2blastmask open `-in` first; dustmasker opens `-out` first, and the interval writer
///   of its directory row writes nothing before a record); on Linux the directory opens and its
///   first read fails;
/// - the paths of the database messages and the search path (`ncbi_paths_on_windows`).
fn expected(program: &str, case: &Case) -> Streams {
    if cfg!(windows) {
        if let Some(directory) = input_directory(program, case) {
            let missing = cases(program)
                .into_iter()
                .find(|row| row.name == "in_missing")
                .unwrap();
            let missing_name = input(&missing).unwrap().to_string();
            let named = |bytes: Vec<u8>| -> Vec<u8> {
                String::from_utf8(bytes)
                    .unwrap()
                    .replace(&format!("`{missing_name}'"), &format!("`{directory}'"))
                    .into_bytes()
            };
            let (out, err, status) = expected_row(program, "in_missing");
            return (named(out), named(err), status);
        }
    }
    expected_row(program, &case.name)
}

/// The value of `-in` in `case`.
fn input(case: &Case) -> Option<&str> {
    let at = case.words.iter().position(|word| word == "-in")?;
    case.words.get(at + 1).map(String::as_str)
}

/// The value of `-in` in `case` when it names a directory.
fn input_directory<'a>(program: &str, case: &'a Case) -> Option<&'a str> {
    input(case).filter(|name| fixtures(program).join(name).is_dir())
}

/// NCBI's stdout, stderr and exit status for the row `name` of `program`, as frozen.
fn expected_row(program: &str, name: &str) -> Streams {
    let file = |extension: &str| {
        fixtures(program)
            .join("ncbi")
            .join(format!("{name}.{extension}"))
    };
    let status = fs::read_to_string(file("rc")).unwrap();
    // The database rows name the probe's directory in NCBI's search path; here it is this one.
    let here = |bytes: Vec<u8>| -> Vec<u8> {
        let text = String::from_utf8_lossy(&bytes);
        let text = if cfg!(windows) {
            ncbi_paths_on_windows(&text, &["@CWD@"])
        } else {
            text.into_owned()
        };
        let dir = fixtures(program);
        text.replace("@CWD@", &dir.to_string_lossy()).into_bytes()
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

/// Whether `case` runs on this platform: a row that writes its output to a device file runs
/// where the file is, `/dev/full` on Linux, `/dev/null` on Unix (`-logfile /dev/null` is NCBI's
/// name for no log file, `CDiagContext::SetLogFile`, on every platform).
fn runs_here(case: &Case) -> bool {
    let out = case
        .words
        .iter()
        .position(|word| word == "-out")
        .and_then(|at| case.words.get(at + 1))
        .map(String::as_str);
    match out {
        Some("/dev/full") => cfg!(target_os = "linux"),
        Some("/dev/null") => cfg!(unix),
        _ => true,
    }
}

/// The rows of `program` that run here (`runs_here`).
fn cases_here(program: &str) -> Vec<Case> {
    cases(program).into_iter().filter(runs_here).collect()
}

/// `rows`, the rows that differ, each with what LOSAT and NCBI wrote (the start of each stream),
/// for the message of a failed test.
fn describe(rows: &[(String, Streams, Streams)]) -> String {
    let start = |bytes: &[u8]| {
        let text = String::from_utf8_lossy(bytes);
        text.chars().take(600).collect::<String>()
    };
    let mut text = String::new();
    for (name, got, want) in rows {
        text.push_str(&format!(
            "\n--- {name}: LOSAT status {}, NCBI {}\nLOSAT stdout: {:?}\nNCBI stdout:  {:?}\n\
             LOSAT stderr: {:?}\nNCBI stderr:  {:?}\n",
            got.2,
            want.2,
            start(&got.0),
            start(&want.0),
            start(&got.1),
            start(&want.1),
        ));
    }
    text
}

/// The rows of `program`'s table whose run differs from NCBI's.
fn rows_that_differ(program: &str) -> Vec<String> {
    let mut differ = Vec::new();
    for case in cases_here(program) {
        let (got, want) = (run(program, &case), expected(program, &case));
        if got != want {
            differ.push((case.name.clone(), got, want));
        }
    }
    if !differ.is_empty() {
        eprintln!("{}", describe(&differ));
    }
    differ.into_iter().map(|(name, _, _)| name).collect()
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
    for case in cases_here(program) {
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
                let (out, err, status) = expected(program, &case);
                (out, err, i32::from(status != 0))
            }
        };
        let got = (out, err, status);
        if got != expected {
            differ.push((case.name.clone(), got, expected));
        }
    }
    assert!(refused > 0);
    if !differ.is_empty() {
        eprintln!("{}", describe(&differ));
    }
    differ.into_iter().map(|(name, _, _)| name).collect()
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
