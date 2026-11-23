//! NCBI BLAST+ 2.17.0's argument reader (and RpsbProc 0.5.1's) on the frozen syntax cases of
//! line RG item K4i (`LOSAT/tests/unit/ncbi_syntax_cases.tsv.gz`: 2,072 cases of every argument
//! error of the 10 programs, 1,892 with two dependency violations, 52 with two requirements; made
//! by remaining-gaps-20261112/K/k4i2/tools/make_syntax_fixture.py, line RG item K4i-2b). The
//! reader that the command line, the library API and LOSAT Web share (`reader::read_argv`) holds
//! NCBI's verdict on each: an argument error is NCBI's report byte for byte in LOSAT's words (the
//! short usage, the line of `=`, NCBI's posts), a help word prints that help, and the words that
//! NCBI's reader takes are taken. The cases are read as the library's start reads them
//! (`Caller::OwnStart`) with `cerr` as the target of the posts: a report is the one that NCBI
//! writes without `-logfile` (with `-logfile`, NCBI's posts are the log's and its standard error
//! points to them: `tests/syntax_cases_cli.rs` runs those cases through the command line's entry,
//! and `api::start`'s tests through a library call), and the versions, which the command line
//! prints before the reader (`tests/syntax_cases_cli.rs`), are refused.

use std::ffi::OsString;
use std::io::Read;

use super::reader::{read_argv, ArgvError, Caller};
use super::Commands;
use crate::blastinput::ncbi_environment::DiagTarget;

/// One frozen case: NCBI's verdict (`kind`), its posts (`report`) and the words after the
/// program's name.
struct Case<'a> {
    name: &'a str,
    program: &'a str,
    kind: &'a str,
    posts: &'a str,
    words: Vec<&'a str>,
}

fn fixture() -> String {
    let gz = include_bytes!("../../tests/unit/ncbi_syntax_cases.tsv.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(&gz[..])
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn cases(text: &str) -> Vec<Case<'_>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut fields = line.split('\t');
            let mut next = || fields.next().unwrap();
            let (name, program, _status, kind, posts) = (next(), next(), next(), next(), next());
            Case {
                name,
                program,
                kind,
                posts,
                words: fields.collect(),
            }
        })
        .collect()
}

/// NCBI's post lines of `compact` (the fixture's form: lines joined by `\n`, `@eX` the line
/// before it again after `Error:  (CArgException::eX) `), each with its line end.
fn posts(compact: &str) -> String {
    let mut lines: Vec<String> = Vec::new();
    for line in compact.split("\\n") {
        let line = match line.strip_prefix('@') {
            Some(code) => {
                let before = lines.last().expect("a line before `@`");
                let message = before.strip_prefix("Error: ").expect("a post");
                format!("Error:  (CArgException::{code}) {message}")
            }
            None => line.to_string(),
        };
        lines.push(line);
    }
    lines.iter().map(|line| format!("{line}\n")).collect()
}

fn read(program: &str, words: &[&str]) -> Result<Commands, ArgvError> {
    let argv: Vec<OsString> = ["losat", program]
        .iter()
        .chain(words)
        .map(OsString::from)
        .collect();
    read_argv(&argv, Caller::OwnStart, &DiagTarget::Stderr)
}

fn printed(program: &str, word: &str) -> String {
    match read(program, &[word]) {
        Err(ArgvError::Printed(text)) => text,
        other => panic!("{program} {word}: {other:?}"),
    }
}

/// What `read` gave, for a failure's line.
fn shown(result: &Result<Commands, ArgvError>) -> String {
    match result {
        Ok(_) => "read".to_string(),
        Err(ArgvError::Report(bytes)) => {
            let text = String::from_utf8_lossy(bytes);
            let posts = text.split_once("=\n\n").map_or(&*text, |(_, posts)| posts);
            format!("report {posts:?}")
        }
        Err(ArgvError::Printed(text)) => format!("printed {:?}", &text[..text.len().min(40)]),
        Err(ArgvError::Refused(text)) => format!("refused {text:?}"),
    }
}

#[test]
fn the_reader_gives_ncbis_verdict_on_every_frozen_syntax_case() {
    let text = fixture();
    let cases = cases(&text);
    assert_eq!(cases.len(), 4016);
    let rule = format!("{}\n\n", "=".repeat(72));
    let mut failures = Vec::new();
    for case in &cases {
        let result = read(case.program, &case.words);
        let ncbi = match (case.kind, &result) {
            ("report", Err(ArgvError::Report(bytes))) => {
                let posts = posts(case.posts);
                let expected = format!(
                    "{}{rule}{}",
                    printed(case.program, "-h"),
                    crate::wording::losat_text(&posts)
                );
                *bytes == expected.as_bytes()
            }
            ("h" | "help" | "help-full", Err(ArgvError::Printed(text))) => {
                *text == printed(case.program, &format!("-{}", case.kind))
            }
            ("version", Err(ArgvError::Refused(text))) => {
                *text
                    == super::version_word_refusal(
                        super::reader::tool(case.program).unwrap(),
                        "-version",
                    )
            }
            ("accepted", Ok(_)) => true,
            _ => false,
        };
        if !ncbi {
            failures.push(format!(
                "{} {} {:?}: NCBI {} {:?}, LOSAT {}",
                case.name,
                case.program,
                case.words,
                case.kind,
                case.posts,
                shown(&result)
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} cases differ:\n{}",
        failures.len(),
        cases.len(),
        failures[..failures.len().min(40)].join("\n")
    );
}
