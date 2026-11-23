//! The command line's entry (`cli::main_entry`) on the frozen syntax cases of line RG item K4i
//! (`tests/unit/ncbi_syntax_cases.tsv.gz`, which `src/cli/syntax_cases_tests.rs` reads through the
//! shared reader) that the reader's callers cannot show: the 98 argument errors with `-logfile`,
//! whose posts NCBI writes to the file of `-logfile` while its standard error says where they are,
//! and the 138 cases with a version word, which the entry prints before the arguments are read.
//! Each runs the `losat` binary in a folder of its own and compares its bytes with NCBI BLAST+
//! 2.17.0's (and RpsbProc 0.5.1's) frozen ones in LOSAT's words (`losat::wording::losat_output_bytes`):
//! `tests/fixtures/syntax_logfile/<program>.err` and `.log` are the standard error and the file of
//! `-logfile log` of `<program> <its inputs> -logfile log -bogus 1` (NCBI's case
//! `<program>.log.unk`, frozen on 2026-11-18 by `remaining-gaps-20261112/K/k4i2/syntax_probe.py`;
//! every `-logfile` case of the fixture has the same standard error, and its posts are the
//! fixture's), `.version` the standard output of `<program> -version` (every version case of the
//! fixture prints it, exit status 0, nothing on the standard error).
//!
//! The version text is LOSAT's where NCBI's names its own build: the ` Package:` line names LOSAT
//! (`ncbi_version::print`, line RG DECISIONS E3-1), and rpsbproc's first line names LOSAT's
//! version where RpsbProc writes `0.5.1` (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`); the other lines
//! are NCBI's in LOSAT's words.

#![cfg(unix)]

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// One frozen case: NCBI's verdict (`kind`), its posts (`posts`) and the words after the
/// program's name.
struct Case {
    name: String,
    program: String,
    kind: String,
    posts: String,
    words: Vec<String>,
}

fn cases() -> Vec<Case> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/unit/ncbi_syntax_cases.tsv.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut fields = line.split('\t').map(str::to_string);
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

/// NCBI's frozen file `<program>.<extension>` of `tests/fixtures/syntax_logfile`, in LOSAT's words.
fn frozen(program: &str, extension: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/syntax_logfile")
        .join(format!("{program}.{extension}"));
    let bytes = std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    losat::wording::losat_output_bytes(&bytes).into_owned()
}

/// A folder of the case `name` alone, empty.
fn folder(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-syntax-cli-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// `losat <program> <words>` in `dir`, in an environment of NCBI's probes (`env -i PATH HOME`).
fn run_in(dir: &Path, case: &Case) -> Output {
    let mut command = Command::new(BIN);
    command
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .current_dir(dir)
        .arg(&case.program)
        .args(&case.words)
        .stdin(Stdio::null());
    // The `--all-features` build's BLASTX diagnostic probe needs its log (as in CI); no case
    // searches.
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        command.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    command.output().unwrap()
}

fn shown(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

#[test]
fn an_argument_error_with_a_logfile_is_ncbis_on_the_standard_error_and_in_the_file() {
    let cases = cases();
    let logged: Vec<&Case> = cases
        .iter()
        .filter(|case| case.kind == "report")
        .filter(|case| case.words.iter().any(|word| word.starts_with("-logfile")))
        .collect();
    assert_eq!(logged.len(), 98);
    let mut failures = Vec::new();
    for case in logged {
        let dir = folder(&case.name);
        let output = run_in(&dir, case);
        let log = std::fs::read(dir.join("log")).unwrap_or_default();
        let want_log = losat::wording::losat_text(&posts(&case.posts)).into_owned();
        // The fixture's posts are NCBI's file of `-logfile` (its case `<program>.log.unk`).
        if case.name.ends_with(".log.unk") {
            assert_eq!(
                shown(&frozen(&case.program, "log")),
                want_log,
                "{}",
                case.name
            );
        }
        let want_err = frozen(&case.program, "err");
        if (output.status.code(), &output.stdout, &output.stderr, &log)
            != (
                Some(1),
                &Vec::new(),
                &want_err,
                &want_log.clone().into_bytes(),
            )
        {
            failures.push(format!(
                "{} {:?}: rc {:?}, stdout {:?}, stderr ends {:?}, log {:?} (NCBI's log {want_log:?})",
                case.name,
                case.words,
                output.status.code(),
                shown(&output.stdout),
                shown(&output.stderr[output.stderr.len().saturating_sub(200)..]),
                shown(&log),
            ));
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn a_version_word_prints_ncbis_version_text_in_losats_words() {
    let cases = cases();
    let versions: Vec<&Case> = cases.iter().filter(|case| case.kind == "version").collect();
    assert_eq!(versions.len(), 138);
    let mut failures = Vec::new();
    for case in versions {
        let dir = folder(&case.name);
        let output = run_in(&dir, case);
        let text = shown(&output.stdout);
        let ncbi = shown(&frozen(&case.program, "version"));
        // The lines compared: all but those that name a build, the package's and rpsbproc's first
        // (its own version).
        let compared_lines = |text: &str| -> Vec<String> {
            text.lines()
                .enumerate()
                .filter(|(at, line)| {
                    let names_a_build =
                        line.starts_with(" Package: ") || (case.program == "rpsbproc" && *at == 0);
                    !names_a_build
                })
                .map(|(_, line)| line.to_string())
                .collect()
        };
        let package = text.lines().nth(1).unwrap_or_default();
        let first = text.lines().next().unwrap_or_default();
        if output.status.code() != Some(0)
            || !output.stderr.is_empty()
            || text.lines().count() != ncbi.lines().count()
            || !text.ends_with('\n')
            || compared_lines(&text) != compared_lines(&ncbi)
            || !package.starts_with(" Package: LOSAT ")
            || first != format!("{}: {}", case.program, losat::wording::LOSAT_VERSION)
        {
            failures.push(format!(
                "{} {:?}: rc {:?}, stdout {text:?}, stderr {:?} (NCBI {ncbi:?})",
                case.name,
                case.words,
                output.status.code(),
                shown(&output.stderr),
            ));
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
