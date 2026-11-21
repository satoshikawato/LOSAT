//! `losat blastn -use_index true` against NCBI BLAST+ 2.17.0 (line RG item D5b-1): each row of
//! `tests/blastn_mb_index/cases.tsv` runs in its own copy of `tests/blastn_mb_index/in`, after
//! its setup steps (LOSAT's `makembindex`, whose index files are NCBI's bytes, and edits of the
//! files), and gives NCBI's standard output, standard error and exit status
//! (`ncbi/<name>.{out,err,rc}`, the row's directory written `@CWD@`), frozen with the local
//! official release and NCBI's makembindex. NCBI is not run here.
//!
//! A row with a refusal is one where NCBI searches with the index (or reads it in a way that
//! leaves the search undefined): LOSAT refuses it with that text, and NCBI's output is kept for
//! the indexed search of line D5b-2.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::masker_fixtures::shell_words;
use crate::ncbi_windows_paths::ncbi_paths_on_windows;

/// `tests/blastn_mb_index`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastn_mb_index")
}

/// One row of the table.
struct Case {
    name: String,
    /// The setup steps, each as shell words.
    setup: Vec<Vec<String>>,
    words: Vec<String>,
    /// LOSAT's refusal (`@CWD@` for the row's directory), for a row that NCBI searches with
    /// the index.
    refusal: Option<String>,
}

/// The rows of the table.
fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut cases = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let setup = match fields[1] {
            "-" => Vec::new(),
            steps => steps.split(" ; ").map(|step| shell_words(step).0).collect(),
        };
        cases.push(Case {
            name: fields[0].to_string(),
            setup,
            words: shell_words(fields[2]).0,
            refusal: fields.get(3).map(|text| text.to_string()),
        });
    }
    cases
}

/// A command of LOSAT with no LOSAT, NCBI or BLAST variable, run in `dir`.
fn command(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .stdin(Stdio::null());
    command
}

/// Copies the files of `from` into `to`.
fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), to.join(entry.file_name())).unwrap();
    }
}

/// Runs the setup step `step` in `dir`.
fn setup(step: &[String], dir: &Path) {
    let file = |at: usize| dir.join(&step[at]);
    let number = |at: usize| step[at].parse::<u64>().unwrap();
    match step[0].as_str() {
        "makembindex" => {
            let output = command(dir).args(step).output().unwrap();
            assert!(output.status.success(), "{step:?}: {output:?}");
        }
        "put" => {
            let mut bytes = fs::read(file(1)).unwrap();
            let at = number(2) as usize;
            bytes[at..at + 4].copy_from_slice(&(number(3) as u32).to_le_bytes());
            fs::write(file(1), bytes).unwrap();
        }
        "cut" => {
            let bytes = fs::read(file(1)).unwrap();
            fs::write(file(1), &bytes[..number(2) as usize]).unwrap();
        }
        "cp" => {
            fs::copy(file(1), file(2)).unwrap();
        }
        "rm" => fs::remove_file(file(1)).unwrap(),
        other => panic!("unknown step {other}"),
    }
}

/// `text` with `@CWD@` read as `dir` (and NCBI's paths as Windows writes them there).
fn here(text: &str, dir: &Path) -> String {
    let text = if cfg!(windows) {
        ncbi_paths_on_windows(text, &["@CWD@"])
    } else {
        text.to_string()
    };
    text.replace("@CWD@", &dir.to_string_lossy())
}

/// `text` without the ends of line that wrap the message strings of an archive (`message "..."`):
/// the writer wraps a string by its length, which the row's directory in a message changes
/// (the wraps themselves are `wording::losat_output_bytes`'s tests').
fn messages_unwrapped(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find("message") {
        let after = at + "message".len();
        let Some(quote) = rest[after..]
            .find(|c: char| c != ' ' && c != '\n')
            .map(|n| after + n)
            .filter(|&quote| rest[quote..].starts_with('"'))
        else {
            out.push_str(&rest[..after]);
            rest = &rest[after..];
            continue;
        };
        out.push_str(&rest[..=quote]);
        rest = &rest[quote + 1..];
        // To the closing quote (a doubled quote is a quote of the string).
        let bytes = rest.as_bytes();
        let mut end = 0;
        while end < bytes.len() && !(bytes[end] == b'"' && bytes.get(end + 1) != Some(&b'"')) {
            end += if bytes[end] == b'"' { 2 } else { 1 };
        }
        out.push_str(&rest[..end].replace('\n', ""));
        rest = &rest[end..];
    }
    out.push_str(rest);
    out
}

/// The rows whose run differs from NCBI's (or from LOSAT's refusal).
fn rows_that_differ() -> Vec<String> {
    let mut differ = Vec::new();
    for case in cases() {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
            "blastn-mb-index-{}-{}",
            case.name,
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        copy_dir(&fixtures().join("in"), &dir);
        for step in &case.setup {
            setup(step, &dir);
        }
        let output = command(&dir)
            .arg("blastn")
            .args(&case.words)
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        let matches = match &case.refusal {
            Some(refusal) => {
                let line = format!(
                    "-use_index true with {} is not supported by LOSAT's BLASTN",
                    here(refusal, &dir)
                );
                output.status.code() == Some(1)
                    && output.stdout.is_empty()
                    && stderr.contains(&line)
            }
            None => {
                let file = |extension: &str| {
                    fs::read(
                        fixtures()
                            .join("ncbi")
                            .join(format!("{}.{extension}", case.name)),
                    )
                    .unwrap()
                };
                let text = |bytes: Vec<u8>| here(&String::from_utf8(bytes).unwrap(), &dir);
                let err = losat::wording::losat_bytes(&file("err")).into_owned();
                let rc: i32 = String::from_utf8(file("rc"))
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap();
                // An archive's messages are in LOSAT's words.
                let want_out = String::from_utf8(
                    losat::wording::losat_output_bytes(text(file("out")).as_bytes()).into_owned(),
                )
                .unwrap();
                let want_err = text(err);
                let same = output.status.code() == Some(rc)
                    && messages_unwrapped(&String::from_utf8_lossy(&output.stdout))
                        == messages_unwrapped(&want_out)
                    && stderr == want_err;
                if !same {
                    eprintln!(
                        "--- {}: LOSAT status {:?}, NCBI {rc}\nLOSAT stderr: {stderr:?}\n\
                         NCBI stderr:  {want_err:?}\nLOSAT stdout: {:?}\nNCBI stdout:  {:?}",
                        case.name,
                        output.status.code(),
                        String::from_utf8_lossy(&output.stdout)
                            .chars()
                            .take(600)
                            .collect::<String>(),
                        want_out.chars().take(600).collect::<String>(),
                    );
                }
                same
            }
        };
        if matches {
            let _ = fs::remove_dir_all(&dir);
        } else {
            if case.refusal.is_some() {
                eprintln!(
                    "--- {}: LOSAT status {:?}, stderr {stderr:?}",
                    case.name,
                    output.status.code()
                );
            }
            differ.push(case.name.clone());
        }
    }
    differ
}

#[test]
fn blastn_use_index_gives_ncbis_bytes_or_refuses_an_indexed_search() {
    let differ = rows_that_differ();
    assert!(differ.is_empty(), "rows that differ: {differ:?}");
}
