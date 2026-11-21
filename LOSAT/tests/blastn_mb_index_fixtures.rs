//! `losat blastn -use_index true` against NCBI BLAST+ 2.17.0 (line RG items D5b-1 and D5b-2):
//! each row of `tests/blastn_mb_index/cases.tsv` runs in its own copy of
//! `tests/blastn_mb_index/in`, after its setup steps (LOSAT's `makembindex`, whose index files are
//! NCBI's bytes, LOSAT's `makeblastdb` of FASTA files made from the genomes of `tests/fasta`
//! (`@TESTS@` in a step is `tests`), and edits of the files), and gives NCBI's standard output,
//! standard error and exit status (`ncbi/<name>.{out,err,rc}`, the row's directory written
//! `@CWD@`), frozen with the local official release and NCBI's makembindex and makeblastdb. NCBI is
//! not run here. A row's first `NAME=value` words are variables of its run
//! (`masker_fixtures::take_variables`).
//!
//! A row with a refusal is one where NCBI searches with an index in a way that LOSAT does not
//! (a mixed search, a two-hit search, `-db_hard_mask`), or reads the index in a way that leaves
//! the search undefined or fails: LOSAT refuses it with that text, after the output of the query
//! batches before (NCBI's output up to there), and NCBI's output is kept for the lines that take
//! these searches.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::masker_fixtures::{shell_words, take_variables};
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
    /// The variables of the row's first `NAME=value` words.
    variables: Vec<(String, String)>,
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
        let mut words = shell_words(fields[2]).0;
        let variables = take_variables(&mut words);
        cases.push(Case {
            name: fields[0].to_string(),
            setup,
            words,
            variables,
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

/// The sequence lines of the FASTA file `bytes` (each with its end of line).
fn sequence_lines(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len());
    for line in bytes.split(|&byte| byte == b'\n') {
        if !line.is_empty() && line[0] != b'>' {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

/// Runs the setup step `step` in `dir`: `makembindex` and `makeblastdb` (LOSAT's), `put FILE
/// OFFSET VALUE` (a little-endian 32-bit word), `cut FILE SIZE`, `cp FROM TO`, `rm FILE`, `cat OUT
/// FROM...` (the bytes of the FROM files one after another) and `join OUT TITLE FROM...` (one FASTA
/// record `>TITLE` with the sequence lines of the FROM files).
fn setup(step: &[String], dir: &Path) {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let step: Vec<String> = step
        .iter()
        .map(|word| word.replace("@TESTS@", &tests.to_string_lossy()))
        .collect();
    let file = |at: usize| dir.join(&step[at]);
    let number = |at: usize| step[at].parse::<u64>().unwrap();
    match step[0].as_str() {
        "makembindex" | "makeblastdb" => {
            let output = command(dir).args(&step).output().unwrap();
            assert!(output.status.success(), "{step:?}: {output:?}");
        }
        "cat" => {
            let mut bytes = Vec::new();
            for from in &step[2..] {
                bytes.extend(fs::read(dir.join(from)).unwrap());
            }
            fs::write(file(1), bytes).unwrap();
        }
        "join" => {
            let mut bytes = format!(">{}\n", step[2]).into_bytes();
            for from in &step[3..] {
                bytes.extend(sequence_lines(&fs::read(dir.join(from)).unwrap()));
            }
            fs::write(file(1), bytes).unwrap();
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
            .envs(case.variables.iter().map(|(name, value)| (name, value)))
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        let matches = match &case.refusal {
            Some(refusal) => {
                let line = format!(
                    "-use_index true with {} is not supported by LOSAT's BLASTN",
                    here(refusal, &dir)
                );
                let ncbi_out = fs::read(fixtures().join("ncbi").join(format!("{}.out", case.name)))
                    .unwrap();
                let ncbi_out = here(&String::from_utf8(ncbi_out).unwrap(), &dir);
                output.status.code() == Some(1)
                    && losat::wording::losat_output_bytes(ncbi_out.as_bytes())
                        .starts_with(&output.stdout)
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
