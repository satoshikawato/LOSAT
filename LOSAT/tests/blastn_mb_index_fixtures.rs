//! `losat blastn -use_index true` against NCBI BLAST+ 2.17.0 (line RG items D5b-1 to D5b-5):
//! each row of `tests/blastn_mb_index/cases.tsv` runs in its own copy of
//! `tests/blastn_mb_index/in`, after its setup steps (LOSAT's `makembindex`, whose index files are
//! NCBI's bytes, LOSAT's `makeblastdb` of FASTA files made from the genomes of `tests/fasta`
//! (`@TESTS@` in a step is `tests`), and edits of the files), and gives NCBI's standard output,
//! standard error and exit status (`ncbi/<name>.{out,err,rc}`, the row's directory written
//! `@CWD@`), frozen with the local official release and NCBI's makembindex and makeblastdb. NCBI is
//! not run here. A row's first `NAME=value` words are variables of its run
//! (`masker_fixtures::take_variables`).
//!
//! A row with a refusal is one where NCBI reads the index in a way that leaves the search
//! undefined or fails:
//! LOSAT refuses it with that text, and NCBI's output is kept for the lines that take these
//! searches. Where NCBI exits with its error, crashes or never ends in a later query batch, LOSAT
//! writes what NCBI wrote before (its output of the batches before) and refuses there; where
//! NCBI's run succeeds, LOSAT writes nothing, or the batches before for a refusal of a later query
//! batch or query chunk (its text names one).
//!
//! A row whose fourth field is `same as NAME` (the setup and the words of row NAME, with another
//! `-num_threads`) is compared with row NAME's NCBI files: NCBI's `-num_threads 1` output, which
//! LOSAT gives at every thread count (Owner 2026-10-13), for a run whose NCBI result at more
//! threads depends on its threads. A row whose fourth field is `plain NAME` (row NAME's search
//! with `-use_index true` (and `-index_name` for an old-style index) and `-db_hard_mask`, its
//! setup with the index made, at any `-num_threads`) is compared with row NAME's NCBI files: NCBI's search without the index, whose
//! output LOSAT's one-hit search with the index gives there, the masked letters matching nothing
//! and the seeds in their hard ranges (the Owner's answer C to Q-HARDIDX, approved exception 24 of
//! PD-LOSAT-NCBI-DEFECTS). A row whose fourth field is `expect NAME` (a search with the index and
//! `-db_hard_mask`) is compared with `exception24/NAME.{out,err,rc}`: LOSAT's output where no NCBI
//! search gives it, the two-hit search (`-window_size`), whose index and scan pair other hits; that
//! file is NCBI's search without the index with the rows that
//! `the_two_hit_index_search_adds_the_index_s_pairs` names, and its rows are the same at every
//! `-num_threads` and on the database whose masked letters are others (hmr).

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
    /// The row whose NCBI files this row is compared with (itself, or that of `same as NAME` or
    /// `plain NAME`), or the name of its `exception24` files (`expect NAME`).
    ncbi: String,
    /// A `plain NAME` row.
    plain: bool,
    /// The directory of the files that the row is compared with: `ncbi`, or `exception24` for
    /// an `expect NAME` row.
    files: &'static str,
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
        let same_as = fields.get(3).and_then(|text| text.strip_prefix("same as "));
        let plain = fields.get(3).and_then(|text| text.strip_prefix("plain "));
        let expect = fields.get(3).and_then(|text| text.strip_prefix("expect "));
        cases.push(Case {
            name: fields[0].to_string(),
            setup,
            words,
            variables,
            refusal: fields
                .get(3)
                .filter(|_| same_as.is_none() && plain.is_none() && expect.is_none())
                .map(|text| text.to_string()),
            ncbi: same_as
                .or(plain)
                .or(expect)
                .unwrap_or(fields[0])
                .to_string(),
            plain: plain.is_some(),
            files: if expect.is_some() {
                "exception24"
            } else {
                "ncbi"
            },
        });
    }
    for case in cases.iter().filter(|case| case.files == "exception24") {
        assert!(
            case.words
                .windows(2)
                .any(|pair| pair == ["-use_index", "true"])
                && case.words.iter().any(|word| word == "-db_hard_mask"),
            "{}: an expect row is a search with the index and -db_hard_mask",
            case.name
        );
    }
    // A `same as` row runs its row's search with another thread count; a `plain` row, with the
    // index (made in its setup, named by `-index_name` for an old-style one) and the hard mask,
    // at any thread count.
    let without_threads = |words: &[String], index: bool| {
        let mut kept = Vec::new();
        let mut words = words.iter();
        while let Some(word) = words.next() {
            if word == "-num_threads" || (index && (word == "-use_index" || word == "-index_name"))
            {
                words.next();
            } else {
                kept.push(word.clone());
            }
        }
        kept
    };
    let without_index = |setup: &[Vec<String>]| -> Vec<Vec<String>> {
        setup
            .iter()
            .filter(|step| step[0] != "makembindex")
            .cloned()
            .collect()
    };
    for case in cases
        .iter()
        .filter(|case| case.ncbi != case.name && case.files == "ncbi")
    {
        let row = cases
            .iter()
            .find(|row| row.name == case.ncbi)
            .unwrap_or_else(|| panic!("{}: no row {}", case.name, case.ncbi));
        let same = if case.plain {
            row.ncbi == row.name
                && without_index(&row.setup) == without_index(&case.setup)
                && case
                    .words
                    .windows(2)
                    .any(|pair| pair == ["-use_index", "true"])
                && case.words.iter().any(|word| word == "-db_hard_mask")
                && !row.words.iter().any(|word| word == "-use_index")
        } else {
            row.ncbi == row.name && row.setup == case.setup
        };
        assert!(
            same && row.variables == case.variables
                && without_threads(&row.words, case.plain)
                    == without_threads(&case.words, case.plain),
            "{}: not row {}'s search with another -num_threads (or with the index)",
            case.name,
            row.name
        );
    }
    cases
}

/// A command of LOSAT with no LOSAT, NCBI or BLAST variable, run in `dir` (the worker log of
/// the `blastx-worker-probe` feature is kept: its threaded searches write to it).
fn command(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI_")
            || text.starts_with("BLAST")
        {
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
                // Where NCBI exits with its error (its search of a later query batch), LOSAT
                // writes NCBI's output of the batches before and refuses there.
                let ncbi = |extension: &str| {
                    fs::read(
                        fixtures()
                            .join(case.files)
                            .join(format!("{}.{extension}", case.ncbi)),
                    )
                    .unwrap()
                };
                let ncbi_rc: i32 = String::from_utf8(ncbi("rc"))
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap();
                let want_out = {
                    let out = here(&String::from_utf8(ncbi("out")).unwrap(), &dir);
                    losat::wording::losat_output_bytes(out.as_bytes()).into_owned()
                };
                // NCBI's run succeeded: LOSAT refuses before the search (nothing written), or a
                // later query batch or query chunk after the batches before. Elsewhere (NCBI's
                // error, its crash or its time limit, 124) what NCBI wrote before.
                let later = refusal.contains("query batch or query chunk");
                let out_ok = match ncbi_rc {
                    0 if later => want_out.starts_with(&output.stdout),
                    0 => output.stdout.is_empty(),
                    _ => output.stdout == want_out,
                };
                output.status.code() == Some(1) && out_ok && stderr.contains(&line)
            }
            None => {
                let file = |extension: &str| {
                    fs::read(
                        fixtures()
                            .join(case.files)
                            .join(format!("{}.{extension}", case.ncbi)),
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

/// The rows of `exception24/hm_i_win.out` (LOSAT's two-hit search with the index and
/// `-db_hard_mask 11`, approved exception 24) are NCBI's search without the index
/// (`ncbi/hm_n_win.out`) and three rows that the index's two-hit search finds and the scan's does
/// not: NCBI's own search with the index and without the mask finds the row of qhm7 in hm4 too
/// (`ncbi/hm_iu_win.out`; hm4's only mask is 700 letters away); the rows of qhm2 in hm0 and of
/// qhm8 in hm5 are stretches that a mask ends (hm0 999..1060, hm5 2304..2562 with 2238..2304), of
/// exact runs a mismatch apart that the index's search pairs within the window (it pairs a hit
/// with one that ends at most the window before its start; NCBI's scan, before its end).
///
/// `exception24/hm_win_stretches.out` is NCBI's own two-hit search with the index of the database
/// of hm's unmasked stretches, each its own sequence (`hm0_1061`: hm0's letters from 1,061), with
/// each query's search space of hm's search with `-db_hard_mask 11` (`-searchsp`), its rows mapped
/// back to hm's sequences and letters (`remaining-gaps-20261112/D/d5/probe/d5b5/st_oracle.py`). Its
/// rows are LOSAT's, the three above included, but for the three rows that cross hm0's 7-letter
/// mask at 1,502..1,508: NCBI's search without the index (and LOSAT's) runs the traceback on the
/// whole sequence with the masked letters as B, across the mask, where the stretches end.
#[test]
fn the_two_hit_index_search_adds_the_index_s_pairs() {
    let read = |path: &str| fs::read_to_string(fixtures().join(path)).unwrap();
    let added = [
        "qhm2\thm0\t98.990\t99\t1\t0\t1152\t1250\t999\t901\t4.13e-47\t178",
        "qhm7\thm4\t97.100\t1000\t29\t0\t1\t1000\t301\t1300\t0.0\t1687",
        "qhm8\thm5\t97.959\t196\t4\t0\t2105\t2300\t2305\t2500\t9.21e-96\t340",
    ];
    let expected = read("exception24/hm_i_win.out");
    let without: Vec<&str> = expected
        .lines()
        .filter(|line| !added.contains(line))
        .collect();
    assert_eq!(
        without,
        read("ncbi/hm_n_win.out").lines().collect::<Vec<_>>()
    );
    assert_eq!(expected.lines().count(), without.len() + added.len());
    assert!(read("ncbi/hm_iu_win.out")
        .lines()
        .any(|line| line == added[1]));
    for extension in ["err", "rc"] {
        assert_eq!(
            read(&format!("exception24/hm_i_win.{extension}")),
            read(&format!("ncbi/hm_n_win.{extension}"))
        );
    }
    let stretches = read("exception24/hm_win_stretches.out");
    let crossing: Vec<&str> = expected
        .lines()
        .filter(|line| !stretches.lines().any(|row| row == *line))
        .collect();
    assert_eq!(
        crossing,
        [
            "qhm0\thm0\t99.222\t900\t7\t0\t261\t1160\t1061\t1960\t0.0\t1637",
            "qhm1\thm0\t97.444\t900\t23\t0\t211\t1110\t1061\t1960\t0.0\t1548",
            "qhm2\thm0\t97.778\t900\t20\t0\t191\t1090\t1960\t1061\t0.0\t1565",
        ]
    );
}

#[test]
fn blastn_use_index_gives_ncbis_bytes_or_refuses_an_indexed_search() {
    let differ = rows_that_differ();
    assert!(differ.is_empty(), "rows that differ: {differ:?}");
}
