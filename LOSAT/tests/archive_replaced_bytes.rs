//! The BLAST archive (outfmt 11, `ARCHIVE_FORMAT` text, `bin`, `xml`) and the search strategy
//! (`-export_search_strategy`) of the five searches with strings that NCBI's serial streams do not
//! write as they are: a byte outside `0x20..=0x7E` in ASN.1 text and binary, a control byte other
//! than TAB, LF and CR in XML. NCBI writes `#` and posts `Error: (803.7) [<program>] <stack
//! trace>` / `Bad char [0x<hex>] in string at line N|byte N` / the string for each byte; the posts
//! of a search strategy request are messages of the batch's archive. Each row of
//! `tests/archive_replaced/cases.tsv` gives NCBI BLAST+ 2.17.0's output file, standard output,
//! standard error, strategy file and log file, frozen in `tests/archive_replaced/ncbi/<name>.*`
//! (probe of line RG, item K4b, 2026-11-14, the local official release; an absent file is
//! empty), and its exit status (`<name>.rc`, 0 when absent). NCBI is not run here. The titles (`é`, DEL, 0xFF), a `-parse_deflines`
//! local id with DEL, the `-db` name `in/dé/ndb` and the `-window_masker_db` paths `in/wé/...`
//! and `in/w\x01/...` come from `in/` and from the database and counts files that the test puts
//! under those names in a directory of its own (the repository holds no such names; the
//! control byte's directory only on Unix). Every row without `-logfile` and
//! `-export_search_strategy` runs through the library API too (`examples/library_command_line`),
//! with the same bytes.
//!
//! NCBI posts each byte when it writes it, and its `cerr` is tied to `cout`: in one file (`2>&1`)
//! a post follows the bytes of the archive that NCBI's serial buffer (`COStreamBuffer`, 4096 bytes
//! that keep the last 80 in ASN.1 text, 1 in XML, none in binary) had written by then. The rows
//! `o_*` write the archive to the standard output, with the byte after one or more flushes (three
//! in text and binary; a title of 1500 bytes, whose runs go through `COStreamBuffer::Write` in
//! binary; XML; three query batches; a `-parse_deflines` id in the Seq-aligns): for them the
//! table gives NCBI's merged stream (`<name>.merged`, NCBI run with `2>&1`) and standard error
//! (`<name>.err`, a run of its own), and the test compares both runs of LOSAT (the merged stream
//! holds the standard output; RG K4m, 2026-11-14). The rows `o_bp_*`, `o_bx_*`, `o_tbn_*` and
//! `o_tbx_*` do the same for the archives of blastp, blastx, tblastn and tblastx, and the rows
//! `k_*` for the search strategy's posts, which NCBI writes at each query batch's start after the
//! batch's reader messages, before its report (DECISIONS K4b-2), with the reports on the standard
//! output; `k_*_pf` stop the first batch with NCBI's exception after its posts, whose exit
//! status `<name>.rc` gives (the status is 0 without the file; RG K4b review L7, 2026-11-14).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::library_examples::example_binary;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/archive_replaced")
}

/// One row of the table.
struct Case {
    name: String,
    program: String,
    env: Vec<(String, String)>,
    words: Vec<String>,
}

/// `%XX` of a cell is the byte XX (the table's words are otherwise printable).
fn decode(cell: &str) -> String {
    let mut bytes = Vec::new();
    let mut rest = cell.as_bytes();
    while let Some((&first, tail)) = rest.split_first() {
        if first == b'%' {
            let hex = std::str::from_utf8(&tail[..2]).unwrap();
            bytes.push(u8::from_str_radix(hex, 16).unwrap());
            rest = &tail[2..];
        } else {
            bytes.push(first);
            rest = tail;
        }
    }
    String::from_utf8(bytes).expect("the table's words are UTF-8")
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            let env = match cells[2] {
                "-" => Vec::new(),
                text => text
                    .split(';')
                    .map(|pair| {
                        let (name, value) = pair.split_once('=').unwrap();
                        (name.to_string(), value.to_string())
                    })
                    .collect(),
            };
            Case {
                name: cells[0].to_string(),
                program: cells[1].to_string(),
                env,
                words: cells[3..].iter().map(|cell| decode(cell)).collect(),
            }
        })
        .collect()
}

/// Whether the row's paths exist on this platform: a directory name with a control byte does
/// not on Windows.
fn runs_here(case: &Case) -> bool {
    !cfg!(windows) || !case.words.iter().any(|word| word.contains('\u{1}'))
}

/// A directory where the rows run: `in/` of the fixtures with the database under `in/dé/` and the
/// WindowMasker counts under `in/wé/` and (Unix) `in/w\x01/`, made for `tag` once per run.
fn run_directory(tag: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("archive_replaced-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let input = dir.join("in");
    fs::create_dir_all(&input).unwrap();
    for entry in fs::read_dir(fixtures().join("in")).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), input.join(entry.file_name())).unwrap();
    }
    let database = input.join("d\u{e9}");
    fs::create_dir_all(&database).unwrap();
    for entry in fs::read_dir(fixtures().join("db")).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), database.join(entry.file_name())).unwrap();
    }
    let mut counts = vec!["w\u{e9}"];
    if cfg!(unix) {
        counts.push("w\u{1}");
    }
    for name in counts {
        fs::create_dir_all(input.join(name)).unwrap();
        fs::copy(
            fixtures().join("wm/counts.obinary"),
            input.join(name).join("counts.obinary"),
        )
        .unwrap();
    }
    fs::create_dir_all(dir.join("files")).unwrap();
    dir
}

/// What a row gives: the files of `-out`, `-export_search_strategy` and `-logfile`, the
/// standard output and error, the exit status, and for the rows with a merged stream the
/// standard output and error of a second run in one file.
#[derive(Debug, PartialEq, Eq)]
struct Streams {
    out: Vec<u8>,
    strat: Vec<u8>,
    log: Vec<u8>,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    merged: Option<Vec<u8>>,
    status: i32,
}

/// `case` run by `binary` in `dir`, and with `merged` once more with the standard output and
/// error in one file (`2>&1`).
fn run(binary: &Path, case: &Case, dir: &Path, merged: bool) -> Streams {
    let file = |ext: &str| format!("files/{}.{ext}", case.name);
    let mut streams = run_once(binary, case, dir, None);
    if merged {
        let path = dir.join(file("merged"));
        let one = fs::File::create(&path).unwrap();
        let again = run_once(binary, case, dir, Some(one));
        assert_eq!(
            again.status, streams.status,
            "{}: the exit status of the second run",
            case.name
        );
        streams.merged = Some(fs::read(&path).unwrap());
    }
    streams
}

/// One run of `case`; with `one`, its standard output and error go to that file.
fn run_once(binary: &Path, case: &Case, dir: &Path, one: Option<fs::File>) -> Streams {
    let file = |ext: &str| format!("files/{}.{ext}", case.name);
    for ext in ["out", "strat", "log"] {
        let _ = fs::remove_file(dir.join(file(ext)));
    }
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        // The worker log stays: a build with the `blastx-worker-probe` feature (the gate's) panics
        // without it.
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI")
            || text.starts_with("BLAST")
            || matches!(
                text.as_ref(),
                "ARCHIVE_FORMAT"
                    | "BATCH_SIZE"
                    | "DIAG_POST_LEVEL"
                    | "PRE_FETCH_SEQS_LIMIT"
                    | "CHUNK_SIZE"
                    | "OVERLAP_CHUNK_SIZE"
            )
        {
            command.env_remove(&name);
        }
    }
    command.env("BLAST_USAGE_REPORT", "false");
    for (name, value) in &case.env {
        command.env(name, value);
    }
    let words = case.words.iter().map(|word| match word.as_str() {
        "{out}" => file("out"),
        "{ss}" => file("strat"),
        "{log}" => file("log"),
        _ => word.clone(),
    });
    command.current_dir(dir).arg(&case.program).args(words);
    if let Some(one) = one {
        command.stdout(one.try_clone().unwrap()).stderr(one);
    }
    let output = command.output().unwrap();
    let read = |ext: &str| fs::read(dir.join(file(ext))).unwrap_or_default();
    Streams {
        out: read("out"),
        strat: read("strat"),
        log: read("log"),
        stdout: output.stdout,
        stderr: output.stderr,
        merged: None,
        status: output.status.code().unwrap(),
    }
}

fn expected(case: &Case) -> Streams {
    let read = |ext: &str| {
        fs::read(fixtures().join("ncbi").join(format!("{}.{ext}", case.name))).unwrap_or_default()
    };
    let merged = fixtures()
        .join("ncbi")
        .join(format!("{}.merged", case.name));
    let status = String::from_utf8(read("rc")).unwrap();
    Streams {
        out: read("out"),
        strat: read("strat"),
        log: read("log"),
        stdout: read("stdout"),
        stderr: read("err"),
        merged: fs::read(merged).ok(),
        status: match status.trim() {
            "" => 0,
            text => text.parse().unwrap(),
        },
    }
}

/// The names of the rows whose run through `binary` differs from NCBI's, printed with the first
/// stream that differs; `library` leaves out the rows with `-logfile` and
/// `-export_search_strategy` (the library refuses both).
fn rows_that_differ(binary: &Path, tag: &str, library: bool) -> Vec<String> {
    let dir = run_directory(tag);
    let mut differ = Vec::new();
    let mut ran = 0;
    for case in cases() {
        let refused_by_library = case
            .words
            .iter()
            .any(|word| word == "-logfile" || word == "-export_search_strategy");
        if !runs_here(&case) || (library && refused_by_library) {
            continue;
        }
        ran += 1;
        let mut want = expected(&case);
        let got = run(binary, &case, &dir, want.merged.is_some());
        if want.merged.is_some() {
            // The merged stream holds the standard output, which the table does not repeat.
            want.stdout = got.stdout.clone();
        }
        if got != want {
            let none = Vec::new();
            let first = [
                ("out", &got.out, &want.out),
                ("strat", &got.strat, &want.strat),
                ("log", &got.log, &want.log),
                ("stdout", &got.stdout, &want.stdout),
                ("stderr", &got.stderr, &want.stderr),
                (
                    "merged stream",
                    got.merged.as_ref().unwrap_or(&none),
                    want.merged.as_ref().unwrap_or(&none),
                ),
            ]
            .into_iter()
            .find(|(_, got, want)| got != want);
            match first {
                Some((stream, got, want)) => eprintln!(
                    "{}: {stream} differs\n--- LOSAT\n{}\n--- NCBI\n{}",
                    case.name,
                    String::from_utf8_lossy(got),
                    String::from_utf8_lossy(want)
                ),
                None => eprintln!(
                    "{}: exit status {} (NCBI {})",
                    case.name, got.status, want.status
                ),
            }
            differ.push(case.name);
        }
    }
    assert!(ran > 30, "{ran} rows ran");
    let _ = fs::remove_dir_all(&dir);
    differ
}

#[test]
fn the_archive_and_the_strategy_give_ncbis_bytes_and_posts() {
    let differ = rows_that_differ(Path::new(env!("CARGO_BIN_EXE_losat")), "cli", false);
    assert!(differ.is_empty(), "rows that differ: {differ:?}");
}

#[test]
fn the_library_gives_the_command_line_s_bytes_and_posts() {
    let example = example_binary("library_command_line");
    let differ = rows_that_differ(&example, "library", true);
    assert!(differ.is_empty(), "rows that differ: {differ:?}");
}
