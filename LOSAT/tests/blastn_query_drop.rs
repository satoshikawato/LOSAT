//! The query records read before the search held in a store, their letters packed as read and
//! dropped once a batch has encoded them (`LOSAT_X_BNQPACK`,
//! `src/algorithm/blastn/blast_engine/x_query_drop.rs`), on whole searches. Until 2026-11-01 the
//! store had a switch of its own, `LOSAT_X_BNQDROP`, and packed only files of two of the task's
//! largest batches or more, not standard input; both were folded into `LOSAT_X_BNQPACK`.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, which keeps every record's letters to the end of the search. It changes when the queries'
//! letters are freed, not what the search computes: after the
//! encoding the search reads the records' lengths only, and its rule keeps the letters where a
//! report shows them. The property under test is that the reports and the messages do not depend
//! on it and that the default run took the new path where the rule says so: with `LOSAT_TIMING`
//! the default run prints `[TIMING] x_bn_q_drop: batches= records= letters=` and
//! `[TIMING] x_bn_q_pack: records= letters= ...`, and `=0` prints no such line. The formats that
//! show no query letter (outfmt 6, 7, 8, 9, 12 and SAM without `SR`)
//! drop the letters of every record, in a search of several batches (`BATCH_SIZE`) with a query
//! without letters and a query without a valid context among them; the formats that show them
//! (outfmt 0, 5, SAM `SQ SR`) and `-query_loc` keep them; a batch split into query chunks
//! (`CHUNK_SIZE`), which cuts its chunks from the records' letters after the encoding, drops them
//! after its chunks. The records are packed whatever the size of the file and from standard input.
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNQPACK", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-qdrop-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

fn letters(path: &Path) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// The queries in `dir`: five records from LC738874 (2,000 letters each, the third without
/// letters, the fourth 50 `N`s, which has no valid context), so `BATCH_SIZE=3000` makes two
/// batches with the two invalid queries in the second; and a database of LC738868 made with
/// `losat makeblastdb`.
fn inputs(dir: &TestDir) -> (PathBuf, PathBuf) {
    let genome = letters(&fasta("LC738874.fasta"));
    let query = dir.0.join("queries.fna");
    std::fs::write(
        &query,
        format!(
            ">q1 first\n{}\n>q2\n{}\n>empty\n>n50\n{}\n>q5 last\n{}\n",
            &genome[..2_000],
            &genome[2_000..4_000],
            "N".repeat(50),
            &genome[6_000..8_000]
        ),
    )
    .expect("write the queries");
    let db = dir.0.join("lc738868");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(fasta("LC738868.fasta"))
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, db)
}

/// The queries of the packing test in `dir`, written 60 letters a line as a FASTA file is: a record
/// with a lower-case stretch, one with runs of `N`, IUPAC letters and `U`s (read as `T`), one without
/// letters, one of other letters only (`NR` repeated, which keeps its letters), one all lower case
/// and a plain one; 7,700 letters, so `BATCH_SIZE=3000` makes several batches of a file of more
/// than two.
fn packing_queries(dir: &TestDir) -> PathBuf {
    let genome = letters(&fasta("LC738874.fasta"));
    let lines = |letters: &str| -> String {
        let bytes = letters.as_bytes();
        bytes
            .chunks(60)
            .map(|line| format!("{}\n", std::str::from_utf8(line).expect("ASCII")))
            .collect()
    };
    let soft = format!(
        "{}{}{}",
        &genome[..500],
        genome[500..800].to_ascii_lowercase(),
        &genome[800..2_000]
    );
    let mut ambiguous = genome[2_000..3_500].to_string();
    ambiguous.replace_range(300..310, "NNNNNNNNNN");
    ambiguous.replace_range(900..904, "RYKM");
    let ambiguous = ambiguous.replacen('T', "U", 40);
    let query = dir.0.join("packing.fna");
    std::fs::write(
        &query,
        format!(
            ">soft lower case\n{}>ambiguous\n{}>empty\n>nr\n{}>lower\n{}>plain\n{}",
            lines(&soft),
            lines(&ambiguous),
            lines(&"NR".repeat(200)),
            lines(&genome[4_000..6_000].to_ascii_lowercase()),
            lines(&genome[6_000..7_800])
        ),
    )
    .expect("write the queries");
    query
}

/// One search of `query` against `subjects` (`-subject` or `-db`, the report on standard output)
/// with the given variables as its only `LOSAT_*` variables (and `BATCH_SIZE`, `CHUNK_SIZE`
/// among `variables`). `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Returns the report, the
/// standard error without its `[TIMING]` lines, and the `x_bn_q_drop` counts (`None` without
/// the line).
fn search(
    task: &str,
    query: &Path,
    subjects: (&str, &Path),
    words: &[&str],
    threads: u32,
    variables: &[(&str, &str)],
) -> (Vec<u8>, String, Option<[u64; 3]>) {
    let (report, messages, counts, _) =
        search_raw(task, query, subjects, words, threads, variables);
    (report, messages, counts)
}

/// `search`, with the whole standard error too.
fn search_raw(
    task: &str,
    query: &Path,
    subjects: (&str, &Path),
    words: &[&str],
    threads: u32,
    variables: &[(&str, &str)],
) -> (Vec<u8>, String, Option<[u64; 3]>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if (name_text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || name == "BATCH_SIZE"
            || name == "CHUNK_SIZE"
            || name == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg(subjects.0)
        .arg(subjects.1)
        .args(["-task", task])
        .args(words)
        .args(["-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {} {words:?} {variables:?} threads {threads}: {stderr}",
        subjects.0
    );
    let counts = timing_counts(&stderr, "x_bn_q_drop");
    let messages = stderr
        .lines()
        .filter(|line| !line.starts_with("[TIMING]"))
        .map(|line| format!("{line}\n"))
        .collect();
    (output.stdout, messages, counts, stderr)
}

/// The first three counts of the `[TIMING] <name>:` line of `stderr` (`None` without the line).
fn timing_counts(stderr: &str, name: &str) -> Option<[u64; 3]> {
    stderr.lines().find_map(|line| {
        let counts = line.strip_prefix(&format!("[TIMING] {name}: "))?;
        let mut values = counts.split(' ').map(|field| {
            field
                .split_once('=')
                .and_then(|(_, value)| value.parse().ok())
                .expect("a count")
        });
        Some([(); 3].map(|_| values.next().expect("three counts")))
    })
}

// Every format that shows no query letter drops the letters of all five records (two batches),
// against a subject file and a database, for megablast and blastn; `LOSAT_X_BNQPACK=0` keeps them.
// The reports and the messages (the queries without letters and without a valid context among
// them) are equal, and the same at four threads (but SAM's, whose header names the threads).
#[test]
fn the_reports_without_query_letters_do_not_depend_on_the_drop() {
    let dir = TestDir::new("letterless");
    let (query, db) = inputs(&dir);
    let subject = fasta("LC738868.fasta");
    for task in ["megablast", "blastn"] {
        // blastn reads a query file of twice its batch or more batch by batch (LOSAT_X_QSTREAM,
        // whose records leave with their batch); megablast reads its records before the search.
        let mut base = vec![("BATCH_SIZE", "3000")];
        if task == "blastn" {
            base.push(("LOSAT_X_QSTREAM", "0"));
        }
        for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
            for format in ["6", "7", "8", "9", "12", "17", "17 SQ"] {
                let words = ["-outfmt", format];
                let what = format!("{task} {} -outfmt {format:?}", subjects.0);
                let run = |threads: u32, extra: &[(&'static str, &'static str)]| {
                    let variables: Vec<(&str, &str)> = base.iter().chain(extra).copied().collect();
                    search(task, &query, subjects, &words, threads, &variables)
                };
                let (default, messages, counts) = run(1, &[TIMING]);
                assert!(!default.is_empty(), "{what}: no report");
                let counts =
                    counts.unwrap_or_else(|| panic!("{what}: no letters dropped: {messages}"));
                assert_eq!(counts, [2, 5, 6_050], "{what}: batches, records, letters");
                let (report, off_messages, counts) = run(1, &[TIMING, SWITCH_OFF]);
                assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
                assert_eq!(off_messages, messages, "{what} {SWITCH_OFF:?}");
                assert_eq!(counts, None, "{what} {SWITCH_OFF:?}");
                // SAM's header names the command line (`@PG CL:`), threads too.
                if !format.starts_with("17") {
                    let (report, threads_messages, counts) = run(4, &[TIMING]);
                    assert_eq!(report, default, "{what} threads 4");
                    assert_eq!(threads_messages, messages, "{what} threads 4");
                    assert_eq!(counts, Some([2, 5, 6_050]), "{what} threads 4");
                }
            }
        }
    }
}

// The formats that show the queries' letters, and `-query_loc`, keep them (no line); the reports
// are those of `LOSAT_X_BNQPACK=0`.
#[test]
fn the_reports_with_query_letters_keep_them() {
    let dir = TestDir::new("letters");
    let (query, _) = inputs(&dir);
    let subject = fasta("LC738868.fasta");
    let subjects = ("-subject", subject.as_path());
    let cases: [&[&str]; 4] = [
        &["-outfmt", "0"],
        &["-outfmt", "5"],
        &["-outfmt", "17 SQ SR"],
        &["-outfmt", "6", "-query_loc", "101-1900"],
    ];
    for words in cases {
        let what = format!("megablast {words:?}");
        let (default, messages, counts) =
            search("megablast", &query, subjects, words, 1, &[TIMING]);
        assert!(!default.is_empty(), "{what}: no report");
        assert_eq!(counts, None, "{what}: letters dropped");
        let (report, off_messages, _) = search(
            "megablast",
            &query,
            subjects,
            words,
            1,
            &[TIMING, SWITCH_OFF],
        );
        assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
        assert_eq!(off_messages, messages, "{what} {SWITCH_OFF:?}");
    }
}

// A batch split into query chunks (`BATCH_SIZE=10000` and `CHUNK_SIZE=2500`: the 6,000 letters of
// three of the records in one batch of two chunks, a file too short for LOSAT_X_QSTREAM) cuts its
// chunks from the records' letters after the encoding: with outfmt 6 it drops the letters after
// its chunks (one batch, three records), with outfmt 0 it keeps them; the reports are those of
// `LOSAT_X_BNQPACK=0`.
#[test]
fn a_batch_of_query_chunks_drops_its_letters_after_them() {
    let dir = TestDir::new("chunks");
    let (_, db) = inputs(&dir);
    let genome = letters(&fasta("LC738874.fasta"));
    let query = dir.0.join("chunked.fna");
    std::fs::write(
        &query,
        format!(
            ">q1\n{}\n>q2\n{}\n>q3\n{}\n",
            &genome[..2_000],
            &genome[2_000..4_000],
            &genome[6_000..8_000]
        ),
    )
    .expect("write the queries");
    let split = [("BATCH_SIZE", "10000"), ("CHUNK_SIZE", "2500")];
    for format in ["0", "6"] {
        let words = ["-outfmt", format];
        let what = format!("blastn -db -outfmt {format}");
        let subjects = ("-db", db.as_path());
        let with = |extra: &[(&'static str, &'static str)]| {
            let variables: Vec<(&str, &str)> = split.iter().chain(extra).copied().collect();
            search("blastn", &query, subjects, &words, 1, &variables)
        };
        let (default, messages, counts) = with(&[TIMING]);
        assert!(!default.is_empty(), "{what}: no report");
        let expected = (format == "6").then_some([1, 3, 6_000]);
        assert_eq!(
            counts, expected,
            "{what}: batches, records, letters: {messages}"
        );
        let (report, off_messages, counts) = with(&[TIMING, SWITCH_OFF]);
        assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
        assert_eq!(off_messages, messages, "{what} {SWITCH_OFF:?}");
        assert_eq!(counts, None, "{what} {SWITCH_OFF:?}");
    }
}

// The letters of the records read before the search packed as they are read (`LOSAT_X_BNQPACK`,
// `src/sequence/packed_letters.rs`), where the store holds them (`BATCH_SIZE=3000`): lower-case
// stretches, runs of `N`, IUPAC letters, `U`s, a record without letters and lines of 60 letters
// come back as the reader gave them. The reports (`-lcase_masking` too) and the messages are those
// of `LOSAT_X_BNQPACK=0`, against a subject file and a database, for megablast and blastn, at one
// and four threads (the pool makes a batch's letters again). Four records pack (the one without
// letters and the one of other letters only keep theirs); without `BATCH_SIZE` the file is shorter
// than two batches and its records pack all the same (no minimum since 2026-11-01).
#[test]
fn the_letters_packed_as_read_give_the_same_reports() {
    let dir = TestDir::new("packing");
    let (_, db) = inputs(&dir);
    let query = packing_queries(&dir);
    let subject = fasta("LC738868.fasta");
    for task in ["megablast", "blastn"] {
        let mut base = vec![("BATCH_SIZE", "3000")];
        if task == "blastn" {
            base.push(("LOSAT_X_QSTREAM", "0"));
        }
        for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
            for words in [&["-outfmt", "6"][..], &["-outfmt", "7", "-lcase_masking"]] {
                let what = format!("{task} {} {words:?}", subjects.0);
                let run = |threads: u32, extra: &[(&'static str, &'static str)]| {
                    let variables: Vec<(&str, &str)> = base.iter().chain(extra).copied().collect();
                    search_raw(task, &query, subjects, words, threads, &variables)
                };
                let (default, messages, _, stderr) = run(1, &[TIMING]);
                assert!(!default.is_empty(), "{what}: no report");
                let packed = timing_counts(&stderr, "x_bn_q_pack")
                    .unwrap_or_else(|| panic!("{what}: nothing packed: {stderr}"));
                assert_eq!(packed[..2], [4, 7_300], "{what}: records, letters");
                let (report, off_messages, _, stderr) = run(1, &[TIMING, SWITCH_OFF]);
                assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
                assert_eq!(off_messages, messages, "{what} {SWITCH_OFF:?}");
                assert_eq!(
                    timing_counts(&stderr, "x_bn_q_pack"),
                    None,
                    "{what} {SWITCH_OFF:?}"
                );
                let (report, threads_messages, _, stderr) = run(4, &[TIMING]);
                assert_eq!(report, default, "{what} threads 4");
                assert_eq!(threads_messages, messages, "{what} threads 4");
                assert!(
                    timing_counts(&stderr, "x_bn_q_pack").is_some(),
                    "{what} threads 4"
                );
            }
        }
        // No `BATCH_SIZE`: the file is far below two batches, and its records pack.
        let one_batch = |extra: &[(&'static str, &'static str)]| {
            let variables: Vec<(&str, &str)> = [TIMING, ("LOSAT_X_QSTREAM", "0")]
                .iter()
                .chain(extra)
                .copied()
                .collect();
            search_raw(
                task,
                &query,
                ("-subject", subject.as_path()),
                &["-outfmt", "6"],
                1,
                &variables,
            )
        };
        let (default, messages, _, stderr) = one_batch(&[]);
        let packed = timing_counts(&stderr, "x_bn_q_pack")
            .unwrap_or_else(|| panic!("{task}: one batch, nothing packed: {stderr}"));
        assert!(packed[0] > 0, "{task}: one batch packed {packed:?}");
        let (report, off_messages, _, stderr) = one_batch(&[SWITCH_OFF]);
        assert_eq!(report, default, "{task}: one batch {SWITCH_OFF:?}");
        assert_eq!(off_messages, messages, "{task}: one batch {SWITCH_OFF:?}");
        assert_eq!(timing_counts(&stderr, "x_bn_q_pack"), None, "{task}");
    }
}

// Standard input (`-query -`, outfmt 6) is packed: the `x_bn_q_pack` line counts records, and the
// report is that of `LOSAT_X_BNQPACK=0`.
#[test]
fn standard_input_is_packed() {
    use std::io::Write;
    use std::process::Stdio;
    let dir = TestDir::new("stdin");
    let query = packing_queries(&dir);
    let text = std::fs::read(&query).expect("read the queries");
    let subject = fasta("LC738868.fasta");
    let run = |extra: &[(&str, &str)]| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            let name_text = name.to_string_lossy();
            if (name_text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
                || matches!(
                    &*name_text,
                    "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE"
                )
            {
                command.env_remove(name);
            }
        }
        let mut child = command
            .args(["blastn", "-task", "megablast", "-query", "-", "-subject"])
            .arg(&subject)
            .args(["-outfmt", "6", "-num_threads", "1"])
            .envs(extra.iter().copied())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("run LOSAT");
        child
            .stdin
            .take()
            .expect("standard input")
            .write_all(&text)
            .expect("write the queries");
        let output = child.wait_with_output().expect("wait for LOSAT");
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(output.status.success(), "{extra:?}: {stderr}");
        (output.stdout, stderr)
    };
    let (default, stderr) = run(&[TIMING]);
    assert!(!default.is_empty(), "no report");
    let packed = timing_counts(&stderr, "x_bn_q_pack")
        .unwrap_or_else(|| panic!("standard input, nothing packed: {stderr}"));
    assert!(packed[0] > 0, "standard input: {packed:?}");
    let (report, off_stderr) = run(&[TIMING, SWITCH_OFF]);
    assert!(report == default, "standard input {SWITCH_OFF:?}");
    assert_eq!(timing_counts(&off_stderr, "x_bn_q_pack"), None);
    // The reader's messages and every other line of standard error are the same.
    let messages = |stderr: &str| -> Vec<String> {
        stderr
            .lines()
            .filter(|line| !line.starts_with("[TIMING]"))
            .map(str::to_owned)
            .collect()
    };
    assert_eq!(messages(&stderr), messages(&off_stderr), "standard input");
}
