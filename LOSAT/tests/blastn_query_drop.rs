//! A query batch's letters dropped once the batch has encoded them (`LOSAT_X_BNQDROP`,
//! `src/algorithm/blastn/blast_engine/x_query_drop.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`. It changes when the queries' letters are freed, not what the search computes: after the
//! encoding the search reads the records' lengths only, and its rule keeps the letters where a
//! report shows them. The property under test is that the reports and the messages do not depend
//! on it and that the default run took the new path where the rule says so: with `LOSAT_TIMING`
//! the default run prints `[TIMING] x_bn_q_drop: batches= records= letters=`, and `=0` prints no
//! such line. The formats that show no query letter (outfmt 6, 7, 8, 9, 12 and SAM without `SR`)
//! drop the letters of every record, in a search of several batches (`BATCH_SIZE`) with a query
//! without letters and a query without a valid context among them; the formats that show them
//! (outfmt 0, 5, SAM `SQ SR`) and `-query_loc` keep them; a batch split into query chunks
//! (`CHUNK_SIZE`) drops its chunks' copies with every format. That those bytes equal NCBI's is the
//! business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNQDROP", "0");
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
/// `losat makedb`.
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
        .arg("makedb")
        .arg("-in")
        .arg(fasta("LC738868.fasta"))
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, db)
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
    let counts = stderr.lines().find_map(|line| {
        let counts = line.strip_prefix("[TIMING] x_bn_q_drop: ")?;
        let mut values = counts.split(' ').map(|field| {
            field
                .split_once('=')
                .and_then(|(_, value)| value.parse().ok())
                .expect("a count")
        });
        Some([(); 3].map(|_| values.next().expect("three counts")))
    });
    let messages = stderr
        .lines()
        .filter(|line| !line.starts_with("[TIMING]"))
        .map(|line| format!("{line}\n"))
        .collect();
    (output.stdout, messages, counts)
}

// Every format that shows no query letter drops the letters of all five records (two batches),
// against a subject file and a database, for megablast and blastn; `LOSAT_X_BNQDROP=0` keeps them.
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
// are those of `LOSAT_X_BNQDROP=0`.
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
// three of the records in one batch of two chunks) drops its chunks' copies of the letters after
// their encoding, also where a report shows the records' letters (outfmt 0, which keeps the
// records'); the reports are those of `LOSAT_X_BNQDROP=0`.
#[test]
fn the_query_chunks_drop_their_letters() {
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
        let [batches, records, _] =
            counts.unwrap_or_else(|| panic!("{what}: no letters dropped: {messages}"));
        // Two chunks, and with outfmt 6 the batch's own three records.
        let (whole_batches, whole_records) = if format == "6" { (1, 3) } else { (0, 0) };
        assert_eq!(batches, 2 + whole_batches, "{what}");
        assert!(records > 3 + whole_records, "{what}: {records} records");
        let (report, off_messages, counts) = with(&[TIMING, SWITCH_OFF]);
        assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
        assert_eq!(off_messages, messages, "{what} {SWITCH_OFF:?}");
        assert_eq!(counts, None, "{what} {SWITCH_OFF:?}");
    }
}

// A query on standard input has no file size for LOSAT_X_QBATCHPAR's default budget, but
// `LOSAT_X_QBATCHPARSLOTS` sets the slots anyway: the slots may then search the batches of a
// `-subject` search at four threads, so the records stay where the slots read them (no line), and
// the search ends with the report of `LOSAT_X_BNQDROP=0`. (A slot that waited on a batch the loop
// had taken out of the store would never end: the run has a deadline.)
#[test]
fn standard_input_with_slots_keeps_the_letters() {
    use std::io::{Read, Write};
    use std::process::Stdio;
    let genome = letters(&fasta("LC738874.fasta"));
    let queries: String = (0..4)
        .map(|index| {
            format!(
                ">s{index}\n{}\n",
                &genome[2_000 * index..2_000 * (index + 1)]
            )
        })
        .collect();
    let subject = fasta("LC738868.fasta");
    let run = |switches: &[(&str, &str)]| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            if (name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
                || name == "BATCH_SIZE"
            {
                command.env_remove(name);
            }
        }
        let mut child = command
            .args(["blastn", "-task", "blastn", "-query", "-", "-subject"])
            .arg(&subject)
            .args(["-outfmt", "6", "-num_threads", "4"])
            .envs([
                ("BATCH_SIZE", "2000"),
                ("LOSAT_X_QBATCHPARSLOTS", "2"),
                TIMING,
            ])
            .envs(switches.iter().copied())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("run LOSAT");
        let mut stdin = child.stdin.take().expect("standard input");
        let text = queries.clone();
        let writer = std::thread::spawn(move || stdin.write_all(text.as_bytes()));
        let mut stdout = child.stdout.take().expect("standard output");
        let reader = std::thread::spawn(move || {
            let mut bytes = Vec::new();
            stdout.read_to_end(&mut bytes).map(|_| bytes)
        });
        let mut stderr = child.stderr.take().expect("standard error");
        let errors = std::thread::spawn(move || {
            let mut text = String::new();
            stderr.read_to_string(&mut text).map(|_| text)
        });
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(120);
        let status = loop {
            if let Some(status) = child.try_wait().expect("wait for LOSAT") {
                break status;
            }
            if std::time::Instant::now() > deadline {
                let _ = child.kill();
                panic!("{switches:?}: the search did not end within 120 s");
            }
            std::thread::sleep(std::time::Duration::from_millis(50));
        };
        writer.join().expect("writer").expect("write the queries");
        let report = reader.join().expect("reader").expect("read the report");
        let messages = errors.join().expect("reader").expect("read standard error");
        assert!(status.success(), "{switches:?}: {messages}");
        (report, messages)
    };
    let (default, messages) = run(&[]);
    assert!(!default.is_empty(), "no report");
    assert!(
        !messages.contains("x_bn_q_drop"),
        "letters dropped where the slots may read them: {messages}"
    );
    let (report, _) = run(&[SWITCH_OFF]);
    assert_eq!(report, default, "{SWITCH_OFF:?}");
}
