//! The reports of each query batch written when the batch loop takes the batch
//! (`LOSAT_X_RPTSTREAM`, `src/algorithm/blastn/blast_engine/x_rpt_stream.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses when LOSAT
//! writes the same reports: each batch's when the loop takes it, or all of them after the loop.
//! So the property under test is that nothing a run that LOSAT does not refuse writes depends on
//! it: for megablast, blastn
//! and dc-megablast, outfmt 0, 6 and 7, at 1 and 4 threads, with small query batches
//! (`BATCH_SIZE`), the searches without the switch and with it give the same exit status,
//! standard output (or `-out` file) and standard error, and the `[X_STATS] RPTSTREAM` line (the arm
//! with `LOSAT_X_STATS`; its lines are taken out of the comparison) shows that the loop wrote the
//! batches. The same holds with `-out` files, with `-query_loc`, `-subject_loc` and `-lcase_masking`, with `-max_target_seqs`,
//! `-max_hsps` and `-subject_besthit`, with the warnings of the reader and of invalid queries,
//! on the error paths that come after the reports of earlier batches (a reader error, a batch of
//! queries without letters, a header file that cannot be read), and with `-db`. Where the switch
//! must not act (a search space that could overflow; scoring options without Karlin-Altschul
//! values; a negative `CHUNK_SIZE`; the debug output) the `[X_STATS]` line names why and the bytes
//! are the same too. The tabular reports with a field list or a delimiter and outfmt 10 and 20 are
//! streamed alike. Where LOSAT refuses a sequence that a later batch reports (a stored title,
//! ambiguity data, a parsed Seq-id, a taxonomy record, a subject's id), the switch writes the
//! reports of the batches before and then the same refusal (drop-in DECISIONS 44): those reports
//! begin the reports of the same search not refused. That the bytes equal NCBI's is the business
//! of the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_RPTSTREAM", "1");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const TASKS: [&str; 3] = ["megablast", "blastn", "dc-megablast"];

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the first record of a FASTA file.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A directory of this test in the temporary directory, removed when dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat_rptstream_{}_{}_{tag}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        TempDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A FASTA record with `title` and `seq` in lines of 61 letters.
fn record(text: &mut Vec<u8>, title: &[u8], seq: &[u8]) {
    text.push(b'>');
    text.extend_from_slice(title);
    text.push(b'\n');
    for line in seq.chunks(61) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
}

/// 16 pieces of LC738874 of 900 to 4,500 letters (every fourth with a run of `N` and lower case
/// letters), one of them with a line of digits (the reader's warning), then a record without
/// letters, one of `N`s only (the warning of an invalid query) and one of two letters, `tail`
/// appended: with a batch size of 4,000 letters about a dozen batches.
fn queries(dir: &TempDir, name: &str, tail: &[u8]) -> PathBuf {
    let genome = letters("LC738874.fasta");
    let mut text = Vec::new();
    for i in 0..16usize {
        let start = i * 11_000 + 137 * i;
        let len = 900 + (i * 1_531) % 3_600;
        let mut seq = genome[start..start + len].to_vec();
        if i % 4 == 0 {
            let at = len / 2;
            seq.splice(at..at, std::iter::repeat_n(b'N', 25));
            for byte in &mut seq[at + 40..at + 140] {
                *byte = byte.to_ascii_lowercase();
            }
        }
        record(
            &mut text,
            format!("piece{i} title words {i}").as_bytes(),
            &seq,
        );
        if i == 6 {
            // A line the reader keeps with a warning for its digits.
            text.extend_from_slice(b"ACGT12ACGT\n");
        }
        if i == 9 {
            text.extend_from_slice(b">without letters\n");
            record(&mut text, b"onlyN", &[b'N'; 600]);
        }
    }
    record(&mut text, b"two", b"AC");
    text.extend_from_slice(tail);
    let path = dir.path(name);
    std::fs::write(&path, text).expect("write the queries");
    path
}

/// One run's exit status, report (standard output, or the `-out` file) and standard error.
#[derive(Debug, PartialEq, Eq)]
struct Run {
    status: Option<i32>,
    report: Vec<u8>,
    stderr: Vec<u8>,
}

/// One search with `args`, the given variables as its only `LOSAT_*` variables and `BATCH_SIZE`
/// 4,000 (`CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` removed, unless in `switches`).
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing (it is set to a file of this test when the environment has
/// none). With `out`, the report is that file's bytes (standard output must then be empty).
fn blastn(args: &[&str], threads: u32, switches: &[(&str, &str)], out: Option<&Path>) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
            || text.starts_with("BLEMIR_")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!("losat_rptstream_{}_worker.log", std::process::id())),
        );
    }
    command
        .arg("blastn")
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "4000")
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(switches.iter().copied());
    if let Some(out) = out {
        let _ = std::fs::remove_file(out);
        command.arg("-out").arg(out);
    }
    let output = command.output().expect("run LOSAT");
    let report = match out {
        Some(out) => {
            assert!(output.stdout.is_empty(), "{args:?}: output beside -out");
            std::fs::read(out).unwrap_or_default()
        }
        None => output.stdout,
    };
    Run {
        status: output.status.code(),
        report,
        stderr: output.stderr,
    }
}

/// The `[X_STATS]` lines taken out of `run`'s standard error: the numbers of the line
/// `[X_STATS] RPTSTREAM batches=B queries=Q off=O ...`, (B, O).
fn take_stats_lines(run: &mut Run) -> (u64, String) {
    let text = String::from_utf8_lossy(&run.stderr).into_owned();
    let lines: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| line.starts_with("[X_STATS] RPTSTREAM "))
        .collect();
    assert_eq!(lines.len(), 1, "one RPTSTREAM line:\n{text}");
    let line = lines[0].trim_end();
    let field = |name: &str| -> String {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
            .to_string()
    };
    let numbers = (field("batches=").parse().unwrap(), field("off="));
    let rest: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_STATS] "))
        .collect();
    run.stderr = rest.concat().into_bytes();
    numbers
}

/// The run without the switch and with it (with `extra` switches in both), and with it and
/// `LOSAT_X_STATS`: the same status, report and standard error (the counters' lines taken out of
/// the last). Returns the off run and the counters (batches written by the loop, why the switch
/// did not act: `none` when it did).
fn same_two(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
) -> (Run, (u64, String)) {
    let off = blastn(args, threads, extra, out);
    let mut with = extra.to_vec();
    with.push(SWITCH);
    let on = blastn(args, threads, &with, out);
    assert!(
        on == off,
        "{what} threads {threads}: on differs\n{on:?}\n{off:?}"
    );
    with.push(STATS);
    let mut counted = blastn(args, threads, &with, out);
    let numbers = take_stats_lines(&mut counted);
    assert!(
        counted == off,
        "{what} threads {threads}: on with statistics differs\n{counted:?}\n{off:?}"
    );
    (off, numbers)
}

/// A search that LOSAT refuses in a batch after the first (drop-in DECISIONS 44): without the
/// switch the refusal comes before every report; with it (and with `LOSAT_X_STATS`, the counters'
/// line taken out) after the reports of the batches before, with the same exit status and
/// standard error ending as without it (the batches' query warnings may come first). Returns the
/// off run, the run with the switch and its counters (the switch acts: `none`).
fn refused_late(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
) -> (Run, Run, (u64, String)) {
    let off = blastn(args, threads, extra, None);
    let mut with = extra.to_vec();
    with.push(SWITCH);
    let on = blastn(args, threads, &with, None);
    with.push(STATS);
    let mut counted = blastn(args, threads, &with, None);
    let numbers = take_stats_lines(&mut counted);
    assert!(
        counted == on,
        "{what} threads {threads}: on with statistics differs\n{counted:?}\n{on:?}"
    );
    assert_ne!(off.status, Some(0), "{what}: the refusal");
    assert_eq!(on.status, off.status, "{what}: the exit status");
    assert!(
        on.stderr.ends_with(&off.stderr),
        "{what}: standard error\n{}\n{}",
        String::from_utf8_lossy(&on.stderr),
        String::from_utf8_lossy(&off.stderr)
    );
    assert_eq!(numbers.1, "none", "{what}: the switch acts");
    (off, on, numbers)
}

/// `on`'s report holds the reports of the batches before the refused one: it begins `reference`
/// (the report of a search of the same queries that is not refused), is shorter and not empty,
/// and holds `queries` queries' sections where `marker` begins each (`None`: not counted). Its
/// standard error is that of the off run `off` after the beginning of `reference_stderr` (the
/// query warnings of the batches before, which a refusal in the search of the batch writes
/// without the switch too: then none).
fn reports_before(
    what: &str,
    off: &Run,
    on: &Run,
    reference: &[u8],
    reference_stderr: &[u8],
    marker: Option<(&str, usize)>,
) {
    assert!(
        !on.report.is_empty()
            && on.report.len() < reference.len()
            && reference.starts_with(&on.report),
        "{what}: the reports of the batches before\n{}\n{}",
        String::from_utf8_lossy(&on.report),
        String::from_utf8_lossy(reference)
    );
    let head = &on.stderr[..on.stderr.len() - off.stderr.len()];
    assert!(
        reference_stderr.starts_with(head),
        "{what}: the warnings of the batches before\n{}\n{}",
        String::from_utf8_lossy(head),
        String::from_utf8_lossy(reference_stderr)
    );
    if let Some((marker, queries)) = marker {
        let text = String::from_utf8_lossy(&on.report);
        assert_eq!(text.matches(marker).count(), queries, "{what}: {text}");
    }
}

/// Every task, outfmt 0, 6 and 7, 1 and 4 threads: the off bytes, with every batch written by
/// the loop and the warnings of the reader and of the invalid queries on standard error.
#[test]
fn every_task_format_and_thread_count_gives_the_off_bytes() {
    let dir = TempDir::new("formats");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    for task in TASKS {
        for outfmt in ["0", "6", "7"] {
            for threads in [1, 4] {
                let args = [
                    "-task",
                    task,
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    outfmt,
                ];
                let what = format!("{task} outfmt {outfmt}");
                let (off, (batches, reason)) = same_two(&what, &args, threads, &[], None);
                assert_eq!(off.status, Some(0), "{what}");
                assert!(!off.report.is_empty(), "{what}: no report");
                let stderr = String::from_utf8_lossy(&off.stderr);
                assert!(
                    stderr.contains("FASTA-Reader: Ignoring invalid residues")
                        && stderr.contains("Sequence contains no data"),
                    "{what}: {stderr}"
                );
                assert_eq!(reason, "none", "{what}");
                assert!(batches >= 8, "{what}: {batches} batches");
            }
        }
    }
}

/// `-out` files: the off bytes. Each task into an `-out` file with one format, so that each format
/// goes once into a file (every task and format at 1 and 4 threads:
/// `every_task_format_and_thread_count_gives_the_off_bytes`).
#[test]
fn out_files_give_the_off_bytes() {
    let dir = TempDir::new("out");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let out = dir.path("report.txt");
    for (task, out_format) in [("megablast", "0"), ("blastn", "6"), ("dc-megablast", "7")] {
        let args = [
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            out_format,
        ];
        let what = format!("{task} outfmt {out_format} -out");
        let (off, (batches, reason)) = same_two(&what, &args, 4, &[], Some(&out));
        assert!(
            !off.report.is_empty() && reason == "none" && batches >= 8,
            "{what}"
        );
    }
}

/// `-query_loc`, `-subject_loc` with `-lcase_masking` (a subject with lower case runs, whose
/// masks the switch makes once from the records as read), and `-max_target_seqs`, `-max_hsps`
/// and `-subject_besthit`: the off bytes, each case at one thread count (the switch writes the
/// same way at both with `-subject`).
#[test]
fn ranges_and_hit_list_options_give_the_off_bytes() {
    let dir = TempDir::new("ranges");
    let query = queries(&dir, "q.fa", b"");
    // LC738868 with lower case letters in runs of 400 (one run in three), so that
    // `-lcase_masking` has subject masks to show, also inside the `-subject_loc` range.
    let mut lower = letters("LC738868.fasta");
    for (position, letter) in lower.iter_mut().enumerate() {
        if (position / 400) % 3 == 1 {
            *letter = letter.to_ascii_lowercase();
        }
    }
    let mut text = Vec::new();
    record(&mut text, b"LC738868 in lower and upper case", &lower);
    let subject = dir.path("lower.fa");
    std::fs::write(&subject, text).unwrap();
    let cases: [(&str, u32, &[&str]); 4] = [
        ("blastn", 1, &["-outfmt", "0", "-query_loc", "100-1500"]),
        ("megablast", 4, &["-outfmt", "6", "-query_loc", "100-1500"]),
        (
            "dc-megablast",
            4,
            &[
                "-outfmt",
                "0",
                "-subject_loc",
                "2000-30000",
                "-lcase_masking",
            ],
        ),
        (
            "blastn",
            1,
            &[
                "-outfmt",
                "7",
                "-max_target_seqs",
                "1",
                "-max_hsps",
                "2",
                "-subject_besthit",
            ],
        ),
    ];
    for (task, threads, extra) in cases {
        let mut args = vec![
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        let what = format!("{task} {extra:?}");
        let (off, (batches, reason)) = same_two(&what, &args, threads, &[], None);
        assert_eq!(off.status, Some(0), "{what}");
        assert!(
            !off.report.is_empty() && reason == "none" && batches >= 8,
            "{what}"
        );
    }
}

/// The errors of a later batch that come after the reports of the batches before (a reader
/// error; a batch of queries without letters): the off bytes, status and standard error, each
/// task at one thread count for one query and at the other for the other.
#[test]
fn errors_after_earlier_batches_give_the_off_bytes() {
    let dir = TempDir::new("errors");
    let subject = fasta("LC738868.fasta");
    let reader_error = queries(&dir, "reader.fa", b">bad\n@@@@\n");
    // A record that ends its batch (4,500 letters), then a batch of records without letters.
    let mut tail = Vec::new();
    record(
        &mut tail,
        b"long",
        &letters("LC738874.fasta")[200_000..204_500],
    );
    tail.extend_from_slice(b">e1\n>e2\n>e3\n");
    let empty_batch = queries(&dir, "empty.fa", &tail);
    for (query, status, error, threads) in [
        (
            &reader_error,
            1,
            "BLAST query error: CFastaReader: Near line",
            [1, 4, 1],
        ),
        (
            &empty_batch,
            3,
            "BLAST engine error: Warning: Sequence contains no data",
            [4, 1, 4],
        ),
    ] {
        let tasks = [("megablast", "7"), ("blastn", "0"), ("dc-megablast", "6")];
        for ((task, outfmt), threads) in tasks.into_iter().zip(threads) {
            let args = [
                "-task",
                task,
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                outfmt,
            ];
            let what = format!("{task} outfmt {outfmt} {query:?}");
            let (off, (batches, reason)) = same_two(&what, &args, threads, &[], None);
            assert_eq!(off.status, Some(status), "{what}");
            assert!(
                String::from_utf8_lossy(&off.stderr).contains(error),
                "{what}: {}",
                String::from_utf8_lossy(&off.stderr)
            );
            assert!(!off.report.is_empty(), "{what}");
            assert!(reason == "none" && batches >= 8, "{what}");
        }
    }
}

/// A database (`-db`, made by `losat makeblastdb` from four genomes): the off bytes, every format,
/// each (task, format) at one thread count.
#[test]
fn a_database_search_gives_the_off_bytes() {
    let dir = TempDir::new("db");
    let query = queries(&dir, "q.fa", b"");
    let mut set = Vec::new();
    for name in ["LC738868", "LC738869", "LC738870", "LC738874"] {
        set.extend(std::fs::read(fasta(&format!("{name}.fasta"))).unwrap());
    }
    std::fs::write(dir.path("set.fna"), set).unwrap();
    let db = dir.path("set").to_str().unwrap().to_string();
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.path("set.fna"))
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(output.status.success(), "makeblastdb");
    for (task, outfmt, threads) in [
        ("megablast", "0", 1),
        ("blastn", "6", 4),
        ("dc-megablast", "7", 1),
        ("megablast", "7", 4),
    ] {
        let args = [
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-db",
            &db,
            "-outfmt",
            outfmt,
        ];
        let what = format!("-db {task} outfmt {outfmt}");
        let (off, (batches, reason)) = same_two(&what, &args, threads, &[], None);
        assert_eq!(off.status, Some(0), "{what}");
        assert!(
            !off.report.is_empty() && reason == "none" && batches >= 8,
            "{what}"
        );
    }
}

/// `losat makeblastdb` of LC738868 and of a part of LC738874 (150,000-190,000, with three runs of `N`, so
/// with ambiguity data) that only the queries of the last batches hit; `patch` changes the files.
fn database(dir: &TempDir, name: &str, patch: impl FnOnce(&str)) -> String {
    let mut part = letters("LC738874.fasta")[150_000..190_000].to_vec();
    for at in [5_000, 17_000, 31_000] {
        part[at..at + 30].fill(b'N');
    }
    let mut text = Vec::new();
    record(&mut text, b"LC738868 first", &letters("LC738868.fasta"));
    record(&mut text, b"LC738874part second", &part);
    let fasta = dir.path(&format!("{name}.fna"));
    std::fs::write(&fasta, text).unwrap();
    let db = dir.path(name).to_str().unwrap().to_string();
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(&fasta)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(output.status.success(), "makeblastdb");
    patch(&db);
    db
}

/// The same database search, intact and refused by LOSAT where a later batch reports a sequence:
/// the second sequence's title holding a byte outside ' ' ..= '~' in its header file, or its
/// ambiguity count raised by one in the sequence file. Both are refused when a batch reports the
/// sequence (its title) or traces it back (its ambiguity data), which only the last batches do:
/// without the switch the refusal comes before every report, with it after the rows of the
/// batches before (drop-in DECISIONS 44; nothing is read before the search). The intact database
/// streams, with its ambiguity data decoded by the traceback, and so do the refused ones for
/// queries that do not reach the second sequence.
#[test]
fn a_database_search_refused_in_a_later_batch_writes_the_batches_before() {
    let dir = TempDir::new("dbrefuse");
    let query = queries(&dir, "q.fa", b"");
    // The pieces before the twelfth do not reach the second sequence (150,000-190,000).
    let text = std::fs::read(&query).unwrap();
    let cut = text
        .windows(8)
        .position(|window| window == b">piece12")
        .unwrap();
    let early = dir.path("early.fa");
    std::fs::write(&early, &text[..cut]).unwrap();
    let args = |query: &Path, db: &str, outfmt: &str| -> Vec<String> {
        [
            "-task",
            "megablast",
            "-query",
            query.to_str().unwrap(),
            "-db",
            db,
            "-outfmt",
            outfmt,
        ]
        .map(String::from)
        .to_vec()
    };
    let run = |query: &Path, db: &str, outfmt: &str, threads: u32| {
        let args = args(query, db, outfmt);
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        same_two(
            &format!("-db {db} outfmt {outfmt}"),
            &args,
            threads,
            &[],
            None,
        )
    };
    let intact = database(&dir, "intact", |_| {});
    let (mut reference, mut reference_stderr) = (Vec::new(), Vec::new());
    for threads in [1, 4] {
        let (off, (batches, reason)) = run(&query, &intact, "6", threads);
        assert_eq!(off.status, Some(0));
        let rows = String::from_utf8_lossy(&off.report).into_owned();
        assert!(
            rows.contains("\tgnl|BL_ORD_ID|1\t") || rows.contains("LC738874part"),
            "{rows}"
        );
        assert!(batches >= 8 && reason == "none");
        (reference, reference_stderr) = (off.report, off.stderr);
    }
    let title = database(&dir, "title", |db| {
        let path = format!("{db}.nhr");
        let mut header = std::fs::read(&path).unwrap();
        let at = header
            .windows(12)
            .position(|window| window == b"LC738874part")
            .expect("the second title");
        header[at + 3] = 0xE9;
        std::fs::write(&path, header).unwrap();
    });
    let ambiguity = database(&dir, "ambiguity", |db| {
        let volume = LOSAT::seqdb::SeqDb::open(db, LOSAT::seqdb::SeqType::Nucleotide).unwrap();
        let index = volume.volumes()[0].index();
        let (at, next) = index.ambiguity_range(1).unwrap().unwrap();
        let (at, next) = (at as usize, next as usize);
        assert!(at < next, "the second sequence has ambiguity data");
        drop(volume);
        let path = format!("{db}.nsq");
        let mut sequences = std::fs::read(&path).unwrap();
        let count = u32::from_be_bytes(sequences[at..at + 4].try_into().unwrap());
        sequences[at..at + 4].copy_from_slice(&(count + 1).to_be_bytes());
        std::fs::write(&path, sequences).unwrap();
    });
    // A title of spaces only: NCBI makes no tabular subject id of it, which LOSAT refuses where a
    // tabular report shows the sequence; outfmt 0 and 17 show it (and stream, with the off bytes).
    let blank = database(&dir, "blank", |db| {
        let path = format!("{db}.nhr");
        let mut header = std::fs::read(&path).unwrap();
        let title = b"LC738874part second";
        let at = header
            .windows(title.len())
            .position(|window| window == title)
            .expect("the second title");
        header[at..at + title.len()].fill(b' ');
        std::fs::write(&path, header).unwrap();
    });
    for (outfmt, threads) in [("0", 4), ("17", 1)] {
        let (off, (batches, reason)) = run(&query, &blank, outfmt, threads);
        assert_eq!(off.status, Some(0), "blank title, outfmt {outfmt}: {off:?}");
        assert!(
            batches >= 8 && reason == "none",
            "blank title, outfmt {outfmt}"
        );
    }
    for (db, refusal) in [
        (&title, "has a stored title with a byte outside"),
        (&blank, "makes no subject id"),
        (
            &ambiguity,
            "the ambiguity data are shorter than their count",
        ),
    ] {
        for threads in [1, 4] {
            let what = format!("-db {db} threads {threads}");
            let tabular = args(&query, db, "6");
            let tabular: Vec<&str> = tabular.iter().map(String::as_str).collect();
            let (off, on, _) = refused_late(&what, &tabular, threads, &[]);
            assert!(
                String::from_utf8_lossy(&off.stderr).contains(refusal),
                "{what}: {}",
                String::from_utf8_lossy(&off.stderr)
            );
            assert!(off.report.is_empty(), "{what}: off, no row");
            reports_before(&what, &off, &on, &reference, &reference_stderr, None);
            // Queries that do not reach the second sequence: no refusal, the off bytes.
            let (off, (batches, reason)) = run(&early, db, "6", threads);
            assert_eq!(off.status, Some(0), "{what}, early queries: {off:?}");
            assert!(batches >= 2 && reason == "none", "{what}, early queries");
            // Outfmt 17 writes its document only at the end: a refused run writes none, with the
            // switch too (the blank title is not refused there).
            if db != &blank {
                let sam = args(&query, db, "17");
                let sam: Vec<&str> = sam.iter().map(String::as_str).collect();
                let what = format!("{what} outfmt 17");
                let (off, on, _) = refused_late(&what, &sam, threads, &[]);
                assert!(
                    off.report.is_empty() && on.report.is_empty(),
                    "{what}: {on:?}"
                );
            }
        }
    }
}

/// Eight queries of 488 random letters (no hits) and then the records of `query`, in
/// `late_<query's file name>`: in batches of four (`BATCH_SIZE=1500`) only the third and later
/// batches report a sequence. The eight alone are in `nohit.fa` (`no_hits`).
fn hits_late(dir: &TempDir, query: &str) -> PathBuf {
    let mut text = std::fs::read(no_hits(dir)).unwrap();
    text.extend(std::fs::read(query).unwrap());
    let name = Path::new(query).file_name().unwrap().to_str().unwrap();
    let late = dir.path(&format!("late_{name}"));
    std::fs::write(&late, text).unwrap();
    late
}

/// The eight queries of 488 random letters that begin `hits_late`'s queries, in `nohit.fa`.
fn no_hits(dir: &TempDir) -> PathBuf {
    let mut text = Vec::new();
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    for at in 0..8 {
        let letters: Vec<u8> = (0..488)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        record(&mut text, format!("nohit{at}").as_bytes(), &letters);
    }
    let path = dir.path("nohit.fa");
    std::fs::write(&path, text).unwrap();
    path
}

/// A database whose header file is missing (`tests/fixtures/outfmt17/db/snucl_nophr`), hit from the
/// third batch on (eight queries without hits first): NCBI's error of the batch that first reads a
/// definition line (`HeaderReadError`, exit status 2) after the reports of the batches before, the
/// off bytes in every kind of format (until 2026-11-03 the check before the search kept such a
/// database's reports after the loop).
#[test]
fn a_header_file_missing_in_a_later_batch_gives_the_off_bytes() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/outfmt17");
    let dir = TempDir::new("nophr");
    let late = hits_late(&dir, fixtures.join("qnucl_valid.fna").to_str().unwrap());
    let db = fixtures.join("db/snucl_nophr");
    for (outfmt, threads) in [
        ("6", 1),
        ("7", 4),
        ("0", 1),
        ("17", 4),
        ("20 qseqid sseqid stitle", 1),
    ] {
        let args = [
            "-task",
            "megablast",
            "-query",
            late.to_str().unwrap(),
            "-db",
            db.to_str().unwrap(),
            "-outfmt",
            outfmt,
        ];
        let what = format!("snucl_nophr outfmt {outfmt}");
        let (off, (batches, reason)) =
            same_two(&what, &args, threads, &[("BATCH_SIZE", "1500")], None);
        assert_eq!(off.status, Some(2), "{what}: {off:?}");
        assert!(batches >= 2 && reason == "none", "{what}");
    }
}

/// Databases with parsed Seq-ids (line DBID; `tests/db_search/ids`, `make_ids.py`), with the 25
/// records of `nucl.fna` (variants of one sequence: every query hits every sequence) as queries
/// in batches of four (`BATCH_SIZE=1500`): the off bytes, written batch by batch, also where
/// LOSAT refuses a report of the first batch (the refusal before every report in both): two
/// sequences that outfmt 0 joins (`njoined`; the tabular formats do not join them), a sequence
/// that NCBI's formatters do not find (`nrefuse_v5`). Two that the check before the search (until
/// 2026-11-03) took for refusals stream with the off bytes: with `taxdb.bti` found, outfmt 0's
/// titles (a taxonomy database with a wrong magic number, which writes its messages with each
/// batch's reports) and the prefetch of a taxonomy database whose records LOSAT refuses where a row
/// looks a name up (`btd_empty`, `custom_tabular_reports_give_the_off_bytes`). Then eight queries without hits come first, so that only
/// the third and later batches report a sequence: a refusal comes after the reports of the first
/// two batches with the switch (drop-in DECISIONS 44; outfmt 7 writes each query's comment lines
/// and outfmt 0 its section, hits or not), and those are the reports of the eight queries alone.
#[test]
fn a_database_with_parsed_seqids_gives_the_off_bytes() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let ids = |name: &str| {
        tests
            .join("db_search/ids")
            .join(name)
            .to_str()
            .unwrap()
            .to_string()
    };
    let taxdb = |variant: &str| {
        tests
            .join("seqdb/taxdb")
            .join(variant)
            .to_str()
            .unwrap()
            .to_string()
    };
    let query = ids("nucl.fna");
    let (bad_magic, btd_empty) = (taxdb("bad_magic"), taxdb("btd_empty"));
    for (db, outfmt, blastdb, threads, refused) in [
        ("nucl_v5", "0", None, 1, false),
        ("nucl_v5", "6", None, 4, false),
        ("nucl_v5", "17", None, 1, false),
        ("nucl_v5", "8", None, 4, false),
        ("nucl_v5", "11", None, 1, false),
        ("nucl_v4", "0", None, 4, false),
        ("nucl_v5", "7", Some(&bad_magic), 4, false),
        ("nucl_v5", "0", Some(&bad_magic), 1, false),
        ("nucl_v5", "17", Some(&bad_magic), 4, false),
        ("nucl_v5", "7", Some(&btd_empty), 1, false),
        ("njoined", "6", None, 4, false),
        ("njoined", "0", None, 1, true),
        ("nrefuse_v5", "6", None, 4, true),
    ] {
        let path = ids(db);
        let args = [
            "-task",
            "megablast",
            "-query",
            &query,
            "-db",
            &path,
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            "3",
        ];
        let mut extra = vec![("BATCH_SIZE", "1500")];
        if let Some(dir) = blastdb {
            extra.push(("BLASTDB", dir.as_str()));
        }
        let what = format!("{db} outfmt {outfmt} {blastdb:?}");
        let (off, (batches, reason)) = same_two(&what, &args, threads, &extra, None);
        assert_eq!(reason, "none", "{what}");
        if refused {
            assert_ne!(off.status, Some(0), "{what}: the refusal");
            // The first batch's (its report call counts as a batch written, though refused).
            assert!(batches <= 1, "{what}: the first batch's");
        } else {
            assert_eq!(off.status, Some(0), "{what}");
            assert!(batches >= 6 && off.report.len() > 1000, "{what}");
        }
    }
    // Eight queries of 488 random letters (no hits) before the 25 records; the reports of the
    // eight alone (each query's comment lines or section).
    let dir = TempDir::new("parsedlate");
    let late = hits_late(&dir, &query);
    let late = late.to_str().unwrap();
    let nohit = no_hits(&dir);
    let nohit = nohit.to_str().unwrap();
    let marker = |outfmt: &str| {
        if outfmt == "0" {
            "Query= "
        } else {
            "# Query: "
        }
    };
    for (db, outfmt, threads, refused) in [
        ("nucl_v5", "7", 4, false),
        ("nrefuse_v5", "7", 1, true),
        ("njoined", "0", 4, true),
    ] {
        let path = ids(db);
        let mut args = [
            "-task",
            "megablast",
            "-query",
            late,
            "-db",
            &path,
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            "3",
        ];
        let what = format!("{db} outfmt {outfmt}, queries without hits first");
        let extra = [("BATCH_SIZE", "1500")];
        if refused {
            let (off, on, (batches, _)) = refused_late(&what, &args, threads, &extra);
            args[3] = nohit;
            let (eight, _) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(eight.status, Some(0), "{what}: {eight:?}");
            reports_before(
                &what,
                &off,
                &on,
                &eight.report,
                &eight.stderr,
                Some((marker(outfmt), 8)),
            );
            assert!(batches >= 2, "{what}");
        } else {
            let (off, (batches, reason)) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(off.status, Some(0), "{what}");
            assert!(batches >= 8 && reason == "none", "{what}");
        }
    }
    // `nambig_v4` (the sequence of the first record of `nucl.fna` and a copy of it with ambiguity
    // codes), intact, and with the second sequence's ambiguity count raised by one: the traceback
    // of the third batch refuses it, after the reports of the first two.
    let source = PathBuf::from(ids("nambig_v4"));
    let copy = |name: &str| -> String {
        for ext in ["nhr", "nin", "njs", "nog", "nsd", "nsi", "nsq"] {
            std::fs::copy(
                source.with_extension(ext),
                dir.path(&format!("{name}.{ext}")),
            )
            .unwrap();
        }
        dir.path(name).to_str().unwrap().to_string()
    };
    let intact = copy("nambig_intact");
    let raised = copy("nambig_raised");
    {
        let db = LOSAT::seqdb::SeqDb::open(&raised, LOSAT::seqdb::SeqType::Nucleotide).unwrap();
        let (at, next) = db.volumes()[0].index().ambiguity_range(1).unwrap().unwrap();
        let (at, next) = (at as usize, next as usize);
        assert!(at < next, "the second sequence has ambiguity data");
        drop(db);
        let path = format!("{raised}.nsq");
        let mut sequences = std::fs::read(&path).unwrap();
        let count = u32::from_be_bytes(sequences[at..at + 4].try_into().unwrap());
        sequences[at..at + 4].copy_from_slice(&(count + 1).to_be_bytes());
        std::fs::write(&path, sequences).unwrap();
    }
    for (db, threads) in [(&intact, 4), (&raised, 1)] {
        let mut args = [
            "-task",
            "megablast",
            "-query",
            late,
            "-db",
            db,
            "-outfmt",
            "7",
            "-max_target_seqs",
            "3",
        ];
        let what = format!("{db} outfmt 7, queries without hits first");
        let extra = [("BATCH_SIZE", "1500")];
        if db == &intact {
            let (off, (batches, reason)) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(off.status, Some(0), "{what}: {off:?}");
            assert!(batches >= 8 && reason == "none", "{what}");
            continue;
        }
        let (off, on, (batches, _)) = refused_late(&what, &args, threads, &extra);
        assert!(
            String::from_utf8_lossy(&off.stderr).contains("ambiguity data"),
            "{what}: {}",
            String::from_utf8_lossy(&off.stderr)
        );
        args[3] = nohit;
        let (eight, _) = same_two(&what, &args, threads, &extra, None);
        assert_eq!(eight.status, Some(0), "{what}: {eight:?}");
        reports_before(
            &what,
            &off,
            &on,
            &eight.report,
            &eight.stderr,
            Some(("# Query: ", 8)),
        );
        assert!(batches >= 2, "{what}");
    }
    // Seq-align ids that NCBI's formatters do not find (`nord_v5`: `gnl|BL_ORD_ID|3`) or find in
    // another sequence first (`nlook_v5`: `lcl|300`), which the third batch reports
    // (`db_search_manifest/` refuses them in one batch).
    for (db, queries, refusal, threads) in [
        ("nord_v5", "nord_qall.fna", "gnl|BL_ORD_ID|3", 1),
        ("nlook_v5", "nlook_qall.fna", "lcl|300", 4),
    ] {
        let queries = hits_late(&dir, &ids(queries));
        let path = ids(db);
        let mut args = [
            "-task",
            "megablast",
            "-query",
            queries.to_str().unwrap(),
            "-db",
            &path,
            "-outfmt",
            "7",
        ];
        let what = format!("{db} outfmt 7, queries without hits first");
        let extra = [("BATCH_SIZE", "1500")];
        let (off, on, (batches, _)) = refused_late(&what, &args, threads, &extra);
        assert!(
            String::from_utf8_lossy(&off.stderr).contains(refusal),
            "{what}: {}",
            String::from_utf8_lossy(&off.stderr)
        );
        args[3] = nohit;
        let (eight, _) = same_two(&what, &args, threads, &extra, None);
        assert_eq!(eight.status, Some(0), "{what}: {eight:?}");
        reports_before(
            &what,
            &off,
            &on,
            &eight.report,
            &eight.stderr,
            Some(("# Query: ", 8)),
        );
        assert!(batches >= 2, "{what}");
    }
    // `-parse_deflines`: the same queries under ids of their own stream; under their own ids (the
    // database's) the third batch's report refuses them (`refuse_id_of_a_query`), after the
    // reports of the first two with the switch.
    let renamed = dir.path("renamed.fa");
    let mut at = 0;
    let text = std::fs::read_to_string(late).unwrap();
    let lines: Vec<String> = text
        .lines()
        .map(|line| {
            if !line.starts_with('>') {
                return line.to_string();
            }
            at += 1;
            format!(">query{at}")
        })
        .collect();
    std::fs::write(&renamed, lines.join("\n") + "\n").unwrap();
    let renamed = renamed.to_str().unwrap();
    for (queries, outfmt, threads) in [
        (renamed, "6", 4),
        (renamed, "0", 1),
        (renamed, "17", 4),
        (late, "7", 4),
        (late, "0", 1),
    ] {
        let path = ids("nucl_v5");
        let mut args = [
            "-task",
            "megablast",
            "-query",
            queries,
            "-db",
            &path,
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            "3",
            "-parse_deflines",
        ];
        let what = format!("-parse_deflines {queries} outfmt {outfmt}");
        let extra = [("BATCH_SIZE", "1500")];
        if queries == late {
            let (off, on, (batches, _)) = refused_late(&what, &args, threads, &extra);
            args[3] = nohit;
            let (eight, _) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(eight.status, Some(0), "{what}: {eight:?}");
            reports_before(
                &what,
                &off,
                &on,
                &eight.report,
                &eight.stderr,
                Some((marker(outfmt), 8)),
            );
            assert!(batches >= 2, "{what}");
        } else {
            let (off, (batches, reason)) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(off.status, Some(0), "{what}: {off:?}");
            assert!(batches >= 8 && off.report.len() > 1000, "{what}");
            assert_eq!(reason, "none", "{what}");
        }
    }
}

/// The two subjects of `subjects_late`: LC738868, hit from the first batch on, and a part of
/// LC738874 (the queries' genome) titled `title`, which only the queries of the last batches hit.
fn subjects_late(dir: &TempDir, name: &str, title: &[u8]) -> PathBuf {
    let mut subjects = Vec::new();
    record(&mut subjects, b"LC738868 first", &letters("LC738868.fasta"));
    record(
        &mut subjects,
        title,
        &letters("LC738874.fasta")[150_000..190_000],
    );
    let subject = dir.path(name);
    std::fs::write(&subject, subjects).unwrap();
    subject
}

/// A subject whose tabular id LOSAT refuses (`Subject_` and a byte in no encoding NCBI
/// recognises), which a later batch reaches: without the switch the refusal comes before every
/// report, with it after the rows of the batches before (drop-in DECISIONS 44), the rows of the
/// same search with the subject titled otherwise. `outfmt` names the format.
fn a_refused_subject_id_comes_after_the_batches_before(dir: &TempDir, query: &Path, outfmt: &str) {
    let refused = subjects_late(dir, "refused.fa", b"Subject_\x81 undecodable");
    let titled = subjects_late(dir, "titled.fa", b"LC738874part second");
    let mut args = [
        "-task",
        "megablast",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        refused.to_str().unwrap(),
        "-outfmt",
        outfmt,
    ];
    let what = format!("Subject_ id, outfmt {outfmt}");
    let (off, on, _) = refused_late(&what, &args, 4, &[]);
    assert!(
        String::from_utf8_lossy(&off.stderr).contains("'Subject_'"),
        "{what}: {}",
        String::from_utf8_lossy(&off.stderr)
    );
    assert!(off.report.is_empty(), "{what}: off, no row");
    // Outfmt 0 and 17 show the subject (and stream, with the off bytes).
    for (outfmt, threads) in [("0", 1), ("17", 4)].into_iter().filter(|_| outfmt == "6") {
        let mut shown = args;
        shown[7] = outfmt;
        let what = format!("Subject_ id, outfmt {outfmt}");
        let (off, (batches, reason)) = same_two(&what, &shown, threads, &[], None);
        assert_eq!(off.status, Some(0), "{what}: {off:?}");
        assert!(batches >= 8 && reason == "none", "{what}");
    }
    args[5] = titled.to_str().unwrap();
    let (reference, (_, reason)) = same_two(&what, &args, 4, &[], None);
    assert!(reference.status == Some(0) && reason == "none", "{what}");
    reports_before(&what, &off, &on, &reference.report, &reference.stderr, None);
}

/// Where the switch must not act, it does not, and the bytes are the off path's: a search space
/// that could overflow (a `-dbsize` near 2^63); scoring options without Karlin-Altschul values; a
/// negative `CHUNK_SIZE` (whose refusal comes where a batch is split); and the debug output
/// (`LOSAT_DEBUG_COORDS`). A subject whose tabular id LOSAT refuses streams, its refusal after
/// the rows of the batches before.
#[test]
fn where_the_switch_must_not_act_the_counters_name_why() {
    let dir = TempDir::new("off");
    let query = queries(&dir, "q.fa", b"");
    a_refused_subject_id_comes_after_the_batches_before(&dir, &query, "6");
    let subject = fasta("LC738868.fasta");
    let plain = [
        "-task",
        "megablast",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    let mut huge = plain.to_vec();
    huge.extend_from_slice(&["-dbsize", "9000000000000000000"]);
    let (_, (batches, reason)) = same_two("-dbsize", &huge, 1, &[], None);
    assert_eq!((batches, reason.as_str()), (0, "search_space"));
    let mut scoring = plain.to_vec();
    scoring[1] = "blastn";
    scoring.extend_from_slice(&[
        "-reward",
        "2",
        "-penalty",
        "-3",
        "-gapopen",
        "1",
        "-gapextend",
        "1",
    ]);
    let (off, (batches, reason)) = same_two("scoring", &scoring, 4, &[], None);
    assert_eq!(off.status, Some(3), "NCBI's Karlin-Altschul error");
    assert_eq!((batches, reason.as_str()), (0, "scoring"));
    let chunks = [("CHUNK_SIZE", "-5"), ("OVERLAP_CHUNK_SIZE", "-10")];
    let (off, (batches, reason)) = same_two("CHUNK_SIZE", &plain, 4, &chunks, None);
    assert_ne!(off.status, Some(0), "LOSAT's refusal");
    assert_eq!((batches, reason.as_str()), (0, "chunk_size"));
    let (off, (batches, reason)) =
        same_two("debug", &plain, 1, &[("LOSAT_DEBUG_COORDS", "1")], None);
    assert!(!off.report.is_empty());
    assert_eq!((batches, reason.as_str()), (0, "debug"));
}

/// The tabular reports with a field list or a delimiter and outfmt 10 and 20 (NCMEM1), written
/// batch by batch: the off bytes. With `-subject` (every task; with a plural name of the taxonomy
/// database, whose formatter warns for each query without `taxdb.bti`), with `-parse_deflines`,
/// with `-db` (`stitle`, `sseq`, `staxids`, `sall*` and the names: the definition lines as the
/// report decodes them, the taxonomy database's lookups and their texts), and on NCBI's errors in
/// the first batch, where outfmt 20's line of the field names comes as without the switch. Where a
/// later batch's report is refused (a subject's id, a stored title, a header or taxonomy record of
/// a database with parsed Seq-ids) the refusal comes after the rows of the batches before
/// (drop-in DECISIONS 44).
#[test]
fn custom_tabular_reports_give_the_off_bytes() {
    let dir = TempDir::new("custom");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let streamed = |what: &str, off: &Run, batches: u64, reason: &str| {
        assert_eq!(off.status, Some(0), "{what}: {off:?}");
        assert!(!off.report.is_empty(), "{what}");
        assert!(
            reason == "none" && batches >= 6,
            "{what}: {batches} {reason}"
        );
    };
    for (at, (outfmt, parse_deflines)) in [
        (
            "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue \
             bitscore qlen slen sseq",
            false,
        ),
        ("7 std qcovs qcovus sstrand staxids sscinames", false),
        ("6 delim=, std qcovs", false),
        ("10", false),
        ("10 qseqid sseqid sseq btop", false),
        ("20 qseqid sseqid qcovs", false),
        ("20 delim=; qseqid sseqid", false),
        ("6 qseqid sseqid sallseqid sacc stitle", true),
    ]
    .into_iter()
    .enumerate()
    {
        let task = TASKS[at % 3];
        let mut args = vec![
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            outfmt,
        ];
        if parse_deflines {
            args.push("-parse_deflines");
        }
        let what = format!("{task} outfmt {outfmt} -parse_deflines {parse_deflines}");
        let (off, (batches, reason)) = same_two(&what, &args, [1, 4][at % 2], &[], None);
        streamed(&what, &off, batches, &reason);
        if outfmt.contains("sscinames") {
            assert!(
                String::from_utf8_lossy(&off.stderr).contains("Taxonomy name lookup"),
                "{what}"
            );
        }
    }
    // A subject whose tabular id LOSAT refuses, which only the last batches report.
    a_refused_subject_id_comes_after_the_batches_before(&dir, &query, "10 qseqid sseqid stitle");
    // NCBI's errors in the first batch: a reader error, a batch of queries without letters. The
    // one batch streams; outfmt 20's line of the field names comes before the error.
    let reader_error = dir.path("bad.fa");
    std::fs::write(&reader_error, b">bad\n@@@@\n").unwrap();
    let empty = dir.path("empty.fa");
    std::fs::write(&empty, b">e1\n>e2\n").unwrap();
    for (bad, status) in [(&reader_error, 1), (&empty, 3)] {
        let args = [
            "-task",
            "megablast",
            "-query",
            bad.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "20 qseqid sseqid",
        ];
        let what = format!("outfmt 20, {bad:?}");
        let (off, (batches, reason)) = same_two(&what, &args, 1, &[], None);
        assert_eq!(off.status, Some(status), "{what}: {off:?}");
        assert!(
            !off.report.is_empty(),
            "{what}: the line of the field names"
        );
        assert_eq!((batches, reason.as_str()), (0, "none"), "{what}");
    }
    // `-db`: two sequences, the second hit by the last batches only; and the same database with
    // the second title holding a byte outside ' ' ..= '~'.
    let intact = database(&dir, "intact", |_| {});
    let title = database(&dir, "title", |db| {
        let path = format!("{db}.nhr");
        let mut header = std::fs::read(&path).unwrap();
        let at = header
            .windows(12)
            .position(|window| window == b"LC738874part")
            .expect("the second title");
        header[at + 3] = 0xE9;
        std::fs::write(&path, header).unwrap();
    });
    for (outfmt, threads) in [
        ("6 qseqid sseqid stitle sseq staxids", 4),
        ("7 qseqid sallseqid salltitles sscinames qcovs", 1),
        ("10 qseqid sseqid slen", 4),
        ("20 qseqid stitle staxid", 1),
    ] {
        let mut args = [
            "-task",
            "megablast",
            "-query",
            query.to_str().unwrap(),
            "-db",
            &intact,
            "-outfmt",
            outfmt,
        ];
        let what = format!("-db {intact} outfmt {outfmt}");
        let (reference, (batches, reason)) = same_two(&what, &args, threads, &[], None);
        streamed(&what, &reference, batches, &reason);
        // Outfmt 7 names the database in each query's comment lines.
        let reference_stderr = reference.stderr;
        let reference = String::from_utf8(reference.report)
            .unwrap()
            .replace(intact.as_str(), title.as_str());
        args[5] = &title;
        let what = format!("-db {title} outfmt {outfmt}");
        let (off, on, _) = refused_late(&what, &args, threads, &[]);
        assert!(off.report.is_empty(), "{what}: off, no row");
        reports_before(
            &what,
            &off,
            &on,
            reference.as_bytes(),
            &reference_stderr,
            None,
        );
    }
    // Databases with parsed Seq-ids and taxids (`tests/db_search/ids`), with taxonomy databases
    // (`tests/seqdb/taxdb`: `good`; `bad_magic`, whose texts of opening the first row that looks a
    // name up writes; `btd_empty`, whose records LOSAT refuses; `offset_mid_past_end`, whose
    // offsets do not increase: the record of taxid 9606 begins past the end of `taxdb.btd`, so its
    // names are not found and no lookup is refused), the queries that hit only from the third
    // batch on, 25 sequences each (the sequences with taxids among them). A refusal comes after
    // the comment lines of the eight queries without hits (outfmt 7), the report of those alone.
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let ids = |name: &str| {
        tests
            .join("db_search/ids")
            .join(name)
            .to_str()
            .unwrap()
            .to_string()
    };
    let taxdb = |variant: &str| {
        tests
            .join("seqdb/taxdb")
            .join(variant)
            .to_str()
            .unwrap()
            .to_string()
    };
    let late = hits_late(&dir, &ids("nucl.fna"));
    let nohit = no_hits(&dir);
    for (db, outfmt, blastdb, threads, refused) in [
        (
            "nucl_v5",
            "6 qseqid sseqid sallseqid sallacc staxids",
            None,
            4,
            false,
        ),
        (
            "nucl_v5",
            "7 qseqid sseqid staxids sscinames scomnames sblastnames sskingdoms",
            Some("good"),
            1,
            false,
        ),
        (
            "nucl_v5",
            "6 qseqid ssciname sscinames",
            Some("bad_magic"),
            4,
            false,
        ),
        ("nucl_v4", "20 qseqid sseqid staxid", Some("good"), 1, false),
        (
            "nucl_v5",
            "10 qseqid ssciname",
            Some("offset_mid_past_end"),
            1,
            false,
        ),
        (
            "nucl_v5",
            "7 qseqid sseqid sscinames",
            Some("btd_empty"),
            4,
            true,
        ),
        ("nrefuse_v5", "7 qseqid sallseqid", None, 4, true),
    ] {
        let path = ids(db);
        let mut args = [
            "-task",
            "megablast",
            "-query",
            late.to_str().unwrap(),
            "-db",
            &path,
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            "25",
        ];
        let blastdb = blastdb.map(taxdb);
        let mut extra = vec![("BATCH_SIZE", "1500")];
        if let Some(dir) = &blastdb {
            extra.push(("BLASTDB", dir.as_str()));
        }
        let what = format!("{db} outfmt {outfmt} {blastdb:?}");
        if refused {
            let (off, on, (batches, _)) = refused_late(&what, &args, threads, &extra);
            assert!(off.report.is_empty(), "{what}: off, no row");
            args[3] = nohit.to_str().unwrap();
            let (eight, _) = same_two(&what, &args, threads, &extra, None);
            assert_eq!(eight.status, Some(0), "{what}: {eight:?}");
            reports_before(
                &what,
                &off,
                &on,
                &eight.report,
                &eight.stderr,
                Some(("# Query: ", 8)),
            );
            assert!(batches >= 2, "{what}");
        } else {
            let (off, (batches, reason)) = same_two(&what, &args, threads, &extra, None);
            streamed(&what, &off, batches, &reason);
        }
    }
}

/// The refusal of NCBI's taxonomy cache (`TaxdbRun::cache_collision`: `tests/db_search/nz_nc`
/// holds taxid 150 beside 99990, the damaged last entry of `tests/seqdb/taxdb/offset_last_past_end`;
/// line DBID L1c) is asked where NCBI fetches a subject, one query per batch: `qx` and `qy`
/// (`nz_qxy.fna`, no hits at `-evalue 1e-3`) are not refused (NCBI looks nothing up, the
/// manifest's `db.nz.cq.bn.*` rows have its bytes), and `qx` then `q2` (`nz_qx2.fna`) are refused
/// where `q2`'s batch prefetches its subjects (or, with `PRE_FETCH_SEQS_LIMIT=0`, where its report
/// fetches them): with the switch after `qx`'s report, as the stream's other refusals of a batch's
/// reported subjects, without it before any report (outfmt 0's prolog is written before the
/// search).
#[test]
fn a_taxonomy_cache_refusal_comes_with_the_first_batch_that_fetches() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let path = |name: &str| tests.join(name).to_str().unwrap().to_string();
    let db = path("db_search/nz_nc");
    let blastdb = path("seqdb/taxdb/offset_last_past_end");
    for (outfmt, threads, limit) in [("7", 1, None), ("0", 4, None), ("7", 4, Some("0"))] {
        let mut extra = vec![("BATCH_SIZE", "1"), ("BLASTDB", blastdb.as_str())];
        if let Some(limit) = limit {
            extra.push(("PRE_FETCH_SEQS_LIMIT", limit));
        }
        let what = format!("outfmt {outfmt}, PRE_FETCH_SEQS_LIMIT {limit:?}");
        let args = |query: &str| -> Vec<String> {
            [
                "-task",
                "blastn",
                "-dust",
                "no",
                "-evalue",
                "1e-3",
                "-query",
                &path(&format!("db_search/{query}")),
                "-db",
                &db,
                "-outfmt",
                outfmt,
            ]
            .map(String::from)
            .to_vec()
        };
        let run = |query: &str| {
            let args = args(query);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            same_two(&format!("{what}, {query}"), &args, threads, &extra, None)
        };
        let (none, (batches, reason)) = run("nz_qxy.fna");
        assert_eq!(none.status, Some(0), "{what}, no hits: {none:?}");
        assert!(
            batches == 2 && reason == "none",
            "{what}, no hits: {batches} {reason}"
        );
        let late = args("nz_qx2.fna");
        let late: Vec<&str> = late.iter().map(String::as_str).collect();
        let (off, on, (batches, _)) = refused_late(&what, &late, threads, &extra);
        assert!(
            String::from_utf8_lossy(&off.stderr).contains("256-slot taxonomy cache"),
            "{what}: {}",
            String::from_utf8_lossy(&off.stderr)
        );
        assert!(batches >= 1, "{what}");
        // Without the switch every report is written after the loop: none before the refusal.
        let marker = if outfmt == "0" {
            "Query= "
        } else {
            "# Query: "
        };
        assert!(
            !String::from_utf8_lossy(&off.report).contains(marker),
            "{what}: off, no report"
        );
        let (qx, _) = run("nz_qx.fna");
        assert_eq!(qx.status, Some(0), "{what}, qx: {qx:?}");
        reports_before(&what, &off, &on, &qx.report, &qx.stderr, Some((marker, 1)));
    }
}
