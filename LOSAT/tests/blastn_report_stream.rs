//! The reports of each query batch written when the batch loop takes the batch
//! (`LOSAT_X_RPTSTREAM`, `src/algorithm/blastn/blast_engine/x_rpt_stream.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses when LOSAT
//! writes the same reports: each batch's when the loop takes it, or all of them after the loop.
//! So the property under test is that nothing a run writes depends on it: for megablast, blastn
//! and dc-megablast, outfmt 0, 6 and 7, at 1 and 4 threads, with small query batches
//! (`BATCH_SIZE`), the searches without the switch and with it give the same exit status,
//! standard output (or `-out` file) and standard error, and the `[X_STATS] RPTSTREAM` line (the arm
//! with `LOSAT_X_STATS`; its lines are taken out of the comparison) shows that the loop wrote the
//! batches. The same holds with `-out` files, with `-query_loc`, `-subject_loc` and `-lcase_masking`, with `-max_target_seqs`,
//! `-max_hsps` and `-subject_besthit`, with the warnings of the reader and of invalid queries,
//! on the error paths that come after the reports of earlier batches (a reader error, a batch of
//! queries without letters), and with `-db`. Where the switch must not act (LOSAT's refusal of a
//! subject's id, which comes before every report without the switch; a search space that could
//! overflow; the debug output) the `[X_STATS]` line names why and the bytes are the same too. That
//! those bytes equal NCBI's is the business of the blastn fixtures.

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
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
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

/// A database (`-db`, made by `losat makedb` from four genomes): the off bytes, every format,
/// each (task, format) at one thread count. The pass over the database's sequences in the pool
/// (`database::x_streaming_refusal`) runs at 1 and 4 threads in
/// `a_database_search_where_the_switch_must_not_act` (the intact database).
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
        .arg("makedb")
        .arg("-in")
        .arg(dir.path("set.fna"))
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(output.status.success(), "makedb");
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

/// `losat makedb` of LC738868 and of a part of LC738874 (150,000-190,000, with three runs of `N`, so
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
        .arg("makedb")
        .arg("-in")
        .arg(&fasta)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(output.status.success(), "makedb");
    patch(&db);
    db
}

/// The same database search where the switch must act and where it must not: the second
/// sequence's title holding a byte outside ' ' ..= '~' in its header file, or its ambiguity count
/// raised by one in the sequence file. Both are refused when a batch reports the sequence (its
/// title) or traces it back (its ambiguity data), which only the last batches do; without the
/// switch the refusal comes before every report, so the switch reads both before the search and
/// does not act (`off=titles`, `off=ambiguity`): the refusal, no row, the same bytes off / on.
/// The intact database streams, with its ambiguity data decoded by the traceback.
#[test]
fn a_database_search_where_the_switch_must_not_act() {
    let dir = TempDir::new("dbrefuse");
    let query = queries(&dir, "q.fa", b"");
    let run = |db: &str, threads: u32| {
        let args = [
            "-task",
            "megablast",
            "-query",
            query.to_str().unwrap(),
            "-db",
            db,
            "-outfmt",
            "6",
        ];
        same_two(&format!("-db {db}"), &args, threads, &[], None)
    };
    let intact = database(&dir, "intact", |_| {});
    for threads in [1, 4] {
        let (off, (batches, reason)) = run(&intact, threads);
        assert_eq!(off.status, Some(0));
        let rows = String::from_utf8_lossy(&off.report).into_owned();
        assert!(
            rows.contains("\tgnl|BL_ORD_ID|1\t") || rows.contains("LC738874part"),
            "{rows}"
        );
        assert!(batches >= 8 && reason == "none");
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
    for (db, refusal, expected) in [
        (&title, "has a stored title with a byte outside", "titles"),
        (
            &ambiguity,
            "the ambiguity data are shorter than their count",
            "ambiguity",
        ),
    ] {
        for threads in [1, 4] {
            let (off, (batches, reason)) = run(db, threads);
            assert_ne!(off.status, Some(0), "{db}: the refusal");
            assert!(
                String::from_utf8_lossy(&off.stderr).contains(refusal),
                "{db}: {}",
                String::from_utf8_lossy(&off.stderr)
            );
            assert!(
                off.report.is_empty(),
                "{db}: the refusal comes before every row"
            );
            assert_eq!((batches, reason.as_str()), (0, expected), "{db}");
        }
    }
}

/// Where the switch must not act, it does not, and the bytes are the off path's: a subject whose
/// tabular id LOSAT refuses (`Subject_` and a byte in no encoding NCBI recognises), which a later
/// batch reaches (without the switch the refusal comes before every report; written batch by
/// batch, the reports of the first batches would come before it); a search space that could
/// overflow (a `-dbsize` near 2^63); scoring options without Karlin-Altschul values; a negative
/// `CHUNK_SIZE` (whose refusal comes where a batch is split); and the debug output
/// (`LOSAT_DEBUG_COORDS`).
#[test]
fn where_the_switch_must_not_act_the_counters_name_why() {
    let dir = TempDir::new("off");
    let query = queries(&dir, "q.fa", b"");
    // The first subject, LC738868, is hit from the first batch on; the second, named
    // `Subject_\x81`, is a part of LC738874 (the queries' genome) that only the queries of the
    // last batches hit.
    let mut subjects = Vec::new();
    record(&mut subjects, b"LC738868 first", &letters("LC738868.fasta"));
    record(
        &mut subjects,
        b"Subject_\x81 undecodable",
        &letters("LC738874.fasta")[150_000..190_000],
    );
    let subject = dir.path("subjects.fa");
    std::fs::write(&subject, subjects).unwrap();
    let args = [
        "-task",
        "megablast",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    let (off, (batches, reason)) = same_two("Subject_ id", &args, 4, &[], None);
    assert_ne!(off.status, Some(0), "the refusal");
    assert!(
        String::from_utf8_lossy(&off.stderr).contains("'Subject_'"),
        "{}",
        String::from_utf8_lossy(&off.stderr)
    );
    assert!(
        off.report.is_empty(),
        "the refusal comes before every report"
    );
    assert_eq!((batches, reason.as_str()), (0, "subject_ids"));
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
