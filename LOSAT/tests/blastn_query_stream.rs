//! The records of each query batch read from the query file when the batch loop (or a
//! LOSAT_X_QBATCHPAR slot) needs them (`LOSAT_X_QSTREAM`,
//! `src/algorithm/blastn/blast_engine/x_query_stream.rs`), on whole searches.
//!
//! The switch is in the strict set since 2026-10-18 (on unless it is `0`),
//! with a rule (off for `-task megablast`, for a file shorter than twice the
//! largest query batch, and where LOSAT_X_QBATCHPAR's slots may read at least half of the file
//! ahead). It chooses when LOSAT reads the same records: all of them before the
//! search, or a first reading that keeps their lengths and a second one batch by batch, which it does only where the reports may be written
//! batch by batch (LOSAT_X_RPTSTREAM on; off, the `[X_STATS]` line says `off=reports`).
//! So the property under test is that nothing a run writes depends on it: for blastn and
//! dc-megablast, outfmt 0, 6 and 7, at 1 and 4 threads, with small query batches (`BATCH_SIZE`)
//! and LOSAT_X_RPTSTREAM on, the searches without the switch and with it give the same exit
//! status, standard output (or `-out` file) and standard error, and the `[X_STATS] QSTREAM` line
//! (the arm with `LOSAT_X_STATS`; its lines are taken out of the comparison) shows that the
//! batches were read again (the tests keep `BATCH_SIZE` small, so the batch part of the rule does
//! not apply). Where the switch does not act (its rule, or LOSAT_X_RPTSTREAM off), `run.rs` reads
//! the records before the search as without the switch (`read_queries`, `Store::off`), so one
//! case each shows why and the same bytes: LOSAT_X_RPTSTREAM off (`off=reports`), `-task
//! megablast` (`off=megablast`) and LOSAT_X_QBATCHPAR slots (`off=slots`). The same holds with
//! `-out` files, with `-query_loc` (skipped records and queries without data), `-subject_loc` and
//! `-lcase_masking`, with `-max_target_seqs`, `-max_hsps` and `-subject_besthit`, with the
//! reader's warnings and errors (late in the file, and in the first record, before any output), a
//! batch of queries without letters, LOSAT's refusals that come before every report
//! (LOSAT_X_RPTSTREAM then keeps its reports after the loop, and for a refusal that comes when
//! they are written the records of every batch are kept for them), and `-db`. Those cases run at
//! one thread count each, 1 or 4: the switch reads records on other threads only for the slots,
//! and its rule reads the thread count only for the slots' window. Standard input and a pipe are
//! read once (the `[X_STATS]` line names why), and the switch is on by default and off by its
//! rule. That those bytes equal NCBI's is the business of the blastn fixtures.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), and this
/// file's own on 2026-10-18 (perf/bn-memory-rules), each `0`, the default before, unless the
/// search sets it: every search here runs as it was written, with the switches off that it does
/// not name (the default runs in `default_bundle.rs` and `the_switch_is_on_by_default`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_QSTREAM", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_QSTREAM", "1");
const RPTSTREAM: (&str, &str) = ("LOSAT_X_RPTSTREAM", "1");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
/// Two slots whatever QBATCHPAR's own choice of the layout (its tuning variable), so that slots run.
const TWO_SLOTS: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
/// The tasks the switch acts for (its rule keeps it off for megablast).
const ACTING: [&str; 2] = ["blastn", "dc-megablast"];

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
            "losat_qstream_{}_{}_{tag}",
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
fn query_text(tail: &[u8]) -> Vec<u8> {
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
    text
}

fn queries(dir: &TempDir, name: &str, tail: &[u8]) -> PathBuf {
    let path = dir.path(name);
    std::fs::write(&path, query_text(tail)).expect("write the queries");
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
/// 4,000 (`CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` removed, unless in `switches`), `stdin` as its
/// standard input (none without).
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing (it is set to a file of this test when the environment has
/// none). With `out`, the report is that file's bytes (standard output must then be empty).
fn blastn_with_input(
    args: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
    out: Option<&Path>,
    stdin: Option<&[u8]>,
) -> Run {
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
            std::env::temp_dir().join(format!("losat_qstream_{}_worker.log", std::process::id())),
        );
    }
    command
        .arg("blastn")
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "4000")
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(switches.iter().copied())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if let Some(out) = out {
        let _ = std::fs::remove_file(out);
        command.arg("-out").arg(out);
    }
    command.stdin(if stdin.is_some() {
        Stdio::piped()
    } else {
        Stdio::null()
    });
    let mut child = command.spawn().expect("run LOSAT");
    if let Some(bytes) = stdin {
        let mut pipe = child.stdin.take().expect("standard input");
        let bytes = bytes.to_vec();
        // Written from another thread: the search reads it while its output fills the pipes.
        std::thread::spawn(move || {
            let _ = pipe.write_all(&bytes);
        });
    }
    let output = child.wait_with_output().expect("wait for LOSAT");
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

fn blastn(args: &[&str], threads: u32, switches: &[(&str, &str)], out: Option<&Path>) -> Run {
    blastn_with_input(args, threads, switches, out, None)
}

/// The numbers of the `[X_STATS] QSTREAM records=R read=N batches=B slots=S ... off=O` line, and
/// `precomputed=` of the `[X_STATS] QBATCHPAR` line where there is one (LOSAT_X_QBATCHPAR on).
#[derive(Debug)]
struct Numbers {
    batches: u64,
    read: u64,
    slots: u64,
    off: String,
    precomputed: Option<u64>,
}

/// Takes the `[X_STATS]` lines out of `run`'s standard error and returns the numbers of the
/// `[X_STATS] QSTREAM` line (and of the `[X_STATS] QBATCHPAR` line).
fn take_stats_lines(run: &mut Run) -> Numbers {
    let text = String::from_utf8_lossy(&run.stderr).into_owned();
    let lines: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| line.starts_with("[X_STATS] QSTREAM "))
        .collect();
    assert_eq!(lines.len(), 1, "one QSTREAM line:\n{text}");
    let line = lines[0].trim_end();
    let field = |name: &str| -> String {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
            .to_string()
    };
    let numbers = Numbers {
        batches: field("batches=").parse().unwrap(),
        read: field("read=").parse().unwrap(),
        slots: field("slots=").parse().unwrap(),
        off: field("off="),
        precomputed: text
            .lines()
            .find(|line| line.starts_with("[X_STATS] QBATCHPAR "))
            .map(|line| {
                line.split_whitespace()
                    .find_map(|word| word.strip_prefix("precomputed="))
                    .and_then(|value| value.parse().ok())
                    .unwrap_or_else(|| panic!("no precomputed= in: {line}"))
            }),
    };
    run.stderr = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_STATS] "))
        .collect::<String>()
        .into_bytes();
    numbers
}

/// The run without the switch and with it (with `extra` switches in both), and with it and
/// `LOSAT_X_STATS`: the same status, report and standard error (the counters' lines taken out of
/// the last). Returns the off run and the counters of the `[X_STATS] QSTREAM` line.
fn same_two_with_input(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
    stdin: Option<&[u8]>,
) -> (Run, Numbers) {
    let off = blastn_with_input(args, threads, extra, out, stdin);
    let mut with = extra.to_vec();
    with.push(SWITCH);
    let on = blastn_with_input(args, threads, &with, out, stdin);
    assert!(
        on == off,
        "{what} threads {threads} {extra:?}: on differs\n{on:?}\n{off:?}"
    );
    with.push(STATS);
    let mut counted = blastn_with_input(args, threads, &with, out, stdin);
    let numbers = take_stats_lines(&mut counted);
    assert!(
        counted == off,
        "{what} threads {threads} {extra:?}: on with statistics differs\n{counted:?}\n{off:?}"
    );
    (off, numbers)
}

fn same_two(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
) -> (Run, Numbers) {
    same_two_with_input(what, args, threads, extra, out, None)
}

/// The switch read every batch again: about a dozen batches, each record read (the 19 input
/// records at least once).
fn read_batch_by_batch(what: &str, numbers: &Numbers) {
    assert!(
        numbers.off == "none" && numbers.batches >= 8 && numbers.read >= 19,
        "{what}: {numbers:?}"
    );
}

/// With LOSAT_X_RPTSTREAM (`reports`) the switch read every batch again, except for megablast (its
/// rule: `off=megablast`); without it the switch did not act (`off=reports`: the records are read
/// before the search).
fn read_as_the_reports_need(what: &str, task: &str, reports: &[(&str, &str)], numbers: &Numbers) {
    if reports.is_empty() || task == "megablast" {
        let why = if reports.is_empty() {
            "reports"
        } else {
            "megablast"
        };
        assert!(
            numbers.off == why && numbers.batches == 0 && numbers.read == 0,
            "{what}: {numbers:?}"
        );
    } else {
        read_batch_by_batch(what, numbers);
    }
}

/// blastn and dc-megablast, outfmt 0, 6 and 7, 1 and 4 threads, the reports per batch
/// (LOSAT_X_RPTSTREAM): the off bytes, with every batch read again, and the warnings of the
/// reader and of the invalid queries on standard error. Where the switch does not act (its rule
/// for megablast; LOSAT_X_RPTSTREAM off) the records are read as without it: one case each.
#[test]
fn every_task_format_thread_count_and_report_mode_gives_the_off_bytes() {
    let dir = TempDir::new("formats");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let mut cases = Vec::new();
    for task in ACTING {
        for outfmt in ["0", "6", "7"] {
            for threads in [1, 4] {
                cases.push((task, outfmt, threads, &[RPTSTREAM][..]));
            }
        }
    }
    cases.push(("megablast", "6", 4, &[RPTSTREAM][..]));
    cases.push(("blastn", "0", 4, &[][..]));
    for (task, outfmt, threads, reports) in cases {
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
        let what = format!("{task} outfmt {outfmt} {reports:?}");
        let (off, numbers) = same_two(&what, &args, threads, reports, None);
        assert_eq!(off.status, Some(0), "{what}");
        assert!(!off.report.is_empty(), "{what}: no report");
        let stderr = String::from_utf8_lossy(&off.stderr);
        assert!(
            stderr.contains("FASTA-Reader: Ignoring invalid residues")
                && stderr.contains("Sequence contains no data"),
            "{what}: {stderr}"
        );
        read_as_the_reports_need(&what, task, reports, &numbers);
        assert_eq!(numbers.slots, 0, "{what}");
    }
}

/// `-out` files (each format once, with a task the switch acts for) and LOSAT_X_QBATCHPAR's slots
/// (two slots of two threads, `LOSAT_X_QBATCHPARSLOTS=2`): the off bytes; the slots did search
/// batches ahead of the loop (`precomputed=`), with the switch kept off by its rule on this small
/// file (`off=slots`), so the records are read as without it, whatever the task (one case).
#[test]
fn out_files_and_qbatchpar_slots_give_the_off_bytes() {
    let dir = TempDir::new("slots");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let out = dir.path("report.txt");
    let args = |task, outfmt| {
        [
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            outfmt,
        ]
    };
    for (task, outfmt) in [("blastn", "0"), ("dc-megablast", "6"), ("blastn", "7")] {
        let what = format!("{task} outfmt {outfmt} -out");
        let (off, numbers) = same_two(&what, &args(task, outfmt), 4, &[RPTSTREAM], Some(&out));
        assert!(!off.report.is_empty(), "{what}");
        read_batch_by_batch(&what, &numbers);
    }
    let what = "blastn outfmt 7 QBATCHPAR";
    let extra = [QBATCHPAR, TWO_SLOTS, RPTSTREAM];
    let (off, numbers) = same_two(what, &args("blastn", "7"), 4, &extra, None);
    assert!(!off.report.is_empty(), "{what}");
    assert!(
        numbers.off == "slots" && numbers.slots == 0 && numbers.read == 0,
        "{what}: {numbers:?}"
    );
    // The slots searched batches ahead of the loop (the run with the counters).
    assert!(
        numbers
            .precomputed
            .is_some_and(|precomputed| precomputed > 0),
        "{what}: {numbers:?}"
    );
}

/// `-query_loc` (a range that skips the records without letters and the two-letter one, and one
/// that makes the two-letter record a query without data), `-subject_loc` with `-lcase_masking`
/// (queries with lower case runs), and `-max_target_seqs`, `-max_hsps` and `-subject_besthit`:
/// the off bytes, with the reports per batch, each case at one thread count.
#[test]
fn ranges_masks_and_hit_list_options_give_the_off_bytes() {
    let dir = TempDir::new("ranges");
    let query = queries(&dir, "q.fa", b"");
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
    let cases: [(&str, u32, &[&str]); 5] = [
        ("blastn", 1, &["-outfmt", "0", "-query_loc", "100-1500"]),
        ("blastn", 4, &["-outfmt", "6", "-query_loc", "3-1000"]),
        (
            "dc-megablast",
            4,
            &["-outfmt", "7", "-query_loc", "100-1500", "-lcase_masking"],
        ),
        (
            "dc-megablast",
            1,
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
            4,
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
        let (off, numbers) = same_two(&what, &args, threads, &[RPTSTREAM], None);
        assert_eq!(off.status, Some(0), "{what}");
        assert!(!off.report.is_empty(), "{what}");
        // The reader's warning of the record with digits (cut with `-query_loc`: the warning is
        // the whole record's, written with its batch).
        assert!(
            String::from_utf8_lossy(&off.stderr)
                .contains("FASTA-Reader: Ignoring invalid residues"),
            "{what}: {}",
            String::from_utf8_lossy(&off.stderr)
        );
        read_batch_by_batch(&what, &numbers);
    }
}

/// The errors of a later batch that come after the reports of the batches before (a reader error
/// late in the file; a batch of queries without letters; with `-query_loc`, a batch of skipped
/// records only), a reader error in the first record (before any output), and LOSAT's refusals,
/// which come before every report (a subject id that a later batch reaches, scoring without
/// Karlin-Altschul values, a negative `CHUNK_SIZE`): the off bytes, status and standard error,
/// with the reports per batch and tasks the switch acts for (each query with each format once).
#[test]
fn errors_and_refusals_give_the_off_bytes() {
    let dir = TempDir::new("errors");
    let subject = fasta("LC738868.fasta");
    let reader_error = queries(&dir, "reader.fa", b">bad\n@@@@\n");
    let mut tail = Vec::new();
    record(
        &mut tail,
        b"long",
        &letters("LC738874.fasta")[200_000..204_500],
    );
    tail.extend_from_slice(b">e1\n>e2\n>e3\n");
    let empty_batch = queries(&dir, "empty.fa", &tail);
    let mut first = b">bad\n@@@@\n".to_vec();
    first.extend_from_slice(&query_text(b""));
    let first_error = dir.path("first.fa");
    std::fs::write(&first_error, first).unwrap();
    for (query, status, error, before) in [
        (
            &reader_error,
            1,
            "BLAST query error: CFastaReader: Near line",
            true,
        ),
        (
            &empty_batch,
            3,
            "BLAST engine error: Warning: Sequence contains no data",
            true,
        ),
        (
            &first_error,
            1,
            "BLAST query error: CFastaReader: Near line",
            false,
        ),
    ] {
        for (task, outfmt, threads) in [
            ("blastn", "0", 1),
            ("blastn", "7", 4),
            ("dc-megablast", "6", 4),
        ] {
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
            let (off, numbers) = same_two(&what, &args, threads, &[RPTSTREAM], None);
            assert_eq!(off.status, Some(status), "{what}");
            assert!(
                String::from_utf8_lossy(&off.stderr).contains(error),
                "{what}: {}",
                String::from_utf8_lossy(&off.stderr)
            );
            if before {
                assert!(!off.report.is_empty(), "{what}");
                read_batch_by_batch(&what, &numbers);
            } else {
                // Only the outfmt 0 prolog comes before the error (written before the first batch
                // is read).
                assert!(
                    outfmt == "0" || off.report.is_empty(),
                    "{what}: {}",
                    String::from_utf8_lossy(&off.report)
                );
                assert_eq!((numbers.off.as_str(), numbers.read), ("none", 0));
            }
        }
    }
    // A batch of skipped records only (`-query_loc` past the end of the last records).
    let mut tail = Vec::new();
    record(
        &mut tail,
        b"long",
        &letters("LC738874.fasta")[200_000..204_500],
    );
    record(&mut tail, b"short1", b"ACGTACGTAC");
    record(&mut tail, b"short2", b"ACGTACGTAC");
    let skipped = queries(&dir, "skipped.fa", &tail);
    let args = [
        "-task",
        "blastn",
        "-query",
        skipped.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
        "-query_loc",
        "100-1500",
    ];
    let (off, numbers) = same_two("skipped", &args, 4, &[RPTSTREAM], None);
    assert_eq!(
        off.status,
        Some(3),
        "{:?}",
        String::from_utf8_lossy(&off.stderr)
    );
    assert!(String::from_utf8_lossy(&off.stderr).contains("Empty CBlastQueryVector"));
    read_batch_by_batch("skipped", &numbers);
    // LOSAT's refusals: a subject id that only the last batches reach, scoring without
    // Karlin-Altschul values, a negative CHUNK_SIZE.
    let query = queries(&dir, "q.fa", b"");
    let mut subjects = Vec::new();
    record(&mut subjects, b"LC738868 first", &letters("LC738868.fasta"));
    record(
        &mut subjects,
        b"Subject_\x81 undecodable",
        &letters("LC738874.fasta")[150_000..190_000],
    );
    let refused = dir.path("subjects.fa");
    std::fs::write(&refused, subjects).unwrap();
    let plain = |subject: &Path| -> Vec<String> {
        [
            "-task",
            "blastn",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "6",
        ]
        .map(String::from)
        .to_vec()
    };
    // dc-megablast: only the queries of the last batches hit the second subject (blastn hits it
    // from the first batch on).
    let mut ids = plain(&refused);
    ids[1] = "dc-megablast".to_string();
    let mut scoring = plain(&subject);
    scoring.extend(
        [
            "-reward",
            "2",
            "-penalty",
            "-3",
            "-gapopen",
            "1",
            "-gapextend",
            "1",
        ]
        .map(String::from),
    );
    let chunks = [("CHUNK_SIZE", "-5"), ("OVERLAP_CHUNK_SIZE", "-10")];
    for (what, args, env) in [
        ("Subject_ id", &ids, &[][..]),
        ("scoring", &scoring, &[][..]),
        ("CHUNK_SIZE", &plain(&subject), &chunks[..]),
    ] {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let mut extra = env.to_vec();
        extra.push(RPTSTREAM);
        let (off, numbers) = same_two(what, &args, 4, &extra, None);
        assert_ne!(off.status, Some(0), "{what}: the refusal");
        assert!(
            off.report.is_empty(),
            "{what}: no report before the refusal"
        );
        // The switch acts and LOSAT_X_RPTSTREAM keeps the reports after the loop. The Subject_ id
        // is refused when they are written, so the records of every batch are kept for them; the
        // other two are refused in the first batch's search.
        assert_eq!(numbers.off, "none", "{what}");
    }
}

/// Standard input (`-query -`) and a pipe (`-query /dev/stdin` with standard input a pipe) are
/// read once, before the search, with the switch too: the off bytes, and the counters' line names
/// why (with the reports per batch: without them the line says `off=reports`, as for a file).
#[test]
fn standard_input_and_pipes_are_read_once() {
    let subject = fasta("LC738868.fasta");
    let text = query_text(b"");
    let mut paths = vec![("-", "standard_input")];
    if cfg!(unix) {
        paths.push(("/dev/stdin", "not_a_file"));
    }
    for (path, reason) in paths {
        let args = [
            "-task",
            "megablast",
            "-query",
            path,
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "6",
        ];
        let (off, numbers) = same_two_with_input(path, &args, 4, &[RPTSTREAM], None, Some(&text));
        assert_eq!(off.status, Some(0), "{path}");
        assert!(!off.report.is_empty(), "{path}");
        assert_eq!(
            (numbers.off.as_str(), numbers.batches, numbers.read),
            (reason, 0, 0),
            "{path}"
        );
    }
}

/// A database (`-db`, made by `losat makedb` from four genomes): the off bytes, every format,
/// each at one thread count, the reports per batch.
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
        ("dc-megablast", "0", 4),
        ("blastn", "6", 1),
        ("dc-megablast", "7", 4),
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
        let (off, numbers) = same_two(&what, &args, threads, &[RPTSTREAM], None);
        assert_eq!(off.status, Some(0), "{what}");
        assert!(!off.report.is_empty(), "{what}");
        read_batch_by_batch(&what, &numbers);
    }
}

/// The switch is on by default (an empty value): with LOSAT_X_RPTSTREAM it reads the batches
/// again (its `[X_STATS]` line with `off=none`); with `LOSAT_X_QSTREAM=0` the records are read
/// before the search (no line). Its rule, where the line says why: `-task megablast`
/// (`off=megablast`); a file shorter than twice the largest batch (`off=batch`; here `BATCH_SIZE`
/// half the file's bytes plus one, and half the file's bytes on); where the slots of
/// LOSAT_X_QBATCHPAR may run, a file shorter than twice what their budget holds (`off=slots`; a
/// budget set by `LOSAT_X_QBATCHPARMEM` for a small file at 2 threads, none at one thread). The
/// bytes are the same.
#[test]
fn the_switch_is_on_by_default() {
    let dir = TempDir::new("default");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let args = |task| {
        [
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "6",
        ]
    };
    let blastn_args = args("blastn");
    let line = |run: &Run| -> String {
        String::from_utf8_lossy(&run.stderr)
            .lines()
            .find(|line| line.starts_with("[X_STATS] QSTREAM "))
            .map(str::to_string)
            .unwrap_or_default()
    };
    let default = ("LOSAT_X_QSTREAM", "");
    let on = blastn(&blastn_args, 4, &[RPTSTREAM, STATS, default], None);
    let read = line(&on);
    assert!(
        read.contains(" off=none") && read.contains(" read=19 "),
        "{read}"
    );
    let switches = &[RPTSTREAM, STATS][..];
    let off = blastn(&blastn_args, 4, switches, None);
    assert_eq!(line(&off), "", "{switches:?}");
    assert_eq!((off.status, &off.report), (on.status, &on.report));
    // The rule: (task, threads, switches, the line's `off=`); each against the switch `0`.
    let half = std::fs::metadata(&query).unwrap().len() / 2;
    let (batch_off, batch_on) = ((half + 1).to_string(), half.to_string());
    let slots = [
        RPTSTREAM,
        STATS,
        default,
        QBATCHPAR,
        ("LOSAT_X_QBATCHPARMEM", "64"),
    ];
    let cases = [
        ("megablast", 4, vec![RPTSTREAM, STATS, default], "megablast"),
        (
            "blastn",
            4,
            vec![RPTSTREAM, STATS, default, ("BATCH_SIZE", &batch_off)],
            "batch",
        ),
        (
            "blastn",
            4,
            vec![RPTSTREAM, STATS, default, ("BATCH_SIZE", &batch_on)],
            "none",
        ),
        ("blastn", 2, slots.to_vec(), "slots"),
        ("blastn", 1, slots.to_vec(), "none"),
    ];
    for (task, threads, switches, why) in cases {
        let what = format!("{task} {threads} {switches:?}");
        let ruled = blastn(&args(task), threads, &switches, None);
        let read = line(&ruled);
        let off_word = format!("off={why}");
        assert!(
            read.split_whitespace().any(|word| word == off_word),
            "{what}: {read}"
        );
        if why != "none" {
            assert!(read.contains(" read=0 "), "{what}: {read}");
        }
        let mut zero = switches.clone();
        zero.push(("LOSAT_X_QSTREAM", "0"));
        let off = blastn(&args(task), threads, &zero, None);
        assert_eq!(
            (ruled.status, &ruled.report),
            (off.status, &off.report),
            "{what}"
        );
    }
}

/// The records of the batches are kept only where a writer reads them after their batch: none
/// with LOSAT_X_RPTSTREAM acting (each batch's records are dropped after its reports).
#[test]
fn the_records_are_kept_only_for_a_writer_after_the_loop() {
    let dir = TempDir::new("kept");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let args = [
        "-task",
        "blastn",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    let kept = |switches: &[(&str, &str)]| -> String {
        let run = blastn(&args, 4, switches, None);
        let stderr = String::from_utf8_lossy(&run.stderr).into_owned();
        stderr
            .lines()
            .find(|line| line.starts_with("[X_STATS] QSTREAM "))
            .and_then(|line| {
                line.split_whitespace()
                    .find_map(|word| word.strip_prefix("kept="))
            })
            .unwrap_or_else(|| panic!("{switches:?}: {stderr}"))
            .to_string()
    };
    assert_eq!(kept(&[SWITCH, RPTSTREAM, STATS]), "0");
    assert_eq!(kept(&[SWITCH, RPTSTREAM, STATS, QBATCHPAR, TWO_SLOTS]), "0");

    // LOSAT_X_RPTSTREAM on but refusing the search (a `-dbsize` whose search space could
    // overflow: `off=search_space`) while the switch acts: the reports are written after the
    // loop from the records it kept, and the bytes are those without the switch.
    let mut huge = args.to_vec();
    huge.extend_from_slice(&["-dbsize", "1200000000000000"]);
    let stats_line = |stderr: &str, name: &str| -> String {
        stderr
            .lines()
            .find(|line| line.starts_with(name))
            .unwrap_or_else(|| panic!("no {name} line: {stderr}"))
            .to_string()
    };
    let refused = blastn(&huge, 4, &[SWITCH, RPTSTREAM, STATS], None);
    let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
    let rpt = stats_line(&stderr, "[X_STATS] RPTSTREAM ");
    assert!(
        rpt.split_whitespace()
            .any(|word| word == "off=search_space"),
        "{rpt}"
    );
    let stream = stats_line(&stderr, "[X_STATS] QSTREAM ");
    assert!(
        stream.split_whitespace().any(|word| word == "off=none"),
        "{stream}"
    );
    let kept_records: u64 = stream
        .split_whitespace()
        .find_map(|word| word.strip_prefix("kept="))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("{stream}"));
    assert!(kept_records > 0, "{stream}");
    let zero = blastn(&huge, 4, &[RPTSTREAM, ("LOSAT_X_QSTREAM", "0")], None);
    assert_eq!(refused.status, Some(0), "{stderr}");
    assert!(!refused.report.is_empty(), "{stderr}");
    let without_stats = |stderr: &str| -> String {
        stderr
            .split_inclusive('\n')
            .filter(|line| !line.starts_with("[X_STATS] "))
            .collect()
    };
    assert_eq!(
        (refused.status, &refused.report, without_stats(&stderr)),
        (
            zero.status,
            &zero.report,
            without_stats(&String::from_utf8_lossy(&zero.stderr))
        )
    );
}
