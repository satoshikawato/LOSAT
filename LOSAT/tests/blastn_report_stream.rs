//! The reports of each query batch written when the batch loop takes the batch
//! (`LOSAT_X_RPTSTREAM`, `src/algorithm/blastn/blast_engine/x_rpt_stream.rs`) on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses when LOSAT
//! writes the same reports: each batch's when the loop takes it, or all of them after the loop.
//! So the property under test is that nothing a run writes depends on it: for megablast, blastn
//! and dc-megablast, outfmt 0, 6 and 7, at 1 and 4 threads, with small query batches
//! (`BATCH_SIZE`), the searches without the switch, with it and with its shadow
//! (`LOSAT_X_RPTSTREAMSHADOW`, which also writes the reports the off path's way in memory and
//! compares the bytes) give the same exit status, standard output (or `-out` file) and standard
//! error (the shadow's line taken out), and the shadow line shows that the loop wrote the
//! batches. The same holds with LOSAT_X_QBATCHPAR slots searching the batches ahead while the
//! loop writes, with `-query_loc`, `-subject_loc` and `-lcase_masking`, with `-max_target_seqs`,
//! `-max_hsps` and `-subject_besthit`, with the warnings of the reader and of invalid queries,
//! on the error paths that come after the reports of earlier batches (a reader error, a batch of
//! queries without letters), and with `-db`. Where the switch must not act (LOSAT's refusal of a
//! subject's id, which comes before every report without the switch; a search space that could
//! overflow; the debug output) the shadow line names why and the bytes are the same too. That
//! those bytes equal NCBI's is the business of the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: (&str, &str) = ("LOSAT_X_RPTSTREAM", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_RPTSTREAMSHADOW", "1");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// The shadow's line `[X_SHADOW] RPTSTREAM batches=B bytes=N checks=C off=O identical` taken out
/// of `run`'s standard error: (B, N, C, O).
fn take_shadow_line(run: &mut Run) -> (u64, u64, u64, String) {
    let text = String::from_utf8_lossy(&run.stderr).into_owned();
    let lines: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| line.starts_with("[X_SHADOW] RPTSTREAM "))
        .collect();
    assert_eq!(lines.len(), 1, "one shadow line:\n{text}");
    let line = lines[0].trim_end();
    assert!(line.ends_with(" identical"), "{line}");
    let field = |name: &str| -> String {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
            .to_string()
    };
    let numbers = (
        field("batches=").parse().unwrap(),
        field("bytes=").parse().unwrap(),
        field("checks=").parse().unwrap(),
        field("off="),
    );
    run.stderr = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_SHADOW] RPTSTREAM "))
        .collect::<String>()
        .into_bytes();
    numbers
}

/// The run without the switch, with it and with the shadow (with `extra` switches in all three):
/// the same status, report and standard error. Returns the off run and the shadow's numbers.
fn same_three(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
) -> (Run, (u64, u64, u64, String)) {
    let off = blastn(args, threads, extra, out);
    let mut with = extra.to_vec();
    with.push(SWITCH);
    let on = blastn(args, threads, &with, out);
    assert!(
        on == off,
        "{what} threads {threads}: on differs\n{on:?}\n{off:?}"
    );
    let mut with = extra.to_vec();
    with.push(SHADOW);
    let mut shadow = blastn(args, threads, &with, out);
    let numbers = take_shadow_line(&mut shadow);
    assert!(
        shadow == off,
        "{what} threads {threads}: shadow differs\n{shadow:?}\n{off:?}"
    );
    (off, numbers)
}

/// Every task, outfmt 0, 6 and 7, 1 and 4 threads: the off bytes, with every batch written by
/// the loop (the shadow compared once, at the end) and the warnings of the reader and of the
/// invalid queries on standard error.
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
                let (off, (batches, bytes, checks, reason)) =
                    same_three(&what, &args, threads, &[], None);
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
                assert_eq!(checks, 1, "{what}");
                // Every byte after the outfmt 0 prolog (written before the search) was compared.
                let written = off.report.len() + off.stderr.len();
                if outfmt == "0" {
                    assert!(bytes > 0 && (bytes as usize) < written, "{what}");
                } else {
                    assert_eq!(bytes as usize, written, "{what}");
                }
            }
        }
    }
}

/// `-out` files, and LOSAT_X_QBATCHPAR's slots (which search the batches ahead while the loop
/// writes): the off bytes; the slots did search batches.
#[test]
fn out_files_and_qbatchpar_slots_give_the_off_bytes() {
    let dir = TempDir::new("slots");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let out = dir.path("report.txt");
    for task in TASKS {
        for outfmt in ["0", "6", "7"] {
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
            let what = format!("{task} outfmt {outfmt} -out");
            let (off, (batches, _, _, reason)) = same_three(&what, &args, 4, &[], Some(&out));
            assert!(
                !off.report.is_empty() && reason == "none" && batches >= 8,
                "{what}"
            );
            let what = format!("{task} outfmt {outfmt} QBATCHPAR");
            let (off, (batches, _, _, reason)) = same_three(&what, &args, 4, &[QBATCHPAR], None);
            assert!(
                !off.report.is_empty() && reason == "none" && batches >= 8,
                "{what}"
            );
            // The slots searched batches ahead of the loop.
            let stats = blastn(&args, 4, &[QBATCHPAR, SWITCH, ("LOSAT_X_STATS", "1")], None);
            let stderr = String::from_utf8_lossy(&stats.stderr);
            let precomputed: u64 = stderr
                .lines()
                .find(|line| line.starts_with("[X_STATS] QBATCHPAR "))
                .and_then(|line| {
                    line.split_whitespace()
                        .find_map(|word| word.strip_prefix("precomputed="))
                })
                .and_then(|value| value.parse().ok())
                .unwrap_or_else(|| panic!("{what}: {stderr}"));
            assert!(precomputed > 0, "{what}: {stderr}");
            assert!(
                stderr.contains("[X_STATS] RPTSTREAM batches=") && stderr.contains(" off=none"),
                "{what}: {stderr}"
            );
        }
    }
}

/// `-query_loc`, `-subject_loc` with `-lcase_masking`, and `-max_target_seqs`, `-max_hsps` and
/// `-subject_besthit`: the off bytes.
#[test]
fn ranges_and_hit_list_options_give_the_off_bytes() {
    let dir = TempDir::new("ranges");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    let cases: [(&str, &[&str]); 4] = [
        ("blastn", &["-outfmt", "0", "-query_loc", "100-1500"]),
        ("megablast", &["-outfmt", "6", "-query_loc", "100-1500"]),
        (
            "dc-megablast",
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
    for (task, extra) in cases {
        let mut args = vec![
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        for threads in [1, 4] {
            let what = format!("{task} {extra:?}");
            let (off, (batches, _, _, reason)) = same_three(&what, &args, threads, &[], None);
            assert_eq!(off.status, Some(0), "{what}");
            assert!(
                !off.report.is_empty() && reason == "none" && batches >= 8,
                "{what}"
            );
        }
    }
}

/// The errors of a later batch that come after the reports of the batches before (a reader
/// error; a batch of queries without letters): the off bytes, status and standard error, the
/// shadow compared where the error is met.
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
    for (query, status, error) in [
        (
            &reader_error,
            1,
            "BLAST query error: CFastaReader: Near line",
        ),
        (
            &empty_batch,
            3,
            "BLAST engine error: Warning: Sequence contains no data",
        ),
    ] {
        for (task, outfmt) in [("megablast", "7"), ("blastn", "0"), ("dc-megablast", "6")] {
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
            for threads in [1, 4] {
                let what = format!("{task} outfmt {outfmt} {query:?}");
                let (off, (batches, _, checks, reason)) =
                    same_three(&what, &args, threads, &[], None);
                assert_eq!(off.status, Some(status), "{what}");
                assert!(
                    String::from_utf8_lossy(&off.stderr).contains(error),
                    "{what}: {}",
                    String::from_utf8_lossy(&off.stderr)
                );
                assert!(!off.report.is_empty(), "{what}");
                assert!(reason == "none" && batches >= 8 && checks == 1, "{what}");
            }
        }
    }
}

/// A database (`-db`, made by `LOSAT makedb` from four genomes): the off bytes, every format.
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
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("makedb")
        .arg("-in")
        .arg(dir.path("set.fna"))
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run LOSAT makedb");
    assert!(output.status.success(), "makedb");
    for (task, outfmt) in [
        ("megablast", "0"),
        ("blastn", "6"),
        ("dc-megablast", "7"),
        ("megablast", "7"),
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
        for threads in [1, 4] {
            let what = format!("-db {task} outfmt {outfmt}");
            let (off, (batches, _, _, reason)) = same_three(&what, &args, threads, &[], None);
            assert_eq!(off.status, Some(0), "{what}");
            assert!(
                !off.report.is_empty() && reason == "none" && batches >= 8,
                "{what}"
            );
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
fn where_the_switch_must_not_act_the_shadow_names_why() {
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
    let (off, (batches, _, checks, reason)) = same_three("Subject_ id", &args, 4, &[], None);
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
    assert_eq!((batches, checks, reason.as_str()), (0, 0, "subject_ids"));
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
    let (_, (batches, _, _, reason)) = same_three("-dbsize", &huge, 1, &[], None);
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
    let (off, (batches, _, _, reason)) = same_three("scoring", &scoring, 4, &[], None);
    assert_eq!(off.status, Some(3), "NCBI's Karlin-Altschul error");
    assert_eq!((batches, reason.as_str()), (0, "scoring"));
    let chunks = [("CHUNK_SIZE", "-5"), ("OVERLAP_CHUNK_SIZE", "-10")];
    let (off, (batches, _, _, reason)) = same_three("CHUNK_SIZE", &plain, 4, &chunks, None);
    assert_ne!(off.status, Some(0), "LOSAT's refusal");
    assert_eq!((batches, reason.as_str()), (0, "chunk_size"));
    let (off, (batches, _, _, reason)) =
        same_three("debug", &plain, 1, &[("LOSAT_DEBUG_COORDS", "1")], None);
    assert!(!off.report.is_empty());
    assert_eq!((batches, reason.as_str()), (0, "debug"));
}
