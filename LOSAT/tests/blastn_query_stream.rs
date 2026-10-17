//! The records of each query batch read from the query file when the batch loop (or a
//! LOSAT_X_QBATCHPAR slot) needs them (`LOSAT_X_QSTREAM`,
//! `src/algorithm/blastn/blast_engine/x_query_stream.rs`), on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses when LOSAT
//! reads the same records: all of them before the search, or a first reading that keeps their
//! lengths and a second one batch by batch. So the property under test is that nothing a run
//! writes depends on it: for megablast, blastn and dc-megablast, outfmt 0, 6 and 7, at 1 and 4
//! threads, with small query batches (`BATCH_SIZE`), with LOSAT_X_RPTSTREAM off (the records of
//! every batch kept for the reports after the loop) and on (dropped after each batch's reports),
//! the searches without the switch, with it and with its shadow (`LOSAT_X_QSTREAMSHADOW`, which
//! also reads the whole file the old way and compares the records of every batch and every slot's
//! copy) give the same exit status, standard output (or `-out` file) and standard error (the
//! shadow's line taken out), and the shadow line shows that the batches were read again. The same
//! holds with LOSAT_X_QBATCHPAR slots copying the records of the batches they search ahead, with
//! `-query_loc` (skipped records and queries without data), `-subject_loc` and `-lcase_masking`,
//! with `-max_target_seqs`, `-max_hsps` and `-subject_besthit`, with the reader's warnings and
//! errors (late in the file, and in the first record, before any output), a batch of queries
//! without letters, LOSAT's refusals that come before every report, and `-db`. Standard input
//! and a pipe are read once (the shadow line names why). That those bytes equal NCBI's is the
//! business of the blastn fixtures.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: (&str, &str) = ("LOSAT_X_QSTREAM", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_QSTREAMSHADOW", "1");
const RPTSTREAM: (&str, &str) = ("LOSAT_X_RPTSTREAM", "1");
const RPTSTREAMSHADOW: (&str, &str) = ("LOSAT_X_RPTSTREAMSHADOW", "1");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
/// Two slots whatever QBATCHPAR's own choice of the layout (its tuning variable), so that slots run.
const TWO_SLOTS: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const TASKS: [&str; 3] = ["megablast", "blastn", "dc-megablast"];
/// The reports after the loop, and each batch's when the loop takes it.
const REPORTS: [&[(&str, &str)]; 2] = [&[], &[RPTSTREAM]];

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
            std::env::temp_dir().join(format!("losat_qstream_{}_worker.log", std::process::id())),
        );
    }
    command
        .arg("blastn")
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "4000")
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

/// The shadow's numbers: `[X_SHADOW] QSTREAM batches=B records=R slot_copies=S off=O identical`.
#[derive(Debug)]
struct Shadow {
    batches: u64,
    records: u64,
    slot_copies: u64,
    off: String,
}

/// Takes the shadow's line out of `run`'s standard error (also the line of LOSAT_X_RPTSTREAM's
/// shadow, which must say `identical`).
fn take_shadow_line(run: &mut Run) -> Shadow {
    let text = String::from_utf8_lossy(&run.stderr).into_owned();
    let lines: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| line.starts_with("[X_SHADOW] QSTREAM "))
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
    let shadow = Shadow {
        batches: field("batches=").parse().unwrap(),
        records: field("records=").parse().unwrap(),
        slot_copies: field("slot_copies=").parse().unwrap(),
        off: field("off="),
    };
    run.stderr = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_SHADOW] QSTREAM "))
        .collect::<String>()
        .into_bytes();
    shadow
}

/// The run without the switch, with it and with the shadow (with `extra` switches in all three):
/// the same status, report and standard error. Returns the off run and the shadow's numbers.
fn same_three_with_input(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
    stdin: Option<&[u8]>,
) -> (Run, Shadow) {
    let off = blastn_with_input(args, threads, extra, out, stdin);
    let mut with = extra.to_vec();
    with.push(SWITCH);
    let on = blastn_with_input(args, threads, &with, out, stdin);
    assert!(
        on == off,
        "{what} threads {threads} {extra:?}: on differs\n{on:?}\n{off:?}"
    );
    let mut with = extra.to_vec();
    with.push(SHADOW);
    let mut shadow = blastn_with_input(args, threads, &with, out, stdin);
    let numbers = take_shadow_line(&mut shadow);
    assert!(
        shadow == off,
        "{what} threads {threads} {extra:?}: shadow differs\n{shadow:?}\n{off:?}"
    );
    (off, numbers)
}

fn same_three(
    what: &str,
    args: &[&str],
    threads: u32,
    extra: &[(&str, &str)],
    out: Option<&Path>,
) -> (Run, Shadow) {
    same_three_with_input(what, args, threads, extra, out, None)
}

/// The shadow read every batch again: about a dozen batches, each record compared (the 19 input
/// records at least once).
fn read_batch_by_batch(what: &str, shadow: &Shadow) {
    assert!(
        shadow.off == "none" && shadow.batches >= 8 && shadow.records >= 19,
        "{what}: {shadow:?}"
    );
}

/// Every task, outfmt 0, 6 and 7, 1 and 4 threads, the reports after the loop and per batch
/// (LOSAT_X_RPTSTREAM): the off bytes, with every batch read again and the warnings of the reader
/// and of the invalid queries on standard error.
#[test]
fn every_task_format_thread_count_and_report_mode_gives_the_off_bytes() {
    let dir = TempDir::new("formats");
    let query = queries(&dir, "q.fa", b"");
    let subject = fasta("LC738868.fasta");
    for task in TASKS {
        for outfmt in ["0", "6", "7"] {
            for threads in [1, 4] {
                for reports in REPORTS {
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
                    let (off, shadow) = same_three(&what, &args, threads, reports, None);
                    assert_eq!(off.status, Some(0), "{what}");
                    assert!(!off.report.is_empty(), "{what}: no report");
                    let stderr = String::from_utf8_lossy(&off.stderr);
                    assert!(
                        stderr.contains("FASTA-Reader: Ignoring invalid residues")
                            && stderr.contains("Sequence contains no data"),
                        "{what}: {stderr}"
                    );
                    read_batch_by_batch(&what, &shadow);
                    assert_eq!(shadow.slot_copies, 0, "{what}");
                }
            }
        }
    }
}

/// `-out` files, LOSAT_X_QBATCHPAR's slots (two slots of two threads, `LOSAT_X_QBATCHPARSLOTS=2`,
/// which copy the records of the batches they search ahead of the loop), and LOSAT_X_RPTSTREAM's
/// shadow (the records of every batch kept for its comparison): the off bytes; the slots did
/// search batches with records read again.
#[test]
fn out_files_qbatchpar_slots_and_the_report_shadow_give_the_off_bytes() {
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
            let (off, shadow) = same_three(&what, &args, 4, &[RPTSTREAM], Some(&out));
            assert!(!off.report.is_empty(), "{what}");
            read_batch_by_batch(&what, &shadow);
            let what = format!("{task} outfmt {outfmt} RPTSTREAMSHADOW");
            let (off, shadow) = same_three(&what, &args, 4, &[RPTSTREAMSHADOW], None);
            assert!(!off.report.is_empty(), "{what}");
            read_batch_by_batch(&what, &shadow);
            let stderr = String::from_utf8_lossy(&off.stderr);
            assert!(
                stderr.contains("[X_SHADOW] RPTSTREAM ") && stderr.contains(" off=none identical"),
                "{what}: {stderr}"
            );
            for reports in REPORTS {
                let mut extra = vec![QBATCHPAR, TWO_SLOTS];
                extra.extend_from_slice(reports);
                let what = format!("{task} outfmt {outfmt} QBATCHPAR {reports:?}");
                let (off, shadow) = same_three(&what, &args, 4, &extra, None);
                assert!(!off.report.is_empty(), "{what}");
                read_batch_by_batch(&what, &shadow);
                assert!(shadow.slot_copies > 0, "{what}: {shadow:?}");
                // The slots searched batches ahead of the loop, with copies of records read again.
                extra.extend_from_slice(&[SWITCH, STATS]);
                let stats = blastn(&args, 4, &extra, None);
                let stderr = String::from_utf8_lossy(&stats.stderr);
                let number = |prefix: &str, name: &str| -> u64 {
                    stderr
                        .lines()
                        .find(|line| line.starts_with(prefix))
                        .and_then(|line| {
                            line.split_whitespace()
                                .find_map(|word| word.strip_prefix(name))
                        })
                        .and_then(|value| value.parse().ok())
                        .unwrap_or_else(|| panic!("{what}: {stderr}"))
                };
                assert!(
                    number("[X_STATS] QBATCHPAR ", "precomputed=") > 0,
                    "{what}: {stderr}"
                );
                assert!(
                    number("[X_STATS] QSTREAM ", "slot_copies=") > 0
                        && number("[X_STATS] QSTREAM ", "read=")
                            == number("[X_STATS] QSTREAM ", "records="),
                    "{what}: {stderr}"
                );
            }
        }
    }
}

/// `-query_loc` (a range that skips the records without letters and the two-letter one, and one
/// that makes the two-letter record a query without data), `-subject_loc` with `-lcase_masking`
/// (queries with lower case runs), and `-max_target_seqs`, `-max_hsps` and `-subject_besthit`:
/// the off bytes, with the reports after the loop and per batch.
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
    let cases: [(&str, &[&str]); 5] = [
        ("blastn", &["-outfmt", "0", "-query_loc", "100-1500"]),
        ("megablast", &["-outfmt", "6", "-query_loc", "3-1000"]),
        (
            "dc-megablast",
            &["-outfmt", "7", "-query_loc", "100-1500", "-lcase_masking"],
        ),
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
            for reports in REPORTS {
                let what = format!("{task} {extra:?}");
                let (off, shadow) = same_three(&what, &args, threads, reports, None);
                assert_eq!(off.status, Some(0), "{what}");
                assert!(!off.report.is_empty(), "{what}");
                read_batch_by_batch(&what, &shadow);
            }
        }
    }
}

/// The errors of a later batch that come after the reports of the batches before (a reader error
/// late in the file; a batch of queries without letters; with `-query_loc`, a batch of skipped
/// records only), a reader error in the first record (before any output), and LOSAT's refusals,
/// which come before every report (a subject id that a later batch reaches, scoring without
/// Karlin-Altschul values, a negative `CHUNK_SIZE`): the off bytes, status and standard error.
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
                for reports in REPORTS {
                    let what = format!("{task} outfmt {outfmt} {query:?}");
                    let (off, shadow) = same_three(&what, &args, threads, reports, None);
                    assert_eq!(off.status, Some(status), "{what}");
                    assert!(
                        String::from_utf8_lossy(&off.stderr).contains(error),
                        "{what}: {}",
                        String::from_utf8_lossy(&off.stderr)
                    );
                    if before {
                        assert!(!off.report.is_empty(), "{what}");
                        read_batch_by_batch(&what, &shadow);
                    } else {
                        // Only the outfmt 0 prolog comes before the error (written before the
                        // first batch is read).
                        assert!(
                            outfmt == "0" || off.report.is_empty(),
                            "{what}: {}",
                            String::from_utf8_lossy(&off.report)
                        );
                        assert_eq!((shadow.off.as_str(), shadow.records), ("none", 0));
                    }
                }
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
        "megablast",
        "-query",
        skipped.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
        "-query_loc",
        "100-1500",
    ];
    for reports in REPORTS {
        let (off, shadow) = same_three("skipped", &args, 4, reports, None);
        assert_eq!(
            off.status,
            Some(3),
            "{:?}",
            String::from_utf8_lossy(&off.stderr)
        );
        assert!(String::from_utf8_lossy(&off.stderr).contains("Empty CBlastQueryVector"));
        read_batch_by_batch("skipped", &shadow);
    }
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
            "megablast",
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
    let ids = plain(&refused);
    let mut scoring = plain(&subject);
    scoring[1] = "blastn".to_string();
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
    for reports in REPORTS {
        for (what, args, env) in [
            ("Subject_ id", &ids, &[][..]),
            ("scoring", &scoring, &[][..]),
            ("CHUNK_SIZE", &plain(&subject), &chunks[..]),
        ] {
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let mut extra = env.to_vec();
            extra.extend_from_slice(reports);
            let (off, shadow) = same_three(what, &args, 4, &extra, None);
            assert_ne!(off.status, Some(0), "{what}: the refusal");
            assert!(
                off.report.is_empty(),
                "{what}: no report before the refusal"
            );
            assert_eq!(shadow.off, "none", "{what}");
        }
    }
}

/// Standard input (`-query -`) and a pipe (`-query /dev/stdin` with standard input a pipe) are
/// read once, before the search, with the switch too: the off bytes, and the shadow line names
/// why.
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
        for reports in REPORTS {
            let (off, shadow) = same_three_with_input(path, &args, 4, reports, None, Some(&text));
            assert_eq!(off.status, Some(0), "{path}");
            assert!(!off.report.is_empty(), "{path}");
            assert_eq!(
                (shadow.off.as_str(), shadow.batches, shadow.records),
                (reason, 0, 0),
                "{path}"
            );
        }
    }
}

/// A database (`-db`, made by `LOSAT makedb` from four genomes): the off bytes, every format,
/// the reports after the loop and per batch.
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
            for reports in REPORTS {
                let what = format!("-db {task} outfmt {outfmt}");
                let (off, shadow) = same_three(&what, &args, threads, reports, None);
                assert_eq!(off.status, Some(0), "{what}");
                assert!(!off.report.is_empty(), "{what}");
                read_batch_by_batch(&what, &shadow);
            }
        }
    }
}
