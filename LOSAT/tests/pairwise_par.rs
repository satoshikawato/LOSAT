//! The pairwise reports (outfmt 0) formatted on the search's pool (`LOSAT_X_PAIRWISEPAR`,
//! `src/report/x_pairwise_par.rs`) on whole searches: BLASTN (the blastn and megablast tasks,
//! `-subject` and `-db`) and TBLASTX (`-subject` and `-db`).
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_PAIRWISEPAR`
//! is `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it on and two or more threads, each
//! query's head, runs of its subjects' alignments and its footer are formatted into buffers by the
//! pool's threads and written in order; `LOSAT_X_PAIRWISEPARSHADOW=1` also writes the report by the
//! serial writer into memory and compares the two.
//!
//! The property under test is that the reports, the messages and the exit status do not depend on
//! the switch or on the thread count: each search is compared with the same LOSAT search with the
//! switch off at one thread. The inputs have a record without letters between queries (a warning
//! before its report), queries with hits and without, several subjects, a query whose alignments
//! are cut into several pieces, and (BLASTN at two threads) more than one window of queries.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const OFF: (&str, &str) = ("LOSAT_X_PAIRWISEPAR", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRWISEPARSHADOW", "1");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const THREADS: [u32; 4] = [1, 2, 4, 16];

fn manifest(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn letters(name: &str) -> Vec<u8> {
    std::fs::read(manifest("tests/fasta").join(name))
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(String, &[u8])]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-pairwisepar-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let path = dir.join(file);
    let mut text = Vec::new();
    for (title, seq) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).expect("write FASTA");
    path
}

/// BLASTN queries: 30 parts of 4 kb of MelaMJNV (titles of several words), a record without
/// letters after the fifth, and a part of AP027131 that has no hits in the subjects; written
/// once (the tests run in parallel).
fn nucleotide_queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let mela = letters("MelaMJNV.fasta");
        let other = letters("AP027131.fasta");
        let mut records: Vec<(String, &[u8])> = Vec::new();
        for part in 0..30 {
            records.push((
                format!("mela part {part}"),
                &mela[part * 4_000..part * 4_000 + 4_000],
            ));
            if part == 4 {
                records.push(("empty record".to_string(), b""));
            }
        }
        records.push(("no hits".to_string(), &other[..3_000]));
        write_fasta("nucleotide_queries.fa", &records)
    })
}

/// Two BLASTN subjects (MejoMJNV and MeenMJNV), written once.
fn nucleotide_subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let mejo = letters("MejoMJNV.fasta");
        let meen = letters("MeenMJNV.fasta");
        write_fasta(
            "nucleotide_subjects.fa",
            &[
                ("mejo genome".to_string(), &mejo[..]),
                ("meen genome".to_string(), &meen[..]),
            ],
        )
    })
}

/// TBLASTX queries: a part of AP027131 with more than 2,048 HSPs (several pieces), a record
/// without letters and a second part; written once.
fn tblastx_queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let a = letters("AP027131.fasta");
        write_fasta(
            "tblastx_queries.fa",
            &[
                ("a".to_string(), &a[..120_000]),
                ("empty record".to_string(), b""),
                ("a second part".to_string(), &a[400_000..430_000]),
            ],
        )
    })
}

/// A TBLASTX subject from AP027133, written once.
fn tblastx_subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let b = letters("AP027133.fasta");
        write_fasta(
            "tblastx_subjects.fa",
            &[("b first part".to_string(), &b[..120_000])],
        )
    })
}

/// A search with the given variables as its only `LOSAT_*` variables that choose a code path
/// (`LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread, and it selects nothing).
fn losat(args: &[&str], threads: u32, switches: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
        .envs(switches.iter().copied())
        .args(args)
        .args(["-outfmt", "0", "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The standard error without the shadow's lines.
fn messages(output: &Output) -> String {
    stderr(output)
        .lines()
        .filter(|line| !line.starts_with("[X_SHADOW]"))
        .map(|line| format!("{line}\n"))
        .collect()
}

/// The `[X_SHADOW] PAIRWISEPAR` lines' counters (queries, chunks, pieces, bytes), one per report.
fn shadow_counts(output: &Output) -> Vec<[u64; 4]> {
    stderr(output)
        .lines()
        .filter(|line| line.starts_with("[X_SHADOW] PAIRWISEPAR "))
        .map(|line| {
            assert!(line.ends_with(" identical"), "{line}");
            ["queries", "chunks", "pieces", "bytes"].map(|name| {
                let key = format!("{name}=");
                line.split_whitespace()
                    .find_map(|field| field.strip_prefix(key.as_str()))
                    .unwrap_or_else(|| panic!("no {name} in {line}"))
                    .parse()
                    .expect("a count")
            })
        })
        .collect()
}

/// Checks the property for one search; returns the shadow's counters at two threads.
fn check(name: &str, args: &[&str]) -> Vec<[u64; 4]> {
    let off = losat(args, 1, &[OFF]);
    assert!(off.status.success(), "{name}: {}", stderr(&off));
    assert!(!off.stdout.is_empty(), "{name}");
    let mut at_two = Vec::new();
    for threads in THREADS {
        let on = losat(args, threads, &[]);
        assert!(
            on.status.success(),
            "{name}, threads {threads}: {}",
            stderr(&on)
        );
        assert_eq!(on.stdout, off.stdout, "{name}, threads {threads}");
        assert_eq!(messages(&on), messages(&off), "{name}, threads {threads}");
        let shadow = losat(args, threads, &[SHADOW]);
        assert!(
            shadow.status.success(),
            "{name}, threads {threads}: {}",
            stderr(&shadow)
        );
        assert_eq!(
            shadow.stdout, off.stdout,
            "{name}, threads {threads}, shadow"
        );
        assert_eq!(
            messages(&shadow),
            messages(&off),
            "{name}, threads {threads}"
        );
        let counts = shadow_counts(&shadow);
        if threads == 1 {
            // One thread: the serial writer, as before.
            assert!(counts.is_empty(), "{name}: {}", stderr(&shadow));
            continue;
        }
        assert!(
            !counts.is_empty(),
            "{name}, threads {threads}: {}",
            stderr(&shadow)
        );
        if args[0] == "blastn" {
            // BLASTN's rows, made on the pool, compared with the serial rows.
            assert!(
                stderr(&shadow).lines().any(|line| {
                    line.starts_with("[X_SHADOW] PAIRWISEPARROWS rows=")
                        && line.ends_with(" identical")
                }),
                "{name}, threads {threads}: {}",
                stderr(&shadow)
            );
        }
        assert_eq!(
            counts.iter().map(|[_, _, _, bytes]| bytes).sum::<u64>() as usize
                + prolog_and_epilog(&off.stdout),
            off.stdout.len(),
            "{name}, threads {threads}: the shadow compared every query's report"
        );
        if threads == 2 {
            at_two = counts;
        }
    }
    let previous = losat(args, 4, &[PREVIOUS]);
    assert!(previous.status.success(), "{name}: {}", stderr(&previous));
    assert_eq!(previous.stdout, off.stdout, "{name}, LOSAT_X_DEFAULT=0");
    at_two
}

/// The bytes before the first query's report (`Query=` after two newlines) and after the last
/// query's footer (the epilog, from the two newlines before `  Database:` or `Lambda`).
fn prolog_and_epilog(report: &[u8]) -> usize {
    let text = String::from_utf8_lossy(report);
    let first = text.find("\n\nQuery=").expect("a query's report");
    let tail = text
        .rfind("\n\n  Database:")
        .or_else(|| text.rfind("\n\nLambda"))
        .expect("the epilog");
    first + (text.len() - tail)
}

fn nucleotide(task: &str) -> Vec<String> {
    vec![
        "blastn".to_string(),
        "-task".to_string(),
        task.to_string(),
        "-query".to_string(),
        nucleotide_queries().display().to_string(),
        "-subject".to_string(),
        nucleotide_subjects().display().to_string(),
    ]
}

fn strings(args: &[String]) -> Vec<&str> {
    args.iter().map(String::as_str).collect()
}

#[test]
fn the_blastn_pairwise_report_does_not_depend_on_the_switch() {
    let args = nucleotide("blastn");
    let counts = check("blastn", &strings(&args));
    // The 32 queries in the reports of the NCBI query batches; a batch's report formatted in more
    // than one chunk at two threads.
    assert_eq!(
        counts.iter().map(|[queries, ..]| queries).sum::<u64>(),
        32,
        "{counts:?}"
    );
    assert!(
        counts.iter().any(|&[_, chunks, ..]| chunks > 1),
        "{counts:?}"
    );
}

#[test]
fn the_megablast_pairwise_report_does_not_depend_on_the_switch() {
    let args = nucleotide("megablast");
    check("megablast", &strings(&args));
}

#[test]
fn fewer_alignments_than_subjects_do_not_depend_on_the_switch() {
    let mut args = nucleotide("blastn");
    args.extend(["-max_target_seqs".to_string(), "1".to_string()]);
    check("blastn -max_target_seqs 1", &strings(&args));
}

#[test]
fn the_blastn_database_report_does_not_depend_on_the_switch() {
    let db = manifest("tests/db_search/nslice");
    let query = manifest("tests/db_search/nslice_q.fna");
    for task in ["blastn", "megablast"] {
        let args = [
            "blastn",
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-db",
            db.to_str().unwrap(),
        ];
        check(&format!("blastn -db {task}"), &args);
    }
}

#[test]
fn the_tblastx_pairwise_report_does_not_depend_on_the_switch() {
    let queries = tblastx_queries().display().to_string();
    let subjects = tblastx_subjects().display().to_string();
    let args = [
        "tblastx",
        "-query",
        &queries,
        "-subject",
        &subjects,
        "-query_gencode",
        "4",
        "-db_gencode",
        "4",
    ];
    let counts = check("tblastx", &args);
    // The first query's alignments with its one subject are cut into several pieces (three
    // heads, two footers and the second query's alignments make 6 pieces).
    let [queries, _, pieces, _] = counts[0];
    assert_eq!((counts.len(), queries), (1, 3), "{counts:?}");
    assert!(pieces > 7, "{counts:?}");
}

#[test]
fn the_tblastx_database_report_does_not_depend_on_the_switch() {
    let db = manifest("tests/db_search/nslice");
    let query = manifest("tests/db_search/nslice_q.fna");
    let args = [
        "tblastx",
        "-query",
        query.to_str().unwrap(),
        "-db",
        db.to_str().unwrap(),
    ];
    check("tblastx -db", &args);
}

#[test]
fn a_failed_write_ends_the_search_as_without_the_switch() {
    if !Path::new("/dev/full").exists() {
        return;
    }
    let mut args = nucleotide("megablast");
    args.extend(["-out".to_string(), "/dev/full".to_string()]);
    let args = strings(&args);
    let off = losat(&args, 4, &[OFF]);
    assert!(!off.status.success(), "{}", stderr(&off));
    for switches in [&[][..], &[SHADOW][..]] {
        let on = losat(&args, 4, switches);
        assert_eq!(on.status.code(), off.status.code(), "{}", stderr(&on));
        assert_eq!(messages(&on), messages(&off));
        assert_eq!(on.stdout, off.stdout);
    }
}
