//! TBLASTX's report stage on the search's pool (`LOSAT_X_PAIRWISEPAR`,
//! `src/algorithm/tblastx/x_rpt_par.rs`) on whole searches. The stage was under a switch of its
//! own, `LOSAT_X_TBXRPTPAR`, until 2026-11-02 (structural line S3), when it was folded into
//! `LOSAT_X_PAIRWISEPAR`, which keeps the whole report stage serial at `0`.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_PAIRWISEPAR`
//! is `0`. With it on and two or more threads, the displayed
//! identities, the outfmt 0 rows and the tabular rows of the final list are computed by the pool's
//! threads (the Seq-align engine counts of outfmt 8 and 12 too). The reach rows (`LOSAT_X_REACH`,
//! `tests/reach_support`) show that the pooled stage ran.
//!
//! The property under test is that the reports, the messages and the exit status do not depend on
//! the switch or on the thread count: the searches are compared with the same LOSAT search with the
//! switch off at one thread. The first query below has more than 2,048 rows (AP027131 x AP027133,
//! 120 kb each, genetic code 4 for both, which a local subject honors: the approved deviation of
//! `docs/release/tblastx_v0.1.0_certification.md`), so the identities, the outfmt 0 rows and the
//! tabular rows are cut into several chunks on the pool. The tabular rows take BLASTN's pieces
//! (512 rows of one query) and windows (2,048 rows per thread) for a query of any number of rows
//! (until 2026-11-02 only a query of more than 1,024 rows was formatted on the pool): a second
//! case below has several queries of 1,024 rows or fewer.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_PAIRWISEPAR", "0");
const THREADS: [u32; 5] = [1, 2, 4, 8, 16];

fn letters(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    std::fs::read(path)
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(&str, &[u8])]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-pwpar-tbx-{}", std::process::id()));
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

/// Two queries from AP027131 (the second with a title of several words) and a record without
/// letters between them, which has a warning before its report; written once.
fn queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let a = letters("AP027131.fasta");
        write_fasta(
            "queries.fa",
            &[
                ("a", &a[..120_000]),
                ("empty record", b""),
                ("a second part", &a[400_000..430_000]),
            ],
        )
    })
}

/// A subject from AP027133 (a title of several words), written once (the tests run in parallel).
fn subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let b = letters("AP027133.fasta");
        write_fasta("subjects.fa", &[("b first part", &b[..120_000])])
    })
}

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path
/// (`LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread, and it selects nothing).
fn tblastx(threads: u32, outfmt: &str, switches: &[(&str, &str)]) -> Output {
    tblastx_of(queries(), threads, outfmt, switches)
}

/// `tblastx` of the queries in `query` against the subject.
fn tblastx_of(query: &Path, threads: u32, outfmt: &str, switches: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects())
        .args(["-query_gencode", "4", "-db_gencode", "4"])
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn check_format(outfmt: &str) {
    let off = tblastx(1, outfmt, &[OFF]);
    assert!(off.status.success(), "{}", stderr(&off));
    assert!(!off.stdout.is_empty(), "the pairs have hits");
    // LOSAT writes no warning with outfmt 8.
    assert!(
        outfmt == "8" || stderr(&off).contains("empty record"),
        "the empty record has a warning: {}",
        stderr(&off)
    );
    if matches!(outfmt, "6" | "7") {
        // The first query's rows are formatted in windows of chunks (more than 1,024 rows).
        let rows = String::from_utf8_lossy(&off.stdout)
            .lines()
            .filter(|line| line.starts_with("a\t"))
            .count();
        assert!(
            rows > 2_048,
            "outfmt {outfmt}: {rows} rows of the first query"
        );
    }
    for threads in THREADS {
        let reach = reach_support::Reach::new("rptpar");
        let on = tblastx(threads, outfmt, &[reach.var()]);
        assert!(on.status.success(), "threads {threads}: {}", stderr(&on));
        assert_eq!(on.stdout, off.stdout, "outfmt {outfmt}, threads {threads}");
        assert_eq!(
            stderr(&on),
            stderr(&off),
            "outfmt {outfmt}, threads {threads}"
        );
        // The pooled stage runs with two threads or more, not with one.
        let sides = reach.sides("LOSAT_X_PAIRWISEPAR");
        assert_eq!(
            sides.0,
            threads > 1,
            "outfmt {outfmt}, threads {threads}: sides {sides:?}"
        );
    }
}

/// Several short queries from AP027131 (each with 1,024 rows or fewer, some with none), written
/// once: their rows were written serially before 2026-11-02.
fn small_queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let a = letters("AP027131.fasta");
        let records: Vec<(String, &[u8])> = (0..8)
            .map(|i| (format!("s{i}"), &a[i * 15_000..i * 15_000 + 9_000]))
            .collect();
        let records: Vec<(&str, &[u8])> = records
            .iter()
            .map(|(title, seq)| (title.as_str(), *seq))
            .collect();
        write_fasta("small_queries.fa", &records)
    })
}

/// The rows of each query of an outfmt 6 / 7 report, in order of appearance.
fn rows_per_query(stdout: &[u8]) -> Vec<(String, usize)> {
    let mut counts: Vec<(String, usize)> = Vec::new();
    for line in String::from_utf8_lossy(stdout).lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let query = line.split('\t').next().unwrap();
        match counts.last_mut() {
            Some((last, n)) if last == query => *n += 1,
            _ => counts.push((query.to_string(), 1)),
        }
    }
    counts
}

fn check_small_queries(outfmt: &str) {
    let off = tblastx_of(small_queries(), 1, outfmt, &[OFF]);
    assert!(off.status.success(), "{}", stderr(&off));
    let counts = rows_per_query(&off.stdout);
    assert!(
        counts.len() > 1,
        "outfmt {outfmt}: {} queries with rows",
        counts.len()
    );
    assert!(
        counts.iter().all(|(_, rows)| *rows <= 1_024),
        "outfmt {outfmt}: rows per query {counts:?}"
    );
    assert!(
        counts.iter().map(|(_, rows)| rows).sum::<usize>() > 100,
        "outfmt {outfmt}: rows per query {counts:?}"
    );
    for threads in [1u32, 2, 4] {
        let reach = reach_support::Reach::new("rptpar-small");
        let on = tblastx_of(small_queries(), threads, outfmt, &[reach.var()]);
        assert!(on.status.success(), "threads {threads}: {}", stderr(&on));
        assert_eq!(on.stdout, off.stdout, "outfmt {outfmt}, threads {threads}");
        assert_eq!(
            stderr(&on),
            stderr(&off),
            "outfmt {outfmt}, threads {threads}"
        );
        let sides = reach.sides("LOSAT_X_PAIRWISEPAR");
        assert_eq!(
            sides.0,
            threads > 1,
            "outfmt {outfmt}, threads {threads}: sides {sides:?}"
        );
    }
}

#[test]
fn the_rows_of_queries_with_few_rows_do_not_depend_on_the_switch() {
    check_small_queries("6");
    check_small_queries("7");
}

#[test]
fn the_tabular_report_does_not_depend_on_the_switch_at_any_thread_count() {
    check_format("6");
}

#[test]
fn the_commented_tabular_report_does_not_depend_on_the_switch() {
    check_format("7");
}

#[test]
fn the_pairwise_report_does_not_depend_on_the_switch() {
    check_format("0");
}

#[test]
fn the_seq_align_report_does_not_depend_on_the_switch() {
    check_format("8");
}
