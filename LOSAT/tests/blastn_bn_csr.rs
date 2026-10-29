//! The chains of a contiguous blastn lookup table laid out as CSR (`LOSAT_X_BNCSR`,
//! `src/algorithm/blastn/x_bn_csr.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It chooses between two LOSAT layouts of the same chains (the runs of a cell list the chain's
//! query offsets in the chain's order), so the property under test is that the report does not
//! depend on it, and that the default took the CSR where its rule says so: `LOSAT_X_REACH` shows
//! the switch's own side, and `LOSAT_TIMING`'s `[TIMING] x_bn_csr: builds= parallel=` line counts
//! the CSR tables and those built on the search pool (all of them at four threads, none at one).
//! The packed head of a newest-first table (`LOSAT_X_BNCSRPACK`, strict set) is on in the default
//! too: its `=0` gives the same report, and the reading-order table never takes it.
//! Two tables: `-task blastn` (word size 11) with a query of 1.27 Mb (a 4^11-cell table, newest
//! first, above the rule), and `-word_size 9` with a query of 9,800 letters (the 4^8-cell table
//! that stands for NCBI's small table, in reading order, above the rule), each against a subject
//! file (`-subject`) and a database (`-db`), at one and four threads.
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_BNCSR";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNCSR", "0");
/// The packed head of a newest-first CSR table (`LOSAT_X_BNCSRPACK`, strict set).
const PACK: &str = "LOSAT_X_BNCSRPACK";
const PACK_OFF: (&str, &str) = ("LOSAT_X_BNCSRPACK", "0");
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
        let dir = std::env::temp_dir().join(format!("losat-bncsr-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// The letters of the first record of a FASTA file.
fn letters(path: &Path) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// A database of `subject` made with `losat makeblastdb` in `dir`.
fn database(dir: &TestDir, subject: &Path, name: &str) -> PathBuf {
    let db = dir.path(name);
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(subject)
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
    db
}

/// One search of `query` against `subjects` (`-subject` or `-db`, `-outfmt 6` on standard output)
/// with `extra` arguments and the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    subjects: (&str, &Path),
    extra: &[&str],
    threads: u32,
    vars: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg(subjects.0)
        .arg(subjects.1)
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {} {vars:?} threads {threads}: {stderr}",
        subjects.0
    );
    (output.stdout, stderr)
}

/// `(builds, parallel)` of `[TIMING] x_bn_csr: builds=<n> parallel=<n> words=<n>` in `stderr`
/// (`None` without the line).
fn csr_builds(stderr: &str) -> Option<(u64, u64)> {
    let line = stderr
        .lines()
        .find_map(|line| line.strip_prefix("[TIMING] x_bn_csr: "))?;
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    Some((field("builds="), field("parallel=")))
}

/// For `query` and `extra` against `subject` as a file and as a database: the report with the
/// switch off at one thread, then off, by default, by default with `LOSAT_TIMING` and with the
/// packed head off (`LOSAT_X_BNCSRPACK=0`) at 1 and 4 threads, all the same bytes; the default takes
/// the switch's side (`LOSAT_X_REACH`), and the packed head's when `packed` (a newest-first table)
/// or its reference side (the reading-order table), builds every CSR table on the pool at four
/// threads and on the search thread at one, and writes nothing else.
fn check(dir: &TestDir, query: &Path, subject: &Path, extra: &[&str], packed: bool) {
    let db = database(dir, subject, "subject_db");
    for subjects in [("-subject", subject), ("-db", db.as_path())] {
        let what = format!("{extra:?} {}", subjects.0);
        let (expected, expected_stderr) = search(query, subjects, extra, 1, &[SWITCH_OFF]);
        assert!(!expected.is_empty(), "{what}: no hits");
        for threads in [1, 4] {
            let reach = reach_support::Reach::new("bncsr");
            let (report, stderr) = search(query, subjects, extra, threads, &[reach.var()]);
            assert!(
                report == expected,
                "{what} threads {threads}: default differs"
            );
            assert_eq!(stderr, expected_stderr, "{what} threads {threads}");
            assert!(
                reach.sides(SWITCH).0,
                "{what} threads {threads}: no table laid out as CSR"
            );
            assert_eq!(
                reach.sides(PACK),
                (packed, !packed),
                "{what} threads {threads}: packed heads"
            );
            let (report, stderr) = search(query, subjects, extra, threads, &[PACK_OFF]);
            assert!(
                report == expected,
                "{what} threads {threads}: packed head off differs"
            );
            assert_eq!(stderr, expected_stderr, "{what} threads {threads} pack off");
            let (report, stderr) = search(query, subjects, extra, threads, &[TIMING]);
            assert!(
                report == expected,
                "{what} threads {threads}: default with LOSAT_TIMING differs"
            );
            let (builds, parallel) = csr_builds(&stderr)
                .unwrap_or_else(|| panic!("{what} threads {threads}: no CSR line:\n{stderr}"));
            assert!(builds > 0, "{what} threads {threads}");
            let on_pool = if threads > 1 { builds } else { 0 };
            assert_eq!(parallel, on_pool, "{what} threads {threads}: pool builds");
            if threads > 1 {
                let (report, stderr) = search(query, subjects, extra, threads, &[SWITCH_OFF]);
                assert!(report == expected, "{what} threads {threads}: off differs");
                assert_eq!(stderr, expected_stderr, "{what} threads {threads} off");
            }
        }
        let (_, stderr) = search(query, subjects, extra, 1, &[SWITCH_OFF, TIMING]);
        assert_eq!(
            csr_builds(&stderr),
            None,
            "{what}: CSR tables with the switch off"
        );
    }
}

// The 4^11-cell table of `-task blastn` (newest first) for one query record of AP027131 followed
// by AP027133 (1.27 Mb, one query batch, 1,739,726 words after dust: above the rule's 1,048,576;
// either genome alone is below it, and as two records they are two batches), against 60 kb of
// AP027132.
#[test]
fn blastn_word_11_large_query_batch() {
    let dir = TestDir::new("w11");
    let query = dir.path("ap131_ap133.fna");
    let text = format!(
        ">ap131_ap133\n{}{}\n",
        letters(&fasta("AP027131.fasta")),
        letters(&fasta("AP027133.fasta"))
    );
    std::fs::write(&query, text).expect("write the query");
    let subject = dir.path("ap132_60k.fna");
    let genome = letters(&fasta("AP027132.fasta"));
    std::fs::write(
        &subject,
        format!(">ap132_60k\n{}\n", &genome[100_000..160_000]),
    )
    .expect("write the subject");
    check(&dir, &query, &subject, &["-task", "blastn"], true);
}

// The 4^8-cell table of `-word_size 9` that stands for NCBI's small table (reading order): 9,800
// letters of LC738874, both strands above the rule's 16,384 positions and below the 21,000 table
// entries where NCBI's choice moves to a megablast-type table, against LC738868.
#[test]
fn blastn_word_9_reading_order_table() {
    let dir = TestDir::new("w9");
    let query = dir.path("lc874_9800.fna");
    let genome = letters(&fasta("LC738874.fasta"));
    std::fs::write(
        &query,
        format!(">lc874_9800\n{}\n", &genome[20_000..29_800]),
    )
    .expect("write the query");
    check(
        &dir,
        &query,
        &fasta("LC738868.fasta"),
        &["-task", "blastn", "-word_size", "9"],
        false,
    );
}
