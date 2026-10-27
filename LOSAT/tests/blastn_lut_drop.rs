//! A query batch's lookup tables dropped when its preliminary stage ends (`LOSAT_X_BNLUTDROP`,
//! `src/algorithm/blastn/blast_engine/x_lut_drop.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`. It changes when the tables are freed, not what the search computes: the traceback reads
//! the HSPs that the preliminary stage kept, and the tables are passed to the scan only (so the
//! compiler shows that the traceback does not read them). The property under test is that the
//! report does not depend on it and that the default run took the new path: with `LOSAT_TIMING`
//! the default run prints `[TIMING] x_bn_lut_drop: batches=` with a count above 0, and `=0`
//! prints no such line. Every task builds its own kind of table (megablast the contiguous
//! megablast table, blastn and blastn-short the small or the direct table, dc-megablast the
//! discontiguous tables), against a subject file (`-subject`) and a database (`-db`), at one and
//! four threads. That those bytes equal NCBI's is the business of the fixture suites.
//!
//! The same holds for how a batch's largest buffers are freed (`LOSAT_X_BNSHRINKFREE`,
//! `src/algorithm/blastn/x_shrink_free.rs`, also in the strict set): the CSR of a table and the
//! batch's BLASTNA query buffer are shrunk to one element before their free. With `LOSAT_TIMING`
//! the default run prints `[TIMING] x_bn_shrink_free: buffers=` with a count above 0, `=0` prints
//! no such line, and the reports are equal.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNLUTDROP", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const TASKS: [&str; 4] = ["megablast", "dc-megablast", "blastn", "blastn-short"];

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-lutdrop-{test}-{}", std::process::id()));
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

/// The query (the first 4,000 letters of LC738874) and a database of LC738868 made with
/// `losat makedb`, in `dir`.
fn inputs(dir: &TestDir) -> (PathBuf, PathBuf) {
    let genome = letters(&fasta("LC738874.fasta"));
    let query = dir.0.join("q4000.fna");
    std::fs::write(&query, format!(">q4000\n{}\n", &genome[..4_000])).expect("write the query");
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

/// One search of `query` against `subjects` (`-subject` or `-db`, `-outfmt 6` on standard output)
/// with the given variables as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays: a
/// binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI)
/// needs that file name as soon as a search has more than one thread, and it selects nothing.
/// Returns the report and the standard error.
fn search(
    task: &str,
    query: &Path,
    subjects: (&str, &Path),
    threads: u32,
    switches: &[(&str, &str)],
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
        .args(["-task", task])
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {} {switches:?} threads {threads}: {stderr}",
        subjects.0
    );
    (output.stdout, stderr)
}

/// The count of `[TIMING] x_bn_lut_drop: batches=` in `stderr` (`None` without the line).
fn dropped_batches(stderr: &str) -> Option<u64> {
    stderr.lines().find_map(|line| {
        line.strip_prefix("[TIMING] x_bn_lut_drop: batches=")
            .map(|count| count.trim().parse().expect("a count"))
    })
}

// The default drops each batch's tables after its preliminary stage and says so with
// LOSAT_TIMING; `LOSAT_X_BNLUTDROP=0` holds them to the end of the batch and prints no line. The
// reports are equal, and the same at four threads.
#[test]
fn the_report_does_not_depend_on_when_the_tables_are_dropped() {
    let dir = TestDir::new("report");
    let (query, db) = inputs(&dir);
    let subject = fasta("LC738868.fasta");
    for task in TASKS {
        for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
            let what = format!("{task} {}", subjects.0);
            let (default, stderr) = search(task, &query, subjects, 1, &[TIMING]);
            assert!(!default.is_empty(), "{what}: no hits");
            assert!(
                dropped_batches(&stderr).is_some_and(|count| count > 0),
                "{what}: no batch dropped its tables: {stderr}"
            );
            let (report, stderr) = search(task, &query, subjects, 1, &[TIMING, SWITCH_OFF]);
            assert_eq!(report, default, "{what} {SWITCH_OFF:?}");
            assert_eq!(dropped_batches(&stderr), None, "{what} {SWITCH_OFF:?}");
            let (report, stderr) = search(task, &query, subjects, 4, &[]);
            assert_eq!(report, default, "{what} threads 4");
            assert_eq!(stderr, "", "{what} threads 4");
        }
    }
}

/// The count of `[TIMING] x_bn_shrink_free: buffers=` in `stderr` (`None` without the line).
fn shrunk_buffers(stderr: &str) -> Option<u64> {
    stderr.lines().find_map(|line| {
        line.strip_prefix("[TIMING] x_bn_shrink_free: buffers=")
            .map(|count| count.trim().parse().expect("a count"))
    })
}

// The default shrinks each batch's query buffer (and its CSR, where the table is built as one)
// before the free and says so with LOSAT_TIMING; `LOSAT_X_BNSHRINKFREE=0` frees them as they are
// and prints no line. The reports are equal.
#[test]
fn the_report_does_not_depend_on_how_the_buffers_are_freed() {
    const SHRINK_OFF: (&str, &str) = ("LOSAT_X_BNSHRINKFREE", "0");
    let dir = TestDir::new("shrink");
    let (query, db) = inputs(&dir);
    let subject = fasta("LC738868.fasta");
    for task in TASKS {
        for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
            let what = format!("{task} {}", subjects.0);
            let (default, stderr) = search(task, &query, subjects, 1, &[TIMING]);
            assert!(!default.is_empty(), "{what}: no hits");
            assert!(
                shrunk_buffers(&stderr).is_some_and(|count| count > 0),
                "{what}: no buffer shrunk before its free: {stderr}"
            );
            let (report, stderr) = search(task, &query, subjects, 1, &[TIMING, SHRINK_OFF]);
            assert_eq!(report, default, "{what} {SHRINK_OFF:?}");
            assert_eq!(shrunk_buffers(&stderr), None, "{what} {SHRINK_OFF:?}");
        }
    }
}
