// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! The run-based one-thread fill of the megablast table (`LOSAT_X_MBFILL`,
//! `src/algorithm/blastn/x_mb_fill.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. It chooses between two LOSAT code
//! paths that build the same table, so the property under test is that the report does not
//! depend on it: the default run, `LOSAT_X_MBFILL=0`, `LOSAT_X_DEFAULT=0` and the shadow run
//! (which builds each one-thread table again with the queued loop and stops on the first
//! difference) give the same bytes, and the default run prints nothing on standard error (the
//! `[X_SHADOW] MBFILL` line is printed only with `LOSAT_X_MBFILLSHADOW`). Two searches build the
//! table on one thread: megablast (the contiguous megablast table) and blastn (word size 11 on a
//! short query: the two-stage table with its cells in the order of indexing) of a query cut from
//! LC738874 with runs of N, against LC738868; at four threads the pool builds the table and the
//! report is the same. That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_MBFILL", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_MBFILLSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-mbfill-{test}-{}", std::process::id()));
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

/// The two queries: the first 60,000 letters of LC738874 with runs of 1 to 7 N every 7,000
/// letters (words stop at them), for megablast; its first 9,000 letters, for blastn.
fn queries(dir: &TestDir) -> [(&'static str, PathBuf); 2] {
    let genome = letters(&fasta("LC738874.fasta"));
    let mut long: Vec<u8> = genome.as_bytes()[..60_000].to_vec();
    for (k, start) in (3_000..60_000).step_by(7_000).enumerate() {
        for letter in &mut long[start..start + 1 + k % 7] {
            *letter = b'N';
        }
    }
    let long = String::from_utf8(long).unwrap();
    let mb = dir.0.join("q60000n.fna");
    std::fs::write(&mb, format!(">q60000n\n{long}\n")).expect("write the query");
    let bn = dir.0.join("q9000.fna");
    std::fs::write(&bn, format!(">q9000\n{}\n", &long[..9_000])).expect("write the query");
    [("megablast", mb), ("blastn", bn)]
}

/// One search of `query` against LC738868 (`-outfmt 6` on standard output) with the given
/// variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing. Returns the report and the standard error.
fn search(task: &str, query: &Path, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-task", task])
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

// The default one-thread search runs the run-based fill; `LOSAT_X_MBFILL=0` and
// `LOSAT_X_DEFAULT=0` run the queued loop, and at four threads the pool fills the table. The
// reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_run_based_fill() {
    let dir = TestDir::new("report");
    for (task, query) in queries(&dir) {
        let (default, default_stderr) = search(task, &query, 1, &[]);
        assert!(!default.is_empty(), "{task}: no hits");
        assert_eq!(default_stderr, "", "{task}");
        for off in [SWITCH_OFF, DEFAULT_OFF] {
            let (report, stderr) = search(task, &query, 1, &[off]);
            assert_eq!(report, default, "{task} {off:?}");
            assert_eq!(stderr, "", "{task} {off:?}");
        }
        let (report, stderr) = search(task, &query, 4, &[]);
        assert_eq!(report, default, "{task} threads 4");
        assert_eq!(stderr, "", "{task} threads 4");
    }
}

// `LOSAT_X_DEFAULT=0` with the variable set to 1 (the run-based fill again) gives the default
// report. The shadow (set to any value, not in the set; it implies the switch) builds each
// one-thread table again with the queued loop and compares every array, prints its summary line
// with a count above 0, and does not change the report.
#[test]
fn the_shadow_compares_the_two_fills_and_the_report_is_unchanged() {
    let dir = TestDir::new("shadow");
    for (task, query) in queries(&dir) {
        let (expected, _) = search(task, &query, 1, &[]);
        let (report, _) = search(task, &query, 1, &[DEFAULT_OFF, ("LOSAT_X_MBFILL", "1")]);
        assert_eq!(report, expected, "{task}");
        for switches in [&[SHADOW][..], &[DEFAULT_OFF, SHADOW][..]] {
            let (report, stderr) = search(task, &query, 1, switches);
            assert_eq!(report, expected, "{task} {switches:?}");
            let line = stderr
                .lines()
                .find(|line| line.starts_with("[X_SHADOW] MBFILL compared="))
                .unwrap_or_else(|| panic!("{task} {switches:?}: no summary line:\n{stderr}"));
            assert!(
                !line.starts_with("[X_SHADOW] MBFILL compared=0 "),
                "{task} {switches:?}: {line}"
            );
        }
    }
}
