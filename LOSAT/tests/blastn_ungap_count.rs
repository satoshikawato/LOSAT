//! The counted loops of the exact ungapped extension (`LOSAT_X_BNUNGAPCOUNT`,
//! `src/algorithm/blastn/x_ungap_count.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. It chooses between two LOSAT code
//! paths with the same result, so the property under test is that the report does not depend
//! on it: the default run, `LOSAT_X_BNUNGAPCOUNT=0`, `LOSAT_X_DEFAULT=0` and the shadow run (which
//! runs both loops on every call and stops on the first difference) give the same bytes, and
//! the default run prints nothing on standard error (the `[X_SHADOW] BNUNGAPCOUNT` line is
//! printed only with `LOSAT_X_BNUNGAPCOUNTSHADOW`). Two searches reach the extension: blastn-short
//! (word size 7: every ungapped extension is the exact one) of short reads cut from LC738874,
//! and blastn (word size 11: the exact pass after the approximate one reaches the reduced cutoff)
//! of the first 9,000 letters of LC738874, both against LC738868. That those bytes equal NCBI's
//! is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNUNGAPCOUNT", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BNUNGAPCOUNTSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-ungapcount-{test}-{}", std::process::id()));
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

/// The two queries: 300 reads of 18 to 30 letters every 900 letters of LC738874, every third one
/// with its middle letter changed (so that extensions also stop at a mismatch), for
/// blastn-short; the first 9,000 letters of LC738874 for blastn.
fn queries(dir: &TestDir) -> [(&'static str, PathBuf); 2] {
    let genome = letters(&fasta("LC738874.fasta"));
    let mut reads = String::new();
    for i in 0..300 {
        let start = i * 900;
        let len = 18 + i % 13;
        let mut read: Vec<u8> = genome.as_bytes()[start..start + len].to_vec();
        if i % 3 == 0 {
            let mid = len / 2;
            read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(read).unwrap()));
    }
    let short = dir.0.join("reads.fna");
    std::fs::write(&short, reads).expect("write the reads");
    let long = dir.0.join("q9000.fna");
    std::fs::write(&long, format!(">q9000\n{}\n", &genome[..9_000])).expect("write the query");
    [("blastn-short", short), ("blastn", long)]
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

// The default search runs the counted loops; `LOSAT_X_BNUNGAPCOUNT=0` and `LOSAT_X_DEFAULT=0` run
// the ported loops. The reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_counted_loops() {
    let dir = TestDir::new("report");
    for (task, query) in queries(&dir) {
        for threads in [1, 4] {
            let (default, default_stderr) = search(task, &query, threads, &[]);
            assert!(!default.is_empty(), "{task} threads {threads}: no hits");
            assert_eq!(default_stderr, "", "{task} threads {threads}");
            for off in [SWITCH_OFF, DEFAULT_OFF] {
                let (report, stderr) = search(task, &query, threads, &[off]);
                assert_eq!(report, default, "{task} {off:?} threads {threads}");
                assert_eq!(stderr, "", "{task} {off:?} threads {threads}");
            }
        }
    }
}

// `LOSAT_X_DEFAULT=0` with the variable set to a value other than 0 runs the counted loops. The
// shadow (set to any value, not in the set; it implies the switch) compares the two loops on
// every call, prints its summary line with a count above 0, and does not change the report.
#[test]
fn the_shadow_compares_the_two_loops_and_the_report_is_unchanged() {
    let dir = TestDir::new("shadow");
    for (task, query) in queries(&dir) {
        let (expected, _) = search(task, &query, 1, &[]);
        let (report, _) = search(task, &query, 1, &[DEFAULT_OFF, ("LOSAT_X_BNUNGAPCOUNT", "1")]);
        assert_eq!(report, expected, "{task}");
        for switches in [&[SHADOW][..], &[DEFAULT_OFF, SHADOW][..]] {
            let (report, stderr) = search(task, &query, 1, switches);
            assert_eq!(report, expected, "{task} {switches:?}");
            let line = stderr
                .lines()
                .find(|line| line.starts_with("[X_SHADOW] BNUNGAPCOUNT compared="))
                .unwrap_or_else(|| panic!("{task} {switches:?}: no summary line:\n{stderr}"));
            assert!(
                !line.starts_with("[X_SHADOW] BNUNGAPCOUNT compared=0 "),
                "{task} {switches:?}: {line}"
            );
        }
    }
}
