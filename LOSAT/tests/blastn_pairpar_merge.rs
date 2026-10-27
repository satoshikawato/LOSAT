//! The merge of a partitioned seed stage putting the stage's hits in key order strip by strip, on
//! the pool in ranges of strips (`LOSAT_X_PAIRPARMERGE`,
//! `src/algorithm/blastn/blast_engine/x_pair_seed_part.rs`), on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes how the hits of a stage are put in the order of their keys, not the order, so the
//! property under test is that the report does not depend on it, and that the default took the
//! switch's path in every stage (`LOSAT_X_REACH`), whatever its number of hits.
//! The three folds that end in the merge, each against a subject file (`-subject`) at two and four
//! threads (and at four with `LOSAT_X_PAIRPARHOME=0`): the one-hit fold of `-task blastn` (LC738874 against LC738868, 17,506 hits in its one
//! stage; and 9,800 letters of LC738874, 412 hits), the one-hit fold of the
//! Na table of `-task blastn-short` (the same 9,800 letters, 42,126 hits) and the two-hit fold of
//! `-task dc-megablast` (Sakai against MG1655, 16,720 hits).
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_PAIRPARMERGE";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_PAIRPARMERGE", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-pmerge-{test}-{}", std::process::id()));
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

/// Letters 20,000 to 29,800 of LC738874 as one record in `dir` (both strands above the 8,000
/// letters of the diagonal hash, so the stage is partitioned).
fn lc874_9800(dir: &TestDir) -> PathBuf {
    let text = std::fs::read_to_string(fasta("LC738874.fasta")).expect("read LC738874");
    let letters: String = text
        .lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect();
    let query = dir.path("lc874_9800.fna");
    std::fs::write(
        &query,
        format!(">lc874_9800\n{}\n", &letters[20_000..29_800]),
    )
    .expect("write the query");
    query
}

/// One search of `query` against `subject` (`-outfmt 6` on standard output) with `extra`
/// arguments and the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    subject: &Path,
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
        .arg("-subject")
        .arg(subject)
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For `query` and `extra` against `subject`: the report with the switch off at four threads,
/// then by default at two and four threads and by default with `LOSAT_X_PAIRPARHOME=0` (the other
/// participant loop before the merge) at four, the same bytes and standard error; the default took
/// the switch's path in its stages and never the sort.
fn check(query: &Path, subject: &Path, extra: &[&str]) {
    let (expected, expected_stderr) = search(query, subject, extra, 4, &[SWITCH_OFF]);
    assert!(!expected.is_empty(), "{extra:?}: no hits");
    for (threads, home) in [
        (2, None),
        (4, None),
        (4, Some(("LOSAT_X_PAIRPARHOME", "0"))),
    ] {
        let reach = reach_support::Reach::new("pmerge");
        let mut vars = vec![reach.var()];
        vars.extend(home);
        let (report, stderr) = search(query, subject, extra, threads, &vars);
        assert!(
            report == expected,
            "{extra:?} threads {threads} {home:?}: default differs"
        );
        assert_eq!(
            stderr, expected_stderr,
            "{extra:?} threads {threads} {home:?}"
        );
        assert_eq!(
            reach.sides(SWITCH),
            (true, false),
            "{extra:?} threads {threads} {home:?}: sides of {SWITCH}"
        );
    }
}

// The one-hit fold of `-task blastn` with 17,506 hits in its stage, and with 412.
#[test]
fn blastn_one_hit_stages() {
    let dir = TestDir::new("bn");
    let subject = fasta("LC738868.fasta");
    check(&fasta("LC738874.fasta"), &subject, &["-task", "blastn"]);
    check(&lc874_9800(&dir), &subject, &["-task", "blastn"]);
}

// The one-hit fold of the Na table of `-task blastn-short` (PAIRPARNA), 42,126 hits.
#[test]
fn blastn_short_na_stage() {
    let dir = TestDir::new("bns");
    check(
        &lc874_9800(&dir),
        &fasta("LC738868.fasta"),
        &["-task", "blastn-short"],
    );
}

// The two-hit fold of `-task dc-megablast` (PAIRPARDISC), 16,720 hits.
#[test]
fn dc_megablast_two_hit_stage() {
    check(
        &fasta("Sakai.fna"),
        &fasta("MG1655.fna"),
        &["-task", "dc-megablast"],
    );
}
