//! dc-megablast's strips stage with `LOSAT_X_DCLATEST` (`src/algorithm/blastn/blast_engine/x_dc_latest.rs`,
//! `x_pair_seed_disc.rs` `seed_stage_latest`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes where the two-hit fold keeps its diagonals' inserts (a map of each diagonal's latest
//! insert instead of NCBI's diagonal hash, per partition of `LOSAT_X_PAIRPARFOLD` or the owner's
//! without partitions), not the decisions any seed makes, so the property under test is that the
//! report does not depend on it, and that the default took the switch's path (`LOSAT_X_REACH`).
//! Against `=0` at one and four threads: 9,800 letters of LC738874 against LC738868 (a query
//! block above the 8,000 letters of the diagonal hash) at the default window (40), above the word
//! length (24), and below it (12, and 4 in the band 2w + 2 < wl where NCBI's chains and LOSAT's
//! differ: below the word length no seed is extended, so no hit), and with the subject's lower case
//! masked
//! (`-lcase_masking`, the stage in the subject's ranges).
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_DCLATEST";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DCLATEST", "0");

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
            std::env::temp_dir().join(format!("losat-dclatest-{test}-{}", std::process::id()));
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

/// The letters of the first record of `name`.
fn letters(name: &str) -> String {
    let text = std::fs::read_to_string(fasta(name)).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// Letters 20,000 to 29,800 of LC738874 as one record in `dir`.
fn lc874_9800(dir: &TestDir) -> PathBuf {
    let query = dir.path("lc874_9800.fna");
    let text = letters("LC738874.fasta");
    std::fs::write(&query, format!(">lc874_9800\n{}\n", &text[20_000..29_800]))
        .expect("write the query");
    query
}

/// LC738868 with every third piece of 2,000 letters in lower case.
fn lc868_lower(dir: &TestDir) -> PathBuf {
    let subject = dir.path("lc868_lower.fna");
    let text: String = letters("LC738868.fasta")
        .as_bytes()
        .chunks(2_000)
        .enumerate()
        .map(|(k, piece)| {
            let piece = String::from_utf8_lossy(piece).into_owned();
            if k % 3 == 1 {
                piece.to_lowercase()
            } else {
                piece
            }
        })
        .collect();
    std::fs::write(&subject, format!(">lc868_lower\n{text}\n")).expect("write the subject");
    subject
}

/// One dc-megablast search of `query` against `subject` (`-outfmt 6` on standard output) with
/// `args` and the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    subject: &Path,
    args: &[&str],
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
        .args(["-task", "dc-megablast"])
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(args)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{args:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For `query` against `subject` with `args`, at one and four threads: the report with the switch
/// off and by default, the same bytes and standard error; the default took the switch's path and
/// never the reference's (in the stage without a pool at one thread, `LOSAT_X_BNDCONE`, and in a
/// partitioned stage at four, `LOSAT_X_PAIRPARFOLD`), and `=0` records nothing for the switch.
/// `hits`: whether the report has hits.
fn check(query: &Path, subject: &Path, args: &[&str], hits: bool) {
    for threads in [1, 4] {
        let at = format!("{args:?} threads {threads}");
        let off_reach = reach_support::Reach::new("dclatest-off");
        let (expected, expected_stderr) = search(
            query,
            subject,
            args,
            threads,
            &[SWITCH_OFF, off_reach.var()],
        );
        assert_eq!(!expected.is_empty(), hits, "{at}: hits");
        assert_eq!(off_reach.sides(SWITCH), (false, false), "{at}: =0");
        let reach = reach_support::Reach::new("dclatest");
        let (report, stderr) = search(query, subject, args, threads, &[reach.var()]);
        assert!(report == expected, "{at}: default differs");
        assert_eq!(stderr, expected_stderr, "{at}");
        assert_eq!(
            reach.sides(SWITCH),
            (true, false),
            "{at}: sides of {SWITCH}"
        );
        assert!(
            reach.sides("LOSAT_X_PAIRPARDISC").0,
            "{at}: no strips stage"
        );
        if threads == 1 {
            assert!(
                reach.sides("LOSAT_X_BNDCONE").0,
                "{at}: no stage without a pool"
            );
        } else {
            assert!(
                reach.sides("LOSAT_X_PAIRPARFOLD").0,
                "{at}: no partitioned stage"
            );
        }
    }
}

// The default window and one above the word length (hits), one below it and one in the band (no
// hit).
#[test]
fn dc_megablast_windows() {
    let dir = TestDir::new("windows");
    let query = lc874_9800(&dir);
    let subject = fasta("LC738868.fasta");
    check(&query, &subject, &[], true);
    check(&query, &subject, &["-window_size", "24"], true);
    check(&query, &subject, &["-window_size", "12"], false);
    check(&query, &subject, &["-window_size", "4"], false);
}

// A subject with lower case masked: the stage scans its ranges.
#[test]
fn dc_megablast_masked_subject() {
    let dir = TestDir::new("masked");
    let query = lc874_9800(&dir);
    let subject = lc868_lower(&dir);
    check(&query, &subject, &["-lcase_masking"], true);
}
