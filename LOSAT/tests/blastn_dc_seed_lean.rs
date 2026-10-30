//! dc-megablast's seed stage with `LOSAT_X_DCSEEDLEAN` (`src/algorithm/blastn/x_disc_lookup_par.rs`
//! `csr_on_pool_steps`, `src/algorithm/blastn/blast_engine/x_pair_seed_disc.rs` `seed_stage_lists`,
//! `src/algorithm/blastn/blast_engine/x_pair_seed_part.rs` `seed_stage_part_lists`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes how a one-template table's CSR is built on the pool (in steps over two sets of
//! buckets) and how the partitioned stage of the strips hands each partition its seeds (the scan
//! into the partitions' lists) and merges them back (every bucket, the folds recording none), not
//! the table, the seeds or their order, so the property under test is that the report does not
//! depend on it, and that the default took the switch's path (`LOSAT_X_REACH`).
//! Against a subject file (`-subject`) at two and four threads, against `=0` at both (and at four with
//! `LOSAT_X_PAIRPARHOME=0`, where the stage cuts the strips as before but folds without recording
//! its buckets): Sakai against MG1655 (a table laid out as CSR, built on the pool, and 16,720 hits
//! in the partitioned stage), and 9,800 letters of LC738874 against LC738868 (a table of chains,
//! below the CSR's size rule, and a partitioned stage).
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_DCSEEDLEAN";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DCSEEDLEAN", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-dclean-{test}-{}", std::process::id()));
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

/// One dc-megablast search of `query` against `subject` (`-outfmt 6` on standard output) with
/// the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(query: &Path, subject: &Path, threads: u32, vars: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For `query` against `subject`: the report with the switch off at two and four threads, then by
/// default at two and four threads and by default with `LOSAT_X_PAIRPARHOME=0` at four, the same
/// bytes and standard error; the default took the switch's path and never the reference's, in a
/// partitioned stage (`LOSAT_X_PAIRPARFOLD`), with the table built on the pool
/// (`LOSAT_X_DISCLUTPAR`) and laid out as CSR (`LOSAT_X_DCCSR`) where `csr`, as chains where not.
fn check(query: &Path, subject: &Path, csr: bool) {
    let (expected, expected_stderr) = search(query, subject, 4, &[SWITCH_OFF]);
    assert!(!expected.is_empty(), "{query:?}: no hits");
    let (off2, off2_stderr) = search(query, subject, 2, &[SWITCH_OFF]);
    assert!(off2 == expected, "{query:?}: =0 at 2 threads differs");
    assert_eq!(off2_stderr, expected_stderr, "{query:?}: =0 at 2 threads");
    for (threads, home) in [
        (2, None),
        (4, None),
        (4, Some(("LOSAT_X_PAIRPARHOME", "0"))),
    ] {
        let reach = reach_support::Reach::new("dclean");
        let mut vars = vec![reach.var()];
        vars.extend(home);
        let (report, stderr) = search(query, subject, threads, &vars);
        let at = format!("{query:?} threads {threads} {home:?}");
        assert!(report == expected, "{at}: default differs");
        assert_eq!(stderr, expected_stderr, "{at}");
        assert_eq!(
            reach.sides(SWITCH),
            (true, false),
            "{at}: sides of {SWITCH}"
        );
        assert!(
            reach.sides("LOSAT_X_PAIRPARFOLD").0,
            "{at}: no partitioned stage"
        );
        assert!(
            reach.sides("LOSAT_X_DISCLUTPAR").0,
            "{at}: no build on the pool"
        );
        if csr {
            assert!(reach.sides("LOSAT_X_DCCSR").0, "{at}: no CSR");
        } else {
            assert_eq!(reach.sides("LOSAT_X_DCCSR"), (false, true), "{at}: CSR");
        }
    }
}

// A table laid out as CSR and built on the pool, and a partitioned stage of 16,720 hits.
#[test]
fn dc_megablast_csr_on_the_pool() {
    check(&fasta("Sakai.fna"), &fasta("MG1655.fna"), true);
}

// A table of chains (below the CSR's size rule) and a partitioned stage.
#[test]
fn dc_megablast_chains() {
    let dir = TestDir::new("chains");
    check(&lc874_9800(&dir), &fasta("LC738868.fasta"), false);
}
