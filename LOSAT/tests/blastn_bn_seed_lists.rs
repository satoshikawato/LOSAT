//! The partitioned seed stages of blastn's tables with `LOSAT_X_BNSEEDLISTS`
//! (`src/algorithm/blastn/blast_engine/x_pair_seed.rs` and `x_pair_seed_na.rs` `seed_stage_of`,
//! `src/algorithm/blastn/blast_engine/x_pair_seed_part.rs` `seed_stage_part_lists`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes how the partitioned stage of the strips hands each partition its seeds (the scan
//! into the partitions' lists) and merges them back (every bucket, the folds recording none), not
//! the seeds or their order, so the property under test is that the report does not depend on
//! it, and that the default took the switch's path (`LOSAT_X_REACH`).
//! Against a subject file (`-subject`) at two and four threads, against `=0` at both (and at four
//! with `LOSAT_X_PAIRPARHOME=0`, where the stage cuts the strips as before but folds without
//! recording its buckets): `-task blastn` and `megablast` of LC738874 against LC738868 (the
//! two-stage table, `LOSAT_X_PAIRPAR`), and `-task blastn-short` of 700 24-mers cut from
//! LC738868 against it (16,800 letters, the direct table, `LOSAT_X_PAIRPARNA`).
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_BNSEEDLISTS";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNSEEDLISTS", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-bnlists-{test}-{}", std::process::id()));
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

/// 700 24-mers of LC738868 at evenly spaced offsets in `dir`, every third one
/// reverse-complemented (16,800 letters: a query block above the 8,000 letters of the diagonal
/// hash, so the stage is partitioned).
fn short_queries(dir: &TestDir) -> PathBuf {
    let text = std::fs::read_to_string(fasta("LC738868.fasta")).expect("read LC738868");
    let letters: Vec<u8> = text
        .lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let step = (letters.len() - 24) / 700;
    let mut out = String::new();
    for i in 0..700 {
        let mut seq = letters[i * step..i * step + 24].to_vec();
        if i % 3 == 1 {
            seq.reverse();
            for b in &mut seq {
                *b = match *b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    other => other,
                };
            }
        }
        out.push_str(&format!(">q{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let query = dir.path("short.fna");
    std::fs::write(&query, out).expect("write the queries");
    query
}

/// One search of `query` against LC738868 (`-outfmt 6` on standard output) with `task_args`
/// and the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    task_args: &[&str],
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
        .args(task_args)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-evalue", "0.001", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task_args:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For `query` with `task_args`: the report with the switch off at two and four threads, then by
/// default at two and four threads and by default with `LOSAT_X_PAIRPARHOME=0` at four, the same
/// bytes and standard error; the default took the switch's path and never the reference's, in a
/// partitioned stage (`LOSAT_X_PAIRPARFOLD`) of `stage`.
fn check(query: &Path, task_args: &[&str], stage: &str) {
    let (expected, expected_stderr) = search(query, task_args, 4, &[SWITCH_OFF]);
    assert!(!expected.is_empty(), "{task_args:?}: no hits");
    let (off2, off2_stderr) = search(query, task_args, 2, &[SWITCH_OFF]);
    assert!(off2 == expected, "{task_args:?}: =0 at 2 threads differs");
    assert_eq!(
        off2_stderr, expected_stderr,
        "{task_args:?}: =0 at 2 threads"
    );
    for (threads, home) in [
        (2, None),
        (4, None),
        (4, Some(("LOSAT_X_PAIRPARHOME", "0"))),
    ] {
        let reach = reach_support::Reach::new("bnlists");
        let mut vars = vec![reach.var()];
        vars.extend(home);
        let (report, stderr) = search(query, task_args, threads, &vars);
        let at = format!("{task_args:?} threads {threads} {home:?}");
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
        assert!(reach.sides(stage).0, "{at}: no stage of {stage}");
    }
}

// The two-stage table's stage (`x_pair_seed.rs`), word size 11.
#[test]
fn blastn_two_stage_table() {
    check(
        &fasta("LC738874.fasta"),
        &["-task", "blastn"],
        "LOSAT_X_PAIRPAR",
    );
}

// The two-stage table's stage, word size 28.
#[test]
fn megablast_two_stage_table() {
    check(
        &fasta("LC738874.fasta"),
        &["-task", "megablast"],
        "LOSAT_X_PAIRPAR",
    );
}

// The direct table's stage (`x_pair_seed_na.rs`), word size 7.
#[test]
fn blastn_short_direct_table() {
    let dir = TestDir::new("short");
    check(
        &short_queries(&dir),
        &["-task", "blastn-short"],
        "LOSAT_X_PAIRPARNA",
    );
}
