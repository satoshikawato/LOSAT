//! The partitioned seed stages of blastn's strips paths in seed pieces with `LOSAT_X_SEEDPIECE`
//! (`src/algorithm/blastn/blast_engine/x_pair_seed_piece.rs` `seed_stage_pieces`, reached from
//! `x_pair_seed_part.rs` `seed_stage_part_lists`: dc-megablast's two-hit stage of
//! `x_pair_seed_disc.rs` and the stages of `x_pair_seed.rs` and `x_pair_seed_na.rs`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`) and carries a number: it is on unless
//! its variable is `0`, with pieces of the number of seeds given (2 or more), else of the built-in
//! 16,384. It changes where a strip's seeds wait for the partitions' folds (pieces of a fixed size
//! whose scan stops at the first pair a piece has no room for and goes on from it, instead of one
//! list per partition and strip), not the seeds, their order or their keys, so the property under
//! test is that the report does not depend on it, and that the default took the switch's path
//! (`LOSAT_X_REACH`). Against a subject file (`-subject`): the report with `=0` at four threads, and
//! by default at two and four threads and with pieces of a few seeds at two and three threads (so
//! that the scans stop inside the strips, inside words' runs, and go on on other threads), the
//! same bytes and standard error. dc-megablast of Sakai against MG1655 (a table laid out as CSR)
//! and of 9,800 letters of LC738874 against LC738868 (a table of chains); `-task blastn` and
//! `megablast` of LC738874 against LC738868 (the two-stage table, `LOSAT_X_PAIRPAR`), and
//! `-task blastn-short` of 700 24-mers cut from LC738868 against it (the direct table,
//! `LOSAT_X_PAIRPARNA`).
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_SEEDPIECE";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-piece-{test}-{}", std::process::id()));
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

/// The letters of the first record of `name`, upper case.
fn letters(name: &str) -> Vec<u8> {
    let text = std::fs::read_to_string(fasta(name)).expect("read the record");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// Letters 20,000 to 29,800 of LC738874 as one record in `dir` (both strands above the 8,000
/// letters of the diagonal hash, so the stage is partitioned).
fn lc874_9800(dir: &TestDir) -> PathBuf {
    let seq = letters("LC738874.fasta");
    let query = dir.path("lc874_9800.fna");
    std::fs::write(
        &query,
        format!(
            ">lc874_9800\n{}\n",
            String::from_utf8_lossy(&seq[20_000..29_800])
        ),
    )
    .expect("write the query");
    query
}

/// 700 24-mers spread over LC738868, every third reverse-complemented (16,800 letters: a query
/// block above the 8,000 letters of the diagonal hash, so the stage is partitioned).
fn short_queries(dir: &TestDir) -> PathBuf {
    let seq = letters("LC738868.fasta");
    let step = (seq.len() - 24) / 700;
    let mut out = String::new();
    for i in 0..700 {
        let mut word = seq[i * step..i * step + 24].to_vec();
        if i % 3 == 1 {
            word.reverse();
            for b in &mut word {
                *b = match *b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    other => other,
                };
            }
        }
        out.push_str(&format!(">q{i}\n{}\n", String::from_utf8(word).unwrap()));
    }
    let query = dir.path("short.fna");
    std::fs::write(&query, out).expect("write the queries");
    query
}

/// One search of `query` against `subject` (`-outfmt 6` on standard output) with `task_args` and
/// the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    subject: &Path,
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
        .arg(subject)
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

/// For `query` against `subject` with `task_args`: the report with the switch off at four
/// threads, then by default at two and four threads and with pieces of `small` seeds at two and
/// three threads, the same bytes and standard error; those runs took the switch's path and never
/// the reference's, in a partitioned stage (`LOSAT_X_PAIRPARFOLD`), and reached the switch's side
/// of each of `stages` (the stage, and the table's layout where the scan of a case depends on it).
fn check(query: &Path, subject: &Path, task_args: &[&str], stages: &[&str], small: &str) {
    let (expected, expected_stderr) = search(query, subject, task_args, 4, &[(SWITCH, "0")]);
    assert!(!expected.is_empty(), "{task_args:?}: no hits");
    for (threads, value) in [(2, None), (4, None), (2, Some(small)), (3, Some(small))] {
        let reach = reach_support::Reach::new("piece");
        let mut vars = vec![reach.var()];
        vars.extend(value.map(|v| (SWITCH, v)));
        let (report, stderr) = search(query, subject, task_args, threads, &vars);
        let at = format!("{task_args:?} threads {threads} {value:?}");
        assert!(report == expected, "{at}: differs from =0");
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
        for stage in stages {
            assert!(reach.sides(stage).0, "{at}: no stage of {stage}");
        }
    }
}

// dc-megablast with a table laid out as CSR (`LOSAT_X_DCCSR`); pieces of 97 seeds.
#[test]
fn dc_megablast_csr() {
    check(
        &fasta("Sakai.fna"),
        &fasta("MG1655.fna"),
        &["-task", "dc-megablast"],
        &["LOSAT_X_PAIRPARDISC", "LOSAT_X_DCCSR"],
        "97",
    );
}

// dc-megablast with a table of chains; pieces of 5 seeds.
#[test]
fn dc_megablast_chains() {
    let dir = TestDir::new("chains");
    check(
        &lc874_9800(&dir),
        &fasta("LC738868.fasta"),
        &["-task", "dc-megablast"],
        &["LOSAT_X_PAIRPARDISC"],
        "5",
    );
}

// The two-stage table's stage (`x_pair_seed.rs`), word size 11; pieces of 5 seeds.
#[test]
fn blastn_two_stage_table() {
    check(
        &fasta("LC738874.fasta"),
        &fasta("LC738868.fasta"),
        &["-task", "blastn"],
        &["LOSAT_X_PAIRPAR", "LOSAT_X_BNCSR"],
        "5",
    );
}

// The two-stage table's stage, word size 28; pieces of 3 seeds.
#[test]
fn megablast_two_stage_table() {
    check(
        &fasta("LC738874.fasta"),
        &fasta("LC738868.fasta"),
        &["-task", "megablast"],
        &["LOSAT_X_PAIRPAR"],
        "3",
    );
}

// The direct table's stage (`x_pair_seed_na.rs`), word size 7; pieces of 5 seeds.
#[test]
fn blastn_short_direct_table() {
    let dir = TestDir::new("short");
    check(
        &short_queries(&dir),
        &fasta("LC738868.fasta"),
        &["-task", "blastn-short"],
        &["LOSAT_X_PAIRPARNA"],
        "5",
    );
}
