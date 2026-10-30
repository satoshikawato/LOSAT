//! Megablast's tables of 4^12 cells laid out as packed CSR (`LOSAT_X_MBCSR`,
//! `src/algorithm/blastn/x_bn_csr.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It chooses between two LOSAT layouts of the same chains, so the property under test is that the
//! report does not depend on it, and that its size rule picks the side it says: `LOSAT_X_REACH`
//! shows the switch's sides. 200 kb of Sakai (both strands 400,000 positions, a 4^12-cell table)
//! against MG1655 (4.6 Mb) is above the rule of 8 letters per position; all of Sakai (11 M
//! positions) against the same subject is below it. Each against a subject file (`-subject`) and
//! a database (`-db`), at one and four threads.
//! That those bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_MBCSR";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_MBCSR", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-mbcsr-{test}-{}", std::process::id()));
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

fn database(dir: &TestDir, subject: &Path) -> PathBuf {
    let db = dir.path("subject_db");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(subject)
        .args(["-dbtype", "nucl", "-out"])
        .arg(&db)
        .output()
        .expect("run makeblastdb");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

fn search(
    query: &Path,
    subjects: (&str, &Path),
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
        .args(["-task", "megablast", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{} {vars:?} threads {threads}: {stderr}",
        subjects.0
    );
    (output.stdout, stderr)
}

fn check(dir: &TestDir, query: &Path, taken: bool) {
    let subject = fasta("MG1655.fna");
    let db = database(dir, &subject);
    for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
        let what = format!("{} taken {taken}", subjects.0);
        let (expected, expected_stderr) = search(query, subjects, 1, &[SWITCH_OFF]);
        assert!(!expected.is_empty(), "{what}: no hits");
        for threads in [1, 4] {
            let reach = reach_support::Reach::new("mbcsr");
            let (report, stderr) = search(query, subjects, threads, &[reach.var()]);
            assert!(
                report == expected,
                "{what} threads {threads}: default differs"
            );
            assert_eq!(stderr, expected_stderr, "{what} threads {threads}");
            assert_eq!(
                reach.sides(SWITCH),
                (taken, !taken),
                "{what} threads {threads}: the size rule's side"
            );
        }
    }
}

#[test]
fn above_the_size_rule_is_csr() {
    let dir = TestDir::new("above");
    let query = dir.path("sakai_200k.fna");
    let genome = letters(&fasta("Sakai.fna"));
    std::fs::write(
        &query,
        format!(">sakai_200k\n{}\n", &genome[1_000_000..1_200_000]),
    )
    .expect("write the query");
    check(&dir, &query, true);
}

#[test]
fn below_the_size_rule_keeps_chains() {
    let dir = TestDir::new("below");
    check(&dir, &fasta("Sakai.fna"), false);
}
