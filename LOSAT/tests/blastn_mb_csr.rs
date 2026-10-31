//! Megablast's tables of 4^12 cells laid out as packed CSR (`LOSAT_X_MBCSR`,
//! `src/algorithm/blastn/x_bn_csr.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It chooses between two LOSAT layouts of the same chains, so the property under test is that the
//! report does not depend on it, and that its size rule picks the side it says: `LOSAT_X_REACH`
//! shows the switch's sides. 200 kb of Sakai (both strands 400,000 positions, a 4^12-cell table)
//! against five copies of MG1655 (23 Mb) is above the rule of 8 letters per position and one per
//! cell (19.9 M); against one copy (4.6 Mb) it is below. Each against a subject file (`-subject`)
//! and a database (`-db`), at one and four threads.
//! That those bytes equal NCBI's is the business of the fixture suites.
//! `-task blastn` of 50 kb against 2 Mb has a packed table of 4^11 cells, whose scan tests each
//! word once (`LOSAT_X_BNPVONCE`); the report does not depend on that either.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_MBCSR";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_MBCSR", "0");
/// The query prefetch of the loop over the offset pairs (`LOSAT_X_MBHITPF`, strict set).
const HIT_PREFETCH_OFF: (&str, &str) = ("LOSAT_X_MBHITPF", "0");

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
    search_task("megablast", query, subjects, threads, vars)
}

fn search_task(
    task: &str,
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
        .args(["-task", task, "-outfmt", "6"])
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

fn check(dir: &TestDir, copies: usize, taken: bool) {
    let query = dir.path("sakai_200k.fna");
    let sakai = letters(&fasta("Sakai.fna"));
    std::fs::write(
        &query,
        format!(">sakai_200k\n{}\n", &sakai[1_000_000..1_200_000]),
    )
    .expect("write the query");
    let query = query.as_path();
    let subject = dir.path("mg1655_copies.fna");
    let genome = letters(&fasta("MG1655.fna"));
    let records: String = (0..copies).map(|i| format!(">mg{i}\n{genome}\n")).collect();
    std::fs::write(&subject, records).expect("write the subjects");
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
            let (report, _) = search(query, subjects, threads, &[HIT_PREFETCH_OFF]);
            assert!(
                report == expected,
                "{what} threads {threads}: hit prefetch off differs"
            );
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
    check(&TestDir::new("above"), 5, true);
}

#[test]
fn below_the_size_rule_keeps_chains() {
    check(&TestDir::new("below"), 1, false);
}

#[test]
fn a_word_is_tested_once_on_a_packed_table_of_4_11_cells() {
    let dir = TestDir::new("pvonce");
    let query = dir.path("sakai_50k.fna");
    let sakai = letters(&fasta("Sakai.fna"));
    std::fs::write(
        &query,
        format!(">sakai_50k\n{}\n", &sakai[1_000_000..1_050_000]),
    )
    .expect("write the query");
    let subject = dir.path("mg1655_2m.fna");
    let genome = letters(&fasta("MG1655.fna"));
    std::fs::write(&subject, format!(">mg\n{}\n", &genome[..2_000_000]))
        .expect("write the subject");
    let subjects = ("-subject", subject.as_path());
    let (expected, _) = search_task("blastn", &query, subjects, 1, &[("LOSAT_X_BNPVONCE", "0")]);
    assert!(!expected.is_empty(), "no hits");
    for threads in [1, 4] {
        let reach = reach_support::Reach::new("pvonce");
        let (report, _) = search_task("blastn", &query, subjects, threads, &[reach.var()]);
        assert!(report == expected, "threads {threads}: default differs");
        assert_eq!(
            reach.sides("LOSAT_X_BNPVONCE"),
            (true, false),
            "threads {threads}: the side"
        );
    }
}
