//! A subject's HSPs grouped by query (`LOSAT_X_BNHITGROUP`, `src/algorithm/blastn/blast_engine/run.rs`
//! `PrelimHitCollector::x_add_grouped`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:116-161
//! ```c
//!       for (index = 0; index < hsp_list->hspcnt; index++) {
//!          Int4 query_index;
//!          hsp = hsp_list->hsp_array[index];
//!          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
//! ```
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes how the lists of one subject are made, not what they hold (the lib test
//! `grouped_hit_lists_equal_the_lists_per_query`), so the property under test is that the report
//! does not depend on it: the default run and `LOSAT_X_BNHITGROUP=0` give the same bytes and print
//! nothing on standard error, with subjects read from FASTA (`-subject`) and from a database
//! (`-db`), at 1 and 4 threads, for blastn-short and megablast, with several queries per subject.
//! That those bytes equal NCBI's is the business of the fixture suites.
//!
//! The queries are 300 reads of 40 to 52 letters every 27 letters of the first 8,100 letters of
//! LC738874, every third one with its middle letter changed; the subjects 150 windows of 4,000
//! letters every 30 letters of the same genome.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNHITGROUP", "0");

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
            std::env::temp_dir().join(format!("losat-hitgroup-{test}-{}", std::process::id()));
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

/// The binary with no `LOSAT_*` variable that chooses a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`,
/// as in CI) needs that file name as soon as a search has more than one thread, and it selects
/// nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
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

/// The reads, the subject windows as FASTA, and their database made by `losat makeblastdb`.
fn inputs(dir: &TestDir) -> (PathBuf, PathBuf, String) {
    let genome = letters(&fasta("LC738874.fasta"));
    let mut reads = String::new();
    for i in 0..300 {
        let start = i * 27;
        let len = 40 + i % 13;
        let mut read: Vec<u8> = genome.as_bytes()[start..start + len].to_vec();
        if i % 3 == 0 {
            let mid = len / 2;
            read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(read).unwrap()));
    }
    let query = dir.0.join("reads.fna");
    std::fs::write(&query, reads).expect("write the reads");
    let mut windows = String::new();
    for j in 0..150 {
        let start = j * 30;
        windows.push_str(&format!(">w{j}\n{}\n", &genome[start..start + 4_000]));
    }
    let set = dir.0.join("windows.fna");
    std::fs::write(&set, windows).expect("write the windows");
    let db = dir.0.join("windows").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, set, db)
}

/// One search (`-outfmt 6` on standard output) with the given variables as its only `LOSAT_*`
/// variables that choose a code path. Returns the report and the standard error.
fn search(
    query: &Path,
    target: &[&str],
    task: &str,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .args(target)
        .args(["-task", task, "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{target:?} {task} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The largest number of queries one subject has rows for.
fn most_queries(report: &[u8]) -> usize {
    let mut per_subject = std::collections::HashMap::<&str, std::collections::HashSet<&str>>::new();
    for line in std::str::from_utf8(report).unwrap().lines() {
        let mut fields = line.split('\t');
        let (query, subject) = (fields.next().unwrap(), fields.next().unwrap());
        per_subject.entry(subject).or_default().insert(query);
    }
    per_subject.values().map(|q| q.len()).max().unwrap_or(0)
}

// The default search groups a subject's HSPs by query; `LOSAT_X_BNHITGROUP=0` makes a list per
// query of the batch. The reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_how_a_subjects_hsps_are_grouped() {
    let dir = TestDir::new("report");
    let (query, set, db) = inputs(&dir);
    let set = set.to_str().unwrap();
    for target in [["-subject", set], ["-db", db.as_str()]] {
        for task in ["blastn-short", "megablast"] {
            for threads in [1, 4] {
                let what = format!("{target:?} {task} threads {threads}");
                let (default, default_stderr) = search(&query, &target, task, threads, &[]);
                assert_eq!(default_stderr, "", "{what}");
                let most = most_queries(&default);
                assert!(most > 10, "{what}: at most {most} queries per subject");
                let (old, old_stderr) = search(&query, &target, task, threads, &[SWITCH_OFF]);
                assert_eq!(old_stderr, "", "{what} =0");
                assert!(old == default, "{what}: the reports differ");
            }
        }
    }
}
