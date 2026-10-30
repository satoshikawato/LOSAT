//! The interval trees of the preliminary stage in one pool shared by the threads
//! (`LOSAT_X_BPTREEPOOL`, `src/algorithm/blastp/blast_engine/x_tree_pool.rs`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3826-3827
//! ```c
//!       tree = Blast_IntervalTreeInit(0, query->length+1,
//!                                     0, subject->length+1);
//! ```
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`. It changes only where a subject's tree lives, never which HSPs it holds or what it answers
//! (the lib tests of `x_tree_pool`), so the property under test is that the report does not depend
//! on it: the default run and `LOSAT_X_BPTREEPOOL=0` give the same bytes and print nothing on
//! standard error, with subjects read from FASTA (`-subject`) and from a database (`-db`), at 1 and
//! 4 threads. That those bytes equal NCBI's is the business of the fixture suites.
//!
//! The queries are 450 windows of 50 letters every 13 letters of the fixture's longest protein
//! (BDT62816.1 of `MejoMJNV.faa`, 6,729 letters); the subjects are that proteome (124 proteins), so
//! the long subject's tree lists at least one HSP per window, and other subjects' trees follow it.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPTREEPOOL", "0");

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
            std::env::temp_dir().join(format!("losat-treepool-{test}-{}", std::process::id()));
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

/// The letters of record `id` of a FASTA file.
fn record(path: &Path, id: &str) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with(&format!(">{id} ")))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// The query windows, and the proteome's database made by `losat makeblastdb`.
fn inputs(dir: &TestDir) -> (PathBuf, String) {
    let protein = record(&fasta("MejoMJNV.faa"), "BDT62816.1");
    assert_eq!(protein.len(), 6_729);
    let mut windows = String::new();
    for i in 0..450 {
        let start = i * 13;
        windows.push_str(&format!(">q{i}\n{}\n", &protein[start..start + 50]));
    }
    let query = dir.0.join("windows.faa");
    std::fs::write(&query, windows).expect("write the windows");
    let db = dir.0.join("proteome").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(fasta("MejoMJNV.faa"))
        .args(["-dbtype", "prot", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, db)
}

/// One blastp search (`-outfmt 6` on standard output) with the given variables as its only
/// `LOSAT_*` variables that choose a code path. Returns the report and the standard error.
fn search(
    query: &Path,
    target: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .args(target)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{target:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

// The reports are equal whichever way the subjects' trees are held, and no run prints a line on
// standard error.
#[test]
fn the_report_does_not_depend_on_where_the_interval_trees_live() {
    let dir = TestDir::new("report");
    let (query, db) = inputs(&dir);
    let proteome = fasta("MejoMJNV.faa");
    let proteome = proteome.to_str().unwrap();
    for target in [["-subject", proteome], ["-db", db.as_str()]] {
        for threads in [1, 4] {
            let what = format!("{target:?} threads {threads}");
            let (default, default_stderr) = search(&query, &target, threads, &[]);
            assert_eq!(default_stderr, "", "{what}");
            assert!(
                default
                    .split(|&b| b == b'\n')
                    .filter(|l| !l.is_empty())
                    .count()
                    >= 450,
                "{what}: fewer rows than windows"
            );
            let (old, old_stderr) = search(&query, &target, threads, &[SWITCH_OFF]);
            assert_eq!(old_stderr, "", "{what} =0");
            assert!(old == default, "{what}: the reports differ");
        }
    }
}
