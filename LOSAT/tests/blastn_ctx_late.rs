//! A lookup hit's context and cutoffs read once the hit is ready for its ungapped extension
//! (`LOSAT_X_BNCTXLATE`, `src/algorithm/blastn/blast_engine/na_scan.rs` `scan_na_lookup` and the
//! fold of LOSAT_X_PAIRPARNA, `x_pair_seed_na.rs`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:900-903
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`. It changes only when the context's rows are read, never which hits are extended or kept,
//! so the property under test is that the report does not depend on it: the default run and
//! `LOSAT_X_BNCTXLATE=0` give the same bytes and print nothing on standard error, for the scans
//! of a table that is not two-stage: blastn-short's word 7 with a query block above 8,000
//! letters (the diagonal hash) at 1 thread and at 4 (the fold of LOSAT_X_PAIRPARNA: one subject,
//! one hit), with a block below it (the diagonal array), in two-hit mode (`-window_size`), and
//! `-task blastn -word_size 8`. That those bytes equal NCBI's is the business of the fixture
//! suites.
//!
//! The queries are reads of 18 to 30 letters every 27 letters of LC738874, every third one with
//! its middle letter changed; the subject is the whole genome (287,061 letters).

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNCTXLATE", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-ctxlate-{test}-{}", std::process::id()));
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

/// `count` reads of 18 to 30 letters every 27 letters of the genome, every third one with its
/// middle letter changed.
fn reads(dir: &TestDir, genome: &str, count: usize) -> PathBuf {
    let mut reads = String::new();
    for i in 0..count {
        let start = i * 27;
        let len = 18 + i % 13;
        let mut read: Vec<u8> = genome.as_bytes()[start..start + len].to_vec();
        if i % 3 == 0 {
            let mid = len / 2;
            read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(read).unwrap()));
    }
    let query = dir.0.join(format!("reads{count}.fna"));
    std::fs::write(&query, reads).expect("write the reads");
    query
}

/// One search of the genome (`-outfmt 6` on standard output) with the given variables as its
/// only `LOSAT_*` variables that choose a code path. Returns the report and the standard error.
fn search(
    query: &Path,
    task: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta("LC738874.fasta"))
        .args(task)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

// The default search reads a hit's context once the hit is ready; `LOSAT_X_BNCTXLATE=0` reads it
// for every lookup hit. The reports are equal and not empty, and no run prints a line on standard
// error.
#[test]
fn the_report_does_not_depend_on_when_the_context_is_read() {
    let dir = TestDir::new("report");
    let genome = letters(&fasta("LC738874.fasta"));
    // 600 reads: 14,400 letters, above the diagonal array's 8,000; 100 reads: 2,400, below it.
    let large = reads(&dir, &genome, 600);
    let small = reads(&dir, &genome, 100);
    let short: &[&str] = &["-task", "blastn-short"];
    let two_hit: &[&str] = &["-task", "blastn-short", "-window_size", "16"];
    let word8: &[&str] = &["-task", "blastn", "-word_size", "8"];
    let cases: [(&Path, &[&str], u32); 6] = [
        (&large, short, 1),
        (&large, short, 4),
        (&small, short, 1),
        (&large, two_hit, 1),
        (&large, two_hit, 4),
        (&large, word8, 1),
    ];
    for (query, task, threads) in cases {
        let what = format!("{query:?} {task:?} threads {threads}");
        let (default, default_stderr) = search(query, task, threads, &[]);
        assert_eq!(default_stderr, "", "{what}");
        assert!(!default.is_empty(), "{what}: no hits");
        let (old, old_stderr) = search(query, task, threads, &[SWITCH_OFF]);
        assert_eq!(old_stderr, "", "{what} =0");
        assert!(old == default, "{what}: the reports differ");
    }
}
