//! The hit lists' arrays grown as lists come (`LOSAT_X_BNHLGROW`, `src/algorithm/blastn/hsp.rs`
//! `HitList::x_room_for_one`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3219-3240
//! ```c
//!     const int kStartValue = 100;
//! ...
//!     if (hit_list->hsplist_current <= 0)
//!        hit_list->hsplist_current = kStartValue;
//!     else
//!        hit_list->hsplist_current =
//!           MIN(2*hit_list->hsplist_current, hit_list->hsplist_max);
//!     hit_list->hsplist_array =
//!        (BlastHSPList**) realloc(hit_list->hsplist_array,
//!                                 hit_list->hsplist_current*sizeof(BlastHSPList*));
//! ```
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`. It changes only how much room a hit list's array holds (the lib test
//! `hit_list_array_grows_as_lists_come_within_ncbi_array_size`), never which lists it keeps, so
//! the property under test is that the report does not depend on it: the default run and
//! `LOSAT_X_BNHLGROW=0` give the same bytes and print nothing on standard error, with subjects
//! read from FASTA (`-subject`) and from a database (`-db`), at 1 and 4 threads, with lists that
//! pass NCBI's first array size (100 subjects) and with lists cut by `-max_target_seqs 50`
//! (the heap of a full hit list). That those bytes equal NCBI's is the business of the fixture
//! suites.
//!
//! The queries are 300 reads of 18 to 30 letters every 27 letters of the first 8,100 letters of
//! LC738874, every third one with its middle letter changed; the subjects 150 windows of 4,000
//! letters every 30 letters of the same genome, so a read lies in up to 134 subjects.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BNHLGROW", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-hlgrow-{test}-{}", std::process::id()));
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

/// The reads, the subject windows as FASTA, and their database made by `losat makedb`.
fn inputs(dir: &TestDir) -> (PathBuf, PathBuf, String) {
    let genome = letters(&fasta("LC738874.fasta"));
    let mut reads = String::new();
    for i in 0..300 {
        let start = i * 27;
        let len = 18 + i % 13;
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
        .arg("makedb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, set, db)
}

/// One blastn-short search (`-outfmt 6` on standard output) with the given variables as its only
/// `LOSAT_*` variables that choose a code path. Returns the report and the standard error.
fn search(
    query: &Path,
    target: &[&str],
    extra: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .args(target)
        .args(["-task", "blastn-short", "-outfmt", "6"])
        .args(extra)
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{target:?} {extra:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The largest number of subjects one query has rows for.
fn most_subjects(report: &[u8]) -> usize {
    let mut per_query = std::collections::HashMap::<&str, std::collections::HashSet<&str>>::new();
    for line in std::str::from_utf8(report).unwrap().lines() {
        let mut fields = line.split('\t');
        let (query, subject) = (fields.next().unwrap(), fields.next().unwrap());
        per_query.entry(query).or_default().insert(subject);
    }
    per_query.values().map(|s| s.len()).max().unwrap_or(0)
}

// The default search grows the hit lists' arrays as lists come; `LOSAT_X_BNHLGROW=0` reserves
// them to NCBI's array size. The reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_how_the_hit_lists_grow() {
    let dir = TestDir::new("report");
    let (query, set, db) = inputs(&dir);
    let set = set.to_str().unwrap();
    for target in [["-subject", set], ["-db", db.as_str()]] {
        for extra in [&[][..], &["-max_target_seqs", "50"][..]] {
            for threads in [1, 4] {
                let what = format!("{target:?} {extra:?} threads {threads}");
                let (default, default_stderr) = search(&query, &target, extra, threads, &[]);
                assert_eq!(default_stderr, "", "{what}");
                let most = most_subjects(&default);
                if extra.is_empty() {
                    assert!(most > 100, "{what}: at most {most} subjects per query");
                } else {
                    assert_eq!(most, 50, "{what}");
                }
                let (old, old_stderr) = search(&query, &target, extra, threads, &[SWITCH_OFF]);
                assert_eq!(old_stderr, "", "{what} =0");
                assert!(old == default, "{what}: the reports differ");
            }
        }
    }
}
