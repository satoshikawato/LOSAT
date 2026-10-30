//! TBLASTX's unit path with the groups of contexts the pool needs (`LOSAT_X_TBXFEWGROUPS`,
//! `src/algorithm/tblastx/x_pair_par.rs`, `x_groups_wanted`, `group_starts`) on whole searches.
//!
//! No NCBI counterpart: the groups are how LOSAT schedules a batch's contexts, and NCBI scans them
//! as one concatenated query. The switch is in the strict set (on unless it is `0`), so the
//! property under test is that the report does not depend on it: at one thread with one group
//! (at most 2^17 letters in the batch), at one thread with the six groups (a longer batch), at four
//! threads with two groups. That reports equal NCBI's is the business of the TBLASTX fixtures.
//! `LOSAT_X_REACH` shows which side each run took.
//!
//! The query is cut from AP027131, the subject from the related AP027133 (genetic code 4 for
//! both). `LOSAT_X_PAIRPARGATE=0` sends the batch to the units whatever its length.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_TBXFEWGROUPS";
const OFF: (&str, &str) = ("LOSAT_X_TBXFEWGROUPS", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_tblastx_few_groups_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The first `n` bases of the one record of a FASTA file, as FASTA text.
fn bases(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut lines = text.lines();
    let header = lines.next().expect("a header");
    let sequence: String = lines.flat_map(str::chars).take(n).collect();
    assert_eq!(sequence.len(), n, "{name} is shorter than {n} bases");
    format!("{header}\n{sequence}\n")
}

/// The query and subject files, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(query_bases: usize, subject_bases: usize) -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("query.fna");
        let subject = dir.join("subject.fna");
        std::fs::write(&query, bases("AP027131.fasta", query_bases)).unwrap();
        std::fs::write(&subject, bases("AP027133.fasta", subject_bases)).unwrap();
        Inputs {
            dir,
            query,
            subject,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The sides of `SWITCH` that a run with `LOSAT_X_REACH=dir` recorded: (new, ref).
fn reached(dir: &Path) -> (bool, bool) {
    let (mut new, mut reference) = (false, false);
    for file in std::fs::read_dir(dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            if f[0] == SWITCH {
                new |= f[1] == "1";
                reference |= f[2] == "1";
            }
        }
    }
    (new, reference)
}

/// `tblastx` of the inputs at `threads` with `env`: the report and the sides of `SWITCH` recorded.
/// `LOSAT_BLASTX_WORKER_LOG` names a file of its own: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn tblastx(
    inputs: &Inputs,
    threads: u32,
    env: &[(&str, &str)],
    args: &[&str],
) -> (Vec<u8>, (bool, bool)) {
    let dir = temp_path("reach");
    std::fs::create_dir_all(&dir).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .env("LOSAT_X_PAIRPARGATE", "0")
        .env("LOSAT_X_REACH", &dir)
        .envs(env.iter().copied())
        .arg("tblastx")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(["-query_gencode", "4", "-db_gencode", "4", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .args(args)
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "{env:?} threads {threads}: {stderr}"
    );
    let sides = reached(&dir);
    let _ = std::fs::remove_dir_all(&dir);
    (output.stdout, sides)
}

/// A batch of 12 kb takes fewer groups at one and four threads, with the report of `=0`; with
/// another matrix too (its tables and walks are per group).
#[test]
fn one_group_gives_the_report_of_the_six() {
    let inputs = Inputs::new(12_000, 12_000);
    for args in [&[][..], &["-matrix", "BLOSUM45"][..]] {
        let (expected, sides) = tblastx(&inputs, 1, &[OFF], args);
        assert!(!expected.is_empty(), "no hits");
        assert_eq!(sides, (false, false), "=0 records nothing");
        let (report, sides) = tblastx(&inputs, 1, &[], args);
        assert_eq!(report, expected, "{args:?}: one group");
        assert_eq!(
            sides,
            (true, false),
            "{args:?}: the side that one thread took"
        );
        let (report, sides) = tblastx(&inputs, 4, &[], args);
        assert_eq!(report, expected, "{args:?}: four threads");
        assert_eq!(
            sides,
            (true, false),
            "{args:?}: the side that four threads took"
        );
    }
}

/// A batch above 2^17 letters (a query of 70 kb) keeps its groups at one thread.
#[test]
fn a_longer_batch_keeps_its_groups() {
    let inputs = Inputs::new(70_000, 6_000);
    let (expected, _) = tblastx(&inputs, 1, &[OFF], &[]);
    assert!(!expected.is_empty(), "no hits");
    let (report, sides) = tblastx(&inputs, 1, &[], &[]);
    assert_eq!(report, expected);
    assert_eq!(sides, (false, true), "the side that a long batch took");
}

/// Records of the first `total` bases of `name`, `pieces` of equal length, as FASTA text.
fn pieces_of(name: &str, total: usize, pieces: usize) -> String {
    let record = bases(name, total);
    let sequence = record.lines().nth(1).expect("the sequence");
    let length = total / pieces;
    (0..pieces)
        .map(|k| format!(">piece{k}\n{}\n", &sequence[k * length..(k + 1) * length]))
        .collect()
}

/// Many queries against many subjects, the gate on (the batches take the units from 5,000 letters
/// at several threads, one group): the report of `=0` (the reference frame loop's at several
/// threads), of one thread, of a chromosome-like subject among short ones. The threads differ in
/// their order of work only.
#[test]
fn queries_against_subjects_give_the_same_report_with_the_gate() {
    let inputs = Inputs::new(12_000, 6_000);
    std::fs::write(&inputs.query, pieces_of("AP027131.fasta", 12_000, 4)).unwrap();
    std::fs::write(&inputs.subject, pieces_of("AP027133.fasta", 12_000, 8)).unwrap();
    let gate = ("LOSAT_X_PAIRPARGATE", "1");
    let (expected, sides) = tblastx(&inputs, 4, &[gate, OFF], &[]);
    assert!(!expected.is_empty(), "no hits");
    assert_eq!(sides, (false, false), "=0 records nothing");
    let (report, sides) = tblastx(&inputs, 4, &[gate], &[]);
    assert_eq!(report, expected, "four threads");
    assert_eq!(sides, (true, false), "one group for the subjects");
    let (report, _) = tblastx(&inputs, 1, &[gate], &[]);
    assert_eq!(report, expected, "one thread");
    let mut mixed = bases("AP027133.fasta", 12_000);
    mixed.push_str(">short\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
    std::fs::write(&inputs.subject, &mixed).unwrap();
    let (expected, _) = tblastx(&inputs, 4, &[gate, OFF], &[]);
    let (report, _) = tblastx(&inputs, 4, &[gate], &[]);
    assert_eq!(report, expected, "one long subject");
}
