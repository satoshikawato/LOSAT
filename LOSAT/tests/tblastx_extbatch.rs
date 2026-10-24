//! The batched two-hit extensions of the TBLASTX unit scan (`LOSAT_X_TBXEXTBATCH`,
//! `src/algorithm/tblastx/x_tbx_extbatch.rs`, module doc) on whole searches.
//!
//! The switch is a value switch of the strict set (`src/utils/xswitch.rs`): on unless it is `0`;
//! its number is the smallest largest-unit query-group span, in Ki letters, of a subject whose
//! units take the path (built-in 256; `2` takes every unit of 2 Ki letters or more). It acts on
//! x86-64 only (on other targets `=2` compares the old path with itself).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-605
//! ```c
//! score = s_BlastAaExtendTwoHit(matrix, subject, query,
//!                               last_hit + wordsize,
//!                               subject_offset, query_offset, ...);
//! ```
//! NCBI extends at once; the path extends later, in batches, with the same cell states and the
//! same saved init HSPs (the module's argument), so the property under test is that the report does
//! not depend on it. That reports equal NCBI's is the business of the fixture comparisons. The input
//! (150,000 bases of MelaMJNV against 150,000 of PemoMJNVA) has frames of about 50,000 letters, below
//! the built-in, so `=2` puts the path on every unit and the default leaves it off: the tests compare
//! `=2` (with the compact records and with today's lists) and the default with `=0`; the reach
//! rows of `=2` show that the path ran (`LOSAT_X_REACH`, `tests/reach_support`).

use std::path::{Path, PathBuf};
use std::process::Command;

#[allow(clippy::duplicate_mod)]
mod reach_support;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const EVERY_UNIT: (&str, &str) = ("LOSAT_X_TBXEXTBATCH", "2");
const OFF: (&str, &str) = ("LOSAT_X_TBXEXTBATCH", "0");
/// Today's `InitHSP` lists instead of the compact records (`run_unit_batched`).
const FULL_LISTS: (&str, &str) = ("LOSAT_X_SEEDCOMPACT", "0");
/// The bases of each piece: frames of about 50,000 letters, so the unit path runs.
const BASES: usize = 150_000;

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
}

/// The letters of a one-record FASTA file.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

/// A FASTA record of `letters` in lines of 70.
fn record(defline: &str, letters: &str) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.as_bytes().chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The query and subject files: the first `BASES` of MelaMJNV and of PemoMJNVA, in a directory of
/// the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-extbatch-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.fna");
        let subject = dir.join("s.fna");
        std::fs::write(
            &query,
            record("mela MelaMJNV piece", &letters("MelaMJNV.fasta")[..BASES]),
        )
        .unwrap();
        std::fs::write(
            &subject,
            record("pemo PemoMJNVA piece", &letters("PemoMJNVA.fasta")[..BASES]),
        )
        .unwrap();
        Self {
            dir,
            query,
            subject,
        }
    }

    /// One TBLASTX search on standard output (`-outfmt 6`); returns the report and the standard
    /// error.
    fn tblastx(&self, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
        let output = losat(switches)
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subject)
            .args(["-outfmt", "6"])
            .args(["-num_threads", &threads.to_string()])
            .output()
            .expect("run LOSAT");
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(
            output.status.success(),
            "{switches:?} threads {threads}: {stderr}"
        );
        (output.stdout, stderr)
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// Asserts that the batched extensions ran in the search that wrote `reach` (the switch acts on
/// x86-64 only).
fn assert_batched(reach: &reach_support::Reach, what: &str) {
    let sides = reach.sides("LOSAT_X_TBXEXTBATCH");
    assert_eq!(
        sides.0,
        cfg!(target_arch = "x86_64"),
        "{what}: sides {sides:?}"
    );
}

#[test]
fn every_unit_batched_gives_the_same_report() {
    let inputs = Inputs::new("every");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        let reach = reach_support::Reach::new("every");
        let (every, _) = inputs.tblastx(threads, &[EVERY_UNIT, reach.var()]);
        assert_eq!(every, off, "threads {threads}: =2 against =0");
        assert_batched(&reach, &format!("threads {threads}: =2"));
        let (default, _) = inputs.tblastx(threads, &[]);
        assert_eq!(default, off, "threads {threads}: default against =0");
        let reach = reach_support::Reach::new("full");
        let (full, _) = inputs.tblastx(threads, &[EVERY_UNIT, FULL_LISTS, reach.var()]);
        assert_eq!(
            full, off,
            "threads {threads}: =2 with full lists against =0"
        );
        assert_batched(&reach, &format!("threads {threads}: =2 with full lists"));
    }
}
