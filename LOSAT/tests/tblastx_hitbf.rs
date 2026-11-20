//! The branch-free hit loop of the TBLASTX unit scan (`LOSAT_X_TBXHITBF`,
//! `src/algorithm/tblastx/x_tbx_hitbf.rs`, module doc) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. It
//! acts on the units of the unit path that run without the probes of `LOSAT_TIMING` /
//! `LOSAT_X_STATS`, in `scan_unit` and in the batched scan of `LOSAT_X_TBXEXTBATCH`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:536-552
//! ```c
//! if (diff >= window) {
//!     diag_array[diag_coord].last_hit =
//!         subject_offset + diag_offset;
//!     continue;
//! }
//! if (diff < wordsize) {
//!     continue;
//! }
//! ```
//! The loop takes these two outcomes without a branch and hands the others to the per-hit body,
//! with the same cell states and the same saved init HSPs (the module's argument), so the property
//! under test is that the report does not depend on it. That reports equal NCBI's is the business of
//! the fixture comparisons. The input (150,000 bases of MelaMJNV against 150,000 of PemoMJNVA) has
//! frames of about 50,000 letters, so the unit path runs; `LOSAT_X_TBXEXTBATCH=2` sends its units to
//! the batched scan. The tests compare the default, the batched scan, the prefetch
//! on every subject and today's lists with `=0`; the reach rows (`LOSAT_X_REACH`,
//! `tests/reach_support`) show that the loop ran.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const OFF: (&str, &str) = ("LOSAT_X_TBXHITBF", "0");
/// The batched scan of LOSAT_X_TBXEXTBATCH on every unit (`x_hit_batched` for the rare hits).
const BATCHED: (&str, &str) = ("LOSAT_X_TBXEXTBATCH", "2");
/// Today's `InitHSP` lists instead of the compact records (`run_unit`, `run_unit_batched`).
const FULL_LISTS: (&str, &str) = ("LOSAT_X_SEEDCOMPACT", "0");
/// The prefetch of LOSAT_X_TBXPREFETCH on every subject (`x_hits_bf` with `PF`, x86-64).
const PREFETCH: (&str, &str) = ("LOSAT_X_TBXPREFETCH", "2");
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
        let dir =
            std::env::temp_dir().join(format!("losat-tblastx-hitbf-{test}-{}", std::process::id()));
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

#[test]
fn the_branch_free_loop_gives_the_same_report() {
    let inputs = Inputs::new("same");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        for switches in [
            &[][..],
            &[FULL_LISTS][..],
            &[PREFETCH][..],
            &[BATCHED][..],
            &[BATCHED, FULL_LISTS][..],
            &[BATCHED, PREFETCH][..],
        ] {
            // The loop ran on the units inside its range (`fits`).
            let reach = reach_support::Reach::new("same");
            let mut arm: Vec<(&str, &str)> = switches.to_vec();
            arm.push(reach.var());
            let (report, _) = inputs.tblastx(threads, &arm);
            assert_eq!(report, off, "threads {threads} {switches:?}: against =0");
            let sides = reach.sides("LOSAT_X_TBXHITBF");
            assert!(
                sides.0,
                "threads {threads} {switches:?}: sides {sides:?}: its path did not run"
            );
            assert!(
                !sides.1,
                "threads {threads} {switches:?}: sides {sides:?}: a subject took the reference path"
            );
        }
        let (batched_off, _) = inputs.tblastx(threads, &[BATCHED, OFF]);
        assert_eq!(batched_off, off, "threads {threads}: batched =0 against =0");
    }
}
