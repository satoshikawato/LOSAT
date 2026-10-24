//! The composition-based redo of tblastn taken match by match, largest first
//! (`LOSAT_X_TBNREDOLPT`, `src/algorithm/tblastn/stage_d_pipeline.rs`, `x_redo_lpt`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. It acts where
//! the redo runs on the search pool (two or more threads, two or more (query, subject) matches),
//! on each chunk with two or more matches to redo; at one thread the redo stays on the calling
//! thread.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3448-3451
//! ```c
//! #pragma omp for schedule(static)
//!         for (b = 0; b < numMatches; ++b) {
//! #pragma omp flush(interrupt)
//!             if (!interrupt) {
//! ```
//! The order changes only which thread redoes which match, and when (the unit tests of
//! `stage_d_pipeline.rs` compare the scope with the `par_iter` path in one process); either way
//! the matches are applied in list order. So the property under test is that the report depends
//! neither on the switch nor on the thread count; that reports equal NCBI's is the business of the
//! fixture comparisons. `LOSAT_X_REACH` shows the side each run took. The input (the first 40
//! proteins of AvCLPV against its genome cut into 42 pieces of 10 kb) has more matches than one
//! chunk holds at 2 threads (64) and at 4 (128); with `-max_target_seqs 1` early termination
//! leaves matches of later chunks out before they start.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_TBNREDOLPT";
const OFF: (&str, &str) = ("LOSAT_X_TBNREDOLPT", "0");
/// The proteins of the query.
const PROTEINS: usize = 40;
/// The bases of each subject piece.
const PIECE: usize = 10_000;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

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

/// The query file (the first `PROTEINS` records of AvCLPV.faa) and the subject file (the AvCLPV
/// genome in pieces of `PIECE` bases), in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastn-redo-lpt-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.faa");
        let text = std::fs::read_to_string(fasta("AvCLPV.faa")).unwrap();
        let end = text
            .match_indices('>')
            .nth(PROTEINS)
            .map_or(text.len(), |(at, _)| at);
        std::fs::write(&query, &text[..end]).unwrap();
        let subject = dir.join("pieces.fna");
        let genome = std::fs::read_to_string(fasta("AvCLPV.fasta")).unwrap();
        let bases: String = genome
            .lines()
            .skip(1)
            .flat_map(|line| line.trim().chars())
            .collect();
        let pieces: String = bases
            .as_bytes()
            .chunks(PIECE)
            .enumerate()
            .map(|(index, piece)| format!(">p{index}\n{}\n", std::str::from_utf8(piece).unwrap()))
            .collect();
        std::fs::write(&subject, pieces).unwrap();
        Self {
            dir,
            query,
            subject,
        }
    }

    /// One tblastn search on standard output (`-outfmt 6`) keeping `targets` subjects per query;
    /// returns the report and the sides of the switch the run took.
    fn tblastn(
        &self,
        targets: u32,
        threads: u32,
        switches: &[(&str, &str)],
    ) -> (Vec<u8>, (bool, bool)) {
        let reach = reach_support::Reach::new("tblastn_redo_lpt");
        let mut variables = switches.to_vec();
        variables.push(reach.var());
        let output = losat(&variables)
            .arg("tblastn")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subject)
            .args(["-max_target_seqs", &targets.to_string()])
            .args(["-outfmt", "6"])
            .args(["-num_threads", &threads.to_string()])
            .output()
            .expect("run LOSAT");
        assert!(
            output.status.success(),
            "targets {targets} {switches:?} threads {threads}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        (output.stdout, reach.sides(SWITCH))
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// At 2 and 4 threads the redo takes the new path, at 1 thread the reference; every report is
/// the one-thread report of `=0`, and `=0` records no side.
#[test]
fn the_redo_largest_first_gives_the_same_report() {
    let inputs = Inputs::new("same");
    for targets in [500, 1] {
        let (reference, _) = inputs.tblastn(targets, 1, &[OFF]);
        assert!(!reference.is_empty(), "targets {targets}: no hits");
        for (threads, side) in [(1, (false, true)), (2, (true, false)), (4, (true, false))] {
            let (off, off_side) = inputs.tblastn(targets, threads, &[OFF]);
            assert_eq!(
                off_side,
                (false, false),
                "targets {targets} threads {threads}: =0's side"
            );
            assert_eq!(
                off, reference,
                "targets {targets} threads {threads}: =0 against one thread"
            );
            let (report, reached) = inputs.tblastn(targets, threads, &[]);
            assert_eq!(
                reached, side,
                "targets {targets} threads {threads}: default's side"
            );
            assert_eq!(
                report, reference,
                "targets {targets} threads {threads}: against =0"
            );
        }
    }
}
