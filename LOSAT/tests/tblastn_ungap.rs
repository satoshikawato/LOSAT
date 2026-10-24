//! The byte-table walks of tblastn's two-hit ungapped extension (`LOSAT_X_TBNUNGAP`,
//! `src/algorithm/tblastn/search_init.rs`, `x_tbn_ungap`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`.
//! It acts on every word finder whose matrix is BLOSUM62 (tblastn's default); another `-matrix`
//! keeps the reference extension.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1088-1096
//! ```c
//! static Int4
//! s_BlastAaExtendTwoHit(Int4 ** matrix,
//!                       const BLAST_SequenceBlk * subject,
//!                       const BLAST_SequenceBlk * query,
//!                       Int4 s_left_off, Int4 s_right_off,
//!                       Int4 q_right_off, Int4 dropoff,
//!                       Int4 * hsp_q, Int4 * hsp_s, Int4 * hsp_len,
//!                       Boolean use_pssm, Int4 word_size,
//!                       Boolean * right_extend, Int4 * s_last_off)
//! ```
//! The walks give the reference's result for every extension (the unit test of `search_init.rs`
//! compares them with the reference in one process), so the property under test is that the report
//! does not depend on the switch; that reports equal NCBI's is the business of the fixture
//! comparisons. `LOSAT_X_REACH` shows the side each run took. The input (the first 40 proteins of
//! AvCLPV against its genome) gives extensions and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_TBNUNGAP";
const OFF: (&str, &str) = ("LOSAT_X_TBNUNGAP", "0");
/// The proteins of the query.
const PROTEINS: usize = 40;

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

/// The query file (the first `PROTEINS` records of AvCLPV.faa), in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-tblastn-ungap-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.faa");
        let text = std::fs::read_to_string(fasta("AvCLPV.faa")).unwrap();
        let end = text
            .match_indices('>')
            .nth(PROTEINS)
            .map_or(text.len(), |(at, _)| at);
        std::fs::write(&query, &text[..end]).unwrap();
        Self { dir, query }
    }

    /// One tblastn search of the AvCLPV genome on standard output (`-outfmt 6`) with `matrix`;
    /// returns the report and the sides of the switch the run took.
    fn tblastn(
        &self,
        matrix: &str,
        threads: u32,
        switches: &[(&str, &str)],
    ) -> (Vec<u8>, (bool, bool)) {
        let reach = reach_support::Reach::new("tblastn_ungap");
        let mut variables = switches.to_vec();
        variables.push(reach.var());
        let output = losat(&variables)
            .arg("tblastn")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(fasta("AvCLPV.fasta"))
            .args(["-matrix", matrix])
            .args(if matrix == "BLOSUM45" {
                // LOSAT's TBLASTN takes BLOSUM45 with these (its message for other values).
                &["-word_size", "2", "-comp_based_stats", "0"][..]
            } else {
                &[][..]
            })
            .args(["-outfmt", "6"])
            .args(["-num_threads", &threads.to_string()])
            .output()
            .expect("run LOSAT");
        assert!(
            output.status.success(),
            "{matrix} {switches:?} threads {threads}: {}",
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

/// BLOSUM62 takes the byte-table walks and BLOSUM45 the reference, at 1 and 2 threads; the
/// default's report is `=0`'s, and `=0` records no side.
#[test]
fn the_byte_table_walks_give_the_same_report() {
    let inputs = Inputs::new("same");
    for (matrix, side) in [("BLOSUM62", (true, false)), ("BLOSUM45", (false, true))] {
        for threads in [1, 2] {
            let (off, off_side) = inputs.tblastn(matrix, threads, &[OFF]);
            assert!(!off.is_empty(), "{matrix} threads {threads}: no hits");
            assert_eq!(
                off_side,
                (false, false),
                "{matrix} threads {threads}: =0's side"
            );
            let (report, reached) = inputs.tblastn(matrix, threads, &[]);
            assert_eq!(reached, side, "{matrix} threads {threads}: default's side");
            assert_eq!(report, off, "{matrix} threads {threads}: against =0");
        }
    }
}
