//! The context of a hit in the TBLASTX unit path from the `XContextTable` of blastp and blastx
//! (`LOSAT_X_TBXCTXTABLE`, `src/algorithm/tblastx/x_pair_par.rs`, module doc) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. It acts on the
//! units of the unit path: every hit that passes the two-hit test takes its context from a table
//! made once per batch instead of the binary search over the batch's frame bases.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
//! The table gives the binary search's index for every offset (the lib test
//! `the_context_table_gives_the_binary_search_index`), so the property under test is that the
//! report does not depend on the switch; that reports equal NCBI's is the business of the fixture
//! comparisons. The inputs (pieces of MelaMJNV against 60,000 bases of PemoMJNVA) have frames of
//! 5,000 letters and more, so the unit path runs: one query (six contexts) and three queries in one
//! batch (eighteen). The reach rows (`LOSAT_X_REACH`, `tests/reach_support`) show that the table's
//! path ran.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const OFF: (&str, &str) = ("LOSAT_X_TBXCTXTABLE", "0");

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

/// The query and subject files in a directory of the test's own: `pieces` pieces of 20,000 bases
/// of MelaMJNV as the queries, the first 60,000 bases of PemoMJNVA as the subject.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(test: &str, pieces: usize) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-ctx-table-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.fna");
        let subject = dir.join("s.fna");
        let mela = letters("MelaMJNV.fasta");
        let queries: String = (0..pieces)
            .map(|p| {
                record(
                    &format!("mela{p} MelaMJNV piece {p}"),
                    &mela[p * 20_000..(p + 1) * 20_000],
                )
            })
            .collect();
        std::fs::write(&query, queries).unwrap();
        std::fs::write(
            &subject,
            record(
                "pemo PemoMJNVA piece",
                &letters("PemoMJNVA.fasta")[..60_000],
            ),
        )
        .unwrap();
        Self {
            dir,
            query,
            subject,
        }
    }

    /// One TBLASTX search on standard output (`-outfmt 6`); returns the report.
    fn tblastx(&self, threads: u32, switches: &[(&str, &str)]) -> Vec<u8> {
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
        assert!(
            output.status.success(),
            "{switches:?} threads {threads}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn the_context_table_gives_the_same_report() {
    for pieces in [1, 3] {
        let inputs = Inputs::new(&format!("same{pieces}"), pieces);
        for threads in [1, 2] {
            let off = inputs.tblastx(threads, &[OFF]);
            assert!(
                !off.is_empty(),
                "{pieces} queries, threads {threads}: no hits"
            );
            let reach = reach_support::Reach::new("ctx");
            let report = inputs.tblastx(threads, &[reach.var()]);
            assert_eq!(
                report, off,
                "{pieces} queries, threads {threads}: against =0"
            );
            let sides = reach.sides("LOSAT_X_TBXCTXTABLE");
            assert!(
                sides.0,
                "{pieces} queries, threads {threads}: sides {sides:?}: its path did not run"
            );
        }
    }
}
