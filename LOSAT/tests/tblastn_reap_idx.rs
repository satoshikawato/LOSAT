//! The reap of contained HSPs after tblastn's composition-based redo with an index by subject start
//! (`LOSAT_X_TBNREAPIDX`, `src/algorithm/tblastn/kappa/x_reap_index.rs`) on a whole search.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. A match whose redo
//! gives `X_REAP_INDEX_FROM` (32) HSPs or more takes the index; shorter lists keep NCBI's loop
//! (`s_HitlistReapContained`, `blast_kappa.c` 224-273). The kept HSPs are the loop's (the unit tests
//! of `x_reap_index.rs` compare them on random lists), so the property under test is that the report
//! does not depend on the switch; that reports equal NCBI's is the business of the fixture
//! comparisons. `LOSAT_X_REACH` shows the side each run took; whether the index removes an HSP is
//! the unit tests' business (their random lists remove more than 10,000). Input: AP027131 proteins
//! BDV02435.1, BDV02433.1, BDV02434.1 and BDV02232.1 (redo lists of 32 or more against the AP027133
//! genome: the index) and BDV02287.1 (a short list: NCBI's loop).

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_TBNREAPIDX";
const OFF: (&str, &str) = ("LOSAT_X_TBNREAPIDX", "0");

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

/// The records of a FASTA fixture with the given identifiers, in that order.
fn records(fixture: &str, ids: &[&str]) -> String {
    let text = std::fs::read_to_string(fasta(fixture)).expect("FASTA");
    let mut out = String::new();
    for id in ids {
        let start = text
            .match_indices('>')
            .map(|(at, _)| at)
            .find(|&at| text[at + 1..].split_whitespace().next() == Some(id))
            .unwrap_or_else(|| panic!("{id} not in {fixture}"));
        let end = text[start + 1..]
            .find('>')
            .map_or(text.len(), |at| start + 1 + at);
        out.push_str(&text[start..end]);
    }
    out
}

/// The query file, in a directory of the test's own.
struct Query {
    dir: PathBuf,
    path: PathBuf,
}

impl Query {
    fn new(text: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-tblastn-reap-idx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.faa");
        std::fs::write(&path, text).unwrap();
        Self { dir, path }
    }
}

impl Drop for Query {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// One search (`-outfmt 6`); returns the report and the sides of the switch the run took.
fn search(query: &Path, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, (bool, bool)) {
    let reach = reach_support::Reach::new("tblastn_reap_idx");
    let mut variables = switches.to_vec();
    variables.push(reach.var());
    let output = losat(&variables)
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta("AP027133.fasta"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (output.stdout, reach.sides(SWITCH))
}

/// The default report equals `=0`'s at 1 and 2 threads, the default takes both the index and NCBI's
/// loop, and `=0` records no side.
#[test]
fn tblastn_reap_with_the_index_gives_the_same_report() {
    let query = Query::new(&records(
        "AP027131.faa",
        &[
            "BDV02435.1",
            "BDV02433.1",
            "BDV02434.1",
            "BDV02232.1",
            "BDV02287.1",
        ],
    ));
    for threads in [1, 2] {
        let (off, off_side) = search(&query.path, threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        assert_eq!(off_side, (false, false), "threads {threads}: =0's side");
        let (report, side) = search(&query.path, threads, &[]);
        assert_eq!(side, (true, true), "threads {threads}: default's side");
        assert_eq!(report, off, "threads {threads}: against =0");
    }
}
