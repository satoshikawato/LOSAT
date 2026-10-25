//! The composition-based redo's kept alignments with an index by subject start
//! (`LOSAT_X_REDOCONTAINIDX`, `src/core/composition_adjustment/x_align_index.rs`) on whole blastp,
//! blastx and tblastn searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. In a redo, a
//! query's list of kept alignments moves into the index once it holds `X_INDEX_FROM` (32)
//! alignments; shorter lists stay NCBI's list. `s_IsContained` and `s_WithDistinctEnds`
//! (`redo_alignment.c` 386-451, 966-997) then look only at the alignments whose subject range can
//! meet the new one; the answers, and the final list with its order, are the list's (the unit
//! tests of `x_align_index.rs` compare them call by call on random lists). So the property under
//! test is that the report does not depend on the switch; that reports equal NCBI's is the business
//! of the fixture comparisons. `LOSAT_X_REACH` shows the side each run took. Inputs, each with a
//! (query, subject) match whose redo keeps more than 32 alignments:
//! - blastp: WSSV proteins YP_009220644.1 (103 HSPs against PajaWSV BDT63508.1), YP_009220531.1
//!   and YP_009220475.1 against PajaWSV BDT63508.1 and BDT63501.1 (short lists: the reference side);
//! - tblastn: AP027131 proteins BDV02435.1, BDV02433.1, BDV02434.1 and BDV02232.1 against the
//!   AP027133 genome (33-40 HSPs each);
//! - blastx: bases 356,000-400,000 of AP027131 against AP027133 proteins BDV03502.1, BDV03503.1,
//!   BDV03665.1 and BDV03795.1.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_REDOCONTAINIDX";
const OFF: (&str, &str) = ("LOSAT_X_REDOCONTAINIDX", "0");

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

/// Bases `from..to` of the single record of a nucleotide fixture, as one record.
fn bases(fixture: &str, from: usize, to: usize) -> String {
    let text = std::fs::read_to_string(fasta(fixture)).expect("FASTA");
    let sequence: String = text
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let mut out = format!(">{fixture}_{from}_{to}\n");
    for line in sequence.as_bytes()[from..to].chunks(70) {
        out.push_str(std::str::from_utf8(line).expect("ASCII"));
        out.push('\n');
    }
    out
}

/// The query and subject files of one case, in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-redo-contain-idx-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self { dir }
    }

    fn write(&self, name: &str, text: &str) -> PathBuf {
        let path = self.dir.join(name);
        std::fs::write(&path, text).unwrap();
        path
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// One search (`-outfmt 6`); returns the report and the sides of the switch the run took.
fn search(
    program: &str,
    query: &Path,
    subject: &Path,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, (bool, bool)) {
    let reach = reach_support::Reach::new("redo_contain_idx");
    let mut variables = switches.to_vec();
    variables.push(reach.var());
    let output = losat(&variables)
        .arg(program)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{program} {switches:?} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (output.stdout, reach.sides(SWITCH))
}

/// The default report equals `=0`'s at 1 and 2 threads, the default takes the index (and the
/// list where `reference` says so), and `=0` records no side.
fn same_report(program: &str, query: &Path, subject: &Path, reference: bool) {
    for threads in [1, 2] {
        let (off, off_side) = search(program, query, subject, threads, &[OFF]);
        assert!(!off.is_empty(), "{program} threads {threads}: no hits");
        assert_eq!(
            off_side,
            (false, false),
            "{program} threads {threads}: =0's side"
        );
        let (report, side) = search(program, query, subject, threads, &[]);
        assert_eq!(
            side,
            (true, reference),
            "{program} threads {threads}: default's side"
        );
        assert_eq!(report, off, "{program} threads {threads}: against =0");
    }
}

#[test]
fn blastp_with_the_index_gives_the_same_report() {
    let inputs = Inputs::new("blastp");
    let query = inputs.write(
        "q.faa",
        &records(
            "WSSV.faa",
            &["YP_009220644.1", "YP_009220531.1", "YP_009220475.1"],
        ),
    );
    let subject = inputs.write(
        "s.faa",
        &records("PajaWSV.faa", &["BDT63508.1", "BDT63501.1"]),
    );
    same_report("blastp", &query, &subject, true);
}

#[test]
fn tblastn_with_the_index_gives_the_same_report() {
    let inputs = Inputs::new("tblastn");
    let query = inputs.write(
        "q.faa",
        &records(
            "AP027131.faa",
            &["BDV02435.1", "BDV02433.1", "BDV02434.1", "BDV02232.1"],
        ),
    );
    same_report("tblastn", &query, &fasta("AP027133.fasta"), false);
}

#[test]
fn blastx_with_the_index_gives_the_same_report() {
    let inputs = Inputs::new("blastx");
    let query = inputs.write("q.fa", &bases("AP027131.fasta", 356_000, 400_000));
    let subject = inputs.write(
        "s.faa",
        &records(
            "AP027133.faa",
            &["BDV03502.1", "BDV03503.1", "BDV03665.1", "BDV03795.1"],
        ),
    );
    same_report("blastx", &query, &subject, true);
}
