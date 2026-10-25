//! The table walk of TBLASTX's preliminary bases of a `-subject` record with a letter other than
//! A, C, G and T (`LOSAT_X_TBXRESOLVE`, `src/algorithm/tblastx/subjects.rs`,
//! `x_resolved_record_frames`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. It acts on each
//! such record while `LOSAT_X_TBXFRAMES` is on; a record of A, C, G and T only copies its letters,
//! and a `-db` sequence is read from its stored codes.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:443-450
//! ```c
//!         if (c  != 0xFF) {
//!             // No ambiguities, so we can do this the easy way.
//!             ncbi2na[i] = c;
//!
//!         } else {
//!             if (b == 0 || b == 0x0F) {
//!                 //gap or N
//!                 ncbi2na[i] = random.GetRand() & 0x3;
//! ```
//! The walk gives the letter-by-letter walk's frames for every record (the unit test of
//! `subjects.rs` compares them in one process), so the property under test is that the report does
//! not depend on the switch; that reports equal NCBI's is the business of the fixture comparisons
//! (`tblastx_regression_fixtures.py`'s `ambig.*` and `input.rna_*`, outfmt0's `tblastx.ambig.*`).
//! `LOSAT_X_REACH` shows the side each run took. The input (the first 6 kb of MjeNMV against the
//! MelaMJNV genome with IUPAC codes, lower-case letters and N blocks, one across a translation
//! window, and the PemoMJNVA genome as it is) gives HSPs in both records.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_TBXRESOLVE";
const OFF: (&str, &str) = ("LOSAT_X_TBXRESOLVE", "0");

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

fn letters(name: &str) -> Vec<u8> {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

fn record(defline: &str, letters: &[u8]) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The query and subject files, in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subjects: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-resolve-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.fna");
        std::fs::write(&query, record("mje6k", &letters("MjeNMV.fasta")[..6000])).unwrap();
        let mut mela = letters("MelaMJNV.fasta");
        let iupac = b"NRYKMSWBDHVnrykmswbdhv";
        for (n, at) in (1_000..mela.len()).step_by(997).enumerate() {
            mela[at] = iupac[n % iupac.len()];
        }
        // Lower-case letters where most HSPs are.
        mela[2_000..8_000].make_ascii_lowercase();
        // N blocks, one across the border of the first translation window (196,608 bases).
        for (start, len) in [(20_000, 300), (196_500, 250)] {
            mela[start..start + len].fill(b'N');
        }
        let subjects = dir.join("s.fna");
        std::fs::write(
            &subjects,
            record("mela MelaMJNV with ambiguity codes", &mela)
                + &record("pemo PemoMJNVA", &letters("PemoMJNVA.fasta")),
        )
        .unwrap();
        Self {
            dir,
            query,
            subjects,
        }
    }

    /// One tblastx search on standard output (`-outfmt 6`); returns the report and the sides of
    /// the switch the run took.
    fn tblastx(&self, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, (bool, bool)) {
        let reach = reach_support::Reach::new("tblastx_resolve");
        let mut variables = switches.to_vec();
        variables.push(reach.var());
        let output = losat(&variables)
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subjects)
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
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The record with ambiguous letters takes the table walk and the other one copies its letters (at
/// 1 and 2 threads); the default's report is `=0`'s, with hits in both records, and `=0` records
/// no side.
#[test]
fn the_table_walk_gives_the_same_report() {
    let inputs = Inputs::new("same");
    for threads in [1, 2] {
        let (off, off_side) = inputs.tblastx(threads, &[OFF]);
        let text = String::from_utf8_lossy(&off);
        for subject in ["mela", "pemo"] {
            assert!(
                text.lines()
                    .any(|line| line.split('\t').nth(1) == Some(subject)),
                "threads {threads}: no hit in {subject}"
            );
        }
        assert_eq!(off_side, (false, false), "threads {threads}: =0's side");
        let (report, reached) = inputs.tblastx(threads, &[]);
        assert_eq!(reached, (true, true), "threads {threads}: default's side");
        assert_eq!(report, off, "threads {threads}: against =0");
    }
}
