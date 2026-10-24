//! The unit scan of TBLASTX in slices of the query group (`LOSAT_X_TBXQSLICE`,
//! `src/algorithm/tblastx/x_tbx_qslice.rs`, module doc) on whole searches.
//!
//! The switch is a value switch of the strict set (`src/utils/xswitch.rs`): on unless it is `0`;
//! its number is the most letters of a slice (`1` gives the built-in). It acts on x86-64 only (on
//! other targets every search runs the unsliced scan).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
//! ```c
//! while (scan_range[1] <= scan_range[2]) {
//!     /* scan the subject sequence for hits */
//!     hits = scansub(lookup_wrap, subject,
//!                               offset_pairs, array_size, scan_range);
//! ```
//! NCBI scans each subject range once with the whole query's table; the path scans it once per
//! slice, with the same hits in the same order on each diagonal (the module's argument), so the
//! property under test is that the report does not depend on it. That reports equal NCBI's is the
//! business of the fixture comparisons. The input (150,000 bases of MelaMJNV against 150,000 of
//! PemoMJNVA) has frames of about 50,000 letters, below the built-in, so the default cuts nothing
//! and `=4096` cuts each group into about 13 slices: the tests compare `=4096` (with the compact
//! records and with today's lists, with the reference table builder, with `LOSAT_X_TBXEXTBATCH` and
//! `LOSAT_X_TBXPREFETCH` on every unit) and the default with `=0`; the reach rows (`LOSAT_X_REACH`,
//! `tests/reach_support`) show that the groups were cut.

use std::path::{Path, PathBuf};
use std::process::Command;

#[allow(clippy::duplicate_mod)]
mod reach_support;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SLICES: (&str, &str) = ("LOSAT_X_TBXQSLICE", "4096");
const OFF: (&str, &str) = ("LOSAT_X_TBXQSLICE", "0");
/// Today's `InitHSP` lists instead of the compact records (`run_unit_sliced`).
const FULL_LISTS: (&str, &str) = ("LOSAT_X_SEEDCOMPACT", "0");
/// The batched extensions on every unit of 2 Ki letters or more.
const BATCHED: (&str, &str) = ("LOSAT_X_TBXEXTBATCH", "2");
/// The prefetch on every table of 2 KiB or more.
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
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-qslice-{test}-{}",
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

#[test]
fn sliced_units_give_the_same_report() {
    let inputs = Inputs::new("slices");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        for switches in [
            &[SLICES][..],
            &[SLICES, FULL_LISTS][..],
            &[SLICES, BATCHED][..],
            &[SLICES, PREFETCH][..],
            &[][..],
        ] {
            // `=4096` cuts each group on x86-64; the default cuts nothing on this input.
            let reach = reach_support::Reach::new("slices");
            let mut arm: Vec<(&str, &str)> = switches.to_vec();
            arm.push(reach.var());
            let (report, _) = inputs.tblastx(threads, &arm);
            assert_eq!(report, off, "threads {threads} {switches:?}: against =0");
            let sides = reach.sides("LOSAT_X_TBXQSLICE");
            assert_eq!(
                sides.0,
                cfg!(target_arch = "x86_64") && switches.contains(&SLICES),
                "threads {threads} {switches:?}: sides {sides:?}"
            );
        }
    }
}
