//! The branch-free hit loop of the TBLASTX unit scan (`LOSAT_X_TBXHITBF`,
//! `src/algorithm/tblastx/x_tbx_hitbf.rs`, module doc) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. It acts on the units of the unit path that run without the
//! probes of `LOSAT_TIMING` / `LOSAT_X_STATS`, in `scan_unit` and in the batched scan of
//! `LOSAT_X_TBXEXTBATCH`; its shadow (`LOSAT_X_TBXHITBFSHADOW`) turns it on and compares it with the
//! loop before the switch on every unit.
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
//! the batched scan. The tests compare the default, `LOSAT_X_DEFAULT=0`, the batched scan, the
//! prefetch on every subject, today's lists and the shadow (which also shows that the loop ran and
//! that no unit fell outside its range) with `=0`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const OFF: (&str, &str) = ("LOSAT_X_TBXHITBF", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBXHITBFSHADOW", "1");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// The lines of `stderr` that start with `prefix`.
fn lines<'a>(stderr: &'a str, prefix: &str) -> Vec<&'a str> {
    stderr
        .lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

/// The number after `name=` in the first line of `stderr` that starts with `prefix`.
fn field(stderr: &str, prefix: &str, name: &str) -> u64 {
    let line = lines(stderr, prefix)
        .first()
        .copied()
        .unwrap_or_else(|| panic!("no {prefix} line in {stderr}"));
    let key = format!("{name}=");
    let at = line
        .find(&key)
        .unwrap_or_else(|| panic!("no {name} in {line}"))
        + key.len();
    line[at..]
        .split(|c: char| !c.is_ascii_digit())
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

#[test]
fn the_branch_free_loop_gives_the_same_report() {
    let inputs = Inputs::new("same");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        for switches in [
            &[][..],
            &[DEFAULT_OFF][..],
            &[FULL_LISTS][..],
            &[PREFETCH][..],
            &[BATCHED][..],
            &[BATCHED, FULL_LISTS][..],
            &[BATCHED, PREFETCH][..],
        ] {
            let (report, _) = inputs.tblastx(threads, switches);
            assert_eq!(report, off, "threads {threads} {switches:?}: against =0");
        }
        let (batched_off, _) = inputs.tblastx(threads, &[BATCHED, OFF]);
        assert_eq!(batched_off, off, "threads {threads}: batched =0 against =0");
    }
}

#[test]
fn the_shadow_compares_every_unit() {
    let inputs = Inputs::new("shadow");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        for switches in [
            &[SHADOW][..],
            &[SHADOW, FULL_LISTS][..],
            &[SHADOW, BATCHED][..],
            &[SHADOW, BATCHED, PREFETCH][..],
        ] {
            let (shadow, stderr) = inputs.tblastx(threads, switches);
            assert_eq!(
                shadow, off,
                "threads {threads} {switches:?}: shadow against =0"
            );
            assert!(
                field(&stderr, "[X_SHADOW] TBXHITBF", "units_compared") > 0,
                "threads {threads} {switches:?}: no unit compared: {stderr}"
            );
            assert!(
                field(&stderr, "[X_SHADOW] TBXHITBF", "records_compared") > 0,
                "threads {threads} {switches:?}: no record compared: {stderr}"
            );
            assert_eq!(
                field(&stderr, "[X_SHADOW] TBXHITBF", "outside"),
                0,
                "threads {threads} {switches:?}: a unit outside the range: {stderr}"
            );
        }
    }
}

#[test]
fn the_switch_off_prints_no_shadow_line() {
    let inputs = Inputs::new("quiet");
    for switches in [&[][..], &[OFF][..]] {
        let (_, stderr) = inputs.tblastx(1, switches);
        assert!(
            lines(&stderr, "[X_SHADOW] TBXHITBF").is_empty(),
            "{switches:?}: {stderr}"
        );
    }
}
