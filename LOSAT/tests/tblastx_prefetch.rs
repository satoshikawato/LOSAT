//! The prefetch of the TBLASTX unit scan (`LOSAT_X_TBXPREFETCH`, `src/algorithm/tblastx/x_pair_par.rs`,
//! module doc) on whole searches.
//!
//! The switch is a value switch of the strict set (`src/utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set; its number is the smallest largest-unit table, in KiB, of
//! a subject whose units scan with the hint (built-in 1,536; `2` takes every table). It acts on
//! x86-64 only. A prefetch has no architectural effect, so the report, and the counts the probes of
//! `LOSAT_TIMING` and `LOSAT_X_STATS` take, do not depend on it.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-516
//! ```c
//!     for (i = 0; i < hits; ++i) {
//!         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ```
//! NCBI's loop over a batch's hits issues no hint; the switch only adds one, so the property under
//! test is that nothing the search reports depends on it. That reports equal NCBI's is the business
//! of the fixture comparisons, which run on the default paths. The input (150,000 bases of MelaMJNV
//! against 150,000 of PemoMJNVA) has unit tables of about 0.4 MiB, below the built-in, so `=2`
//! puts the hint on every unit and the default leaves it off: the tests compare `=2` with `=0` and
//! with the default, through the unit functions without a probe, with `XStatsProbe` and with
//! `XTimeProbe`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const EVERY_TABLE: (&str, &str) = ("LOSAT_X_TBXPREFETCH", "2");
const OFF: (&str, &str) = ("LOSAT_X_TBXPREFETCH", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
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
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-prefetch-{test}-{}",
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
fn the_hint_on_every_table_gives_the_same_report_and_probe_counts() {
    let inputs = Inputs::new("every");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastx(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        let (every, _) = inputs.tblastx(threads, &[EVERY_TABLE]);
        assert_eq!(every, off, "threads {threads}: =2 against =0");
        let (default, _) = inputs.tblastx(threads, &[]);
        assert_eq!(default, off, "threads {threads}: default against =0");

        // XStatsProbe: the diagonal statistics count every hit, so they are the same.
        let (off_stats, off_err) = inputs.tblastx(threads, &[OFF, STATS]);
        let (every_stats, every_err) = inputs.tblastx(threads, &[EVERY_TABLE, STATS]);
        assert_eq!(
            (&off_stats, &every_stats),
            (&off, &off),
            "threads {threads}"
        );
        let seed = lines(&off_err, "[X_STATS] TBXSEED");
        assert!(!seed.is_empty(), "threads {threads}: no TBXSEED line");
        assert_eq!(
            lines(&every_err, "[X_STATS] TBXSEED"),
            seed,
            "threads {threads}"
        );

        // XTimeProbe: the unit path ran, with the same hits and extensions.
        let (off_timed, off_err) = inputs.tblastx(threads, &[OFF, TIMING]);
        let (every_timed, every_err) = inputs.tblastx(threads, &[EVERY_TABLE, TIMING]);
        assert_eq!(
            (&off_timed, &every_timed),
            (&off, &off),
            "threads {threads}"
        );
        assert!(
            field(&every_err, "[TIMING] x_pairpar:", "units") > 0,
            "threads {threads}: the unit path did not run"
        );
        for name in ["units", "hits", "extends"] {
            assert_eq!(
                field(&every_err, "[TIMING] x_unit:", name),
                field(&off_err, "[TIMING] x_unit:", name),
                "threads {threads}: x_unit {name}"
            );
        }
    }
}
