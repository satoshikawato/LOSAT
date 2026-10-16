//! `LOSAT_X_PAIRPARGATE` on whole runs of the binary: the unit path of `LOSAT_X_PAIRPAR` runs only
//! for a batch where it is expected to be faster than the reference frame loop
//! (`units_gain` in `src/algorithm/tblastx/x_pair_par.rs`): a reference diagonal table of 2^18
//! cells or more, or several threads and one subject with more than one and a half thread's share
//! of the bases, or a lookup threshold below the default (13). Other batches run the reference
//! loop ("fallback: cost").
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-504
//! ```c
//! while (scan_range[1] <= scan_range[2]) {
//!     /* scan the subject sequence for hits */
//!     hits = scansub(lookup_wrap, subject,
//!                               offset_pairs, array_size, scan_range);
//! ```
//! Both paths give NCBI's list of each subject (module doc of x_pair_par.rs), so the switch may not
//! change a report. The tests compare the reports with and without it, and read
//! `[X_STATS] PAIRPAR ... fallback_cost=` and `units=` to see which path ran.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path, then `switches`.
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
    for (name, value) in switches {
        command.env(name, value);
    }
    command
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-tblastx-unit-gate-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
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

/// A file of one record, the first `length` bases of `name`.
fn query(dir: &Path, name: &str, length: usize) -> PathBuf {
    let file = dir.join(format!("q{length}.fna"));
    std::fs::write(&file, record("q", &letters(name)[..length])).unwrap();
    file
}

/// 40 pieces of 1,500 bases of the MelaMJNV and PemoMJNVA genomes (many short subjects).
fn pieces(dir: &Path) -> PathBuf {
    let mela = letters("MelaMJNV.fasta");
    let pemo = letters("PemoMJNVA.fasta");
    let mut set = String::new();
    for piece in 0..40 {
        let source = if piece % 2 == 0 { &mela } else { &pemo };
        let start = (piece * 7_919) % (source.len() - 1_500);
        set += &record(&format!("piece{piece}"), &source[start..start + 1_500]);
    }
    let file = dir.join("pieces.fna");
    std::fs::write(&file, set).unwrap();
    file
}

/// The report (`-outfmt 6`) and the counters of `[X_STATS] PAIRPAR`: (units, fallback_cost).
fn tblastx(
    query: &Path,
    subject: &Path,
    threads: u32,
    extra: &[&str],
    switches: &[(&str, &str)],
) -> (Vec<u8>, u64, u64) {
    let mut all = vec![("LOSAT_X_STATS", "1")];
    all.extend_from_slice(switches);
    let output = losat(&all)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .args(extra)
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] PAIRPAR "))
        .unwrap_or_else(|| panic!("no PAIRPAR counters:\n{stderr}"));
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .unwrap()
    };
    (output.stdout, field("units="), field("fallback_cost="))
}

const GATE: (&str, &str) = ("LOSAT_X_PAIRPARGATE", "1");

/// The reports with and without the gate are the same; returns (units, fallback_cost) with it.
fn with_and_without(query: &Path, subject: &Path, threads: u32, extra: &[&str]) -> (u64, u64) {
    let (reference, units_without, cost_without) = tblastx(query, subject, threads, extra, &[]);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 5,
        "the search finds hits"
    );
    assert!(units_without > 0, "without the gate the units run");
    assert_eq!(cost_without, 0);
    let (report, units, cost) = tblastx(query, subject, threads, extra, &[GATE]);
    assert!(
        report == reference,
        "threads {threads} {}: the report differs with the gate",
        query.display()
    );
    (units, cost)
}

// A short query (diagonal table of 2^14 cells) against many short subjects: the reference loop at
// one and at several threads.
#[test]
fn a_short_query_against_short_subjects_takes_the_reference_loop() {
    let dir = temp_dir("short");
    let q = query(&dir, "MjeNMV.fasta", 6_000);
    let set = pieces(&dir);
    for threads in [1, 2, 4] {
        assert_eq!(
            with_and_without(&q, &set, threads, &[]),
            (0, 1),
            "threads {threads}"
        );
    }
}

// Below the default threshold the gate does not act: the same search keeps the units.
#[test]
fn a_threshold_below_the_default_keeps_the_units() {
    let dir = temp_dir("threshold");
    let q = query(&dir, "MjeNMV.fasta", 6_000);
    let set = pieces(&dir);
    let (units, cost) = with_and_without(&q, &set, 1, &["-threshold", "12"]);
    assert!(units > 0);
    assert_eq!(cost, 0);
}

// The same query against one genome: the reference loop at one thread, the units (which split the
// subject over the pool) at several.
#[test]
fn one_long_subject_keeps_the_units_at_several_threads() {
    let dir = temp_dir("one");
    let q = query(&dir, "MjeNMV.fasta", 6_000);
    let genome = fasta("PemoMJNVA.fasta");
    assert_eq!(with_and_without(&q, &genome, 1, &[]), (0, 1));
    for threads in [2, 4] {
        let (units, cost) = with_and_without(&q, &genome, threads, &[]);
        assert!(units > 0, "threads {threads}");
        assert_eq!(cost, 0, "threads {threads}");
    }
}

// A query of 70 kb (a table of 2^18 cells) keeps the units at one thread.
#[test]
fn a_long_query_keeps_the_units() {
    let dir = temp_dir("long");
    let q = query(&dir, "MelaMJNV.fasta", 70_000);
    let set = pieces(&dir);
    let (units, cost) = with_and_without(&q, &set, 1, &[]);
    assert!(units > 0);
    assert_eq!(cost, 0);
}
