//! `LOSAT_X_TBXLUT` (in the strict set: on unless `=0`) on whole runs of the binary: the full
//! lookup table of a TBLASTX query batch, the table of the reference frame loop, filled in place
//! from per-cell counts (`x_tbx_build_ncbi_lookup` in `src/algorithm/tblastx/lookup/backbone.rs`,
//! `x_build_direct` in `x_lut_direct.rs`) instead of being copied out of one chain per cell.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:356-360
//! ```c
//!             for (j=0; j <lookup->thin_backbone[i][1]; j++)
//!                 dest[j] = lookup->thin_backbone[i][j + 2];
//!             /* done with this chain- free it */
//!             sfree(lookup->thin_backbone[i]);
//! ```
//! The table must be the same, so the report must be: the tests run the reference frame loop
//! (`LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1`, whose gate sends every batch to the frame loop at any
//! threshold and thread count) and compare the report with the switch (the default), whose in-place
//! build the run's reach file shows (`LOSAT_X_REACH`, `tests/reach_support`), against the report with
//! `LOSAT_X_TBXLUT=0`, at the default threshold (13), at 11 and at `+inf` (NCBI's INT_MIN: every
//! word enters every cell it can), at one and four threads; and on the default path, where the unit
//! path builds its own tables and the switch acts only on a batch the guard or
//! `LOSAT_X_PAIRPARGATE` sends back.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

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
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
        "losat-tblastx-ref-lut-{tag}-{}",
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

/// Writes `records` (defline, letters) to `dir/name`.
fn write(dir: &Path, name: &str, records: &[(&str, &[u8])]) -> PathBuf {
    let path = dir.join(name);
    let text: String = records.iter().map(|(d, l)| record(d, l)).collect();
    std::fs::write(&path, text).unwrap();
    path
}

/// The report (outfmt 6) and stderr of one run.
fn tblastx(
    query: &Path,
    subject: &Path,
    threshold: Option<&str>,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = losat(switches);
    command
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()]);
    if let Some(threshold) = threshold {
        command.args(["-threshold", threshold]);
    }
    command.env("LOSAT_X_STATS", "1");
    let output = command.output().expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threshold {threshold:?} threads {threads}: {stderr}"
    );
    if switches.contains(&REFERENCE) {
        // The chunk-switch gate sends every batch to the reference frame loop: no unit runs.
        let [units, fallback_gate] = pair_par(&stderr);
        assert_eq!(
            (units, fallback_gate >= 1),
            (0, true),
            "{switches:?} threshold {threshold:?} threads {threads}: the frame loop:\n{stderr}"
        );
    }
    (output.stdout, stderr)
}

/// `[X_STATS] PAIRPAR`: units, fallback_gate.
fn pair_par(stderr: &str) -> [u64; 2] {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] PAIRPAR "))
        .unwrap_or_else(|| panic!("no [X_STATS] PAIRPAR line in: {stderr}"));
    ["units", "fallback_gate"].map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("counter")
    })
}

const REFERENCE: (&str, &str) = ("LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_TBXLUT", "0");

/// Whether the run of `reach` built the full table in place (`LOSAT_X_TBXLUT` took its new side).
fn built_in_place(reach: &reach_support::Reach) -> bool {
    reach.sides("LOSAT_X_TBXLUT").0
}

/// Two queries in one batch (1,500 and 2,000 bases of MjeNMV) against two subjects (40 kb of
/// MelaMJNV, 30 kb of PemoMJNVA), threshold 13 and 11, one and four threads.
#[test]
fn the_reference_loop_reports_the_same_with_the_table_filled_in_place() {
    let dir = temp_dir("thresholds");
    let mje = letters("MjeNMV.fasta");
    let query = write(
        &dir,
        "q.fna",
        &[("mje_a", &mje[..1_500]), ("mje_b", &mje[20_000..22_000])],
    );
    let subject = write(
        &dir,
        "s.fna",
        &[
            ("mela", &letters("MelaMJNV.fasta")[..40_000]),
            ("pemo", &letters("PemoMJNVA.fasta")[..30_000]),
        ],
    );
    for threshold in [None, Some("11")] {
        for threads in [1, 4] {
            let (base, _) = tblastx(
                &query,
                &subject,
                threshold,
                threads,
                &[REFERENCE, SWITCH_OFF],
            );
            assert!(
                !base.is_empty(),
                "threshold {threshold:?}: no hits to compare"
            );
            let reach = reach_support::Reach::new("tbxlut-thresholds");
            let (with, _) = tblastx(
                &query,
                &subject,
                threshold,
                threads,
                &[REFERENCE, reach.var()],
            );
            assert_eq!(
                base, with,
                "LOSAT_X_TBXLUT, threshold {threshold:?}, threads {threads}"
            );
            assert!(
                built_in_place(&reach),
                "threshold {threshold:?}, threads {threads}: the table was not built in place"
            );
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// `-threshold +inf` (every word in every cell): a 450-base query against 3 kb, one thread.
#[test]
fn every_word_threshold_gives_the_same_report() {
    let dir = temp_dir("inf");
    let query = write(&dir, "q.fna", &[("mje", &letters("MjeNMV.fasta")[..450])]);
    let subject = write(
        &dir,
        "s.fna",
        &[("mela", &letters("MelaMJNV.fasta")[..3_000])],
    );
    let (base, _) = tblastx(&query, &subject, Some("+inf"), 1, &[REFERENCE, SWITCH_OFF]);
    let reach = reach_support::Reach::new("tbxlut-inf");
    let (with, _) = tblastx(&query, &subject, Some("+inf"), 1, &[REFERENCE, reach.var()]);
    assert_eq!(base, with, "LOSAT_X_TBXLUT at -threshold +inf");
    assert!(
        built_in_place(&reach),
        "-threshold +inf: the table was not built in place"
    );
    std::fs::remove_dir_all(&dir).ok();
}

/// The default path (the unit path of the strict set builds its own tables; at one thread the gate
/// sends this short query's batch to the reference loop, where the switch acts): the report is the
/// one of `LOSAT_X_TBXLUT=0`.
#[test]
fn the_default_path_reports_the_same_with_the_switch() {
    let dir = temp_dir("default");
    let query = write(&dir, "q.fna", &[("mje", &letters("MjeNMV.fasta")[..3_000])]);
    let subject = write(
        &dir,
        "s.fna",
        &[("mela", &letters("MelaMJNV.fasta")[..40_000])],
    );
    for threads in [1, 4] {
        let (base, _) = tblastx(&query, &subject, None, threads, &[SWITCH_OFF]);
        let (with, _) = tblastx(&query, &subject, None, threads, &[]);
        assert_eq!(base, with, "default path, threads {threads}");
    }
    std::fs::remove_dir_all(&dir).ok();
}
