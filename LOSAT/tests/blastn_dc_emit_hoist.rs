//! dc-megablast's strips scan of a table laid out as CSR with `LOSAT_X_DCEMITHOIST`
//! (`src/algorithm/blastn/blast_engine/x_pair_seed_disc.rs`: `scan_strip_disc_csr`'s `HOIST`
//! instance, `DcPairs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes how the scan makes each word's seeds from the word's run of query offsets (what
//! `seed_of_pair` reads for every pair, read once per strip or word), not which seeds it makes or
//! their order, so the property under test is that the report does not depend on it, and that the
//! default took the hoisted CSR emit where the strips scan reads a CSR (`[TIMING] x_dc_emit_hoist:
//! strips=`) and not when `=0` or where it walks chains. MG1655 (4.6 Mb, one table of about 9.2
//! million words, above the CSR's size rule) against AP027131: at one thread (the strips stage on
//! that thread, `LOSAT_X_BNDCONE`, into one list per strip), and at four (`LOSAT_X_PAIRPARFOLD`'s
//! partitions: the seed pieces of `LOSAT_X_SEEDPIECE`, pieces of 2 seeds that stop inside a word's
//! run, and with `LOSAT_X_SEEDPIECE=0` the partitions' lists). LC738874 against LC738868 (a table
//! of about 580,000 words, below the size rule: chains) at one thread.
//! This test compares the default with `=0` only and runs no NCBI BLAST+. That the `=0` bytes
//! equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DCEMITHOIST", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One dc-megablast search of `query` against `subject` (`-outfmt 6` on standard output) with the
/// given variables as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built
/// with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file
/// name as soon as a search has more than one thread, and it selects nothing. Returns `(stdout,
/// stderr)`.
fn search(query: &str, subject: &str, threads: u32, vars: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .args(["-task", "dc-megablast"])
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subject))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{query} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The `strips=` count on `[TIMING] x_dc_emit_hoist:` in `stderr`.
fn hoist_strips(stderr: &str) -> u64 {
    let line = stderr
        .lines()
        .find_map(|line| line.strip_prefix("[TIMING] x_dc_emit_hoist:"))
        .unwrap_or_else(|| panic!("no x_dc_emit_hoist timing line in:\n{stderr}"));
    line.split_whitespace()
        .find_map(|field| field.strip_prefix("strips="))
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("strips= in {line}"))
}

// The default and `=0` give the same report; the default's strips scans read the CSR through the
// hoisted emit (`strips > 0`), `=0` records `strips == 0`, and the stage that runs the scan is the
// one named: the strips stage (`LOSAT_X_PAIRPARDISC`) at one thread, its seed pieces
// (`LOSAT_X_SEEDPIECE`) at four, or the partitions' lists (`LOSAT_X_PAIRPARHOME`'s home loop with
// `LOSAT_X_SEEDPIECE=0`: no piece).
#[test]
fn dc_megablast_csr_scan_with_the_hoisted_emit() {
    let runs: [(u32, &[(&str, &str)], &str); 4] = [
        (1, &[], "LOSAT_X_PAIRPARDISC"),
        (4, &[], "LOSAT_X_SEEDPIECE"),
        (4, &[("LOSAT_X_SEEDPIECE", "2")], "LOSAT_X_SEEDPIECE"),
        (4, &[("LOSAT_X_SEEDPIECE", "0")], "LOSAT_X_PAIRPARHOME"),
    ];
    for (threads, vars, stage) in runs {
        let at = format!("threads {threads} {vars:?}");
        let mut off_vars = vars.to_vec();
        off_vars.extend([SWITCH_OFF, ("LOSAT_TIMING", "1")]);
        let (expected, expected_stderr) =
            search("MG1655.fna", "AP027131.fasta", threads, &off_vars);
        assert!(!expected.is_empty(), "{at}: no hits");
        assert_eq!(hoist_strips(&expected_stderr), 0, "{at}: =0 strips");
        let reach = reach_support::Reach::new("dcemithoist");
        let mut on_vars = vars.to_vec();
        on_vars.extend([reach.var(), ("LOSAT_TIMING", "1")]);
        let (report, stderr) = search("MG1655.fna", "AP027131.fasta", threads, &on_vars);
        assert!(report == expected, "{at}: default differs from =0");
        assert!(hoist_strips(&stderr) > 0, "{at}: strips == 0");
        assert!(reach.sides("LOSAT_X_DCCSR").0, "{at}: no CSR table");
        assert!(reach.sides(stage).0, "{at}: {stage} did not run");
        let pieces = threads > 1 && !vars.contains(&("LOSAT_X_SEEDPIECE", "0"));
        assert_eq!(
            reach.sides("LOSAT_X_SEEDPIECE").0,
            pieces,
            "{at}: seed pieces"
        );
    }
}

// Below the CSR's size rule the strips scan walks the chains: `strips == 0`, and the same report
// as `=0`.
#[test]
fn dc_megablast_chains_scan_takes_the_reference_side() {
    let (expected, expected_stderr) = search(
        "LC738874.fasta",
        "LC738868.fasta",
        1,
        &[SWITCH_OFF, ("LOSAT_TIMING", "1")],
    );
    assert!(!expected.is_empty(), "no hits");
    assert_eq!(hoist_strips(&expected_stderr), 0, "=0 strips");
    let reach = reach_support::Reach::new("dcemithoist-chains");
    let (report, stderr) = search(
        "LC738874.fasta",
        "LC738868.fasta",
        1,
        &[reach.var(), ("LOSAT_TIMING", "1")],
    );
    assert!(report == expected, "default differs from =0");
    assert_eq!(hoist_strips(&stderr), 0, "chains strips");
    assert_eq!(reach.sides("LOSAT_X_DCCSR"), (false, true), "a CSR table");
    assert!(reach.sides("LOSAT_X_PAIRPARDISC").0, "no strips stage");
}
