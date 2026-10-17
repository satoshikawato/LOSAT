//! The private counters and presence-vector words of `LOSAT_X_PAIRPAR`'s megablast table fill
//! (`LOSAT_X_MBLUTPAD`, `src/algorithm/blastn/x_mb_lookup_par.rs`) on whole searches.
//!
//! The switch is in the strict set and on by default (since 2026-10-17): `LOSAT_X_MBLUTPAD=0` or
//! `LOSAT_X_DEFAULT=0` turns it off. It chooses between two LOSAT code paths with the same tables:
//! the pool fill counting into the shared `helper_array` / `pv_array`, or into private copies added
//! in at the end. So the property under test is that the report does not depend on it, and that
//! the padded fill really runs where the switch says: with `LOSAT_TIMING`, `[TIMING]
//! x_mb_lookup_build` counts `padded=` builds equal to its `parallel=` builds (more than 0) at 4
//! threads by default and with the switch set, 0 with it off (`=0`, or `LOSAT_X_DEFAULT=0`) and 0
//! at one thread (no pool); with `LOSAT_X_PAIRPARSHADOW`, the shadow builds the second table of
//! its comparison with the padded fill and says so (`[X_SHADOW] MBLUTPAD tables=` more than 0),
//! and says nothing with the switch off. That the bytes equal NCBI's is the business of the blastn
//! fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH: (&str, &str) = ("LOSAT_X_MBLUTPAD", "1");
const OFF: (&str, &str) = ("LOSAT_X_MBLUTPAD", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One search, `-outfmt 6` on standard output, with the given variables as its only `LOSAT_*`
/// variables. `LOSAT_BLASTX_WORKER_LOG` stays, or is set to a file of this test when the
/// environment has none: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one thread,
/// and it selects nothing. Returns the report and the standard error.
fn blastn(
    task: &str,
    query: &str,
    subject: &str,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!("losat_mblutpad_worker_{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subject))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The value of `name` (`name=VALUE`) in the line of `stderr` that starts with `prefix`.
fn field(stderr: &str, prefix: &str, name: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(prefix))
        .unwrap_or_else(|| panic!("no {prefix} line:\n{stderr}"));
    line.split(|c: char| c.is_whitespace() || c == '(' || c == ')')
        .find_map(|word| word.strip_prefix(name))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

const BUILD: &str = "[TIMING] x_mb_lookup_build: ";

/// For one search: the same report with the switch off, by default, with it set, with
/// `LOSAT_X_DEFAULT=0`, with `LOSAT_TIMING`, with the shadow, and at one thread; the padded fill
/// counted where it ran.
fn check(task: &str, query: &str, subject: &str) {
    let (reference, stderr) = blastn(task, query, subject, 4, &[OFF]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert_eq!(stderr, "", "{task}");
    for switches in [&[][..], &[SWITCH], &[DEFAULT_OFF]] {
        let (report, stderr) = blastn(task, query, subject, 4, switches);
        assert!(
            report == reference,
            "{task} {switches:?}: the report differs"
        );
        assert_eq!(stderr, "", "{task} {switches:?}: the switch writes nothing");
    }

    let (report, stderr) = blastn(task, query, subject, 4, &[OFF, TIMING]);
    assert!(report == reference, "{task}");
    let parallel = field(&stderr, BUILD, "parallel=");
    assert!(parallel > 0, "{task}: no table on the pool:\n{stderr}");
    assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");
    // `LOSAT_X_DEFAULT=0` also turns `LOSAT_X_PAIRPAR` off, so no table is filled on the pool.
    let (report, stderr) = blastn(task, query, subject, 4, &[DEFAULT_OFF, TIMING]);
    assert!(report == reference, "{task}");
    assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
    assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");

    for switches in [&[TIMING][..], &[SWITCH, TIMING]] {
        let (report, stderr) = blastn(task, query, subject, 4, switches);
        assert!(report == reference, "{task} {switches:?}");
        assert_eq!(field(&stderr, BUILD, "parallel="), parallel, "{task}");
        assert_eq!(
            field(&stderr, BUILD, "padded="),
            parallel,
            "{task} {switches:?}: every table on the pool padded:\n{stderr}"
        );
    }

    let (report, stderr) = blastn(task, query, subject, 1, &[TIMING]);
    assert!(report == reference, "{task} at one thread");
    assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
    assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");

    let (report, stderr) = blastn(task, query, subject, 4, &[SHADOW]);
    assert!(report == reference, "{task} with the shadow");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] MBLUTPAD "))
        .unwrap_or_else(|| panic!("{task}: no shadow line:\n{stderr}"));
    assert!(line.ends_with(" diffs=0 identical"), "{line}");
    assert!(
        field(&stderr, "[X_SHADOW] MBLUTPAD ", "tables=") > 0,
        "{task}: the shadow compared no padded table: {line}"
    );
    assert!(
        field(&stderr, "[X_SHADOW] MBLUTPAD ", "cells=") > 0,
        "{line}"
    );
    let (_, stderr) = blastn(task, query, subject, 4, &[OFF, SHADOW]);
    assert!(
        !stderr.contains("[X_SHADOW] MBLUTPAD"),
        "{task}: the line needs the switch:\n{stderr}"
    );
}

#[test]
fn megablast_gives_the_same_report_and_pads_its_table_fill() {
    check("megablast", "LC738868.fasta", "LC738868.fasta");
}

#[test]
fn blastn_gives_the_same_report_and_pads_its_table_fill() {
    check("blastn", "LC738868.fasta", "LC738869.fasta");
}
