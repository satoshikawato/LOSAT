//! The private counters and presence-vector words of `LOSAT_X_PAIRPAR`'s megablast table fill
//! (`LOSAT_X_MBLUTPAD`, `src/algorithm/blastn/x_mb_lookup_par.rs`) on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set), and it acts only with
//! `LOSAT_X_MBLUTONCE=0` (the read-once fill, on by default since 2026-10-18, always counts into
//! the private copies), so every search here sets `LOSAT_X_MBLUTONCE=0`. It chooses between
//! two LOSAT code paths with the same tables: the pool fill counting into the shared
//! `helper_array` / `pv_array`, or into private copies added in at the end. So the property under
//! test is that the report does not depend on it, and that the padded fill really runs where the
//! switch says: with `LOSAT_TIMING`, `[TIMING] x_mb_lookup_build` counts `padded=` builds equal to
//! its `parallel=` builds (more than 0) at 4 threads with the switch, 0 without it and 0 at one
//! thread (no pool); with `LOSAT_X_PAIRPARSHADOW`, the shadow builds the second table of its
//! comparison with the padded fill and says so (`[X_SHADOW] MBLUTPAD tables=` more than 0). That
//! the bytes equal NCBI's is the business of the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH: (&str, &str) = ("LOSAT_X_MBLUTPAD", "1");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One search, `-outfmt 6` on standard output, with `LOSAT_X_MBLUTONCE=0` and the given variables
/// as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays, or is set to a file of this test when the
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
        .env("LOSAT_X_MBLUTONCE", "0")
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

/// For one search: the same report without the switch, with it, with it and `LOSAT_TIMING`, with
/// it and the shadow, and with it at one thread; the padded fill counted where it ran.
fn check(task: &str, query: &str, subject: &str) {
    let (reference, stderr) = blastn(task, query, subject, 4, &[]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert_eq!(stderr, "", "{task}");
    let (report, stderr) = blastn(task, query, subject, 4, &[SWITCH]);
    assert!(
        report == reference,
        "{task}: the report differs with the switch"
    );
    assert_eq!(stderr, "", "{task}: the switch writes nothing");

    let (report, stderr) = blastn(task, query, subject, 4, &[TIMING]);
    assert!(report == reference, "{task}");
    let parallel = field(&stderr, BUILD, "parallel=");
    assert!(parallel > 0, "{task}: no table on the pool:\n{stderr}");
    assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");

    let (report, stderr) = blastn(task, query, subject, 4, &[SWITCH, TIMING]);
    assert!(report == reference, "{task}");
    assert_eq!(field(&stderr, BUILD, "parallel="), parallel, "{task}");
    assert_eq!(
        field(&stderr, BUILD, "padded="),
        parallel,
        "{task}: every table on the pool padded:\n{stderr}"
    );

    let (report, stderr) = blastn(task, query, subject, 1, &[SWITCH, TIMING]);
    assert!(report == reference, "{task} at one thread");
    assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
    assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");

    let (report, stderr) = blastn(task, query, subject, 4, &[SWITCH, SHADOW]);
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
    let (_, stderr) = blastn(task, query, subject, 4, &[SHADOW]);
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
