//! The read-once pool fill of the megablast lookup table (`LOSAT_X_MBLUTONCE`,
//! `src/algorithm/blastn/x_mb_lookup_once.rs`) on whole searches.
//!
//! The switch is in the strict set since 2026-10-18 (on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set). Where `LOSAT_X_PAIRPAR` fills the contiguous megablast table on
//! a pool of two threads or more it chooses between two LOSAT code paths with the same tables: each
//! worker reads every word (`x_mb_lookup_par.rs`), or each position is read once and each word
//! applied by the owner of its cell. So the property under test is that the report does not depend
//! on it, and that the read-once fill really runs where the switch says: with `LOSAT_TIMING`,
//! `[TIMING] x_mb_lookup_build` counts `once=` builds equal to its `parallel=` builds (more than
//! 0) at 4 threads with the switch (also with `LOSAT_X_MBLUTPAD`, which it takes precedence over),
//! with `[TIMING] x_mb_lookup_once` giving as many fills, and so without the variable (the
//! default); 0 with `LOSAT_X_MBLUTONCE=0`, with `LOSAT_X_DEFAULT=0`
//! (which turns `LOSAT_X_PAIRPAR` off) and at one thread (no pool); with `LOSAT_X_PAIRPARSHADOW`,
//! the shadow builds the second table of its comparison with the read-once fill and says so
//! (`[X_SHADOW] MBLUTONCE tables=` more than 0, `diffs=0 identical`), and says nothing with
//! `LOSAT_X_MBLUTONCE=0`. That the bytes equal NCBI's is the business of the blastn fixtures.
//!
//! The prefault round of the read-once fill (`LOSAT_X_MBLUTPREFAULT`, in the strict set since
//! 2026-10-19) runs where the first chunk has at least 8 positions per page of the table: here
//! megablast of LC738868 (306 kb, a table of 4^12 cells, 16,385 pages) at 8 threads (chunks of
//! 262,144 positions), counted by `prefault_rounds=` / `prefault_touches=` of `[TIMING]
//! x_mb_lookup_once`; not at 4 threads (131,072 positions), nor with `LOSAT_X_MBLUTPREFAULT=0`; the
//! report is the same either way, and the shadow's read-once table equals the loop's with it.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH: (&str, &str) = ("LOSAT_X_MBLUTONCE", "1");
const OFF: (&str, &str) = ("LOSAT_X_MBLUTONCE", "0");
const PAD: (&str, &str) = ("LOSAT_X_MBLUTPAD", "1");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARSHADOW", "1");
const PREFAULT_OFF: (&str, &str) = ("LOSAT_X_MBLUTPREFAULT", "0");

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
            std::env::temp_dir().join(format!("losat_mblutonce_worker_{}.log", std::process::id())),
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
const ONCE: &str = "[TIMING] x_mb_lookup_once: ";
const SHADOW_LINE: &str = "[X_SHADOW] MBLUTONCE ";

/// For one search: the same report with the switch off, with it, by default, with it and
/// `LOSAT_X_DEFAULT=0`, with `LOSAT_TIMING`, with the shadow, and at one thread; the read-once fill
/// counted where it ran.
fn check(task: &str, query: &str, subject: &str) {
    let (reference, stderr) = blastn(task, query, subject, 4, &[OFF]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert_eq!(stderr, "", "{task}");
    for switches in [&[SWITCH][..], &[], &[SWITCH, DEFAULT_OFF], &[SWITCH, PAD]] {
        let (report, stderr) = blastn(task, query, subject, 4, switches);
        assert!(
            report == reference,
            "{task} {switches:?}: the report differs"
        );
        assert_eq!(stderr, "", "{task} {switches:?}: the switch writes nothing");
    }

    // With the switch off: the fill of `x_mb_lookup_par` makes every pool build.
    let (report, stderr) = blastn(task, query, subject, 4, &[OFF, TIMING]);
    assert!(report == reference, "{task}");
    let parallel = field(&stderr, BUILD, "parallel=");
    assert!(parallel > 0, "{task}: no table on the pool:\n{stderr}");
    assert_eq!(field(&stderr, BUILD, "once="), 0, "{task}: {stderr}");
    assert!(!stderr.contains(ONCE), "{task}: {stderr}");
    // `LOSAT_X_DEFAULT=0` turns `LOSAT_X_PAIRPAR` off, so no table is filled on the pool.
    let (report, stderr) = blastn(task, query, subject, 4, &[SWITCH, DEFAULT_OFF, TIMING]);
    assert!(report == reference, "{task}");
    assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
    assert_eq!(field(&stderr, BUILD, "once="), 0, "{task}: {stderr}");

    for switches in [&[SWITCH, TIMING][..], &[TIMING], &[SWITCH, PAD, TIMING]] {
        let (report, stderr) = blastn(task, query, subject, 4, switches);
        assert!(report == reference, "{task} {switches:?}");
        assert_eq!(field(&stderr, BUILD, "parallel="), parallel, "{task}");
        assert_eq!(
            field(&stderr, BUILD, "once="),
            parallel,
            "{task} {switches:?}: every table on the pool read once:\n{stderr}"
        );
        assert_eq!(field(&stderr, BUILD, "padded="), 0, "{task}: {stderr}");
        assert_eq!(field(&stderr, ONCE, "fills="), parallel, "{task}: {stderr}");
        assert!(
            field(&stderr, ONCE, "chunks=") >= parallel,
            "{task}: {stderr}"
        );
    }

    let (report, stderr) = blastn(task, query, subject, 1, &[SWITCH, TIMING]);
    assert!(report == reference, "{task} at one thread");
    assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
    assert_eq!(field(&stderr, BUILD, "once="), 0, "{task}: {stderr}");
    assert!(!stderr.contains(ONCE), "{task}: {stderr}");

    let (report, stderr) = blastn(task, query, subject, 4, &[SWITCH, SHADOW]);
    assert!(report == reference, "{task} with the shadow");
    let line = stderr
        .lines()
        .find(|line| line.starts_with(SHADOW_LINE))
        .unwrap_or_else(|| panic!("{task}: no shadow line:\n{stderr}"));
    assert!(line.ends_with(" diffs=0 identical"), "{line}");
    assert!(
        field(&stderr, SHADOW_LINE, "tables=") > 0,
        "{task}: the shadow compared no read-once table: {line}"
    );
    assert!(field(&stderr, SHADOW_LINE, "cells=") > 0, "{line}");
    assert!(
        !stderr.contains("[X_SHADOW] MBLUTPAD"),
        "{task}: no padded build:\n{stderr}"
    );
    let (_, stderr) = blastn(task, query, subject, 4, &[OFF, SHADOW]);
    assert!(
        !stderr.contains(SHADOW_LINE),
        "{task}: the line needs the switch:\n{stderr}"
    );
}

#[test]
fn megablast_gives_the_same_report_and_reads_its_table_once() {
    check("megablast", "LC738868.fasta", "LC738868.fasta");
}

#[test]
fn blastn_gives_the_same_report_and_reads_its_table_once() {
    check("blastn", "LC738868.fasta", "LC738869.fasta");
}

#[test]
fn megablast_prefaults_the_table_where_the_first_chunk_is_dense() {
    let (task, query) = ("megablast", "LC738868.fasta");
    let (reference, _) = blastn(task, query, query, 8, &[OFF]);
    assert!(!reference.is_empty());
    let (report, stderr) = blastn(task, query, query, 8, &[TIMING]);
    assert!(report == reference, "the report differs with the round");
    let fills = field(&stderr, ONCE, "fills=");
    assert!(fills > 0, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), fills, "{stderr}");
    let touches = field(&stderr, ONCE, "prefault_touches=");
    assert!(
        touches > 16_000 && touches <= 16_385,
        "nearly every page of the 4^12 table: {stderr}"
    );
    let (report, stderr) = blastn(task, query, query, 8, &[PREFAULT_OFF, TIMING]);
    assert!(report == reference, "the report differs without the round");
    assert_eq!(field(&stderr, ONCE, "fills="), fills, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), 0, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_touches="), 0, "{stderr}");
    // 4 threads: a first chunk of 131,072 positions, below 8 per page of 16,385 pages.
    let (report, stderr) = blastn(task, query, query, 4, &[TIMING]);
    assert!(report == reference, "the report differs at 4 threads");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), 0, "{stderr}");
    let (report, stderr) = blastn(task, query, query, 8, &[SHADOW]);
    assert!(report == reference, "the report differs with the shadow");
    let line = stderr
        .lines()
        .find(|line| line.starts_with(SHADOW_LINE))
        .unwrap_or_else(|| panic!("no shadow line:\n{stderr}"));
    assert!(line.ends_with(" diffs=0 identical"), "{line}");
    assert!(field(&stderr, SHADOW_LINE, "tables=") > 0, "{line}");
}
