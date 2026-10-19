//! The megablast hashtable kept from batch to batch of blastn's batch loop (`LOSAT_X_MBLUTKEEP`,
//! `src/algorithm/blastn/x_mb_lut_keep.rs`) on whole searches.
//!
//! The switch is in the strict set since 2026-10-19 (on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set). It chooses between two LOSAT code paths with the same tables: each batch
//! allocates a zeroed table and frees it, or the loop keeps the largest table and the next build of
//! that size zeroes and fills it. So the property under test is that the report does not depend on
//! it, and that the kept path really runs where the switch says: with `LOSAT_TIMING`, `[TIMING]
//! x_mb_lut_keep` counts `kept=` builds (at least one) and `fresh=` builds that add up to the
//! search's builds (`[TIMING] x_mb_lookup_build: builds=`), at 4 threads (zeroing and fill on the
//! pool) and at one thread; the prefault round of `LOSAT_X_MBLUTPREFAULT` runs only on the fresh
//! tables (`prefault_rounds=` of `[TIMING] x_mb_lookup_once` equals `fresh=` where every first
//! chunk is dense enough, and the builds without the switch); with `LOSAT_X_MBLUTKEEPSHADOW` every
//! build on a kept table is compared with a fresh build (`[X_SHADOW] MBLUTKEEP tables=` = `kept=`,
//! `diffs=0 identical`); no line with `LOSAT_X_MBLUTKEEP=0` or `LOSAT_X_DEFAULT=0`, kept builds
//! with the switch set and `LOSAT_X_DEFAULT=0` (the original fill loop); and no table is kept while
//! `LOSAT_X_QBATCHPAR`'s slots may search batches. That the bytes equal NCBI's is the business of
//! the blastn fixtures.
//!
//! The query file holds four genomes of 287-306 kb (LC738868, LC738869, LC738870, LC738873) and
//! `BATCH_SIZE=4000` puts each in a batch of its own (a batch holds whole records): four tables of
//! the same width per search, so the last three builds are on the kept table. The subject is
//! LC738868.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const OFF: (&str, &str) = ("LOSAT_X_MBLUTKEEP", "0");
const SWITCH: (&str, &str) = ("LOSAT_X_MBLUTKEEP", "1");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_MBLUTKEEPSHADOW", "1");
const BATCH: (&str, &str) = ("BATCH_SIZE", "4000");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
const SLOTS2: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");

const BUILD: &str = "[TIMING] x_mb_lookup_build: ";
const ONCE: &str = "[TIMING] x_mb_lookup_once: ";
const KEEP: &str = "[TIMING] x_mb_lut_keep: ";
const SHADOW_LINE: &str = "[X_SHADOW] MBLUTKEEP ";
const QBATCHPAR_LINE: &str = "[X_STATS] QBATCHPAR ";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The four genomes in one query file (written once per test binary).
fn queries() -> &'static Path {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        let mut text = Vec::new();
        for name in [
            "LC738868.fasta",
            "LC738869.fasta",
            "LC738870.fasta",
            "LC738873.fasta",
        ] {
            text.extend(std::fs::read(fasta(name)).expect("read fixture"));
            if text.last() != Some(&b'\n') {
                text.push(b'\n');
            }
        }
        let path =
            std::env::temp_dir().join(format!("losat_mblutkeep_q_{}.fasta", std::process::id()));
        std::fs::write(&path, text).expect("write the query file");
        path
    })
}

/// One search, `-outfmt 6` on standard output, with the given variables as its only `LOSAT_*`
/// variables (and batch-size variables). `LOSAT_BLASTX_WORKER_LOG` stays, or is set to a file of
/// this test when the environment has none: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing. Returns the report and the standard error.
fn blastn(task: &str, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
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
            std::env::temp_dir().join(format!("losat_mblutkeep_worker_{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(queries())
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
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
        .and_then(|value| value.trim_end_matches('s').parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

/// For one task: the same report with the switch off, by default, with it, with it and
/// `LOSAT_X_DEFAULT=0`, with `LOSAT_TIMING` and with the shadow at 4 threads and one; the kept
/// builds counted where they ran.
fn check(task: &str) {
    let (reference, stderr) = blastn(task, 4, &[OFF, BATCH]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert_eq!(stderr, "", "{task}");
    for switches in [
        &[BATCH][..],
        &[SWITCH, BATCH],
        &[SWITCH, DEFAULT_OFF, BATCH],
    ] {
        let (report, stderr) = blastn(task, 4, switches);
        assert!(
            report == reference,
            "{task} {switches:?}: the report differs"
        );
        assert_eq!(stderr, "", "{task} {switches:?}: the switch writes nothing");
    }

    for threads in [4, 1] {
        let (report, stderr) = blastn(task, threads, &[BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        let builds = field(&stderr, BUILD, "builds=");
        let kept = field(&stderr, KEEP, "kept=");
        assert_eq!(
            builds, 4,
            "{task} threads {threads}: one table per genome:\n{stderr}"
        );
        assert_eq!(
            kept + field(&stderr, KEEP, "fresh="),
            builds,
            "{task} threads {threads}: every build of the loop through the keep:\n{stderr}"
        );
        assert_eq!(kept, 3, "{task} threads {threads}: {stderr}");

        let (report, stderr) = blastn(task, threads, &[OFF, BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        assert!(!stderr.contains(KEEP), "{task}: {stderr}");
        let (report, stderr) = blastn(task, threads, &[DEFAULT_OFF, BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        assert!(!stderr.contains(KEEP), "{task}: {stderr}");
        // Set, the switch keeps the table with `LOSAT_X_DEFAULT=0` too (where the original loop
        // fills it: `LOSAT_X_PAIRPAR` is off).
        let (report, stderr) = blastn(task, threads, &[SWITCH, DEFAULT_OFF, BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
        assert_eq!(field(&stderr, KEEP, "kept="), kept, "{task}: {stderr}");

        let (report, stderr) = blastn(task, threads, &[BATCH, SHADOW]);
        assert!(
            report == reference,
            "{task} threads {threads} with the shadow"
        );
        let line = stderr
            .lines()
            .find(|line| line.starts_with(SHADOW_LINE))
            .unwrap_or_else(|| panic!("{task}: no shadow line:\n{stderr}"));
        assert!(line.ends_with(" diffs=0 identical"), "{line}");
        assert_eq!(field(&stderr, SHADOW_LINE, "tables="), kept, "{line}");
        assert!(field(&stderr, SHADOW_LINE, "cells=") > 0, "{line}");
        let (_, stderr) = blastn(task, threads, &[OFF, BATCH, SHADOW]);
        assert!(
            stderr.contains("[X_SHADOW] MBLUTKEEP tables=0 "),
            "{task}: nothing compared without the switch:\n{stderr}"
        );
    }

    // With two slots searching the batches after the first, the loop keeps no table: it searches
    // the first batch while the slots may still start, and none after.
    let (report, stderr) = blastn(task, 4, &[BATCH, TIMING, STATS, QBATCHPAR, SLOTS2]);
    assert!(report == reference, "{task} with slots");
    assert!(
        field(&stderr, QBATCHPAR_LINE, "precomputed=") > 0,
        "{task}: the slots searched no batch:\n{stderr}"
    );
    assert!(!stderr.contains(KEEP), "{task}: {stderr}");
}

#[test]
fn megablast_gives_the_same_report_on_a_kept_table() {
    check("megablast");
}

#[test]
fn blastn_gives_the_same_report_on_a_kept_table() {
    check("blastn");
}

// The prefault round of `LOSAT_X_MBLUTPREFAULT` runs on the fresh table only: at 8 threads every
// first chunk of these genomes is dense enough for it (306 kb, 262,144 positions per chunk against
// 131,072 for 4^12 cells), so it runs in every build without the switch and in the fresh ones with
// it.
#[test]
fn megablast_runs_the_prefault_round_on_fresh_tables_only() {
    let (reference, stderr) = blastn("megablast", 8, &[OFF, BATCH, TIMING]);
    let builds = field(&stderr, BUILD, "builds=");
    assert_eq!(field(&stderr, ONCE, "fills="), builds, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), builds, "{stderr}");
    let (report, stderr) = blastn("megablast", 8, &[BATCH, TIMING]);
    assert!(report == reference);
    let fresh = field(&stderr, KEEP, "fresh=");
    assert_eq!(fresh, 1, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), fresh, "{stderr}");
}
