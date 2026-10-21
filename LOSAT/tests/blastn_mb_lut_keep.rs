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
//! The query file holds four genomes of 292-306 kb (LC738868 306,008 letters, LC738869 300,002,
//! LC738870 294,144, LC738873 291,934) and `BATCH_SIZE=4000` puts each in a batch of its own (a
//! batch holds whole records): four tables of 4^12 cells per search, so the last three builds are
//! on the kept table. More searches: LC738868 followed by four 100 kb pieces of LC738869 (tables of
//! 4^11 cells: they use the start of the kept 4^12 table, no second table is allocated), LC738868
//! split into query chunks by `CHUNK_SIZE` (the chunks' tables go through the keep), and 16 pieces
//! of 20 kb of LC738868 in batches of one with two slots and corrupted predictions
//! (`LOSAT_X_QBATCHPARMISS`), where the slots stop partway and the loop keeps its table from then
//! on. The subject is LC738868.

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

/// The letters of a one-record FASTA fixture.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .expect("read fixture")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// A query file of these records in the temporary directory, named after this test process.
fn query_file(name: &str, records: &[(String, &str)]) -> PathBuf {
    let text: String = records
        .iter()
        .map(|(title, seq)| format!(">{title}\n{seq}\n"))
        .collect();
    let path = std::env::temp_dir().join(format!(
        "losat_mblutkeep_{name}_{}.fasta",
        std::process::id()
    ));
    std::fs::write(&path, text).expect("write the query file");
    path
}

/// The four genomes in one query file (written once per test binary).
fn genomes() -> &'static Path {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        let names = ["LC738868", "LC738869", "LC738870", "LC738873"];
        let seqs: Vec<String> = names
            .iter()
            .map(|name| letters(&format!("{name}.fasta")))
            .collect();
        let records: Vec<(String, &str)> = names
            .iter()
            .zip(&seqs)
            .map(|(name, seq)| (name.to_string(), seq.as_str()))
            .collect();
        query_file("genomes", &records)
    })
}

/// One search, `-outfmt 6` on standard output, with the given variables as its only `LOSAT_*`
/// variables (and batch-size variables). `LOSAT_BLASTX_WORKER_LOG` stays, or is set to a file of
/// this test when the environment has none: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing. Returns the report and the standard error.
fn blastn(task: &str, query: &Path, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
        .arg(query)
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
    let (reference, stderr) = blastn(task, genomes(), 4, &[OFF, BATCH]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert_eq!(stderr, "", "{task}");
    for switches in [
        &[BATCH][..],
        &[SWITCH, BATCH],
        &[SWITCH, DEFAULT_OFF, BATCH],
    ] {
        let (report, stderr) = blastn(task, genomes(), 4, switches);
        assert!(
            report == reference,
            "{task} {switches:?}: the report differs"
        );
        assert_eq!(stderr, "", "{task} {switches:?}: the switch writes nothing");
    }

    for threads in [4, 1] {
        let (report, stderr) = blastn(task, genomes(), threads, &[BATCH, TIMING]);
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

        let (report, stderr) = blastn(task, genomes(), threads, &[OFF, BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        assert!(!stderr.contains(KEEP), "{task}: {stderr}");
        let (report, stderr) = blastn(task, genomes(), threads, &[DEFAULT_OFF, BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        assert!(!stderr.contains(KEEP), "{task}: {stderr}");
        // Set, the switch keeps the table with `LOSAT_X_DEFAULT=0` too (where the original loop
        // fills it: `LOSAT_X_PAIRPAR` is off).
        let (report, stderr) = blastn(
            task,
            genomes(),
            threads,
            &[SWITCH, DEFAULT_OFF, BATCH, TIMING],
        );
        assert!(report == reference, "{task} threads {threads}");
        assert_eq!(field(&stderr, BUILD, "parallel="), 0, "{task}: {stderr}");
        assert_eq!(field(&stderr, KEEP, "kept="), kept, "{task}: {stderr}");

        let (report, stderr) = blastn(task, genomes(), threads, &[BATCH, SHADOW]);
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
        let (_, stderr) = blastn(task, genomes(), threads, &[OFF, BATCH, SHADOW]);
        assert!(
            stderr.contains("[X_SHADOW] MBLUTKEEP tables=0 "),
            "{task}: nothing compared without the switch:\n{stderr}"
        );
    }

    // With two slots searching the batches after the first, the loop keeps no table: it searches
    // the first batch while the slots may still start, and none after.
    let (report, stderr) = blastn(
        task,
        genomes(),
        4,
        &[BATCH, TIMING, STATS, QBATCHPAR, SLOTS2],
    );
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
// first chunk of these genomes is dense enough for it (chunks of 262,144 positions, which each
// genome fills with its two strands, against 131,072 for 4^12 cells), so it runs in every build
// without the switch and in the fresh one with it.
#[test]
fn megablast_runs_the_prefault_round_on_fresh_tables_only() {
    let (reference, stderr) = blastn("megablast", genomes(), 8, &[OFF, BATCH, TIMING]);
    let builds = field(&stderr, BUILD, "builds=");
    assert_eq!(field(&stderr, ONCE, "fills="), builds, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), builds, "{stderr}");
    let (report, stderr) = blastn("megablast", genomes(), 8, &[BATCH, TIMING]);
    assert!(report == reference);
    let fresh = field(&stderr, KEEP, "fresh=");
    assert_eq!(fresh, 1, "{stderr}");
    assert_eq!(field(&stderr, ONCE, "prefault_rounds="), fresh, "{stderr}");
}

/// The same report with the switch off and on, the shadow's tables identical, and the keep's
/// counts (`kept=`, `fresh=`) as given.
fn check_counts(task: &str, query: &Path, env: &[(&str, &str)], kept: u64, fresh: u64) -> String {
    let (reference, _) = blastn(task, query, 4, &[env, &[OFF]].concat());
    assert!(!reference.is_empty(), "{task}: no hits");
    let (report, stderr) = blastn(task, query, 4, &[env, &[TIMING, STATS]].concat());
    assert!(report == reference, "{task} {env:?}: the report differs");
    assert_eq!(
        field(&stderr, KEEP, "kept="),
        kept,
        "{task} {env:?}:\n{stderr}"
    );
    assert_eq!(
        field(&stderr, KEEP, "fresh="),
        fresh,
        "{task} {env:?}:\n{stderr}"
    );
    let (report, shadow) = blastn(task, query, 4, &[env, &[SHADOW]].concat());
    assert!(report == reference, "{task} {env:?} with the shadow");
    assert_eq!(field(&shadow, SHADOW_LINE, "tables="), kept, "{shadow}");
    assert!(shadow.contains(" diffs=0 identical"), "{shadow}");
    stderr
}

// Megablast of LC738868 (a table of 4^12 cells) and then four 100 kb pieces of LC738869 (4^11
// cells each): the pieces' tables are the start of the kept table, so the search allocates one
// table.
#[test]
fn shorter_tables_are_the_start_of_the_kept_table() {
    let (long, short) = (letters("LC738868.fasta"), letters("LC738869.fasta"));
    let mut records = vec![("LC738868".to_string(), long.as_str())];
    for i in 0..4 {
        records.push((
            format!("LC738869_{i}"),
            &short[i * 65_000..i * 65_000 + 100_000],
        ));
    }
    let query = query_file("mixed", &records);
    check_counts("megablast", &query, &[BATCH], 4, 1);
}

// LC738868 against itself split into query chunks of 100,000 letters (`CHUNK_SIZE`): the batch
// builds no table, each of its three chunks builds one, the last two on the kept table.
#[test]
fn the_chunks_of_a_split_batch_use_the_kept_table() {
    let query = fasta("LC738868.fasta");
    check_counts("megablast", &query, &[("CHUNK_SIZE", "100000")], 2, 1);
}

// 16 pieces of 20 kb of LC738868 in batches of one, two slots, every second prediction corrupted:
// the slots stop partway (`off=` other than `none`, some batches taken from them), and the loop's
// own searches after that keep its table.
#[test]
fn the_loop_keeps_its_table_once_the_slots_stop() {
    let genome = letters("LC738868.fasta");
    let records: Vec<(String, &str)> = (0..16)
        .map(|i| {
            (
                format!("LC738868_{i}"),
                &genome[i * 19_000..i * 19_000 + 20_000],
            )
        })
        .collect();
    let query = query_file("pieces", &records);
    let env = [
        ("BATCH_SIZE", "1"),
        QBATCHPAR,
        SLOTS2,
        ("LOSAT_X_QBATCHPARMISS", "2"),
    ];
    let (reference, _) = blastn("megablast", &query, 4, &[OFF, ("BATCH_SIZE", "1")]);
    let (report, stderr) = blastn(
        "megablast",
        &query,
        4,
        &[&env[..], &[TIMING, STATS]].concat(),
    );
    assert!(report == reference, "the report differs");
    assert!(
        field(&stderr, QBATCHPAR_LINE, "precomputed=") > 0,
        "{stderr}"
    );
    let line = stderr
        .lines()
        .find(|line| line.starts_with(QBATCHPAR_LINE))
        .unwrap();
    assert!(
        !line.ends_with(" off=none"),
        "the slots never stopped: {line}"
    );
    assert!(field(&stderr, KEEP, "kept=") > 0, "{stderr}");
}
