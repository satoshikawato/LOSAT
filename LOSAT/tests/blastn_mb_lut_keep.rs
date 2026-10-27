//! The megablast hashtable kept from batch to batch of blastn's batch loop (`LOSAT_X_MBLUTKEEP`,
//! `src/algorithm/blastn/x_mb_lut_keep.rs`) on whole searches.
//!
//! The switch is in the strict set since 2026-10-19 (on unless it is `0`). It chooses between two
//! LOSAT code paths with the same tables: each batch allocates a zeroed table and frees it, or the
//! loop keeps the largest table and the next build of that size zeroes and fills it. So the
//! property under test is that the report does not depend on it, and that the kept path really runs
//! where the switch says: with `LOSAT_TIMING`, `[TIMING] x_mb_lut_keep` counts `kept=` builds (at
//! least one) and `fresh=` builds that add up to the search's builds (`[TIMING] x_mb_lookup_build:
//! builds=`), at 4 threads (zeroing and fill on the pool) and at one thread; the prefault round of
//! `LOSAT_X_MBLUTPREFAULT` runs only on the fresh tables (`prefault_rounds=` of `[TIMING]
//! x_mb_lookup_once` equals `fresh=` where every first chunk is dense enough, and the builds
//! without the switch); no line with
//! `LOSAT_X_MBLUTKEEP=0`, kept builds with `LOSAT_X_MBLUTONCE=0` and `LOSAT_X_MBFILL=0` (the
//! original fill loops); and no table is kept while `LOSAT_X_QBATCHPAR`'s slots may search
//! batches. That the bytes equal NCBI's is the business of the blastn fixtures.
//!
//! The query file holds four genomes of 292-306 kb (LC738868 306,008 letters, LC738869 300,002,
//! LC738870 294,144, LC738873 291,934) and `BATCH_SIZE=4000` puts each in a batch of its own (a
//! batch holds whole records): the megablast task builds four tables of 4^12 cells per search, so
//! the last three builds are on the kept table. blastn searches four pieces of 120,000 letters
//! instead, one from the middle of each genome: a blastn search of the genomes took 2.5 CPU s at
//! one thread and 5 at 4, most of it LC738868's alignment with itself; the pieces take 0.8 and 1.2.
//! Each piece still gets blastn's table of the genomes, 4^11 cells: 211,000-230,000 positions on
//! its two strands, above the 180,000 from which word size 11 takes 4^11 cells
//! (`LUT_WIDTH_11_THRESHOLD_10`), filled in two chunks at 4 threads. Four pieces of 5,000 letters
//! (9,478-9,920 positions, below 12,000) give blastn's two-stage tables (a lookup word of 8, built
//! in `ascending_cells` order), which go through the keep too. More searches: LC738868 followed by
//! four 100 kb pieces of LC738869 (tables of 4^11 cells: they use the start of the kept 4^12
//! table, no second table is allocated), LC738868 split into query chunks by `CHUNK_SIZE` (the
//! chunks' tables go through the keep), and 16 pieces of 20 kb of LC738868 in batches of one with
//! two slots and corrupted predictions (`LOSAT_X_QBATCHPARMISS`), where the slots stop partway and
//! the loop keeps its table from then on. The subject is LC738868.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const OFF: (&str, &str) = ("LOSAT_X_MBLUTKEEP", "0");
const SWITCH: (&str, &str) = ("LOSAT_X_MBLUTKEEP", "1");
const ONCE_OFF: (&str, &str) = ("LOSAT_X_MBLUTONCE", "0");
const FILL_OFF: (&str, &str) = ("LOSAT_X_MBFILL", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const BATCH: (&str, &str) = ("BATCH_SIZE", "4000");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
const SLOTS2: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");

const BUILD: &str = "[TIMING] x_mb_lookup_build: ";
const ONCE: &str = "[TIMING] x_mb_lookup_once: ";
const KEEP: &str = "[TIMING] x_mb_lut_keep: ";
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

/// Four pieces of `len` letters, one from the middle of each genome, in one query file.
fn middles(len: usize) -> PathBuf {
    let names = ["LC738868", "LC738869", "LC738870", "LC738873"];
    let seqs: Vec<String> = names
        .iter()
        .map(|name| {
            let genome = letters(&format!("{name}.fasta"));
            let start = genome.len() / 2 - len / 2;
            genome[start..start + len].to_string()
        })
        .collect();
    let records: Vec<(String, &str)> = names
        .iter()
        .zip(&seqs)
        .map(|(name, seq)| (format!("{name}_middle"), seq.as_str()))
        .collect();
    query_file(&format!("middles{len}"), &records)
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

/// For one task and four records: the same report with the switch off (at 4 threads with
/// `LOSAT_TIMING`, which prints no keep line), by default, and with `LOSAT_TIMING` by default, off
/// (at one thread) and on the original fill loops, at 4 threads and one; the kept builds counted
/// where they ran. No run sets `LOSAT_X_MBLUTKEEP=1` alone: unset and `1` are the same switch
/// state.
fn check(task: &str, query: &Path) {
    let (reference, stderr) = blastn(task, query, 4, &[OFF, BATCH, TIMING]);
    assert!(!reference.is_empty(), "{task}: no hits");
    assert!(!stderr.contains(KEEP), "{task}: {stderr}");
    let (report, stderr) = blastn(task, query, 4, &[BATCH]);
    assert!(report == reference, "{task}: the report differs");
    assert_eq!(stderr, "", "{task}: the switch writes nothing");

    for threads in [4, 1] {
        let (report, stderr) = blastn(task, query, threads, &[BATCH, TIMING]);
        assert!(report == reference, "{task} threads {threads}");
        let builds = field(&stderr, BUILD, "builds=");
        let kept = field(&stderr, KEEP, "kept=");
        assert_eq!(
            builds, 4,
            "{task} threads {threads}: one table per record:\n{stderr}"
        );
        assert_eq!(
            kept + field(&stderr, KEEP, "fresh="),
            builds,
            "{task} threads {threads}: every build of the loop through the keep:\n{stderr}"
        );
        assert_eq!(kept, 3, "{task} threads {threads}: {stderr}");

        // At 4 threads the search with the switch off is `reference` itself.
        if threads == 1 {
            let (report, stderr) = blastn(task, query, threads, &[OFF, BATCH, TIMING]);
            assert!(report == reference, "{task} threads {threads}");
            assert!(!stderr.contains(KEEP), "{task}: {stderr}");
        }
        // Set, the switch keeps the table on the original fill loops too.
        let (report, stderr) = blastn(
            task,
            query,
            threads,
            &[SWITCH, ONCE_OFF, FILL_OFF, BATCH, TIMING],
        );
        assert!(report == reference, "{task} threads {threads}");
        assert_eq!(field(&stderr, KEEP, "kept="), kept, "{task}: {stderr}");
    }

    // With two slots searching the batches after the first, the loop keeps no table: it searches
    // the first batch while the slots may still start, and none after.
    let (report, stderr) = blastn(task, query, 4, &[BATCH, TIMING, STATS, QBATCHPAR, SLOTS2]);
    assert!(report == reference, "{task} with slots");
    assert!(
        field(&stderr, QBATCHPAR_LINE, "precomputed=") > 0,
        "{task}: the slots searched no batch:\n{stderr}"
    );
    assert!(!stderr.contains(KEEP), "{task}: {stderr}");
}

#[test]
fn megablast_gives_the_same_report_on_a_kept_table() {
    check("megablast", genomes());
}

#[test]
fn blastn_gives_the_same_report_on_a_kept_table() {
    check("blastn", &middles(120_000));
}

// blastn's two-stage tables (a lookup word of 8 below the word of 11, built in `ascending_cells`
// order) on four pieces of 5,000 letters.
#[test]
fn blastn_two_stage_tables_go_through_the_keep() {
    check("blastn", &middles(5_000));
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

/// The same report with the switch off and on, and the keep's counts (`kept=`, `fresh=`) as
/// given.
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

// dc-megablast (discontiguous tables of 4^11 cells) and blastn (megablast tables of 4^11 cells) of
// four records in batches of one: LC738868 and LC738869 joined into one record, LC738870, LC738873
// and LC738871 joined, LC738875. A joined record has at least a quarter as many positions as
// cells and its table is laid out as CSR (`LOSAT_X_DCCSR`, `LOSAT_X_BNCSR`), a genome alone has
// fewer and gets chains. The first build (CSR) allocates its `offsets`, which the keep takes back
// when the batch drops the table; every build after it is on that table: the chain heads or
// hashtable (`kept=`) and the second CSR's `offsets` (`csr=`). The report equals the one with the
// switch off.
#[test]
fn csr_and_chain_builds_share_the_kept_table() {
    let seqs: Vec<String> = [
        "LC738868", "LC738869", "LC738870", "LC738873", "LC738871", "LC738875",
    ]
    .iter()
    .map(|name| letters(&format!("{name}.fasta")))
    .collect();
    let (first, second) = (seqs[0].clone() + &seqs[1], seqs[3].clone() + &seqs[4]);
    let records = [
        ("LC738868_LC738869".to_string(), first.as_str()),
        ("LC738870".to_string(), seqs[2].as_str()),
        ("LC738873_LC738871".to_string(), second.as_str()),
        ("LC738875".to_string(), seqs[5].as_str()),
    ];
    let query = query_file("csr_chains", &records);
    for task in ["dc-megablast", "blastn"] {
        let stderr = check_counts(task, &query, &[BATCH], 2, 1);
        assert_eq!(field(&stderr, KEEP, "csr="), 1, "{task}:\n{stderr}");
    }
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
