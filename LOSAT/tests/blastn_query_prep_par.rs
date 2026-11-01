//! The batch preparation of blastn on the search pool (`LOSAT_X_QPREPPAR`,
//! `src/algorithm/blastn/blast_engine/x_query_prep_par.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses between
//! two LOSAT code paths with the same result: the BLASTNA encoding of a query batch and the
//! ungapped Karlin blocks of its contexts made one after the other, or as ordered parallel maps on
//! the search pool. So the property under test is that the report does not depend on it: for
//! megablast, blastn, dc-megablast and blastn-short (at `-evalue 10`) with several queries
//! (pieces of LC738874 with microsatellites, runs of `N`, lower case letters, a poly-A query, a
//! query shorter than a triplet, one of 20 letters and one of 300 `N` followed by 400 letters) at
//! 1 and 4 threads, the searches without the switch and with it give the same bytes. The same
//! holds with several query batches (`BATCH_SIZE`) and when a batch is split into query chunks
//! (`CHUNK_SIZE`). That the steps run on the pool's
//! threads is checked by the unit test `steps_run_on_the_pool_threads` in `x_query_prep_par.rs`;
//! here the sides of the switch's decision (`LOSAT_X_REACH`) show that the pool path ran at 4
//! threads and not at one. That those bytes equal NCBI's is the business of the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

mod reach_support;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_QPREPPAR", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the first record of a FASTA file.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A file name in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_qpreppar_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A query file removed when dropped: 24 pieces of LC738874 of 900 to 4,500 letters, every third
/// with a microsatellite, every fourth with a run of `N`s and lower case letters, then a poly-A
/// query (`polyA`), a query of two letters (`two`), a short one of 20 letters (`short`) and one of
/// 300 `N` followed by 400 letters (`mostlyN`): 28 queries.
struct Queries(PathBuf);

impl Queries {
    fn new() -> Self {
        let genome = letters("LC738874.fasta");
        let mut text = Vec::new();
        let mut put = |title: String, mut seq: Vec<u8>| {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            seq.push(b'\n');
            for line in seq.chunks(61) {
                text.extend_from_slice(line);
            }
            if !text.ends_with(b"\n") {
                text.push(b'\n');
            }
        };
        for i in 0..24usize {
            let start = i * 11_000 + 137 * i;
            let len = 900 + (i * 1_531) % 3_600;
            let mut seq = genome[start..start + len].to_vec();
            if i % 3 == 0 {
                let at = len / 3;
                let unit: &[u8] = [b"CA".as_slice(), b"GAA", b"T", b"ACGGT"][(i / 3) % 4];
                let run: Vec<u8> = unit.iter().copied().cycle().take(180).collect();
                seq.splice(at..at, run);
            }
            if i % 4 == 0 {
                let at = len / 2;
                seq.splice(at..at, std::iter::repeat_n(b'N', 25));
                for byte in &mut seq[at + 40..at + 140] {
                    *byte = byte.to_ascii_lowercase();
                }
            }
            put(format!("piece{i}"), seq);
        }
        put("polyA".into(), vec![b'A'; 700]);
        put("two".into(), b"AC".to_vec());
        put("short".into(), b"ACGTACGTACGTTTGACAGT".to_vec());
        let mut mostly_n = vec![b'N'; 300];
        mostly_n.extend_from_slice(&genome[5_000..5_400]);
        put("mostlyN".into(), mostly_n);
        let path = temp_path("queries.fa");
        std::fs::write(&path, text).expect("write the queries");
        Queries(path)
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// One search of the queries against LC738868, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables (and `BATCH_SIZE` / `CHUNK_SIZE`) that choose a code
/// path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as
/// in CI) needs that file name as soon as a search has more than one thread, and it selects
/// nothing. Returns the report and the standard error.
fn blastn(
    queries: &Queries,
    extra: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(&queries.0)
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For one task: the search without the switch at one thread, then off at 4 threads and on at 1
/// and 4 threads, all with the same bytes and no standard error; the switch's own path runs only at
/// 4 threads (the sides of its decision).
fn check_task(task: &str, extra: &[&str]) {
    let queries = Queries::new();
    let mut args = vec!["-task", task];
    args.extend_from_slice(extra);
    let (expected, expected_stderr) = blastn(&queries, &args, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [1, 4] {
        // At one thread the search without the switch is `expected` itself.
        if threads > 1 {
            let (report, stderr) = blastn(&queries, &args, threads, &[]);
            assert!(report == expected, "{task} threads {threads}: off differs");
            assert_eq!(stderr, "", "{task} threads {threads}");
        }
        let reach = reach_support::Reach::new("qpreppar");
        let (report, stderr) = blastn(&queries, &args, threads, &[SWITCH, reach.var()]);
        assert!(report == expected, "{task} threads {threads}: on differs");
        assert_eq!(stderr, "", "{task} threads {threads}");
        let sides = reach.sides("LOSAT_X_QPREPPAR");
        assert_eq!(sides.0, threads > 1, "{task} threads {threads}: {sides:?}");
    }
}

type Var = (&'static str, &'static str);

/// The searches of `env` (`BATCH_SIZE` and `CHUNK_SIZE` values), megablast: the report without
/// the switch at one thread, then off, and with the switch, at 4 threads, all the same bytes and no standard error. Returns whether the switch's own path ran (its decision's
/// new side).
fn check_env(env: &[Var]) -> bool {
    let queries = Queries::new();
    let args = ["-task", "megablast"];
    let with = |extra: &[Var]| -> Vec<Var> { env.iter().chain(extra).copied().collect() };
    let (expected, stderr) = blastn(&queries, &args, 1, &with(&[]));
    assert!(!expected.is_empty(), "{env:?}: no hits");
    assert_eq!(stderr, "", "{env:?}");
    let (report, stderr) = blastn(&queries, &args, 4, &with(&[]));
    assert!(report == expected, "{env:?}: off differs");
    assert_eq!(stderr, "", "{env:?}");
    let reach = reach_support::Reach::new("qpreppar-env");
    let mut on: Vec<(&str, &str)> = with(&[SWITCH]);
    on.push(reach.var());
    let (report, stderr) = blastn(&queries, &args, 4, &on);
    assert!(report == expected, "{env:?}: on differs");
    assert_eq!(stderr, "", "{env:?}");
    reach.sides("LOSAT_X_QPREPPAR").0
}

// Several batches (a batch size of 10,000 letters): each batch's steps run on the pool, and the
// report is the same.
#[test]
fn batches_prepared_in_several_batches() {
    assert!(check_env(&[("BATCH_SIZE", "10000")]));
}

// A batch that is split into query chunks (`CHUNK_SIZE` below the batch's letters): each part
// of a split query is searched as a batch of its own.
#[test]
fn split_query_chunks_prepared_on_the_pool() {
    assert!(check_env(&[
        ("BATCH_SIZE", "1000000"),
        ("CHUNK_SIZE", "7000")
    ],));
}

#[test]
fn megablast_batch_prepared_on_the_pool() {
    check_task("megablast", &[]);
}

#[test]
fn blastn_batch_prepared_on_the_pool() {
    check_task("blastn", &[]);
}

#[test]
fn dc_megablast_batch_prepared_on_the_pool() {
    check_task("dc-megablast", &[]);
}

// blastn-short (word size 7, the short-query scores); `check_task` asserts that it finds hits.
// `-evalue 10` instead of the task's 1000: a search writes 53,000 rows in 1.5 CPU s instead of
// 359,000 in 6-8 s. Every row still reads its context's ungapped block: blastn-short's gap costs
// (5/2 at 1/-3) are not in the table, so the gapped block of a context is made from its ungapped
// block. At 10 every query but `two` (no hits at either E-value) has hits on both strands.
#[test]
fn blastn_short_batch_prepared_on_the_pool() {
    check_task("blastn-short", &["-evalue", "10"]);
}

// Lower-case masking and another DUST window and level do not touch the encoding or the blocks.
#[test]
fn lowercase_masking_and_dust_parameters() {
    check_task("megablast", &["-dust", "10 32 2", "-lcase_masking"]);
}

// `-dust no` and other scores (other ungapped Karlin blocks): the same bytes with the switch, whose
// own path runs at 4 threads only. (There is no single-strand case: `-strand plus` and `minus`
// build their buffer serially, `build_strand_query_blastna_buffer`, whatever the switch says; its
// unit tests are in `blast_engine/run.rs`, the whole searches in `tests/blastn_strand.rs`.)
#[test]
fn dust_no_and_other_scores() {
    let queries = Queries::new();
    for args in [
        &["-task", "megablast", "-dust", "no"][..],
        &[
            "-task",
            "blastn",
            "-reward",
            "2",
            "-penalty",
            "-3",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
        ][..],
        &[
            "-task",
            "megablast",
            "-reward",
            "1",
            "-penalty",
            "-3",
            "-gapopen",
            "2",
            "-gapextend",
            "2",
        ][..],
    ] {
        let (expected, _) = blastn(&queries, args, 1, &[]);
        for threads in [1, 4] {
            let reach = reach_support::Reach::new("qpreppar-scores");
            let (report, _) = blastn(&queries, args, threads, &[SWITCH, reach.var()]);
            assert!(report == expected, "{args:?} threads {threads}: on differs");
            let sides = reach.sides("LOSAT_X_QPREPPAR");
            assert_eq!(
                sides.0,
                threads > 1,
                "{args:?} threads {threads}: {sides:?}"
            );
        }
    }
}

// With LOSAT_X_STATS the two steps print their summed wall times and the number of fallbacks at
// exit (`encode_s= karlin_s= fallbacks=0`), with or without the switch, when a batch was prepared.
#[test]
fn stats_line_is_printed() {
    let queries = Queries::new();
    let args = ["-task", "megablast"];
    let (expected, _) = blastn(&queries, &args, 1, &[]);
    for switches in [
        &[("LOSAT_X_STATS", "1")][..],
        &[("LOSAT_X_STATS", "1"), SWITCH][..],
    ] {
        let (report, stderr) = blastn(&queries, &args, 4, switches);
        assert!(report == expected, "{switches:?}: the report differs");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[X_STATS] QPREPPAR "))
            .unwrap_or_else(|| panic!("{switches:?}: no stats line:\n{stderr}"));
        let numbers: Vec<f64> = line["[X_STATS] QPREPPAR ".len()..]
            .split_whitespace()
            .map(|word| {
                word.split_once('=')
                    .and_then(|(_, value)| value.parse().ok())
                    .unwrap_or_else(|| panic!("{line}"))
            })
            .collect();
        assert_eq!(numbers.len(), 3, "{line}");
        assert_eq!(numbers[2], 0.0, "{line}");
    }
}
