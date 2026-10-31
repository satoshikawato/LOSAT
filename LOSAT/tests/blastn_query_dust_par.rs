//! The DUST masks of a query batch on the search pool (`LOSAT_X_QPREPPAR`,
//! `src/algorithm/blastn/x_query_dust_par.rs`) on whole searches. The masks had a switch of their
//! own, `LOSAT_X_QDUSTPAR`, until 2026-11-01, when it was folded into `LOSAT_X_QPREPPAR` (same
//! rule: a pool of two or more threads and two or more queries).
//!
//! The switch is in the strict set (on unless it is `0`). It chooses between
//! two LOSAT code paths with the same result: the queries of a batch masked one after the other,
//! or the per-query map split over the search pool and collected in query order. So the property under
//! test is that the report does not depend on it: for megablast, blastn and dc-megablast with
//! several queries (pieces of LC738874 with microsatellites, runs of `N`, a poly-A query and a
//! query shorter than a triplet) at 1 and 4 threads, the searches without the switch and with it
//! give the same bytes, and `-dust no` makes no masks (no `[X_STATS] QDUSTPAR` line). The
//! same holds with several query batches (`BATCH_SIZE`) and when a batch is split into query
//! chunks (`CHUNK_SIZE`, which reaches `chunk_query_masks`). That the masks are made on the
//! pool's threads is checked by the unit test `masks_are_made_on_the_pool_threads` in
//! `x_query_dust_par.rs`; here the `pool=` count of the `[X_STATS] QDUSTPAR` line (the DUST calls
//! whose masks were made on the pool) shows that the pool path ran at 4 threads and not at one
//! (QPREPPAR's reach marks also fire for the encoding and the Karlin blocks of the same searches,
//! so they cannot tell the masks apart). That those bytes equal NCBI's is the business of the
//! blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

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
        "losat_qdustpar_{}_{}_{tag}",
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

/// For one task: the search without the switch at one thread, then off and on at 1 and 4 threads,
/// all with the same bytes and no standard error but the stats lines; the masks are made on the
/// pool only at 4 threads (`pool=` of the stats line).
fn check_task(task: &str, extra: &[&str]) {
    let queries = Queries::new();
    let mut args = vec!["-task", task];
    args.extend_from_slice(extra);
    let (expected, expected_stderr) = blastn(&queries, &args, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [1, 4] {
        let (report, stderr) = blastn(&queries, &args, threads, &[]);
        assert!(report == expected, "{task} threads {threads}: off differs");
        assert_eq!(stderr, "", "{task} threads {threads}");
        let (report, stderr) = blastn(&queries, &args, threads, &[SWITCH, STATS]);
        assert!(report == expected, "{task} threads {threads}: on differs");
        assert_eq!(
            without_stats(&stderr),
            "",
            "{task} threads {threads}: {stderr}"
        );
        let pool = dust_pool_calls(&stderr);
        assert_eq!(pool > 0, threads > 1, "{task} threads {threads}: {stderr}");
    }
}

type Var = (&'static str, &'static str);

/// The searches of `env` (`BATCH_SIZE` and `CHUNK_SIZE` values), megablast: the report without
/// the switch at one thread, then off, and with the switch, at 4 threads, all the same bytes and no
/// standard error but the stats lines. Returns whether the masks were made on the pool (`pool=`).
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
    let (report, stderr) = blastn(&queries, &args, 4, &with(&[SWITCH, STATS]));
    assert!(report == expected, "{env:?}: on differs");
    assert_eq!(without_stats(&stderr), "", "{env:?}: {stderr}");
    dust_pool_calls(&stderr) > 0
}

const STATS: Var = ("LOSAT_X_STATS", "1");

/// Standard error without the `[X_STATS]` lines.
fn without_stats(stderr: &str) -> String {
    stderr
        .lines()
        .filter(|line| !line.starts_with("[X_STATS]"))
        .map(|line| format!("{line}\n"))
        .collect()
}

/// The `pool=` of the `[X_STATS] QDUSTPAR` line: the DUST calls whose masks were made on the
/// search pool (the switch's own path); 0 without the line.
fn dust_pool_calls(stderr: &str) -> u64 {
    stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] QDUSTPAR "))
        .and_then(|line| {
            line.split_whitespace()
                .find_map(|field| field.strip_prefix("pool="))
        })
        .map(|value| value.parse().expect("pool= is a count"))
        .unwrap_or(0)
}

// Several batches (a batch size of 10,000 letters): each batch's DUST masks are made on the pool,
// and the report is the same.
#[test]
fn queries_masked_in_several_batches() {
    assert!(check_env(&[("BATCH_SIZE", "10000")]));
}

// A batch that is split into query chunks (`CHUNK_SIZE` below the batch's letters): the parts of
// the split queries are masked by `chunk_query_masks`.
#[test]
fn split_query_chunks_masked_on_the_pool() {
    assert!(check_env(&[
        ("BATCH_SIZE", "1000000"),
        ("CHUNK_SIZE", "7000")
    ],));
}

#[test]
fn megablast_queries_masked_on_the_pool() {
    check_task("megablast", &[]);
}

#[test]
fn blastn_queries_masked_on_the_pool() {
    check_task("blastn", &[]);
}

#[test]
fn dc_megablast_queries_masked_on_the_pool() {
    check_task("dc-megablast", &[]);
}

// Another window and level, and the lower-case masking that merges with the DUST masks.
#[test]
fn dust_parameters_and_lowercase_masking() {
    check_task("megablast", &["-dust", "10 32 2", "-lcase_masking"]);
}

// `-dust no` makes no masks: the report is the same, and no DUST mask is made (QPREPPAR's other
// marks fire at 4 threads, so the proof is the absence of the `[X_STATS] QDUSTPAR` line, which a
// masked batch prints; the same search with DUST on prints it).
#[test]
fn dust_no_is_not_acted_on() {
    let queries = Queries::new();
    let args = ["-task", "megablast", "-dust", "no"];
    let (expected, _) = blastn(&queries, &args, 1, &[]);
    assert!(!expected.is_empty());
    for threads in [1, 4] {
        let (report, stderr) = blastn(&queries, &args, threads, &[SWITCH, ("LOSAT_X_STATS", "1")]);
        assert!(report == expected, "threads {threads}: on differs");
        assert!(!stderr.contains("QDUSTPAR"), "threads {threads}: {stderr}");
    }
    let (_, stderr) = blastn(
        &queries,
        &["-task", "megablast"],
        4,
        &[SWITCH, ("LOSAT_X_STATS", "1")],
    );
    assert!(stderr.contains("[X_STATS] QDUSTPAR "), "{stderr}");
}

// With LOSAT_TIMING (or LOSAT_X_STATS) the DUST step prints its summed wall time at exit, with or
// without the switch, when a batch was masked.
#[test]
fn timing_line_is_printed() {
    let queries = Queries::new();
    let args = ["-task", "megablast"];
    let (expected, _) = blastn(&queries, &args, 1, &[]);
    for switches in [
        &[("LOSAT_TIMING", "1")][..],
        &[("LOSAT_TIMING", "1"), SWITCH][..],
    ] {
        let (report, stderr) = blastn(&queries, &args, 4, switches);
        assert!(report == expected, "{switches:?}: the report differs");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] query_dust: "))
            .unwrap_or_else(|| panic!("{switches:?}: no timing line:\n{stderr}"));
        assert!(line.ends_with('s'), "{line}");
    }
    let (report, stderr) = blastn(&queries, &args, 4, &[("LOSAT_X_STATS", "1"), SWITCH]);
    assert!(report == expected, "the report differs with LOSAT_X_STATS");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] QDUSTPAR dust_s="))
        .unwrap_or_else(|| panic!("no stats line:\n{stderr}"));
    let mut fields = line["[X_STATS] QDUSTPAR ".len()..].split_whitespace();
    assert!(
        fields
            .next()
            .and_then(|field| field.strip_prefix("dust_s="))
            .is_some_and(|value| value.parse::<f64>().is_ok()),
        "{line}"
    );
    assert!(
        fields
            .next()
            .and_then(|field| field.strip_prefix("pool="))
            .is_some_and(|value| value.parse::<u64>().is_ok_and(|calls| calls > 0)),
        "{line}"
    );
    assert!(!stderr.contains("[TIMING] query_dust"), "{stderr}");
    // `-dust no` masks no batch: neither line is printed.
    let args = ["-task", "megablast", "-dust", "no"];
    let (_, stderr) = blastn(
        &queries,
        &args,
        4,
        &[("LOSAT_TIMING", "1"), ("LOSAT_X_STATS", "1"), SWITCH],
    );
    assert!(
        !stderr.contains("query_dust") && !stderr.contains("QDUSTPAR"),
        "{stderr}"
    );
}
