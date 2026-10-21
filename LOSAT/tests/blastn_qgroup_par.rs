//! The gapped stages of a blastn subject run per query on the search pool (`LOSAT_X_QGROUPPAR`,
//! `src/algorithm/blastn/blast_engine/x_qgroup_par.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses between two LOSAT code paths
//! with the same result: the preliminary gapped loop and the traceback of the subject run on the
//! owner thread over all queries in one order, or per query on the pool. So the property under
//! test is that the report does not depend on it: for megablast, blastn, dc-megablast and
//! blastn-short with several queries (pieces of LC738874 with microsatellites, runs of `N`, lower
//! case letters, a poly-A query and two very short ones) against LC738868 at 1 and 4 threads, the
//! searches with `LOSAT_X_QGROUPPAR=0`, by default and with the shadow (`LOSAT_X_QGROUPPARSHADOW`,
//! which also runs the owner loop and compares every HSP of both stages) give the same bytes. The
//! same holds in `LOSAT_X_QBATCHPAR` slots of two threads and for batches split into query chunks.
//! With `LOSAT_X_STATS` the counters show which form ran: per query at 4 threads, the owner loop
//! where one query holds most of the HSPs. That the bytes equal NCBI's is the business of the
//! blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_QGROUPPAR", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_QGROUPPARSHADOW", "1");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

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
        "losat_qgrouppar_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A query file removed when dropped.
struct Queries(PathBuf);

impl Queries {
    fn write(records: &[(String, Vec<u8>)]) -> Self {
        let mut text = Vec::new();
        for (title, seq) in records {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            for line in seq.chunks(61) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
            if seq.is_empty() {
                text.push(b'\n');
            }
        }
        let path = temp_path("queries.fa");
        std::fs::write(&path, text).expect("write the queries");
        Queries(path)
    }

    /// 24 pieces of LC738874 of 900 to 4,500 letters, every third with a microsatellite, every
    /// fourth with a run of `N`s and lower case letters, then a poly-A query, one of 20 letters
    /// and one of 300 `N` followed by 400 letters: 27 queries.
    fn pieces() -> Self {
        Self::pieces_of("LC738874.fasta", 0)
    }

    /// The same kind of 27 queries cut from the subject itself (LC738868, 3,000 letters further
    /// on, which keeps a repeat region of the subject from holding most of the HSPs), so that
    /// megablast and dc-megablast find HSPs for most of them.
    fn subject_pieces() -> Self {
        Self::pieces_of("LC738868.fasta", 3_000)
    }

    fn pieces_of(name: &str, shift: usize) -> Self {
        let genome = letters(name);
        let mut records = Vec::new();
        for i in 0..24usize {
            let start = i * 11_000 + 137 * i + shift;
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
            records.push((format!("piece{i}"), seq));
        }
        records.push(("polyA".into(), vec![b'A'; 700]));
        records.push(("short".into(), b"ACGTACGTACGTTTGACAGT".to_vec()));
        let mut mostly_n = vec![b'N'; 300];
        mostly_n.extend_from_slice(&genome[5_000..5_400]);
        records.push(("mostlyN".into(), mostly_n));
        Self::write(&records)
    }

    /// One query of 120,000 letters of LC738874 and two of 1,000: the long query holds most of the
    /// HSPs, so the owner loop stays.
    fn one_long() -> Self {
        let genome = letters("LC738874.fasta");
        Self::write(&[
            ("long".into(), genome[10_000..130_000].to_vec()),
            ("small1".into(), genome[150_000..151_000].to_vec()),
            ("small2".into(), genome[170_000..171_000].to_vec()),
        ])
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// One search of the queries against LC738868, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables (and `BATCH_SIZE` / `CHUNK_SIZE`) that choose a code
/// path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing. Returns the report and the standard error.
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

/// A number of the `[X_SHADOW] QGROUPPAR` or `[X_STATS] QGROUPPAR` line of `stderr`.
fn field(stderr: &str, line_start: &str, name: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(line_start))
        .unwrap_or_else(|| panic!("no {line_start} line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|word| word.strip_prefix(name))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

/// The shadow's line says `diffs=0 identical`; returns its `(prelim, traceback, hsps)`.
fn shadow_counts(stderr: &str) -> (u64, u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] QGROUPPAR "))
        .unwrap_or_else(|| panic!("no shadow line:\n{stderr}"));
    assert!(line.ends_with(" diffs=0 identical"), "{line}");
    let at = "[X_SHADOW] QGROUPPAR ";
    (
        field(stderr, at, "prelim="),
        field(stderr, at, "traceback="),
        field(stderr, at, "hsps="),
    )
}

/// For one task: the search with the switch off at one thread, then off, by default and with the
/// shadow at 1 and 4 threads, all with the same bytes; the default writes nothing to standard
/// error, the shadow one line, which counts stages compared only at 4 threads. Returns the shadow's
/// counts at 4 threads.
fn check_task(
    queries: &Queries,
    task: &str,
    extra: &[&str],
    env: &[(&str, &str)],
) -> (u64, u64, u64) {
    let mut args = vec!["-task", task];
    args.extend_from_slice(extra);
    let with = |more: &[(&'static str, &'static str)]| -> Vec<(&str, &str)> {
        env.iter().chain(more).copied().collect()
    };
    let (expected, expected_stderr) = blastn(queries, &args, 1, &with(&[OFF]));
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    let mut counts = (0, 0, 0);
    for threads in [1, 4] {
        let (report, stderr) = blastn(queries, &args, threads, &with(&[OFF]));
        assert!(report == expected, "{task} threads {threads}: off differs");
        assert_eq!(stderr, "", "{task} threads {threads}");
        let (report, stderr) = blastn(queries, &args, threads, &with(&[]));
        assert!(
            report == expected,
            "{task} threads {threads}: default differs"
        );
        assert_eq!(stderr, "", "{task} threads {threads}");
        let (report, stderr) = blastn(queries, &args, threads, &with(&[SHADOW]));
        assert!(
            report == expected,
            "{task} threads {threads}: shadow differs"
        );
        let shadow = shadow_counts(&stderr);
        if threads == 1 {
            assert_eq!(shadow, (0, 0, 0), "{task}");
        } else {
            counts = shadow;
        }
    }
    counts
}

#[test]
fn megablast_stages_per_query() {
    let (prelim, traceback, hsps) = check_task(&Queries::subject_pieces(), "megablast", &[], &[]);
    assert!(
        prelim > 0 && traceback > 0 && hsps > 0,
        "{prelim} {traceback} {hsps}"
    );
}

#[test]
fn blastn_stages_per_query() {
    let (prelim, traceback, hsps) = check_task(&Queries::pieces(), "blastn", &[], &[]);
    assert!(
        prelim > 0 && traceback > 0 && hsps > 0,
        "{prelim} {traceback} {hsps}"
    );
}

#[test]
fn dc_megablast_stages_per_query() {
    let (prelim, traceback, _) = check_task(&Queries::subject_pieces(), "dc-megablast", &[], &[]);
    assert!(prelim > 0 && traceback > 0, "{prelim} {traceback}");
}

#[test]
fn blastn_short_stages_per_query() {
    check_task(&Queries::pieces(), "blastn-short", &[], &[]);
}

// Lower-case masking and other scores (other cutoffs per context).
#[test]
fn masking_and_other_scores() {
    let queries = Queries::pieces();
    check_task(
        &queries,
        "megablast",
        &["-lcase_masking", "-dust", "no"],
        &[],
    );
    check_task(
        &queries,
        "blastn",
        &[
            "-reward",
            "2",
            "-penalty",
            "-3",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
        ],
        &[],
    );
}

// Several batches (`BATCH_SIZE`), searched in LOSAT_X_QBATCHPAR slots of two threads (each with a
// private pool, on which the stages run per query), and a batch split into query chunks.
#[test]
fn slots_and_split_batches() {
    let queries = Queries::pieces();
    let (_, traceback, _) = check_task(
        &queries,
        "blastn",
        &[],
        &[
            ("BATCH_SIZE", "20000"),
            ("LOSAT_X_QBATCHPAR", "1"),
            ("LOSAT_X_QBATCHPARSLOTS", "2"),
        ],
    );
    assert!(traceback > 0, "{traceback}");
    check_task(
        &queries,
        "megablast",
        &[],
        &[("BATCH_SIZE", "1000000"), ("CHUNK_SIZE", "7000")],
    );
}

// One query with most of the HSPs keeps the owner loop (counted as `owner=`), and the report is
// the same; the pieces run per query at 4 threads and on the owner at 1 thread (no pool).
#[test]
fn stats_show_which_form_ran() {
    let stats = "[X_STATS] QGROUPPAR ";
    let long = Queries::one_long();
    let args = ["-task", "megablast"];
    let (expected, _) = blastn(&long, &args, 1, &[OFF]);
    let (report, stderr) = blastn(&long, &args, 4, &[STATS]);
    assert!(report == expected, "one long query: the report differs");
    assert!(field(&stderr, stats, "owner=") > 0, "{stderr}");
    assert_eq!(field(&stderr, stats, "traceback="), 0, "{stderr}");

    let pieces = Queries::pieces();
    let args = ["-task", "blastn"];
    let (_, stderr) = blastn(&pieces, &args, 4, &[STATS]);
    assert!(field(&stderr, stats, "prelim=") > 0, "{stderr}");
    assert!(field(&stderr, stats, "traceback_queries=") > 1, "{stderr}");
    // At one thread no stage is planned, so the line is not printed.
    let (_, stderr) = blastn(&pieces, &args, 1, &[STATS]);
    assert!(!stderr.contains(stats), "{stderr}");
}
