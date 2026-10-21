//! blastp's two-hit test with `LOSAT_X_BPDEFER` (`src/algorithm/blastp/x_seed_fused.rs`) on whole
//! searches: the candidates of a lookup cell's AVX2 batches are called after the cell's last batch
//! of eight instead of at the end of each batch.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`) since 2026-10-17: on unless
//! `LOSAT_X_BPDEFER=0`, or `LOSAT_X_DEFAULT=0` and it is not set. It acts only where
//! `LOSAT_X_BPDIAG`'s AVX2 batch runs (the fused seed loop of `LOSAT_X_BPSEED`, both in the set, a
//! CPU with AVX2, not `LOSAT_X_BPDIAG_SIMD=0`). Each search states the default it runs with: the
//! reference loop (`LOSAT_X_DEFAULT=0`), the default (no variable: the deferred calls where the CPU
//! has AVX2), the default with the switch off (`LOSAT_X_BPDEFER=0`: the AVX2 batch with the calls
//! per batch), and `LOSAT_X_DEFAULT=0` with `LOSAT_X_BPSEED` and `LOSAT_X_BPDIAG` set, with and
//! without `LOSAT_X_BPDEFER` (the runs they gave before the switches joined the set). The tests
//! compare the reports with the reference and with the batches at one and four threads on a
//! generated pair (long lookup cells from repeats, related proteins for two-hit extensions) and
//! the fixture pair, check from the
//! `[TIMING] x_bpdiag: mode=<m> cells=<n> batched_hits=<h> candidates=<k>` lines that the deferred
//! path ran and tested the same cells, hits and candidates as the batches, and run the shadow of
//! `LOSAT_X_BPSEEDSHADOW` (the deferred loop against the reference for every subject).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
//! ```c
//!             diag_coord = (query_offset - subject_offset) & diag_mask;
//!             /* If the reset bit is set, an extension just happened. */
//!             if (diag_array[diag_coord].flag) {
//!                 /* If we've already extended past this hit, skip it. */
//!                 if ((Int4) (subject_offset + diag_offset) <
//!                     diag_array[diag_coord].last_hit) {
//!                     continue;
//!                 }
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The mode and counts of one `[TIMING] x_bpdiag` line.
#[derive(Clone, Debug)]
struct DiagCounts {
    mode: String,
    cells: u64,
    hits: u64,
    candidates: u64,
}

/// One blastp run of `query` against `subject` with the given variables as its only `LOSAT_*`
/// variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing; without the caller's value
/// a file of the tests is named. `BATCH_SIZE` is the test's own. The run must exit with 0.
/// Returns the standard output, the `[TIMING] x_bpdiag` lines' counts, and the standard error.
fn blastp(
    query: &Path,
    subject: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, Vec<DiagCounts>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpdefer-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .env("LOSAT_TIMING", "1")
        .env("BATCH_SIZE", "100000")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    let counts = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bpdiag: mode="))
        .map(|rest| {
            let (mode, rest) = rest.split_once(" cells=").expect("cells");
            let (cells, rest) = rest.split_once(" batched_hits=").expect("hits");
            let (hits, candidates) = rest.split_once(" candidates=").expect("candidates");
            DiagCounts {
                mode: mode.to_string(),
                cells: cells.parse().expect("a count of cells"),
                hits: hits.parse().expect("a count of hits"),
                candidates: candidates.parse().expect("a count of candidates"),
            }
        })
        .collect();
    (output.stdout, counts, stderr)
}

fn next(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 33) as u32
}

/// A generated query and subject set in a directory of the calling test (`test`; the tests run in
/// parallel, so they must not share files): random proteins, a protein with
/// long repeats (lookup cells far above three hits, so batches of eight and the overflow array),
/// and subjects that are copies of the queries with every sixth letter changed (two hits on a
/// diagonal, extensions with and without a reset of the diagonal), plus unrelated proteins.
fn generated_pair(test: &str) -> (PathBuf, PathBuf) {
    const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("blastp_diag_defer")
        .join(test);
    std::fs::create_dir_all(&dir).expect("test directory");
    let mut seed = 0xd1a9_u64;
    let mut protein = |n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| AMINO[next(&mut seed) as usize % AMINO.len()])
            .collect()
    };
    let mut queries: Vec<Vec<u8>> = (0..12).map(|i| protein(120 + 37 * i)).collect();
    let mut repeats = Vec::new();
    for _ in 0..12 {
        repeats.extend_from_slice(b"KLSEELLKKAEE");
    }
    repeats.extend(protein(60));
    for _ in 0..10 {
        repeats.extend_from_slice(b"GSGSPGTPGS");
    }
    queries.push(repeats);
    let mut subjects: Vec<Vec<u8>> = queries
        .iter()
        .map(|q| {
            q.iter()
                .enumerate()
                .map(|(i, &a)| if i % 6 == 5 { AMINO[(i * 7) % 20] } else { a })
                .collect()
        })
        .collect();
    subjects.extend((0..8).map(|i| protein(90 + 61 * i)));
    let write = |name: &str, seqs: &[Vec<u8>], tag: &str| {
        let mut text = String::new();
        for (i, seq) in seqs.iter().enumerate() {
            text.push_str(&format!(">{tag}{i}\n{}\n", String::from_utf8_lossy(seq)));
        }
        let path = dir.join(name);
        std::fs::write(&path, text).expect("write FASTA");
        path
    };
    (
        write("query.faa", &queries, "q"),
        write("subject.faa", &subjects, "s"),
    )
}

/// The reference loop: every switch of the strict set off.
const REFERENCE: &[(&str, &str)] = &[("LOSAT_X_DEFAULT", "0")];
/// The default: the strict set, with the deferred calls.
const DEFAULT: &[(&str, &str)] = &[];
/// The default without the deferred calls: the AVX2 batch with its calls per batch.
const DEFAULT_NO_DEFER: &[(&str, &str)] = &[("LOSAT_X_BPDEFER", "0")];
/// The batches of `LOSAT_X_BPDIAG` without the strict set (as before either joined it).
const BATCHES: &[(&str, &str)] = &[
    ("LOSAT_X_DEFAULT", "0"),
    ("LOSAT_X_BPSEED", "1"),
    ("LOSAT_X_BPDIAG", "1"),
];
/// The batches with the deferred calls without the strict set (the switch as it ran opt-in).
const DEFERRED: &[(&str, &str)] = &[
    ("LOSAT_X_DEFAULT", "0"),
    ("LOSAT_X_BPSEED", "1"),
    ("LOSAT_X_BPDIAG", "1"),
    ("LOSAT_X_BPDEFER", "1"),
];

fn avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// The one `[TIMING] x_bpdiag` line of a search.
fn one_line(counts: &[DiagCounts], what: &str) -> DiagCounts {
    assert_eq!(counts.len(), 1, "{what}: {counts:?}");
    counts[0].clone()
}

#[test]
fn deferred_calls_give_the_reports_and_counts_of_the_batches() {
    let (query, subject) = generated_pair("reports");
    let pairs = [
        (query, subject),
        (fixture("MeenMJNV.faa"), fixture("MejoMJNV.faa")),
    ];
    for (query, subject) in &pairs {
        for threads in [1, 4] {
            let what = format!("{} threads {threads}", query.display());
            let (reference, reference_counts, _) =
                blastp(query, subject, threads, REFERENCE, &["-outfmt", "6"]);
            assert!(!reference.is_empty(), "{what}: the pair has a report");
            assert!(
                reference_counts.is_empty(),
                "{what}: no line in the reference"
            );
            let (batched, batched_counts, _) =
                blastp(query, subject, threads, BATCHES, &["-outfmt", "6"]);
            assert!(batched == reference, "{what}: the batches differ");
            let b = one_line(&batched_counts, &what);
            // The runs with the deferred calls (`true`) or the calls per batch, by default and
            // as they ran opt-in.
            let runs = [(DEFAULT, true), (DEFAULT_NO_DEFER, false), (DEFERRED, true)];
            for (vars, deferred_calls) in runs {
                let what = format!("{what} {vars:?}");
                let (rows, counts, _) = blastp(query, subject, threads, vars, &["-outfmt", "6"]);
                assert!(rows == reference, "{what}: the report differs");
                let d = one_line(&counts, &what);
                if avx2() {
                    assert_eq!(b.mode, "avx2", "{what}");
                    let expected = if deferred_calls { "avx2-defer" } else { "avx2" };
                    assert_eq!(d.mode, expected, "{what}");
                    assert!(d.cells > 0 && d.candidates > 0, "{what}: {d:?}");
                } else {
                    // Without AVX2 there is no AVX2 batch to defer: the per-hit loop runs.
                    assert_eq!((b.mode.as_str(), d.mode.as_str()), ("per-hit", "per-hit"));
                }
                // The same cells, hits and candidates: only the time of the calls moved.
                assert_eq!(
                    (b.cells, b.hits, b.candidates),
                    (d.cells, d.hits, d.candidates),
                    "{what}: counts"
                );
            }
        }
    }
    // The pairwise report too, on the fixture pair, by default and as the switch ran opt-in.
    let (query, subject) = &pairs[1];
    let (reference, _, _) = blastp(query, subject, 1, REFERENCE, &["-outfmt", "0"]);
    for vars in [DEFAULT, DEFERRED] {
        let (report, _, _) = blastp(query, subject, 1, vars, &["-outfmt", "0"]);
        assert!(report == reference, "outfmt 0 {vars:?}: the report differs");
    }
}

#[test]
fn the_switch_acts_only_on_the_avx2_batches() {
    let (query, subject) = (fixture("MeenMJNV.faa"), fixture("MejoMJNV.faa"));
    let (reference, _, _) = blastp(&query, &subject, 1, REFERENCE, &["-outfmt", "6"]);
    // Without LOSAT_X_BPDIAG: the per-hit fused loop, no line. By default with the switch off,
    // and as the switch ran opt-in.
    let without_diag: [&[(&str, &str)]; 2] = [
        &[("LOSAT_X_BPDIAG", "0")],
        &[
            ("LOSAT_X_DEFAULT", "0"),
            ("LOSAT_X_BPSEED", "1"),
            ("LOSAT_X_BPDEFER", "1"),
        ],
    ];
    for vars in without_diag {
        let (rows, counts, _) = blastp(&query, &subject, 1, vars, &["-outfmt", "6"]);
        assert!(rows == reference, "without BPDIAG {vars:?}: rows differ");
        assert!(counts.is_empty(), "without BPDIAG {vars:?}: {counts:?}");
    }
    // With the scalar batch (LOSAT_X_BPDIAG_SIMD=0): the scalar batch, its calls per batch. By
    // default, and as the switch ran opt-in.
    let mut scalar = DEFERRED.to_vec();
    scalar.push(("LOSAT_X_BPDIAG_SIMD", "0"));
    let scalar_runs: [&[(&str, &str)]; 2] = [&[("LOSAT_X_BPDIAG_SIMD", "0")], &scalar];
    for vars in scalar_runs {
        let (rows, counts, _) = blastp(&query, &subject, 1, vars, &["-outfmt", "6"]);
        assert!(rows == reference, "scalar batch {vars:?}: rows differ");
        assert_eq!(one_line(&counts, "scalar batch").mode, "scalar");
    }
}

#[test]
fn the_shadow_compares_the_deferred_loop_with_the_reference() {
    let (query, subject) = generated_pair("shadow");
    let (reference, _, _) = blastp(&query, &subject, 1, REFERENCE, &["-outfmt", "6"]);
    // The shadow by default (the switch on), and as it ran opt-in.
    let shadows: [&[(&str, &str)]; 2] = [
        &[("LOSAT_X_BPSEEDSHADOW", "1")],
        &[
            ("LOSAT_X_DEFAULT", "0"),
            ("LOSAT_X_BPSEEDSHADOW", "1"),
            ("LOSAT_X_BPDIAG", "1"),
            ("LOSAT_X_BPDEFER", "1"),
        ],
    ];
    for vars in shadows {
        for threads in [1, 4] {
            let what = format!("{vars:?} threads {threads}");
            let (rows, counts, stderr) = blastp(&query, &subject, threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "{what}: rows differ");
            assert!(
                stderr
                    .lines()
                    .any(|line| line.starts_with("[X_SHADOW] BPSEED subjects=")
                        && line.contains(" identical")),
                "{what}: no shadow line: {stderr}"
            );
            let expected = if avx2() { "avx2-defer" } else { "per-hit" };
            assert_eq!(one_line(&counts, &what).mode, expected);
        }
    }
}
