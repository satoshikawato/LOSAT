//! blastp's two-hit test in batches with `LOSAT_X_BPDIAG` (`src/algorithm/blastp/x_seed_fused.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`) since 2026-10-15: on unless
//! `LOSAT_X_BPDIAG=0`, or `LOSAT_X_DEFAULT=0` and it is not set. It acts in the fused seed loop of
//! `LOSAT_X_BPSEED` (also in the set). With it on, the hits of one lookup cell are tested in
//! batches of eight in AVX2 registers where the CPU has AVX2; without AVX2 the per-hit loop runs,
//! and `LOSAT_X_BPDIAG_SIMD=0` forces the scalar batch. These tests compare the reports with the
//! reference loop (`LOSAT_X_DEFAULT=0`) and the per-hit fused loop at one and four threads, on a
//! generated pair (long lookup cells from repeats, related proteins for two-hit extensions) and on
//! the fixture pair, run the shadow of `LOSAT_X_BPSEEDSHADOW` (the batched loop against the
//! reference for every subject), and check from the
//! `[TIMING] x_bpdiag: mode=<m> cells=<n> batched_hits=<h> candidates=<k>` lines which path ran and
//! that cells were tested in batches, not only that the report is equal. Each search states the
//! default it runs with: no variable (the default), `=0`, and `LOSAT_X_DEFAULT=0` with the
//! variables set (the run they gave before the switch joined the set).
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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
            std::env::temp_dir().join("losat-bpdiag-worker.log"),
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

/// Checks one search's `[TIMING] x_bpdiag` line: `LOSAT_X_BPDIAG_SIMD=0` runs the scalar batch,
/// otherwise the AVX2 batch where this CPU has AVX2 and the per-hit loop (nothing counted) where
/// it has not. A batched run has tested cells and passed candidates on.
fn check_counts(counts: &[DiagCounts], simd: &str, what: &str) {
    assert_eq!(counts.len(), 1, "{what}: {counts:?}");
    let c = &counts[0];
    #[cfg(target_arch = "x86_64")]
    let avx2 = std::arch::is_x86_feature_detected!("avx2");
    #[cfg(not(target_arch = "x86_64"))]
    let avx2 = false;
    let expected = if simd == "0" {
        "scalar"
    } else if avx2 {
        "avx2"
    } else {
        "per-hit"
    };
    assert_eq!(c.mode, expected, "{what}: {c:?}");
    if expected == "per-hit" {
        assert!(
            c.cells == 0 && c.hits == 0 && c.candidates == 0,
            "{what}: {c:?}"
        );
    } else {
        assert!(
            c.cells > 0 && c.hits >= c.cells && c.candidates > 0 && c.candidates <= c.hits,
            "{what}: {c:?}"
        );
    }
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
        .join("blastp_diag_batch")
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

#[test]
fn batched_cells_give_the_reports_of_the_per_hit_loop() {
    let (query, subject) = generated_pair("per_hit_loop");
    let pairs = [
        (query, subject),
        (fixture("MeenMJNV.faa"), fixture("MejoMJNV.faa")),
    ];
    // The fused per-hit loop: the switch off in the default, and the run of `LOSAT_X_BPSEED`
    // alone before either joined the set.
    let per_hit: [&[(&str, &str)]; 2] = [
        &[("LOSAT_X_BPDIAG", "0")],
        &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_BPSEED", "1")],
    ];
    // The batches, with the `LOSAT_X_BPDIAG_SIMD` value each runs with: the default (the CPU's
    // path, the AVX2 batch or the per-hit loop without AVX2) and its scalar batch, and the run the
    // switch gave before it joined the set.
    let batched: [(&[(&str, &str)], &str); 4] = [
        (&[], "1"),
        (&[("LOSAT_X_BPDIAG_SIMD", "0")], "0"),
        (
            &[
                ("LOSAT_X_DEFAULT", "0"),
                ("LOSAT_X_BPSEED", "1"),
                ("LOSAT_X_BPDIAG", "1"),
            ],
            "1",
        ),
        (
            &[
                ("LOSAT_X_DEFAULT", "0"),
                ("LOSAT_X_BPSEED", "1"),
                ("LOSAT_X_BPDIAG", "1"),
                ("LOSAT_X_BPDIAG_SIMD", "0"),
            ],
            "0",
        ),
    ];
    for (query, subject) in &pairs {
        for threads in [1, 4] {
            let (reference, reference_counts, _) =
                blastp(query, subject, threads, REFERENCE, &["-outfmt", "6"]);
            assert!(!reference.is_empty(), "the pair has a report");
            assert!(reference_counts.is_empty(), "no line with the switch off");
            for vars in per_hit {
                let (fused, fused_counts, _) =
                    blastp(query, subject, threads, vars, &["-outfmt", "6"]);
                assert!(
                    fused == reference,
                    "{vars:?} threads {threads}: the fused loop differs"
                );
                assert!(
                    fused_counts.is_empty(),
                    "{vars:?}: no line with the switch off"
                );
            }
            for (vars, simd) in batched {
                let (rows, counts, _) = blastp(query, subject, threads, vars, &["-outfmt", "6"]);
                let what = format!("{} threads {threads} {vars:?}", query.display());
                assert!(rows == reference, "{what}: rows differ");
                check_counts(&counts, simd, &what);
            }
            // The switch on without the fused loop: the reference loop, nothing tested in batches.
            let (rows, counts, _) = blastp(
                query,
                subject,
                threads,
                &[("LOSAT_X_BPSEED", "0")],
                &["-outfmt", "6"],
            );
            assert!(
                rows == reference,
                "LOSAT_X_BPSEED=0 threads {threads}: rows differ"
            );
            assert!(
                counts.len() == 1 && counts[0].cells == 0 && counts[0].hits == 0,
                "LOSAT_X_BPSEED=0 threads {threads}: {counts:?}"
            );
        }
    }
    // The pairwise report too, on the fixture pair: the default against the reference.
    let (query, subject) = &pairs[1];
    let (reference, _, _) = blastp(query, subject, 1, REFERENCE, &["-outfmt", "0"]);
    let (report, counts, _) = blastp(query, subject, 1, &[], &["-outfmt", "0"]);
    assert!(report == reference, "outfmt 0: the report differs");
    check_counts(&counts, "1", "outfmt 0");
}

#[test]
fn the_shadow_compares_the_batched_loop_with_the_reference() {
    let (query, subject) = generated_pair("shadow");
    let (reference, _, _) = blastp(&query, &subject, 1, REFERENCE, &["-outfmt", "6"]);
    // The shadow with the switch by default, and as it ran before the switch joined the set.
    let defaults: [&[(&str, &str)]; 2] =
        [&[], &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_BPDIAG", "1")]];
    for default in defaults {
        for simd in ["1", "0"] {
            let mut vars = default.to_vec();
            vars.extend([("LOSAT_X_BPSEEDSHADOW", "1"), ("LOSAT_X_BPDIAG_SIMD", simd)]);
            let what = format!("shadow {vars:?}");
            let (rows, counts, stderr) = blastp(&query, &subject, 1, &vars, &["-outfmt", "6"]);
            assert!(rows == reference, "{what}: rows differ");
            assert!(
                stderr
                    .lines()
                    .any(|line| line.starts_with("[X_SHADOW] BPSEED subjects=")
                        && line.contains(" identical")),
                "{what}: no shadow line: {stderr}"
            );
            check_counts(&counts, simd, &what);
        }
    }
}
