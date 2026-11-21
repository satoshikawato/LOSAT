//! blastp's two-hit test with `LOSAT_X_BPDEFER` and `LOSAT_X_BPMASKED`
//! (`src/algorithm/blastp/x_seed_fused.rs`) on whole searches: the candidates of a lookup cell's
//! AVX2 batches are called after the cell's last batch of eight instead of at the end of each batch
//! (`LOSAT_X_BPDEFER`), and with `LOSAT_X_BPMASKED` the batches take their entry coordinates from
//! scalar loads of the hits and a cell's last batch of fewer than eight hits is one masked AVX2
//! batch whose candidates join the buffer (`x_cell_avx2_defer_masked`, the masked batch of
//! `LOSAT_X_BXDIAG` in blastx).
//!
//! Both switches are in the strict set (`src/utils/xswitch.rs`), `LOSAT_X_BPDEFER` since
//! 2026-10-17 and `LOSAT_X_BPMASKED` since 2026-10-24: on unless `=0`. `LOSAT_X_BPDEFER` acts only
//! where `LOSAT_X_BPDIAG`'s AVX2 batch runs (the fused seed loop, a CPU with AVX2, not
//! `LOSAT_X_BPDIAG_SIMD=0`), and `LOSAT_X_BPMASKED` only where
//! `LOSAT_X_BPDEFER`'s deferred calls run. Each search states the default it runs with: the
//! per-hit fused loop (`LOSAT_X_BPDIAG=0`), the default (no variable: the masked batch
//! where the CPU has AVX2), the default with `LOSAT_X_BPMASKED=0` (the deferred calls), and the
//! default with `LOSAT_X_BPDEFER=0` (the AVX2 batch with the calls per batch). The tests
//! compare the reports with the per-hit loop and with the batches on a generated pair (long lookup
//! cells from repeats, related proteins for two-hit extensions) at one and four threads and on the
//! fixture pair at four, check from the
//! `[TIMING] x_bpdiag: mode=<m> cells=<n> batched_hits=<h> candidates=<k>` lines that the masked
//! and the deferred paths ran and tested the same cells, hits and candidates as the batches.
//!
//! `LOSAT_X_BPCELL128` (strict set, class platform) runs the batches in 128-bit vectors where the
//! AVX2 ones do not run (the WebAssembly build, x86-64 with SSSE3, SSE4.1 and POPCNT but no AVX2;
//! `=2` on any x86-64 CPU with them): the calls at once, a cell's last hits one after the other,
//! on the tables `LOSAT_X_BPMASKEDMIN` gives the masked batch (the others keep the per-hit loop).
//! `the_128_bit_batches_give_the_reports_of_the_reference` compares it with `=2`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552

use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::{worker_log_if_unset, TempFile};

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
) -> (Vec<u8>, Vec<DiagCounts>, String, Option<TempFile>) {
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
    let log = worker_log_if_unset(&mut command, "blastp_diag_defer_worker");
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
    (output.stdout, counts, stderr, log)
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

/// The reference: the per-hit fused loop, without the batches and without their line.
const REFERENCE: &[(&str, &str)] = &[("LOSAT_X_BPDIAG", "0")];
/// The default: the strict set, with the deferred calls and the masked batch.
const DEFAULT: &[(&str, &str)] = &[];
/// The default without the masked batch: the deferred calls of `LOSAT_X_BPDEFER`.
const DEFERRED: &[(&str, &str)] = &[("LOSAT_X_BPMASKED", "0")];
/// The default without the deferred calls: the AVX2 batch with its calls per batch.
const BATCHES: &[(&str, &str)] = &[("LOSAT_X_BPDEFER", "0")];

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

/// Whether this CPU has what `LOSAT_X_BPCELL128`'s 128-bit batches need (SSSE3, SSE4.1 and
/// POPCNT): they run by default without AVX2, and with `=2` in place of the AVX2 batches.
fn cell128() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("sse4.1")
            && std::arch::is_x86_feature_detected!("popcnt")
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

/// The (thread, table) choices of `LOSAT_X_BPMASKEDMIN`'s rule from a search's
/// `[TIMING] x_bpmasked: min_hits_per_cell=<n> tables_masked=<a> tables_deferred=<b> ...` lines:
/// the tables that took the masked batch and those it gave the deferred calls (`None` without a
/// line: the search did not choose the masked batch).
fn masked_tables(stderr: &str) -> Option<(u64, u64)> {
    let mut found = None;
    for line in stderr.lines() {
        let Some(rest) = line.strip_prefix("[TIMING] x_bpmasked: ") else {
            continue;
        };
        let field = |name: &str| -> u64 {
            rest.split(' ')
                .find_map(|f| f.strip_prefix(name))
                .expect("a field of the line")
                .parse()
                .expect("a count of tables")
        };
        let (masked, deferred) = found.unwrap_or((0, 0));
        found = Some((
            masked + field("tables_masked="),
            deferred + field("tables_deferred="),
        ));
    }
    found
}

#[test]
fn deferred_and_masked_calls_give_the_reports_and_counts_of_the_batches() {
    let (query, subject) = generated_pair("reports");
    // The generated pair (0.01 CPU s a search) at one and four threads with tabular rows; the
    // fixture pair (7-12 CPU s a search) at four threads with the pairwise report (outfmt 0), which
    // shows the alignments too and, with 124 subjects (under the 250 alignments shown), every HSP
    // of the rows. `tests/blastp_diag_batch.rs` relies on this test for the fixture pair's run of
    // `BATCHES` and its pairwise reports of the per-hit loop and the default.
    let pairs = [
        (query, subject, &[1, 4][..], "6"),
        (
            fixture("MeenMJNV.faa"),
            fixture("MejoMJNV.faa"),
            &[4][..],
            "0",
        ),
    ];
    for (query, subject, thread_counts, outfmt) in &pairs {
        let extra = ["-outfmt", *outfmt];
        for &threads in *thread_counts {
            let what = format!("{} threads {threads}", query.display());
            let (reference, reference_counts, _, _log) =
                blastp(query, subject, threads, REFERENCE, &extra);
            // Rows, or a pairwise report with a subject's alignments (it has a header and says
            // "No hits found" for a query without them).
            let has_hits = if *outfmt == "0" {
                reference.windows(2).any(|pair| pair == b"\n>")
            } else {
                !reference.is_empty()
            };
            assert!(has_hits, "{what}: the pair has a report");
            assert!(
                reference_counts.is_empty(),
                "{what}: no line in the reference"
            );
            let (batched, batched_counts, _, _log) =
                blastp(query, subject, threads, BATCHES, &extra);
            assert!(batched == reference, "{what}: the batches differ");
            let b = one_line(&batched_counts, &what);
            // The run with the deferred calls (LOSAT_X_BPMASKED=0).
            let what_d = format!("{what} {DEFERRED:?}");
            let (report, counts, _, _log) = blastp(query, subject, threads, DEFERRED, &extra);
            assert!(report == reference, "{what_d}: the report differs");
            let d = one_line(&counts, &what_d);
            // The run with the masked batch, by default.
            let what = format!("{what} {DEFAULT:?}");
            let (report, counts, stderr, _log) = blastp(query, subject, threads, DEFAULT, &extra);
            assert!(report == reference, "{what}: the report differs");
            let m = one_line(&counts, &what);
            if avx2() {
                assert_eq!(b.mode, "avx2", "{what}");
                assert_eq!(d.mode, "avx2-defer", "{what_d}");
                assert_eq!(m.mode, "avx2-defer-masked", "{what}");
                assert!(m.cells > 0 && m.candidates > 0, "{what}: {m:?}");
                // The query sets have enough hits per cell for the masked batch on every table.
                let tables = masked_tables(&stderr);
                assert!(
                    matches!(tables, Some((n, 0)) if n > 0),
                    "{what}: {tables:?}"
                );
            } else {
                // Without AVX2 there is no AVX2 batch to defer: LOSAT_X_BPCELL128's batches run
                // where the CPU has their instructions, otherwise the per-hit loop.
                let expected = if cell128() { "w128" } else { "per-hit" };
                let modes = (b.mode.as_str(), d.mode.as_str(), m.mode.as_str());
                assert_eq!(modes, (expected, expected, expected));
            }
            // The same cells, hits and candidates: only the time of the calls and the form of the
            // batches changed.
            for (c, what) in [(&d, &what_d), (&m, &what)] {
                assert_eq!(
                    (b.cells, b.hits, b.candidates),
                    (c.cells, c.hits, c.candidates),
                    "{what}: counts"
                );
            }
        }
    }
}

#[test]
fn the_switch_acts_only_on_the_avx2_batches() {
    // These runs' cell tests (the per-hit loop, the scalar batch) are chosen once per run from the
    // variables and the CPU (`x_diag_mode`), not from the input: the generated pair.
    let (query, subject) = generated_pair("only_avx2");
    let (reference, _, _, _log) = blastp(&query, &subject, 1, REFERENCE, &["-outfmt", "6"]);
    assert!(!reference.is_empty(), "the pair has a report");
    // Without LOSAT_X_BPDIAG (the reference): the per-hit fused loop, no line.
    let (_, counts, _, _log) = blastp(&query, &subject, 1, REFERENCE, &["-outfmt", "6"]);
    assert!(counts.is_empty(), "without BPDIAG: {counts:?}");
    // With the scalar batch (LOSAT_X_BPDIAG_SIMD=0): the scalar batch, its calls per batch.
    let vars: &[(&str, &str)] = &[("LOSAT_X_BPDIAG_SIMD", "0")];
    let (rows, counts, _, _log) = blastp(&query, &subject, 1, vars, &["-outfmt", "6"]);
    assert!(rows == reference, "scalar batch {vars:?}: rows differ");
    assert_eq!(one_line(&counts, "scalar batch").mode, "scalar");
}

/// `LOSAT_X_BPMASKEDMIN`'s rule on whole searches: a short query (about one hit per used lookup
/// cell) keeps the deferred calls by default and takes the masked batch with
/// `LOSAT_X_BPMASKEDMIN=0`, with the per-hit loop's report and the batches' counts either way.
#[test]
fn the_rule_gives_tables_with_few_hits_per_cell_the_deferred_calls() {
    let (query, subject) = generated_pair("rule");
    // The first 60 letters of the first query (a random protein, its copy among the subjects).
    let text = std::fs::read_to_string(&query).expect("read the query set");
    let first = text.lines().nth(1).expect("the first protein");
    let short = query.with_file_name("short.faa");
    std::fs::write(&short, format!(">short\n{}\n", &first[..60])).expect("write FASTA");
    const FORCED: &[(&str, &str)] = &[("LOSAT_X_BPMASKEDMIN", "0")];
    for threads in [1, 4, 8, 16] {
        let what = format!("short query threads {threads}");
        let (reference, _, _, _log) =
            blastp(&short, &subject, threads, REFERENCE, &["-outfmt", "6"]);
        assert!(
            !reference.is_empty(),
            "{what}: the short query has a report"
        );
        let (deferred, deferred_counts, stderr, _log) =
            blastp(&short, &subject, threads, DEFERRED, &["-outfmt", "6"]);
        assert!(
            deferred == reference,
            "{what} {DEFERRED:?}: the report differs"
        );
        let d = one_line(&deferred_counts, &what);
        // With the switch off the rule is never consulted, except by LOSAT_X_BPCELL128's batches
        // (a CPU without AVX2), which LOSAT_X_BPMASKED does not turn off.
        let tables = masked_tables(&stderr);
        if !avx2() && cell128() {
            assert!(
                matches!(tables, Some((0, n)) if n > 0),
                "{what} {DEFERRED:?}: {tables:?}"
            );
        } else {
            assert_eq!(tables, None, "{what} {DEFERRED:?}: {tables:?}");
        }
        for (vars, rule_masks) in [(DEFAULT, false), (FORCED, true)] {
            let what = format!("{what} {vars:?}");
            let (rows, counts, stderr, _log) =
                blastp(&short, &subject, threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "{what}: the report differs");
            let c = one_line(&counts, &what);
            // LOSAT_X_BPCELL128's batches count only the tables they take; the per-hit loop of
            // the others counts nothing (the deferred calls count every table).
            if avx2() || !cell128() || !rule_masks {
                assert_eq!(
                    (c.cells, c.hits, c.candidates),
                    (d.cells, d.hits, d.candidates),
                    "{what}: counts"
                );
            }
            let tables = masked_tables(&stderr);
            if !avx2() && !cell128() {
                assert_eq!(tables, None, "{what}: no masked batch without AVX2");
            } else if rule_masks {
                assert!(
                    matches!(tables, Some((n, 0)) if n > 0),
                    "{what}: {tables:?}"
                );
            } else {
                assert!(
                    matches!(tables, Some((0, n)) if n > 0),
                    "{what}: {tables:?}"
                );
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_BPCELL128): the 128-bit batches (`=2`, in place of the AVX2 batches on this
// CPU) give the per-hit loop's reports with LOSAT_X_BPMASKEDMIN's rule and on every table, and on
// every table test the cells, hits and candidates of the scalar batch (`LOSAT_X_BPDIAG_SIMD=0`,
// which counts every table too). The generated pair's queries have enough hits per used cell for
// the batches on every table; the short query (about one hit per cell) keeps the per-hit loop
// under the rule.
#[test]
fn the_128_bit_batches_give_the_reports_of_the_reference() {
    if !cell128() {
        eprintln!("no SSSE3, SSE4.1 and POPCNT on this CPU: the 128-bit batches are not compared");
        return;
    }
    let (query, subject) = generated_pair("cell128");
    let text = std::fs::read_to_string(&query).expect("read the query set");
    let first = text.lines().nth(1).expect("the first protein");
    let short = query.with_file_name("short.faa");
    std::fs::write(&short, format!(">short\n{}\n", &first[..60])).expect("write FASTA");
    const SCALAR: &[(&str, &str)] = &[("LOSAT_X_BPDIAG_SIMD", "0")];
    const RULE: &[(&str, &str)] = &[("LOSAT_X_BPCELL128", "2")];
    const EVERY: &[(&str, &str)] = &[("LOSAT_X_BPCELL128", "2"), ("LOSAT_X_BPMASKEDMIN", "0")];
    for (q, threads, dense) in [(&query, 1, true), (&query, 4, true), (&short, 1, false)] {
        let what = format!("{} threads {threads}", q.display());
        let (reference, _, _, _log) = blastp(q, &subject, threads, REFERENCE, &["-outfmt", "6"]);
        assert!(!reference.is_empty(), "{what}: the pair has a report");
        let (rows, counts, _, _log) = blastp(q, &subject, threads, SCALAR, &["-outfmt", "6"]);
        assert!(rows == reference, "{what} {SCALAR:?}: rows differ");
        let scalar = one_line(&counts, &what);
        for (vars, every) in [(RULE, false), (EVERY, true)] {
            let what = format!("{what} {vars:?}");
            let (rows, counts, stderr, _log) =
                blastp(q, &subject, threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "{what}: rows differ");
            let c = one_line(&counts, &what);
            assert_eq!(c.mode, "w128", "{what}");
            let tables = masked_tables(&stderr);
            if every || dense {
                assert!(
                    matches!(tables, Some((n, 0)) if n > 0),
                    "{what}: {tables:?}"
                );
                assert_eq!(
                    (c.cells, c.hits, c.candidates),
                    (scalar.cells, scalar.hits, scalar.candidates),
                    "{what}: counts"
                );
                assert!(c.cells > 0 && c.candidates > 0, "{what}: {c:?}");
            } else {
                assert!(
                    matches!(tables, Some((0, n)) if n > 0),
                    "{what}: {tables:?}"
                );
                assert_eq!((c.cells, c.hits, c.candidates), (0, 0, 0), "{what}");
            }
        }
    }
}
