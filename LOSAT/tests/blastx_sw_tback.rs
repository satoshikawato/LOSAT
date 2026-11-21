//! BLASTX `-use_sw_tback` (line RG item X4b-2): what LOSAT decides of NCBI's Smith-Waterman redo of
//! a match (`algorithm::blastx::kappa`, `runtime::kappa_batch`), on X4's crafted proteins
//! back-translated into nucleotide queries (`tests/fasta/x4b2`, `tests/db_search/x4b2/make_x4b2.py`).
//! With `-searchsp 1000 -evalue 1e4 -sum_stats false` the e-value of a score of 0 passes the
//! threshold, so a window of two or more preliminary HSPs makes NCBI's round loop reach a score of
//! 0: tks_W's redo fails against sfail (W against D), tks_M's never ends against X4's endless
//! subjects (M against M); `pre20` (the first 20 letters of X4's qok1) keeps one match against
//! sok1. NCBI's one-thread redo keeps the status of the batch's last redone match: a failed last
//! one fails the batch's traceback (`LOSAT engine error: Traceback failed`, exit status 3) after
//! the reports of the batches before, unless the batch's first query holds a redone match, where
//! NCBI's clean-up crashes and LOSAT refuses the batch (exit status 1). In e-value mode NCBI tests
//! the heap that a window's query context (six per query) indexes in its array of the batch's
//! heaps (one per query): a context below the batch's query count reads another query's heap,
//! which LOSAT reproduces (`q12.fna`: two queries with hits, ten without); at or past it NCBI
//! reads past the array, where LOSAT tests an empty heap (approved exception 20 of
//! PD-LOSAT-NCBI-DEFECTS): NCBI's result for the batch followed by five queries without a hit
//! per query (`*_pad.fna`).
//! `-ungapped -comp_based_stats 0 -use_sw_tback` crashes NCBI at the redo's set-up of every batch
//! it searches. That the bytes equal NCBI's is the business of the fixture rows `db.x4b2.*`; these
//! tests pin LOSAT's decisions and their independence of `-num_threads`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1465-1474

use std::path::{Path, PathBuf};
use std::process::Command;

fn input(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(name)
}

struct Run {
    stdout: String,
    stderr: String,
    code: Option<i32>,
}

/// `losat blastx ARGS -num_threads THREADS` with `env` as the only variables that choose a code
/// path or NCBI's batches.
fn blastx(args: &[&str], threads: usize, env: &[(&str, &str)]) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        // `LOSAT_BLASTX_WORKER_LOG` is the log file of the build with `--all-features` (its worker
        // probe), not a code path.
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .envs(env.iter().copied())
        .arg("blastx")
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// `blastx -query QUERY -subject SUBJECT -searchsp 1000 -evalue 1e4 -sum_stats false
/// -use_sw_tback -outfmt OUTFMT`: the score of 0 significant, e-value mode.
fn sw_search(
    query: &str,
    subject: &str,
    outfmt: &str,
    threads: usize,
    env: &[(&str, &str)],
) -> Run {
    let (query, subject) = (input(query), input(subject));
    blastx(
        &[
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-searchsp",
            "1000",
            "-evalue",
            "1e4",
            "-sum_stats",
            "false",
            "-use_sw_tback",
            "-outfmt",
            outfmt,
        ],
        threads,
        env,
    )
}

/// The query ids of the rows of a tabular report.
fn row_queries(stdout: &str) -> Vec<&str> {
    stdout
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').next())
        .collect()
}

const TRACEBACK_FAILED: &str = "LOSAT engine error: Traceback failed\n";
/// The part of LOSAT's refusal of a failed batch whose first query holds a redone match.
const NCBI_CRASHES: &str = "with which NCBI BLAST+ 2.17.0's blastx crashes";
/// The part of LOSAT's refusal of a redo that NCBI never ends.
const NEVER_ENDS: &str = "does not end in NCBI BLAST+ 2.17.0";
const OK_FAIL: &str = "fasta/x4b2/ok_fail.faa";

// tks_W's only redone match (on sfail, the stream's last) fails: the batch's traceback fails
// (NCBI 2.17.0: probes f.w.okfail.* of item X4b-2, rows db.x4b2.fail.w.*).
#[test]
fn a_failed_only_match_fails_the_batch() {
    for threads in [1, 4] {
        for outfmt in ["6", "7"] {
            let run = sw_search("fasta/x4b2/w.fna", OK_FAIL, outfmt, threads, &[]);
            assert_eq!(run.code, Some(3), "n{threads} {outfmt}: {}", run.stderr);
            assert_eq!(run.stdout, "", "n{threads} {outfmt}");
            assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads} {outfmt}");
        }
        // The batch's first query has no hit, tks_W is its second: the same failure.
        let run = sw_search("fasta/x4b2/nohit_w12.fna", OK_FAIL, "6", threads, &[]);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// pre20's match is kept in the first query's heap, tks_W's (context 6, below the batch's twelve
// queries) fails last: NCBI's clean-up crashes and LOSAT refuses the batch. Split by BATCH_SIZE
// into [pre20 and eleven queries without a hit] and [tks_W], the first batch is reported and the
// second fails (exit 3).
#[test]
fn a_failed_last_match_after_a_kept_first_query_is_refused() {
    for threads in [1, 4] {
        for outfmt in ["6", "0"] {
            let run = sw_search("fasta/x4b2/pre20_w12.fna", OK_FAIL, outfmt, threads, &[]);
            assert_eq!(run.code, Some(1), "n{threads} {outfmt}: {}", run.stderr);
            assert!(run.stderr.contains(NCBI_CRASHES), "{}", run.stderr);
            if outfmt == "0" {
                assert!(
                    run.stdout.starts_with(concat!(
                        "BLASTX (LOSAT ",
                        env!("CARGO_PKG_VERSION"),
                        ")"
                    )),
                    "n{threads}"
                );
                assert!(!run.stdout.contains("Query="), "n{threads}");
            } else {
                assert_eq!(run.stdout, "", "n{threads}");
            }
        }
        let batches = [("BATCH_SIZE", "1644")];
        let run = sw_search(
            "fasta/x4b2/pre20_pad_w.fna",
            OK_FAIL,
            "6",
            threads,
            &batches,
        );
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(row_queries(&run.stdout), ["pre20_nt"], "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// A window whose query context is at or past the batch's query count, where NCBI's heap test reads
// past its array of the batch's heaps: LOSAT tests an empty heap there (approved exception 20 of
// PD-LOSAT-NCBI-DEFECTS), so a batch gives the rows of the same batch followed by five queries
// without a hit per query (`*_pad.fna`, whose NCBI runs are the rows `db.x4b2.past.*.pad.*`), at
// every thread count.
#[test]
fn a_heap_test_past_the_batch_reads_an_empty_heap() {
    let subject = input(OK_FAIL);
    let bslice = input("db_search/bslice.faa");
    // (query, padded query, subject, as-is batch size, padded batch size): tks_W alone (contexts 1
    // to 5 past one heap), q12 one query per batch (each with five pads in a batch of 1720 letters).
    let cases = [
        (
            "fasta/x4b2/w.fna",
            "fasta/x4b2/w_pad.fna",
            &subject,
            None,
            None,
        ),
        (
            "fasta/x4b2/q12.fna",
            "fasta/x4b2/q12_b1_pad.fna",
            &bslice,
            Some("1"),
            Some("1720"),
        ),
    ];
    for (query, padded, subject, size, padded_size) in cases {
        let run = |query: &str, size: Option<&str>, threads: usize| {
            let query = input(query);
            let env: Vec<(&str, &str)> =
                size.map(|size| ("BATCH_SIZE", size)).into_iter().collect();
            blastx(
                &[
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-sum_stats",
                    "false",
                    "-use_sw_tback",
                    "-outfmt",
                    "6 qseqid sseqid qstart qend sstart send evalue score qframe",
                ],
                threads,
                &env,
            )
        };
        let expected = run(padded, padded_size, 1);
        assert_eq!(expected.code, Some(0), "{padded}: {}", expected.stderr);
        assert!(
            expected.stdout.lines().any(|line| !line.ends_with("\t1")),
            "{padded}: a row in a frame other than +1"
        );
        for threads in [1, 4, 8] {
            let as_is = run(query, size, threads);
            assert_eq!(as_is.code, Some(0), "n{threads} {query}: {}", as_is.stderr);
            assert_eq!(as_is.stdout, expected.stdout, "n{threads} {query}");
        }
    }
    for threads in [1, 4] {
        // tks_W in frame -1 alone (context 3): the empty heap takes the window, whose redo fails,
        // as in NCBI's run with the five queries after it (exit 3).
        let run = sw_search("fasta/x4b2/w_rc.fna", OK_FAIL, "6", threads, &[]);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
        // tks_W as the second of two queries (context 6) and as the thirteenth of thirteen
        // (context 72): the window is taken and fails last after pre20's kept match, where NCBI's
        // clean-up crashes (with the queries without a hit after it too): LOSAT refuses the batch.
        for query in ["fasta/x4b2/pre20_w.fna", "fasta/x4b2/pre20_pad_w.fna"] {
            let run = sw_search(query, OK_FAIL, "6", threads, &[]);
            assert_eq!(run.code, Some(1), "n{threads} {query}: {}", run.stderr);
            assert!(run.stderr.contains(NCBI_CRASHES), "{query}: {}", run.stderr);
            assert_eq!(run.stdout, "", "n{threads} {query}");
        }
    }
}

// tks_M's redo never ends in NCBI: LOSAT refuses it.
#[test]
fn the_endless_redo_is_refused() {
    for threads in [1, 4] {
        let run = sw_search(
            "fasta/x4b2/m.fna",
            "fasta/x4/sw_endless_subjects.faa",
            "6",
            threads,
            &[],
        );
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, "", "n{threads}");
        assert!(run.stderr.contains(NEVER_ENDS), "{}", run.stderr);
        assert!(run.stderr.contains("LOSAT's BLASTX"), "{}", run.stderr);
    }
}

// `-ungapped -comp_based_stats 0 -use_sw_tback`: NCBI crashes at the redo's set-up of every batch
// it searches, a batch without hits too; LOSAT refuses at the first one (probes u.*).
#[test]
fn the_ungapped_search_is_refused_at_its_first_batch() {
    let subject = input(OK_FAIL);
    for query in ["fasta/x4b2/w.fna", "db_search/bslice_nohit.fna"] {
        let query = input(query);
        for threads in [1, 4] {
            let run = blastx(
                &[
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-ungapped",
                    "-comp_based_stats",
                    "0",
                    "-use_sw_tback",
                    "-outfmt",
                    "6",
                ],
                threads,
                &[],
            );
            assert_eq!(run.code, Some(1), "n{threads} {query:?}: {}", run.stderr);
            assert_eq!(run.stdout, "", "n{threads} {query:?}");
            assert!(
                run.stderr.contains(
                    "-use_sw_tback with -ungapped, on which NCBI BLAST+ 2.17.0's blastx crashes"
                ),
                "{}",
                run.stderr
            );
        }
    }
}

// The e-value mode on twelve queries, whose windows' contexts (below twelve) read the batch's
// heaps, and the default mode (sum statistics): the result of one thread at four, with every
// composition mode.
#[test]
fn the_redo_gives_the_one_thread_result_at_four_threads() {
    let (query, subject) = (input("fasta/x4b2/q12.fna"), input("db_search/bslice.faa"));
    for mode in [&["-sum_stats", "false"][..], &[][..]] {
        for cbs in ["0", "1", "2", "3"] {
            let mut args = vec![
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-comp_based_stats",
                cbs,
                "-use_sw_tback",
                "-outfmt",
                "6",
            ];
            args.extend_from_slice(mode);
            let one = blastx(&args, 1, &[]);
            assert_eq!(one.code, Some(0), "{mode:?} CBS {cbs}: {}", one.stderr);
            assert!(!one.stdout.is_empty(), "{mode:?} CBS {cbs}");
            let four = blastx(&args, 4, &[]);
            assert_eq!(four.code, Some(0), "{mode:?} CBS {cbs}: {}", four.stderr);
            assert_eq!(four.stdout, one.stdout, "{mode:?} CBS {cbs}");
        }
    }
}
