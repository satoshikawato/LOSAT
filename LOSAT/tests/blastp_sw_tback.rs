//! BLASTP `-use_sw_tback` (line RG item X4): what LOSAT decides of the status of NCBI's
//! Smith-Waterman redo (`blast_engine::redo_run::one_thread_redo`), on crafted inputs in
//! `tests/fasta/x4`. With `-searchsp 1000 -evalue 1e4` the e-value of a score of 0 passes the
//! threshold, so the redo of a match with two or more preliminary HSPs fails when the first query
//! letter scores the subject's first letter 0 or less (`smith_waterman_redo_window`): `tks_W`
//! against `sfail` (tks_M with D first, twice: two HSPs, W against D) in `sw_subjects.faa`. The
//! other queries (`qok1` and `qok2`, a protein of PajaWSV with W first) have one HSP each, with that
//! protein (`sok1`), and their redo succeeds; in `sw_endless_subjects.faa`, `send` (tks_M twice)
//! gives tks_M the redo that NCBI never ends (M against M). NCBI's one-thread redo keeps the status of the batch's last redone
//! match: a failed match before another is dropped, a failed last one fails the batch's traceback
//! (`LOSAT engine error: Traceback failed`, exit status 3) after the reports of the batches
//! before, unless the batch's first query holds a redone match: NCBI's clean-up then reads that
//! query's freed heap and crashes (SIGSEGV), and LOSAT refuses the batch after the same reports
//! (exit status 1; PD-LOSAT-NCBI-DEFECTS 1.53). `qnohit` (W 87 times) has no hit (the subjects
//! have no W). That the bytes equal NCBI's is the business of the fixture rows `db.x4.sw.fail.*`
//! (frozen on ws); these tests pin LOSAT's decisions and their independence of `-num_threads`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:251-253
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3796-3802 (the arrays of
//! every query's heap are freed, then the first query's heap is popped:
//! `BlastCompo_HeapPop`, compo_heap.c:444-457, reads `heapArray[1]` when its count is above 0)

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

struct Run {
    stdout: String,
    stderr: String,
    code: Option<i32>,
}

/// `losat blastp -query QUERY -subject tests/fasta/x4/sw_subjects.faa -use_sw_tback -searchsp 1000
/// -evalue 1e4 -outfmt OUTFMT -num_threads THREADS` with `env` as the only variables that choose a
/// code path or NCBI's batches.
fn sw_search(query: &str, outfmt: &str, threads: usize, env: &[(&str, &str)]) -> Run {
    sw_search_in("x4/sw_subjects.faa", query, outfmt, threads, env)
}

/// [`sw_search`] against the subjects `subjects`.
fn sw_search_in(
    subjects: &str,
    query: &str,
    outfmt: &str,
    threads: usize,
    env: &[(&str, &str)],
) -> Run {
    let options = ["-searchsp", "1000", "-evalue", "1e4"];
    sw_search_with(subjects, query, &options, outfmt, threads, env)
}

/// `losat blastp -query QUERY -subject SUBJECTS -use_sw_tback OPTIONS -outfmt OUTFMT -num_threads
/// THREADS` with `env` as the only variables that choose a code path or NCBI's batches.
fn sw_search_with(
    subjects: &str,
    query: &str,
    options: &[&str],
    outfmt: &str,
    threads: usize,
    env: &[(&str, &str)],
) -> Run {
    let output = losat_command()
        .envs(env.iter().copied())
        .arg("blastp")
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subjects))
        .arg("-use_sw_tback")
        .args(options)
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// `losat` without the variables that choose a code path or NCBI's batches.
fn losat_command() -> Command {
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
    command
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
const NCBI_CRASHES: &str = "with which NCBI BLAST+ 2.17.0's blastp crashes";

// Case a of review 1: the batch's only match fails (two HSPs, W against D), so its traceback fails.
#[test]
fn a_failed_only_match_fails_the_batch() {
    for threads in [1, 4] {
        let run = sw_search("rps/tks_w.faa", "6", threads, &[]);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, "", "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// Case b: the failed match of the batch's first query is dropped, and the second query's redone
// match resets the status: the batch completes with the second query's rows.
#[test]
fn b_a_failed_match_before_a_redone_one_is_dropped() {
    let alone = sw_search("x4/sw_ok.faa", "6", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    assert!(row_queries(&alone.stdout).iter().all(|&q| q == "qok1"));
    for threads in [1, 4] {
        let run = sw_search("x4/sw_fail_first.faa", "6", threads, &[]);
        assert_eq!(run.code, Some(0), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert_eq!(run.stderr, "", "n{threads}");
    }
}

// Case c, the mirror of b: the failing query is the batch's last, and the batch's first query
// (qok1) holds a redone match, on which NCBI's clean-up crashes: LOSAT refuses the batch before any
// report (the outfmt 0 prolog is written).
#[test]
fn c_a_failed_last_match_after_a_kept_first_query_is_refused() {
    for threads in [1, 4] {
        for outfmt in ["6", "7", "0"] {
            let run = sw_search("x4/sw_fail_last.faa", outfmt, threads, &[]);
            assert_eq!(
                run.code,
                Some(1),
                "n{threads} outfmt {outfmt}: {}",
                run.stderr
            );
            assert!(
                run.stderr.contains(NCBI_CRASHES),
                "n{threads} outfmt {outfmt}: {}",
                run.stderr
            );
            if outfmt == "0" {
                assert!(run.stdout.starts_with("BLASTP "), "n{threads}");
                assert!(!run.stdout.contains("Query="), "n{threads}");
            } else {
                assert_eq!(run.stdout, "", "n{threads} outfmt {outfmt}");
            }
        }
    }
}

// Case d: one query per batch (`BATCH_SIZE=1`), the failure in batch 2 of 3 (whose first query,
// tks_W, holds no redone match): batch 1's report is written, without the epilog, then the error;
// batch 3 is not searched. In one batch (the default
// batch size) the third query's match is the last and the run completes without the failed match.
#[test]
fn d_a_failure_in_batch_2_of_3_comes_after_batch_1() {
    let one_per_batch = [("BATCH_SIZE", "1")];
    let alone = sw_search("x4/sw_ok.faa", "6", 1, &[]);
    let alone7 = sw_search("x4/sw_ok.faa", "7", 1, &[]);
    let alone0 = sw_search("x4/sw_ok.faa", "0", 1, &[]);
    for threads in [1, 4] {
        let run = sw_search("x4/sw_fail_mid.faa", "6", threads, &one_per_batch);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
        // Outfmt 7: batch 1's lines without the last one (`# BLAST processed 1 queries`, the
        // epilog).
        let commented = sw_search("x4/sw_fail_mid.faa", "7", threads, &one_per_batch);
        assert_eq!(commented.code, Some(3), "n{threads}");
        let epilog = alone7
            .stdout
            .rfind("# BLAST processed")
            .expect("outfmt 7's last line");
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        // Outfmt 0: the prolog and batch 1's report, a part of the query's report alone.
        let pairwise = sw_search("x4/sw_fail_mid.faa", "0", threads, &one_per_batch);
        assert_eq!(pairwise.code, Some(3), "n{threads}");
        assert!(pairwise.stdout.contains("Query= qok1"), "n{threads}");
        assert!(!pairwise.stdout.contains("Query= tks_W"), "n{threads}");
        assert!(alone0.stdout.starts_with(&pairwise.stdout), "n{threads}");
        assert_eq!(pairwise.stderr, TRACEBACK_FAILED, "n{threads}");
        let one_batch = sw_search("x4/sw_fail_mid.faa", "6", threads, &[]);
        assert_eq!(one_batch.code, Some(0), "n{threads}: {}", one_batch.stderr);
        let queries = row_queries(&one_batch.stdout);
        assert!(
            queries.contains(&"qok1") && queries.contains(&"qok2"),
            "n{threads}"
        );
        assert!(!queries.contains(&"tks_W"), "n{threads}");
    }
}

// The redo that NCBI never ends (tks_M against send: two HSPs, M against M) in batch 2 of 2
// (`BATCH_SIZE=1`): LOSAT refuses it after batch 1's report (without the epilog), as NCBI's search
// of batch 2 never ends after NCBI has printed batch 1 (PD-LOSAT-NCBI-DEFECTS 1.53); in one batch
// nothing comes before the refusal.
#[test]
fn e_the_endless_redo_is_refused_after_the_batches_before() {
    let subjects = "x4/sw_endless_subjects.faa";
    let alone = sw_search_in(subjects, "x4/sw_ok.faa", "6", 1, &[]);
    let alone7 = sw_search_in(subjects, "x4/sw_ok.faa", "7", 1, &[]);
    let epilog = alone7
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line");
    for threads in [1, 4] {
        let run = sw_search_in(
            subjects,
            "x4/sw_endless_last.faa",
            "6",
            threads,
            &[("BATCH_SIZE", "1")],
        );
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert!(
            run.stderr.contains("does not end in NCBI BLAST+ 2.17.0"),
            "n{threads}: {}",
            run.stderr
        );
        let commented = sw_search_in(
            subjects,
            "x4/sw_endless_last.faa",
            "7",
            threads,
            &[("BATCH_SIZE", "1")],
        );
        assert_eq!(commented.code, Some(1), "n{threads}");
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        let one_batch = sw_search_in(subjects, "x4/sw_endless_last.faa", "6", threads, &[]);
        assert_eq!(one_batch.code, Some(1), "n{threads}");
        assert_eq!(one_batch.stdout, "", "n{threads}");
    }
}

// The batch's first query holds no redone match, so NCBI's clean-up pops an empty heap and the
// batch ends with NCBI's error, as in case a: f, the first query has no hit (qnohit); g, its only
// match fails (tks_W, then qok1's redone match, then tks_W2's failed match last). A redone match of
// a later query (qok1) does not matter: NCBI pops the first query's heap alone.
#[test]
fn f_g_a_first_query_without_a_redone_match_ends_in_traceback_failed() {
    for query in [
        "x4/sw_fail_last_nohit_first.faa",
        "x4/sw_fail_last_failed_first.faa",
    ] {
        for threads in [1, 4] {
            let run = sw_search(query, "6", threads, &[]);
            assert_eq!(run.code, Some(3), "{query} n{threads}: {}", run.stderr);
            assert_eq!(run.stdout, "", "{query} n{threads}");
            assert_eq!(run.stderr, TRACEBACK_FAILED, "{query} n{threads}");
        }
    }
}

// Case h: batch 1 is [qok2, qnohit] and batch 2 [qok1, tks_W] (`BATCH_SIZE=200`: 124 + 87, then
// 124 + 80 letters): batch 1's report is written, without the epilog, then the refusal of batch 2,
// on whose clean-up NCBI crashes.
#[test]
fn h_the_refusal_comes_after_the_batches_before() {
    let batches = [("BATCH_SIZE", "200")];
    let alone = sw_search("x4/sw_ok.faa", "6", 1, &[]);
    for threads in [1, 4] {
        let run = sw_search("x4/sw_crash_batch2.faa", "6", threads, &batches);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(
            run.stderr.contains(NCBI_CRASHES),
            "n{threads}: {}",
            run.stderr
        );
        // qok2 is qok1's sequence: the same rows under its id.
        assert!(!alone.stdout.is_empty());
        assert_eq!(
            run.stdout,
            alone.stdout.replace("qok1\t", "qok2\t"),
            "n{threads}"
        );
        let pairwise = sw_search("x4/sw_crash_batch2.faa", "0", threads, &batches);
        assert_eq!(pairwise.code, Some(1), "n{threads}: {}", pairwise.stderr);
        assert!(pairwise.stdout.contains("Query= qok2"), "n{threads}");
        assert!(pairwise.stdout.contains("Query= qnohit"), "n{threads}");
        assert!(!pairwise.stdout.contains("Query= qok1"), "n{threads}");
        assert!(pairwise.stderr.contains(NCBI_CRASHES), "n{threads}");
    }
}

// Case x5 (line RG item X5, review R3 M1): unified P-values without `-searchsp`, the batches
// [qok1], [tks_W] and [qshort, allx] (`BATCH_SIZE=80`: 124, 80, then 40 + 50 letters). Batch 3 is
// refused before any redo (a query with letters to score followed by one without a search space),
// but NCBI stops before it: batch 2's traceback fails (case d). LOSAT's search of all batches
// meets the refusal of batch 3 first and searches batches 1 and 2 again, where the redo fails;
// it then writes batch 1's report alone, then batch 2's error, as NCBI does (the build before the
// fix wrote nothing on stdout here).
#[test]
fn x5_a_failure_in_the_redo_comes_before_a_later_batchs_refusal() {
    let unified = ["-evalue", "1e4", "-comp_based_stats", "2u"];
    let batches = [("BATCH_SIZE", "80")];
    let run = |query: &str, outfmt: &str, threads: usize, env: &[(&str, &str)]| {
        sw_search_with("x4/sw_subjects.faa", query, &unified, outfmt, threads, env)
    };
    let alone = run("x4/sw_ok.faa", "6", 1, &[]);
    let alone7 = run("x4/sw_ok.faa", "7", 1, &[]);
    let alone0 = run("x4/sw_ok.faa", "0", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    for threads in [1, 4] {
        let tabular = run("x5/sw_order.faa", "6", threads, &batches);
        assert_eq!(tabular.code, Some(3), "n{threads}: {}", tabular.stderr);
        assert_eq!(tabular.stdout, alone.stdout, "n{threads}");
        assert_eq!(tabular.stderr, TRACEBACK_FAILED, "n{threads}");
        let commented = run("x5/sw_order.faa", "7", threads, &batches);
        assert_eq!(commented.code, Some(3), "n{threads}");
        let epilog = alone7
            .stdout
            .rfind("# BLAST processed")
            .expect("outfmt 7's last line");
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        let pairwise = run("x5/sw_order.faa", "0", threads, &batches);
        assert_eq!(pairwise.code, Some(3), "n{threads}");
        assert!(pairwise.stdout.contains("Query= qok1"), "n{threads}");
        assert!(!pairwise.stdout.contains("Query= tks_W"), "n{threads}");
        assert!(alone0.stdout.starts_with(&pairwise.stdout), "n{threads}");
        assert_eq!(pairwise.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// The rows `db.x4excl.*` (`tests/db_search_manifest/blastp.x4excl.tsv`, line RG item X4-excl): a
// `-negative_seqidlist` that keeps a subject of a version 5 database (`tests/db_search/x4excl`)
// without any of its definition lines, so NCBI's fetch excludes that subject's match at the redo,
// which sets the batch's status back to 0 unless the redo ends early at it
// (`redo_run::query_fetch_excluded_reached`). The rows' bytes are NCBI's at one thread, which
// `db_search_fixtures.py check` compares on ws; this test is what CI runs of them. Each row runs
// at `-num_threads` 1, 2, 4 and 8 with every redo path (`LOSAT_X_BPREDOPAR` unset and `=0`) and
// must exit with the row's status and, at every count and path, write the bytes of the first run;
// a tabular row (outfmt 6) of kind `ncbi` must also write the frozen NCBI stdout. Rows
// `last.kept.*` and `last.q2.*` keep two matches or more, so above one thread they reach a redo
// with more than one job: a redo path there that does not read the places of the excluded lists
// (`RedoInputs::fetch_excluded`) fails or refuses them, where NCBI exits 0.
#[test]
fn x4excl_rows_end_alike_at_every_thread_count_and_redo_path() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let manifest = std::fs::read_to_string(root.join("tests/db_search_manifest/blastp.x4excl.tsv"))
        .expect("read blastp.x4excl.tsv");
    let mut rows = 0;
    for line in manifest.lines() {
        if line.starts_with('#') || line.starts_with("row_id\t") || line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let (row_id, kind, db, query, extra, outfmt, rc) = (
            cells[0], cells[2], cells[3], cells[4], cells[5], cells[6], cells[7],
        );
        assert!(row_id.starts_with("db.x4excl."), "{row_id}");
        rows += 1;
        let run = |threads: usize, env: &[(&str, &str)]| {
            let output = losat_command()
                // The paths of the row are relative to LOSAT/, where the fixture tool runs both
                // programs; no `.ncbirc` is in that directory.
                .current_dir(root)
                .env("HOME", root)
                .env("USERPROFILE", root)
                .env_remove("BLASTDB")
                .envs(env.iter().copied())
                .args(["blastp", "-db", db, "-query", query])
                .args(extra.split_whitespace())
                .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
                .output()
                .expect("run LOSAT");
            (output.status.code(), output.stdout, output.stderr)
        };
        let first = run(1, &[]);
        let expected_code: i32 = rc.parse().expect("rc");
        assert_eq!(
            first.0,
            Some(expected_code),
            "{row_id}: {}",
            String::from_utf8_lossy(&first.2)
        );
        if kind == "ncbi" && outfmt == "6" {
            let frozen = std::fs::read(root.join(format!("tests/fixtures/db_search/{row_id}.out")))
                .expect("read the frozen stdout");
            assert_eq!(first.1, frozen, "{row_id}: the frozen NCBI stdout");
        }
        for threads in [1, 2, 4, 8] {
            for env in [&[][..], &[("LOSAT_X_BPREDOPAR", "0")][..]] {
                let again = run(threads, env);
                assert_eq!(
                    again,
                    first,
                    "{row_id} n{threads} {env:?}: {}",
                    String::from_utf8_lossy(&again.2)
                );
            }
        }
    }
    assert_eq!(rows, 13, "the rows of blastp.x4excl.tsv");
}
