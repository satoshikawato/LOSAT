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
//! (`BLAST engine error: Traceback failed`, exit status 3) after the reports of the batches
//! before. That the bytes equal NCBI's is the business of the fixture rows `db.x4.sw.fail.*`
//! (frozen on ws); these tests pin LOSAT's decisions and their independence of `-num_threads`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3774-3778
//! ```c
//!     for (i = 0; i < actual_num_threads; ++i) {
//!         if (status_code_tld[i] != 0) {
//!             status_code = status_code_tld[i];
//!         }
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:250-252
//! ```c++
//!     if (status) {
//!         NCBI_THROW(CBlastException, eCoreBlastError, "Traceback failed");
//!     }
//! ```

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
        .envs(env.iter().copied())
        .arg("blastp")
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subjects))
        .args(["-use_sw_tback", "-searchsp", "1000", "-evalue", "1e4"])
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// The query ids of the rows of a tabular report.
fn row_queries(stdout: &str) -> Vec<&str> {
    stdout
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').next())
        .collect()
}

const TRACEBACK_FAILED: &str = "BLAST engine error: Traceback failed\n";

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

// Case c, the mirror of b: the failing query is the batch's last, so the batch fails before any
// report (the outfmt 0 prolog is written).
#[test]
fn c_a_failed_last_match_fails_the_batch() {
    for threads in [1, 4] {
        let run = sw_search("x4/sw_fail_last.faa", "6", threads, &[]);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, "", "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
        let pairwise = sw_search("x4/sw_fail_last.faa", "0", threads, &[]);
        assert_eq!(pairwise.code, Some(3), "n{threads}: {}", pairwise.stderr);
        assert!(pairwise.stdout.starts_with("BLASTP "), "n{threads}");
        assert!(!pairwise.stdout.contains("Query="), "n{threads}");
        assert_eq!(pairwise.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// Case d: one query per batch (`BATCH_SIZE=1`), the failure in batch 2 of 3: batch 1's report is
// written, without the epilog, then the error; batch 3 is not searched. In one batch (the default
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
// of batch 2 never ends after NCBI has printed batch 1 (PD-LOSAT-NCBI-DEFECTS 1.52); in one batch
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
