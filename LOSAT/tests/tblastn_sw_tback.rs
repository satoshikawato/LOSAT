//! TBLASTN `-use_sw_tback` (line RG item X4b-1): what LOSAT decides of the status of NCBI's
//! Smith-Waterman redo of a match (`algorithm::tblastn::match_redo`), on X4's crafted proteins
//! back-translated into nucleotide subjects (`tests/fasta/x4b`). With `-searchsp 1000 -evalue 1e4
//! -sum_stats false` the e-value of a score of 0 passes the threshold, so a window of two or more
//! preliminary HSPs makes NCBI's round loop reach a score of 0: the redo of the match fails when
//! the window's first cell scores 0 or less (`tks_W` against `sfail_rc`, X4's sfail in frame -1: W
//! against D), and never ends when it scores above 0 (`tks_M` against `send_nt`: M against M).
//! `qpre20` (the first 20 letters of X4's qok1) has one window of one HSP against `sok1_nt` (X4's
//! sok1 in frame +1), whose redo succeeds; `qpre40` (its first 40 letters) also fails against
//! `sfail_rc`. NCBI's one-thread redo keeps the status of the batch's last redone match: a failed
//! last one fails the batch's traceback (NCBI's `BLAST engine error: Traceback failed`, LOSAT's
//! `LOSAT engine error: Traceback failed`, exit status 3) after the reports of the batches before,
//! unless the batch's first query holds a redone match:
//! NCBI's clean-up then reads that query's freed heap and crashes, and LOSAT refuses the batch
//! after the same reports (exit status 1; PD-LOSAT-NCBI-DEFECTS 1.73). That the bytes equal NCBI's
//! is the business of the fixture rows `db.x4b.*`; these tests pin LOSAT's decisions and their
//! independence of `-num_threads`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:251-253

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

/// `losat tblastn ARGS -num_threads THREADS` with `env` as the only variables that choose a code
/// path or NCBI's batches.
fn tblastn(args: &[&str], threads: usize, env: &[(&str, &str)]) -> Run {
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
        .arg("tblastn")
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

/// `tblastn -query QUERY -subject SUBJECT -searchsp 1000 -evalue 1e4 -sum_stats false
/// -use_sw_tback -outfmt OUTFMT`: the score of 0 significant.
fn sw_search(
    query: &str,
    subject: &str,
    outfmt: &str,
    threads: usize,
    env: &[(&str, &str)],
) -> Run {
    let (query, subject) = (input(query), input(subject));
    tblastn(
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
const NCBI_CRASHES: &str = "with which NCBI BLAST+ 2.17.0's tblastn crashes";
/// The part of LOSAT's refusal of a redo that NCBI never ends.
const NEVER_ENDS: &str = "does not end in NCBI BLAST+ 2.17.0";
const FAIL_RC: &str = "fasta/x4b/ok1_fail_rc.fna";

// The batch's only redone match fails: its traceback fails (NCBI 2.17.0: probe tbn.fail.a.rc of
// line RG X4b's survey, rows db.x4b.fail.a.*).
#[test]
fn a_failed_only_match_fails_the_batch() {
    for threads in [1, 4] {
        let run = sw_search(
            "fasta/rps/tks_w.faa",
            "fasta/x4b/t_sw_subjects_rc.fna",
            "6",
            threads,
            &[],
        );
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, "", "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
    }
}

// qpre20's match succeeds, tks_W's fails last while qpre20's heap holds a match: NCBI's clean-up
// crashes and LOSAT refuses the batch (before any report; the outfmt 0 prolog is written).
// One query per batch: batch 1 (qpre20) is reported, batch 2 (tks_W, whose first query holds no
// match) ends in NCBI's error.
#[test]
fn a_failed_last_match_after_a_kept_first_query_is_refused() {
    let alone = sw_search("fasta/x4b/q_pre20.faa", FAIL_RC, "6", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert_eq!(row_queries(&alone.stdout), ["qpre20"]);
    let alone0 = sw_search("fasta/x4b/q_pre20.faa", FAIL_RC, "0", 1, &[]);
    for threads in [1, 4] {
        for outfmt in ["6", "0"] {
            let run = sw_search("fasta/x4b/q_pre20_w.faa", FAIL_RC, outfmt, threads, &[]);
            assert_eq!(run.code, Some(1), "n{threads} {outfmt}: {}", run.stderr);
            assert!(run.stderr.contains(NCBI_CRASHES), "{}", run.stderr);
            if outfmt == "0" {
                assert!(run.stdout.starts_with("TBLASTN "), "n{threads}");
                assert!(!run.stdout.contains("Query="), "n{threads}");
            } else {
                assert_eq!(run.stdout, "", "n{threads}");
            }
        }
        let batches = [("BATCH_SIZE", "1")];
        let run = sw_search("fasta/x4b/q_pre20_w.faa", FAIL_RC, "6", threads, &batches);
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads}");
        let pairwise = sw_search("fasta/x4b/q_pre20_w.faa", FAIL_RC, "0", threads, &batches);
        assert_eq!(pairwise.code, Some(3), "n{threads}");
        assert!(pairwise.stdout.contains("Query= qpre20"), "n{threads}");
        assert!(!pairwise.stdout.contains("Query= tks_W"), "n{threads}");
        assert!(alone0.stdout.starts_with(&pairwise.stdout), "n{threads}");
        // Batch 2 (qpre40) fails last after its own kept match: refused after batch 1's report.
        let crash = sw_search(
            "fasta/x4b/q_pre20_pre40.faa",
            FAIL_RC,
            "6",
            threads,
            &batches,
        );
        assert_eq!(crash.code, Some(1), "n{threads}: {}", crash.stderr);
        assert_eq!(crash.stdout, alone.stdout, "n{threads}");
        assert!(crash.stderr.contains(NCBI_CRASHES), "{}", crash.stderr);
    }
}

// tks_M's redo never ends in NCBI: LOSAT refuses it, after batch 1's report with one query per
// batch.
#[test]
fn the_endless_redo_is_refused_after_the_batches_before() {
    let subjects = "fasta/x4b/t_sw_endless_subjects.fna";
    let alone = sw_search("fasta/x4b/q_pre20.faa", subjects, "6", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert_eq!(row_queries(&alone.stdout), ["qpre20"]);
    for threads in [1, 4] {
        let run = sw_search("fasta/rps/tks_m.faa", subjects, "6", threads, &[]);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, "", "n{threads}");
        assert!(run.stderr.contains(NEVER_ENDS), "{}", run.stderr);
        assert!(run.stderr.contains("LOSAT's TBLASTN"), "{}", run.stderr);
        let batches = [("BATCH_SIZE", "1")];
        let run = sw_search("fasta/x4b/q_pre20_m.faa", subjects, "6", threads, &batches);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert!(run.stderr.contains(NEVER_ENDS), "{}", run.stderr);
    }
}

// The default mode with `-db`: sum statistics with SeqDB's shortest length 10, whose cutoff score
// is `BLAST_Cutoffs`' floor of 1 (32 scaled with composition-based statistics). NCBI 2.17.0 reports
// 3 alignments in both modes (rows db.x4b.cut.*); a cutoff below 0 made a score of 0 significant.
#[test]
fn the_default_mode_with_a_database_reads_the_cutoff_of_1() {
    let (query, db) = (input("fasta/x4b/tq6_160.faa"), input("db_search/tslice"));
    for cbs in ["0", "2"] {
        let args = [
            "-query",
            query.to_str().unwrap(),
            "-db",
            db.to_str().unwrap(),
            "-comp_based_stats",
            cbs,
            "-use_sw_tback",
            "-outfmt",
            "6",
        ];
        let one = tblastn(&args, 1, &[]);
        assert_eq!(one.code, Some(0), "CBS {cbs}: {}", one.stderr);
        assert_eq!(row_queries(&one.stdout).len(), 3, "CBS {cbs}");
        let four = tblastn(&args, 4, &[]);
        assert_eq!(four.stdout, one.stdout, "CBS {cbs}");
    }
}

// A PSSM query's Smith-Waterman redo is refused at its first match, after the prolog (NCBI 2.17.0
// fails, crashes, does not end or writes rows that differ from run to run, by the input;
// `match_redo::sw_pssm_refusal`); a search without hits ends as NCBI's (exit 0, rows
// db.x4b.pssm.*).
#[test]
fn a_pssm_query_is_refused_at_its_first_match() {
    let pssm = input("fixtures/db_search/db.ps.asn.rnd.6.asn");
    for (subject, outfmt) in [("db_search/tslice.fna", "6"), ("db_search/tslice.fna", "0")] {
        let subject = input(subject);
        let run = tblastn(
            &[
                "-in_pssm",
                pssm.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-use_sw_tback",
                "-outfmt",
                outfmt,
            ],
            1,
            &[],
        );
        assert_eq!(run.code, Some(1), "{}", run.stderr);
        assert!(
            run.stderr
                .contains("-use_sw_tback with -in_pssm on a search with hits"),
            "{}",
            run.stderr
        );
        if outfmt == "6" {
            assert_eq!(run.stdout, "");
        } else {
            assert!(
                run.stdout.starts_with(concat!(
                    "PSITBLASTN (LOSAT ",
                    env!("CARGO_PKG_VERSION"),
                    ")"
                )),
                "{}",
                run.stdout
            );
            assert!(!run.stdout.contains("Query="), "{}", run.stdout);
        }
    }
    let none = input("fasta/x4b/polycg.fna");
    for threads in [1, 4] {
        let run = tblastn(
            &[
                "-in_pssm",
                pssm.to_str().unwrap(),
                "-subject",
                none.to_str().unwrap(),
                "-use_sw_tback",
                "-outfmt",
                "6",
            ],
            threads,
            &[],
        );
        assert_eq!(run.code, Some(0), "{}", run.stderr);
        assert_eq!(run.stdout, "");
    }
}

// IDENTITY, which NCBI searches without composition-based statistics, has no frequency ratios for
// the redo's set-up: every batch's traceback fails before the redo reads a match, a query without
// hits too (`LOSAT engine error: Traceback failed`, exit status 3; NCBI 2.17.0: rows
// db.x4b.ident.*). With one query per batch the first batch fails: nothing is reported.
#[test]
fn identity_fails_every_batch_at_the_redo_set_up() {
    let subject = input("db_search/tslice.fna");
    for (query, env) in [
        ("db_search/tslice_q.faa", &[][..]),
        ("db_search/tslice_nohit.faa", &[][..]),
        ("db_search/tslice_q.faa", &[("BATCH_SIZE", "1")][..]),
    ] {
        let query = input(query);
        for threads in [1, 4] {
            let run = tblastn(
                &[
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-matrix",
                    "IDENTITY",
                    "-use_sw_tback",
                    "-outfmt",
                    "6",
                ],
                threads,
                env,
            );
            assert_eq!(run.code, Some(3), "n{threads} {env:?}: {}", run.stderr);
            assert_eq!(run.stdout, "", "n{threads} {env:?}");
            assert_eq!(run.stderr, TRACEBACK_FAILED, "n{threads} {env:?}");
        }
    }
}
