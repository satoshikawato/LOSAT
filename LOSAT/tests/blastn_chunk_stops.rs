//! BLASTN query-chunk stops and refusal order (line RG item X-chunk):
//! - `CQuerySplitter::Split` (`split_query_cxx.cpp:234-249,876-880`): when any query chunk holds
//!   no query (positive pair with evened chunk size below or equal to overlap at L <= overlap;
//!   negative pair whose wrapped chunk size covers the batch; or negative pair whose chunk start
//!   reaches or passes the batch end before the last chunk), NCBI stops with `eNullPtr` (exit 3)
//!   before chunk masking, chunk setup (`e_nF`, `avg2p32`), database reads (`nonsq`), or empty-DB
//!   completion (`e_nA`, `sl_tfn_all.txt`), and after earlier batches' reports and the prologs
//!   that come with them (`0`, `15`, `16`, `20`; `0` bytes for `5`, `6`, `7`, `10`, `17` on batch
//!   1). LOSAT refuses there with `ErrorAfterBatchReports`, matching streamed and off paths.
//! - Negative pair with gaps where every chunk holds a query (PD-LOSAT-NCBI-DEFECTS D13): refused
//!   at the first valid chunk with subjects in hand, after that chunk's setup checks and sequence
//!   file reads; without subjects in hand or without a valid chunk (`gapvalid`), reported as NCBI.
//! - Approved Exception 1 (`CHUNK_SIZE` 101-150 with default overlap 100 and L >= 101): every
//!   chunk holds a query, and LOSAT searches each chunk once (exit 0) at 1 and 4 threads.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

// Matches `algorithm::blastn::query_split::{empty_chunk_refusal, negative_split_refusal}`
// (`algorithm` is `pub(crate)`, so integration tests assert the rendered message directly).
const EMPTY_CHUNK: &str = "a query chunk without a query is not supported by LOSAT's BLASTN";
const NEGATIVE_GAPS: &str =
    "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this \
query batch with NCBI BLAST+'s chunk ranges wrapped, which is not supported by LOSAT's BLASTN";

struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastn-chunk-stops-{tag}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn manifest_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn run_blastn(args: &[&str], env: &[(&str, &str)], cwd: Option<&Path>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
    cmd.arg("blastn").args(args);
    cmd.env("LOSAT_BLASTX_WORKER_LOG", "/dev/null");
    cmd.env_remove("CHUNK_SIZE");
    cmd.env_remove("OVERLAP_CHUNK_SIZE");
    cmd.env_remove("BATCH_SIZE");
    cmd.env_remove("BLAST_USAGE_REPORT");
    cmd.env_remove("NCBI_CONFIG_PATH");
    cmd.env_remove("NCBI_DONT_USE_NCBIRC");
    cmd.env_remove("NCBI_DONT_USE_LOCAL_CONFIG");
    cmd.env_remove("WINDOW_MASKER_PATH");
    cmd.env_remove("NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH");
    cmd.env_remove("LOSAT_X_RPTSTREAM");
    cmd.env_remove("LOSAT_X_STATS");
    for (k, v) in env {
        cmd.env(k, v);
    }
    if let Some(dir) = cwd {
        cmd.current_dir(dir);
    }
    cmd.output().expect("run losat blastn")
}

/// Checks that `LOSAT_X_RPTSTREAM=0` and `=1` produce identical exit status, stdout, and stderr,
/// and returns the off-path output.
fn run_both_streams(args: &[&str], env: &[(&str, &str)]) -> Output {
    let mut off_env = env.to_vec();
    off_env.push(("LOSAT_X_RPTSTREAM", "0"));
    let off = run_blastn(args, &off_env, None);

    let mut on_env = env.to_vec();
    on_env.push(("LOSAT_X_RPTSTREAM", "1"));
    let on = run_blastn(args, &on_env, None);

    assert_eq!(
        off.status.code(),
        on.status.code(),
        "stream exit mismatch for args={args:?} env={env:?}"
    );
    assert_eq!(
        off.stdout, on.stdout,
        "stream stdout mismatch for args={args:?} env={env:?}"
    );
    assert_eq!(
        off.stderr, on.stderr,
        "stream stderr mismatch for args={args:?} env={env:?}"
    );
    off
}

#[test]
fn empty_chunk_refusal_matches_across_formats_streams_and_threads() {
    let db = manifest_path("tests/db_search/taxfilter/tfn");
    let subj = manifest_path("tests/db_search/taxfilter/tfn.fna");
    let q_first = manifest_path("tests/fasta/xchunk/q_first.fna");
    let q_under1 = manifest_path("tests/fasta/xchunk/q_under1.fna");
    let q_under = manifest_path("tests/fasta/xchunk/q_under.fna");
    let q_neg_empty = manifest_path("tests/fasta/xchunk/q_neg_empty.fna");
    let q_neg201 = manifest_path("tests/fasta/xchunk/q_neg201.fna");
    let db_s = db.to_str().unwrap();
    let subj_s = subj.to_str().unwrap();
    let q_first_s = q_first.to_str().unwrap();
    let q_under1_s = q_under1.to_str().unwrap();
    let q_under_s = q_under.to_str().unwrap();
    let q_neg_empty_s = q_neg_empty.to_str().unwrap();
    let q_neg201_s = q_neg201.to_str().unwrap();

    // Batch 1 empty chunk (CHUNK_SIZE=140 on L=90): prolog bytes for 0, 15 (20 B), 16 (275 B),
    // 20 (87 B header); 0 bytes for 5, 6, 7, 10, 17; identical across streams and 1/4 threads.
    for (flag, target) in [("-db", db_s), ("-subject", subj_s)] {
        for (outfmt, min_bytes, exact_bytes) in [
            ("0", 100usize, None),
            ("5", 0, Some(0)),
            ("6", 0, Some(0)),
            ("7", 0, Some(0)),
            ("10", 0, Some(0)),
            ("15", 20, Some(20)),
            ("16", 275, Some(275)),
            ("17", 0, Some(0)),
            ("20", 87, Some(87)),
        ] {
            let out1 = run_both_streams(
                &[
                    "-task",
                    "blastn",
                    flag,
                    target,
                    "-query",
                    q_under1_s,
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    "1",
                ],
                &[("CHUNK_SIZE", "140")],
            );
            assert_eq!(out1.status.code(), Some(1), "{flag} outfmt={outfmt}");
            let err = String::from_utf8_lossy(&out1.stderr);
            assert!(
                err.contains(EMPTY_CHUNK),
                "{flag} outfmt={outfmt} stderr={err}"
            );
            if let Some(exact) = exact_bytes {
                assert_eq!(
                    out1.stdout.len(),
                    exact,
                    "{flag} outfmt={outfmt} stdout len"
                );
            } else {
                assert!(out1.stdout.len() >= min_bytes, "{flag} outfmt={outfmt}");
            }
            let out4 = run_blastn(
                &[
                    "-task",
                    "blastn",
                    flag,
                    target,
                    "-query",
                    q_under1_s,
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    "4",
                ],
                &[("CHUNK_SIZE", "140")],
                None,
            );
            assert_eq!(out1.status.code(), out4.status.code());
            assert_eq!(out1.stdout, out4.stdout, "{flag} outfmt={outfmt} t1 vs t4");
            assert_eq!(out1.stderr, out4.stderr, "{flag} outfmt={outfmt} t1 vs t4");
        }
    }

    // Batch 2 empty chunk across all three shapes (evened < overlap, negative wrap, unwrapped
    // negative at L=201): writes batch 1's report without epilog, then refuses.
    for (label, query_s, env) in [
        (
            "under",
            q_under_s,
            &[("BATCH_SIZE", "1"), ("CHUNK_SIZE", "140")][..],
        ),
        (
            "neg_wrap",
            q_neg_empty_s,
            &[
                ("BATCH_SIZE", "1"),
                ("CHUNK_SIZE", "-100"),
                ("OVERLAP_CHUNK_SIZE", "-200"),
            ][..],
        ),
        (
            "neg201",
            q_neg201_s,
            &[
                ("BATCH_SIZE", "1"),
                ("CHUNK_SIZE", "-132"),
                ("OVERLAP_CHUNK_SIZE", "-200"),
            ][..],
        ),
    ] {
        for outfmt in ["0", "5", "6", "7", "10", "15", "16", "17", "20"] {
            let out = run_both_streams(
                &[
                    "-task", "blastn", "-db", db_s, "-query", query_s, "-outfmt", outfmt,
                ],
                env,
            );
            assert_eq!(out.status.code(), Some(1), "{label} outfmt={outfmt}");
            let err = String::from_utf8_lossy(&out.stderr);
            assert!(
                err.contains(EMPTY_CHUNK),
                "{label} outfmt={outfmt} stderr={err}"
            );
            assert!(
                !out.stdout.is_empty(),
                "{label} outfmt={outfmt} must keep batch 1 report"
            );
            if outfmt == "6" || outfmt == "10" || outfmt == "20" {
                let first_only = run_blastn(
                    &[
                        "-task", "blastn", "-db", db_s, "-query", q_first_s, "-outfmt", outfmt,
                    ],
                    &[],
                    None,
                );
                assert_eq!(first_only.status.code(), Some(0));
                assert_eq!(
                    out.stdout, first_only.stdout,
                    "{label} outfmt={outfmt} batch 1 prefix"
                );
            }
            if outfmt == "7" {
                let text = String::from_utf8_lossy(&out.stdout);
                assert!(text.contains("# Query: q1"), "{label}: {text}");
                assert!(!text.contains("# BLAST processed"), "{label}: {text}");
            }
        }
    }

    // Batch 2 -subject outfmt 15 and 16: NCBI holds back the previous query, so only the 20 B /
    // 275 B prolog is written before batch 2's stop.
    for (outfmt, expected_len) in [("15", 20usize), ("16", 275usize)] {
        let out = run_both_streams(
            &[
                "-task", "blastn", "-subject", subj_s, "-query", q_under_s, "-outfmt", outfmt,
            ],
            &[("BATCH_SIZE", "1"), ("CHUNK_SIZE", "140")],
        );
        assert_eq!(out.status.code(), Some(1));
        assert_eq!(out.stdout.len(), expected_len, "subject b2 outfmt={outfmt}");
    }

    // File sink (-out file) on batch 1 writes the same prolog bytes as stdout.
    let dir = TempDir::new("filesink");
    for outfmt in ["0", "15", "16", "20"] {
        let out_file = dir.path(&format!("b1.{outfmt}.out"));
        let stdout_run = run_blastn(
            &[
                "-task", "blastn", "-db", db_s, "-query", q_under1_s, "-outfmt", outfmt,
            ],
            &[("CHUNK_SIZE", "140")],
            None,
        );
        let file_run = run_blastn(
            &[
                "-task",
                "blastn",
                "-db",
                db_s,
                "-query",
                q_under1_s,
                "-outfmt",
                outfmt,
                "-out",
                out_file.to_str().unwrap(),
            ],
            &[("CHUNK_SIZE", "140")],
            None,
        );
        assert_eq!(file_run.status.code(), Some(1));
        assert!(file_run.stdout.is_empty());
        assert_eq!(
            std::fs::read(&out_file).unwrap(),
            stdout_run.stdout,
            "file sink outfmt={outfmt}"
        );
    }
}

#[test]
fn resplit_exception_one_boundary_and_negative_gaps_match_across_tasks() {
    let db = manifest_path("tests/db_search/taxfilter/tfn");
    let q_l100 = manifest_path("tests/fasta/xchunk/q_l100.fna");
    let q_l101 = manifest_path("tests/fasta/xchunk/q_l101.fna");
    let q_neg_gap = manifest_path("tests/fasta/xchunk/q_neg_gap.fna");
    let q_first = manifest_path("tests/fasta/xchunk/q_first.fna");
    let db_s = db.to_str().unwrap();
    let q_l100_s = q_l100.to_str().unwrap();
    let q_l101_s = q_l101.to_str().unwrap();
    let q_neg_gap_s = q_neg_gap.to_str().unwrap();
    let q_first_s = q_first.to_str().unwrap();

    // CHUNK_SIZE=150 (default overlap 100): L=100 leaves chunk 1 empty (refused), while L=101 has
    // two non-empty chunks (Approved Exception 1: searched once, exit 0 at 1 and 4 threads).
    for task in ["blastn", "megablast", "dc-megablast", "blastn-short"] {
        let stop = run_both_streams(
            &[
                "-task", task, "-db", db_s, "-query", q_l100_s, "-outfmt", "7",
            ],
            &[("BATCH_SIZE", "1"), ("CHUNK_SIZE", "150")],
        );
        assert_eq!(stop.status.code(), Some(1), "task={task} L=100");
        assert!(String::from_utf8_lossy(&stop.stderr).contains(EMPTY_CHUNK));

        let ok1 = run_both_streams(
            &[
                "-task",
                task,
                "-db",
                db_s,
                "-query",
                q_l101_s,
                "-outfmt",
                "7",
                "-num_threads",
                "1",
            ],
            &[("BATCH_SIZE", "1"), ("CHUNK_SIZE", "150")],
        );
        assert_eq!(ok1.status.code(), Some(0), "task={task} L=101");
        let ok4 = run_blastn(
            &[
                "-task",
                task,
                "-db",
                db_s,
                "-query",
                q_l101_s,
                "-outfmt",
                "7",
                "-num_threads",
                "4",
            ],
            &[("BATCH_SIZE", "1"), ("CHUNK_SIZE", "150")],
            None,
        );
        assert_eq!(ok1.stdout, ok4.stdout, "task={task} L=101 t1 vs t4");

        // Negative pair with gaps on batch 2: writes batch 1's report and refuses with D13 text.
        let gap = run_both_streams(
            &[
                "-task",
                task,
                "-db",
                db_s,
                "-query",
                q_neg_gap_s,
                "-outfmt",
                "6",
            ],
            &[
                ("BATCH_SIZE", "1"),
                ("CHUNK_SIZE", "-100"),
                ("OVERLAP_CHUNK_SIZE", "-200"),
            ],
        );
        assert_eq!(gap.status.code(), Some(1), "task={task} neg_gap");
        assert!(String::from_utf8_lossy(&gap.stderr).contains(NEGATIVE_GAPS));
        let anchor = run_blastn(
            &[
                "-task", task, "-db", db_s, "-query", q_first_s, "-outfmt", "6",
            ],
            &[],
            None,
        );
        assert_eq!(gap.stdout, anchor.stdout, "task={task} neg_gap batch 1");
    }
}

#[test]
fn window_masker_registry_loads_follow_ncbis_split_step_order() {
    // With a malformed `.ncbirc` and `-window_masker_taxid`:
    // - startup posts the registry error twice;
    // - batch 1's query masking posts it once (3 total);
    // - an empty-chunk stop (`CQuerySplitter::Split`) happens before any chunk is masked, so it
    //   stays at 3 posts;
    // - a negative-pair gap refusal happens inside chunk 0 after chunk 0's mask load, so it posts
    //   a 4th time.
    let dir = TempDir::new("wm");
    let db = manifest_path("tests/db_search/taxfilter/tfn");
    let q_under1 = manifest_path("tests/fasta/xchunk/q_under1.fna");
    let q_neg_gap1 = manifest_path("tests/fasta/xchunk/q_neg_gap1.fna");
    let counts_src = manifest_path("tests/windowmasker/counts/u4.obinary");
    let wm_dir = dir.path("wm");
    std::fs::create_dir_all(wm_dir.join("9606")).unwrap();
    std::fs::copy(&counts_src, wm_dir.join("9606/wmasker.obinary")).unwrap();
    std::fs::write(
        dir.path(".ncbirc"),
        format!(
            "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = {}\nmalformed line without equals\n",
            wm_dir.display()
        ),
    )
    .unwrap();

    let count_posts = |stderr: &[u8]| -> usize {
        String::from_utf8_lossy(stderr)
            .matches("CNcbiRegistry: Syntax error")
            .count()
    };

    let empty = run_blastn(
        &[
            "-task",
            "blastn",
            "-db",
            db.to_str().unwrap(),
            "-query",
            q_under1.to_str().unwrap(),
            "-window_masker_taxid",
            "9606",
            "-outfmt",
            "7",
        ],
        &[("CHUNK_SIZE", "140")],
        Some(&dir.0),
    );
    assert_eq!(empty.status.code(), Some(1));
    let empty_err = String::from_utf8_lossy(&empty.stderr);
    assert!(empty_err.contains(EMPTY_CHUNK), "{empty_err}");
    assert_eq!(
        count_posts(&empty.stderr),
        3,
        "empty chunk stops before chunk 0 masking: {empty_err}"
    );

    let gap = run_blastn(
        &[
            "-task",
            "blastn",
            "-db",
            db.to_str().unwrap(),
            "-query",
            q_neg_gap1.to_str().unwrap(),
            "-window_masker_taxid",
            "9606",
            "-outfmt",
            "7",
        ],
        &[("CHUNK_SIZE", "-100"), ("OVERLAP_CHUNK_SIZE", "-200")],
        Some(&dir.0),
    );
    assert_eq!(gap.status.code(), Some(1));
    let gap_err = String::from_utf8_lossy(&gap.stderr);
    assert!(gap_err.contains(NEGATIVE_GAPS), "{gap_err}");
    assert_eq!(
        count_posts(&gap.stderr),
        4,
        "negative gap refuses after chunk 0 masking: {gap_err}"
    );

    // `q_gapvalid` (no chunk has a valid context, only the gap does):
    // - on `e_nF` (OID 0 read fails in `BLAST_GapAlignSetUp` before context check), chunk 0 stops
    //   before chunk 1 is masked (4 posts total, matching NCBI 2.17.0);
    // - on `avg2p32` (average subject length check is reached only after context check), both
    //   chunks are masked and skipped, and traceback setup fails with `Traceback failed` (5 posts).
    let q_gapvalid = manifest_path("tests/fasta/xchunk/q_gapvalid.fna");
    let db_enf = manifest_path("tests/db_search/aliaskeys/e_nF");
    let db_avg = manifest_path("tests/db_search/dbcounts/avg2p32");

    let gv_enf = run_blastn(
        &[
            "-task",
            "blastn",
            "-db",
            db_enf.to_str().unwrap(),
            "-query",
            q_gapvalid.to_str().unwrap(),
            "-window_masker_taxid",
            "9606",
            "-outfmt",
            "7",
        ],
        &[("CHUNK_SIZE", "-100"), ("OVERLAP_CHUNK_SIZE", "-200")],
        Some(&dir.0),
    );
    assert_eq!(gv_enf.status.code(), Some(3));
    let gv_enf_err = String::from_utf8_lossy(&gv_enf.stderr);
    assert!(
        gv_enf_err.contains("Database memory map file error"),
        "{gv_enf_err}"
    );
    assert_eq!(
        count_posts(&gv_enf.stderr),
        4,
        "gapvalid e_nF stops in chunk 0 OID 0 read before chunk 1 masking: {gv_enf_err}"
    );

    let gv_avg = run_blastn(
        &[
            "-task",
            "blastn",
            "-db",
            db_avg.to_str().unwrap(),
            "-query",
            q_gapvalid.to_str().unwrap(),
            "-window_masker_taxid",
            "9606",
            "-outfmt",
            "7",
        ],
        &[("CHUNK_SIZE", "-100"), ("OVERLAP_CHUNK_SIZE", "-200")],
        Some(&dir.0),
    );
    assert_eq!(gv_avg.status.code(), Some(3));
    let gv_avg_err = String::from_utf8_lossy(&gv_avg.stderr);
    assert!(gv_avg_err.contains("Traceback failed"), "{gv_avg_err}");
    assert_eq!(
        count_posts(&gv_avg.stderr),
        5,
        "gapvalid avg2p32 masks both chunks before traceback setup fails: {gv_avg_err}"
    );

    // Refused split duplicate-id batch (`q_dup300.fna`, `CHUNK_SIZE=150 OVERLAP_CHUNK_SIZE=10`):
    // the refused duplicate-id batch loads chunk 0's registry as NCBI does, so it posts 4 times
    // before emitting NCBI's error on `avg2p32`, `nonsq`, and `e_nF`, and posts 4 times before
    // LOSAT's duplicate-id refusal on `tfn`.
    let q_dup300 = manifest_path("tests/fasta/xchunk/q_dup300.fna");
    let db_nonsq = manifest_path("tests/db_search/arch_tbn_snucl_nonsq");
    for (tag, db_path, expected_rc, expected_err) in [
        (
            "avg2p32",
            &db_avg,
            3,
            "The average subject length is too short",
        ),
        ("nonsq", &db_nonsq, 3, "Database memory map file error"),
        ("e_nF", &db_enf, 3, "Database memory map file error"),
        ("tfn", &db, 1, "two query records have the Seq-id lcl|dup"),
    ] {
        let out = run_blastn(
            &[
                "-task",
                "blastn",
                "-parse_deflines",
                "-db",
                db_path.to_str().unwrap(),
                "-query",
                q_dup300.to_str().unwrap(),
                "-window_masker_taxid",
                "9606",
                "-outfmt",
                "7",
            ],
            &[("CHUNK_SIZE", "150"), ("OVERLAP_CHUNK_SIZE", "10")],
            Some(&dir.0),
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(
            out.status.code(),
            Some(expected_rc),
            "dup300 {tag}: {stderr}"
        );
        assert!(stderr.contains(expected_err), "dup300 {tag}: {stderr}");
        assert_eq!(
            count_posts(&out.stderr),
            4,
            "dup300 {tag} loads chunk 0 registry as NCBI does: {stderr}"
        );
    }
}
