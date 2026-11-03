//! BLASTX's `-query_loc` and `-subject_loc` on the command line (line SX, sessions SXd1 and SXd2).
//!
//! That the reports equal NCBI's is the business of the fixtures: the `sxd.blastx.*` and
//! `sxd2.blastx.*` rows of `tests/outfmt0_manifest.tsv`, the `blastx.*` rows of
//! `tests/fixtures/range_regression/` and the `db.bxloc.*` rows of `tests/db_search_manifest.tsv`.
//! Here: what a run writes with a range does not depend on the threads or on the path its query
//! batches take (`LOSAT_X_BXFLOW`, the serial batches of `-parse_deflines`), also when a later
//! batch holds only skipped records (NCBI's `Empty CBlastQueryVector` after the reports of the
//! batches before); and that a range goes with every format but those NCBI's formatting handler
//! refuses.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fixture(path: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(path)
        .to_string_lossy()
        .into_owned()
}

fn temp_path(name: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat-blastx-range-{}-{}-{name}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

struct Run {
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    code: Option<i32>,
}

/// One `losat blastx` run with `args` and `env` as its only variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` is a file of its own: a binary built with the `blastx-worker-probe`
/// feature (`cargo test --all-features`) needs it as soon as a search has more than one thread.
fn blastx(args: &[&str], env: &[(&str, &str)]) -> Run {
    let log = temp_path("worker.log");
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastx")
        .args(args)
        .env("LOSAT_BLASTX_WORKER_LOG", &log)
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&log);
    Run {
        stdout: output.stdout,
        stderr: output.stderr,
        code: output.status.code(),
    }
}

/// Every ranged search below writes, at two to four threads, with and without the flow of query
/// batches, what it writes at one thread: standard output, standard error and the exit status.
/// The query batches (`BATCH_SIZE`) hold skipped records (their interval starts more than one
/// letter past their end), records without letters in the range and, in the failing search, a
/// last batch of skipped records only.
#[cfg(feature = "parallel")]
#[test]
fn ranged_runs_do_not_depend_on_the_threads_or_the_batch_flow() {
    let many = fixture("tests/fasta/outfmt0/sxd_bx_many_q.fna");
    let mje = fixture("tests/fasta/MjeNMV.faa");
    let unsearched = fixture("tests/fasta/outfmt0/x16_bx_q_unsearched.fna");
    let bslice = fixture("tests/db_search/bslice.faa");
    let pd_query = fixture("tests/fasta/outfmt0/pd_sub_q.fna");
    let pd_subject = fixture("tests/fasta/outfmt0/pd_sub_p.faa");
    let searches: [(&str, Vec<&str>, &str, bool); 7] = [
        (
            "many records, outfmt 7",
            vec![
                "-query",
                &many,
                "-subject",
                &mje,
                "-query_loc",
                "2001-6000",
                "-outfmt",
                "7",
            ],
            "2500",
            false,
        ),
        (
            "many records, minus strand, outfmt 0",
            vec![
                "-query",
                &many,
                "-subject",
                &mje,
                "-query_loc",
                "1500-6000",
                "-strand",
                "minus",
                "-subject_loc",
                "20-400",
            ],
            "3000",
            false,
        ),
        (
            "a last batch of skipped records, outfmt 0",
            vec![
                "-query",
                &unsearched,
                "-subject",
                &bslice,
                "-query_loc",
                "400-1000",
            ],
            "500",
            true,
        ),
        (
            "a last batch of skipped records, outfmt 7",
            vec![
                "-query",
                &unsearched,
                "-subject",
                &bslice,
                "-query_loc",
                "400-1000",
                "-outfmt",
                "7",
            ],
            "500",
            true,
        ),
        (
            "many records, outfmt 11 (an archive per batch)",
            vec![
                "-query",
                &many,
                "-subject",
                &mje,
                "-query_loc",
                "2001-6000",
                "-outfmt",
                "11",
            ],
            "2500",
            false,
        ),
        (
            "a last batch of skipped records, outfmt 5",
            vec![
                "-query",
                &unsearched,
                "-subject",
                &bslice,
                "-query_loc",
                "400-1000",
                "-outfmt",
                "5",
            ],
            "500",
            true,
        ),
        (
            "-parse_deflines (serial batches)",
            vec![
                "-query",
                &pd_query,
                "-subject",
                &pd_subject,
                "-parse_deflines",
                "-query_loc",
                "31-600",
                "-outfmt",
                "6",
            ],
            "1000",
            false,
        ),
    ];
    for (what, args, batch_size, fails) in &searches {
        let reference = blastx(args, &[("BATCH_SIZE", batch_size)]);
        if *fails {
            assert_eq!(reference.code, Some(3), "{what}");
            assert!(
                String::from_utf8_lossy(&reference.stderr).contains("Empty CBlastQueryVector"),
                "{what}"
            );
        } else {
            assert_eq!(
                reference.code,
                Some(0),
                "{what}: {}",
                String::from_utf8_lossy(&reference.stderr)
            );
        }
        assert!(!reference.stdout.is_empty(), "{what}: no report");
        for threads in ["2", "3", "4"] {
            for flow in ["0", "1"] {
                let mut threaded = args.clone();
                threaded.extend(["-num_threads", threads]);
                let run = blastx(
                    &threaded,
                    &[("BATCH_SIZE", batch_size), ("LOSAT_X_BXFLOW", flow)],
                );
                let context = format!("{what}, {threads} threads, LOSAT_X_BXFLOW={flow}");
                assert_eq!(run.code, reference.code, "{context}");
                assert!(run.stdout == reference.stdout, "{context}: standard output");
                assert_eq!(
                    String::from_utf8_lossy(&run.stderr),
                    String::from_utf8_lossy(&reference.stderr),
                    "{context}"
                );
            }
        }
    }
}

/// SAM, AIRR and FASTA stop with NCBI's errors of BLASTX's formatting handler, which comes after
/// the ranges (an empty `-out` file, as NCBI's); the other formats search the ranges (session
/// SXd2: their bytes are the `sxd2.blastx.*` rows of outfmt0_manifest.tsv, outfmt 13 and 14 with
/// their side files).
#[test]
fn a_range_goes_with_every_format_but_ncbis_errors() {
    let query = fixture("tests/db_search/bslice_q.fna");
    let subject = fixture("tests/db_search/bslice.faa");
    for format in ["5", "8", "9", "11", "12", "15", "16", "17", "19", "21"] {
        for range in ["-query_loc", "-subject_loc"] {
            let out = temp_path("out");
            let out_text = out.to_string_lossy().into_owned();
            let run = blastx(
                &[
                    "-query", &query, "-subject", &subject, range, "1-100", "-outfmt", format,
                    "-out", &out_text,
                ],
                &[],
            );
            let stderr = String::from_utf8_lossy(&run.stderr);
            let ncbi_error = match format {
                "17" => Some("SAM format is only applicable to blastn"),
                "19" => Some("AIRR rearrangement format is only applicable to igblastn"),
                "21" => Some("FASTA output format is only applicable to magicblast"),
                _ => None,
            };
            assert!(run.stdout.is_empty(), "{format} {range}");
            let written = std::fs::read(&out).ok();
            let _ = std::fs::remove_file(&out);
            if let Some(ncbi_error) = ncbi_error {
                assert_eq!(run.code, Some(1), "{format} {range}: {stderr}");
                assert!(stderr.contains(ncbi_error), "{format} {range}: {stderr}");
                // The handler of `-query` and `-out` comes before the formatting handler.
                assert_eq!(written, Some(Vec::new()), "{range}");
            } else {
                assert_eq!(run.code, Some(0), "{format} {range}: {stderr}");
                assert!(
                    written.is_some_and(|bytes| !bytes.is_empty()),
                    "{format} {range}: no report"
                );
            }
        }
    }
}
