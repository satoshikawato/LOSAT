//! `losat deltablast`'s refusals (series PR, item PS13a): each runs where LOSAT's DELTABLAST stops, with its text,
//! and writes to the standard output only what NCBI BLAST+ 2.17.0 writes before the same step (nothing for all of
//! them here). NCBI runs the commands of the first group to their end (the rows of
//! `tests/db_search_manifest/deltablast.tsv` freeze its bytes for the cases that LOSAT writes); with those of the
//! second group (PS13a fix round 1) NCBI crashes or stops with a C++ exception. The refusal of a Seq-id that is
//! not local or general is a unit test of `rps_search::cdd_input`.

use std::path::{Path, PathBuf};
use std::process::Command;

const RPS: [&str; 2] = ["-rpsdb", "tests/db_search/delta/deltatest"];

struct Run {
    rc: i32,
    stdout: Vec<u8>,
    stderr: String,
}

fn run_query(query: &str, args: &[&str]) -> Run {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args([
            "deltablast",
            "-query",
            &format!("tests/db_search/delta/{query}"),
            "-db",
            "tests/db_search/delta/dbs",
        ])
        .args(args)
        .output()
        .expect("runs losat");
    Run {
        rc: output.status.code().unwrap_or(-1),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// Asserts the run refused with `text` and wrote nothing to the standard output.
fn assert_refused(run: &Run, what: &str, text: &str) {
    assert_ne!(run.rc, 0, "{what}");
    assert!(run.stderr.contains(text), "{what}: {}", run.stderr);
    assert!(
        run.stderr
            .contains("is not supported by LOSAT's DELTABLAST"),
        "{what}: {}",
        run.stderr
    );
    assert!(
        run.stdout.is_empty(),
        "{what}: stdout {:?}",
        String::from_utf8_lossy(&run.stdout)
    );
}

fn refused(extra: &[&str], text: &str) {
    let mut args = RPS.to_vec();
    args.extend_from_slice(extra);
    assert_refused(
        &run_query("q_serpin.faa", &args),
        &format!("{extra:?}"),
        text,
    );
}

#[test]
fn show_domain_hits_is_refused() {
    refused(&["-show_domain_hits"], "-show_domain_hits");
}

#[test]
fn query_loc_is_refused() {
    refused(&["-query_loc", "10-200"], "query_loc");
}

#[test]
fn outfmts_other_than_0_6_7_are_refused() {
    for format in ["5", "8", "11", "15"] {
        refused(&["-outfmt", format], &format!("output format {format}"));
    }
}

#[test]
fn use_sw_tback_is_refused_as_psiblasts_pssm_search_refuses_it() {
    refused(&["-use_sw_tback"], "use_sw_tback");
}

#[test]
fn a_domain_database_of_several_volumes_is_refused() {
    let run = run_query("q_serpin.faa", &["-rpsdb", "tests/db_search/delta/delta2"]);
    assert_refused(&run, "delta2", "several volumes");
}

#[test]
fn a_subject_set_larger_than_max_target_seqs_is_refused_before_any_output() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args([
            "deltablast",
            "-query",
            "tests/db_search/delta/q_serpin.faa",
            "-subject",
            "tests/db_search/delta/subj_mix.faa",
            "-rpsdb",
            "tests/db_search/delta/deltatest",
            "-max_target_seqs",
            "2",
        ])
        .output()
        .expect("runs losat");
    let run = Run {
        rc: output.status.code().unwrap_or(-1),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    };
    assert_refused(&run, "subject set", "-max_target_seqs");
}

// --- PS13a fix round 1: NCBI 2.17.0 probes of `PR/ps13/fix1_probe` (x2, xs, none1) ---

#[test]
fn a_query_of_x_with_other_queries_in_a_one_round_batch_is_refused() {
    // NCBI's deltablast ends with SIGSEGV (exit 139) in the batch's domain search, before any output, whether
    // the X query comes after a query with letters (probe x2) or before it (probe xs).
    for query in ["q_x2.faa", "q_xs.faa"] {
        let mut args = RPS.to_vec();
        args.extend_from_slice(&["-outfmt", "6"]);
        assert_refused(
            &run_query(query, &args),
            query,
            "no Karlin-Altschul parameters",
        );
    }
}

struct TempDir(PathBuf);

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn save_last_with_one_round_and_a_first_query_without_alignments_is_refused_after_the_cdd_pssm_files(
) {
    // NCBI writes the batch's CDD PSSM files and its (empty) outfmt 6 report, then stops with an exception for a
    // null pointer (exit 255) where it would compute the last PSSM from the query's alignments (probe none1).
    let temp = TempDir(
        std::env::temp_dir().join(format!("losat-deltablast-none1-{}", std::process::id())),
    );
    std::fs::create_dir_all(&temp.0).unwrap();
    let ascii = temp.0.join("a");
    let asn = temp.0.join("s");
    let mut args = RPS.to_vec();
    args.extend_from_slice(&["-outfmt", "6", "-save_pssm_after_last_round"]);
    args.extend_from_slice(&[
        "-out_ascii_pssm",
        ascii.to_str().unwrap(),
        "-out_pssm",
        asn.to_str().unwrap(),
    ]);
    let run = run_query("q_none.faa", &args);
    assert_refused(&run, "none1", "-save_pssm_after_last_round with one round");
    for file in [&ascii, &asn] {
        assert!(
            std::fs::metadata(file).is_ok_and(|meta| meta.len() > 0),
            "{} is written",
            file.display()
        );
    }
}

#[cfg(target_os = "linux")]
#[test]
fn a_pssm_file_that_is_the_standard_output_is_refused_in_deltablasts_name() {
    for option in ["-out_pssm", "-out_ascii_pssm"] {
        let mut args = RPS.to_vec();
        args.extend_from_slice(&[option, "/dev/stdout"]);
        let text = format!("an {option} file that is also the standard output (/dev/stdout)");
        assert_refused(&run_query("q_serpin.faa", &args), option, &text);
    }
}
