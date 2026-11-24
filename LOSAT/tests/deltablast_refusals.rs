//! `losat deltablast`'s refusals (series PR, item PS13a): each runs where LOSAT's DELTABLAST stops, with its text, and
//! writes nothing to the standard output. NCBI BLAST+ 2.17.0 runs all of them (the rows of
//! `tests/db_search_manifest/deltablast.tsv` freeze its bytes for the cases that LOSAT writes); the refusals of a
//! tie of two hits and of a Seq-id that is not local or general are the unit tests of `deltablast::cdd_input`.

use std::path::Path;
use std::process::Command;

fn run(args: &[&str]) -> (i32, String, String) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args(["deltablast", "-query", "tests/db_search/delta/q_serpin.faa", "-db", "tests/db_search/delta/dbs"])
        .args(args)
        .output()
        .expect("runs losat");
    (
        output.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

const RPS: [&str; 2] = ["-rpsdb", "tests/db_search/delta/deltatest"];

fn refused(extra: &[&str], text: &str) {
    let mut args = RPS.to_vec();
    args.extend_from_slice(extra);
    let (rc, _, err) = run(&args);
    assert_ne!(rc, 0, "{extra:?}");
    assert!(err.contains(text), "{extra:?}: {err}");
    assert!(err.contains("LOSAT's DELTABLAST"), "{extra:?}: {err}");
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
    let (rc, _, err) = run(&["-rpsdb", "tests/db_search/delta/delta2"]);
    assert_ne!(rc, 0);
    assert!(err.contains("several volumes") && err.contains("LOSAT's DELTABLAST"), "{err}");
}

#[test]
fn a_subject_set_larger_than_max_target_seqs_is_refused() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args([
            "deltablast", "-query", "tests/db_search/delta/q_serpin.faa", "-subject",
            "tests/db_search/delta/subj_mix.faa", "-rpsdb", "tests/db_search/delta/deltatest",
            "-max_target_seqs", "2",
        ])
        .output()
        .expect("runs losat");
    assert!(!output.status.success());
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(err.contains("-max_target_seqs") && err.contains("LOSAT's DELTABLAST"), "{err}");
}
