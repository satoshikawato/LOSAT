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

// --- PR #1664 R2: error order (c1..c5), short/zero-count delta files, A3, D7 ---

fn run_raw(args: &[&str], envs: &[(&str, &str)]) -> Run {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
    cmd.current_dir(root)
        .env_remove("BLASTDB")
        .arg("deltablast")
        .args(args);
    for &(k, v) in envs {
        cmd.env(k, v);
    }
    let output = cmd.output().expect("runs losat");
    Run {
        rc: output.status.code().unwrap_or(-1),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

fn copy_db_without(src_dir: &Path, stem: &str, dst_dir: &Path, skip_ext: &str) -> PathBuf {
    std::fs::create_dir_all(dst_dir).unwrap();
    for entry in std::fs::read_dir(src_dir).unwrap() {
        let path = entry.unwrap().path();
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if let Some(ext) = name.strip_prefix(&format!("{stem}.")) {
            if ext != skip_ext {
                std::fs::copy(&path, dst_dir.join(name)).unwrap();
            }
        }
    }
    dst_dir.join(stem)
}

#[test]
fn first_batch_stops_check_rpsdb_and_domain_options_first() {
    let temp = TempDir(
        std::env::temp_dir().join(format!("losat-deltablast-order-{}", std::process::id())),
    );
    std::fs::create_dir_all(&temp.0).unwrap();
    let qe2 = temp.0.join("qe2.fa");
    let qe = temp.0.join("qe.fa");
    let q1 = temp.0.join("q1.fa");
    std::fs::write(&qe2, ">q1\nMKVLAAGIVALLLAAGCSSSKEETLKKAMDELTANG\n>empty\n\n").unwrap();
    std::fs::write(&qe, ">empty\n\n").unwrap();
    std::fs::write(&q1, ">q1\nMKVLAAGIVALLLAAGCSSSKEETLKKAMDELTANG\n").unwrap();
    let delta_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/delta");
    let nopsq = copy_db_without(&delta_dir, "dbs", &temp.0.join("nopsq"), "psq");

    // c1: multi-query batch with empty query + missing -rpsdb -> -rpsdb error (exit 2)
    let c1 = run_raw(
        &[
            "-db",
            "tests/db_search/delta/dbs",
            "-query",
            qe2.to_str().unwrap(),
            "-rpsdb",
            "nosuch",
        ],
        &[],
    );
    assert_eq!(c1.rc, 2, "c1: {}", c1.stderr);
    assert!(
        c1.stderr.contains("No alias or index file found for protein database [nosuch]"),
        "c1: {}",
        c1.stderr
    );
    assert!(!c1.stdout.is_empty(), "c1 should write outfmt 0 prolog");

    // c2: single empty query + missing -rpsdb -> -rpsdb error (exit 2)
    let c2 = run_raw(
        &[
            "-db",
            "tests/db_search/delta/dbs",
            "-query",
            qe.to_str().unwrap(),
            "-rpsdb",
            "nosuch",
        ],
        &[],
    );
    assert_eq!(c2.rc, 2, "c2: {}", c2.stderr);
    assert!(
        c2.stderr.contains("No alias or index file found for protein database [nosuch]"),
        "c2: {}",
        c2.stderr
    );
    assert!(!c2.stdout.is_empty(), "c2 should write outfmt 0 prolog");

    // c3: invalid CHUNK_SIZE + missing -rpsdb -> -rpsdb error (exit 2)
    let c3 = run_raw(
        &[
            "-db",
            "tests/db_search/delta/dbs",
            "-query",
            qe2.to_str().unwrap(),
            "-rpsdb",
            "nosuch",
        ],
        &[("CHUNK_SIZE", "abc")],
    );
    assert_eq!(c3.rc, 2, "c3: {}", c3.stderr);
    assert!(
        c3.stderr.contains("No alias or index file found for protein database [nosuch]"),
        "c3: {}",
        c3.stderr
    );

    // c4: missing main -db .psq + missing -rpsdb -> -rpsdb error (exit 2)
    let c4 = run_raw(
        &[
            "-db",
            nopsq.to_str().unwrap(),
            "-query",
            q1.to_str().unwrap(),
            "-rpsdb",
            "nosuch",
        ],
        &[],
    );
    assert_eq!(c4.rc, 2, "c4: {}", c4.stderr);
    assert!(
        c4.stderr.contains("No alias or index file found for protein database [nosuch]"),
        "c4: {}",
        c4.stderr
    );

    // c4 with valid -rpsdb -> main -db .psq error (exit 2) after domain search
    let c4_valid = run_raw(
        &[
            "-db",
            nopsq.to_str().unwrap(),
            "-query",
            q1.to_str().unwrap(),
            "-rpsdb",
            "tests/db_search/delta/deltatest",
        ],
        &[],
    );
    assert_eq!(c4_valid.rc, 2, "c4_valid: {}", c4_valid.stderr);
    assert!(
        c4_valid.stderr.contains("dbs.psq) not found"),
        "c4_valid: {}",
        c4_valid.stderr
    );

    // c5: batch with empty query + valid -rpsdb + -domain_inclusion_ethresh 0 -> options error (exit 1)
    let c5 = run_raw(
        &[
            "-db",
            "tests/db_search/delta/dbs",
            "-query",
            qe2.to_str().unwrap(),
            "-rpsdb",
            "tests/db_search/delta/deltatest",
            "-domain_inclusion_ethresh",
            "0",
        ],
        &[],
    );
    assert_eq!(c5.rc, 1, "c5: {}", c5.stderr);
    assert!(
        c5.stderr
            .contains("expect value or cutoff score must be greater than zero"),
        "c5: {}",
        c5.stderr
    );
    assert!(!c5.stdout.is_empty(), "c5 should write outfmt 0 prolog");
}

#[test]
fn short_and_corrupt_wcounts_and_obsr_match_ncbi_open_and_refuse_on_hit() {
    let temp = TempDir(
        std::env::temp_dir().join(format!("losat-deltablast-cdfiles-{}", std::process::id())),
    );
    let delta_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/delta");
    let orig_wcounts = std::fs::read(delta_dir.join("deltatest.wcounts")).unwrap();
    let orig_obsr = std::fs::read(delta_dir.join("deltatest.obsr")).unwrap();

    // m4: 4-byte valid magic .wcounts -> no-hit query exits 0; hit query refuses cleanly.
    let m4 = copy_db_without(&delta_dir, "deltatest", &temp.0.join("m4"), "");
    std::fs::write(temp.0.join("m4/deltatest.wcounts"), &orig_wcounts[..4]).unwrap();
    let run = run_query("q_none.faa", &["-rpsdb", m4.to_str().unwrap(), "-outfmt", "6"]);
    assert_eq!(run.rc, 0, "m4 no-hit: {}", run.stderr);
    let run = run_query(
        "q_serpin.faa",
        &["-rpsdb", m4.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_refused(&run, "m4 hit", "shorter than its header");

    // h20: 20-byte truncated header .wcounts -> no-hit exits 0; hit refuses.
    let h20 = copy_db_without(&delta_dir, "deltatest", &temp.0.join("h20"), "");
    std::fs::write(temp.0.join("h20/deltatest.wcounts"), &orig_wcounts[..20]).unwrap();
    let run = run_query("q_none.faa", &["-rpsdb", h20.to_str().unwrap(), "-outfmt", "6"]);
    assert_eq!(run.rc, 0, "h20 no-hit: {}", run.stderr);
    let run = run_query(
        "q_serpin.faa",
        &["-rpsdb", h20.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_refused(&run, "h20 hit", "shorter than its header");

    // o4: 4-byte valid magic .obsr -> no-hit exits 0; hit refuses.
    let o4 = copy_db_without(&delta_dir, "deltatest", &temp.0.join("o4"), "");
    std::fs::write(temp.0.join("o4/deltatest.obsr"), &orig_obsr[..4]).unwrap();
    let run = run_query("q_none.faa", &["-rpsdb", o4.to_str().unwrap(), "-outfmt", "6"]);
    assert_eq!(run.rc, 0, "o4 no-hit: {}", run.stderr);
    let run = run_query(
        "q_serpin.faa",
        &["-rpsdb", o4.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_refused(&run, "o4 hit", "shorter than its header");

    // bad4: 4-byte wrong magic .wcounts -> NCBI corrupt profile error (exit 3).
    let bad4 = copy_db_without(&delta_dir, "deltatest", &temp.0.join("bad4"), "");
    std::fs::write(temp.0.join("bad4/deltatest.wcounts"), [0u8; 4]).unwrap();
    let run = run_query(
        "q_none.faa",
        &["-rpsdb", bad4.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_eq!(run.rc, 3, "bad4: {}", run.stderr);
    assert!(
        run.stderr
            .contains("is either corrupt or constructed for an incompatible architecture"),
        "bad4: {}",
        run.stderr
    );

    // e0: empty .wcounts -> refused before NCBI crashes.
    let e0 = copy_db_without(&delta_dir, "deltatest", &temp.0.join("e0"), "");
    std::fs::write(temp.0.join("e0/deltatest.wcounts"), []).unwrap();
    let run = run_query("q_none.faa", &["-rpsdb", e0.to_str().unwrap(), "-outfmt", "6"]);
    assert_refused(&run, "e0", "shorter than its 4-byte magic number");

    // Zero-count .wcounts column: refused on a non-X query residue, succeeds on an X query residue.
    let zc = copy_db_without(&delta_dir, "deltatest", &temp.0.join("zc"), "");
    let mut zbytes = orig_wcounts.clone();
    let num_profiles = i32::from_le_bytes(zbytes[4..8].try_into().unwrap()) as usize;
    let data_start = 8 + 4 * (num_profiles + 1);
    // Zero out every column of every profile's data so any hit hits a zero-count column.
    for b in &mut zbytes[data_start..] {
        *b = 0;
    }
    std::fs::write(temp.0.join("zc/deltatest.wcounts"), &zbytes).unwrap();
    let run = run_query(
        "q_serpin.faa",
        &["-rpsdb", zc.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_refused(&run, "zero-count non-X", "column of zero counts");
}

#[test]
fn domain_search_refusal_and_subject_limit_by_num_descriptions_name_deltablast() {
    let temp = TempDir(
        std::env::temp_dir().join(format!("losat-deltablast-a3d7-{}", std::process::id())),
    );
    let delta_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/delta");
    // A3: truncate .loo of a copy of deltatest so domain search refuses in DELTABLAST's name.
    let badloo = copy_db_without(&delta_dir, "deltatest", &temp.0.join("badloo"), "");
    std::fs::write(temp.0.join("badloo/deltatest.loo"), [0u8; 4]).unwrap();
    let run = run_query(
        "q_serpin.faa",
        &["-rpsdb", badloo.to_str().unwrap(), "-outfmt", "6"],
    );
    assert_refused(&run, "A3 domain refusal", "DELTABLAST");
    assert!(
        !run.stderr.contains("RPSBLAST"),
        "A3 must not name RPSBLAST: {}",
        run.stderr
    );

    // D7: -subject refusal when the limit comes from -num_descriptions / -num_alignments.
    let run = run_raw(
        &[
            "-query",
            "tests/db_search/delta/q_serpin.faa",
            "-subject",
            "tests/db_search/delta/subj_mix.faa",
            "-rpsdb",
            "tests/db_search/delta/deltatest",
            "-num_descriptions",
            "1",
            "-num_alignments",
            "1",
        ],
        &[],
    );
    assert_refused(&run, "D7 num_descriptions/num_alignments", "-num_descriptions");
}
