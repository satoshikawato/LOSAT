//! The two-hit seed stage of dc-megablast on the search pool (`LOSAT_X_PAIRPARDISC`,
//! `src/algorithm/blastn/blast_engine/x_pair_seed_disc.rs`) on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses between
//! two LOSAT code paths with the same result: the subject scan cut into strips that the pool
//! threads scan, and the diagonal-hash fold done by the owner thread in the scan order, or the
//! scan and fold of the reference code. So the property under test is that the report does not
//! depend on it: for every discontiguous template (word size 11 and 12; coding, optimal, and
//! both templates; lengths 16, 18 and 21) the search with the switch at 1, 2 and 4 threads
//! gives the bytes of the search without it, and the shadow run
//! (`LOSAT_X_PAIRPARDISCSHADOW`, which computes the stage on copies and compares the hit list
//! and the whole diagonal hash with the reference code for every chunk) gives them too and
//! compares at least one chunk. That those bytes equal NCBI's is the business of the dc-megablast
//! fixtures (`tests/blastn_regression_fixtures.py`, `dc.*`).
//!
//! The input is one subject (the stage takes a search of one subject by the owner of a pool of
//! two threads or more) and one query of 287 kb, whose query block is longer than 8000 letters
//! (the stage takes the diagonal hash only).

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH: (&str, &str) = ("LOSAT_X_PAIRPARDISC", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARDISCSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One search of LC738874 (query) against LC738868 (subject), `-outfmt 6` on standard output,
/// with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing. Returns the report and the standard error.
fn blastn(task_args: &[&str], threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(task_args)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task_args:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For one template: the search without the switch at one thread, then with the switch at 1, 2
/// and 4 threads (the same bytes, nothing on standard error), then the shadow at 4 threads (the
/// same bytes, and its summary line with at least one chunk compared).
fn check_template(template: &[&str]) {
    let mut task_args = vec!["-task", "dc-megablast"];
    task_args.extend_from_slice(template);
    let (expected, expected_stderr) = blastn(&task_args, 1, &[]);
    assert!(!expected.is_empty(), "{template:?}: no hits");
    assert_eq!(expected_stderr, "", "{template:?}");
    for threads in [1, 2, 4] {
        let (report, stderr) = blastn(&task_args, threads, &[SWITCH]);
        assert!(
            report == expected,
            "{template:?} threads {threads}: the report differs with {SWITCH:?}"
        );
        assert_eq!(stderr, "", "{template:?} threads {threads}");
    }
    let (report, stderr) = blastn(&task_args, 4, &[SHADOW]);
    assert!(
        report == expected,
        "{template:?}: the report differs with {SHADOW:?}"
    );
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_PAIRPARDISCSHADOW] blastn: seed chunks equal="))
        .unwrap_or_else(|| panic!("{template:?}: no summary line:\n{stderr}"));
    assert!(
        !line.starts_with("[X_PAIRPARDISCSHADOW] blastn: seed chunks equal=0 "),
        "{template:?}: {line}"
    );
    assert_eq!(stderr.lines().count(), 1, "{template:?}: {stderr}");
}

// The default template of dc-megablast (word size 11, coding, length 18) and the other one-
// template coding cases: each length picks its own scan kernel (11-of-18 and 11-of-21 coding
// have their own, the rest the general one).
#[test]
fn coding_templates_give_the_same_report() {
    check_template(&[]);
    check_template(&[
        "-word_size",
        "11",
        "-template_type",
        "coding",
        "-template_length",
        "16",
    ]);
    check_template(&[
        "-word_size",
        "11",
        "-template_type",
        "coding",
        "-template_length",
        "21",
    ]);
    check_template(&[
        "-word_size",
        "12",
        "-template_type",
        "coding",
        "-template_length",
        "21",
    ]);
}

#[test]
fn optimal_templates_give_the_same_report() {
    check_template(&[
        "-word_size",
        "11",
        "-template_type",
        "optimal",
        "-template_length",
        "21",
    ]);
    check_template(&[
        "-word_size",
        "12",
        "-template_type",
        "optimal",
        "-template_length",
        "16",
    ]);
}

// Two templates (coding_and_optimal): per subject offset the first template's chain, then the
// second's (the two-template scan kernel).
#[test]
fn two_templates_give_the_same_report() {
    check_template(&[
        "-word_size",
        "11",
        "-template_type",
        "coding_and_optimal",
        "-template_length",
        "16",
    ]);
    check_template(&[
        "-word_size",
        "11",
        "-template_type",
        "coding_and_optimal",
        "-template_length",
        "18",
    ]);
    check_template(&[
        "-word_size",
        "12",
        "-template_type",
        "coding_and_optimal",
        "-template_length",
        "21",
    ]);
}

// A contiguous table in two-hit mode is not taken: the CLI has no way to ask for it
// (`-window_size` is an unported option), and a megablast search (one hit) with the switch set
// is the search without it (the unit test `scope_is_the_two_hit_discontiguous_stage` checks
// that the stage refuses a contiguous table with a window).
#[test]
fn contiguous_tables_are_not_taken() {
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-task", "megablast", "-window_size", "40", "-outfmt", "6"])
        .output()
        .expect("run LOSAT");
    assert!(!output.status.success());
    let task_args = ["-task", "megablast"];
    let (expected, _) = blastn(&task_args, 4, &[]);
    let (report, stderr) = blastn(&task_args, 4, &[SWITCH]);
    assert!(report == expected, "megablast with {SWITCH:?}");
    assert_eq!(stderr, "");
    let (report, stderr) = blastn(&task_args, 4, &[SHADOW]);
    assert!(report == expected, "megablast with {SHADOW:?}");
    assert_eq!(
        stderr.trim_end(),
        "[X_PAIRPARDISCSHADOW] blastn: seed chunks equal=0 (ungapped hits 0, hash cells 0)"
    );
}
