//! The two-hit seed stage of dc-megablast for the diagonal array on the search pool
//! (`LOSAT_X_PAIRPARDISCARR`, `src/algorithm/blastn/blast_engine/x_pair_seed_disc.rs`) on whole
//! searches.
//!
//! NCBI keeps the diagonals of a query block of at most 8000 letters (both strands) in an array
//! (`kQueryLenForHashTable`), so a query of at most about 4000 bases takes the array, and
//! `LOSAT_X_PAIRPARDISC` (the diagonal hash only) never takes its chunks. The switch (in the
//! strict set: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set) chooses between
//! two LOSAT code paths with the same result: the subject scan cut into strips that the pool
//! threads scan, and the diagonal-array fold done by the owner thread in the scan order, or the
//! scan and fold of the reference code. So the property under test is that the report does not
//! depend on it: for discontiguous templates of both word sizes, one and two templates, the
//! search with the switch on at 1, 2, 4 and 5 threads gives the bytes of the search with it at
//! `0`, and the shadow run (`LOSAT_X_PAIRPARDISCARRSHADOW`, which computes the stage on copies
//! and compares the hit list and both arrays with the reference code for every chunk) gives
//! them too and compares at least one chunk. That those bytes equal NCBI's is the business of
//! the dc-megablast fixtures (`tests/blastn_regression_fixtures.py`, `dc.*`).
//!
//! The inputs are one subject (the stage takes a search of one subject by the owner of a pool
//! of two threads or more) and queries cut from LC738874 into files of this test: one query of
//! 3900 bases (7801 letters with both strands and the sentinel), and two queries of 1900 bases
//! (four contexts in one array).

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_PAIRPARDISCARR", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARDISCARRSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A query file of this test, removed when dropped.
struct QueryFile(PathBuf);

impl Drop for QueryFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// A query file of this test (one per test and process: the tests of this file run in
/// parallel) in Cargo's temporary directory for integration tests, with the pieces
/// `(start, length)` of LC738874 as records `q<start>`.
fn query_file(test: &str, pieces: &[(usize, usize)]) -> QueryFile {
    let text = std::fs::read_to_string(fasta("LC738874.fasta")).expect("read LC738874");
    let seq: String = text
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect();
    let mut out = String::new();
    for &(start, len) in pieces {
        out.push_str(&format!(">q{start}\n"));
        for line in seq.as_bytes()[start..start + len].chunks(70) {
            out.push_str(std::str::from_utf8(line).expect("ascii"));
            out.push('\n');
        }
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "blastn_pairpar_disc_arr_{}_{test}.fa",
        std::process::id()
    ));
    std::fs::write(&path, out).expect("write the query file");
    QueryFile(path)
}

/// One search of `query` against LC738868, `-outfmt 6` on standard output, with every `LOSAT_*`
/// variable dropped but the given ones (the strict set at its defaults) and
/// `LOSAT_BLASTX_WORKER_LOG`, which a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs as soon as a search has more than one thread
/// and which selects nothing. Returns the report and the standard error.
fn blastn(
    query: &Path,
    task_args: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-task", "dc-megablast"])
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

/// The `chunks=` and `seeds=` counts of the stage's `LOSAT_TIMING` line
/// (`[TIMING] x_pairpardiscarr_seed: chunks=N strips=M seeds=S ...`) in `stderr`.
fn stage_counts(stderr: &str) -> (u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] x_pairpardiscarr_seed: chunks="))
        .unwrap_or_else(|| panic!("no x_pairpardiscarr_seed line:\n{stderr}"));
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    (field("chunks="), field("seeds="))
}

/// For one template: the search with the switch at `0` at one thread, then by default at 1, 2,
/// 4 and 5 threads (the same bytes, nothing on standard error), then the shadow at 4 threads
/// (the same bytes, and its summary line with at least one chunk and one hit compared).
fn check_template(query: &Path, template: &[&str]) {
    let (expected, expected_stderr) = blastn(query, template, 1, &[SWITCH_OFF]);
    assert!(!expected.is_empty(), "{template:?}: no hits");
    assert_eq!(expected_stderr, "", "{template:?}");
    for threads in [1, 2, 4, 5] {
        let (report, stderr) = blastn(query, template, threads, &[]);
        assert!(
            report == expected,
            "{template:?} threads {threads}: the report differs from {SWITCH_OFF:?}"
        );
        assert_eq!(stderr, "", "{template:?} threads {threads}");
    }
    let (report, stderr) = blastn(query, template, 4, &[SHADOW]);
    assert!(
        report == expected,
        "{template:?}: the report differs with {SHADOW:?}"
    );
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_PAIRPARDISCARRSHADOW] blastn: seed chunks equal="))
        .unwrap_or_else(|| panic!("{template:?}: no summary line:\n{stderr}"));
    assert!(
        !line.starts_with("[X_PAIRPARDISCARRSHADOW] blastn: seed chunks equal=0 ")
            && !line.contains("(ungapped hits 0,"),
        "{template:?}: {line}"
    );
    assert_eq!(stderr.lines().count(), 1, "{template:?}: {stderr}");
}

// The byte comparisons above would also pass if the stage silently stopped running: with
// LOSAT_TIMING the stage prints how many chunks it took; with the switch at `0` or with
// `LOSAT_X_DEFAULT=0` it takes none and prints no line.
#[test]
fn the_stage_runs_by_default_and_not_when_switched_off() {
    let query = query_file("runs", &[(180_000, 3900)]);
    let (expected, _) = blastn(&query.0, &[], 1, &[SWITCH_OFF]);
    assert!(!expected.is_empty());
    let (report, stderr) = blastn(&query.0, &[], 4, &[("LOSAT_TIMING", "1")]);
    assert!(report == expected, "the report differs by default");
    let (chunks, seeds) = stage_counts(&stderr);
    assert!(chunks > 0 && seeds > 0, "the stage took nothing:\n{stderr}");
    for off in [
        &[SWITCH_OFF, ("LOSAT_TIMING", "1")][..],
        &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_TIMING", "1")][..],
    ] {
        let (report, stderr) = blastn(&query.0, &[], 4, off);
        assert!(report == expected, "the report differs with {off:?}");
        assert!(
            !stderr.contains("[TIMING] x_pairpardiscarr_seed:"),
            "{off:?}:\n{stderr}"
        );
    }
}

// The default template of dc-megablast (word size 11, coding, length 18), word size 12, the
// optimal template, and two templates (coding_and_optimal: per subject offset the first
// template's chain, then the second's).
#[test]
fn templates_give_the_same_report() {
    let query = query_file("templates", &[(180_000, 3900)]);
    check_template(&query.0, &[]);
    check_template(
        &query.0,
        &[
            "-word_size",
            "12",
            "-template_type",
            "coding",
            "-template_length",
            "21",
        ],
    );
    check_template(
        &query.0,
        &[
            "-word_size",
            "11",
            "-template_type",
            "optimal",
            "-template_length",
            "16",
        ],
    );
    check_template(
        &query.0,
        &[
            "-word_size",
            "11",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "18",
        ],
    );
}

// Two queries in one array: four contexts, each with its own cutoffs and masks, folded into one
// `hit_level_array` whose diagonals the contexts share.
#[test]
fn two_queries_in_one_array_give_the_same_report() {
    let query = query_file("two", &[(165_000, 1900), (45_000, 1900)]);
    check_template(&query.0, &[]);
}

// A query block longer than 8000 letters keeps the diagonal hash: the stage takes no chunk (the
// hash stage, LOSAT_X_PAIRPARDISC, takes it).
#[test]
fn the_diagonal_hash_is_not_taken() {
    let query = query_file("hash", &[(180_000, 4100)]);
    let (expected, _) = blastn(&query.0, &[], 1, &[SWITCH_OFF]);
    let (report, stderr) = blastn(&query.0, &[], 4, &[SHADOW, ("LOSAT_TIMING", "1")]);
    assert!(report == expected);
    assert!(
        stderr.contains(
            "[X_PAIRPARDISCARRSHADOW] blastn: seed chunks equal=0 (ungapped hits 0, array cells 0)"
        ),
        "{stderr}"
    );
    let (report, stderr) = blastn(&query.0, &[], 4, &[("LOSAT_TIMING", "1")]);
    assert!(report == expected);
    assert_eq!(stage_counts(&stderr).0, 0, "{stderr}");
}
