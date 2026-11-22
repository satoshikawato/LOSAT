//! TBLASTX on a query batch whose average query length is 0, against a subject that the database
//! adds no protein letter to (one subject, or other subjects of fewer than 3 bases in all): NCBI
//! BLAST+ 2.17.0 crashes there (SIGSEGV in its linking, exit status 139, no output) and LOSAT
//! writes the linked result, approved exception 20 of PD-LOSAT-NCBI-DEFECTS (version 1.68; the
//! Owner's decision of 2026-11-21, 「O2b-Q1=A」, line NP item O2b).
//!
//! NCBI's average query length of a batch is the end of the last context in the batch's
//! concatenated query (the letters and one separator after each context with letters) divided by
//! all six contexts of every query; with `-strand plus` or `minus` the contexts of the other strand
//! have no letters. So a 1,200-nt query with 300 queries of 3 nt on one strand, or with 1,300 on
//! both, averages 0. Its expected length is then 0, and when the database's length divided by 3 is
//! the subject's, the big-gap cutoff of the linking is one more than `INT_MIN`, and its `Int4` sums
//! with `cutoff[1]` overflow. LOSAT computes the same cutoff and adds those sums in 64 bits. With a
//! subject of up to 60,002 bases the linking keeps big gaps only, so each pass takes the chain of
//! the most HSPs and, among those, of the highest score; with a longer one (`long_subject_minus`)
//! each pass chooses between that chain and the small-gap chain by their e-values, as NCBI does.
//! NCBI has no output to compare with: the expected files in
//! `tests/fixtures/tblastx_avg0_one_subject` are LOSAT's (line NP item O2b, 2026-11-21), which the
//! default linking kernel and the literal port (`LOSAT_LINKING_LEGACY`) both give, at 1 and 4
//! threads. The nearest runs that NCBI completes (a second subject of 12 N and `-dbsize 3000`, a
//! big-gap cutoff of -4) are the cases `*avg0_pad_subject` of `tblastx_regression_fixtures.py`,
//! frozen from NCBI: each of their alignments is one of the first three files' rows, but they link
//! otherwise, so the sets and their e-values differ.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1023-1076
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:868,942-943

use std::path::{Path, PathBuf};
use std::process::Command;

/// (expected file, query and subject of `tests/fixtures/tblastx_regression/inputs`, arguments
/// after them)
const CASES: [(&str, &str, &str, &[&str]); 5] = [
    (
        "strand_plus.out",
        "avg0_query.fa",
        "avg0_one_subject.fa",
        &["-strand", "plus"],
    ),
    (
        "strand_minus.out",
        "avg0_query.fa",
        "avg0_one_subject.fa",
        &["-strand", "minus"],
    ),
    (
        "both_strands.out",
        "avg0_both_query.fa",
        "avg0_one_subject.fa",
        &[],
    ),
    // A second subject of 2 nt: no codon, so the database's length divided by 3 stays Sa0mut's.
    (
        "second_subject_2nt_plus.out",
        "avg0_query.fa",
        "avg0_two_subject.fa",
        &["-strand", "plus"],
    ),
    // Sa0mut followed by 57,003 N: a search space past 8 x 50^2, so the small-gap chains compete.
    (
        "long_subject_minus.out",
        "avg0_query.fa",
        "avg0_long_subject.fa",
        &["-strand", "minus"],
    ),
];

fn engine(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

/// `losat tblastx` of `query` against `subject` with `-outfmt 6`, with no `LOSAT_*`, `NCBI*`,
/// `BATCH_SIZE` or `BLASTDB` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with the
/// `blastx-worker-probe` feature needs that file name once a search has more than one thread; it
/// selects nothing) and, for the literal port, `LOSAT_LINKING_LEGACY`.
fn tblastx(
    query: &str,
    subject: &str,
    args: &[&str],
    threads: &str,
    literal_port: bool,
) -> Vec<u8> {
    let inputs = engine("tests/fixtures/tblastx_regression/inputs");
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(text.as_ref(), "BATCH_SIZE" | "BLASTDB")
        {
            command.env_remove(name);
        }
    }
    if literal_port {
        command.env("LOSAT_LINKING_LEGACY", "1");
    }
    let worker_log = std::env::temp_dir().join(format!(
        "losat-tblastx-avg0-one-subject-{}.log",
        std::process::id()
    ));
    let output = command
        .env("LOSAT_BLASTX_WORKER_LOG", &worker_log)
        .arg("tblastx")
        .arg("-query")
        .arg(inputs.join(query))
        .arg("-subject")
        .arg(inputs.join(subject))
        .args(args)
        .args(["-outfmt", "6", "-num_threads", threads])
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&worker_log);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success() && stderr.is_empty(),
        "{query} {subject} {args:?} threads {threads}: {stderr}"
    );
    output.stdout
}

/// The number of the first line where `actual` and `expected` differ, with both line counts.
fn first_difference(actual: &[u8], expected: &[u8]) -> String {
    let actual_lines: Vec<&[u8]> = actual.split(|&byte| byte == b'\n').collect();
    let expected_lines: Vec<&[u8]> = expected.split(|&byte| byte == b'\n').collect();
    let line = actual_lines
        .iter()
        .zip(&expected_lines)
        .position(|(a, b)| a != b)
        .unwrap_or(actual_lines.len().min(expected_lines.len()));
    format!(
        "{} lines, expected {}; first difference at line {}",
        actual_lines.len() - 1,
        expected_lines.len() - 1,
        line + 1
    )
}

#[test]
fn a_batch_averaging_0_links_against_one_subject_as_both_kernels_do() {
    for (expected, query, subject, args) in CASES {
        let expected_bytes =
            std::fs::read(engine("tests/fixtures/tblastx_avg0_one_subject").join(expected))
                .expect("expected output");
        assert!(!expected_bytes.is_empty(), "{expected}");
        for threads in ["1", "4"] {
            for literal_port in [false, true] {
                let actual = tblastx(query, subject, args, threads, literal_port);
                assert!(
                    actual == expected_bytes,
                    "{expected}: threads {threads}, literal port {literal_port}: {}",
                    first_difference(&actual, &expected_bytes)
                );
            }
        }
    }
}
