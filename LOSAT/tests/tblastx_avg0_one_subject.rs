//! TBLASTX on a query batch whose average query length is 0, against one subject: NCBI BLAST+
//! 2.17.0 crashes there (SIGSEGV in its linking, exit status 139, no output) and LOSAT writes the
//! linked result, approved exception 20 of PD-LOSAT-NCBI-DEFECTS (version 1.68; the Owner's
//! decision of 2026-11-21, 「O2b-Q1=A」, line NP item O2b).
//!
//! NCBI divides the letters of a batch by all six contexts of every query (with `-strand plus` or
//! `minus` the contexts of the other strand count with length 0), so a 1,200-nt query with 300
//! queries of 3 nt on one strand, or with 1,300 on both, averages 0. Its expected length is then 0,
//! and with one subject the database is as long as the subject, so the big-gap cutoff of the
//! linking is `INT_MIN + 1`, and its `Int4` sums `score - cutoff[1]` overflow. LOSAT computes the
//! same cutoff and adds those sums in 64 bits, so each pass takes the chain of the most HSPs and,
//! among those, of the highest score. NCBI has no output to compare with: the expected files in
//! `tests/fixtures/tblastx_avg0_one_subject` are LOSAT's (dev ff07c6a7d, 2026-11-21), which the
//! default linking kernel and the literal port (`LOSAT_LINKING_LEGACY`) both give, at 1 and 4
//! threads. The nearest runs that NCBI completes (a second subject of 12 N and `-dbsize 3000`, a
//! big-gap cutoff of -4) are the cases `*avg0_pad_subject` of `tblastx_regression_fixtures.py`,
//! frozen from NCBI; they link otherwise (3 of 15, 3 of 20 and 6 of 35 of NCBI's rows are rows
//! here).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1023-1076
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:753,868

use std::path::{Path, PathBuf};
use std::process::Command;

/// (expected file, query of `tests/fixtures/tblastx_regression/inputs`, arguments after it)
const CASES: [(&str, &str, &[&str]); 3] = [
    ("strand_plus.out", "avg0_query.fa", &["-strand", "plus"]),
    ("strand_minus.out", "avg0_query.fa", &["-strand", "minus"]),
    ("both_strands.out", "avg0_both_query.fa", &[]),
];

fn engine(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

/// `losat tblastx` of `query` against the one subject Sa0mut with `-outfmt 6`, with no `LOSAT_*`,
/// `NCBI*`, `BATCH_SIZE` or `BLASTDB` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with
/// the `blastx-worker-probe` feature needs that file name once a search has more than one thread;
/// it selects nothing) and, for the literal port, `LOSAT_LINKING_LEGACY`.
fn tblastx(query: &str, args: &[&str], threads: &str, literal_port: bool) -> Vec<u8> {
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
        .arg(inputs.join("avg0_one_subject.fa"))
        .args(args)
        .args(["-outfmt", "6", "-num_threads", threads])
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&worker_log);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success() && stderr.is_empty(),
        "{query} {args:?} threads {threads}: {stderr}"
    );
    output.stdout
}

#[test]
fn a_batch_averaging_0_links_against_one_subject_as_both_kernels_do() {
    for (expected, query, args) in CASES {
        let expected_bytes =
            std::fs::read(engine("tests/fixtures/tblastx_avg0_one_subject").join(expected))
                .expect("expected output");
        assert!(!expected_bytes.is_empty(), "{expected}");
        for threads in ["1", "4"] {
            for literal_port in [false, true] {
                assert!(
                    tblastx(query, args, threads, literal_port) == expected_bytes,
                    "{expected}: threads {threads}, literal port {literal_port}"
                );
            }
        }
    }
}
