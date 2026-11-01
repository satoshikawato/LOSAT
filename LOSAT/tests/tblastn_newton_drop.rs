//! The Newton results computed ahead and not taken, dropped at the next prefetch
//! (`LOSAT_X_NEWTONDROP`, in the strict set; `src/core/composition_adjustment/x_newton_lanes.rs`
//! `prefetch`), on whole TBLASTN searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1143-1143
//! ```c
//!     for (window_index = 0;  window_index < nWindows;  window_index++) {
//! ```
//! TBLASTN's redo walks a match's windows once, in order, and solves the Newton problems of the
//! coming windows ahead (Newton lanes); a result is taken only by the call with the same input
//! bits, and a call whose result is missing solves its problem itself with the same bits. The
//! results of earlier windows not taken when the next windows are solved ahead (HSPs found
//! contained in the redo) are never asked for again: with the switch they are dropped then, with
//! `=0` they stay in the thread's store until it holds 512. These tests compare the reports and
//! read the `[X_STATS] NEWTONLANES` and `[X_STATS] NEWTONDROP` lines (`LOSAT_X_STATS=1`): the same
//! results taken in both, and only the default drops any.
//!
//! The queries are the 115 proteins of `tests/fasta/PemoMJNVA.faa`; the subject is
//! `tests/fasta/MelaMJNV.fasta` (63 of the 2,999 results computed ahead are not taken).

use std::path::Path;
use std::process::Command;

/// One TBLASTN run (`-outfmt 6` on standard output) at `threads` with `variables` and
/// `LOSAT_X_STATS` as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build with
/// the worker probe (`--all-features`) writes to that file as soon as a pool starts, so the
/// caller's value is kept, or a file of the tests is named. Returns the report and the standard
/// error.
fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_newton_drop");
        std::fs::create_dir_all(&dir).expect("work directory");
        command.env(
            WORKER_LOG,
            dir.join(format!("worker.{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(fasta.join("PemoMJNVA.faa"))
        .arg("-subject")
        .arg(fasta.join("MelaMJNV.fasta"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The values of `fields` on the one line of `stderr` that starts with `head`.
fn counts<const K: usize>(stderr: &str, head: &str, fields: [&str; K]) -> [u64; K] {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with(head))
        .collect();
    assert_eq!(lines.len(), 1, "one {head} line per process:\n{stderr}");
    fields.map(|field| {
        lines[0]
            .split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {}", lines[0]))
            .parse()
            .expect("a number")
    })
}

/// The default's report is `=0`'s at one and four threads, with the same results computed ahead and
/// taken; the default drops some of the results not taken, `=0` none.
#[test]
fn dropping_results_not_taken_keeps_the_report() {
    for threads in [1, 4] {
        let (drop, drop_err) = tblastn(threads, &[]);
        let (keep, keep_err) = tblastn(threads, &[("LOSAT_X_NEWTONDROP", "0")]);
        assert!(!drop.is_empty(), "the search finds hits");
        assert_eq!(drop, keep, "{threads} threads: the reports differ");
        let lanes = ["prefetched", "hits", "misses", "unused"];
        let [prefetched, hits, misses, unused] = counts(&drop_err, "[X_STATS] NEWTONLANES ", lanes);
        assert_eq!(
            counts(&keep_err, "[X_STATS] NEWTONLANES ", lanes),
            [prefetched, hits, misses, unused],
            "{threads} threads: the results computed ahead and taken differ"
        );
        assert!(unused > 0, "{threads} threads: some results are not taken");
        let [dropped] = counts(&drop_err, "[X_STATS] NEWTONDROP ", ["dropped"]);
        assert!(
            dropped > 0 && dropped <= unused,
            "{threads} threads: {dropped} dropped of {unused}"
        );
        assert_eq!(
            counts(&keep_err, "[X_STATS] NEWTONDROP ", ["dropped"]),
            [0],
            "=0 drops none"
        );
    }
}
