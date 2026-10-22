//! Which TBLASTN query batches are searched side by side (`LOSAT_X_TBNBATCH` and
//! `LOSAT_X_TBNBATCHRULE`) and which at the same time on pools of their own
//! (`LOSAT_X_TBNOWNPOOLS`), all in the strict set; `src/algorithm/tblastn/stage_d_pipeline.rs`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-289
//! ```c++
//!             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!
//!                 query =  input->GetNextSeqBatch(*scope);
//!                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
//!                 ...
//!                     CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
//! ```
//! NCBI searches each query batch on its own, one after another. With `LOSAT_X_TBNBATCH`, two or
//! more threads, the subject plans shared (four or more batches) and at least one batch per
//! thread, LOSAT searches the batches side by side, each as the one-thread search. With
//! `LOSAT_X_TBNBATCHRULE` as well (the default), only the batches of the full waves (a multiple of
//! the threads) go side by side, and only for a small subject set or fewer than two subjects per
//! thread; the rest are searched one after another with the pool inside each. With
//! `LOSAT_X_TBNBATCHRULE=0` every batch goes side by side once there is one per thread. With
//! `LOSAT_X_TBNOWNPOOLS` as well (the default), the batches the rule leaves one after another go
//! at the same time, each on a pool of its own; `LOSAT_X_TBNOWNPOOLS=0` keeps them one after
//! another. None of these may change the report, so the tests compare it with the one-thread search and read the
//! `[losat-thread-stage] ... stage=query_batches` line (`LOSAT_WASI_THREADS_DEBUG=1`) to check how
//! many batches went side by side.
//!
//! `tests/fasta/PemoMJNVA.faa` (115 proteins, 72,314 residues) is four query batches as NCBI's
//! reader cuts them (20,000 residues each); the subject is `tests/fasta/MelaMJNV.fasta`.

use std::path::Path;
use std::process::Command;

/// The query batches of `PemoMJNVA.faa`.
const BATCHES: usize = 4;

/// One TBLASTN run of `PemoMJNVA.faa` against `MelaMJNV.fasta` (`-outfmt 6` on standard output)
/// with the given variables as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a
/// build with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so
/// the caller's value is kept, or a file of the tests is named. Returns the report and the
/// standard error.
fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    let query = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
    tblastn_query(&query, threads, variables)
}

/// `tblastn` with another query file.
fn tblastn_query(query: &Path, threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_batch_rule");
        std::fs::create_dir_all(&dir).expect("work directory");
        command.env(WORKER_LOG, dir.join("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta.join("MelaMJNV.fasta"))
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .env("LOSAT_WASI_THREADS_DEBUG", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The `side_by_side` field of the one `query_batches` line, or `None` without the line (the
/// batches left to the loop of one batch after another before the switch is read).
fn side_by_side(stderr: &str) -> Option<usize> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| {
            line.starts_with("[losat-thread-stage] program=tblastn stage=query_batches ")
        })
        .collect();
    assert!(lines.len() <= 1, "one line per search:\n{stderr}");
    let line = lines.first()?;
    assert!(
        line.contains(&format!(" batches={BATCHES} ")),
        "four query batches: {line}"
    );
    Some(
        line.split_whitespace()
            .find_map(|field| field.strip_prefix("side_by_side="))
            .unwrap_or_else(|| panic!("no side_by_side in {line}"))
            .parse()
            .expect("number of batches"),
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-289 (above). Every batch is
// searched as NCBI searches it, side by side or not, and the reports keep the batch order, so
// each arm must give the one-thread report.
#[test]
fn the_report_is_the_one_thread_report_with_and_without_the_rule() {
    let (expected, _) = tblastn(1, &[]);
    assert!(!expected.is_empty(), "the search must report hits");
    for threads in [2, 3, 4, 8] {
        // The default has the rule: the full waves only (one small subject).
        let (report, stderr) = tblastn(threads, &[]);
        assert_eq!(report, expected, "default, {threads} threads");
        assert_eq!(
            side_by_side(&stderr),
            Some(BATCHES / threads * threads),
            "default, {threads} threads"
        );
        // Without the rule: every batch side by side once there is one per thread.
        let (report, stderr) = tblastn(threads, &[("LOSAT_X_TBNBATCHRULE", "0")]);
        assert_eq!(report, expected, "rule off, {threads} threads");
        let all = if BATCHES >= threads { BATCHES } else { 0 };
        assert_eq!(
            side_by_side(&stderr),
            Some(all),
            "rule off, {threads} threads"
        );
    }
}

// The rule acts only with `LOSAT_X_TBNBATCH`: with it off, alone or through `LOSAT_X_DEFAULT=0`,
// the batches are searched one after another and the report is unchanged.
#[test]
fn the_rule_does_nothing_without_the_side_by_side_switch() {
    let (expected, _) = tblastn(1, &[]);
    for variables in [
        &[("LOSAT_X_TBNBATCH", "0")][..],
        &[("LOSAT_X_DEFAULT", "0")],
    ] {
        let (report, stderr) = tblastn(3, variables);
        assert_eq!(report, expected, "{variables:?}");
        assert_eq!(side_by_side(&stderr), None, "{variables:?}");
    }
}

/// The `batches` and `own_pools` fields of the one `query_batches` line, and the number of
/// `[losat-thread-pool]` lines (one per search pool started).
fn own_pools(stderr: &str) -> (usize, usize, usize) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[losat-thread-stage] program=tblastn stage=query_batches "))
        .unwrap_or_else(|| panic!("no query_batches line:\n{stderr}"));
    let field = |name: &str| -> usize {
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("number")
    };
    let pools = stderr
        .lines()
        .filter(|line| line.starts_with("[losat-thread-pool] program=tblastn "))
        .count();
    (field("batches="), field("own_pools="), pools)
}

// EXPERIMENT (LOSAT_X_TBNOWNPOOLS, in the strict set): the batches the rule leaves one after another (two or
// more, fewer than the threads, one small subject here) are searched at the same time, each on a
// pool of its own with its share of the threads. NCBI reference (598d8ae6):
// c++/src/app/blast/tblastn_app.cpp:275-289 (above): each batch is still searched as NCBI
// searches it and the reports keep the batch order, so every arm must give the one-thread
// report.
#[test]
fn own_pools_keep_the_one_thread_report() {
    let (expected, _) = tblastn(1, &[]);
    for threads in [2, 3, 5, 8, 16] {
        let (report, stderr) = tblastn(threads, &[]);
        assert_eq!(report, expected, "own pools, {threads} threads");
        let full = BATCHES / threads * threads;
        let left = BATCHES - full;
        let own = if (2..threads).contains(&left) {
            left
        } else {
            0
        };
        assert_eq!(
            own_pools(&stderr),
            (BATCHES, own, 1 + own),
            "own pools, {threads} threads"
        );
        // With `LOSAT_X_TBNOWNPOOLS=0` no batch goes to a pool of its own.
        let (report, stderr) = tblastn(threads, &[("LOSAT_X_TBNOWNPOOLS", "0")]);
        assert_eq!(report, expected, "own pools off, {threads} threads");
        assert_eq!(
            own_pools(&stderr),
            (BATCHES, 0, 1),
            "own pools off, {threads} threads"
        );
    }
    // The rule decides which batches are left: without it, no pools of their own.
    let (report, stderr) = tblastn(8, &[("LOSAT_X_TBNBATCHRULE", "0")]);
    assert_eq!(report, expected, "own pools without the rule");
    assert_eq!(
        own_pools(&stderr),
        (BATCHES, 0, 1),
        "own pools without the rule"
    );
}

// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): two query batches (fewer than four, so the subject plans are
// not shared and no batch goes side by side) on pools of their own at four threads.
#[test]
fn own_pools_without_shared_subject_plans_keep_the_one_thread_report() {
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let text = std::fs::read_to_string(fasta.join("PemoMJNVA.faa")).expect("query");
    // The records up to 25,000 residues: two batches of NCBI's reader (20,000 residues each).
    let mut kept = String::new();
    let mut residues = 0;
    for record in text.split_inclusive('>').filter(|r| *r != ">") {
        let record = format!(">{}", record.trim_end_matches('>'));
        residues += record
            .lines()
            .skip(1)
            .map(|line| line.trim().len())
            .sum::<usize>();
        kept.push_str(&record);
        if residues > 25_000 {
            break;
        }
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_batch_rule");
    std::fs::create_dir_all(&dir).expect("work directory");
    let query = dir.join("two_batches.faa");
    std::fs::write(&query, kept).expect("query file");
    let (expected, _) = tblastn_query(&query, 1, &[]);
    assert!(!expected.is_empty(), "the search must report hits");
    let (report, stderr) = tblastn_query(&query, 4, &[]);
    assert_eq!(report, expected, "own pools, two batches");
    assert_eq!(own_pools(&stderr), (2, 2, 3), "own pools, two batches");
}
