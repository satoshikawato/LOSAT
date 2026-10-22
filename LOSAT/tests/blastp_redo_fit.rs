//! blastp's redone HSP lists with `LOSAT_X_BPREDOFIT` (`src/algorithm/blastp/x_redo_fit.rs`) on
//! whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. With it on, the
//! composition-based redo builds each redone HSP list with room for exactly its alignments. Every
//! test states the default it runs with: the reference is `LOSAT_X_DEFAULT=0` (every switch off)
//! and the searches with the switch on are the default (with `LOSAT_X_BPQBLOCK=0` where one line
//! per search is counted) and `LOSAT_X_DEFAULT=0` with `LOSAT_X_BPREDOFIT=1`, the run the
//! variable gave before it joined the set; `LOSAT_X_BPREDOFIT=0` gives no line. These tests
//! compare the reports, alone and with the block search of `LOSAT_X_BPQBLOCK`, at one and four
//! threads, and check from the
//! `[TIMING] x_bpredofit: lists=<n> hsps=<m>` lines that lists were built that way, not only that
//! the report is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-356
//! ```c
//! for (align = *alignments; NULL != align; align = align->next) {
//!     GapEditScript * editScript = align->context;
//!     align->context = NULL;
//!     status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
//! }
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPREDOFIT", "0");
const BLOCKS_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");

/// One blastp run of the fixture pair (two related nimavirus proteomes) with the given variables
/// as its only `LOSAT_*` variables that choose a code path, and `extra` arguments.
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the tests is named.
/// `BATCH_SIZE` is the test's own. The run must exit with 0. Returns the standard output and the
/// (lists, HSPs) that the `[TIMING] x_bpredofit` lines report, one per search of a query set.
fn blastp(threads: u32, vars: &[(&str, &str)], extra: &[&str]) -> (Vec<u8>, Vec<(usize, usize)>) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpredofit-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
        .arg("-subject")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .env("LOSAT_TIMING", "1")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    let fitted = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bpredofit: lists="))
        .map(|counts| {
            let (lists, hsps) = counts.split_once(" hsps=").expect("lists and HSPs");
            (
                lists.parse().expect("a count of lists"),
                hsps.parse().expect("a count of HSPs"),
            )
        })
        .collect();
    (output.stdout, fitted)
}

#[test]
fn fitted_redone_lists_give_the_reports_of_pushed_lists() {
    let batch = ("BATCH_SIZE", "3000");
    for outfmt in ["6", "0"] {
        for threads in [1, 4] {
            let (reference, reference_fitted) =
                blastp(threads, &[batch, DEFAULT_OFF], &["-outfmt", outfmt]);
            assert!(!reference.is_empty(), "the fixture pair has a report");
            assert!(
                reference_fitted.is_empty(),
                "no line with LOSAT_X_DEFAULT=0"
            );
            let (report, fitted) = blastp(
                threads,
                &[batch, BLOCKS_OFF, SWITCH_OFF],
                &["-outfmt", outfmt],
            );
            assert!(report == reference, "the switch off: the report differs");
            assert!(fitted.is_empty(), "no line with the switch off: {fitted:?}");
            for vars in [
                &[batch, BLOCKS_OFF][..],
                &[batch, DEFAULT_OFF, ("LOSAT_X_BPREDOFIT", "1")][..],
            ] {
                let (report, fitted) = blastp(threads, vars, &["-outfmt", outfmt]);
                assert!(
                    report == reference,
                    "{vars:?} outfmt {outfmt} threads {threads}: the report differs"
                );
                assert_eq!(
                    fitted.len(),
                    1,
                    "{vars:?} outfmt {outfmt} threads {threads}: {fitted:?}"
                );
                let (lists, hsps) = fitted[0];
                assert!(
                    lists > 0 && hsps >= lists,
                    "{vars:?} outfmt {outfmt} threads {threads}: {fitted:?}"
                );
            }
        }
    }
}

#[test]
fn fitted_redone_lists_in_query_blocks_give_the_rows_of_the_search_of_all_queries() {
    let batch = ("BATCH_SIZE", "3000");
    for threads in [1, 4] {
        let (reference, _) = blastp(threads, &[batch, DEFAULT_OFF], &["-outfmt", "6"]);
        // The default blocks (50,000 letters; the fixture's queries have 77,448), one batch per
        // block, and the run the two variables gave before they joined the set.
        for vars in [
            &[batch][..],
            &[batch, ("LOSAT_X_BPQBLOCK", "1")][..],
            &[
                batch,
                DEFAULT_OFF,
                ("LOSAT_X_BPQBLOCK", "1"),
                ("LOSAT_X_BPREDOFIT", "1"),
            ][..],
        ] {
            let (rows, fitted) = blastp(threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "threads {threads} {vars:?}: rows differ");
            assert!(
                fitted.len() > 1,
                "threads {threads} {vars:?}: one line per block, {fitted:?}"
            );
            assert!(
                fitted.iter().map(|(lists, _)| lists).sum::<usize>() > 0,
                "threads {threads} {vars:?}: no list fitted: {fitted:?}"
            );
        }
    }
}
