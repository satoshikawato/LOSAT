// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! blastp's preliminary hit lists with `LOSAT_X_BPPFREE` (`src/algorithm/blastp/x_prelim_free.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. With it on, a reaped
//! preliminary HSP list gives back its spare room and the serial redo (one thread) frees each
//! query's preliminary list after its redo. Every test states the default it runs with: the
//! reference is `LOSAT_X_DEFAULT=0` (every switch off) and the searches with the switch on are the
//! default (with `LOSAT_X_BPQBLOCK=0` where one line per search is counted) and
//! `LOSAT_X_DEFAULT=0` with `LOSAT_X_BPPFREE=1`, the run the variable gave before it joined the
//! set; `LOSAT_X_BPPFREE=0` gives no line. These tests compare the reports, alone and with the
//! block search of `LOSAT_X_BPQBLOCK`, and check from the `[TIMING] x_bppfree: freed_lists=<n>`
//! lines that the serial redo freed lists, not only that the report is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3734-3738
//! ```c
//! match_loop_cleanup:
//!                 if (seqSrc) {
//!                     localMatch = Blast_HSPListFree(localMatch);
//!                 } else {
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One blastp run of the fixture pair (two related nimavirus proteomes) with the given variables
/// as its only `LOSAT_*` variables that choose a code path, and `extra` arguments.
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the tests is named.
/// `BATCH_SIZE` is the test's own. The run must exit with 0. Returns the standard output and the
/// numbers of lists that the `[TIMING] x_bppfree` lines report, one per search of a query set.
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPPFREE", "0");
const BLOCKS_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");

fn blastp(threads: u32, vars: &[(&str, &str)], extra: &[&str]) -> (Vec<u8>, Vec<usize>) {
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
            std::env::temp_dir().join("losat-bppfree-worker.log"),
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
    let freed = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bppfree: freed_lists="))
        .map(|count| count.parse().expect("a count of lists"))
        .collect();
    (output.stdout, freed)
}

#[test]
fn freed_lists_give_the_reports_of_the_kept_lists() {
    let batch = ("BATCH_SIZE", "3000");
    for outfmt in ["6", "0"] {
        for threads in [1, 4] {
            let (reference, reference_freed) =
                blastp(threads, &[batch, DEFAULT_OFF], &["-outfmt", outfmt]);
            assert!(!reference.is_empty(), "the fixture pair has a report");
            assert!(reference_freed.is_empty(), "no line with LOSAT_X_DEFAULT=0");
            let (report, freed) = blastp(
                threads,
                &[batch, BLOCKS_OFF, SWITCH_OFF],
                &["-outfmt", outfmt],
            );
            assert!(report == reference, "the switch off: the report differs");
            assert!(freed.is_empty(), "no line with the switch off: {freed:?}");
            for vars in [
                &[batch, BLOCKS_OFF][..],
                &[batch, DEFAULT_OFF, ("LOSAT_X_BPPFREE", "1")][..],
            ] {
                let (report, freed) = blastp(threads, vars, &["-outfmt", outfmt]);
                assert!(
                    report == reference,
                    "{vars:?} outfmt {outfmt} threads {threads}: the report differs"
                );
                assert_eq!(
                    freed.len(),
                    1,
                    "{vars:?} outfmt {outfmt} threads {threads}: {freed:?}"
                );
                if threads == 1 {
                    assert!(
                        freed[0] > 0,
                        "{vars:?} outfmt {outfmt}: the serial redo freed no list"
                    );
                }
            }
        }
    }
}

#[test]
fn freed_lists_in_query_blocks_give_the_rows_of_the_search_of_all_queries() {
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
                ("LOSAT_X_BPPFREE", "1"),
            ][..],
        ] {
            let (rows, freed) = blastp(threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "threads {threads} {vars:?}: rows differ");
            assert!(
                freed.len() > 1,
                "threads {threads} {vars:?}: one line per block, {freed:?}"
            );
            if threads == 1 {
                assert!(
                    freed.iter().sum::<usize>() > 0,
                    "{vars:?}: the blocks' serial redo freed no list: {freed:?}"
                );
            }
        }
    }
}
