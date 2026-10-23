//! blastp's preliminary hit lists with `LOSAT_X_BPPCOMPACT` (`src/algorithm/blastp/x_prelim_compact.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`. With it on, the subject
//! scan keeps each preliminary HSP as a 64-byte element and the redo reads the hit lists rebuilt
//! from them. Every test states the default it runs with: the reference is the search of all
//! queries with the switch off (`LOSAT_X_BPQBLOCK=0 LOSAT_X_BPPCOMPACT=0`), and the search with
//! the switch on is the default (with `LOSAT_X_BPQBLOCK=0` where one line per search is
//! counted); `LOSAT_X_BPPCOMPACT=0` gives no line. These tests compare the reports: with one and four
//! threads, with `LOSAT_X_BPPFREE`, in the blocks of `LOSAT_X_BPQBLOCK` and with query chunks
//! (`CHUNK_SIZE`); they check from the `[TIMING] x_bppcompact: expanded_lists=<n>` lines that
//! compact lists were expanded, and that the shadow mode (`LOSAT_X_BPPCOMPACT=shadow`) compared
//! lists (`[X_SHADOW] BPPCOMPACT compared=<n>`; it panics on a difference).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1523-1539
//! ```c
//!          Blast_HSPListGetEvalues(program_number, query_info,
//!                                  stat_length, hsp_list,
//!                                  gapped_calculation, FALSE,
//!                                  sbp, 0, 1.0);
//!          status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPPCOMPACT", "0");
const BLOCKS_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");

/// One blastp run of the fixture pair (two related nimavirus proteomes) with the given variables
/// as its only `LOSAT_*` and query-batch variables, and `extra` arguments.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the tests is named.
/// The run must exit with 0. Returns the standard output, the numbers of lists that the
/// `[TIMING] x_bppcompact` lines report (one per search of a query set) and the number of lists
/// on the `[X_SHADOW] BPPCOMPACT` line, if any.
fn blastp(
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, Vec<usize>, Option<usize>) {
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
            std::env::temp_dir().join("losat-bppcompact-worker.log"),
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
    let expanded = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bppcompact: expanded_lists="))
        .map(|count| count.parse().expect("a count of lists"))
        .collect();
    let compared = stderr.lines().find_map(|line| {
        line.strip_prefix("[X_SHADOW] BPPCOMPACT compared=")
            .and_then(|rest| rest.split_whitespace().next())
            .map(|count| count.parse().expect("a count of lists"))
    });
    (output.stdout, expanded, compared)
}

#[test]
fn compact_lists_give_the_reports_of_the_full_lists() {
    let batch = ("BATCH_SIZE", "3000");
    for outfmt in ["6", "0"] {
        for threads in [1, 4] {
            let (reference, reference_expanded, _) = blastp(
                threads,
                &[batch, BLOCKS_OFF, SWITCH_OFF],
                &["-outfmt", outfmt],
            );
            assert!(!reference.is_empty(), "the fixture pair has a report");
            assert!(
                reference_expanded.is_empty(),
                "no line with the switch off: {reference_expanded:?}"
            );
            for vars in [
                &[batch, BLOCKS_OFF][..],
                &[batch, BLOCKS_OFF, ("LOSAT_X_BPPFREE", "1")][..],
            ] {
                let (report, expanded, _) = blastp(threads, vars, &["-outfmt", outfmt]);
                assert!(
                    report == reference,
                    "{vars:?} outfmt {outfmt} threads {threads}: the report differs"
                );
                assert!(
                    expanded.len() == 1 && expanded[0] > 0,
                    "{vars:?} outfmt {outfmt} threads {threads}: {expanded:?}"
                );
            }
        }
    }
}

#[test]
fn compact_lists_in_query_blocks_and_chunks_give_the_same_rows() {
    let batch = ("BATCH_SIZE", "3000");
    let chunks = [("CHUNK_SIZE", "1000"), ("OVERLAP_CHUNK_SIZE", "100")];
    for threads in [1, 4] {
        let (reference, _, _) =
            blastp(threads, &[batch, BLOCKS_OFF, SWITCH_OFF], &["-outfmt", "6"]);
        // The default blocks (50,000 letters; the fixture's queries have 77,448), one batch per
        // block.
        for vars in [&[batch][..], &[batch, ("LOSAT_X_BPQBLOCK", "1")][..]] {
            let (rows, expanded, _) = blastp(threads, vars, &["-outfmt", "6"]);
            assert!(
                rows == reference,
                "threads {threads} {vars:?}: block rows differ"
            );
            assert!(
                expanded.len() > 1,
                "threads {threads} {vars:?}: one line per block, {expanded:?}"
            );
        }
        let (chunk_reference, _, _) = blastp(
            threads,
            &[batch, chunks[0], chunks[1], BLOCKS_OFF, SWITCH_OFF],
            &["-outfmt", "6"],
        );
        let vars = &[batch, chunks[0], chunks[1], BLOCKS_OFF][..];
        let (chunk_rows, chunk_expanded, _) = blastp(threads, vars, &["-outfmt", "6"]);
        assert!(
            chunk_rows == chunk_reference,
            "threads {threads} {vars:?}: rows with query chunks differ"
        );
        assert!(
            chunk_expanded.len() == 1 && chunk_expanded[0] > 0,
            "threads {threads} {vars:?}: {chunk_expanded:?}"
        );
    }
}

#[test]
fn shadow_mode_compares_every_list() {
    let batch = ("BATCH_SIZE", "3000");
    let chunks = [("CHUNK_SIZE", "1000"), ("OVERLAP_CHUNK_SIZE", "100")];
    let shadow = ("LOSAT_X_BPPCOMPACT", "shadow");
    for threads in [1, 4] {
        let (reference, _, _) =
            blastp(threads, &[batch, BLOCKS_OFF, SWITCH_OFF], &["-outfmt", "6"]);
        // The shadow runs with every other switch on (one search of all queries).
        let (rows, _, compared) = blastp(threads, &[batch, BLOCKS_OFF, shadow], &["-outfmt", "6"]);
        assert!(rows == reference, "threads {threads}: shadow rows differ");
        assert!(
            compared.is_some_and(|lists| lists > 0),
            "threads {threads}: {compared:?}"
        );
        let (_, _, chunk_compared) = blastp(
            threads,
            &[batch, chunks[0], chunks[1], BLOCKS_OFF, shadow],
            &["-outfmt", "6"],
        );
        assert!(
            chunk_compared.is_some_and(|lists| lists > compared.unwrap()),
            "threads {threads}: query chunks add their lists, {chunk_compared:?}"
        );
    }
}
