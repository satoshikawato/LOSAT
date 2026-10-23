//! blastp's preliminary hit lists with `LOSAT_X_BPPCOMPACT` (`src/algorithm/blastp/x_prelim_compact.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`. With it on, the subject
//! scan keeps each preliminary HSP as a 64-byte element and the redo reads the hit lists rebuilt
//! from them. Every test states the default it runs with: the reference is the search of all
//! queries with the switch off (`LOSAT_X_BPQBLOCK=0 LOSAT_X_BPPCOMPACT=0`), and the search with
//! the switch on is the default (with `LOSAT_X_BPQBLOCK=0` where one line per search is
//! counted); `LOSAT_X_BPPCOMPACT=0` gives no line. These tests compare the outfmt 6 reports
//! (every report is written from the redo's results, which read the expanded lists), one search
//! for each place where these searches expand the lists: the one-thread redo (each query's list
//! when it reaches it), the windows of the redo at four threads (`LOSAT_X_BPREDOWINEXPAND`, in
//! the default blocks of `LOSAT_X_BPQBLOCK`) and the end of a query chunk's search
//! (`CHUNK_SIZE`). They check from the `[TIMING] x_bppcompact: expanded_lists=<n>` lines that
//! compact lists were expanded, and that the shadow mode (`LOSAT_X_BPPCOMPACT=shadow`) compared
//! lists (`[X_SHADOW] BPPCOMPACT compared=<n>`; it panics on a difference).
//!
//! Other files run the compact lists (the default) beside other settings:
//! `blastp_prelim_free.rs` compares `LOSAT_X_BPPFREE=0` with the default at one and four
//! threads, in outfmt 6 and 0 and in blocks; `blastp_redo_par.rs` compares
//! `LOSAT_X_BPREDOWINEXPAND=0` (every list expanded before a redo of more than one thread) with
//! the windows; `blastp_query_blocks.rs` runs blocks of 1 and 20,000 letters at one and four
//! threads; `blastp_hsp_tie_order.rs` runs a four-thread search of one query against one subject
//! (one match, no windows), whose list is expanded before the redo.
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
    let (reference, reference_expanded, _) =
        blastp(1, &[batch, BLOCKS_OFF, SWITCH_OFF], &["-outfmt", "6"]);
    assert!(!reference.is_empty(), "the fixture pair has a report");
    assert!(
        reference_expanded.is_empty(),
        "no line with the switch off: {reference_expanded:?}"
    );
    // One thread: the redo expands each query's list when it reaches it, in one search of all
    // queries (one line).
    let (rows, expanded, _) = blastp(1, &[batch, BLOCKS_OFF], &["-outfmt", "6"]);
    assert!(rows == reference, "threads 1: the report differs");
    assert!(
        expanded.len() == 1 && expanded[0] > 0,
        "threads 1: {expanded:?}"
    );
    // Four threads: the windows of the redo expand the lists as they reach them, in the default
    // blocks (50,000 letters; the fixture's queries have 77,448), one line per block. The report
    // is the one-thread search's.
    let (rows, expanded, _) = blastp(4, &[batch], &["-outfmt", "6"]);
    assert!(
        rows == reference,
        "threads 4, default blocks: the report differs"
    );
    assert!(
        expanded.len() > 1 && expanded.iter().all(|&lists| lists > 0),
        "threads 4, default blocks: one line per block, {expanded:?}"
    );
}

#[test]
fn compact_lists_of_query_chunks_give_the_same_rows() {
    let batch = ("BATCH_SIZE", "3000");
    let chunks = [("CHUNK_SIZE", "1000"), ("OVERLAP_CHUNK_SIZE", "100")];
    // A query batch of 1,800 letters or more is split into chunks: here every batch of 3,000
    // letters but the last (6 queries, 1,343 letters). Each chunk's search expands its lists
    // before the chunk merge; at four threads the windows of the redo then read the merged lists
    // beside the compact lists of the last batch's queries.
    let (reference, _, _) = blastp(
        4,
        &[batch, chunks[0], chunks[1], BLOCKS_OFF, SWITCH_OFF],
        &["-outfmt", "6"],
    );
    let vars = &[batch, chunks[0], chunks[1], BLOCKS_OFF][..];
    let (rows, expanded, _) = blastp(4, vars, &["-outfmt", "6"]);
    assert!(rows == reference, "{vars:?}: rows with query chunks differ");
    assert!(
        expanded.len() == 1 && expanded[0] > 0,
        "{vars:?}: {expanded:?}"
    );
}

#[test]
fn shadow_mode_compares_every_list() {
    let batch = ("BATCH_SIZE", "3000");
    let chunks = [("CHUNK_SIZE", "1000"), ("OVERLAP_CHUNK_SIZE", "100")];
    let shadow = ("LOSAT_X_BPPCOMPACT", "shadow");
    // The shadow mode compares the lists before the redo at any thread count; four threads.
    let (reference, _, _) = blastp(4, &[batch, BLOCKS_OFF, SWITCH_OFF], &["-outfmt", "6"]);
    // The shadow runs with every other switch on (one search of all queries).
    let (rows, _, compared) = blastp(4, &[batch, BLOCKS_OFF, shadow], &["-outfmt", "6"]);
    assert!(rows == reference, "the shadow rows differ");
    assert!(compared.is_some_and(|lists| lists > 0), "{compared:?}");
    let (_, _, chunk_compared) = blastp(
        4,
        &[batch, chunks[0], chunks[1], BLOCKS_OFF, shadow],
        &["-outfmt", "6"],
    );
    assert!(
        chunk_compared.is_some_and(|lists| lists > compared.unwrap()),
        "query chunks add their lists, {chunk_compared:?}"
    );
}
