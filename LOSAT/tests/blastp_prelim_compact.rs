//! blastp's preliminary hit lists as compact elements (`src/algorithm/blastp/x_prelim_compact.rs`)
//! on whole searches.
//!
//! The subject scan keeps each preliminary HSP as a 64-byte element and the redo reads the hit
//! lists rebuilt from them. Every test states the default it runs with: the reference is the
//! one-thread search of all queries (`LOSAT_X_BPQBLOCK=0`, where one line per search is counted,
//! and blocks of 1,000 letters where one line per block is). These tests compare the outfmt 6
//! reports (every report is written from the redo's results, which read the expanded lists) with
//! the reference (that the reports equal NCBI's is the business of the frozen fixtures), one search
//! for each place where these searches expand the lists: the one-thread redo (each query's list
//! when it reaches it), the windows of the redo at four threads (`LOSAT_X_BPREDOWINEXPAND`, in
//! blocks of `LOSAT_X_BPQBLOCK`) and the end of a query chunk's search
//! (`CHUNK_SIZE`). They check from the `[TIMING] x_bppcompact: expanded_lists=<n>` lines that
//! compact lists were expanded, and from the `[TIMING] x_bpredopar:` lines that the redo ran in
//! windows.
//!
//! The searches run five of the fixture's 101 queries (`queries`), one of them the fixture's only
//! query without a preliminary list. The compact lists read no size that the five change: a list
//! is expanded where a search reaches it whatever its length, and no preliminary list is full
//! (the fixture has 124 subjects; the rule of a full list is the unit test
//! `hit_list_helpers_keep_the_same_lists_on_both_types` in `x_prelim_compact.rs`). A default
//! window of the redo (64 jobs per thread) holds all the jobs of a search here, so the search with
//! query chunks runs in windows of one job per thread (`SMALL_WINDOWS`), several of them;
//! the expansion ahead over several windows of compact lists runs on all 101 queries in
//! `blastp_redo_par.rs`, which compares `LOSAT_X_BPREDOWINEXPAND=0` (every list expanded before a
//! redo of more than one thread) with the windows, small and default. The default blocks of
//! 50,000 letters run in the default test of `blastp_query_blocks.rs`.
//!
//! Other files run the compact lists beside other settings:
//! `blastp_prelim_free.rs` compares the reports at one and four threads, in outfmt 6 and 0 and in
//! blocks; `blastp_query_blocks.rs` runs blocks of 1 and 1,000
//! letters at one and four threads; `blastp_hsp_tie_order.rs` runs a four-thread search of one
//! query against one subject (one match, no windows), whose list is expanded before the redo.
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

const BLOCKS_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");
/// The batches of `queries`: {BDT62581.1}, {BDT62582.1, BDT62583.1}, {BDT62584.1} and
/// {BDT62553.1} (a batch takes records while it holds fewer letters).
const BATCH: (&str, &str) = ("BATCH_SIZE", "500");
/// Blocks of several batches of `queries` (a block takes batches while it holds fewer letters):
/// {BDT62581.1 to BDT62583.1} (1,595 letters) and {BDT62584.1, BDT62553.1}.
const TWO_BLOCKS: (&str, &str) = ("LOSAT_X_BPQBLOCK", "1000");
/// A batch of `queries` is split into chunks when it holds at least twice `CHUNK_SIZE -
/// OVERLAP_CHUNK_SIZE` letters (700): {BDT62582.1, BDT62583.1} (897 letters) into two, the first
/// of which holds parts of both queries, and {BDT62584.1} (944) into two; {BDT62581.1} (698) and
/// {BDT62553.1} are not split.
const CHUNKS: [(&str, &str); 2] = [("CHUNK_SIZE", "450"), ("OVERLAP_CHUNK_SIZE", "100")];
/// Windows of the redo of one job per thread (`LOSAT_X_BPREDOPARWINDOW`, as in
/// `blastp_redo_par.rs`; the default holds 64 jobs per thread).
const SMALL_WINDOWS: (&str, &str) = ("LOSAT_X_BPREDOPARWINDOW", "1");

/// The queries of the searches, in a directory of the test `tag` (the tests run at once; the test
/// removes it): four consecutive records of `MeenMJNV.faa`, BDT62581.1 (698 letters, one
/// subject), BDT62582.1 (268, another), BDT62583.1 (629, no rows but a preliminary list) and
/// BDT62584.1 (944, five subjects), then BDT62553.1 (174), which has no preliminary list. Each
/// of the four has one; one search costs 0.1-0.2 CPU s (one of all 101 queries 3-9 s).
fn queries(tag: &str) -> PathBuf {
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let records: Vec<String> = text
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| format!(">{record}"))
        .collect();
    let first = records
        .iter()
        .position(|record| record.starts_with(">BDT62581.1 "))
        .expect("the record is in the fixture");
    let no_list = records
        .iter()
        .find(|record| record.starts_with(">BDT62553.1 "))
        .expect("the record is in the fixture");
    let dir = std::env::temp_dir().join(format!("losat-bppcompact-q-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("queries.faa");
    std::fs::write(&path, records[first..first + 4].concat() + no_list).unwrap();
    path
}

/// One blastp run of `query` against the subject proteome (the fixture pair is two related
/// nimavirus proteomes) with the given variables as its only `LOSAT_*` and query-batch variables,
/// and `extra` arguments.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the tests is named.
/// The run must exit with 0. Returns the standard output, the numbers of lists that the
/// `[TIMING] x_bppcompact` lines report (one per search of a query set), and the windows of each
/// `[TIMING] x_bpredopar:` line (one per search whose redo ran in windows).
fn blastp(
    query: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, Vec<usize>, Vec<usize>) {
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
        .arg(query)
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
    let redo_windows = stderr
        .lines()
        .filter(|line| line.starts_with("[TIMING] x_bpredopar: "))
        .map(|line| {
            line.split(' ')
                .find_map(|field| field.strip_prefix("windows="))
                .and_then(|count| count.parse().ok())
                .unwrap_or_else(|| panic!("windows= in {line}"))
        })
        .collect();
    (output.stdout, expanded, redo_windows)
}

#[test]
fn compact_lists_give_the_reports_of_the_one_thread_search() {
    let query = queries("full");
    // One thread: the redo expands each query's list when it reaches it, in one search of all
    // queries (one line, a list per query that has one; the loop passes BDT62553.1 by).
    let (reference, expanded, _) = blastp(&query, 1, &[BATCH, BLOCKS_OFF], &["-outfmt", "6"]);
    assert!(!reference.is_empty(), "the queries have a report");
    assert_eq!(expanded, [4], "threads 1");
    // Four threads: the window of the redo expands the lists as it reaches them, in the two blocks
    // of `TWO_BLOCKS` (one window each), one line per block. The report is the one-thread
    // search's.
    let (rows, expanded, redo_windows) = blastp(&query, 4, &[BATCH, TWO_BLOCKS], &["-outfmt", "6"]);
    assert!(
        rows == reference,
        "threads 4, two blocks: the report differs"
    );
    assert_eq!(
        expanded,
        [3, 1],
        "threads 4, two blocks: one line per block"
    );
    assert_eq!(
        redo_windows,
        [1, 1],
        "threads 4, two blocks: the redo in windows"
    );
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

#[test]
fn compact_lists_of_query_chunks_give_the_same_rows() {
    let query = queries("chunks");
    // Each chunk's search (`CHUNKS`: two of each split batch) expands its lists, one per query
    // part (five), before the chunk merge; at four threads the windows of the redo
    // (`SMALL_WINDOWS`, several) then read the merged lists beside the compact list of
    // BDT62581.1, which they expand: six in all. The reference expands every list before the
    // redo (`LOSAT_X_BPREDOWINEXPAND=0`) instead of window by window.
    let (reference, _, _) = blastp(
        &query,
        4,
        &[
            BATCH,
            CHUNKS[0],
            CHUNKS[1],
            SMALL_WINDOWS,
            BLOCKS_OFF,
            ("LOSAT_X_BPREDOWINEXPAND", "0"),
        ],
        &["-outfmt", "6"],
    );
    assert!(!reference.is_empty(), "the queries have a report");
    let vars = &[BATCH, CHUNKS[0], CHUNKS[1], SMALL_WINDOWS, BLOCKS_OFF][..];
    let (rows, expanded, redo_windows) = blastp(&query, 4, vars, &["-outfmt", "6"]);
    assert!(rows == reference, "{vars:?}: rows with query chunks differ");
    assert_eq!(expanded, [6], "{vars:?}");
    assert!(
        redo_windows.len() == 1 && redo_windows[0] > 1,
        "{vars:?}: several windows, {redo_windows:?}"
    );
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}
