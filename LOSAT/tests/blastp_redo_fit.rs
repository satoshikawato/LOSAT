//! blastp's redone HSP lists with `LOSAT_X_BPREDOFIT` (`src/algorithm/blastp/x_redo_fit.rs`) on
//! whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`. With it on, the
//! composition-based redo builds each redone HSP list with room for exactly its alignments. Every
//! test states the default it runs with: the reference is the search of all queries with the
//! switch off (`LOSAT_X_BPQBLOCK=0 LOSAT_X_BPREDOFIT=0`), and the search with the switch on is the
//! default (with `LOSAT_X_BPQBLOCK=0` where one line per search is counted);
//! `LOSAT_X_BPREDOFIT=0` gives no line. These tests
//! compare the reports, alone and with the block search of `LOSAT_X_BPQBLOCK`, at one and four
//! threads, and check from the
//! `[TIMING] x_bpredofit: lists=<n> hsps=<m>` lines that lists were built that way, not only that
//! the report is equal.
//! The switch acts on each query's lists, so three queries of the fixture are enough (`queries`),
//! each its own query batch (`BATCH`); the block test sets the blocks' letters so that there are
//! two and three blocks. `blastp_redo_par.rs` runs all the fixture's queries (batches of 3,000
//! letters) at one and four threads with the switch on by default: a query without a preliminary
//! list and batches of many queries.
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

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPREDOFIT", "0");
const BLOCKS_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");
/// Each query its own batch: the three queries of `queries` have 629, 944 and 918 letters.
const BATCH: (&str, &str) = ("BATCH_SIZE", "500");

/// Three records of the fixture's query proteome (`MeenMJNV.faa`) from BDT62583.1 on, written to a
/// directory of the test `tag`: BDT62583.1 (629 letters), which has preliminary matches but none
/// within the E-value after the redo (no hits), then BDT62584.1 (944) and BDT62585.1 (918), which
/// have hits.
fn queries(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpredofit-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let start = text
        .find(">BDT62583.1 ")
        .expect("the record is in the fixture");
    let mut end = start;
    for _ in 0..3 {
        end = text[end + 1..]
            .find("\n>")
            .map_or(text.len(), |at| end + 1 + at + 1);
    }
    let path = dir.join("queries.faa");
    std::fs::write(&path, &text[start..end]).unwrap();
    path
}

/// One blastp run of `query` against the fixture's subject proteome (`MejoMJNV.faa`, a nimavirus
/// related to the queries') with the given variables as its only `LOSAT_*` variables that choose
/// a code path, and `extra` arguments.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the tests is named.
/// `BATCH_SIZE` is the test's own. The run must exit with 0. Returns the standard output and the
/// (lists, HSPs) that the `[TIMING] x_bpredofit` lines report, one per search of a query set.
fn blastp(
    query: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, Vec<(usize, usize)>) {
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
    let query = queries("reports");
    for outfmt in ["6", "0"] {
        for threads in [1, 4] {
            let (reference, reference_fitted) = blastp(
                &query,
                threads,
                &[BATCH, BLOCKS_OFF, SWITCH_OFF],
                &["-outfmt", outfmt],
            );
            assert!(!reference.is_empty(), "the queries have a report");
            assert!(
                reference_fitted.is_empty(),
                "no line with the switch off: {reference_fitted:?}"
            );
            let vars = &[BATCH, BLOCKS_OFF][..];
            let (report, fitted) = blastp(&query, threads, vars, &["-outfmt", outfmt]);
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

#[test]
fn fitted_redone_lists_in_query_blocks_give_the_rows_of_the_search_of_all_queries() {
    let query = queries("blocks");
    for threads in [1, 4] {
        let (reference, _) = blastp(
            &query,
            threads,
            &[BATCH, BLOCKS_OFF, SWITCH_OFF],
            &["-outfmt", "6"],
        );
        // Blocks of 1,500 letters (BDT62583.1 and BDT62584.1, two batches, then BDT62585.1), as the
        // default blocks of 50,000 letters hold several batches; and one batch per block.
        for (vars, blocks) in [
            (&[BATCH, ("LOSAT_X_BPQBLOCK", "1500")][..], 2),
            (&[BATCH, ("LOSAT_X_BPQBLOCK", "1")][..], 3),
        ] {
            let (rows, fitted) = blastp(&query, threads, vars, &["-outfmt", "6"]);
            assert!(rows == reference, "threads {threads} {vars:?}: rows differ");
            assert_eq!(
                fitted.len(),
                blocks,
                "threads {threads} {vars:?}: one line per block, {fitted:?}"
            );
            assert!(
                fitted.iter().map(|(lists, _)| lists).sum::<usize>() > 0,
                "threads {threads} {vars:?}: no list fitted: {fitted:?}"
            );
        }
    }
}
