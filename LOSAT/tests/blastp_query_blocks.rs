//! The blastp search in blocks of whole query batches with `LOSAT_X_BPQBLOCK`
//! (`src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switch is not in the strict set (`src/utils/xswitch.rs`): it is off unless its variable
//! is a positive number of query letters. With it on, a search whose only report is tabular
//! format 6 searches its queries block by block and writes each block's rows before it searches
//! the next block. These tests compare the rows with the search of all queries (the variable
//! unset) and check, from the `[TIMING] blastp_total` lines that each search of a query set
//! prints with `LOSAT_TIMING`, that the block path ran (more searches than the reference) or did
//! not (the same number), not only that the report is equal. `BATCH_SIZE` (NCBI's query batch
//! size) makes small batches, so that the fixture pair has many batches and blocks.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-287
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             ...
//!                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!             ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const THREADS: [u32; 2] = [1, 4];

/// One blastp run of the fixture pair (two related nimavirus proteomes) with the given
/// variables as its only `LOSAT_*` variables that choose a code path, and `extra` arguments.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing; without the caller's value a file of the tests is
/// named. `BATCH_SIZE` is the test's own. The run must exit
/// with 0. Returns the standard output and the number of `[TIMING] blastp_total` lines.
fn blastp(threads: u32, vars: &[(&str, &str)], extra: &[&str]) -> (Vec<u8>, usize) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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
            std::env::temp_dir().join("losat-bpqblock-worker.log"),
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
    let searches = stderr
        .lines()
        .filter(|line| line.starts_with("[TIMING] blastp_total:"))
        .count();
    (output.stdout, searches)
}

#[test]
fn blocks_write_the_rows_of_the_search_of_all_queries() {
    for batch_size in [None, Some("3000")] {
        let batch: Vec<(&str, &str)> = batch_size
            .map(|size| ("BATCH_SIZE", size))
            .into_iter()
            .collect();
        for threads in THREADS {
            let (reference, reference_searches) = blastp(threads, &batch, &["-outfmt", "6"]);
            assert!(!reference.is_empty(), "the fixture pair has rows");
            for letters in ["1", "20000"] {
                let mut vars = batch.clone();
                vars.push(("LOSAT_X_BPQBLOCK", letters));
                let (rows, searches) = blastp(threads, &vars, &["-outfmt", "6"]);
                assert!(
                    rows == reference,
                    "BATCH_SIZE {batch_size:?} threads {threads} blocks of {letters}: rows differ"
                );
                assert!(
                    searches > reference_searches,
                    "BATCH_SIZE {batch_size:?} threads {threads} blocks of {letters}: \
                     {searches} searches, the reference {reference_searches}"
                );
            }
        }
    }
}

// blastp-fast's compressed lookup table: each block builds its own (exactness argument for the
// compressed scan in `x_query_blocks.rs`).
#[test]
fn blastp_fast_blocks_write_the_rows_of_the_search_of_all_queries() {
    let fast = ["-task", "blastp-fast", "-outfmt", "6"];
    for batch_size in [None, Some("3000")] {
        let batch: Vec<(&str, &str)> = batch_size
            .map(|size| ("BATCH_SIZE", size))
            .into_iter()
            .collect();
        for threads in THREADS {
            let (reference, reference_searches) = blastp(threads, &batch, &fast);
            assert!(!reference.is_empty(), "the fixture pair has rows");
            for letters in ["1", "20000"] {
                let mut vars = batch.clone();
                vars.push(("LOSAT_X_BPQBLOCK", letters));
                let (rows, searches) = blastp(threads, &vars, &fast);
                assert!(
                    rows == reference,
                    "blastp-fast BATCH_SIZE {batch_size:?} threads {threads} blocks of {letters}: \
                     rows differ"
                );
                assert!(
                    searches > reference_searches,
                    "blastp-fast BATCH_SIZE {batch_size:?} threads {threads} blocks of {letters}: \
                     {searches} searches, the reference {reference_searches}"
                );
            }
        }
    }
}

#[test]
fn blocks_write_the_rows_to_one_output_file() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblock-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("rows.tsv");
    let path_text = path.to_str().unwrap();
    let custom = "6 qseqid sseqid pident length evalue bitscore qlen slen";
    let (reference, _) = blastp(1, &[("BATCH_SIZE", "3000")], &["-outfmt", custom]);
    let (stdout, searches) = blastp(
        1,
        &[("BATCH_SIZE", "3000"), ("LOSAT_X_BPQBLOCK", "1")],
        &["-outfmt", custom, "-out", path_text],
    );
    assert!(stdout.is_empty());
    assert!(searches > 1);
    assert!(
        std::fs::read(&path).unwrap() == reference,
        "the -out file differs"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn other_reports_and_one_block_take_the_search_of_all_queries() {
    let batch = ("BATCH_SIZE", "3000");
    for (letters, outfmt) in [("1", "7"), ("1", "0"), ("100000000", "6")] {
        let (reference, reference_searches) = blastp(1, &[batch], &["-outfmt", outfmt]);
        let (report, searches) = blastp(
            1,
            &[batch, ("LOSAT_X_BPQBLOCK", letters)],
            &["-outfmt", outfmt],
        );
        assert!(
            report == reference,
            "outfmt {outfmt} blocks of {letters}: the report differs"
        );
        assert_eq!(
            searches, reference_searches,
            "outfmt {outfmt} blocks of {letters}: not the search of all queries"
        );
    }
}
