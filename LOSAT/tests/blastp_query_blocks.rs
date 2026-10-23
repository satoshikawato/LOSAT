//! The blastp search in blocks of whole query batches with `LOSAT_X_BPQBLOCK`
//! (`src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, with blocks of
//! 50,000 query letters, unless its variable is `0`; a positive number in the variable sets the
//! size of a block. With it on, a search
//! whose only report is tabular format 6 searches its queries block by block and writes each
//! block's rows before it searches the next block (with `LOSAT_X_BPQBLOCKFMT`, also outfmt 0 and
//! 7: `blastp_query_blocks_fmt.rs`). Every test states the default it runs with:
//! these tests compare the rows with the search of all queries (`LOSAT_X_BPQBLOCK=0`) and check,
//! from the `[TIMING] blastp_total` lines that each search of a query set prints with
//! `LOSAT_TIMING`, that the block path ran (more searches than the reference) or did not (the
//! same number), not only that the report is equal. `BATCH_SIZE` (NCBI's query batch size) makes
//! small batches, so that the fixture pair has many batches and blocks. The same
//! comparisons run with `-db` on a database of the subject proteome made with `losat makedb`.
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

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const THREADS: [u32; 2] = [1, 4];
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");

/// One blastp run of the fixture pair (two related nimavirus proteomes) with the given
/// variables as its only `LOSAT_*` variables that choose a code path, and `extra` arguments.
fn blastp(threads: u32, vars: &[(&str, &str)], extra: &[&str]) -> (Vec<u8>, usize) {
    let subject = fixture("MejoMJNV.faa");
    blastp_against(
        &[OsStr::new("-subject"), subject.as_os_str()],
        threads,
        vars,
        extra,
    )
}

/// `blastp` with the subject arguments `subject` (`-subject FILE` or `-db NAME`).
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing; without the caller's value a file of the tests is
/// named. `BATCH_SIZE` is the test's own. The run must exit
/// with 0. Returns the standard output and the number of `[TIMING] blastp_total` lines.
fn blastp_against(
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, usize) {
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
            std::env::temp_dir().join("losat-bpqblock-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
        .args(subject)
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
            let mut off = batch.clone();
            off.push(SWITCH_OFF);
            let (reference, reference_searches) = blastp(threads, &off, &["-outfmt", "6"]);
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
            let mut off = batch.clone();
            off.push(SWITCH_OFF);
            let (reference, reference_searches) = blastp(threads, &off, &fast);
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
    let (reference, _) = blastp(
        1,
        &[("BATCH_SIZE", "3000"), SWITCH_OFF],
        &["-outfmt", custom],
    );
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

// outfmt 0 and 7 with `LOSAT_X_BPQBLOCKFMT=0` (with it on, `blastp_query_blocks_fmt.rs`), and one
// block, take the search of all queries.
#[test]
fn other_reports_and_one_block_take_the_search_of_all_queries() {
    let batch = ("BATCH_SIZE", "3000");
    let fmt_off = ("LOSAT_X_BPQBLOCKFMT", "0");
    for (letters, outfmt) in [("1", "7"), ("1", "0"), ("100000000", "6")] {
        let (reference, reference_searches) = blastp(1, &[batch, SWITCH_OFF], &["-outfmt", outfmt]);
        let (report, searches) = blastp(
            1,
            &[batch, ("LOSAT_X_BPQBLOCK", letters), fmt_off],
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

// The default: with no variable the block is 50,000 query letters (the fixture's queries have
// 77,448), `LOSAT_X_BPQBLOCK=0` searches all queries together, and a value that is not a
// positive number is the built-in size. The same with `-db`.
#[test]
fn the_default_searches_in_blocks_of_50000_letters() {
    let batch = ("BATCH_SIZE", "3000");
    let fmt6 = ["-outfmt", "6"];
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("default");
    for subject in [
        [OsStr::new("-subject"), fasta.as_os_str()],
        [OsStr::new("-db"), db.as_os_str()],
    ] {
        for threads in THREADS {
            let what = format!("{:?} threads {threads}", subject[0]);
            let (reference, reference_searches) =
                blastp_against(&subject, threads, &[batch, SWITCH_OFF], &fmt6);
            assert!(!reference.is_empty(), "the fixture pair has rows");
            for vars in [
                &[batch][..],
                &[batch, ("LOSAT_X_BPQBLOCK", "abc")][..],
                &[batch, ("LOSAT_X_BPQBLOCK", "50000")][..],
            ] {
                let (rows, searches) = blastp_against(&subject, threads, vars, &fmt6);
                assert!(rows == reference, "{what} {vars:?}: rows differ");
                assert_eq!(
                    searches,
                    reference_searches + 1,
                    "{what} {vars:?}: two blocks of at least 50,000 letters"
                );
            }
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}

/// A protein database of the subject proteome, made with `losat makedb` in a directory of the
/// test (`tag` names it); returns the database name.
fn subject_database(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpqblock-db-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

// `-db`: the database's subjects are searched one after the other on one diagonal table
// (exactness argument in `x_query_blocks.rs`), and the blocks are searched only when no
// definition line of the database can fail a batch (`x_titles_allow_blocks`).
#[test]
fn database_blocks_write_the_rows_of_the_search_of_all_queries() {
    let db = subject_database("rows");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    for task in [
        &["-outfmt", "6"][..],
        &["-task", "blastp-fast", "-outfmt", "6"],
    ] {
        for batch_size in [None, Some("3000")] {
            let batch: Vec<(&str, &str)> = batch_size
                .map(|size| ("BATCH_SIZE", size))
                .into_iter()
                .collect();
            for threads in THREADS {
                let mut off = batch.clone();
                off.push(SWITCH_OFF);
                let (reference, reference_searches) = blastp_against(&subject, threads, &off, task);
                assert!(!reference.is_empty(), "the fixture pair has rows");
                for letters in ["1", "20000"] {
                    let mut vars = batch.clone();
                    vars.push(("LOSAT_X_BPQBLOCK", letters));
                    let (rows, searches) = blastp_against(&subject, threads, &vars, task);
                    assert!(
                        rows == reference,
                        "-db {task:?} BATCH_SIZE {batch_size:?} threads {threads} blocks of \
                         {letters}: rows differ"
                    );
                    assert!(
                        searches > reference_searches,
                        "-db {task:?} BATCH_SIZE {batch_size:?} threads {threads} blocks of \
                         {letters}: {searches} searches, the reference {reference_searches}"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}

// The exact blastp switches together with `-db`, in blocks and in the search of all queries: the
// rows of the search with the exact blastp switches and the block search off. The default runs
// them all, in blocks of 50,000 letters (two blocks here); `LOSAT_X_BPPCOMPACT=shadow`
// compares every compact preliminary list of the database search with the off path's list (it
// panics on a difference).
#[test]
fn database_search_with_the_exact_blastp_experiments() {
    const EXACT_OFF: [(&str, &str); 8] = [
        ("LOSAT_X_BPSEED", "0"),
        ("LOSAT_X_BPUNGAP", "0"),
        ("LOSAT_X_BPGAP", "0"),
        ("LOSAT_X_BPLUT", "0"),
        ("LOSAT_X_BPPFREE", "0"),
        ("LOSAT_X_BPPCOMPACT", "0"),
        ("LOSAT_X_BPREDOFIT", "0"),
        ("LOSAT_X_BPQBLOCK", "0"),
    ];
    let db = subject_database("exact");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let batch = ("BATCH_SIZE", "3000");
    let outfmt = ["-outfmt", "6"];
    let runs: Vec<(Vec<(&str, &str)>, bool)> = vec![
        (vec![batch], true),
        (vec![batch, SWITCH_OFF], false),
        (vec![batch, ("LOSAT_X_BPQBLOCK", "1")], true),
    ];
    let mut exact_off = EXACT_OFF.to_vec();
    exact_off.push(batch);
    for threads in THREADS {
        let (reference, reference_searches) =
            blastp_against(&subject, threads, &exact_off, &outfmt);
        for (vars, blocks) in &runs {
            let (rows, searches) = blastp_against(&subject, threads, vars, &outfmt);
            assert!(
                rows == reference,
                "-db threads {threads} {vars:?}: rows differ"
            );
            assert_eq!(
                searches > reference_searches,
                *blocks,
                "-db threads {threads} {vars:?}: {searches} searches, the reference \
                 {reference_searches}"
            );
        }
        for vars in [
            &[batch, ("LOSAT_X_BPPCOMPACT", "shadow")][..],
            &[batch, SWITCH_OFF, ("LOSAT_X_BPPCOMPACT", "shadow")][..],
        ] {
            let (rows, _) = blastp_against(&subject, threads, vars, &outfmt);
            assert!(
                rows == reference,
                "-db threads {threads} {vars:?}: rows differ"
            );
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}
