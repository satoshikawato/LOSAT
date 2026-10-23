// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! The blastp search in blocks of whole query batches with `LOSAT_X_BPQBLOCK`
//! (`src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, with blocks of
//! 50,000 query letters, unless its variable is `0`, or `LOSAT_X_DEFAULT=0` and its variable is
//! not set; a positive number in the variable sets the size of a block. With it on, a search
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

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const BATCH: (&str, &str) = ("BATCH_SIZE", "3000");
const EXACT: &[&str] = &["-outfmt", "6"];
const FAST: &[&str] = &["-task", "blastp-fast", "-outfmt", "6"];

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
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with.
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

/// Each case (the task's arguments, the threads, the block letters): its rows in blocks are the
/// rows of the search of all queries at the same thread count, from more searches. `BATCH_SIZE`
/// 3000 makes 21 query batches.
fn assert_block_rows(subject: &[&OsStr], cases: &[(&[&str], u32, &str)]) {
    for &(task, threads, letters) in cases {
        let what = format!(
            "{:?} {task:?} threads {threads} blocks of {letters}",
            subject[0]
        );
        let (reference, reference_searches) =
            blastp_against(subject, threads, &[BATCH, SWITCH_OFF], task);
        assert!(!reference.is_empty(), "the fixture pair has rows");
        let (rows, searches) = blastp_against(
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", letters)],
            task,
        );
        assert!(rows == reference, "{what}: rows differ");
        assert!(
            searches > reference_searches,
            "{what}: {searches} searches, the reference {reference_searches}"
        );
    }
}

// The rows of the block search are the rows of the search of all queries. One case per path (each
// search of the fixture pair costs CI seconds): each lookup table (exact blastp here, blastp-fast's
// compressed table below), one thread and four (the pool encodes the subjects for all blocks,
// `x_bp_encode_block_subjects`, and computes the contexts' Karlin blocks), blocks of one batch (`1`)
// and of several (`20000`). Each combination of the three runs once: with `-subject` here and
// below, and with `-db` in `database_blocks_write_the_rows_of_the_search_of_all_queries`. NCBI's
// default batches (10,000 letters) take the same paths with fewer blocks: these queries split none
// of them (a split needs a batch of at least 19,800 letters, `calculate_num_chunks`).
#[test]
fn blocks_write_the_rows_of_the_search_of_all_queries() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    assert_block_rows(&subject, &[(EXACT, 1, "1"), (EXACT, 4, "20000")]);
}

// blastp-fast's compressed lookup table: each block builds its own (exactness argument for the
// compressed scan in `x_query_blocks.rs`), after the test of its batches' words
// (`x_block_words_fit`, in groups on the pool at four threads).
#[test]
fn blastp_fast_blocks_write_the_rows_of_the_search_of_all_queries() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    assert_block_rows(&subject, &[(FAST, 1, "20000"), (FAST, 4, "1")]);
}

#[test]
fn blocks_write_the_rows_to_one_output_file() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblock-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("rows.tsv");
    let path_text = path.to_str().unwrap();
    let custom = "6 qseqid sseqid pident length evalue bitscore qlen slen";
    let (reference, _) = blastp(1, &[BATCH, SWITCH_OFF], &["-outfmt", custom]);
    let (stdout, searches) = blastp(
        1,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
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

// outfmt 7 with `LOSAT_X_BPQBLOCKFMT=0` (with it on, `blastp_query_blocks_fmt.rs`, which also has
// outfmt 0 with it `0`), and one block, take the search of all queries.
#[test]
fn other_reports_and_one_block_take_the_search_of_all_queries() {
    let fmt_off = ("LOSAT_X_BPQBLOCKFMT", "0");
    for (letters, outfmt) in [("1", "7"), ("100000000", "6")] {
        let (reference, reference_searches) = blastp(1, &[BATCH, SWITCH_OFF], &["-outfmt", outfmt]);
        let (report, searches) = blastp(
            1,
            &[BATCH, ("LOSAT_X_BPQBLOCK", letters), fmt_off],
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
// 77,448: two blocks), and `LOSAT_X_BPQBLOCK=0` and `LOSAT_X_DEFAULT=0` search all queries
// together. How a value switch reads its variable (a value that is not a positive number is the
// built-in size, `LOSAT_X_DEFAULT=0` leaves a variable that is set on) is `state` and `given`,
// unit-tested in `src/utils/xswitch.rs`; the default with `-db`, and a variable set under
// `LOSAT_X_DEFAULT=0`, run in `database_search_with_the_exact_blastp_experiments`.
#[test]
fn the_default_searches_in_blocks_of_50000_letters() {
    let threads = 4;
    let (reference, reference_searches) = blastp(threads, &[BATCH, SWITCH_OFF], EXACT);
    assert!(!reference.is_empty(), "the fixture pair has rows");
    let (rows, searches) = blastp(threads, &[BATCH, DEFAULT_OFF], EXACT);
    assert!(rows == reference, "LOSAT_X_DEFAULT=0: rows differ");
    assert_eq!(searches, reference_searches, "LOSAT_X_DEFAULT=0");
    let (rows, searches) = blastp(threads, &[BATCH], EXACT);
    assert!(rows == reference, "the default: rows differ");
    assert_eq!(
        searches,
        reference_searches + 1,
        "the default: two blocks of at least 50,000 letters"
    );
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
// definition line of the database can fail a batch (`x_titles_allow_blocks`, on the pool at four
// threads). The combinations of lookup table, threads and block letters that `-subject` does not
// run above.
#[test]
fn database_blocks_write_the_rows_of_the_search_of_all_queries() {
    let db = subject_database("rows");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    assert_block_rows(
        &subject,
        &[
            (EXACT, 1, "20000"),
            (EXACT, 4, "1"),
            (FAST, 1, "1"),
            (FAST, 4, "20000"),
        ],
    );
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}

// The exact blastp switches together with `-db`, in blocks and in the search of all queries: the
// rows of the search with every switch of the strict set off (`LOSAT_X_DEFAULT=0`). The default
// runs them all, in blocks of 50,000 letters (two blocks here); `LOSAT_X_DEFAULT=0` with the
// variables set gives the runs they gave before they joined the set. `LOSAT_X_BPPCOMPACT=shadow`
// compares every compact preliminary list of the database search with the off path's list (it
// panics on a difference). One case per path: the runs at four threads, and at one thread the
// experiments in blocks and the shadow (whose `-db` subject loop at one thread is sequential,
// `database.rs` `search_blocks`). The experiments' thread gates are left to their own tests; the
// default in blocks at one thread with `-db`: `database_blocks_write_the_rows_...` above.
#[test]
fn database_search_with_the_exact_blastp_experiments() {
    const EXPERIMENTS: [(&str, &str); 7] = [
        ("LOSAT_X_BPSEED", "1"),
        ("LOSAT_X_BPUNGAP", "1"),
        ("LOSAT_X_BPGAP", "1"),
        ("LOSAT_X_BPLUT", "1"),
        ("LOSAT_X_BPPFREE", "1"),
        ("LOSAT_X_BPPCOMPACT", "1"),
        ("LOSAT_X_BPREDOFIT", "1"),
    ];
    let db = subject_database("exact");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let experiments = |blocks: Option<&'static str>| {
        let mut vars = EXPERIMENTS.to_vec();
        vars.extend([BATCH, DEFAULT_OFF]);
        vars.extend(blocks.map(|letters| ("LOSAT_X_BPQBLOCK", letters)));
        vars
    };
    let shadow = ("LOSAT_X_BPPCOMPACT", "shadow");
    /// A run's variables and whether the blocks are searched (`None`: a shadow run, whose searches
    /// are not counted).
    type Run<'a> = (Vec<(&'a str, &'a str)>, Option<bool>);
    let cases: [(u32, Vec<Run>); 2] = [
        (
            4,
            vec![
                (vec![BATCH], Some(true)),
                (vec![BATCH, SWITCH_OFF], Some(false)),
                (experiments(None), Some(false)),
                (vec![BATCH, shadow], None),
                (vec![BATCH, DEFAULT_OFF, shadow], None),
            ],
        ),
        (
            1,
            vec![
                (experiments(Some("1")), Some(true)),
                (vec![BATCH, shadow], None),
            ],
        ),
    ];
    for (threads, runs) in &cases {
        let threads = *threads;
        let (reference, reference_searches) =
            blastp_against(&subject, threads, &[BATCH, DEFAULT_OFF], EXACT);
        for (vars, blocks) in runs {
            let (rows, searches) = blastp_against(&subject, threads, vars, EXACT);
            assert!(
                rows == reference,
                "-db threads {threads} {vars:?}: rows differ"
            );
            if let Some(blocks) = blocks {
                assert_eq!(
                    searches > reference_searches,
                    *blocks,
                    "-db threads {threads} {vars:?}: {searches} searches, the reference \
                     {reference_searches}"
                );
            }
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}
