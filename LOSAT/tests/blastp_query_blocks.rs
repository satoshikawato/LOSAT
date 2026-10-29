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
//! `LOSAT_TIMING`, that the block path ran (a search per block, more than the reference's) or did
//! not (the same number), not only that the report is equal. `BATCH_SIZE` (NCBI's query batch
//! size) makes small batches, so that a few queries make several batches and blocks. The same
//! comparisons run with `-db` on a database of the subject proteome made with `losat makedb`.
//!
//! The searches run four of the fixture's 101 queries (`queries`), except the default test's:
//! blocks of 50,000 letters need all of them. The block plan reads no other size whose outcome the
//! four change (it needs more than one block, the `Int8` search space of the database and the
//! longest query, the diagonal table's window, the batches' split sizes and, for blastp-fast, their
//! words, `x_block_words_fit`). Three rules of other default switches that read sizes are not
//! reached by the four: a build ahead that lacks the room (`LOSAT_X_BPLKPAHEAD`:
//! `blastp_lookup_ahead.rs`), several redo windows that end inside queries (`LOSAT_X_BPREDOPAR`:
//! the default test, `blastp_redo_par.rs`) and the DP band past lane 1024 (`LOSAT_X_DPWINDOW`: the
//! default test, the reach gate's fixtures).
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
/// The batches of `queries`: {BDT62581.1}, {BDT62582.1, BDT62583.1}, {BDT62584.1} (a batch takes
/// records while it holds fewer letters).
const BATCH: (&str, &str) = ("BATCH_SIZE", "500");
/// `LOSAT_X_BPQBLOCK` for blocks of several batches of `queries` (a block takes batches while it
/// holds fewer letters): {BDT62581.1 to BDT62583.1} (1,595 letters) and {BDT62584.1}. `1` makes a
/// block per batch, three.
const TWO_BLOCKS: &str = "1000";
const EXACT: &[&str] = &["-outfmt", "6"];
const FAST: &[&str] = &["-task", "blastp-fast", "-outfmt", "6"];

/// The queries of the searches here but the default test's, in a directory of the test `tag` (the
/// tests run at once; the test removes it): four consecutive records of `MeenMJNV.faa`, BDT62581.1
/// (698 letters, one subject), BDT62582.1 (268, another), BDT62583.1 (629, no hits) and BDT62584.1
/// (944, five subjects), three batches of `BATCH`. One search costs 0.1-0.25 CPU s (one of all 101
/// queries 2-8 s).
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
    let dir = std::env::temp_dir().join(format!("losat-bpqblock-q-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("queries.faa");
    std::fs::write(&path, records[first..first + 4].concat()).unwrap();
    path
}

/// The blocks of `queries` in blocks of `letters`.
fn blocks(letters: &str) -> usize {
    match letters {
        "1" => 3,
        TWO_BLOCKS => 2,
        _ => unreachable!("blocks of {letters}"),
    }
}

/// One blastp run of `query` against the subject proteome (the fixture pair is two related
/// nimavirus proteomes) with the given variables as its only `LOSAT_*` variables that choose a code
/// path, and `extra` arguments.
fn blastp(query: &Path, threads: u32, vars: &[(&str, &str)], extra: &[&str]) -> (Vec<u8>, usize) {
    let subject = fixture("MejoMJNV.faa");
    blastp_against(
        query,
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
    query: &Path,
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
        .arg(query)
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

/// Each case (the task's arguments, the threads, the block letters): its rows in blocks of
/// `queries` are the rows of the search of all of them at the same thread count, from one search
/// per block (the reference: one).
fn assert_block_rows(query: &Path, subject: &[&OsStr], cases: &[(&[&str], u32, &str)]) {
    for &(task, threads, letters) in cases {
        let what = format!(
            "{:?} {task:?} threads {threads} blocks of {letters}",
            subject[0]
        );
        let (reference, reference_searches) =
            blastp_against(query, subject, threads, &[BATCH, SWITCH_OFF], task);
        assert!(!reference.is_empty(), "the fixture pair has rows");
        assert_eq!(reference_searches, 1, "{what}: the reference");
        let (rows, searches) = blastp_against(
            query,
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", letters)],
            task,
        );
        assert!(rows == reference, "{what}: rows differ");
        assert_eq!(searches, blocks(letters), "{what}: searches");
    }
}

// The rows of the block search are the rows of the search of all queries. One case per path (CI
// time): each lookup table (exact blastp here, blastp-fast's compressed table below), one thread
// and four (the pool encodes the subjects for all blocks, `x_bp_encode_block_subjects`, and
// computes the contexts' Karlin blocks), blocks of one batch (`1`) and of several (`TWO_BLOCKS`).
// Each combination of the three runs once: with `-subject` here and below, and with `-db` in
// `database_blocks_write_the_rows_of_the_search_of_all_queries`. NCBI's default batches (10,000
// letters) take the same paths with fewer blocks (the default test: the fixture's queries split
// none of them, a split needs a batch of at least 19,800 letters, `calculate_num_chunks`).
#[test]
fn blocks_write_the_rows_of_the_search_of_all_queries() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let query = queries("rows");
    assert_block_rows(&query, &subject, &[(EXACT, 1, "1"), (EXACT, 4, TWO_BLOCKS)]);
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

// blastp-fast's compressed lookup table: each block builds its own (exactness argument for the
// compressed scan in `x_query_blocks.rs`), after the test of its batches' words
// (`x_block_words_fit`, in groups on the pool at four threads: three groups of one batch here; one
// counter for every batch at one thread).
#[test]
fn blastp_fast_blocks_write_the_rows_of_the_search_of_all_queries() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let query = queries("fast");
    assert_block_rows(&query, &subject, &[(FAST, 1, TWO_BLOCKS), (FAST, 4, "1")]);
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

#[test]
fn blocks_write_the_rows_to_one_output_file() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblock-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("rows.tsv");
    let path_text = path.to_str().unwrap();
    let custom = "6 qseqid sseqid pident length evalue bitscore qlen slen";
    let query = queries("out");
    let (reference, _) = blastp(&query, 1, &[BATCH, SWITCH_OFF], &["-outfmt", custom]);
    let (stdout, searches) = blastp(
        &query,
        1,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
        &["-outfmt", custom, "-out", path_text],
    );
    assert!(stdout.is_empty());
    assert_eq!(searches, blocks("1"));
    assert!(
        std::fs::read(&path).unwrap() == reference,
        "the -out file differs"
    );
    std::fs::remove_dir_all(&dir).unwrap();
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

// outfmt 7 with `LOSAT_X_BPQBLOCKFMT=0` (with it on, `blastp_query_blocks_fmt.rs`, which also has
// outfmt 0 with it `0`), and one block, take the search of all queries.
#[test]
fn other_reports_and_one_block_take_the_search_of_all_queries() {
    let fmt_off = ("LOSAT_X_BPQBLOCKFMT", "0");
    let query = queries("other");
    for (letters, outfmt) in [("1", "7"), ("100000000", "6")] {
        let (reference, reference_searches) =
            blastp(&query, 1, &[BATCH, SWITCH_OFF], &["-outfmt", outfmt]);
        let (report, searches) = blastp(
            &query,
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
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

// The default: with no variable the block is 50,000 query letters (the fixture's queries have
// 77,448: two blocks), and `LOSAT_X_BPQBLOCK=0` searches all queries together. How a value
// switch reads its variable (a value that is not a positive number is the built-in size) is
// `state` and `given`, unit-tested in `src/utils/xswitch.rs`. The only test of all 101 queries
// here, in `BATCH_SIZE` 3000 (21 batches). `x_bp_qblock` is one value, read by the block plan for
// `-subject` and `-db` alike (`-db` in blocks: below).
#[test]
fn the_default_searches_in_blocks_of_50000_letters() {
    let threads = 4;
    let all = fixture("MeenMJNV.faa");
    let batch = ("BATCH_SIZE", "3000");
    let (reference, reference_searches) = blastp(&all, threads, &[batch, SWITCH_OFF], EXACT);
    assert!(!reference.is_empty(), "the fixture pair has rows");
    let (rows, searches) = blastp(&all, threads, &[batch], EXACT);
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
    let query = queries("db-rows");
    assert_block_rows(
        &query,
        &subject,
        &[
            (EXACT, 1, TWO_BLOCKS),
            (EXACT, 4, "1"),
            (FAST, 1, "1"),
            (FAST, 4, TWO_BLOCKS),
        ],
    );
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}

// The exact blastp search with `-db`, in blocks and in the search of all queries: the rows of the
// search with `LOSAT_X_BPGAP` and the block search off (`EXACT_OFF`) against the default's. The
// runs at the default's block size are in `TWO_BLOCKS` here (the default size needs all 101
// queries: `the_default_searches_in_blocks_of_50000_letters`). One case per path: the runs at four
// threads, and at one thread the default in blocks of one letter (whose `-db` subject loop at one
// thread is sequential, `database.rs` `search_blocks`). The experiments' thread gates are left to
// their own tests; the default in blocks at one thread with `-db`:
// `database_blocks_write_the_rows_...` above.
#[test]
fn database_search_with_the_exact_blastp_experiments() {
    const EXACT_OFF: [(&str, &str); 3] = [BATCH, ("LOSAT_X_BPGAP", "0"), SWITCH_OFF];
    let db = subject_database("exact");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let query = queries("exact");
    let two_blocks = ("LOSAT_X_BPQBLOCK", TWO_BLOCKS);
    /// A run's variables and its searches.
    type Run<'a> = (Vec<(&'a str, &'a str)>, usize);
    let cases: [(u32, Vec<Run>); 2] = [
        (
            4,
            vec![
                (vec![BATCH, two_blocks], blocks(TWO_BLOCKS)),
                (vec![BATCH, SWITCH_OFF], 1),
            ],
        ),
        (
            1,
            vec![(vec![BATCH, ("LOSAT_X_BPQBLOCK", "1")], blocks("1"))],
        ),
    ];
    for (threads, runs) in &cases {
        let threads = *threads;
        let (reference, reference_searches) =
            blastp_against(&query, &subject, threads, &EXACT_OFF, EXACT);
        assert_eq!(
            reference_searches, 1,
            "-db threads {threads}: the reference"
        );
        for (vars, expected) in runs {
            let (rows, searches) = blastp_against(&query, &subject, threads, vars, EXACT);
            assert!(
                rows == reference,
                "-db threads {threads} {vars:?}: rows differ"
            );
            assert_eq!(
                searches, *expected,
                "-db threads {threads} {vars:?}: searches"
            );
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(query.parent().unwrap()).unwrap();
}
