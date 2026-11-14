//! The preliminary stage of blastp's search of all queries in blocks of whole query batches with
//! `LOSAT_X_BPQBLOCK` (`src/algorithm/blastp/blast_engine.rs`: `x_prelim_block_plan` and the block
//! loop of `run_call`, item SC8b-2).
//!
//! The reports that the block search (`blastp_query_blocks.rs`) does not write come from the search
//! of all queries: outfmt 10, 18 and 20, a custom outfmt 7 that reads more of `-db`'s definition
//! lines, `-query_loc`, and others. That search runs its preliminary stage block by block, each
//! block a call with its own lookup table, and joins the blocks' hit lists before the redo and the
//! reports, which run over all queries as after one table. These tests compare the standard output,
//! the standard error (without its `[TIMING]` lines) and the exit status with the search of all
//! queries with one table (`LOSAT_X_BPQBLOCK=0`), and read from `LOSAT_X_REACH` that the blocks ran
//! (the switch's own side; the reference side is recorded too, by the block search's plan, which
//! does not take these reports) or did not. `BATCH_SIZE` (NCBI's query batch size) makes small
//! batches, so that a few queries make several batches and blocks.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             ...
//!                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//!                 results = lcl_blast.Run();
//! ```

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;
mod temp_support;
use temp_support::{worker_log_if_unset, TempFile, TempPath};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");
/// The batches of `queries(4)`: {BDT62581.1}, {BDT62582.1, BDT62583.1}, {BDT62584.1} (a batch takes
/// records while it holds fewer letters).
const BATCH: (&str, &str) = ("BATCH_SIZE", "500");
/// `LOSAT_X_BPQBLOCK` for blocks of several batches of `queries(4)`: {BDT62581.1 to BDT62583.1}
/// (1,595 letters) and {BDT62584.1}. `1` makes a block per batch, three.
const TWO_BLOCKS: &str = "1000";

/// One run.
struct Run {
    code: Option<i32>,
    stdout: Vec<u8>,
    /// Standard error without the `[TIMING]` lines.
    stderr: String,
    _log: Option<TempFile>,
}

/// `count` consecutive records of `MeenMJNV.faa` from BDT62581.1 (698 letters, one subject),
/// BDT62582.1 (268, another), BDT62583.1 (629, no hits), BDT62584.1 (944, five subjects), ... in a
/// directory of the test `tag` (removed with the returned guard).
fn queries(tag: &str, count: usize) -> (PathBuf, TempPath) {
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
    queries_of(tag, &records[first..first + count].concat())
}

/// A query file of `text` in a directory of the test `tag`.
fn queries_of(tag: &str, text: &str) -> (PathBuf, TempPath) {
    let dir = TempPath::new(
        std::env::temp_dir().join(format!("losat-bpprelim-q-{tag}-{}", std::process::id())),
    );
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("queries.faa");
    std::fs::write(&path, text).unwrap();
    (path, dir)
}

/// A protein database of the subject proteome, made with `losat makeblastdb` in a directory of the
/// test (`tag` names it); returns the database name and the directory's guard.
fn subject_database(tag: &str) -> (PathBuf, TempPath) {
    let dir = TempPath::new(
        std::env::temp_dir().join(format!("losat-bpprelim-db-{tag}-{}", std::process::id())),
    );
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (name, dir)
}

/// `blastp -query <query> <subject...> -evalue 0.001 -num_threads <threads> <extra...>`, run in
/// `LOSAT/`, with the given variables as its only `LOSAT_*` and batch variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs
/// that file name as soon as a search has more than one thread, and it selects nothing.
fn blastp(
    query: &Path,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
            || text.starts_with("NCBI")
            || text == "BLASTDB"
        {
            command.env_remove(name);
        }
    }
    let log = worker_log_if_unset(&mut command, "blastp_prelim_blocks_worker");
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
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
    let stderr = String::from_utf8_lossy(&output.stderr)
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[TIMING]"))
        .collect();
    Run {
        code: output.status.code(),
        stdout: output.stdout,
        stderr,
        _log: log,
    }
}

/// Runs the case with one table (`LOSAT_X_BPQBLOCK=0`) and with blocks of `letters`
/// (`LOSAT_X_BPQBLOCK=<letters>`, `None`: the switch's default of 50,000 letters), both with
/// `vars`; asserts that both write the same bytes and end the same way, and returns the second run
/// and its reach directory (`Reach::sides("LOSAT_X_BPQBLOCK")`: (the blocks ran, one table ran);
/// a run writes it when it ends with 0).
fn same_as_one_table(
    query: &Path,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    letters: Option<&str>,
    extra: &[&str],
) -> (reach_support::Reach, Run) {
    let what = format!(
        "{:?} {extra:?} {vars:?} threads {threads} blocks {letters:?}",
        subject[0]
    );
    let mut off = vars.to_vec();
    off.push(SWITCH_OFF);
    let reference = blastp(query, subject, threads, &off, extra);
    let reach = reach_support::Reach::new("bpprelim");
    let mut on = vars.to_vec();
    on.push(reach.var());
    if let Some(letters) = letters {
        on.push(("LOSAT_X_BPQBLOCK", letters));
    }
    let run = blastp(query, subject, threads, &on, extra);
    assert_eq!(run.code, reference.code, "{what}: {}", run.stderr);
    assert!(run.stdout == reference.stdout, "{what}: the report differs");
    assert_eq!(
        run.stderr, reference.stderr,
        "{what}: standard error differs"
    );
    (reach, run)
}

// The reports that only the search of all queries writes, with `-subject`, in blocks of the
// preliminary stage: each lookup table (exact blastp's words of 3 and of 4, blastp-fast's
// compressed table), one thread and four, blocks of one batch (`1`) and of several (`TWO_BLOCKS`);
// outfmt 10 (CSV), 18 (the organism report) and 20 (CSV with a header line). One case per path (CI
// time).
#[test]
fn reports_of_the_search_of_all_queries_with_blocks_of_the_preliminary_stage() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let (query, _dir) = queries("subject", 4);
    let cases: [(&[&str], u32, &str); 6] = [
        (&["-outfmt", "10"], 1, "1"),
        (&["-outfmt", "18"], 4, TWO_BLOCKS),
        (&["-outfmt", "20 qseqid sseqid evalue bitscore"], 4, "1"),
        (&["-task", "blastp-fast", "-outfmt", "10"], 4, TWO_BLOCKS),
        (&["-task", "blastp-fast", "-outfmt", "18"], 1, "1"),
        (&["-word_size", "4", "-outfmt", "10"], 4, TWO_BLOCKS),
    ];
    for (extra, threads, letters) in cases {
        let (reach, run) =
            same_as_one_table(&query, &subject, threads, &[BATCH], Some(letters), extra);
        assert_eq!(run.code, Some(0), "{extra:?}: {}", run.stderr);
        assert!(
            !run.stdout.is_empty(),
            "{extra:?}: the fixture pair has rows"
        );
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert!(
            sides.0,
            "{extra:?} threads {threads}: the blocks did not run: {sides:?}"
        );
    }
}

// `-db`: each block reads the database (`database::search_blocks`); the lists name the sequences by
// their ordinal ids, which the call of all queries ranks after the join (`ListedSubjects::x_list`).
// Exact blastp and blastp-fast, one thread and four, and a custom outfmt 7 whose `sallseqid` reads
// more of the definition lines than the block search takes.
#[test]
fn database_reports_with_blocks_of_the_preliminary_stage() {
    let (db, _db_dir) = subject_database("rows");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let (query, _dir) = queries("db", 4);
    let cases: [(&[&str], u32, &str); 3] = [
        (&["-outfmt", "10"], 4, TWO_BLOCKS),
        (&["-task", "blastp-fast", "-outfmt", "10"], 1, "1"),
        (
            &["-outfmt", "7 qseqid sallseqid evalue bitscore"],
            1,
            TWO_BLOCKS,
        ),
    ];
    for (extra, threads, letters) in cases {
        let (reach, run) =
            same_as_one_table(&query, &subject, threads, &[BATCH], Some(letters), extra);
        assert_eq!(run.code, Some(0), "{extra:?}: {}", run.stderr);
        assert!(
            !run.stdout.is_empty(),
            "{extra:?}: the fixture pair has rows"
        );
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert!(
            sides.0,
            "{extra:?} threads {threads}: the blocks did not run: {sides:?}"
        );
    }
}

// `-query_loc` (which the block search does not take): each query is its interval, and the batches
// count the whole records (`seq_range::QueryInput::batches`); a block's input is the records of its
// batches (`QueryInput::slice`). Outfmt 6, one thread and four.
#[test]
fn query_loc_takes_blocks_of_the_preliminary_stage() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let (query, _dir) = queries("loc", 4);
    for threads in [1, 4] {
        let extra = ["-query_loc", "10-300", "-outfmt", "6"];
        let (reach, run) =
            same_as_one_table(&query, &subject, threads, &[BATCH], Some("1"), &extra);
        assert_eq!(run.code, Some(0), "threads {threads}: {}", run.stderr);
        assert!(!run.stdout.is_empty(), "the fixture pair has rows");
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert!(
            sides.0,
            "threads {threads}: the blocks did not run: {sides:?}"
        );
    }
}

// `-word_size 4` at the default block of 50,000 letters: sixteen queries (8,709 letters) are one
// block of letters, but their table's words pass the most of a block (`x_aa_one_block_splits`), so
// the blocks close on their words (`x_word_blocks`, counted on the pool at four threads) and the
// preliminary stage runs in more than one block.
#[test]
fn word_size_4_blocks_close_on_their_words() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let (query, _dir) = queries("w4", 16);
    for threads in [1, 4] {
        let extra = ["-word_size", "4", "-outfmt", "10"];
        let (reach, run) = same_as_one_table(&query, &subject, threads, &[BATCH], None, &extra);
        assert_eq!(run.code, Some(0), "threads {threads}: {}", run.stderr);
        assert!(!run.stdout.is_empty(), "the fixture pair has rows");
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert!(
            sides.0,
            "threads {threads}: the blocks did not run: {sides:?}"
        );
    }
}

// A first batch whose queries are all invalid (two all-X records of 450 letters fill a batch of
// `BATCH`): NCBI does not search it, which is a property of the batch, and a block holds whole
// batches; its queries are reported as before (outfmt 10 writes nothing for them). And one block
// (a block size above the queries' letters): the search of all queries with one table, its
// reference side only.
#[test]
fn an_unsearched_first_batch_and_one_block() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let real = std::fs::read_to_string(queries("unsrch-real", 4).0).unwrap();
    let mut text = String::new();
    for index in 0..2 {
        text.push_str(&format!(">allx{index} all X\n{}\n", "X".repeat(450)));
    }
    text.push_str(&real);
    let (query, _dir) = queries_of("unsrch", &text);
    for threads in [1, 4] {
        let extra = ["-outfmt", "10"];
        let (reach, run) =
            same_as_one_table(&query, &subject, threads, &[BATCH], Some("1"), &extra);
        assert_eq!(run.code, Some(0), "threads {threads}: {}", run.stderr);
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert!(
            sides.0,
            "threads {threads}: the blocks did not run: {sides:?}"
        );
        let (reach, _) = same_as_one_table(
            &query,
            &subject,
            threads,
            &[BATCH],
            Some("100000000"),
            &extra,
        );
        assert_eq!(
            reach.sides("LOSAT_X_BPQBLOCK"),
            (false, true),
            "threads {threads}: one block"
        );
    }
}

// A database that the search cannot read: its errors come where one table gives them. Without a
// sequence file (`arch_sprot_nopsq`) the first block's read fails (`PreliminaryReadError`), and
// the search of all queries makes the error with its own batches and reports
// (`preliminary_search_error`); a database without sequences (`tests/seqdb/db/empty`) fails in the
// set-up of the search of all queries, before the first block. Outfmt 10 and 18, one thread and
// four.
#[test]
fn a_database_that_cannot_be_read_fails_as_with_one_table() {
    let (query, _dir) = queries("dberr", 4);
    for (db, extra) in [
        ("tests/db_search/arch_sprot_nopsq", ["-outfmt", "10"]),
        ("tests/db_search/arch_sprot_nopsq", ["-outfmt", "18"]),
        ("tests/seqdb/db/empty", ["-outfmt", "10"]),
    ] {
        for threads in [1, 4] {
            let subject = [OsStr::new("-db"), OsStr::new(db)];
            let (_, run) =
                same_as_one_table(&query, &subject, threads, &[BATCH], Some("1"), &extra);
            assert_eq!(
                run.code,
                Some(3),
                "{db} {extra:?} threads {threads}: {}",
                run.stderr
            );
            assert!(
                run.stderr
                    .ends_with("BLAST engine error: Database memory map file error\n"),
                "{db} {extra:?} threads {threads}: {}",
                run.stderr
            );
        }
    }
}

// NCBI's frozen bytes of a report of the search of all queries with several query batches
// (`db.of18.bp.batch.of18`: outfmt 18, `BATCH_SIZE=100`, `tests/db_search_manifest/blastp.tsv`):
// written in blocks of one batch each (`LOSAT_X_BPQBLOCK=1`), at one thread and four.
#[test]
fn the_frozen_ncbi_rows_of_several_batches_are_written_in_blocks() {
    let manifest = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search_manifest/blastp.tsv"),
    )
    .unwrap();
    let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
    let columns: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let column = |name: &str| columns.iter().position(|c| *c == name).unwrap();
    let (id, db, query, extra, format, rc) = (
        column("row_id"),
        column("db"),
        column("query"),
        column("extra_args"),
        column("outfmt"),
        column("rc"),
    );
    let frozen = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search");
    let mut rows = 0;
    for line in lines.filter(|line| line.starts_with("db.of18.bp.batch.")) {
        let cells: Vec<&str> = line.split('\t').collect();
        let row = cells[id];
        let expected_err =
            std::fs::read_to_string(frozen.join(format!("{row}.err"))).unwrap_or_default();
        let expected_out = std::fs::read(frozen.join(format!("{row}.out"))).unwrap();
        // `NAME=VALUE` tokens are variables of the run, the others arguments.
        let (variables, arguments): (Vec<&str>, Vec<&str>) = cells[extra]
            .split_whitespace()
            .partition(|token| token.contains('='));
        let variables: Vec<(&str, &str)> = variables
            .iter()
            .map(|token| token.split_once('=').unwrap())
            .collect();
        assert!(
            variables.iter().any(|(name, _)| *name == "BATCH_SIZE"),
            "{row}: several batches"
        );
        let mut args = arguments.clone();
        args.extend(["-outfmt", cells[format]]);
        for threads in [1, 4] {
            let subject = [OsStr::new("-db"), OsStr::new(cells[db])];
            let reach = reach_support::Reach::new("bpprelim-frozen");
            let mut vars = variables.clone();
            vars.extend([("LOSAT_X_BPQBLOCK", "1"), reach.var()]);
            let what = format!("{row} threads {threads}");
            let run = blastp_args(Path::new(cells[query]), &subject, threads, &vars, &args);
            assert_eq!(
                run.code.map(|code| code.to_string()).as_deref(),
                Some(cells[rc]),
                "{what}: {}",
                run.stderr
            );
            assert!(
                run.stdout == expected_out,
                "{what}: the bytes of NCBI differ"
            );
            assert_eq!(run.stderr, expected_err, "{what}: standard error differs");
            assert!(
                reach.sides("LOSAT_X_BPQBLOCK").0,
                "{what}: the blocks did not run"
            );
        }
        rows += 1;
    }
    assert_eq!(rows, 1);
}

/// `blastp` without the `-evalue 0.001` it adds (a frozen row's arguments are its own).
fn blastp_args(
    query: &Path,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    args: &[&str],
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text.starts_with("NCBI")
            || text.starts_with("BLAST")
        {
            command.env_remove(name);
        }
    }
    let log = worker_log_if_unset(&mut command, "blastp_prelim_blocks_worker");
    let home = TempPath::new(std::env::temp_dir().join(format!(
        "losat-bpprelim-home-{}-{threads}",
        std::process::id()
    )));
    std::fs::create_dir_all(&home).unwrap();
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("HOME", &*home)
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .args(subject)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("LOSAT_TIMING", "1")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr)
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[TIMING]"))
        .collect();
    Run {
        code: output.status.code(),
        stdout: output.stdout,
        stderr,
        _log: log,
    }
}
