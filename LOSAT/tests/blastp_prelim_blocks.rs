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
/// `TWO_BLOCKS` for `-query_loc 10-300` (291 letters a query, 1,164 in all, one block at 1,000):
/// {BDT62581.1 to BDT62583.1} (873 letters) and {BDT62584.1}.
const TWO_BLOCKS_LOC: &str = "500";

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
    subject_database_with(tag, &[])
}

/// `subject_database` made with the further makeblastdb arguments `extra` (`-parse_seqids`).
fn subject_database_with(tag: &str, extra: &[&str]) -> (PathBuf, TempPath) {
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
        .args(extra)
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
// outfmt 10 (CSV), 18 (the organism report) and 20 (CSV with a header line), and the ordinary
// traceback of `-comp_based_stats 0` (which takes each list's query from the list's place among
// all queries, so it does not check the join's rebase: `the_join_moves_each_list_to_its_query`
// does). One case per path (CI time).
#[test]
fn reports_of_the_search_of_all_queries_with_blocks_of_the_preliminary_stage() {
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let (query, _dir) = queries("subject", 4);
    let cases: [(&[&str], u32, &str); 7] = [
        (&["-comp_based_stats", "0", "-outfmt", "10"], 1, TWO_BLOCKS),
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

// The join of the blocks' hit lists moves each list and its HSPs to their query among all queries
// (`query_split::rebase_hit_list`): the list's index and each HSP's query index, context and query
// length. The queries with rows are BDT62581.1 (698 letters) and BDT62582.1 (268) in the first block
// of `TWO_BLOCKS` (the second block of `1`), and BDT62584.1 (944) in the last, and the ordinary
// traceback (`-comp_based_stats 0`) reads both:
// - without the Gumbel block (`OLD_FSC`) each HSP gets the e-value of its list's query's search
//   space (`list.query_index`, `blast_engine/traceback.rs`), and the queries' spaces differ: a list
//   left with its index in its block reads another query's;
// - `-qcov_hsp_perc 50` keeps the HSPs that cover half of their query's length (the HSP's
//   `query_length`): BDT62582.1's HSP of 221 letters and BDT62584.1's of 376 to 446 letters are
//   kept or dropped on the length of the query they name.
// Blocks of one batch and of several, one thread and four, `-subject` and `-db`.
#[test]
fn the_join_moves_each_list_to_its_query() {
    let fasta = fixture("MejoMJNV.faa");
    let fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let (db, _db_dir) = subject_database("join");
    let db = [OsStr::new("-db"), db.as_os_str()];
    let (query, _dir) = queries("join", 4);
    let old_fsc: &[(&str, &str)] = &[BATCH, ("OLD_FSC", "1")];
    let space: &[&str] = &["-comp_based_stats", "0", "-outfmt", "10"];
    let coverage: &[&str] = &[
        "-comp_based_stats",
        "0",
        "-qcov_hsp_perc",
        "50",
        "-outfmt",
        "10",
    ];
    let cases: [(&[&OsStr], u32, &[(&str, &str)], &str, &[&str]); 8] = [
        (&fasta, 1, old_fsc, TWO_BLOCKS, space),
        (&fasta, 4, old_fsc, "1", space),
        (&db, 1, old_fsc, "1", space),
        (&db, 4, old_fsc, TWO_BLOCKS, space),
        (&fasta, 1, &[BATCH], "1", coverage),
        (&fasta, 4, &[BATCH], TWO_BLOCKS, coverage),
        (&db, 1, &[BATCH], TWO_BLOCKS, coverage),
        (&db, 4, &[BATCH], "1", coverage),
    ];
    for (subject, threads, vars, letters, extra) in cases {
        let what = format!(
            "{:?} {extra:?} {vars:?} threads {threads} blocks {letters}",
            subject[0]
        );
        let (reach, run) = same_as_one_table(&query, subject, threads, vars, Some(letters), extra);
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(!run.stdout.is_empty(), "{what}: the fixture pair has rows");
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert_eq!(
            sides,
            (true, true),
            "{what}: the blocks of the preliminary stage did not run"
        );
    }
}

/// One case of `cases_with_blocks`: the subject arguments, the threads, the variables beside
/// `BATCH`, the blocks (`LOSAT_X_BPQBLOCK`) and the arguments.
type BlockCase<'a> = (
    &'a [&'a OsStr],
    u32,
    &'a [(&'a str, &'a str)],
    &'a str,
    &'a [&'a str],
);

/// Runs each case against one table (`same_as_one_table`, with `BATCH`), and asserts that it ends
/// with 0, writes rows and runs the preliminary stage of the search of all queries in blocks: both
/// sides of the switch, the block search's plan refused (its reference side) and the blocks of the
/// preliminary stage (its new side); a run of the block search records the new side only.
fn cases_with_blocks(query: &Path, cases: &[BlockCase<'_>]) {
    for &(subject, threads, vars, letters, extra) in cases {
        let mut all = vec![BATCH];
        all.extend_from_slice(vars);
        let (reach, run) = same_as_one_table(query, subject, threads, &all, Some(letters), extra);
        let what = format!(
            "{:?} {extra:?} {vars:?} threads {threads} blocks {letters}",
            subject[0]
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(!run.stdout.is_empty(), "{what}: the fixture pair has rows");
        let sides = reach.sides("LOSAT_X_BPQBLOCK");
        assert_eq!(
            sides,
            (true, true),
            "{what}: the blocks of the preliminary stage did not run"
        );
    }
}

// The hit-saving options that act inside a block's call: `-culling_limit` (each block's call has
// its own culling writer, `culling::PrelimCulling`, which equals one writer of all queries because
// the culling keeps one tree per query; the limit 1 drops most of the fixture pair's 27 rows) and
// `-max_target_seqs 1` (the preliminary hit list's size). One thread and four, `-subject` and
// `-db`, blocks of one batch and of several.
#[test]
fn culling_and_hit_list_sizes_with_blocks() {
    let fasta = fixture("MejoMJNV.faa");
    let fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let (db, _db_dir) = subject_database("cull");
    let db = [OsStr::new("-db"), db.as_os_str()];
    let (query, _dir) = queries("cull", 4);
    let culling: &[&str] = &["-culling_limit", "1", "-outfmt", "10"];
    let max_targets: &[&str] = &["-max_target_seqs", "1", "-outfmt", "10"];
    cases_with_blocks(
        &query,
        &[
            (&fasta, 1, &[], TWO_BLOCKS, culling),
            (&fasta, 4, &[], "1", culling),
            (&db, 1, &[], "1", culling),
            (&db, 4, &[], TWO_BLOCKS, culling),
            (&fasta, 4, &[], TWO_BLOCKS, max_targets),
            (&db, 1, &[], "1", max_targets),
        ],
    );
}

// Batches that the call of all queries treats as a whole, inside blocks other than the first:
// - a batch that NCBI splits into query chunks (`CHUNK_SIZE=800`: the batches {BDT62582.1,
//   BDT62583.1}, 897 letters, and {BDT62584.1}, 944): the block scans the batch, and the call of
//   all queries searches its chunks after the join (`search_split_query_batches`), whose lists
//   replace the block's;
// - a batch of invalid queries only (an all-X record of 600 letters after BDT62581.1), which NCBI
//   does not search: with blocks of one batch it is a block of its own, of invalid queries only;
// - an invalid-only first batch with a database whose sequence file cannot be read
//   (`arch_sprot_nopsq`): the unsearched batch is reported, then the read error of the first block
//   that searches, made by the call of all queries (`preliminary_search_error`), exit 3.
#[test]
fn split_and_unsearched_batches_with_blocks() {
    let fasta = fixture("MejoMJNV.faa");
    let fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let (db, _db_dir) = subject_database("split");
    let db = [OsStr::new("-db"), db.as_os_str()];
    let outfmt: &[&str] = &["-outfmt", "10"];
    let split: &[(&str, &str)] = &[("CHUNK_SIZE", "800")];
    let (query, _dir) = queries("split", 4);
    cases_with_blocks(
        &query,
        &[
            (&fasta, 1, split, TWO_BLOCKS, outfmt),
            (&fasta, 4, split, "1", outfmt),
            (&db, 1, split, "1", outfmt),
            (&db, 4, split, TWO_BLOCKS, outfmt),
        ],
    );
    let real = std::fs::read_to_string(&query).unwrap();
    let first_end = real[1..].find('>').unwrap() + 1;
    let all_x = format!(">allx all X\n{}\n", "X".repeat(600));
    let text = format!("{}{all_x}{}", &real[..first_end], &real[first_end..]);
    let (mid, _mid_dir) = queries_of("split-mid", &text);
    cases_with_blocks(
        &mid,
        &[
            (&fasta, 1, &[], "1", outfmt),
            (&fasta, 4, &[], TWO_BLOCKS, outfmt),
            (&db, 4, &[], "1", outfmt),
            (&db, 1, &[], TWO_BLOCKS, outfmt),
        ],
    );
    let (first, _first_dir) = queries_of("split-first", &format!("{all_x}{real}"));
    let nopsq = [
        OsStr::new("-db"),
        OsStr::new("tests/db_search/arch_sprot_nopsq"),
    ];
    for threads in [1, 4] {
        let (_, run) = same_as_one_table(&first, &nopsq, threads, &[BATCH], Some("1"), outfmt);
        assert_eq!(run.code, Some(3), "threads {threads}: {}", run.stderr);
        assert!(
            run.stderr
                .ends_with("BLAST engine error: Database memory map file error\n"),
            "threads {threads}: {}",
            run.stderr
        );
    }
}

/// `text` with the letters at offsets `from..to` of every sequence (each record's letters on one
/// line) in lower case.
fn lower_case(text: &str, from: usize, to: usize) -> String {
    text.lines()
        .map(|line| {
            if line.starts_with('>') || line.len() <= from {
                format!("{line}\n")
            } else {
                let to = to.min(line.len());
                format!(
                    "{}{}{}\n",
                    &line[..from],
                    line[from..to].to_ascii_lowercase(),
                    &line[to..]
                )
            }
        })
        .collect()
}

// `-lcase_masking`, lower-case letters in the queries (offsets 40 to 100) and in the subjects (20 to
// 60): the queries' masks are in the frames each block takes from the call of all queries, and the
// subjects' scan ranges are made once by that call for its blocks. Outfmt 10, and 7 and 0 with
// `-query_loc` (without it the block search writes them); one thread and four.
#[test]
fn lcase_masking_with_blocks() {
    let (query, _dir) = queries("lcase", 4);
    let text = std::fs::read_to_string(&query).unwrap();
    let joined: String = text
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| {
            let (title, letters) = record.split_once('\n').unwrap();
            format!(">{title}\n{}\n", letters.replace('\n', ""))
        })
        .collect();
    let (query, _lc_dir) = queries_of("lcase-q", &lower_case(&joined, 40, 100));
    let subjects = std::fs::read_to_string(fixture("MejoMJNV.faa")).unwrap();
    let subjects: String = subjects
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| {
            let (title, letters) = record.split_once('\n').unwrap();
            format!(">{title}\n{}\n", letters.replace('\n', ""))
        })
        .collect();
    let (subject, _subject_dir) = queries_of("lcase-s", &lower_case(&subjects, 20, 60));
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    cases_with_blocks(
        &query,
        &[
            (
                &subject,
                1,
                &[],
                TWO_BLOCKS,
                &["-lcase_masking", "-outfmt", "10"],
            ),
            (
                &subject,
                4,
                &[],
                "1",
                &["-lcase_masking", "-query_loc", "10-300", "-outfmt", "7"],
            ),
            (
                &subject,
                4,
                &[],
                TWO_BLOCKS_LOC,
                &["-lcase_masking", "-query_loc", "10-300", "-outfmt", "0"],
            ),
        ],
    );
}

// `-seqidlist` and `-negative_seqidlist` (a database made with `-parse_seqids`): the call of all
// queries ranks the listed subjects after the join (`ListedSubjects`); the list keeps subjects with
// rows of the first block's queries and of the last block's (BDT62772.1, BDT62777.1, BDT62780.1),
// the negative list the others. One thread and four.
#[test]
fn seqidlist_with_blocks() {
    let (db, db_dir) = subject_database_with("ids", &["-parse_seqids"]);
    let db = [OsStr::new("-db"), db.as_os_str()];
    let ids = db_dir.join("ids.txt");
    std::fs::write(&ids, "BDT62772.1\nBDT62777.1\nBDT62780.1\nBDT62750.1\n").unwrap();
    let ids = ids.to_str().unwrap();
    let (query, _dir) = queries("ids", 4);
    cases_with_blocks(
        &query,
        &[
            (
                &db,
                1,
                &[],
                TWO_BLOCKS,
                &["-seqidlist", ids, "-outfmt", "10"],
            ),
            (&db, 4, &[], "1", &["-seqidlist", ids, "-outfmt", "10"]),
            (
                &db,
                4,
                &[],
                TWO_BLOCKS,
                &["-negative_seqidlist", ids, "-outfmt", "10"],
            ),
        ],
    );
}

// The reports that send a `-query_loc` run to the search of all queries (the block search does not
// take `-query_loc`): outfmt 0, 5 (XML) and 11 (the archive); and the lookup tables of another
// matrix (PAM70) and of the word sizes 2 (exact blastp's direct table) and 5 (the compressed table
// of `-task blastp`), outfmt 10. One thread and four, `-subject`.
#[test]
fn query_loc_reports_and_other_tables_with_blocks() {
    let fasta = fixture("MejoMJNV.faa");
    let fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let (query, _dir) = queries("tables", 4);
    cases_with_blocks(
        &query,
        &[
            (
                &fasta,
                1,
                &[],
                "1",
                &["-query_loc", "10-300", "-outfmt", "0"],
            ),
            (
                &fasta,
                4,
                &[],
                TWO_BLOCKS_LOC,
                &["-query_loc", "10-300", "-outfmt", "5"],
            ),
            (
                &fasta,
                4,
                &[],
                "1",
                &["-query_loc", "10-300", "-outfmt", "11"],
            ),
            (
                &fasta,
                1,
                &[],
                TWO_BLOCKS,
                &["-matrix", "PAM70", "-outfmt", "10"],
            ),
            (&fasta, 4, &[], "1", &["-word_size", "2", "-outfmt", "10"]),
            (&fasta, 1, &[], "1", &["-word_size", "5", "-outfmt", "10"]),
        ],
    );
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
