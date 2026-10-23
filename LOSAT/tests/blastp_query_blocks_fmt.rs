// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! The blastp search in blocks of whole query batches (`LOSAT_X_BPQBLOCK`) for outfmt 0 and 7
//! with `LOSAT_X_BPQBLOCKFMT`, outfmt 8, 9 and 12 with `LOSAT_X_BPQBLOCKSA`, outfmt 5 with
//! `LOSAT_X_BPQBLOCKXML`, and outfmt 13 to 16 with `LOSAT_X_BPQBLOCKXML2`
//! (`src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switches are in the strict set (`src/utils/xswitch.rs`): on by default, off with `0` or
//! with `LOSAT_X_DEFAULT=0` and the variable not set. With them on, a search whose one report is
//! outfmt 0, 5, 6, 7, 8, 9 or 12 searches its queries block by block and writes each block's report
//! before it searches the next block: the outfmt 0 prolog and the XML header with the first block,
//! the epilog, the XML end and outfmt 7's count of the queries with the last. A run that the
//! block search could not report as the search of all queries does (an XML name NCBI cannot
//! write, a `-db` definition line that does not read) takes the search of all queries. These tests compare the standard output and the standard
//! error (without its `[TIMING]` lines) with the search of all queries (`LOSAT_X_BPQBLOCK=0`), and
//! count the `[TIMING] blastp_total` lines (one per search of a query set) to tell whether the
//! blocks ran. `BATCH_SIZE` (NCBI's query batch size) makes small batches, so that the fixture
//! pair has many batches and blocks. The shadow (`LOSAT_X_BPQBLOCKFMTSHADOW`) compares the block
//! search's bytes with the search of all queries inside the run.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-295
//! ```c++
//!         formatter.PrintProlog();
//!         ...
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!         ...
//!         formatter.PrintEpilog(opt);
//! ```

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const THREADS: [u32; 2] = [1, 4];
const BATCH: (&str, &str) = ("BATCH_SIZE", "3000");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");
const FMT_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCKFMT", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BPQBLOCKFMTSHADOW", "1");
const SA_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCKSA", "0");
const XML_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCKXML", "0");
const XML2_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCKXML2", "0");

/// The family switch of an outfmt of the Seq-annot and XML reports, set to `0`.
fn family_off(spec: &str) -> (&'static str, &'static str) {
    match spec {
        "5" => XML_OFF,
        "13" | "14" | "15" | "16" => XML2_OFF,
        "8" | "9" | "11" | "12" => SA_OFF,
        _ => FMT_OFF,
    }
}

/// One blastp run.
struct Run {
    success: bool,
    /// The exit status (`None` when a signal ended the process).
    code: Option<i32>,
    stdout: Vec<u8>,
    /// Standard error without the `[TIMING]` and `[X_SHADOW]` lines.
    stderr: String,
    /// The `[TIMING] blastp_total` lines: one per search of a query set.
    searches: usize,
    /// The `[X_SHADOW] BPQBLOCKFMT` lines.
    shadow: Vec<String>,
    /// The `[TIMING] x_bppfree` lines: one per search of a query set, also of a run that failed
    /// (the `blastp_total` line comes only from a run that ended well).
    frees: usize,
}

/// `blastp -query <query> <subject...> -evalue 0.001 -num_threads <threads> <extra...>` with the
/// given variables as its only `LOSAT_*` variables that choose a code path (`LOSAT_X_DEFAULT` and
/// `BATCH_SIZE` of the test's environment are dropped too). `LOSAT_BLASTX_WORKER_LOG` stays: a
/// binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs
/// that file name as soon as a search has more than one thread, and it selects nothing.
fn blastp(
    query: &Path,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&OsStr],
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "ARCHIVE_FORMAT"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
            || text.starts_with("NCBI")
            || text == "BLASTDB"
        {
            command.env_remove(name);
        }
    }
    // No `~/.ncbirc` names a taxonomy database or a BLASTDB (outfmt 15 and 16 look for one).
    let home = std::env::temp_dir().join(format!("losat-bpqblockfmt-home-{}", std::process::id()));
    std::fs::create_dir_all(&home).expect("empty home");
    command.env("HOME", &home).env("USERPROFILE", &home);
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpqblockfmt-worker.log"),
        );
    }
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .args(subject)
        .args(if extra.contains(&OsStr::new("-evalue")) {
            [].as_slice()
        } else {
            ["-evalue", "0.001"].as_slice()
        })
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .env("LOSAT_TIMING", "1")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let all = String::from_utf8_lossy(&output.stderr).into_owned();
    let searches = all
        .lines()
        .filter(|line| line.starts_with("[TIMING] blastp_total:"))
        .count();
    let frees = all
        .lines()
        .filter(|line| line.starts_with("[TIMING] x_bppfree:"))
        .count();
    let shadow = all
        .lines()
        .filter(|line| line.starts_with("[X_SHADOW] BPQBLOCKFMT"))
        .map(str::to_string)
        .collect();
    let stderr = all
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[TIMING]") && !line.starts_with("[X_SHADOW]"))
        .collect();
    Run {
        success: output.status.success(),
        code: output.status.code(),
        stdout: output.stdout,
        stderr,
        searches,
        shadow,
        frees,
    }
}

/// The subject of the first-hit row of the latest query in `rows` (outfmt 6 of the search of
/// all queries, in the order of the queries): a subject that no query of block 0 hits.
/// Returns the subject and the number (from 1) of its first-hit query in `query`'s file. The
/// refusal of a title that is made for that subject can only come after the earlier blocks'
/// bytes, so a check that is missing from the block plan changes the run's output.
fn subject_first_hit_late(query: &Path, rows: &str) -> (String, usize) {
    let text = std::fs::read_to_string(query).unwrap();
    let ids: Vec<&str> = text
        .lines()
        .filter_map(|line| line.strip_prefix('>'))
        .map(|line| line.split_whitespace().next().unwrap())
        .collect();
    let mut first: Vec<(&str, usize)> = Vec::new();
    for row in rows.lines() {
        let mut cells = row.split('\t');
        let (id, subject) = (cells.next().unwrap(), cells.next().unwrap());
        if first.iter().all(|(seen, _)| *seen != subject) {
            let number = ids.iter().position(|name| *name == id).unwrap() + 1;
            first.push((subject, number));
        }
    }
    let (subject, number) = first.into_iter().max_by_key(|(_, number)| *number).unwrap();
    // Block 0 holds the first five records at one thread (one batch of `BATCH`, 3,268 letters)
    // and about 26 at four (20,000 letters).
    assert!(number > 30, "the latest first hit is by record {number}");
    (subject.to_string(), number)
}

/// The (outfmt, threads) runs of one case of a test.
type Runs = &'static [(&'static str, u32)];

fn outfmt(spec: &str) -> [OsString; 2] {
    [OsString::from("-outfmt"), OsString::from(spec)]
}

fn os(args: &[OsString]) -> Vec<&OsStr> {
    args.iter().map(OsString::as_os_str).collect()
}

/// A protein database of the subject proteome, made with `losat makedb` in a directory of the
/// test (`tag` names it); returns the database name.
fn subject_database(tag: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-bpqblockfmt-db-{tag}-{}", std::process::id()));
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

/// Asserts that `run` wrote what `reference` wrote, and that it searched in blocks (more searches
/// than the reference) when `blocks` is set, else in one search.
fn assert_same(run: &Run, reference: &Run, blocks: bool, what: &str) {
    assert!(
        reference.success,
        "{what}: the reference failed: {}",
        reference.stderr
    );
    assert!(run.success, "{what}: {}", run.stderr);
    assert!(run.stdout == reference.stdout, "{what}: the report differs");
    assert_eq!(
        run.stderr, reference.stderr,
        "{what}: standard error differs"
    );
    assert_eq!(
        run.searches > reference.searches,
        blocks,
        "{what}: {} searches, the reference {}",
        run.searches,
        reference.searches
    );
    if !blocks {
        assert_eq!(run.searches, reference.searches, "{what}");
    }
}

// outfmt 0, 7 and a custom outfmt 7 (whose fields need the alignment) with `-subject` and `-db`:
// the blocks (one per query batch at one thread, of 20,000 letters at four) write the bytes of the
// search of all queries, and with `LOSAT_X_BPQBLOCKFMT=0` these formats take the search of all
// queries. One case per path (each search of the fixture pair costs CI seconds): outfmt 0 with
// each subject source, its serial writer at one thread (`-db`, whose reported sequences are
// numbered per block) and the writer on the pool at four (`LOSAT_X_PAIRWISEPAR`); outfmt 7 with
// `-db` (with `-subject`: the default test below); the custom fields; the switch off once (with
// outfmt 7: `blastp_query_blocks.rs`, which also has outfmt 6).
#[test]
fn outfmt_0_and_7_blocks_write_the_report_of_the_search_of_all_queries() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("reports");
    let with_fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let with_db = [OsStr::new("-db"), db.as_os_str()];
    // (subject, outfmt, threads, block letters, also with LOSAT_X_BPQBLOCKFMT=0)
    let cases = [
        (&with_fasta, "0", 4, "20000", true),
        (&with_db, "0", 1, "1", false),
        (&with_db, "7", 4, "20000", false),
        (
            &with_fasta,
            "7 qseqid sseqid pident evalue bitscore btop",
            1,
            "1",
            false,
        ),
    ];
    for (subject, spec, threads, letters, switch_off) in cases {
        let fmt = outfmt(spec);
        let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
        let reference = blastp(&query, subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        assert!(
            !reference.stdout.is_empty(),
            "{what}: the fixture pair has a report"
        );
        let run = blastp(
            &query,
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", letters)],
            &os(&fmt),
        );
        assert_same(
            &run,
            &reference,
            true,
            &format!("{what} blocks of {letters}"),
        );
        if switch_off {
            let run = blastp(
                &query,
                subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1"), FMT_OFF],
                &os(&fmt),
            );
            assert_same(
                &run,
                &reference,
                false,
                &format!("{what} LOSAT_X_BPQBLOCKFMT=0"),
            );
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}

// The default: with no variable outfmt 0 and 7 are searched in blocks of 50,000 query letters (the
// fixture's queries have 77,448: two blocks); `LOSAT_X_DEFAULT=0` searches all queries together,
// and leaves both variables on when they are set.
#[test]
fn the_default_writes_outfmt_0_and_7_in_blocks() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    for spec in ["0", "7"] {
        let fmt = outfmt(spec);
        let reference = blastp(&query, &subject, 4, &[BATCH, SWITCH_OFF], &os(&fmt));
        let run = blastp(&query, &subject, 4, &[BATCH], &os(&fmt));
        assert_same(&run, &reference, true, &format!("outfmt {spec} default"));
        assert_eq!(
            run.searches,
            reference.searches + 1,
            "outfmt {spec}: two blocks"
        );
        let run = blastp(&query, &subject, 4, &[BATCH, DEFAULT_OFF], &os(&fmt));
        assert_same(
            &run,
            &reference,
            false,
            &format!("outfmt {spec} LOSAT_X_DEFAULT=0"),
        );
        let run = blastp(
            &query,
            &subject,
            4,
            &[
                BATCH,
                DEFAULT_OFF,
                ("LOSAT_X_BPQBLOCK", "50000"),
                ("LOSAT_X_BPQBLOCKFMT", "1"),
            ],
            &os(&fmt),
        );
        assert_same(
            &run,
            &reference,
            true,
            &format!("outfmt {spec} LOSAT_X_DEFAULT=0 with the variables set"),
        );
    }
}

// outfmt 5, 8, 9, 11, 12, 15 and 16 with `-subject` and `-db`: the blocks (one per query batch at
// one thread, of 20,000 letters at four) write the bytes of the search of all queries (BER and the
// XML are compared as bytes), and with the family's switch at `0` these formats take the search of
// all queries. One case per path (each search of the fixture pair costs CI seconds): every format
// once, each family with both subject sources and both thread counts; outfmt 15 and 16 with both
// sources, as JSON joins the documents and a `-subject` run holds each document until the next
// (across a block's end); each family's switch off once.
#[test]
fn outfmt_5_8_9_12_15_16_blocks_write_the_report_of_the_search_of_all_queries() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("sareports");
    let with_fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let with_db = [OsStr::new("-db"), db.as_os_str()];
    // (subject, outfmt, threads, block letters, also with the family's switch off)
    let cases = [
        (&with_fasta, "5", 1, "1", false),
        (&with_db, "5", 4, "20000", true),
        (&with_fasta, "8", 4, "20000", true),
        (&with_db, "9", 1, "1", false),
        (&with_db, "11", 1, "1", false),
        (&with_fasta, "12", 1, "1", false),
        (&with_fasta, "15", 1, "1", false),
        (&with_db, "15", 4, "20000", false),
        (&with_fasta, "16", 4, "20000", false),
        (&with_db, "16", 1, "1", true),
    ];
    for (subject, spec, threads, letters, switch_off) in cases {
        let fmt = outfmt(spec);
        let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
        let reference = blastp(&query, subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        assert!(
            !reference.stdout.is_empty(),
            "{what}: the fixture pair has a report"
        );
        let run = blastp(
            &query,
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", letters)],
            &os(&fmt),
        );
        assert_same(
            &run,
            &reference,
            true,
            &format!("{what} blocks of {letters}"),
        );
        if switch_off {
            let run = blastp(
                &query,
                subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1"), family_off(spec)],
                &os(&fmt),
            );
            assert_same(&run, &reference, false, &format!("{what} family off"));
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
}

// The default: with no variable outfmt 5, 8, 9, 12, 15 and 16 are searched in blocks of 50,000
// query letters (two blocks); `LOSAT_X_DEFAULT=0` searches all queries together, and with the
// variables set (`LOSAT_X_BPQBLOCK` and the family's switch) the blocks are back. One format of
// each family (XML 5, Seq-annot 12, XML2 15); every format's blocks with the family's switch at
// its default: above.
#[test]
fn the_default_writes_outfmt_5_12_and_15_in_blocks() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    for spec in ["5", "12", "15"] {
        let fmt = outfmt(spec);
        let family = match spec {
            "5" => ("LOSAT_X_BPQBLOCKXML", "1"),
            "15" | "16" => ("LOSAT_X_BPQBLOCKXML2", "1"),
            _ => ("LOSAT_X_BPQBLOCKSA", "1"),
        };
        let reference = blastp(&query, &subject, 4, &[BATCH, SWITCH_OFF], &os(&fmt));
        let run = blastp(&query, &subject, 4, &[BATCH], &os(&fmt));
        assert_same(&run, &reference, true, &format!("outfmt {spec} default"));
        assert_eq!(
            run.searches,
            reference.searches + 1,
            "outfmt {spec}: two blocks"
        );
        let run = blastp(&query, &subject, 4, &[BATCH, DEFAULT_OFF], &os(&fmt));
        assert_same(
            &run,
            &reference,
            false,
            &format!("outfmt {spec} LOSAT_X_DEFAULT=0"),
        );
        let run = blastp(
            &query,
            &subject,
            4,
            &[BATCH, DEFAULT_OFF, ("LOSAT_X_BPQBLOCK", "50000"), family],
            &os(&fmt),
        );
        assert_same(
            &run,
            &reference,
            true,
            &format!("outfmt {spec} LOSAT_X_DEFAULT=0 with the variables set"),
        );
    }
}

/// A file of the test's directory `dir` holding `text`.
fn write_file(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

// Queries with warnings (O read as X, a record without letters) in several blocks: each query's
// warnings come before its report, in the same bytes, in every kind of writer: outfmt 0's serial
// writer (one thread) and the one on the pool (four), the tabular writer (7; 6 is the same
// function), the XML, JSON and XML2 messages (5, 15, 16), and a Seq-annot format (8), which writes
// none.
#[test]
fn query_warnings_keep_their_place_before_each_report() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-warn-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let queries = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    // The records again, with a record without letters and a copy of one query with O letters
    // after every 20th record.
    let mut text = String::new();
    let mut number = 0;
    for record in queries.split('>').filter(|record| !record.is_empty()) {
        text.push('>');
        text.push_str(record);
        number += 1;
        if number % 20 == 0 {
            text.push_str(&format!(
                ">empty{number} no letters\n>withO{number} O letters\n"
            ));
            let sequence: String = record.lines().skip(1).collect();
            text.push_str(&sequence.replacen('L', "O", 3));
            text.push('\n');
        }
    }
    assert!(number > 40);
    let query = write_file(&dir, "queries.faa", &text);
    let fasta = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    for (spec, threads) in [
        ("0", 1),
        ("0", 4),
        ("7", 1),
        ("5", 4),
        ("8", 1),
        ("11", 4),
        ("15", 4),
        ("16", 1),
    ] {
        let fmt = outfmt(spec);
        let what = format!("outfmt {spec} threads {threads}");
        let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        // Outfmt 0 and 7 write a query's warnings to standard error, the XML report in the
        // Iteration_message of the query; the Seq-annot format (8) and the archive (11) write none
        // (as NCBI), so there the equality below only shows that the blocks add none.
        match spec {
            "5" | "15" | "16" => {
                let report = String::from_utf8_lossy(&reference.stdout);
                assert!(
                    report.contains("O characters replaced")
                        && report.contains("Sequence contains no data"),
                    "{what}: the warnings in the Iteration_message"
                );
                assert!(reference.stderr.is_empty(), "{what}: {}", reference.stderr);
            }
            "8" | "11" => {
                assert!(reference.stderr.is_empty(), "{what}: {}", reference.stderr);
            }
            _ => assert!(
                reference.stderr.contains("O characters replaced")
                    && reference.stderr.contains("Sequence contains no data"),
                "{what}: the queries have warnings: {}",
                reference.stderr
            ),
        }
        let run = blastp(
            &query,
            &subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        assert_same(&run, &reference, true, &what);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// A subject title that LOSAT refuses where the outfmt 0 report shows the subject (a definition line
// of punctuation that NCBI's x_CleanAndCompress reads past the end of), on a subject with hits: the
// search of all queries refuses it before it writes anything, so with outfmt 0 the blocks are not
// searched and the run fails as without them; the tabular formats show the subject's id only, and
// are searched in blocks.
#[test]
fn a_refused_outfmt_0_title_leaves_the_search_of_all_queries() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-title-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let rows = blastp(
        &query,
        &[OsStr::new("-subject"), fasta.as_os_str()],
        1,
        &[BATCH, SWITCH_OFF],
        &os(&outfmt("6")),
    );
    let first_row = String::from_utf8(rows.stdout).unwrap();
    let hit = first_row.split('\t').nth(1).unwrap().to_string();
    let subjects = std::fs::read_to_string(&fasta).unwrap();
    let defline = subjects
        .lines()
        .find(|line| line.starts_with(&format!(">{hit} ")))
        .unwrap()
        .to_string();
    let fasta = write_file(
        &dir,
        "subjects.faa",
        &subjects.replacen(&defline, ">, ,", 1),
    );
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    for (spec, blocks) in [("0", false), ("7", true)] {
        let fmt = outfmt(spec);
        let reference = blastp(&query, &subject, 4, &[BATCH, SWITCH_OFF], &os(&fmt));
        let run = blastp(
            &query,
            &subject,
            4,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        if blocks {
            assert_same(&run, &reference, true, &format!("outfmt {spec}"));
            continue;
        }
        assert!(
            !reference.success,
            "outfmt {spec}: the reference refuses the title"
        );
        assert!(
            reference
                .stderr
                .contains("x_CleanAndCompress reads past the end"),
            "outfmt {spec}: {}",
            reference.stderr
        );
        assert!(!run.success, "outfmt {spec}");
        assert!(
            run.stdout == reference.stdout,
            "outfmt {spec}: the report differs"
        );
        assert_eq!(run.stderr, reference.stderr, "outfmt {spec}");
        assert_eq!(
            run.searches, reference.searches,
            "outfmt {spec}: not in blocks"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// The report written to a `-out` file, and the shadow: with `LOSAT_X_BPQBLOCKFMTSHADOW` each block
// search compares its output and its warnings with the search of all queries written into memory
// and says so in one line (it panics on a difference). The shadow implies the family's switch,
// which is `0` here (outfmt 6 has none of its own; `LOSAT_X_BPQBLOCKFMT=0` does not stop its
// blocks). One case per family (each search of the fixture pair costs CI seconds), with both
// subject sources and both thread counts; outfmt 6 has the shadow only here (its `-out` file:
// `blastp_query_blocks.rs`).
#[test]
fn the_out_file_and_the_shadow() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-out-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("shadow");
    let with_fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let with_db = [OsStr::new("-db"), db.as_os_str()];
    // (subject, outfmt, threads, also to a `-out` file)
    let cases = [
        (&with_fasta, "0", 1, true),
        (&with_db, "6", 4, false),
        (&with_db, "9", 4, true),
        (&with_fasta, "5", 1, true),
        (&with_fasta, "15", 4, true),
    ];
    for (subject, spec, threads, out) in cases {
        let fmt = outfmt(spec);
        let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
        let reference = blastp(&query, subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        if out {
            let path = dir.join(format!("report{spec}.txt"));
            let mut extra = fmt.to_vec();
            extra.extend([OsString::from("-out"), path.clone().into_os_string()]);
            let run = blastp(
                &query,
                subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
                &os(&extra),
            );
            assert!(
                run.success && run.stdout.is_empty(),
                "{what}: {}",
                run.stderr
            );
            assert!(run.searches > reference.searches, "{what}: blocks");
            assert!(
                std::fs::read(&path).unwrap() == reference.stdout,
                "{what}: the -out file differs"
            );
        }
        let run = blastp(
            &query,
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1"), family_off(spec), SHADOW],
            &os(&fmt),
        );
        assert!(run.success, "{what} shadow: {}", run.stderr);
        assert!(
            run.stdout == reference.stdout,
            "{what} shadow: the report differs"
        );
        assert_eq!(run.shadow.len(), 1, "{what} shadow: {:?}", run.shadow);
        assert!(
            run.shadow[0].starts_with("[X_SHADOW] BPQBLOCKFMT blocks=")
                && run.shadow[0].ends_with(" identical"),
            "{what} shadow: {}",
            run.shadow[0]
        );
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A protein of `length` random letters (a query or a subject that has no hit).
fn random_protein(state: &mut u64, length: usize) -> String {
    const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    (0..length)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            LETTERS[(*state % 20) as usize] as char
        })
        .collect()
}

/// Asserts that `run` ended as `reference` did (the same status, standard output and standard
/// error), in blocks (more searches) when `blocks` is set, else in one search.
fn assert_same_outcome(run: &Run, reference: &Run, blocks: bool, what: &str) {
    assert_eq!(run.success, reference.success, "{what}: the status differs");
    assert!(run.stdout == reference.stdout, "{what}: the output differs");
    assert_eq!(
        run.stderr, reference.stderr,
        "{what}: standard error differs"
    );
    if blocks {
        assert!(
            run.searches > reference.searches,
            "{what}: {} searches, the reference {}",
            run.searches,
            reference.searches
        );
    } else {
        assert_eq!(run.searches, reference.searches, "{what}: not in blocks");
    }
}

// Queries without hits at the edges of the blocks: the first query (the XML header is written from
// it), a run of 30 queries that fills whole batches (outfmt 8, 9 and 12 write nothing for them; the
// XML writes Iterations with "No hits found"), the last query (the XML end after a no-hit
// Iteration), and a query of X only (invalid). One case per path (CI time): outfmt 5, outfmt 8
// with `-db` and 12 with `-subject`, outfmt 15 and 16 holding the documents of a `-subject` run
// across the blocks of no-hit queries, outfmt 16 with `-db`; one and four threads.
#[test]
fn queries_without_hits_at_the_edges_of_blocks() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-nohit-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let real = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let real: Vec<&str> = real
        .split('>')
        .filter(|record| !record.is_empty())
        .collect();
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut text = String::new();
    let mut none = |text: &mut String, name: &str| {
        text.push_str(&format!(
            ">{name} no hits\n{}\n",
            random_protein(&mut state, 300)
        ));
    };
    none(&mut text, "first");
    for record in &real[..20] {
        text.push('>');
        text.push_str(record);
    }
    for index in 0..30 {
        none(&mut text, &format!("run{index}"));
    }
    for record in &real[20..40] {
        text.push('>');
        text.push_str(record);
    }
    text.push_str(">allx all X\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    none(&mut text, "last");
    let query = write_file(&dir, "queries.faa", &text);
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("nohit");
    let with_fasta = [OsStr::new("-subject"), fasta.as_os_str()];
    let with_db = [OsStr::new("-db"), db.as_os_str()];
    // (subject, outfmt, threads)
    let cases = [
        (&with_fasta, "5", 1),
        (&with_db, "8", 1),
        (&with_fasta, "12", 4),
        (&with_fasta, "15", 4),
        (&with_fasta, "16", 1),
        (&with_db, "16", 4),
    ];
    for (subject, spec, threads) in cases {
        let fmt = outfmt(spec);
        let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
        let reference = blastp(&query, subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        assert!(reference.success, "{what}: {}", reference.stderr);
        if matches!(spec, "5" | "15" | "16") {
            let xml = String::from_utf8_lossy(&reference.stdout);
            assert!(xml.contains("No hits found"), "{what}: no-hit Iterations");
            // The element of a hit in XML (5, 16) and JSON (15).
            assert!(
                xml.contains("<Hit>") || xml.contains("\"description\""),
                "{what}: hits too"
            );
        }
        let run = blastp(
            &query,
            subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        assert_same(&run, &reference, true, &what);
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

// Refusals that the block search could meet in a later block, with a name outside US-ASCII:
// - outfmt 5, 15 and 16, a query title in the 60th record: the search of all queries refuses
//   before it writes a byte, so the blocks are not used and the run fails as without them (outfmt
//   8, 9 and 12 write no title and are searched in blocks);
// - outfmt 5, 15 and 16, a `-subject` title on a subject with hits: refused the same way, as is the
//   `>, ,` title; on a subject without hits both complete, and the blocks are not used either (the
//   check names every subject);
// - outfmt 8, 9 and 12 take the same subjects in blocks.
// One case per path (each search of the fixture pair costs CI seconds): each refusal with outfmt 5
// and with outfmt 15 or 16 (the XML families share the checks, whose thread counts are unit tests
// in `x_query_blocks.rs`), the subject without hits and the `>, ,` title with one XML format each,
// outfmt 9 with the query title and 8 with a subject title.
#[test]
fn xml_names_that_ncbi_cannot_write_leave_the_search_of_all_queries() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-names-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    // (a) the 60th query's title.
    let queries = std::fs::read_to_string(&query).unwrap();
    let mut records: Vec<String> = queries
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| format!(">{record}"))
        .collect();
    let (header, rest) = records[59].split_once('\n').unwrap();
    records[59] = format!(
        "{}caf\u{e9}\n{rest}",
        &header[..header.find(' ').unwrap() + 1]
    );
    let cafe_query = write_file(&dir, "queries_cafe.faa", &records.concat());
    for (spec, threads) in [("5", 1), ("16", 4)] {
        let fmt = outfmt(spec);
        let what = format!("a query title, outfmt {spec}, threads {threads}");
        let reference = blastp(
            &cafe_query,
            &subject,
            threads,
            &[BATCH, SWITCH_OFF],
            &os(&fmt),
        );
        assert!(!reference.success && reference.stdout.is_empty(), "{what}");
        let run = blastp(
            &cafe_query,
            &subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        assert_same_outcome(&run, &reference, false, &what);
    }
    let (spec, threads) = ("9", 4);
    let fmt = outfmt(spec);
    let what = format!("a query title, outfmt {spec}, threads {threads}");
    let reference = blastp(
        &cafe_query,
        &subject,
        threads,
        &[BATCH, SWITCH_OFF],
        &os(&fmt),
    );
    let run = blastp(
        &cafe_query,
        &subject,
        threads,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
        &os(&fmt),
    );
    assert_same(&run, &reference, true, &what);
    // (b), (c), (d): subjects with a title outside US-ASCII.
    let rows = blastp(&query, &subject, 1, &[BATCH, SWITCH_OFF], &os(&outfmt("6")));
    let rows = String::from_utf8(rows.stdout).unwrap();
    let (hit, _) = subject_first_hit_late(&query, &rows);
    let subjects = std::fs::read_to_string(&fasta).unwrap();
    let defline = subjects
        .lines()
        .find(|line| line.starts_with(&format!(">{hit} ")))
        .unwrap()
        .to_string();
    let mut state = 0x1234_5678_9abc_def1u64;
    let with_hit = write_file(
        &dir,
        "subjects_hit.faa",
        &subjects.replacen(&defline, &format!(">{hit} caf\u{e9}"), 1),
    );
    let without_hit = write_file(
        &dir,
        "subjects_nohit.faa",
        &format!(
            "{subjects}>nohit caf\u{e9}\n{}\n",
            random_protein(&mut state, 200)
        ),
    );
    let punctuation = write_file(
        &dir,
        "subjects_punct.faa",
        &subjects.replacen(&defline, ">, ,", 1),
    );
    // (subject file, name, whether the XML report completes, its (outfmt, threads), the
    // Seq-annot (outfmt, threads))
    let cases: [(&PathBuf, &str, bool, Runs, Runs); 3] = [
        (
            &with_hit,
            "a title on a subject with hits",
            false,
            &[("5", 4), ("15", 1)],
            &[("8", 1)],
        ),
        (
            &without_hit,
            "a title on a subject without hits",
            true,
            &[("16", 1)],
            &[],
        ),
        (&punctuation, "the >, , title", false, &[("15", 4)], &[]),
    ];
    for (file, name, xml_succeeds, xml, seqalign) in cases {
        let subject = [OsStr::new("-subject"), file.as_os_str()];
        for &(spec, threads) in xml {
            let fmt = outfmt(spec);
            let what = format!("{name}, outfmt {spec}, threads {threads}");
            let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
            assert_eq!(
                reference.success, xml_succeeds,
                "{what}: {}",
                reference.stderr
            );
            if !xml_succeeds {
                assert!(
                    reference.stdout.is_empty(),
                    "{what}: no byte before the refusal"
                );
            }
            let run = blastp(
                &query,
                &subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
                &os(&fmt),
            );
            assert_same_outcome(&run, &reference, false, &what);
        }
        for &(spec, threads) in seqalign {
            let fmt = outfmt(spec);
            let what = format!("{name}, outfmt {spec}, threads {threads}");
            let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
            let run = blastp(
                &query,
                &subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
                &os(&fmt),
            );
            assert_same(&run, &reference, true, &what);
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// `-db`: one stored title with a Latin-1 byte, on a subject with hits (refused by outfmt 5, 8, 15
// and 16) and on one without (all complete): the block search is not used, and the run ends as the
// search of all queries does. One case per path (each search of the fixture pair costs CI
// seconds): with hits outfmt 5 and 8 (one check, `x_titles_read_allow_blocks`) and 15 (the titles
// with the taxids, `x_titles_taxids_allow_blocks`), without hits 8 and 16; both checks are unit
// tests at one and four threads (`database.rs`).
#[test]
fn a_database_title_outside_the_printable_range_leaves_the_search_of_all_queries() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-dbt-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let rows = blastp(
        &query,
        &[OsStr::new("-subject"), fasta.as_os_str()],
        1,
        &[BATCH, SWITCH_OFF],
        &os(&outfmt("6")),
    );
    let rows = String::from_utf8(rows.stdout).unwrap();
    let (hit, _) = subject_first_hit_late(&query, &rows);
    // `makedb` stores a byte outside ' ' ..= '~' as '#', so the title is made with four letters
    // that the test then replaces by a Latin-1 one in the header file (the same length).
    let subjects = std::fs::read_to_string(&fasta).unwrap();
    let defline = subjects
        .lines()
        .find(|line| line.starts_with(&format!(">{hit} ")))
        .unwrap()
        .to_string();
    let mut state = 0x0fed_cba9_8765_4321u64;
    let with_hit = subjects.replacen(&defline, &format!(">{hit} cafX"), 1);
    let without_hit = format!(
        "{subjects}>nohit cafY\n{}\n",
        random_protein(&mut state, 200)
    );
    // (database text, marker, name, whether the search of all queries completes, its (outfmt,
    // threads))
    let cases: [(&String, &str, &str, bool, Runs); 2] = [
        (
            &with_hit,
            "cafX",
            "a title on a subject with hits",
            false,
            &[("5", 1), ("8", 4), ("15", 1)],
        ),
        (
            &without_hit,
            "cafY",
            "a title on a subject without hits",
            true,
            &[("8", 1), ("16", 4)],
        ),
    ];
    for (text, marker, name, completes, runs) in cases {
        let fasta = dir.join(format!("{}.faa", name.replace(' ', "_")));
        std::fs::write(&fasta, text).unwrap();
        let db = dir.join(name.replace(' ', "_"));
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("makedb")
            .arg("-in")
            .arg(&fasta)
            .args(["-dbtype", "prot"])
            .arg("-out")
            .arg(&db)
            .output()
            .expect("run LOSAT makedb");
        assert!(
            output.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let header = db.with_extension("phr");
        let stored = std::fs::read(&header).unwrap();
        let at = stored
            .windows(4)
            .position(|window| window == marker.as_bytes())
            .expect("the title is stored once");
        let mut latin = stored.clone();
        latin[at + 3] = 0xe9;
        std::fs::write(&header, latin).unwrap();
        let subject = [OsStr::new("-db"), db.as_os_str()];
        for &(spec, threads) in runs {
            let fmt = outfmt(spec);
            let what = format!("{name}, outfmt {spec}, threads {threads}");
            let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
            assert_eq!(reference.success, completes, "{what}: {}", reference.stderr);
            let run = blastp(
                &query,
                &subject,
                threads,
                &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
                &os(&fmt),
            );
            assert_same_outcome(&run, &reference, false, &what);
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// outfmt 5 with a `-db` path in a directory whose name is outside US-ASCII: the database name is
// written as `BlastOutput_db`, which NCBI cannot write, so the search of all queries refuses before
// its first byte; the block search is not used and the run ends the same way (outfmt 16 too; outfmt
// 8 writes no name and is searched in blocks). One case per path (CI time); the check's thread
// counts are unit tests in `x_query_blocks.rs`.
#[cfg(unix)]
#[test]
fn a_database_name_outside_us_ascii_leaves_the_search_of_all_queries() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-dbn-{}", std::process::id()));
    let named = dir.join("d\u{e9}");
    std::fs::create_dir_all(&named).unwrap();
    let db = named.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let query = fixture("MeenMJNV.faa");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    // The same database under an ASCII name takes the blocks with outfmt 5, so the plan reaches
    // the name check (and the run leaves the blocks there, not before); outfmt 16's: the `-db`
    // cases of `outfmt_5_8_9_12_15_16_blocks_write_the_report_of_the_search_of_all_queries`.
    let plain = subject_database("ascii");
    let plain_subject = [OsStr::new("-db"), plain.as_os_str()];
    let (spec, threads) = ("5", 1);
    let fmt = outfmt(spec);
    let what = format!("an ASCII database name, outfmt {spec}, threads {threads}");
    let reference = blastp(
        &query,
        &plain_subject,
        threads,
        &[BATCH, SWITCH_OFF],
        &os(&fmt),
    );
    let run = blastp(
        &query,
        &plain_subject,
        threads,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
        &os(&fmt),
    );
    assert_same(&run, &reference, true, &what);
    std::fs::remove_dir_all(plain.parent().unwrap()).unwrap();
    for (spec, threads) in [("5", 4), ("16", 1), ("8", 4)] {
        let fmt = outfmt(spec);
        let what = format!("outfmt {spec} threads {threads}");
        let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        let run = blastp(
            &query,
            &subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        if spec != "8" {
            assert!(!reference.success && reference.stdout.is_empty(), "{what}");
            assert_same_outcome(&run, &reference, false, &what);
        } else {
            assert_same(&run, &reference, true, &what);
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// The eight NCBI-frozen rows `db.sa.bpqb.{5,8,9,12,15,16}`, `db.sa.bpqb.taxdb.16` and
// `db.sa.bpqb.taxdb.lastoff.15` (`tests/db_search_manifest.tsv`; the last two have a taxonomy
// database on `BLASTDB`): their 30
// queries have 11,940 letters, which is one block at the default `LOSAT_X_BPQBLOCK` (the fixture
// script drops `LOSAT_*`, so CI's db_search check never leaves the search of all queries for
// them). With `LOSAT_X_BPQBLOCK=1` the block search runs two searches (one per NCBI batch) and
// must write the bytes that NCBI froze, standard error included, at one and four threads.
#[test]
fn the_frozen_ncbi_rows_of_two_query_batches_are_written_by_two_blocks() {
    let manifest = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search_manifest.tsv"),
    )
    .unwrap();
    let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
    let columns: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let column = |name: &str| columns.iter().position(|c| *c == name).unwrap();
    let (id, db, query, extra, format) = (
        column("row_id"),
        column("db"),
        column("query"),
        column("extra_args"),
        column("outfmt"),
    );
    let frozen = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search");
    let mut rows = 0;
    for line in lines.filter(|line| line.starts_with("db.sa.bpqb.")) {
        let cells: Vec<&str> = line.split('\t').collect();
        // Outfmt 13 and 14 write side files: `the_frozen_ncbi_rows_of_outfmt_13_and_14_...` below.
        if matches!(cells[format], "13" | "14") {
            continue;
        }
        let row = cells[id];
        let mut expected_err = String::new();
        if let Ok(text) = std::fs::read_to_string(frozen.join(format!("{row}.err"))) {
            expected_err = text;
        }
        let expected_out = std::fs::read(frozen.join(format!("{row}.out"))).unwrap();
        for threads in THREADS {
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
            // No `~/.ncbirc` names a taxonomy database or a BLASTDB.
            let home =
                std::env::temp_dir().join(format!("losat-bpqblockfmt-home-{}", std::process::id()));
            std::fs::create_dir_all(&home).unwrap();
            command.env("HOME", &home).env("USERPROFILE", &home);
            if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
                command.env(
                    "LOSAT_BLASTX_WORKER_LOG",
                    std::env::temp_dir().join("losat-bpqblockfmt-worker.log"),
                );
            }
            // A `BLASTDB=DIR` token is the environment variable of the run (a `taxdb.bti`'s
            // directory, relative to LOSAT/, where the run starts), not an argument.
            let (variables, arguments): (Vec<&str>, Vec<&str>) = cells[extra]
                .split_whitespace()
                .partition(|token| token.starts_with("BLASTDB="));
            for variable in variables {
                command.env("BLASTDB", variable.trim_start_matches("BLASTDB="));
            }
            let output = command
                .current_dir(env!("CARGO_MANIFEST_DIR"))
                .arg("blastp")
                .args(["-db", cells[db], "-query", cells[query]])
                .args(arguments)
                .args(["-outfmt", cells[format]])
                .args(["-num_threads", &threads.to_string()])
                .env("LOSAT_TIMING", "1")
                .env("LOSAT_X_BPQBLOCK", "1")
                .output()
                .expect("run LOSAT");
            let all = String::from_utf8_lossy(&output.stderr).into_owned();
            let searches = all
                .lines()
                .filter(|line| line.starts_with("[TIMING] blastp_total:"))
                .count();
            let stderr: String = all
                .split_inclusive('\n')
                .filter(|line| !line.starts_with("[TIMING]"))
                .collect();
            let what = format!("{row} threads {threads}");
            assert!(output.status.success(), "{what}: {stderr}");
            assert!(searches >= 2, "{what}: {searches} searches, not in blocks");
            assert!(
                output.stdout == expected_out,
                "{what}: the bytes of NCBI differ"
            );
            assert_eq!(stderr, expected_err, "{what}: standard error differs");
        }
        rows += 1;
    }
    assert_eq!(rows, 8);
}

// `-db` with a stored title of punctuation (`, ,`) on a subject with hits: outfmt 6 cannot make the
// tabular id NCBI writes from it and refuses, so it leaves the blocks (`x_titles_allow_blocks`);
// outfmt 8, 9, 12 and 5 write no id, so they run in blocks (`x_titles_read_allow_blocks`) and write
// the bytes of the search of all queries, as do 15 and 16 (`x_titles_taxids_allow_blocks`). One
// case per path (CI time); the checks' thread counts are unit tests in `database.rs`.
#[test]
fn a_punctuation_database_title_stops_outfmt_6_blocks_but_not_the_others() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-punct-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let rows = blastp(
        &query,
        &[OsStr::new("-subject"), fasta.as_os_str()],
        1,
        &[BATCH, SWITCH_OFF],
        &os(&outfmt("6")),
    );
    let rows = String::from_utf8(rows.stdout).unwrap();
    let (hit, _) = subject_first_hit_late(&query, &rows);
    let subjects = std::fs::read_to_string(&fasta).unwrap();
    let defline = subjects
        .lines()
        .find(|line| line.starts_with(&format!(">{hit} ")))
        .unwrap()
        .to_string();
    // `makedb` stores the whole definition line as the title, so the title is `, ,` itself.
    let source = write_file(
        &dir,
        "subjects.faa",
        &subjects.replacen(&defline, ">, ,", 1),
    );
    let db = dir.join("punct");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(&source)
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let threads = 1;
    let fmt = outfmt("6");
    let what = format!("outfmt 6, threads {threads}");
    let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
    assert!(
        !reference.success && reference.stdout.is_empty(),
        "{what}: the title of punctuation is refused: {}",
        reference.stderr
    );
    let run = blastp(
        &query,
        &subject,
        threads,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
        &os(&fmt),
    );
    assert_same_outcome(&run, &reference, false, &what);
    for (spec, threads) in [("5", 4), ("9", 1), ("16", 4)] {
        let fmt = outfmt(spec);
        let what = format!("outfmt {spec}, threads {threads}");
        let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
        assert!(reference.success, "{what}: {}", reference.stderr);
        let run = blastp(
            &query,
            &subject,
            threads,
            &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
            &os(&fmt),
        );
        assert_same(&run, &reference, true, &what);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// Outfmt 15 and 16 with `-db` taxids (`tests/db_search/paja_taxid`: taxid 10090 in every definition
// line but BDT63547.1's 9606, a subject that only the second NCBI query batch of `bpqb_q.faa` hits)
// and a taxonomy database on `BLASTDB`. The taxonomy database's refusals come before the first byte
// in the search of all queries, so a name NCBI cannot write for 9606 (which a later block would
// refuse after the first block's bytes) and an unreadable `taxdb.btd` leave the block search; its
// messages on standard error (the open, the volume's lookup once per taxid, each `Hit`'s lookup) are
// written by the blocks as by the search of all queries, at the `Hit` that first needs them.
#[test]
fn taxonomy_names_for_outfmt_15_and_16_keep_the_bytes_and_the_refusals() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let query = manifest_dir.join("tests/db_search/bpqb_q.faa");
    let db = manifest_dir.join("tests/db_search/paja_taxid");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let extra = ["-max_target_seqs", "2", "-max_hsps", "1", "-evalue", "0.1"];
    let taxdb = |variant: &str| manifest_dir.join("tests/seqdb").join(variant);
    // The conservative fallback: the same queries without BDT62609.1, the one query that hits
    // BDT63547.1 (taxid 9606), so that no subject with hits carries 9606.
    let dir = std::env::temp_dir().join(format!("losat-bpqblockxml2-tax-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let text = std::fs::read_to_string(&query).unwrap();
    let mut kept = String::new();
    let mut skipping = false;
    for line in text.split_inclusive('\n') {
        if line.starts_with('>') {
            skipping = line.starts_with(">BDT62609.1 ");
        }
        if !skipping {
            kept.push_str(line);
        }
    }
    let query_without_9606 = write_file(&dir, "q_without_9606.faa", &kept);
    assert!(kept.len() < text.len() && kept.matches('>').count() == 29);
    // (taxonomy database variant, whether the search of all queries completes, whether it writes
    // to standard error, whether the queries are the reduced ones)
    let cases = [
        (None, true, false, false),
        (Some("taxdb/good"), true, false, false),
        // 9606's scientific name has UTF-8 bytes: refused, before the first byte.
        (Some("taxdb/high_bytes_utf8"), false, false, false),
        // An empty `taxdb.btd`: refused.
        (Some("taxdb/btd_empty"), false, false, false),
        // Messages only: the magic number (once, at the first `Hit`), and 10090's offset past the
        // end of the file (the volume's lookup once for the whole run, each `Hit`'s lookup).
        (Some("taxdb/bad_magic"), true, true, false),
        (Some("taxdb_xml2/last_10090_past_end"), true, true, false),
        // 9606's name cannot be written, but no subject with hits has it: the search of all
        // queries completes. The block search checks every taxid of the database, so it is not
        // taken (the check is conservative; a narrowing of it to the subjects with hits would
        // need the hits of the blocks to come first).
        (Some("taxdb/high_bytes_utf8"), true, false, true),
    ];
    for spec in ["15", "16"] {
        let mut fmt = outfmt(spec).to_vec();
        fmt.extend(extra.map(OsString::from));
        for threads in [1, 2, 4] {
            for (variant, completes, messages, reduced) in cases {
                let query = if reduced { &query_without_9606 } else { &query };
                // The block search is taken where the search of all queries completes and the
                // conservative check passes.
                let blocks = completes && !reduced;
                let what = format!(
                    "outfmt {spec}, taxdb {variant:?}, reduced {reduced}, threads {threads}"
                );
                let mut vars: Vec<(&str, String)> = Vec::new();
                if let Some(variant) = variant {
                    vars.push(("BLASTDB", taxdb(variant).to_str().unwrap().to_string()));
                }
                let run_with = |extra_vars: &[(&str, &str)]| {
                    let mut all: Vec<(&str, &str)> =
                        vars.iter().map(|(n, v)| (*n, v.as_str())).collect();
                    all.extend_from_slice(extra_vars);
                    blastp(query, &subject, threads, &all, &os(&fmt))
                };
                let reference = run_with(&[SWITCH_OFF]);
                assert_eq!(reference.success, completes, "{what}: {}", reference.stderr);
                // (LOSAT writes NCBI's "Examining 5 or more matches" line for `-max_target_seqs 2`.)
                assert_eq!(
                    reference.stderr.lines().count() > 1,
                    messages || !completes,
                    "{what}: {}",
                    reference.stderr
                );
                if !completes {
                    assert!(
                        reference.stdout.is_empty(),
                        "{what}: no byte before the refusal"
                    );
                }
                let run = run_with(&[("LOSAT_X_BPQBLOCK", "1")]);
                assert_same_outcome(&run, &reference, blocks, &what);
                if blocks {
                    let run = run_with(&[("LOSAT_X_BPQBLOCK", "1"), XML2_OFF]);
                    assert_same_outcome(&run, &reference, false, &format!("{what}, switch off"));
                    // The shadow's second search starts the taxonomy database over, so that it
                    // writes the open and lookup messages the blocks wrote.
                    let run = run_with(&[("LOSAT_X_BPQBLOCK", "1"), SHADOW]);
                    assert!(run.success, "{what} shadow: {}", run.stderr);
                    assert!(run.stdout == reference.stdout, "{what} shadow: stdout");
                    assert_eq!(run.shadow.len(), 1, "{what} shadow: {:?}", run.shadow);
                    assert!(
                        run.shadow[0].ends_with(" identical"),
                        "{what}: {}",
                        run.shadow[0]
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// Queries of a batch that NCBI did not search (every query in it is invalid): outfmt 15 and 16
// write them with no `query-masking` and a message that ends with "No hits found" (a property of
// the query's batch, and a block holds whole batches; outfmt 5 is compared too). Eight all-X records between
// real ones fill whole batches of 3,000 letters, so they sit in blocks of their own and in a later
// block than the first query; `-seg yes` gives the real queries a mask to list.
#[test]
fn queries_of_an_unsearched_batch_are_written_by_the_block_that_holds_them() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-unsrch-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let real = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let real: Vec<&str> = real
        .split('>')
        .filter(|record| !record.is_empty())
        .collect();
    let mut text = String::new();
    for record in &real[..6] {
        text.push('>');
        text.push_str(record);
    }
    for index in 0..8 {
        text.push_str(&format!(">allx{index} all X\n{}\n", "X".repeat(450)));
    }
    for record in &real[6..12] {
        text.push('>');
        text.push_str(record);
    }
    let query = write_file(&dir, "queries.faa", &text);
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("unsearched");
    for subject in [
        [OsStr::new("-subject"), fasta.as_os_str()],
        [OsStr::new("-db"), db.as_os_str()],
    ] {
        for spec in ["5", "15", "16"] {
            let mut fmt = outfmt(spec).to_vec();
            fmt.extend(["-seg", "yes"].map(OsString::from));
            for threads in THREADS {
                let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
                let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
                assert!(reference.success, "{what}: {}", reference.stderr);
                let report = String::from_utf8_lossy(&reference.stdout);
                // (Outfmt 5 lists no masks.)
                assert_eq!(
                    report.contains("query-masking") || report.contains("query_masking"),
                    spec != "5",
                    "{what}: the searched queries list a mask"
                );
                let run = blastp(
                    &query,
                    &subject,
                    threads,
                    &[BATCH, ("LOSAT_X_BPQBLOCK", "1")],
                    &os(&fmt),
                );
                assert_same(&run, &reference, true, &what);
            }
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

// The BLAST archive (outfmt 11) in blocks, on the fixture pair of the frozen `db.sa.bpqb.*` rows
// (`tests/db_search/paja`, 30 queries of 11,940 letters): `db.sa.arch.bp.paja.11` (two NCBI batches,
// two archives) and `db.sa.arch.bp.pajab.11` (`BATCH_SIZE=3000`, several). With `LOSAT_X_BPQBLOCK=80`
// every NCBI batch is a block of its own: the blocks' bytes equal the search of all queries
// (`LOSAT_X_BPQBLOCK=0`) and NCBI's frozen bytes (standard error included, the option warning of
// `-max_target_seqs 2` in the first archive), the shadow reports `identical`, `LOSAT_X_BPQBLOCKSA=0`
// keeps the search of all queries, and `ARCHIVE_FORMAT=bin` gives the same documents as ASN.1
// binary in blocks and out of them.
#[test]
fn the_archive_of_two_or_more_batches_is_written_by_blocks_as_by_the_search_of_all_queries() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let frozen = engine.join("tests/fixtures/db_search");
    // Relative to LOSAT/ (the helper runs there): the archive writes the `-db` name as typed.
    let query = PathBuf::from("tests/db_search/bpqb_q.faa");
    let database = PathBuf::from("tests/db_search/paja");
    let extra: Vec<OsString> = ["-evalue", "0.1", "-max_target_seqs", "2", "-max_hsps", "1"]
        .iter()
        .map(OsString::from)
        .collect();
    let mut fmt = outfmt("11").to_vec();
    fmt.extend(extra.iter().cloned());
    let subject = [OsStr::new("-db"), database.as_os_str()];
    let rows = [
        ("db.sa.arch.bp.paja.11", None),
        ("db.sa.arch.bp.pajab.11", Some(("BATCH_SIZE", "3000"))),
    ];
    for (row, batch) in rows {
        let want = std::fs::read(frozen.join(format!("{row}.out"))).unwrap();
        let want_err = std::fs::read_to_string(frozen.join(format!("{row}.err"))).unwrap();
        for threads in THREADS {
            let what = format!("{row} threads {threads}");
            let mut vars: Vec<(&str, &str)> = batch.into_iter().collect();
            // The -evalue of the helper is the fixture's: the row's own arguments follow.
            let reference_vars: Vec<(&str, &str)> =
                vars.iter().copied().chain([SWITCH_OFF]).collect();
            let reference = blastp(&query, &subject, threads, &reference_vars, &os(&fmt));
            assert!(reference.success, "{what}: {}", reference.stderr);
            assert!(
                reference.stdout == want,
                "{what}: the search of all queries is not NCBI's"
            );
            assert_eq!(reference.stderr, want_err, "{what}: standard error");
            vars.push(("LOSAT_X_BPQBLOCK", "80"));
            let blocks = blastp(&query, &subject, threads, &vars, &os(&fmt));
            assert_same(&blocks, &reference, true, &format!("{what} blocks"));
            assert!(
                blocks.stdout == want,
                "{what}: blocks differ from NCBI's bytes"
            );
            vars.push(SHADOW);
            let shadow = blastp(&query, &subject, threads, &vars, &os(&fmt));
            assert!(shadow.success && shadow.stdout == want, "{what} shadow");
            assert_eq!(shadow.shadow.len(), 1, "{what}: {:?}", shadow.shadow);
            assert!(
                shadow.shadow[0].ends_with(" identical"),
                "{what}: {}",
                shadow.shadow[0]
            );
            assert_eq!(shadow.stderr, want_err, "{what} shadow: standard error");
            vars.pop();
            vars.push(SA_OFF);
            let off = blastp(&query, &subject, threads, &vars, &os(&fmt));
            assert_same(&off, &reference, false, &format!("{what} family off"));
            // The binary archive: blocks and search of all queries agree.
            vars.pop();
            vars.push(("ARCHIVE_FORMAT", "bin"));
            let binary = blastp(&query, &subject, threads, &vars, &os(&fmt));
            vars.retain(|(name, _)| *name != "LOSAT_X_BPQBLOCK");
            let binary_reference = blastp(&query, &subject, threads, &vars, &os(&fmt));
            assert_same(&binary, &binary_reference, true, &format!("{what} binary"));
        }
    }
}

// ---- outfmt 13 and 14 (OF-BPQB3): side files and a master file ----

/// What a run with `-out` left in its directory: the file name and bytes of every regular file, a
/// directory as `<directory>` and a link as `-> target` (a link to `/dev/full` must not be read).
type Files = std::collections::BTreeMap<String, Vec<u8>>;

/// A run of outfmt `spec` (13 or 14) whose `-out` is `<dir>/rep.json` (13) or `<dir>/rep.xml` (14),
/// in a fresh directory `dir` that `prepare` may fill first; the side files are `dir/rep_<n>.json`
/// or `.xml` and the master is `dir/rep.json` / `dir/rep.xml`.
#[allow(clippy::too_many_arguments)]
fn files_run(
    query: &Path,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    spec: &str,
    dir: &Path,
    extra: &[&str],
    prepare: &dyn Fn(&Path),
) -> (Run, Files) {
    let _ = std::fs::remove_dir_all(dir);
    std::fs::create_dir_all(dir).unwrap();
    prepare(dir);
    let extension = if spec == "13" { "json" } else { "xml" };
    let mut args = outfmt(spec).to_vec();
    args.extend([
        OsString::from("-out"),
        dir.join(format!("rep.{extension}")).into_os_string(),
    ]);
    args.extend(extra.iter().map(OsString::from));
    let run = blastp(query, subject, threads, vars, &os(&args));
    let mut files = Files::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().to_string_lossy().into_owned();
        let kind = entry.path().symlink_metadata().unwrap().file_type();
        let bytes = if kind.is_symlink() {
            format!("-> {}", entry.path().read_link().unwrap().display()).into_bytes()
        } else if kind.is_dir() {
            b"<directory>".to_vec()
        } else {
            std::fs::read(entry.path()).unwrap()
        };
        files.insert(name, bytes);
    }
    (run, files)
}

/// Asserts that two runs left the same files, the same status and the same standard error, and
/// that `first` searched in blocks (more searches) when `blocks` is set, else in one search.
fn assert_same_files(
    first: &(Run, Files),
    reference: &(Run, Files),
    blocks: Option<bool>,
    what: &str,
) {
    assert_eq!(first.0.code, reference.0.code, "{what}: the status differs");
    assert_eq!(
        first.0.stderr, reference.0.stderr,
        "{what}: standard error differs"
    );
    let names = |files: &Files| files.keys().cloned().collect::<Vec<_>>();
    assert_eq!(
        names(&first.1),
        names(&reference.1),
        "{what}: the files differ"
    );
    for (name, bytes) in &first.1 {
        assert!(
            *bytes == reference.1[name],
            "{what}: {name} differs ({} bytes, the reference {})",
            bytes.len(),
            reference.1[name].len()
        );
    }
    // The searches of a query set: `blastp_total` lines, or `x_bppfree` lines when the run failed.
    let searches = |run: &Run| run.searches.max(run.frees);
    match blocks {
        Some(true) => assert!(
            searches(&first.0) > searches(&reference.0),
            "{what}: {} searches, the reference {}",
            searches(&first.0),
            searches(&reference.0)
        ),
        Some(false) => assert_eq!(
            searches(&first.0),
            searches(&reference.0),
            "{what}: not in blocks"
        ),
        None => {}
    }
}

/// `blocks=<n>` and `side_files=<n>` of the one shadow line of `run`, which must end in `identical`.
fn shadow_blocks_and_side_files(run: &Run, what: &str) -> (usize, usize) {
    assert_eq!(run.shadow.len(), 1, "{what}: {:?}", run.shadow);
    let line = &run.shadow[0];
    assert!(line.ends_with(" identical"), "{what}: {line}");
    let value = |key: &str| -> usize {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(key))
            .unwrap_or_else(|| panic!("{what}: no {key} in {line}"))
            .parse()
            .unwrap()
    };
    (value("blocks="), value("side_files="))
}

/// The reports of 13 and 14 are written serially after each block's search, so the thread count does
/// not reach the side-file code: the cases that look at the files alone run at one thread.
const ONE_THREAD: [u32; 1] = [1];
/// The side file of a report in the second block (queries 27-30) of `bpqb_q.faa` at `LOSAT_X_BPQBLOCK=80`.
const SECOND_BLOCK: [usize; 1] = [28];
const BLOCK_80: (&str, &str) = ("LOSAT_X_BPQBLOCK", "80");
const NO_PREPARE: &dyn Fn(&Path) = &|_| {};
/// The options of the frozen `db.sa.bpqb.*` rows.
const PAJA_OPTIONS: [&str; 6] = ["-evalue", "0.1", "-max_target_seqs", "2", "-max_hsps", "1"];

fn test_dir(tag: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-bpqblockfiles-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn paja() -> (PathBuf, PathBuf) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    (root.join("bpqb_q.faa"), root.join("paja"))
}

// Outfmt 13 and 14 with `-db`: the 30 queries of `bpqb_q.faa` are two NCBI batches, so with
// `LOSAT_X_BPQBLOCK=80` two blocks (queries 1-26, 27-30) write the master and the 30 side files
// that the search of all queries writes (`LOSAT_X_BPQBLOCKXML2=0`); the shadow run writes the blocks'
// files and says so (two blocks, 30 side files), with the side files compared with the reference's.
// The reports are written serially after each block's search, so the thread count does not reach
// the side-file code; both thread counts are run here once, for the files to show it.
#[test]
fn outfmt_13_and_14_blocks_write_the_side_files_and_master_of_the_search_of_all_queries() {
    let (query, db) = paja();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let dir = test_dir("db");
    for spec in ["13", "14"] {
        for threads in THREADS {
            let what = format!("outfmt {spec} threads {threads}");
            let run = |vars: &[(&str, &str)], at: &str| {
                files_run(
                    &query,
                    &subject,
                    threads,
                    vars,
                    spec,
                    &dir.join(at),
                    &PAJA_OPTIONS,
                    NO_PREPARE,
                )
            };
            let reference = run(&[BLOCK_80, XML2_OFF], "ref");
            assert!(reference.0.success, "{what}: {}", reference.0.stderr);
            assert_eq!(
                reference.1.len(),
                31,
                "{what}: the master and 30 side files"
            );
            // The default switch, with `LOSAT_X_DEFAULT=0` (the variable set: on).
            let shadow = run(&[BLOCK_80, SHADOW], "shadow");
            assert_same_files(&shadow, &reference, None, &format!("{what} shadow"));
            let (count, sides) = shadow_blocks_and_side_files(&shadow.0, &what);
            assert_eq!((count, sides), (2, 30), "{what}: blocks and side files");
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// Outfmt 14 with `-subject` (a FASTA subject file): the last report of a block is held until the
// next report arrives, which is in the next block, and is written as the next number. The shadow run
// writes the blocks' files and compares them with the reference's. Outfmt 13 holds a report through
// the same `Xml2Document`, `Xml2Carry` and `FsSideFiles` (the db test above and the cases below run
// both), and the search runs at one thread (the reports are written serially).
#[test]
fn outfmt_13_and_14_with_subject_hold_a_report_across_blocks() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), fasta.as_os_str()];
    let dir = test_dir("subject");
    for spec in ["14"] {
        for threads in ONE_THREAD {
            let what = format!("outfmt {spec} threads {threads}");
            let run = |vars: &[(&str, &str)], at: &str| {
                files_run(
                    &query,
                    &subject,
                    threads,
                    vars,
                    spec,
                    &dir.join(at),
                    &[],
                    NO_PREPARE,
                )
            };
            let reference = run(&[BATCH, ("LOSAT_X_BPQBLOCK", "1"), XML2_OFF], "ref");
            assert!(reference.0.success, "{what}: {}", reference.0.stderr);
            assert!(
                reference.1.len() > 10,
                "{what}: {} files",
                reference.1.len()
            );
            let shadow = run(&[BATCH, ("LOSAT_X_BPQBLOCK", "1"), SHADOW], "shadow");
            assert_same_files(&shadow, &reference, None, &format!("{what} shadow"));
            let (count, sides) = shadow_blocks_and_side_files(&shadow.0, &what);
            assert!(count >= 2, "{what}: {count} blocks");
            assert_eq!(sides, reference.1.len() - 1, "{what}: side files");
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// Outfmt 13 and 14 with taxids and a taxonomy database: the names in the side files and the
// messages on standard error are those of the search of all queries; a database that the search of
// all queries refuses (9606's name has UTF-8 bytes) leaves the blocks and fails before any file.
// Both formats are run (names render differently in JSON and XML), at one thread; the shadow run
// writes the blocks' files, and only a refused database has a separate blocks run, which no
// shadow run repeats.
#[test]
fn outfmt_13_and_14_with_taxids_and_a_taxonomy_database_keep_the_files_and_messages() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let query = manifest_dir.join("tests/db_search/bpqb_q.faa");
    let db = manifest_dir.join("tests/db_search/paja_taxid");
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let dir = test_dir("taxdb");
    // (taxonomy database variant, whether the search of all queries completes, whether it writes
    // messages to standard error)
    let cases = [
        (None, true, false),
        (Some("taxdb/good"), true, false),
        (Some("taxdb/bad_magic"), true, true),
        (Some("taxdb_xml2/last_10090_past_end"), true, true),
        (Some("taxdb/high_bytes_utf8"), false, false),
    ];
    for spec in ["13", "14"] {
        for threads in ONE_THREAD {
            for (variant, completes, messages) in cases {
                let what = format!("outfmt {spec}, taxdb {variant:?}, threads {threads}");
                let blastdb = variant
                    .map(|variant| manifest_dir.join("tests/seqdb").join(variant))
                    .map(|path| path.to_str().unwrap().to_string());
                let mut vars: Vec<(&str, &str)> = vec![BLOCK_80];
                if let Some(path) = &blastdb {
                    vars.push(("BLASTDB", path));
                }
                let run = |extra_vars: &[(&str, &str)], at: &str| {
                    let mut all = vars.clone();
                    all.extend_from_slice(extra_vars);
                    files_run(
                        &query,
                        &subject,
                        threads,
                        &all,
                        spec,
                        &dir.join(at),
                        &PAJA_OPTIONS,
                        NO_PREPARE,
                    )
                };
                let reference = run(&[XML2_OFF], "ref");
                assert_eq!(
                    reference.0.success, completes,
                    "{what}: {}",
                    reference.0.stderr
                );
                assert_eq!(
                    reference.0.stderr.lines().count() > 1,
                    messages || !completes,
                    "{what}: {}",
                    reference.0.stderr
                );
                if completes {
                    assert_eq!(reference.1.len(), 31, "{what}");
                    let shadow = run(&[SHADOW], "shadow");
                    assert_same_files(&shadow, &reference, None, &format!("{what} shadow"));
                    let (count, sides) = shadow_blocks_and_side_files(&shadow.0, &what);
                    assert_eq!((count, sides), (2, 30), "{what}");
                } else {
                    // The refused database leaves the blocks (the check is made before the first
                    // block), so that no side file is written.
                    let blocks = run(&[], "blocks");
                    assert_same_files(&blocks, &reference, Some(false), &what);
                    assert!(
                        reference.1.values().all(Vec::is_empty),
                        "{what}: the refused run wrote a file"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `bpqb_q.faa` with an O letter in queries `with_o` (numbers from 1), whose "O characters
/// replaced" warning goes to standard error before the query's report.
fn queries_with_warnings(dir: &Path, with_o: &[usize]) -> PathBuf {
    let text = std::fs::read_to_string(paja().0).unwrap();
    let mut out = String::new();
    let mut number = 0;
    for record in text.split('>').filter(|record| !record.is_empty()) {
        number += 1;
        let (head, body) = record.split_once('\n').unwrap();
        let body = if with_o.contains(&number) {
            body.replacen('L', "O", 2)
        } else {
            body.to_string()
        };
        out.push('>');
        out.push_str(head);
        out.push('\n');
        out.push_str(&body);
    }
    write_file(dir, "q_warn.faa", &out)
}

// A side file that cannot be opened (its name is a directory), in the first block, at the first and
// last query of the second block and at its middle: NCBI stops at that report with "Cannot open
// output file" (status 1), the side files before it stay and the master stays empty. The blocks
// stop at the same report with the same files and the same standard error, where queries after
// the failing one (and one in the first block) have warnings. Outfmt 13 only, at one thread, at the
// first block's first report and at a later block's: outfmt 14 opens its side files through the
// same `FsSideFiles`, and the reports are written serially.
#[test]
fn an_unopenable_side_file_stops_the_blocks_where_it_stops_the_search_of_all_queries() {
    let dir = test_dir("open");
    let query = queries_with_warnings(&dir, &[5, 29, 30]);
    let (_, db) = paja();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    for spec in ["13"] {
        let extension = if spec == "13" { "json" } else { "xml" };
        for threads in ONE_THREAD {
            for n in [1usize, 28] {
                let what = format!("outfmt {spec}, side file {n}, threads {threads}");
                let prepare = |at: &Path| {
                    std::fs::create_dir(at.join(format!("rep_{n}.{extension}"))).unwrap();
                };
                let run = |vars: &[(&str, &str)], at: &str| {
                    files_run(
                        &query,
                        &subject,
                        threads,
                        vars,
                        spec,
                        &dir.join(at),
                        &PAJA_OPTIONS,
                        &prepare,
                    )
                };
                let reference = run(&[BLOCK_80, XML2_OFF], "ref");
                assert_eq!(reference.0.code, Some(1), "{what}: {}", reference.0.stderr);
                assert!(
                    reference
                        .0
                        .stderr
                        .contains("Command line argument error: Cannot open output file"),
                    "{what}: {}",
                    reference.0.stderr
                );
                // The side files before n, the directory and an empty master.
                assert_eq!(reference.1.len(), n + 1, "{what}");
                assert!(
                    reference.1[&format!("rep.{extension}")].is_empty(),
                    "{what}: the master is not empty"
                );
                assert_eq!(
                    reference.1[&format!("rep_{n}.{extension}")],
                    b"<directory>",
                    "{what}"
                );
                let blocks = run(&[BLOCK_80], "blocks");
                // A failure in the first block (26 queries) is reached by its search alone.
                assert_same_files(&blocks, &reference, (n > 26).then_some(true), &what);
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// A side file that is a link to `/dev/full`, in the second block: the write fails
// after the file was opened. LOSAT stops the run in every program with "BLAST failed to write
// output" and status 6 (DECISIONS 35: NCBI ignores the error); the blocks do the same at the same
// report with the same files. Outfmt 14 only, at one thread, in the second block (the first block's
// stop is that of the test above, and outfmt 13 writes through the same `FsSideFiles`).
#[cfg(unix)]
#[test]
fn a_side_file_that_cannot_be_written_stops_the_blocks_where_it_stops_the_search_of_all_queries() {
    let dir = test_dir("full");
    let query = queries_with_warnings(&dir, &[5, 29, 30]);
    let (_, db) = paja();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    for spec in ["14"] {
        let extension = if spec == "13" { "json" } else { "xml" };
        for threads in ONE_THREAD {
            for n in SECOND_BLOCK {
                let what = format!("outfmt {spec}, side file {n}, threads {threads}");
                let prepare = |at: &Path| {
                    std::os::unix::fs::symlink(
                        "/dev/full",
                        at.join(format!("rep_{n}.{extension}")),
                    )
                    .unwrap();
                };
                let run = |vars: &[(&str, &str)], at: &str| {
                    files_run(
                        &query,
                        &subject,
                        threads,
                        vars,
                        spec,
                        &dir.join(at),
                        &PAJA_OPTIONS,
                        &prepare,
                    )
                };
                let reference = run(&[BLOCK_80, XML2_OFF], "ref");
                assert_eq!(reference.0.code, Some(6), "{what}: {}", reference.0.stderr);
                assert!(
                    reference.0.stderr.contains("BLAST failed to write output"),
                    "{what}: {}",
                    reference.0.stderr
                );
                assert_eq!(reference.1.len(), n + 1, "{what}");
                assert!(
                    reference.1[&format!("rep.{extension}")].is_empty(),
                    "{what}"
                );
                let blocks = run(&[BLOCK_80], "blocks");
                // A failure in the first block (26 queries) is reached by its search alone.
                assert_same_files(&blocks, &reference, (n > 26).then_some(true), &what);
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// Queries without a hit and an empty query file: the blocks write what the search of all queries
// writes (a side file per query, even one without hits; for an empty query file the plan is not
// taken and the master lists nothing). Both formats and both inputs, at one thread; the shadow run
// writes the blocks' files and compares them with the reference's.
#[test]
fn outfmt_13_and_14_without_hits_and_without_queries_are_unchanged() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let (_, db) = paja();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let dir = test_dir("nohit");
    for spec in ["13", "14"] {
        for name in ["nohit.faa", "empty.faa"] {
            let query = root.join(name);
            for threads in ONE_THREAD {
                let what = format!("outfmt {spec}, {name}, threads {threads}");
                let run = |vars: &[(&str, &str)], at: &str| {
                    files_run(
                        &query,
                        &subject,
                        threads,
                        vars,
                        spec,
                        &dir.join(at),
                        &PAJA_OPTIONS,
                        NO_PREPARE,
                    )
                };
                let reference = run(&[BLOCK_80, XML2_OFF], "ref");
                let shadow = run(&[BLOCK_80, SHADOW], "shadow");
                assert_same_files(&shadow, &reference, None, &format!("{what} shadow"));
                // One block or no query: the search of all queries, in one search, so that the
                // plan of blocks is not taken and the shadow has no line.
                assert!(shadow.0.shadow.is_empty(), "{what}: {:?}", shadow.0.shadow);
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// The frozen NCBI rows `db.sa.bpqb.{13,14,taxdb.13,taxdb.14}` (`tests/db_search_manifest.tsv`; the
// last two have a taxonomy database on `BLASTDB`): 30 queries are one block at the default
// `LOSAT_X_BPQBLOCK`, so with `LOSAT_X_BPQBLOCK=1` the two blocks (one per NCBI batch) must write the
// master and the 30 side files that NCBI froze, and its standard error, at one thread (the reports
// are written serially, so the thread count does not reach the files).
#[test]
fn the_frozen_ncbi_rows_of_outfmt_13_and_14_are_written_by_two_blocks() {
    let manifest = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search_manifest.tsv"),
    )
    .unwrap();
    let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
    let columns: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let column = |name: &str| columns.iter().position(|c| *c == name).unwrap();
    let (id, db, query, extra, format) = (
        column("row_id"),
        column("db"),
        column("query"),
        column("extra_args"),
        column("outfmt"),
    );
    let frozen = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search");
    let dir = test_dir("frozen");
    let mut rows = 0;
    for line in lines.filter(|line| line.starts_with("db.sa.bpqb.")) {
        let cells: Vec<&str> = line.split('\t').collect();
        if !matches!(cells[format], "13" | "14") {
            continue;
        }
        let (row, spec) = (cells[id], cells[format]);
        let extension = if spec == "13" { "json" } else { "xml" };
        let expected_err =
            std::fs::read_to_string(frozen.join(format!("{row}.err"))).unwrap_or_default();
        let (variables, arguments): (Vec<&str>, Vec<&str>) = cells[extra]
            .split_whitespace()
            .partition(|token| token.starts_with("BLASTDB="));
        // The side files are `<row>.out_<n>.<ext>` next to the master `<row>.out`.
        let mut expected: Vec<(String, Vec<u8>)> = vec![(
            format!("{row}.out"),
            std::fs::read(frozen.join(format!("{row}.out"))).unwrap(),
        )];
        for n in 1.. {
            let name = format!("{row}.out_{n}.{extension}");
            match std::fs::read(frozen.join(&name)) {
                Ok(bytes) => expected.push((name, bytes)),
                Err(_) => break,
            }
        }
        assert_eq!(expected.len(), 31, "{row}");
        for threads in ONE_THREAD {
            let what = format!("{row} threads {threads}");
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            let out = dir.join(format!("{row}.out"));
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
            let home =
                std::env::temp_dir().join(format!("losat-bpqblockfmt-home-{}", std::process::id()));
            std::fs::create_dir_all(&home).unwrap();
            command.env("HOME", &home).env("USERPROFILE", &home);
            if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
                command.env(
                    "LOSAT_BLASTX_WORKER_LOG",
                    std::env::temp_dir().join("losat-bpqblockfmt-worker.log"),
                );
            }
            for variable in &variables {
                command.env("BLASTDB", variable.trim_start_matches("BLASTDB="));
            }
            let output = command
                .current_dir(env!("CARGO_MANIFEST_DIR"))
                .arg("blastp")
                .args(["-db", cells[db], "-query", cells[query]])
                .args(&arguments)
                .args(["-outfmt", spec, "-out"])
                .arg(&out)
                .args(["-num_threads", &threads.to_string()])
                .env("LOSAT_TIMING", "1")
                .env("LOSAT_X_BPQBLOCK", "1")
                .output()
                .expect("run LOSAT");
            let all = String::from_utf8_lossy(&output.stderr).into_owned();
            let searches = all
                .lines()
                .filter(|line| line.starts_with("[TIMING] blastp_total:"))
                .count();
            let stderr: String = all
                .split_inclusive('\n')
                .filter(|line| !line.starts_with("[TIMING]"))
                .collect();
            assert!(output.status.success(), "{what}: {stderr}");
            assert!(searches >= 2, "{what}: {searches} searches, not in blocks");
            assert!(output.stdout.is_empty(), "{what}: standard output");
            assert_eq!(stderr, expected_err, "{what}: standard error differs");
            for (name, bytes) in &expected {
                // `<row>.out_<n>.<ext>` is written as `<row>.out_<n>.<ext>` next to `-out`.
                let written = std::fs::read(dir.join(name)).unwrap_or_default();
                assert!(written == *bytes, "{what}: {name} differs from NCBI's");
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 4);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The queries are read before the first block, so a side file that has the name of the query file
// (`rep_1.json` for `-out rep.json`) truncates it only after it was read: the blocks write what the
// search of all queries writes.
#[test]
fn a_side_file_named_like_the_query_file_changes_nothing() {
    let (query, db) = paja();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let dir = test_dir("alias");
    let text = std::fs::read(&query).unwrap();
    let prepare = |at: &Path| std::fs::write(at.join("rep_1.json"), &text).unwrap();
    let run = |vars: &[(&str, &str)], at: &str| {
        let at = dir.join(at);
        let (run, files) = files_run(
            &at.join("rep_1.json"),
            &subject,
            1,
            vars,
            "13",
            &at,
            &PAJA_OPTIONS,
            &prepare,
        );
        (run, files)
    };
    let reference = run(&[BLOCK_80, XML2_OFF], "ref");
    assert!(reference.0.success, "{}", reference.0.stderr);
    let blocks = run(&[BLOCK_80], "blocks");
    assert_same_files(&blocks, &reference, Some(true), "side file 1 is the query");
    std::fs::remove_dir_all(&dir).unwrap();
}
