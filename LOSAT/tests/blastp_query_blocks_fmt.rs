//! The blastp search in blocks of whole query batches (`LOSAT_X_BPQBLOCK`) for outfmt 0 and 7
//! with `LOSAT_X_BPQBLOCKFMT` (`src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! Both switches are in the strict set (`src/utils/xswitch.rs`): on by default, off with `0` or
//! with `LOSAT_X_DEFAULT=0` and the variable not set. With both on, a search whose one report is
//! outfmt 0, 6 or 7 searches its queries block by block and writes each block's report before it
//! searches the next block: the outfmt 0 prolog with the first block, the epilog and outfmt 7's
//! count of the queries with the last. These tests compare the standard output and the standard
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

/// One blastp run.
struct Run {
    success: bool,
    stdout: Vec<u8>,
    /// Standard error without the `[TIMING]` and `[X_SHADOW]` lines.
    stderr: String,
    /// The `[TIMING] blastp_total` lines: one per search of a query set.
    searches: usize,
    /// The `[X_SHADOW] BPQBLOCKFMT` lines.
    shadow: Vec<String>,
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
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpqblockfmt-worker.log"),
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
    let all = String::from_utf8_lossy(&output.stderr).into_owned();
    let searches = all
        .lines()
        .filter(|line| line.starts_with("[TIMING] blastp_total:"))
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
        stdout: output.stdout,
        stderr,
        searches,
        shadow,
    }
}

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
// queries.
#[test]
fn outfmt_0_and_7_blocks_write_the_report_of_the_search_of_all_queries() {
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("reports");
    for subject in [
        [OsStr::new("-subject"), fasta.as_os_str()],
        [OsStr::new("-db"), db.as_os_str()],
    ] {
        for spec in ["0", "7", "7 qseqid sseqid pident evalue bitscore btop"] {
            let fmt = outfmt(spec);
            for (threads, letters) in [(1, "1"), (4, "20000")] {
                let what = format!("{:?} outfmt {spec} threads {threads}", subject[0]);
                let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
                assert!(
                    !reference.stdout.is_empty(),
                    "{what}: the fixture pair has a report"
                );
                let run = blastp(
                    &query,
                    &subject,
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
                if threads == 4 {
                    let run = blastp(
                        &query,
                        &subject,
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

/// A file of the test's directory `dir` holding `text`.
fn write_file(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

// Queries with warnings (O read as X, a record without letters) in several blocks: each query's
// warnings come before its report, in the same bytes, in every format.
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
    for spec in ["0", "7", "6"] {
        let fmt = outfmt(spec);
        for threads in THREADS {
            let what = format!("outfmt {spec} threads {threads}");
            let reference = blastp(&query, &subject, threads, &[BATCH, SWITCH_OFF], &os(&fmt));
            assert!(
                reference.stderr.contains("Warning:"),
                "{what}: the queries have warnings: {}",
                reference.stderr
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
// and says so in one line (it panics on a difference).
#[test]
fn the_out_file_and_the_shadow() {
    let dir = std::env::temp_dir().join(format!("losat-bpqblockfmt-out-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = fixture("MeenMJNV.faa");
    let fasta = fixture("MejoMJNV.faa");
    let db = subject_database("shadow");
    for subject in [
        [OsStr::new("-subject"), fasta.as_os_str()],
        [OsStr::new("-db"), db.as_os_str()],
    ] {
        for spec in ["0", "6", "7"] {
            let fmt = outfmt(spec);
            let what = format!("{:?} outfmt {spec}", subject[0]);
            let reference = blastp(&query, &subject, 4, &[BATCH, SWITCH_OFF], &os(&fmt));
            let path = dir.join(format!("report{spec}.txt"));
            let mut extra = fmt.to_vec();
            extra.extend([OsString::from("-out"), path.clone().into_os_string()]);
            let run = blastp(
                &query,
                &subject,
                4,
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
            for threads in THREADS {
                let run = blastp(
                    &query,
                    &subject,
                    threads,
                    &[BATCH, ("LOSAT_X_BPQBLOCK", "1"), FMT_OFF, SHADOW],
                    &os(&fmt),
                );
                assert!(
                    run.success,
                    "{what} shadow threads {threads}: {}",
                    run.stderr
                );
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
        }
    }
    std::fs::remove_dir_all(db.parent().unwrap()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}
