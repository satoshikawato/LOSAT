//! blastp's `-db` blocks cut by the database size and the thread count with `LOSAT_X_BPDBCUT`
//! (`src/algorithm/blastp/blast_engine/database/x_db_cut.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`, on by default since 2026-10-17; the
//! searches here set it to `0` unless they name it, `PREVIOUS_DEFAULT`), and it acts only where a
//! `-db` search runs its rounds of blocks on the pool (more than one thread). With it on, the
//! blocks of the preliminary stage hold about a quarter of the database's sequence bytes per
//! thread (at least 4 KiB, at most the 1 MiB of the search without it), each starting at the
//! one-thread order's diagonal offset of its first sequence, and the results are merged in ordinal
//! order. These tests compare the reports with the switch off and on (blastp and blastp-fast, one
//! search of all queries and the block search of `LOSAT_X_BPQBLOCK`); the reach rows of the run
//! with the switch on (`LOSAT_X_REACH`, `tests/reach_support`) show that the cut ran (and, at one
//! thread, that the blocks keep 1 MiB). Each task and case runs at one thread count: the thread
//! count sets only the number of blocks (`blocks_of_the_cut`), the blocks per round and the number
//! of workers, and the report is written from the merged hit lists whatever the outfmt.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413,1553-1554
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ...
//!          /* Save the results. */
//!          status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use LOSAT::seqdb::{SeqDb, SeqType};

mod reach_support;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// What one run printed: the report and standard error.
struct Run {
    code: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
}

/// A protein database of the fixture's subject proteome (`MejoMJNV.faa`, 124 sequences, about 90 KB
/// of sequence: one block of 1 MiB, about 4 x threads blocks of the cut: 9 at 2 threads, 19 at 4), made with
/// `losat makeblastdb` in a directory of the test (`label`, which no other test uses); returns the
/// database name.
fn subject_database(label: &str) -> PathBuf {
    subject_database_with(label, &[])
}

/// `subject_database` with `extra` arguments for `losat makeblastdb` (`-max_file_sz` for volumes).
fn subject_database_with(label: &str, extra: &[&str]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpdbcut-{label}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .args(extra)
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// One blastp run of the fixture pair (two related nimavirus proteomes) with `-task task` against
/// the database `db`, with `vars` as its only `LOSAT_*` variables that choose a code path,
/// `BATCH_SIZE` 3000 (several query batches), `extra` arguments and `-outfmt outfmt`.
/// `LOSAT_X_BPQBLOCK` is 0 (one search of all queries) unless `vars` sets it, so that the counts do
/// not depend on the strict set's default. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with
/// the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name
/// as soon as a search has more than one thread, and it selects nothing; without the caller's value
/// a file of the tests is named. The run is returned whatever its exit code.
fn blastp_any_exit(
    task: &str,
    db: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
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
            std::env::temp_dir().join("losat-bpdbcut-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .args(["-task", task])
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
        .arg("-db")
        .arg(db)
        .args(["-evalue", "0.001", "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .env("BATCH_SIZE", "3000")
        .env("LOSAT_X_BPQBLOCK", "0")
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    Run {
        code: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// `blastp_any_exit`, which must exit with 0.
fn blastp(
    task: &str,
    db: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
) -> Run {
    let run = blastp_any_exit(task, db, threads, vars, extra, outfmt);
    assert!(
        run.code == Some(0),
        "{task} {vars:?} {extra:?} threads {threads}: {:?} {}",
        run.code,
        run.stderr
    );
    run
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

const ON: (&str, &str) = ("LOSAT_X_BPDBCUT", "1");

/// The switch off and on give the same report and standard error, the report of `one_thread` (the
/// same search at one thread without the switch). Returns the reach sides of the run with the
/// switch on (`LOSAT_X_BPDBCUT`: `(the cut ran, the blocks of 1 MiB ran)`) and of
/// `LOSAT_X_BPQBLOCK`.
fn same_reports(
    task: &str,
    db: &Path,
    threads: u32,
    base: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
    one_thread: &Run,
) -> ((bool, bool), (bool, bool)) {
    let what = format!("{task} {base:?} {extra:?} threads {threads} outfmt {outfmt}");
    assert!(
        !one_thread.stdout.is_empty(),
        "{what}: the fixture pair has rows"
    );
    let off = blastp(task, db, threads, base, extra, outfmt);
    assert!(
        off.stdout == one_thread.stdout,
        "{what}: the report differs from one thread's without the switch"
    );
    let reach = reach_support::Reach::new("dbcut");
    let on = blastp(
        task,
        db,
        threads,
        &with(base, &[ON, reach.var()]),
        extra,
        outfmt,
    );
    assert!(
        on.stdout == one_thread.stdout,
        "{what}: the report differs with the switch"
    );
    assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
    (
        reach.sides("LOSAT_X_BPDBCUT"),
        reach.sides("LOSAT_X_BPQBLOCK"),
    )
}

/// The blocks of the cut of the fixture database at `threads` threads: its sequence bytes
/// (letters, one sentinel per sequence and one more: its `.psq` file) over 4 x threads, at least
/// 4 KiB, as `SeqDb::next_block` cuts the sequences.
fn blocks_of_the_cut(db: &Path, threads: u64) -> u64 {
    let seqdb = SeqDb::open(db.to_str().unwrap(), SeqType::Protein).unwrap();
    let bytes = seqdb.total_length() + u64::from(seqdb.num_seqs()) + 1;
    assert_eq!(
        bytes,
        std::fs::metadata(db.with_extension("psq")).unwrap().len()
    );
    let cut = (bytes / (4 * threads)).clamp(4096, 1 << 20);
    let (mut blocks, mut next) = (0, 0);
    while let Some(range) = seqdb.next_block(next, cut) {
        blocks += 1;
        next = range.end_ordinal();
    }
    blocks
}

// blastp at 2 threads and blastp-fast at 4, one search of all queries.
#[test]
fn search_of_all_queries_reports_with_the_cut_are_the_reports_without_it() {
    let db = subject_database("all");
    for (task, threads) in [("blastp", 2u64), ("blastp-fast", 4)] {
        let blocks = blocks_of_the_cut(&db, threads);
        assert!(blocks > 2 * threads, "{threads} threads: {blocks} blocks");
        let one_thread = blastp(task, &db, 1, &[], &[], "6");
        let (cut, _) = same_reports(task, &db, threads as u32, &[], &[], "6", &one_thread);
        assert!(
            cut.0,
            "{task} threads {threads}: the cut did not run: {cut:?}"
        );
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// The block search (`LOSAT_X_BPQBLOCK`, every block of queries a search of the database). Every
// outfmt's search reaches the cut through `database::search_blocks`; outfmt 0 with the cut on runs
// in `pairwise_par.rs`.
#[test]
fn block_search_reports_with_the_cut_are_the_reports_without_it() {
    let db = subject_database("blocks");
    let base = [("LOSAT_X_BPQBLOCK", "1")];
    let one_thread = blastp("blastp-fast", &db, 1, &base, &[], "6");
    let (cut, blocks) = same_reports("blastp-fast", &db, 4, &base, &[], "6", &one_thread);
    assert!(cut.0, "block search: the cut did not run: {cut:?}");
    assert!(
        blocks.0,
        "block search: the block search did not run: {blocks:?}"
    );
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// At one thread the blocks keep 1 MiB: the run with the switch on takes the reference side (the
// cut's rule does not hold) and writes the report of the run without it.
#[test]
fn the_cut_does_not_act_at_one_thread() {
    let db = subject_database("one");
    let off = blastp("blastp", &db, 1, &[], &[], "6");
    assert!(!off.stdout.is_empty(), "the fixture pair has rows");
    let reach = reach_support::Reach::new("dbcut-one");
    let on = blastp("blastp", &db, 1, &[ON, reach.var()], &[], "6");
    assert!(
        on.stdout == off.stdout,
        "the report differs with the switch"
    );
    assert_eq!(on.stderr, off.stderr, "standard error differs");
    let sides = reach.sides("LOSAT_X_BPDBCUT");
    assert_eq!(sides, (false, true), "{sides:?}");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A damaged database (its `.psq` cut short by 500 bytes) through the command line: the search fails
// with the error of the one-thread order, so standard error and the exit code are the same at one
// thread and, with the switch off and on, at 4 threads. The unit test compares
// the texts of the errors; this one also pins their type, which the message of the error maps to.
#[test]
fn a_damaged_database_fails_with_the_same_error_with_the_cut() {
    let db = subject_database("damaged");
    let psq = db.with_extension("psq");
    let bytes = std::fs::read(&psq).unwrap();
    std::fs::write(&psq, &bytes[..bytes.len() - 500]).unwrap();
    let one_thread = blastp_any_exit("blastp", &db, 1, &[], &[], "6");
    assert!(
        one_thread.code.is_some_and(|code| code != 0),
        "the damaged database must fail: {:?} {}",
        one_thread.code,
        one_thread.stderr
    );
    assert!(!one_thread.stderr.is_empty());
    for (what, vars) in [("off", vec![]), ("on", vec![ON])] {
        let run = blastp_any_exit("blastp", &db, 4, &vars, &[], "6");
        assert_eq!(run.code, one_thread.code, "{what}: exit code differs");
        assert_eq!(
            run.stderr, one_thread.stderr,
            "{what}: standard error differs from one thread's"
        );
        assert!(run.stdout == one_thread.stdout, "{what}: report differs");
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A database of several volumes (`-max_file_sz 20000`) with the real kernel: the cut ends blocks at
// volume ends, and 4 threads with the switch on give one thread's report.
#[test]
fn a_database_of_volumes_reports_with_the_cut_are_the_reports_without_it() {
    let db = subject_database_with("volumes", &["-max_file_sz", "20000"]);
    let volumes = (0..100)
        .take_while(|volume| {
            db.with_file_name(format!("MejoMJNV.{volume:02}.psq"))
                .exists()
        })
        .count();
    assert!(volumes >= 3, "{volumes} volumes");
    let one_thread = blastp("blastp", &db, 1, &[], &[], "6");
    assert!(!one_thread.stdout.is_empty());
    let reach = reach_support::Reach::new("dbcut-volumes");
    let on = blastp("blastp", &db, 4, &[ON, reach.var()], &[], "6");
    assert!(
        on.stdout == one_thread.stdout,
        "the report differs with the switch"
    );
    let sides = reach.sides("LOSAT_X_BPDBCUT");
    assert!(sides.0, "the cut did not run: {sides:?}");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}
