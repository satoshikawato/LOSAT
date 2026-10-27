//! What blastp holds during its composition-based redo: the search's store of SEG-masked subject
//! ranges emptied after a window of the redo in windows once it holds the switch's number of
//! entries (`LOSAT_X_BPSEGWIN`, `src/algorithm/blastp/x_redo_par.rs` `x_bp_seg_win`), and the lookup table dropped before the
//! redo where no later block's build takes its vectors (`LOSAT_X_BPLKDROP`,
//! `src/algorithm/blastp/blast_engine.rs` `x_lookup_after_prelim`), and the vectors of a block's
//! table taken by the next block's own build where no table is built ahead (`LOSAT_X_BPLKREUSE`,
//! `src/algorithm/blastp/x_query_blocks.rs` `x_bp_aa_block_spare`), on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
//! ```c
//!                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//!                                              subject_maybe_biased);
//! ```
//! Both switches are in the strict set (on unless `0`) and change only what is kept in memory: a
//! masked range is a function of its subject, range and SEG parameters, so a window that masks a
//! range again gets the same bytes (the lib test
//! `windows_with_the_seg_store_emptied_after_each_equal_windows_with_one_store`), and nothing reads
//! the lookup table after the preliminary stage. So the property under test is that the report,
//! the messages and the exit code do not depend on them, with subjects read from FASTA (`-subject`)
//! and from a database (`-db`), at 1 and 2 threads. That reports equal NCBI's is the business of
//! the BLASTP fixtures. `LOSAT_X_REACH` shows the side each run took: the store is emptied where
//! the redo runs in windows (`LOSAT_X_BPREDOPAR`: more than one thread) and the store holds the
//! switch's number of entries after a window (built in 4,096; these searches hold fewer, and
//! `LOSAT_X_BPSEGWIN=2` empties it), and kept otherwise;
//! the table is dropped by a search of all queries and by the last blocks of the block search
//! (`LOSAT_X_BPQBLOCK`), and kept by a block whose table's vectors a later block's build takes
//! (more than one thread: the build ahead of the block after the next, `LOSAT_X_BPLKPAHEAD`; one
//! thread: the next block's own build, `LOSAT_X_BPLKREUSE`, which records its side per block
//! search: its own at one thread, the other where the tables are built ahead).
//!
//! The query is the first 30 proteins of AP027131 (9,948 residues: one query batch), the subjects
//! the first 150 proteins of the related AP027133. The block search reads the first 100 proteins
//! (30,926 residues), with `LOSAT_X_BPQBLOCK=1`: one block per query batch.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SEG_WIN: &str = "LOSAT_X_BPSEGWIN";
const LK_DROP: &str = "LOSAT_X_BPLKDROP";
const LK_REUSE: &str = "LOSAT_X_BPLKREUSE";
const BLOCKS: (&str, &str) = ("LOSAT_X_BPQBLOCK", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_redo_memory_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Records `from..to` (0-based) of a FASTA file, as text.
fn records(name: &str, from: usize, to: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut out = String::new();
    let mut record = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            record += 1;
        }
        if record > from && record <= to {
            out.push_str(line);
            out.push('\n');
        }
    }
    assert!(record >= to, "{name} has fewer than {to} records");
    out
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    /// The first 100 proteins, for the block search.
    blocks_query: PathBuf,
    subject: PathBuf,
    db: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("query.faa");
        let blocks_query = dir.join("blocks_query.faa");
        let subject = dir.join("subject.faa");
        std::fs::write(&query, records("AP027131.faa", 0, 30)).unwrap();
        std::fs::write(&blocks_query, records("AP027131.faa", 0, 100)).unwrap();
        std::fs::write(&subject, records("AP027133.faa", 0, 150)).unwrap();
        let db = dir.join("db");
        let made = losat(&[], &dir.join("worker.log"))
            .arg("makedb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run LOSAT makedb");
        assert!(
            made.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&made.stderr)
        );
        Inputs {
            dir,
            query,
            blocks_query,
            subject,
            db,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The binary with `env` as its only variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// names `worker_log`, in a directory the caller removes: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn losat(env: &[(&str, &str)], worker_log: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", worker_log)
        // `LOSAT_X_BPQSLOTS=0`: the block loop of `LOSAT_X_BPQBLOCK`, one block at a time (the switch searches
        // blastp-fast's blocks two at a time, of other sizes, with no table built ahead), unless the run names it.
        .env("LOSAT_X_BPQSLOTS", "0")
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
    /// The sides each switch's run recorded: (new, ref) of `SEG_WIN`, `LK_DROP` and `LK_REUSE`.
    seg_win: (bool, bool),
    lk_drop: (bool, bool),
    lk_reuse: (bool, bool),
}

/// `blastp` with `env` and `LOSAT_X_REACH` set to a directory of its own (removed after the run,
/// with the worker log).
fn blastp(args: &[&str], env: &[(&str, &str)]) -> Run {
    let dir = temp_path("reach");
    std::fs::create_dir_all(&dir).unwrap();
    let worker_log = temp_path("worker");
    std::fs::create_dir_all(&worker_log).unwrap();
    let output = losat(env, &worker_log.join("worker.log"))
        .env("LOSAT_X_REACH", &dir)
        .arg("blastp")
        .args(args)
        .output()
        .expect("run LOSAT");
    let (mut seg_win, mut lk_drop, mut lk_reuse) = ((false, false), (false, false), (false, false));
    for file in std::fs::read_dir(&dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            let sides = match f[0] {
                SEG_WIN => &mut seg_win,
                LK_DROP => &mut lk_drop,
                LK_REUSE => &mut lk_reuse,
                _ => continue,
            };
            sides.0 |= f[1] == "1";
            sides.1 |= f[2] == "1";
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&worker_log);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
        seg_win,
        lk_drop,
        lk_reuse,
    }
}

/// Default, `LOSAT_X_BPSEGWIN=2` and `=0` give the same report and messages for exact blastp on
/// `-subject` and `-db` at 1 and 2 threads. The default keeps the store (fewer than 4,096
/// entries); `=2` empties it after a window at 2 threads (the redo in windows) and not at 1; `=0`
/// records nothing.
#[test]
fn blastp_reports_do_not_depend_on_emptying_the_seg_store_after_each_window() {
    let inputs = Inputs::new();
    let query = inputs.query.to_str().unwrap();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    for (source, target) in [("-subject", subject), ("-db", db)] {
        for threads in ["1", "2"] {
            let case = format!("{source} n{threads}");
            let args = [
                "-query",
                query,
                source,
                target,
                "-outfmt",
                "6",
                "-evalue",
                "10",
                "-num_threads",
                threads,
            ];
            let on = blastp(&args, &[]);
            assert_eq!(on.code, Some(0), "{case}: {}", on.stderr);
            assert!(!on.stdout.is_empty(), "{case}: no hits");
            assert_eq!(on.seg_win, (false, true), "{case}: default's side");
            let windows = threads == "2";
            for (value, side) in [("2", (windows, !windows)), ("0", (false, false))] {
                let arm = blastp(&args, &[(SEG_WIN, value)]);
                assert_eq!(arm.code, Some(0), "{case} ={value}: {}", arm.stderr);
                assert!(arm.stdout == on.stdout, "{case} ={value}: report differs");
                assert_eq!(arm.stderr, on.stderr, "{case} ={value}: messages");
                assert_eq!(arm.seg_win, side, "{case} ={value}: side");
            }
        }
    }
}

/// Default and `LOSAT_X_BPLKDROP=0` give the same report and messages for exact blastp and
/// blastp-fast on `-subject` and `-db`: in one search of all queries at 1 and 2 threads, and in
/// blocks at 2 threads. The table is dropped before the redo by every search except the blocks
/// whose table's vectors a later block's build takes (all but the last two), and `=0` records
/// nothing.
#[test]
fn blastp_reports_do_not_depend_on_dropping_the_lookup_table_before_the_redo() {
    let inputs = Inputs::new();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    for task in ["blastp", "blastp-fast"] {
        for (source, target) in [("-subject", subject), ("-db", db)] {
            for (threads, blocks) in [("1", false), ("2", false), ("2", true)] {
                let query = if blocks {
                    inputs.blocks_query.to_str().unwrap()
                } else {
                    inputs.query.to_str().unwrap()
                };
                let case = format!("{task} {source} n{threads} blocks {blocks}");
                let args = [
                    "-task",
                    task,
                    "-query",
                    query,
                    source,
                    target,
                    "-outfmt",
                    "6",
                    "-evalue",
                    "10",
                    "-num_threads",
                    threads,
                ];
                let env: &[(&str, &str)] = if blocks { &[BLOCKS] } else { &[] };
                let on = blastp(&args, env);
                assert_eq!(on.code, Some(0), "{case}: {}", on.stderr);
                assert!(!on.stdout.is_empty(), "{case}: no hits");
                assert_eq!(on.lk_drop, (true, blocks), "{case}: default's side");
                let off = blastp(&args, &[env, &[(LK_DROP, "0")]].concat());
                assert_eq!(off.code, Some(0), "{case} off: {}", off.stderr);
                assert!(off.stdout == on.stdout, "{case} off: report differs");
                assert_eq!(off.stderr, on.stderr, "{case} off: messages");
                assert_eq!(off.lk_drop, (false, false), "{case} off: side");
            }
        }
    }
}

/// Default and `LOSAT_X_BPLKREUSE=0` give the same report and messages for exact blastp on
/// `-subject` and `-db` in blocks, at 1 and 2 threads. At 1 thread each block's build takes the
/// vectors of the block before: every block but the last keeps its table to its search's end and
/// the last drops it before the redo (`LOSAT_X_BPLKDROP` records both sides; with `=0` every block
/// drops it). At 2 threads the tables are built ahead (`LOSAT_X_BPLKPAHEAD`), and the switch
/// records its other side. `=0` records nothing.
#[test]
fn blastp_reports_do_not_depend_on_reusing_the_vectors_of_the_block_before() {
    let inputs = Inputs::new();
    let query = inputs.blocks_query.to_str().unwrap();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    for (source, target) in [("-subject", subject), ("-db", db)] {
        for threads in ["1", "2"] {
            let case = format!("{source} n{threads}");
            let args = [
                "-query",
                query,
                source,
                target,
                "-outfmt",
                "6",
                "-evalue",
                "10",
                "-num_threads",
                threads,
            ];
            let one = threads == "1";
            let on = blastp(&args, &[BLOCKS]);
            assert_eq!(on.code, Some(0), "{case}: {}", on.stderr);
            assert!(!on.stdout.is_empty(), "{case}: no hits");
            assert_eq!(on.lk_reuse, (one, !one), "{case}: default's side");
            assert_eq!(on.lk_drop, (true, true), "{case}: LK_DROP's sides");
            let off = blastp(&args, &[BLOCKS, (LK_REUSE, "0")]);
            assert_eq!(off.code, Some(0), "{case} off: {}", off.stderr);
            assert!(off.stdout == on.stdout, "{case} off: report differs");
            assert_eq!(off.stderr, on.stderr, "{case} off: messages");
            assert_eq!(off.lk_reuse, (false, false), "{case} off: side");
            assert_eq!(off.lk_drop, (true, !one), "{case} off: LK_DROP's sides");
        }
    }
}
