//! blastp-fast's fused loop on the flat compressed table with each cell's hits staged before the
//! two-hit test and the scan two offsets per step (`LOSAT_X_BPFSTAGE`,
//! `src/algorithm/blastp/x_seed_fused_clut.rs` `x_scan_flat_stage`) on whole searches, and its
//! form with the scan's steps and the copy out of line (`LOSAT_X_BPFLEAN`, `x_scan_flat_lean`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518
//! ```c
//!         for (i = 0; i < hits; ++i) {
//!             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ```
//! Both switches are in the strict set (on unless they are `0`). The loops make the candidates of
//! the loop without them and leave the same diagonal table (the lib test
//! `the_flat_loop_is_the_fused_loop_on_the_reference_table`), so the property under test is that
//! the report does not depend on them, for subjects read from FASTA (`-subject`) and from a database
//! (`-db`), at 1 and 2 threads. That reports equal NCBI's is the business of the BLASTP fixtures.
//! `LOSAT_X_REACH` shows the side each run took: the loop runs where LOSAT_X_BPFSEED's loop runs
//! (the dense query's table, and every table with `LOSAT_X_BPCLUTSCANBFMIN=0`; every subject here
//! keeps its sentinel); the sparse query's table, every table with `LOSAT_X_BPCLUTSCANBF=0`, and
//! the reference compressed table (`LOSAT_X_BPCLUT=0`) keep their loops, and with
//! `LOSAT_X_BPFSEED=0` the switch does not act. `LOSAT_X_BPFLEAN` takes the sides of
//! `LOSAT_X_BPFSTAGE` and does not act with `LOSAT_X_BPFSTAGE=0`; with `LOSAT_X_BPFLEAN=0` the
//! staged loop runs. Exact blastp builds no compressed table and records no side.
//!
//! The dense query is the first 30 proteins of AP027131 (9,948 residues; its table uses 30.9 % of
//! the cells), the sparse one its first protein (455 residues; 1.66 %); the subjects are the first
//! 150 proteins of the related AP027133.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_BPFSTAGE";
const OFF: (&str, &str) = ("LOSAT_X_BPFSTAGE", "0");
const LEAN: &str = "LOSAT_X_BPFLEAN";
const LEAN_OFF: (&str, &str) = ("LOSAT_X_BPFLEAN", "0");
const SEED_OFF: (&str, &str) = ("LOSAT_X_BPFSEED", "0");
const BRANCH_LIST: (&str, &str) = ("LOSAT_X_BPCLUTSCANBF", "0");
const EVERY_TABLE: (&str, &str) = ("LOSAT_X_BPCLUTSCANBFMIN", "0");
const REFERENCE_TABLE: (&str, &str) = ("LOSAT_X_BPCLUT", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_fast_stage_{}_{}_{tag}",
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
    dense: PathBuf,
    sparse: PathBuf,
    subject: PathBuf,
    db: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let dense = dir.join("dense.faa");
        let sparse = dir.join("sparse.faa");
        let subject = dir.join("subject.faa");
        std::fs::write(&dense, records("AP027131.faa", 0, 30)).unwrap();
        std::fs::write(&sparse, records("AP027131.faa", 0, 1)).unwrap();
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
            dense,
            sparse,
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
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
    /// The sides of the switch the run recorded: (new, ref).
    reached: (bool, bool),
    /// The sides of `LOSAT_X_BPFLEAN` the run recorded.
    lean: (bool, bool),
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
    let (mut new, mut reference) = (false, false);
    let mut lean = (false, false);
    for file in std::fs::read_dir(&dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            if f[0] == SWITCH {
                new |= f[1] == "1";
                reference |= f[2] == "1";
            }
            if f[0] == LEAN {
                lean.0 |= f[1] == "1";
                lean.1 |= f[2] == "1";
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&worker_log);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
        reached: (new, reference),
        lean,
    }
}

/// Default, `=0` and `LOSAT_X_BPFLEAN=0` give the same report and messages for blastp-fast on
/// `-subject` and `-db` at 1 and 2 threads, and so do the other lists of the scan, the reference
/// compressed table and `LOSAT_X_BPFSEED=0` on one case each; a flat table with the branch-free
/// list takes the loop (the lean one by default), the other tables keep theirs, and `=0` records
/// nothing.
#[test]
fn blastp_fast_reports_do_not_depend_on_the_staged_loop() {
    let inputs = Inputs::new();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    for (query, dense) in [(&inputs.dense, true), (&inputs.sparse, false)] {
        let query = query.to_str().unwrap();
        for (source, target) in [("-subject", subject), ("-db", db)] {
            for threads in ["1", "2"] {
                let case = format!("dense {dense} {source} n{threads}");
                let args = [
                    "-task",
                    "blastp-fast",
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
                assert_eq!(on.reached, (dense, !dense), "{case}: default's side");
                assert_eq!(on.lean, (dense, !dense), "{case}: default's BPFLEAN side");
                // The tables' other choices on one case each: the sparse table with the branch-free
                // list, the dense one with the list with the branch and the reference table.
                // (label, variable, the sides of LOSAT_X_BPFSTAGE and of LOSAT_X_BPFLEAN)
                let first = source == "-subject" && threads == "1";
                let mut arms = vec![("off", OFF, (false, false), (false, false))];
                // The sparse table does not take the loop: its LOSAT_X_BPFLEAN=0 on one case.
                if dense || first {
                    arms.push(("lean off", LEAN_OFF, (dense, !dense), (false, false)));
                }
                if first && !dense {
                    arms.push(("every table", EVERY_TABLE, (true, false), (true, false)));
                }
                if first && dense {
                    arms.push(("branch list", BRANCH_LIST, (false, true), (false, true)));
                    arms.push((
                        "reference table",
                        REFERENCE_TABLE,
                        (false, true),
                        (false, true),
                    ));
                    arms.push(("BPFSEED off", SEED_OFF, (false, false), (false, false)));
                }
                for (label, env, reached, lean) in arms {
                    let run = blastp(&args, &[env]);
                    assert_eq!(run.code, Some(0), "{case} {label}: {}", run.stderr);
                    assert!(run.stdout == on.stdout, "{case} {label}: report differs");
                    assert_eq!(run.stderr, on.stderr, "{case} {label}: messages");
                    assert_eq!(run.reached, reached, "{case} {label}: side");
                    assert_eq!(run.lean, lean, "{case} {label}: BPFLEAN side");
                }
            }
        }
    }
}

/// Exact blastp builds no compressed table: the same report with `=0`, and no side recorded.
#[test]
fn exact_blastp_records_no_side() {
    let inputs = Inputs::new();
    let args = [
        "-query",
        inputs.sparse.to_str().unwrap(),
        "-subject",
        inputs.subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    let on = blastp(&args, &[]);
    assert_eq!(on.code, Some(0), "{}", on.stderr);
    assert!(!on.stdout.is_empty(), "no hits");
    assert_eq!(on.reached, (false, false), "default's side");
    assert_eq!(on.lean, (false, false), "default's BPFLEAN side");
    let off = blastp(&args, &[OFF]);
    assert_eq!(off.code, Some(0), "{}", off.stderr);
    assert!(off.stdout == on.stdout, "report differs");
    assert_eq!(off.reached, (false, false), "=0 records a side");
    assert_eq!(off.lean, (false, false), "=0 records a BPFLEAN side");
}
