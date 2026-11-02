//! blastp-fast's flat compressed scan lists the subject offsets without a branch on the presence
//! bit on every query's table (`LOSAT_X_BPCLUTSCANBF`,
//! `src/algorithm/tblastx/lookup/x_clut_direct.rs` `for_each_cell`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341
//! ```c
//!        index = preshift + compressed_char;
//!        preshift = (Int4)((((Int8)( index )) * recip) >> 32);
//!
//!        /* if there are hits */
//!        if (PV_TEST(pv, index, pv_array_bts)) {
//! ```
//! The switch is in the strict set (on unless it is `0`). The branch-free list holds the offsets
//! and indices of the loop with the branch, in the same order (the lib test
//! `flat_scans_give_the_reference_scans_hits` compares both with the reference scan), so the
//! property under test is that the report does not depend on it, for subjects read from FASTA
//! (`-subject`) and from a database (`-db`), at 1 and 2 threads. That reports equal NCBI's is the
//! business of the BLASTP fixtures. Every flat compressed table takes the list (until 2026-11-02 a
//! table with fewer than 500 used cells per 10,000 kept the branch natively); only
//! `LOSAT_X_BPCLUTSCANBF=0` and the reference compressed table (`LOSAT_X_BPCLUT=0`) keep the loop
//! with the branch. Exact blastp builds no compressed table, and the switch does not act.
//!
//! The loops after the list (blastp-fast's fused loop with the branch-free two-hit test and the
//! staged copy, `src/algorithm/blastp/x_seed_fused_clut.rs`; always taken, natively, on a table
//! with the branch-free list) are part of the default run here, and the report must not change
//! with the loop with the branch. The subjects include records shorter than a word, which take
//! the branch-free test's loop without the staged copy.
//!
//! With `LOSAT_X_STATS` (and the `xstats` feature, as in `cargo test --all-features`) the
//! `[X_STATS] BPFLAT lean= fseed= bf= branch=` line counts the subjects each loop on the flat table
//! scanned, so the default arm is shown to run `x_scan_flat_lean` (the proteins) and
//! `x_scan_flat_fseed_loop` (the records shorter than a word, natively), and the `=0` arm the loop with
//! the branch alone; the reference table prints no such line; a build without the feature prints
//! none. `LOSAT_X_STATS` changes neither the report nor the other messages. A build without the
//! `parallel` feature runs the one-thread cases only.
//!
//! The dense query is the first 30 proteins of AP027131 (9,948 residues; its table uses 30.9 % of
//! the cells), the sparse one its first protein (455 residues; 1.66 %, a table that kept the
//! branch under the old rule); the subjects are the first 150 proteins of the related AP027133.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_BPCLUTSCANBF", "0");
const REFERENCE_TABLE: (&str, &str) = ("LOSAT_X_BPCLUT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_fast_scan_bf_{}_{}_{tag}",
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
        // Records of 1, 3 and 4 residues (shorter than a word of 5) follow the proteins.
        let mut subjects = records("AP027133.faa", 0, 150);
        subjects.push_str(">short1\nM\n>short3\nMKV\n>short4\nMKVL\n");
        std::fs::write(&subject, subjects).unwrap();
        let db = dir.join("db");
        let made = losat(&[], &dir.join("worker.log"))
            .arg("makeblastdb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run LOSAT makeblastdb");
        assert!(
            made.status.success(),
            "makeblastdb: {}",
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
}

impl Run {
    /// Standard error without the `[X_STATS]` lines.
    fn messages(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[X_STATS]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// The number after `name=` in the `[X_STATS] BPFLAT` line, or `None` without that line.
    fn flat(&self, name: &str) -> Option<u64> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with("[X_STATS] BPFLAT "))?;
        let key = format!(" {name}=");
        let at = line
            .find(&key)
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            + key.len();
        Some(
            line[at..]
                .split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap(),
        )
    }
}

/// `blastp` with `env` (the worker log is in a directory of its own, removed after the run).
fn blastp(args: &[&str], env: &[(&str, &str)]) -> Run {
    let worker_log = temp_path("worker");
    std::fs::create_dir_all(&worker_log).unwrap();
    let output = losat(env, &worker_log.join("worker.log"))
        .arg("blastp")
        .args(args)
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_dir_all(&worker_log);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// Default and `=0` (and, on one case, the reference compressed table) give the same report and
/// messages for blastp-fast on `-subject` and `-db` at 1 and 2 threads, for the dense and the
/// sparse query: the default takes the branch-free list on both tables.
#[test]
fn blastp_fast_reports_do_not_depend_on_the_branch_free_list() {
    let inputs = Inputs::new();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    // The `xstats` feature counts the subjects; a build without it prints no BPFLAT line.
    let counted = cfg!(feature = "xstats");
    for (query, dense) in [(&inputs.dense, true), (&inputs.sparse, false)] {
        let query = query.to_str().unwrap();
        for (source, target) in [("-subject", subject), ("-db", db)] {
            for threads in ["1", "2"] {
                // A build without the `parallel` feature refuses two threads.
                if threads == "2" && !cfg!(feature = "parallel") {
                    continue;
                }
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
                let mut arms = vec![("off", OFF)];
                // The reference compressed table, on one case: it builds no flat table and keeps
                // its loop.
                if dense && source == "-subject" && threads == "1" {
                    arms.push(("reference table", REFERENCE_TABLE));
                }
                for (label, env) in arms {
                    let run = blastp(&args, &[env]);
                    assert_eq!(run.code, Some(0), "{case} {label}: {}", run.stderr);
                    assert!(run.stdout == on.stdout, "{case} {label}: report differs");
                    assert_eq!(run.stderr, on.stderr, "{case} {label}: messages");
                }
                // The loops that ran, counted per subject. A build without the feature prints no
                // line, with `LOSAT_X_STATS` or without.
                let on_stats = blastp(&args, &[STATS]);
                assert_eq!(on_stats.code, Some(0), "{case} stats: {}", on_stats.stderr);
                assert!(on_stats.stdout == on.stdout, "{case} stats: report differs");
                assert_eq!(on_stats.messages(), on.stderr, "{case} stats: messages");
                let off_stats = blastp(&args, &[STATS, OFF]);
                assert_eq!(
                    off_stats.code,
                    Some(0),
                    "{case} off stats: {}",
                    off_stats.stderr
                );
                assert!(
                    off_stats.stdout == on.stdout,
                    "{case} off stats: report differs"
                );
                assert_eq!(
                    off_stats.messages(),
                    on.stderr,
                    "{case} off stats: messages"
                );
                let reference_stats = blastp(&args, &[STATS, REFERENCE_TABLE]);
                assert_eq!(
                    reference_stats.code,
                    Some(0),
                    "{case}: {}",
                    reference_stats.stderr
                );
                assert!(
                    reference_stats.stdout == on.stdout,
                    "{case} ref stats: report differs"
                );
                // The reference table runs none of the flat table's loops.
                assert_eq!(
                    reference_stats.flat("lean"),
                    None,
                    "{case}: {}",
                    reference_stats.stderr
                );
                if counted {
                    let loops = ["lean", "fseed", "bf", "branch"];
                    let on_n = loops.map(|name| on_stats.flat(name));
                    let off_n = loops.map(|name| off_stats.flat(name));
                    let msg = format!("{case}: {}", on_stats.stderr);
                    // Each subject is scanned once per search (150 proteins and the three short
                    // records). Default: the proteins take the list two offsets per step, the
                    // records of 1, 3 and 4 residues the loop for subjects shorter than a word
                    // (natively; wasm32 runs `bf` on every subject), none the loop with the branch.
                    assert_eq!(on_n[0], Some(150), "{msg}");
                    assert_eq!(on_n[1], Some(3), "{msg}");
                    assert_eq!(on_n[2], Some(0), "{msg}");
                    assert_eq!(on_n[3], Some(0), "{msg}");
                    // `=0`: the loop with the branch alone, on every subject.
                    let msg = format!("{case} off: {}", off_stats.stderr);
                    assert_eq!(off_n[3], Some(153), "{msg}");
                    assert_eq!(off_n[..3], [Some(0); 3], "{msg}");
                } else {
                    assert_eq!(on_stats.flat("lean"), None, "{case}: {}", on_stats.stderr);
                    assert_eq!(
                        off_stats.flat("branch"),
                        None,
                        "{case}: {}",
                        off_stats.stderr
                    );
                }
            }
        }
    }
}

/// Exact blastp builds no compressed table: the same report with `=0`.
#[test]
fn exact_blastp_does_not_depend_on_the_switch() {
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
    let off = blastp(&args, &[OFF]);
    assert_eq!(off.code, Some(0), "{}", off.stderr);
    assert!(off.stdout == on.stdout, "report differs");
}
