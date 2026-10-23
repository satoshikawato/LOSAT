//! blastp-fast's flat compressed scan lists the subject offsets without a branch on the presence
//! bit on dense tables (`LOSAT_X_BPSCANBF`, `src/algorithm/tblastx/lookup/x_clut_direct.rs`
//! `for_each_cell`, `XFlatCompressedLookup::x_scan_bf_pays`) on whole searches.
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
//! business of the BLASTP fixtures. `LOSAT_X_REACH` shows the side each run took: a query block's
//! table with at least `LOSAT_X_BPSCANBFMIN` used cells per 10,000 (built-in 500) takes the
//! branch-free list, a sparser one keeps the branch, and `LOSAT_X_BPSCANBFMIN=0` takes the list
//! on every table. Exact blastp's standard table always takes its branch-free scan.
//!
//! The dense query is the first 30 proteins of AP027131 (9,948 residues), the sparse one its
//! first protein (455 residues); the subjects are the first 150 proteins of the related AP027133.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_BPSCANBF";
const OFF: (&str, &str) = ("LOSAT_X_BPSCANBF", "0");
const EVERY_TABLE: (&str, &str) = ("LOSAT_X_BPSCANBFMIN", "0");

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
        std::fs::write(&subject, records("AP027133.faa", 0, 150)).unwrap();
        let db = dir.join("db");
        let made = losat(&[])
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
/// names a file of its own: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one thread,
/// and it selects nothing.
fn losat(env: &[(&str, &str)]) -> Command {
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
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
    /// The sides of the switch the run recorded: (new, ref).
    reached: (bool, bool),
}

/// `blastp` with `env` and `LOSAT_X_REACH` set to a directory of its own.
fn blastp(args: &[&str], env: &[(&str, &str)]) -> Run {
    let dir = temp_path("reach");
    std::fs::create_dir_all(&dir).unwrap();
    let output = losat(env)
        .env("LOSAT_X_REACH", &dir)
        .arg("blastp")
        .args(args)
        .output()
        .expect("run LOSAT");
    let (mut new, mut reference) = (false, false);
    for file in std::fs::read_dir(&dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            if f[0] == SWITCH {
                new |= f[1] == "1";
                reference |= f[2] == "1";
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
        reached: (new, reference),
    }
}

/// Default, `=0` and `LOSAT_X_BPSCANBFMIN=0` give the same report and messages for blastp-fast
/// on `-subject` and `-db` at 1 and 2 threads; the dense query's table takes the branch-free list
/// and the sparse one's keeps the branch unless every table takes the list; `=0` records nothing.
#[test]
fn blastp_fast_reports_do_not_depend_on_the_branch_free_list() {
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
                for (label, env, reached) in [
                    ("off", OFF, (false, false)),
                    ("every table", EVERY_TABLE, (true, false)),
                ] {
                    let run = blastp(&args, &[env]);
                    assert_eq!(run.code, Some(0), "{case} {label}: {}", run.stderr);
                    assert!(run.stdout == on.stdout, "{case} {label}: report differs");
                    assert_eq!(run.stderr, on.stderr, "{case} {label}: messages");
                    assert_eq!(run.reached, reached, "{case} {label}: side");
                }
            }
        }
    }
}

/// Exact blastp's fused loop on the standard table always takes its branch-free scan.
#[test]
fn exact_blastp_takes_the_branch_free_scan() {
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
    assert_eq!(on.reached, (true, false), "default's side");
    let off = blastp(&args, &[OFF]);
    assert_eq!(off.code, Some(0), "{}", off.stderr);
    assert!(off.stdout == on.stdout, "report differs");
    assert_eq!(off.reached, (false, false), "=0 records a side");
}
