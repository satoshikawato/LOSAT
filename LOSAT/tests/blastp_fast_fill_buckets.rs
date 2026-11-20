//! blastp-fast's flat compressed table filled by buckets of cells (`LOSAT_X_BPFFILL`,
//! `src/algorithm/tblastx/lookup/x_clut_direct.rs` `x_try_build_flat_by_buckets`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
//! ```c
//!     for (loc = location; loc; loc = loc->next){
//!         Int4 from = loc->ssr->left;
//!         Int4 to = loc->ssr->right - lookup->word_length + 1;
//!
//!         for (offset = from; offset <= to; offset++){
//!             s_CompressedAddWordHits(&info, query->sequence, offset);
//!         }
//!     }
//! ```
//! The switch is in the strict set (on unless it is `0`). The bucket fill builds the table of the
//! count-and-fill builder slot for slot (the lib test `x_bucket_fill_is_the_count_and_fill_build`),
//! so the property under test is that the report does not depend on it, for subjects read from
//! FASTA (`-subject`) and from a database (`-db`). That reports equal NCBI's is the business of the
//! BLASTP fixtures. `LOSAT_X_REACH` shows the side each run took: at one thread (and at any thread
//! count in a build without the `parallel` feature) the serial builder runs and takes the bucket
//! fill (its other side needs a query of 2^22 letters or more given to one table, or a table past
//! NCBI's overflow capacity); at two threads the pool builds the table and no side is recorded, nor
//! with the reference compressed table (`LOSAT_X_BPCLUT=0`) or with exact blastp, which builds no
//! compressed table.
//!
//! The dense query is the first 30 proteins of AP027131 (9,948 residues; its table uses 30.9 % of
//! the cells, cells of more than 5 offsets among them), the sparse one its first protein (455
//! residues; 1.66 %); the subjects are the first 150 proteins of the related AP027133.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_BPFFILL";
const OFF: (&str, &str) = ("LOSAT_X_BPFFILL", "0");
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
        "losat_blastp_fast_fill_buckets_{}_{}_{tag}",
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
    /// The sides of the switch the run recorded: (new, ref).
    reached: (bool, bool),
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
    let _ = std::fs::remove_dir_all(&worker_log);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
        reached: (new, reference),
    }
}

/// Default and `=0` give the same report and messages for blastp-fast on `-subject` and `-db` at 1
/// and 2 threads, and so does the reference compressed table on one case; the serial builder (one
/// thread, or a build without threads) takes the bucket fill, the pool's build (two threads) and the
/// reference table record no side, and `=0` records nothing.
#[test]
fn blastp_fast_reports_do_not_depend_on_the_bucket_fill() {
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
                let serial = threads == "1" || !cfg!(feature = "parallel");
                assert_eq!(on.reached, (serial, false), "{case}: default's side");
                let mut arms = vec![("off", OFF)];
                if dense && source == "-subject" && threads == "1" {
                    arms.push(("reference table", REFERENCE_TABLE));
                }
                for (label, env) in arms {
                    let run = blastp(&args, &[env]);
                    assert_eq!(run.code, Some(0), "{case} {label}: {}", run.stderr);
                    assert!(run.stdout == on.stdout, "{case} {label}: report differs");
                    assert_eq!(run.stderr, on.stderr, "{case} {label}: messages");
                    assert_eq!(run.reached, (false, false), "{case} {label}: side");
                }
            }
        }
    }
}

/// Exact blastp builds no compressed table: the same report with `=0`, and no side recorded.
#[test]
fn exact_blastp_records_no_bucket_fill_side() {
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
    let off = blastp(&args, &[OFF]);
    assert_eq!(off.code, Some(0), "{}", off.stderr);
    assert!(off.stdout == on.stdout, "report differs");
    assert_eq!(off.reached, (false, false), "=0 records a side");
}
