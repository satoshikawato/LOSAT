//! BLASTP's two-hit ungapped extension inlined into the fused seed loop with sentinel-stopped walks
//! (`LOSAT_X_BPUNGAPSENT`, `src/algorithm/blastp/x_ungapped_sent.rs`, module header) on whole
//! searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155
//! ```c
//! static Int4 s_BlastAaExtendTwoHit(...)
//! {
//!     ...
//!     left_score = s_BlastAaExtendLeft(matrix, subject, query, ...);
//!     ...
//!         right_score = s_BlastAaExtendRight(matrix, subject, query, ...);
//! ```
//! The switch is in the strict set (on unless it is `0`).
//! The path gives the extension's result of the closure it replaces (the module's argument), so
//! the property under test is that the report does not depend on it, for subjects read from FASTA
//! (`-subject`) and from a database (`-db`), whose buffers carry the sentinels in different ways.
//! That reports equal NCBI's is the business of the BLASTP fixtures. The shadow
//! (`LOSAT_X_BPUNGAPSENTSHADOW`) runs both extensions on every candidate the path takes and panics
//! on a difference; its line must show compared extensions. blastp-fast (word size 5, the
//! compressed table) stays on the closure without the path and prints no line.
//!
//! The queries are the first 30 proteins of AP027131 (a query block of 30 contexts), the subjects
//! the first 150 proteins of the related AP027133, so there are long walks and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_BPUNGAPSENT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BPUNGAPSENTSHADOW", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_ungap_sent_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The first `n` records of a FASTA file, as text.
fn head(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut out = String::new();
    let mut records = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            records += 1;
            if records > n {
                break;
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    assert_eq!(records, n + 1, "{name} has more than {n} records");
    out
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
    db: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("query.faa");
        let subject = dir.join("subject.faa");
        std::fs::write(&query, head("AP027131.faa", 30)).unwrap();
        std::fs::write(&subject, head("AP027133.faa", 150)).unwrap();
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
            query,
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
}

impl Run {
    /// Standard error without the diagnostic lines (`[TIMING]`, `[X_SHADOW]`).
    fn messages(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[TIMING]") && !line.starts_with("[X_SHADOW]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// The number after `name=` in the `[X_SHADOW] BPUNGAPSENT` line, or `None` without that line.
    fn shadow(&self, name: &str) -> Option<u64> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with("[X_SHADOW] BPUNGAPSENT "))?;
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

fn blastp(args: &[&str], env: &[(&str, &str)]) -> Run {
    let output = losat(env)
        .arg("blastp")
        .args(args)
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// Default, `=0` and the shadow give the same report, on `-subject` and `-db`
/// at 1 and 2 threads; the shadow compares extensions.
#[test]
fn reports_do_not_depend_on_the_sentinel_extension() {
    let inputs = Inputs::new();
    let query = inputs.query.to_str().unwrap();
    let subject = inputs.subject.to_str().unwrap();
    let db = inputs.db.to_str().unwrap();
    for (source, target) in [("-subject", subject), ("-db", db)] {
        for threads in ["1", "2"] {
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
            assert_eq!(on.code, Some(0), "{source} n{threads}: {}", on.stderr);
            assert!(!on.stdout.is_empty(), "{source} n{threads}: no hits");
            for (label, env) in [("off", &[OFF][..]), ("shadow", &[SHADOW][..])] {
                let run = blastp(&args, env);
                assert_eq!(
                    run.code,
                    Some(0),
                    "{source} n{threads} {label}: {}",
                    run.stderr
                );
                assert!(
                    run.stdout == on.stdout,
                    "{source} n{threads} {label}: report differs"
                );
                assert_eq!(run.messages(), on.messages(), "{source} n{threads} {label}");
                if label == "shadow" {
                    let compared = run.shadow("compared").expect("shadow line");
                    assert!(compared > 0, "{source} n{threads}: nothing compared");
                    assert!(run.shadow("blocks").unwrap() > 0);
                } else {
                    assert_eq!(
                        run.shadow("compared"),
                        None,
                        "{label} printed a shadow line"
                    );
                }
            }
        }
    }
}

/// blastp-fast (word size 5) keeps the closure without the path: the same report, no shadow line.
#[test]
fn blastp_fast_keeps_the_closure_without_the_path() {
    let inputs = Inputs::new();
    let args = [
        "-query",
        inputs.query.to_str().unwrap(),
        "-subject",
        inputs.subject.to_str().unwrap(),
        "-task",
        "blastp-fast",
        "-outfmt",
        "6",
    ];
    let on = blastp(&args, &[]);
    assert_eq!(on.code, Some(0), "{}", on.stderr);
    let shadow = blastp(&args, &[SHADOW]);
    assert_eq!(shadow.code, Some(0), "{}", shadow.stderr);
    assert!(shadow.stdout == on.stdout, "report differs");
    assert_eq!(
        shadow.shadow("compared"),
        None,
        "shadow line on blastp-fast"
    );
}
