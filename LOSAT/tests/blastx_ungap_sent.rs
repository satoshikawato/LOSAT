// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! BLASTX's two-hit ungapped extension with sentinel-stopped walks in the seed loops
//! (`LOSAT_X_BXUNGAPSENT`, `src/algorithm/blastx/x_ungap_sent.rs`, module header) on whole
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
//! The switch is in the strict set (on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set).
//! The path gives the reference extension's result (the module's argument), so the property under
//! test is that the report does not depend on it, for subjects read from FASTA (`-subject`) and
//! from a database (`-db`, read in blocks and, with `LOSAT_X_BXDBSTREAM=0`, from the store), whose
//! buffers carry the sentinel before the residues in different ways, in the fused loop and in the
//! direct loop (`LOSAT_X_BXDIAG=0`, with and without `LOSAT_X_BXHITBF`). That reports equal NCBI's
//! is the business of the BLASTX fixtures. The shadow (`LOSAT_X_BXUNGAPSENTSHADOW`) runs both
//! extensions on every candidate the path takes and panics on a difference; its line must show
//! compared extensions. `-word_size 5` (the compressed table) and `-window_size 0` (one-hit) stay
//! off the path and print no line.
//!
//! The query is the first 30,000 bases of AP027131 (genetic code 4), the subjects the first 150
//! proteins of the related AP027133, so there are long walks and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_BXUNGAPSENT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BXUNGAPSENTSHADOW", "1");
const DIRECT: (&str, &str) = ("LOSAT_X_BXDIAG", "0");
const NO_HITBF: (&str, &str) = ("LOSAT_X_BXHITBF", "0");
const NO_STREAM: (&str, &str) = ("LOSAT_X_BXDBSTREAM", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_ungap_sent_{}_{}_{tag}",
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

/// The first `n` bases of the one record of a FASTA file, as FASTA text.
fn bases(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut lines = text.lines();
    let header = lines.next().expect("a header");
    let sequence: String = lines.flat_map(str::chars).take(n).collect();
    assert_eq!(sequence.len(), n, "{name} is shorter than {n} bases");
    format!("{header}\n{sequence}\n")
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
        let query = dir.join("query.fa");
        let subject = dir.join("subject.faa");
        std::fs::write(&query, bases("AP027131.fasta", 30_000)).unwrap();
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

    /// The number after `name=` in the `[X_SHADOW] BXUNGAPSENT` line, or `None` without that line.
    fn shadow(&self, name: &str) -> Option<u64> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with("[X_SHADOW] BXUNGAPSENT "))?;
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

fn blastx(args: &[&str], env: &[(&str, &str)]) -> Run {
    let output = losat(env)
        .arg("blastx")
        .args(args)
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// Default, `=0`, the shadow (fused and direct loops) and `LOSAT_X_DEFAULT=0` give the same
/// report, on `-subject` and `-db` at 1 and 2 threads; the shadow compares extensions.
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
                "-query_gencode",
                "4",
                "-outfmt",
                "6",
                "-evalue",
                "10",
                "-num_threads",
                threads,
            ];
            let on = blastx(&args, &[]);
            assert_eq!(on.code, Some(0), "{source} n{threads}: {}", on.stderr);
            assert!(!on.stdout.is_empty(), "{source} n{threads}: no hits");
            for (label, env) in [
                ("off", &[OFF][..]),
                ("shadow", &[SHADOW][..]),
                ("shadow direct", &[SHADOW, DIRECT][..]),
                ("shadow direct no-bf", &[SHADOW, DIRECT, NO_HITBF][..]),
                ("shadow store", &[SHADOW, NO_STREAM][..]),
                ("LOSAT_X_DEFAULT=0", &[DEFAULT_OFF][..]),
            ] {
                let run = blastx(&args, env);
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
                if label.starts_with("shadow") {
                    let compared = run.shadow("compared").expect("shadow line");
                    assert!(
                        compared > 0,
                        "{source} n{threads} {label}: nothing compared"
                    );
                    assert!(run.shadow("subjects").unwrap() > 0);
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

/// Word size 5 (the compressed table) and the one-hit extension (`-window_size 0`) stay off the
/// path: the same report, no shadow line.
#[test]
fn word_size_5_and_one_hit_stay_off_the_path() {
    let inputs = Inputs::new();
    let query = inputs.query.to_str().unwrap();
    let subject = inputs.subject.to_str().unwrap();
    for extra in [&["-word_size", "5"][..], &["-window_size", "0"][..]] {
        let mut args = vec![
            "-query",
            query,
            "-subject",
            subject,
            "-query_gencode",
            "4",
            "-outfmt",
            "6",
        ];
        args.extend_from_slice(extra);
        let on = blastx(&args, &[]);
        assert_eq!(on.code, Some(0), "{extra:?}: {}", on.stderr);
        let shadow = blastx(&args, &[SHADOW]);
        assert_eq!(shadow.code, Some(0), "{extra:?}: {}", shadow.stderr);
        assert!(shadow.stdout == on.stdout, "{extra:?}: report differs");
        assert_eq!(
            shadow.shadow("compared"),
            None,
            "{extra:?}: shadow line off the path"
        );
    }
}
