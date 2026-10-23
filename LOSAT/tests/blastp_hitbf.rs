//! The branch-free two-hit test of BLASTP's per-hit seed loop (`LOSAT_X_BPHITBF`,
//! `src/algorithm/blastp/x_seed_fused.rs`, module header) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
//! ```c
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//!
//!                 /* If the difference is less than the wordsize (i.e. last
//!                    hit and this hit overlap), give up */
//!
//!                 if (diff < wordsize) {
//!                     continue;
//!                 }
//! ```
//! The switch is in the strict set (on unless it is `0`).
//! The per-hit loop takes these two outcomes without a branch where its range holds and hands the
//! others to the reference's statements, with the same entries and the same calls of the
//! extension (the module's argument), so the property under test is that the report does not
//! depend on it. That reports equal NCBI's is the business of the BLASTP fixtures. The searches set
//! `LOSAT_X_BPDIAG=0`, so that the fused loop tests the hits one after the other (the default where
//! the CPU has no AVX2: aarch64, WebAssembly), besides the default. The shadow
//! (`LOSAT_X_BPSEEDSHADOW`) compares the fused loop with the reference on every subject and counts
//! the subjects the branch-free loop scanned: all of them at a two-hit window above the word size,
//! none at a window equal to the word size, where the loop keeps its statements, and none of
//! blastp-fast's (its loops on the compressed table keep their tests).
//!
//! The queries are the first 30 proteins of AP027131, the subjects the first 150 proteins of the
//! related AP027133, so there are extensions, flag-1 entries and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The per-hit loop of the fused loop (`LOSAT_X_BPDIAG`'s batches off).
const PER_HIT: (&str, &str) = ("LOSAT_X_BPDIAG", "0");
const OFF: (&str, &str) = ("LOSAT_X_BPHITBF", "0");
const ON: (&str, &str) = ("LOSAT_X_BPHITBF", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_BPSEEDSHADOW", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_hitbf_{}_{}_{tag}",
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

    /// The number after `name=` in the `[X_SHADOW] <switch>` line, or `None` without that line.
    fn shadow(&self, switch: &str, name: &str) -> Option<u64> {
        let prefix = format!("[X_SHADOW] {switch} ");
        let line = self.stderr.lines().find(|line| line.starts_with(&prefix))?;
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

/// The arguments of one search: `-subject` or `-db`, task, outfmt, threads, window.
fn search(
    inputs: &Inputs,
    db: bool,
    task: &str,
    outfmt: &str,
    threads: &str,
    window: &str,
) -> Vec<String> {
    let target = if db {
        ["-db", inputs.db.to_str().unwrap()]
    } else {
        ["-subject", inputs.subject.to_str().unwrap()]
    };
    [
        "-query",
        inputs.query.to_str().unwrap(),
        target[0],
        target[1],
        "-task",
        task,
        "-outfmt",
        outfmt,
        "-num_threads",
        threads,
        "-window_size",
        window,
    ]
    .map(String::from)
    .to_vec()
}

#[test]
fn the_branch_free_test_gives_the_reports_of_the_loop_without_it() {
    let inputs = Inputs::new();
    for db in [false, true] {
        for (task, outfmt, threads, window) in [
            ("blastp", "6", "1", "40"),
            ("blastp", "7", "3", "40"),
            ("blastp", "6", "1", "200"),
            ("blastp-fast", "6", "1", "40"),
            ("blastp-fast", "7", "3", "200"),
        ] {
            let args = search(&inputs, db, task, outfmt, threads, window);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let what = format!("db={db} {task} outfmt={outfmt} threads={threads} window={window}");
            let reference = blastp(&args, &[PER_HIT, OFF]);
            assert_eq!(reference.code, Some(0), "{what}: {}", reference.stderr);
            assert!(!reference.stdout.is_empty(), "{what}: no report");
            // The per-hit loop with the test; the default (the AVX2 batch where the CPU has it, else
            // the per-hit loop).
            for (arm, run, against) in [
                ("on", blastp(&args, &[PER_HIT]), &reference),
                ("named", blastp(&args, &[PER_HIT, ON]), &reference),
                ("default", blastp(&args, &[]), &reference),
                ("default off", blastp(&args, &[OFF]), &reference),
            ] {
                assert_eq!(run.code, against.code, "{what} {arm}: {}", run.stderr);
                assert!(run.stdout == against.stdout, "{what} {arm}: report differs");
                assert_eq!(run.messages(), against.messages(), "{what} {arm}");
            }
        }
    }
}

#[test]
fn the_shadow_compares_the_branch_free_loop_on_every_subject_in_its_range() {
    let inputs = Inputs::new();
    // A two-hit window above the word size (3): every subject; a window equal to the word size, and
    // blastp-fast (the compressed table's loops): none.
    for (task, window, all) in [
        ("blastp", "40", true),
        ("blastp", "3", false),
        ("blastp-fast", "40", false),
    ] {
        let what = format!("{task} window={window}");
        let args = search(&inputs, false, task, "6", "1", window);
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastp(&args, &[PER_HIT, OFF]);
        let run = blastp(&args, &[PER_HIT, SHADOW]);
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        let compared = run
            .shadow("BPSEED", "subjects")
            .unwrap_or_else(|| panic!("{what}: no BPSEED line: {}", run.stderr));
        let branch_free = run
            .shadow("BPHITBF", "subjects_compared")
            .unwrap_or_else(|| panic!("{what}: no BPHITBF line: {}", run.stderr));
        assert!(compared >= 150, "{what}: {compared} subjects");
        assert_eq!(branch_free, if all { compared } else { 0 }, "{what}");
        // Without the switch the shadow compares the loop without the test and has no BPHITBF
        // line.
        let off = blastp(&args, &[PER_HIT, SHADOW, OFF]);
        assert_eq!(off.code, Some(0), "{what}: {}", off.stderr);
        assert!(off.stdout == reference.stdout, "{what}: off differs");
        assert_eq!(off.shadow("BPSEED", "subjects"), Some(compared));
        assert_eq!(off.shadow("BPHITBF", "subjects_compared"), None);
    }
}
