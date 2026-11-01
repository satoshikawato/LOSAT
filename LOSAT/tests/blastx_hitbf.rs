//! The branch-free two-hit test of BLASTX's direct seed loop (`LOSAT_X_BXHITBF`,
//! `src/algorithm/blastx/x_seed_scan.rs`, module header) on whole searches.
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
//! The direct loop takes these two outcomes without a branch where its range holds and hands the
//! others to the statements of the reference word finder, with the same cell states and the same
//! init HSPs (the module's argument), so the property under test is that the report does not
//! depend on it. That reports equal NCBI's is the business of the BLASTX fixtures. The searches set
//! `LOSAT_X_BXDIAG=0`, so that the direct loop scans every subject (the default where the CPU has no
//! AVX2: aarch64, WebAssembly), besides the default. At the one-hit window and at a window equal to
//! the word size the loop keeps its statements, so the reports must not depend on the switch there
//! either.
//!
//! The queries are 6 pieces of AP027131 of 1,500 to 3,900 letters each; the subjects are the first
//! 150 proteins of AP027131, so there are extensions, flag-1 cells and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The direct loop on every subject (`LOSAT_X_BXDIAG`'s fused loop off).
const DIRECT: (&str, &str) = ("LOSAT_X_BXDIAG", "0");
const OFF: (&str, &str) = ("LOSAT_X_BXHITBF", "0");
const ON: (&str, &str) = ("LOSAT_X_BXHITBF", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_hitbf_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The records of a FASTA file as (title, letters).
fn records(name: &str) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(title) = line.strip_prefix('>') {
            records.push((title.to_string(), Vec::new()));
        } else if let Some(record) = records.last_mut() {
            record.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    records
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    let mut text = Vec::new();
    for (title, letters) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in letters.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).expect("write FASTA");
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
        let genome = &records("AP027131.fasta")[0].1;
        let queries: Vec<(String, Vec<u8>)> = (0..6)
            .map(|index| {
                let start = index * 80_000;
                let length = 1_500 + (index * 479) % 2_400;
                (
                    format!("q{index} piece {start}"),
                    genome[start..start + length].to_vec(),
                )
            })
            .collect();
        let mut proteins = records("AP027131.faa");
        proteins.truncate(150);
        let query = dir.join("query.fna");
        let subject = dir.join("subject.faa");
        write_fasta(&query, &queries);
        write_fasta(&subject, &proteins);
        let db = dir.join("db");
        let made = losat(&[])
            .arg("makeblastdb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run losat makeblastdb");
        assert!(
            made.status.success(),
            "makeblastdb: {}",
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

/// The arguments of one search: `-subject` or `-db`, outfmt, threads, window.
fn search(inputs: &Inputs, db: bool, outfmt: &str, threads: &str, window: &str) -> Vec<String> {
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
        for (outfmt, threads, window) in [("6", "1", "40"), ("7", "3", "40"), ("6", "1", "200")] {
            let args = search(&inputs, db, outfmt, threads, window);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let what = format!("db={db} outfmt={outfmt} threads={threads} window={window}");
            let reference = blastx(&args, &[DIRECT, OFF]);
            assert_eq!(reference.code, Some(0), "{what}: {}", reference.stderr);
            assert!(!reference.stdout.is_empty(), "{what}: no report");
            // The direct loop with the test; the default (the fused loop where it pays, the direct
            // loop elsewhere).
            for (arm, run, against) in [
                ("on", blastx(&args, &[DIRECT]), &reference),
                ("named", blastx(&args, &[DIRECT, ON]), &reference),
                ("default", blastx(&args, &[]), &reference),
            ] {
                assert_eq!(run.code, against.code, "{what} {arm}: {}", run.stderr);
                assert!(run.stdout == against.stdout, "{what} {arm}: report differs");
                assert_eq!(run.messages(), against.messages(), "{what} {arm}");
            }
        }
    }
}

#[test]
fn the_reports_do_not_depend_on_the_switch_where_the_loop_keeps_its_statements() {
    let inputs = Inputs::new();
    // A window equal to the word size (3) and the one-hit window: the branch-free test does not act.
    for window in ["3", "0"] {
        let args = search(&inputs, false, "6", "1", window);
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastx(&args, &[DIRECT, OFF]);
        assert_eq!(
            reference.code,
            Some(0),
            "window={window}: {}",
            reference.stderr
        );
        let run = blastx(&args, &[DIRECT, ON]);
        assert_eq!(run.code, Some(0), "window={window}: {}", run.stderr);
        assert!(
            run.stdout == reference.stdout,
            "window={window}: report differs"
        );
        assert_eq!(run.messages(), reference.messages(), "window={window}");
    }
}
