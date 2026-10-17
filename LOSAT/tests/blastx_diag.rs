//! BLASTX's two-hit seed stage with the subject scan and the diagonal test fused, and the hits of
//! one lookup cell tested in batches of eight (`LOSAT_X_BXDIAG`, `seed/x_seed_diag.rs`), on whole
//! searches.
//!
//! The switch is off unless set to any value (it is not in the strict set); its value `per-hit`
//! or `scalar` forces a cell test, and `LOSAT_X_BXDIAGSHADOW` runs the reference beside it for
//! every subject and stops at a difference. So the properties under test are that the report,
//! standard error without the diagnostic lines and the exit status equal those of the search
//! without the switch, with every cell test, with `-subject` and `-db`, outfmt 6 and 7, one and
//! three threads and a two-hit window other than the default; that the shadow compared subjects;
//! and that the timing line names the cell test that ran. That the bytes equal NCBI's is the
//! business of the BLASTX fixtures.
//!
//! The queries are 12 pieces of AP027131 of 1,500 to 3,900 letters each; the subjects are the
//! first 150 proteins of AP027131, so the pieces have many hits and the lookup cells many
//! entries.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_diag_{}_{}_{tag}",
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
        let queries: Vec<(String, Vec<u8>)> = (0..12)
            .map(|index| {
                let start = index * 40_000;
                let length = 1_500 + (index * 373) % 2_400;
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// The searches of the test: (`-subject` or `-db`, outfmt, threads, window).
fn searches(inputs: &Inputs) -> Vec<Vec<String>> {
    let mut searches = Vec::new();
    for target in [
        ["-subject", inputs.subject.to_str().unwrap()],
        ["-db", inputs.db.to_str().unwrap()],
    ] {
        for (outfmt, threads, window) in [("6", "1", "40"), ("7", "3", "40"), ("6", "1", "200")] {
            searches.push(
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
                .to_vec(),
            );
        }
    }
    searches
}

#[test]
fn the_fused_seed_loop_gives_the_reports_of_the_reference() {
    let inputs = Inputs::new();
    for args in searches(&inputs) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastx(&args, &[]);
        assert_eq!(reference.code, Some(0), "{args:?}: {}", reference.stderr);
        assert!(
            reference.stdout.split(|&b| b == b'\n').count() > 50,
            "{args:?}: too few rows"
        );
        for value in ["1", "per-hit", "scalar"] {
            let run = blastx(&args, &[("LOSAT_X_BXDIAG", value)]);
            let what = format!("{args:?} LOSAT_X_BXDIAG={value}");
            assert_eq!(run.code, reference.code, "{what}: {}", run.stderr);
            assert!(run.stdout == reference.stdout, "{what}: report differs");
            assert_eq!(run.messages(), reference.messages(), "{what}");
        }
    }
}

#[test]
fn the_shadow_compares_every_subject_and_the_timing_line_names_the_cell_test() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[]);
    for value in ["1", "per-hit", "scalar"] {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXDIAG", value),
                ("LOSAT_X_BXDIAGSHADOW", "1"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{value}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{value}: report differs");
        // Every query batch scans every subject; the shadow compares each scan.
        let compared = run
            .stderr
            .lines()
            .find_map(|line| line.strip_prefix("[X_SHADOW] BXDIAG subjects_compared="))
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("{value}: no shadow line: {}", run.stderr));
        assert!(compared >= 150, "{value}: {compared} subjects compared");
        let line = run
            .stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] blastx_diag: "))
            .unwrap_or_else(|| panic!("{value}: no timing line: {}", run.stderr));
        let mode = line
            .split_whitespace()
            .find_map(|field| field.strip_prefix("mode="))
            .unwrap();
        match value {
            "per-hit" => assert_eq!(mode, "per-hit", "{line}"),
            "scalar" => assert_eq!(mode, "scalar", "{line}"),
            _ => {
                #[cfg(target_arch = "x86_64")]
                let expected = if std::arch::is_x86_feature_detected!("avx2") {
                    "avx2"
                } else {
                    "per-hit"
                };
                #[cfg(not(target_arch = "x86_64"))]
                let expected = "per-hit";
                assert_eq!(mode, expected, "{line}");
            }
        }
        if mode != "per-hit" {
            assert!(!line.contains(" cells=0 "), "{line}");
        }
    }
}
