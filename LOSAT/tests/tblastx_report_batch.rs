//! TBLASTX's final order made for each query batch when it ends (`LOSAT_X_TBXRPTBATCH`,
//! `src/algorithm/tblastx/x_rpt_batch.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_TBXRPTBATCH`
//! is `0`. With it on, the batch loop puts each batch's HSPs in report order when the batch ends
//! and the write after the last batch keeps that order; off, the write orders every batch's HSPs
//! at once. The property under test is that the reports, the messages and the exit status do not
//! depend on the switch: searches of several query batches (`BATCH_SIZE` makes every query its own
//! batch; without it the short queries share a batch) are compared with the switch off, with the
//! options that act on the order (`-max_target_seqs`, `-max_hsps`, `-culling_limit`,
//! `-qcov_hsp_perc`), the report formats, `-subject` and `-db`, at 1 and 4 threads.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const OFF: (&str, &str) = ("LOSAT_X_TBXRPTBATCH", "0");

fn letters(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    std::fs::read(path)
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// This test's temporary directory.
fn dir() -> &'static Path {
    static DIR: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    DIR.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("losat-tbxrptbatch-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temporary directory");
        dir
    })
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(&str, &[u8])]) -> PathBuf {
    let path = dir().join(file);
    let mut text = Vec::new();
    for (title, seq) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).expect("write FASTA");
    path
}

/// Five queries from MjeNMV (a record without letters among them, which has a warning before its
/// report; two short ones that share a batch without `BATCH_SIZE`), written once (the tests run
/// in parallel).
fn queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let q = letters("MjeNMV.fasta");
        write_fasta(
            "queries.fa",
            &[
                ("q1 first part", &q[..12_000]),
                ("empty record", b""),
                ("q2", &q[100_000..103_000]),
                ("q3", &q[150_000..154_000]),
                ("q4 last part", &q[200_000..215_000]),
            ],
        )
    })
}

/// The MelaMJNV and PemoMJNVA genomes, as a FASTA file and as a database made by `losat makedb`,
/// made once.
fn subjects() -> &'static (PathBuf, String) {
    static SUBJECTS: std::sync::OnceLock<(PathBuf, String)> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let fasta = write_fasta(
            "subjects.fa",
            &[
                ("mela MelaMJNV", &letters("MelaMJNV.fasta")),
                ("pemo PemoMJNVA", &letters("PemoMJNVA.fasta")),
            ],
        );
        let db = dir().join("subjects").to_str().expect("path").to_string();
        let output = losat()
            .arg("makedb")
            .arg("-in")
            .arg(&fasta)
            .args(["-dbtype", "nucl", "-out", &db])
            .output()
            .expect("run losat makedb");
        assert!(
            output.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        (fasta, db)
    })
}

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature needs it as soon as a search has
/// more than one thread, and it selects nothing) and without `BATCH_SIZE`.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env_remove("BATCH_SIZE");
    command
}

/// One search; `db` chooses `-db` over `-subject`.
fn tblastx(db: bool, threads: u32, args: &[&str], env: &[(&str, &str)]) -> Output {
    let (fasta, database) = subjects();
    let mut command = losat();
    command.envs(env.iter().copied()).arg("tblastx").arg("-query").arg(queries());
    if db {
        command.args(["-db", database]);
    } else {
        command.arg("-subject").arg(fasta);
    }
    command
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn check(args: &[&str]) {
    for db in [false, true] {
        for batch_size in [None, Some("4000")] {
            let env = |off: bool| -> Vec<(&str, &str)> {
                let mut env = Vec::new();
                if off {
                    env.push(OFF);
                }
                if let Some(size) = batch_size {
                    env.push(("BATCH_SIZE", size));
                }
                env
            };
            let off = tblastx(db, 1, args, &env(true));
            assert!(
                off.status.success(),
                "{args:?} db {db}: {}",
                String::from_utf8_lossy(&off.stderr)
            );
            assert!(!off.stdout.is_empty(), "{args:?} db {db}: the queries have hits");
            for threads in [1, 4] {
                let on = tblastx(db, threads, args, &env(false));
                let case = format!("{args:?} db {db} BATCH_SIZE {batch_size:?} threads {threads}");
                assert_eq!(on.status.code(), off.status.code(), "{case}");
                assert_eq!(
                    String::from_utf8_lossy(&on.stderr),
                    String::from_utf8_lossy(&off.stderr),
                    "{case}"
                );
                assert!(on.stdout == off.stdout, "{case}: reports differ");
            }
        }
    }
}

#[test]
fn tabular_reports_do_not_depend_on_the_switch() {
    check(&["-outfmt", "6"]);
    check(&["-outfmt", "7"]);
}

#[test]
fn hit_list_options_do_not_depend_on_the_switch() {
    check(&["-outfmt", "6", "-max_target_seqs", "1", "-max_hsps", "3"]);
    check(&["-outfmt", "6", "-qcov_hsp_perc", "2"]);
    check(&["-outfmt", "6", "-culling_limit", "1"]);
}

#[test]
fn pairwise_and_structured_reports_do_not_depend_on_the_switch() {
    check(&["-outfmt", "0", "-max_target_seqs", "1"]);
    check(&["-outfmt", "5"]);
    check(&["-outfmt", "11"]);
}
