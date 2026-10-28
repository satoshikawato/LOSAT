//! The subject search's scratches kept from one query batch to the next of a blastn search
//! (`LOSAT_X_BNSCRATCHKEEP`, `src/algorithm/blastn/blast_engine/x_scratch_keep.rs`) on whole
//! searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses between two LOSAT code paths
//! with the same scratch: each query batch makes its threads' scratches and frees them, or the
//! search keeps them and the next batch takes them renewed. So the property under test is that the
//! report does not depend on it, and that kept scratches really serve the later batches: with
//! `LOSAT_TIMING`, `[TIMING] x_scratch_keep` counts `renewed=` scratches (taken from the keep) and
//! `fresh=` ones (made new), and with `LOSAT_X_BNSCRATCHKEEP=0` none is renewed. That the bytes
//! equal NCBI's is the business of the blastn fixtures.
//!
//! The query file holds four pieces of 20,000 letters, one from the middle of each of four genomes
//! (LC738868, LC738869, LC738870, LC738873), and `BATCH_SIZE=30000` puts each in a batch of its
//! own (a batch holds whole records): four batches per search, of which the last three take the
//! first one's scratch at one thread. Every task (megablast, blastn, dc-megablast, blastn-short)
//! searches LC738868 (`-subject`) at one and four threads and with two `LOSAT_X_QBATCHPAR` slots
//! (whose batches take their scratches from the same keep), and megablast and blastn search a
//! database of three genomes made by `losat makedb` (`-db`, whose drivers make one scratch per
//! worker of a block).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const OFF: (&str, &str) = ("LOSAT_X_BNSCRATCHKEEP", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const BATCH: (&str, &str) = ("BATCH_SIZE", "30000");
const QBATCHPAR: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
const SLOTS2: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");

const KEEP: &str = "[TIMING] x_scratch_keep: ";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of a one-record FASTA fixture.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .expect("read fixture")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// A file of this test process in the temporary directory.
fn temp_file(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("losat_scratchkeep_{name}_{}", std::process::id()))
}

/// Four pieces of 20,000 letters, one from the middle of each genome, in one query file (written
/// once per test binary).
fn pieces() -> &'static Path {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        let mut text = String::new();
        for name in ["LC738868", "LC738869", "LC738870", "LC738873"] {
            let genome = letters(&format!("{name}.fasta"));
            let start = genome.len() / 2 - 10_000;
            text.push_str(&format!(
                ">{name}_middle\n{}\n",
                &genome[start..start + 20_000]
            ));
        }
        let path = temp_file("pieces.fasta");
        std::fs::write(&path, text).expect("write the query file");
        path
    })
}

/// The binary with no `LOSAT_*` variable that chooses a code path and no batch-size variable.
/// `LOSAT_BLASTX_WORKER_LOG` stays, or is set to a file of this test when the environment has none:
/// a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI)
/// needs that file name as soon as a search has more than one thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env("LOSAT_BLASTX_WORKER_LOG", temp_file("worker.log"));
    }
    command
}

/// A database of three genomes made by `losat makedb` (once per test binary).
fn database() -> &'static str {
    static DB: OnceLock<String> = OnceLock::new();
    DB.get_or_init(|| {
        let set = temp_file("set.fna");
        let mut bytes = Vec::new();
        for name in ["LC738868", "LC738869", "LC738870"] {
            bytes.extend(std::fs::read(fasta(&format!("{name}.fasta"))).unwrap());
        }
        std::fs::write(&set, bytes).unwrap();
        let db = temp_file("set").to_str().unwrap().to_string();
        let output = losat()
            .arg("makedb")
            .arg("-in")
            .arg(&set)
            .args(["-dbtype", "nucl", "-out", &db])
            .output()
            .expect("run losat makedb");
        assert!(
            output.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        db
    })
}

/// One search of the pieces in batches of one record, `-outfmt 6` on standard output, against
/// LC738868 (`-subject`) or the database (`-db`), with the given variables as its only `LOSAT_*`
/// variables. Returns the report and the standard error.
fn blastn(task: &str, db: bool, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = losat();
    command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(pieces());
    if db {
        command.args(["-db", database()]);
    } else {
        command.arg("-subject").arg(fasta("LC738868.fasta"));
    }
    let output = command
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .env(BATCH.0, BATCH.1)
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} db {db} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The value of `name` (`name=VALUE`) in the `[TIMING] x_scratch_keep` line of `stderr`.
fn field(stderr: &str, name: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(KEEP))
        .unwrap_or_else(|| panic!("no {KEEP} line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|word| word.strip_prefix(name))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

/// For one task and subject source: the same report with the switch off (one thread, where every
/// batch makes its scratch) and by default at one and four threads, and with two slots where the
/// subjects are a `-subject` file; at one thread every batch after the first takes the first one's
/// scratch.
fn check(task: &str, db: bool) {
    let (reference, stderr) = blastn(task, db, 1, &[OFF, TIMING]);
    assert!(!reference.is_empty(), "{task} db {db}: no hits");
    assert_eq!(field(&stderr, "renewed="), 0, "{task} db {db}: {stderr}");
    assert!(field(&stderr, "fresh=") >= 4, "{task} db {db}: {stderr}");

    let (report, stderr) = blastn(task, db, 1, &[TIMING]);
    assert!(report == reference, "{task} db {db} threads 1");
    assert_eq!(
        field(&stderr, "fresh="),
        1,
        "{task} db {db} threads 1: one scratch for the search:\n{stderr}"
    );
    assert!(
        field(&stderr, "renewed=") >= 3,
        "{task} db {db} threads 1: the later batches take it:\n{stderr}"
    );

    let (report, stderr) = blastn(task, db, 4, &[TIMING]);
    assert!(report == reference, "{task} db {db} threads 4");
    assert!(field(&stderr, "renewed=") >= 1, "{task} db {db}: {stderr}");

    if !db {
        let (report, stderr) = blastn(task, db, 4, &[QBATCHPAR, SLOTS2, TIMING]);
        assert!(report == reference, "{task} slots");
        assert!(field(&stderr, "renewed=") >= 1, "{task} slots: {stderr}");
    }
}

#[test]
fn megablast_reports_the_same_with_kept_scratches() {
    check("megablast", false);
    check("megablast", true);
}

#[test]
fn blastn_reports_the_same_with_kept_scratches() {
    check("blastn", false);
    check("blastn", true);
}

#[test]
fn dc_megablast_reports_the_same_with_kept_scratches() {
    check("dc-megablast", false);
}

#[test]
fn blastn_short_reports_the_same_with_kept_scratches() {
    check("blastn-short", false);
}
