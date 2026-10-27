//! blastp-fast's block search with `LOSAT_X_BPQSLOTS` (`src/algorithm/blastp/x_query_slots.rs`):
//! several query blocks searched at a time on the search's pool, their reports written in block
//! order.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): every run here names it, `0` (the block
//! loop) or `1` (two blocks in flight), and states `LOSAT_X_BPQBLOCK` itself. These tests compare
//! the reports, standard error (without the `[TIMING]` and `[X_STATS]` lines) and exit codes with
//! the switch off and on, with `-subject` and `-db`, for the reports the slots take, and read the
//! `[X_STATS] BPQSLOTS` line (`LOSAT_X_STATS`, printed when a count is not zero): `searches` (block
//! searches that took the slots) and `blocks` (their blocks). That the reports come in block order
//! whatever order the blocks finish in, and the failure of a block, are unit tests of
//! `x_query_slots.rs`. Each test writes its files to a directory of its own.
//!
//! The queries are records of `MeenMJNV.faa` against `MejoMJNV.faa` (two related nimavirus
//! proteomes): the four records BDT62581.1 to BDT62584.1 (2,539 letters, a batch each at
//! `BATCH_SIZE` 250, four blocks at `LOSAT_X_BPQBLOCK=1`; a run costs 0.05-0.45 CPU s, one on all
//! 101 records 2-8).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!     ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//! ```

use std::collections::HashMap;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of the test `tag` (removed first).
fn test_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpqslots-{tag}-{}", std::process::id()));
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The queries of a search (a file of records of `MeenMJNV.faa`) and its `BATCH_SIZE`.
struct Queries {
    file: PathBuf,
    batch_size: &'static str,
}

/// The four records BDT62581.1 to BDT62584.1 of `MeenMJNV.faa`, a batch each at `BATCH_SIZE` 250.
fn four_queries(dir: &Path) -> Queries {
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let records: Vec<String> = text
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| format!(">{record}"))
        .collect();
    let position = |accession: &str| {
        records
            .iter()
            .position(|record| record.starts_with(&format!(">{accession} ")))
            .expect("the record is in the fixture")
    };
    let file = dir.join("four.faa");
    std::fs::write(
        &file,
        records[position("BDT62581.1")..=position("BDT62584.1")].concat(),
    )
    .unwrap();
    Queries {
        file,
        batch_size: "250",
    }
}

/// What one run printed and how it ended.
#[derive(Debug)]
struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

impl Run {
    /// Standard error without the `[TIMING]` and `[X_STATS]` lines.
    fn stderr_compared(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[TIMING]") && !line.starts_with("[X_STATS]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// The counters of the `[X_STATS] BPQSLOTS` line: empty when the run printed none.
    fn slots(&self) -> HashMap<String, u64> {
        let prefix = "[X_STATS] BPQSLOTS ";
        let Some(line) = self.stderr.lines().find(|line| line.starts_with(prefix)) else {
            return HashMap::new();
        };
        line[prefix.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
    }

    /// The searches of query sets that ran to their reports (`[TIMING] blastp_total` lines): the
    /// blocks of a block search.
    fn searches(&self) -> u64 {
        self.stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] blastp_total:"))
            .count() as u64
    }
}

/// One blastp run (`task`) of `queries` against `subject` (`-subject FILE` or `-db NAME`) with
/// `slots` as `LOSAT_X_BPQSLOTS` (`0` or `1`), `LOSAT_X_BPQBLOCK=1`, the queries' `BATCH_SIZE`, `LOSAT_TIMING`,
/// `LOSAT_X_STATS` and the arguments `extra`; no other `LOSAT_*` variable that chooses a path.
/// `LOSAT_BLASTX_WORKER_LOG` stays (a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread); without the caller's value a file in `dir` is named.
fn blastp(
    dir: &Path,
    task: &str,
    queries: &Queries,
    subject: &[&OsStr],
    threads: u32,
    slots: &str,
    extra: &[&str],
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
    }
    command
        .arg("blastp")
        .args(["-task", task])
        .arg("-query")
        .arg(&queries.file)
        .args(subject)
        .args(["-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra);
    if !extra.contains(&"-outfmt") {
        command.args(["-outfmt", "6"]);
    }
    let output = command
        .env("BATCH_SIZE", queries.batch_size)
        .env("LOSAT_X_BPQBLOCK", "1")
        .env("LOSAT_X_BPQSLOTS", slots)
        .env("LOSAT_TIMING", "1")
        .env("LOSAT_X_STATS", "1")
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// The switch off and on (`slots`) give the same report, standard error and exit code; returns the
/// run with the switch on.
#[allow(clippy::too_many_arguments)]
fn same_runs(
    dir: &Path,
    task: &str,
    queries: &Queries,
    subject: &[&OsStr],
    threads: u32,
    slots: &str,
    extra: &[&str],
) -> Run {
    let off = blastp(dir, task, queries, subject, threads, "0", extra);
    let on = blastp(dir, task, queries, subject, threads, slots, extra);
    let what = format!("{task} slots {slots} threads {threads} {extra:?}");
    assert!(on.stdout == off.stdout, "{what}: the report");
    assert_eq!(
        on.stderr_compared(),
        off.stderr_compared(),
        "{what}: standard error"
    );
    assert_eq!(on.code, off.code, "{what}: exit code");
    assert!(
        off.slots().is_empty(),
        "{what}: the block loop counts nothing"
    );
    on
}

/// The `-db` of `MejoMJNV.faa` in `dir`.
fn database(dir: &Path) -> PathBuf {
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

// The reports the slots take, `-subject`: outfmt 6, 7 (its count of the queries in the last block),
// 0 (the prolog with the first block, the epilog with the last), 5 (the iteration numbers), 8, 9 and
// 12 (a Seq-annot per query); at
// four threads and at two (the pool's least). Every block of the search counted.
#[test]
fn block_reports_with_the_slots_are_the_block_loops() {
    let dir = test_dir("subject");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let four = four_queries(&dir);
    for outfmt in ["6", "7", "0", "5", "8", "9", "12"] {
        for threads in [4, 2] {
            let on = same_runs(
                &dir,
                "blastp-fast",
                &four,
                &subject,
                threads,
                "1",
                &["-outfmt", outfmt],
            );
            assert!(!on.stdout.is_empty(), "the fixture pair has a report");
            assert_eq!(on.searches(), 4, "outfmt {outfmt}: four blocks");
            let counts = on.slots();
            assert_eq!(
                counts.get("searches"),
                Some(&1),
                "outfmt {outfmt}: {counts:?}"
            );
            assert_eq!(
                counts.get("blocks"),
                Some(&4),
                "outfmt {outfmt}: {counts:?}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

// `-db`: the database search in blocks; outfmt 6, 0, 7 and 5 at two and four threads.
#[test]
fn database_block_reports_with_the_slots_are_the_block_loops() {
    let dir = test_dir("db");
    let name = database(&dir);
    let subject = [OsStr::new("-db"), name.as_os_str()];
    let four = four_queries(&dir);
    for (outfmt, threads) in [("6", 2), ("0", 4), ("7", 4), ("5", 2)] {
        let on = same_runs(
            &dir,
            "blastp-fast",
            &four,
            &subject,
            threads,
            "1",
            &["-outfmt", outfmt],
        );
        assert!(!on.stdout.is_empty(), "the fixture pair has a report");
        assert_eq!(on.slots().get("blocks"), Some(&4), "outfmt {outfmt}");
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Where the slots do not act, the switch on counts nothing and the reports are the block loop's: one
// thread, exact blastp, outfmt 11 (the archive's state moves from document to document). A report
// that cannot be written (`-out /dev/full`) fails alike: the same message and exit code.
#[test]
fn the_slots_do_not_act_where_the_block_loop_is_kept() {
    let dir = test_dir("kept");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let four = four_queries(&dir);
    for (task, threads, extra) in [
        ("blastp-fast", 1, &[][..]),
        ("blastp", 4, &[][..]),
        ("blastp-fast", 4, &["-outfmt", "11"][..]),
    ] {
        let on = same_runs(&dir, task, &four, &subject, threads, "1", extra);
        assert_eq!(on.searches(), 4, "{task} threads {threads} {extra:?}");
        assert!(
            on.slots().is_empty(),
            "{task} threads {threads} {extra:?}: {:?}",
            on.slots()
        );
    }
    if Path::new("/dev/full").exists() {
        let on = same_runs(
            &dir,
            "blastp-fast",
            &four,
            &subject,
            4,
            "1",
            &["-out", "/dev/full"],
        );
        assert_ne!(on.code, Some(0), "a report that cannot be written fails");
    }
    std::fs::remove_dir_all(&dir).ok();
}
