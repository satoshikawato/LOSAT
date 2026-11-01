//! blastp-fast's lookup build with `LOSAT_X_BPLKPFLAT` (`src/algorithm/tblastx/lookup/x_clut_par.rs`,
//! `src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`) with `LOSAT_X_BPCLUT`: it is on by default,
//! `LOSAT_X_BPLKPFLAT=0` turns it off. It acts
//! only where the search's pool builds the flat compressed table of `LOSAT_X_BPCLUT` (on by default as
//! well, `LOSAT_X_BPCLUT=0` turns it off), with more than one thread. With it on, the build's two passes run with one latch
//! each, a block search's table takes the vectors of the previous block's table, and with
//! `LOSAT_X_BPQBLOCK` the contexts' Karlin blocks are tasks of the build. These tests compare the
//! reports with the switch off (`LOSAT_X_BPLKPFLAT=0`) and on, with `-subject` and `-db`, and read
//! the `[X_STATS] BPLKPFLAT tables= recycled= recycled_hits=` line of the run with the switch on
//! (printed with `LOSAT_X_STATS`, when the pool built a table), which shows that the path ran: the
//! tables built and the builds that took the vectors of an earlier table. Every run but the
//! default-set test states `LOSAT_X_BPQBLOCK` and `LOSAT_X_BPCLUT` itself, so that the tests do not
//! depend on the defaults of the strict set.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
//! ```c
//!     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//!                                     query, locations);
//!     s_CompressedLookupFinalize(lookup);
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

/// What one run printed: the report and standard error.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
}

impl Run {
    /// The counters of the `[X_STATS] <name>` line (`LOSAT_X_STATS`; printed only when a counter is
    /// not zero): empty when the run printed none.
    fn stats(&self, name: &str) -> HashMap<String, u64> {
        let prefix = format!("[X_STATS] {name} ");
        let Some(line) = self.stderr.lines().find(|line| line.starts_with(&prefix)) else {
            return HashMap::new();
        };
        line[prefix.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
    }

    /// `stats(name)` where the line is expected: fails with the run's stderr when it is missing.
    fn printed(&self, name: &str) -> HashMap<String, u64> {
        let counts = self.stats(name);
        assert!(
            !counts.is_empty(),
            "no [X_STATS] {name} line; stderr:\n{}",
            self.stderr
        );
        counts
    }

    /// The searches of query sets that ran to their reports (`[TIMING] blastp_total` lines, with
    /// `LOSAT_TIMING`): the blocks of a block search; a query chunk of a split batch ends after its
    /// preliminary stage and prints none.
    fn searches(&self) -> u64 {
        self.stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] blastp_total:"))
            .count() as u64
    }
}

/// The `recycled_hits` counter: builds whose `hits` are an earlier table's buffer, which happens only
/// where that buffer has the room; at least one, and at most `recycled`.
fn assert_hits_recycled(counts: &HashMap<String, u64>, what: &str) {
    assert!(
        (1..=counts["recycled"]).contains(&counts["recycled_hits"]),
        "{what}: {counts:?}"
    );
}

/// One blastp-fast run of the fixture pair (two related nimavirus proteomes) against `subject`
/// (`-subject FILE` or `-db NAME`) with `vars` as its only `LOSAT_*` variables that choose a code
/// path, `BATCH_SIZE` 3000 (several query batches) and `-outfmt outfmt`. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects nothing;
/// without the caller's value a file of the tests is named. The run must exit with 0.
fn blastp_fast(subject: &[&OsStr], threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
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
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bplkpflat-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .args(["-task", "blastp-fast"])
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
        .args(subject)
        .args(["-evalue", "0.001", "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "3000")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    Run {
        stdout: output.stdout,
        stderr,
    }
}

fn subject_file() -> PathBuf {
    fixture("MejoMJNV.faa")
}

/// A protein database of the subject proteome, made with `losat makeblastdb` in a directory of the
/// test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bplkpflat-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(subject_file())
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

const ON: (&str, &str) = ("LOSAT_X_BPLKPFLAT", "1");
const OFF: (&str, &str) = ("LOSAT_X_BPLKPFLAT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

/// The switch off and on give the same report; returns the run with the switch on (the callers read
/// its `[X_STATS] BPLKPFLAT` line: `printed` where it is expected, `stats` where it may be
/// missing).
fn same_reports(subject: &[&OsStr], threads: u32, base: &[(&str, &str)], outfmt: &str) -> Run {
    let off = blastp_fast(subject, threads, &with(base, &[OFF]), outfmt);
    assert!(!off.stdout.is_empty(), "the fixture pair has rows");
    let on = blastp_fast(subject, threads, &with(base, &[ON, STATS]), outfmt);
    assert!(
        on.stdout == off.stdout,
        "{base:?} threads {threads} outfmt {outfmt}: the report differs with the switch"
    );
    on
}

// The block search (every block a search with its own table) at 2 threads, with blocks of one
// batch: every build is the pool's and every build after a block's first takes the vectors of the
// table before it. The thread count only sets how the build is cut into
// chunks, stripes and tasks (the unit tests of `x_clut_par.rs` build at 1, 2, 3 and 7 threads); the
// other tests where the pool builds run at 4 threads. Blocks of 20,000 letters run in the database
// test, which takes the same lookup region per block.
#[test]
fn block_search_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", "1")];
    let counts = same_reports(&subject, 2, &base, "6").printed("BPLKPFLAT");
    let what = format!("{counts:?}");
    assert!(counts["tables"] > 1, "{what}");
    assert_eq!(counts["recycled"], counts["tables"] - 1, "{what}");
    assert_hits_recycled(&counts, &what);
}

// The search of all queries: one build, by the pool; no earlier table. `LOSAT_X_BPQBLOCK=0`, and a
// block larger than the 101 queries (77,448 letters, above the built-in 50,000), which is one block
// with the block switch on (`x_query_block_plan` takes the search of all queries for it).
// The pairwise report is written in blocks, so only the block size makes the search of all queries.
#[test]
fn search_of_all_queries_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for letters in ["100000000", "0"] {
        let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", letters)];
        let counts = same_reports(&subject, 4, &base, "0").printed("BPLKPFLAT");
        let what = format!("blocks {letters}: {counts:?}");
        assert_eq!(counts["tables"], 1, "{what}");
        assert_eq!(counts["recycled"], 0, "{what}");
        assert_eq!(counts["recycled_hits"], 0, "{what}");
    }
}

// Where the pool does not build the flat table (one thread), the switch changes nothing and the run
// prints no `[X_STATS] BPLKPFLAT` line. The other case, `LOSAT_X_BPCLUT=0` at 4 threads, runs in
// `the_default_and_the_switches_off_give_one_report` (the same searches, but the switch on by
// default).
#[test]
fn the_switch_does_not_act_without_the_pools_build() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [("LOSAT_X_BPQBLOCK", "1"), ("LOSAT_X_BPCLUT", "1")];
    let counts = same_reports(&subject, 1, &base, "6").stats("BPLKPFLAT");
    assert!(counts.is_empty(), "one thread: {counts:?}");
}

// Split query batches inside the blocks (`CHUNK_SIZE` 1000, `OVERLAP_CHUNK_SIZE` 100): each query
// chunk of a split batch is a search of its own inside its block's search (its preliminary stage
// only), whose table the pool builds and which drops that table when it ends. A chunk takes no
// spare vectors, so the block's table keeps them and the next block reuses them, with the switch:
// the pool builds more tables than there are blocks (the chunks' too), and `recycled` is the blocks
// after the first.
#[test]
fn query_chunks_of_split_batches_take_no_spare() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [
        ("LOSAT_X_BPCLUT", "1"),
        ("LOSAT_X_BPLKPFLAT", "0"),
        ("LOSAT_X_BPQBLOCK", "1"),
        ("CHUNK_SIZE", "1000"),
        ("OVERLAP_CHUNK_SIZE", "100"),
        ("LOSAT_TIMING", "1"),
    ];
    let off = blastp_fast(&subject, 4, &base, "6");
    assert!(!off.stdout.is_empty(), "the fixture pair has rows");
    let blocks = off.searches();
    assert!(blocks > 1, "{blocks} blocks");
    let extra = [ON, STATS];
    let run = blastp_fast(&subject, 4, &with(&base, &extra), "6");
    assert!(run.stdout == off.stdout, "{extra:?}: the report differs");
    assert_eq!(run.searches(), blocks, "{extra:?}");
    let counts = run.printed("BPLKPFLAT");
    let what = format!("{extra:?}: {blocks} blocks, {counts:?}");
    assert!(counts["tables"] > blocks, "{what}");
    assert_eq!(counts["recycled"], blocks - 1, "{what}");
    assert_hits_recycled(&counts, &what);
}

// `-db`: the database search in blocks runs the same lookup region per block. This test also holds
// the blocks of 20,000 letters (several batches per table) for `-subject`.
#[test]
fn database_block_search_reports_with_the_switch_are_the_reports_without_it() {
    let db = subject_database();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", "20000")];
    let counts = same_reports(&subject, 4, &base, "6").printed("BPLKPFLAT");
    assert!(counts["tables"] > 1, "{counts:?}");
    assert_eq!(counts["recycled"], counts["tables"] - 1, "{counts:?}");
    assert_hits_recycled(&counts, "-db");
    std::fs::remove_dir_all(db.parent().unwrap()).ok();
}

// Both switches are on by default: the report with no variable, with either `=0` and with both
// `=0` is one report, in blocks of one batch at 4 threads. `LOSAT_X_BPLKPFLAT` acts only where the
// pool builds the flat table of `LOSAT_X_BPCLUT` and the switch is on, so with no other variable
// the run counts tables exactly when both are on by default, and counts nothing with
// `LOSAT_X_BPCLUT=0`.
#[test]
fn the_default_and_the_switches_off_give_one_report() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let block = ("LOSAT_X_BPQBLOCK", "1");
    let default = blastp_fast(&subject, 4, &[block], "6");
    assert!(!default.stdout.is_empty(), "the fixture pair has rows");
    for vars in [
        vec![block, ("LOSAT_X_BPCLUT", "0")],
        vec![block, OFF],
        vec![block, ("LOSAT_X_BPCLUT", "0"), OFF],
    ] {
        let run = blastp_fast(&subject, 4, &vars, "6");
        assert!(
            run.stdout == default.stdout,
            "{vars:?}: the report differs from the default's"
        );
    }
    let on = blastp_fast(&subject, 4, &[block, STATS], "6");
    assert!(on.stdout == default.stdout, "the run with LOSAT_X_STATS");
    let counts = on.printed("BPLKPFLAT");
    assert!(
        counts["tables"] > 1,
        "no variable: the default builds the flat table on the pool: {counts:?}"
    );
    let off = blastp_fast(&subject, 4, &[block, ("LOSAT_X_BPCLUT", "0"), STATS], "6");
    assert!(off.stdout == default.stdout, "the report, LOSAT_X_BPCLUT=0");
    let counts = off.stats("BPLKPFLAT");
    assert!(
        counts.is_empty(),
        "LOSAT_X_BPCLUT=0 builds no flat table: {counts:?}"
    );
}
