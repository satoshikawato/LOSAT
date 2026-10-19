//! blastp-fast's lookup build with `LOSAT_X_BPLKPFLAT` (`src/algorithm/tblastx/lookup/x_clut_par.rs`,
//! `src/algorithm/blastp/x_query_blocks.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`) with `LOSAT_X_BPCLUT`: it is on by default,
//! `LOSAT_X_BPLKPFLAT=0` turns it off, and `LOSAT_X_DEFAULT=0` leaves it off unless it is set. It acts
//! only where the search's pool builds the flat compressed table of `LOSAT_X_BPCLUT` (on by default as
//! well, `LOSAT_X_BPCLUT=0` turns it off), with more than one thread. With it on, the build's two passes run with one latch
//! each, a block search's table takes the vectors of the previous block's table, and with
//! `LOSAT_X_BPQBLOCK` the contexts' Karlin blocks are tasks of the build. These tests compare the
//! reports with the switch off (`LOSAT_X_BPLKPFLAT=0`) and on, with `-subject` and `-db`, and run its shadow mode
//! (`LOSAT_X_BPLKPFLATSHADOW`), whose line at exit shows that the path ran (tables compared, builds
//! on the vectors of an earlier table, contexts computed by the tasks); a difference stops the run.
//! Every run but the default-set test states `LOSAT_X_BPQBLOCK` and `LOSAT_X_BPCLUT` itself, so that
//! the tests do not depend on the defaults of the strict set.
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
    /// The counters of the `[X_SHADOW] <name>` line, which must end with `identical`.
    fn shadow(&self, name: &str) -> HashMap<String, u64> {
        let prefix = format!("[X_SHADOW] {name} ");
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with(&prefix))
            .unwrap_or_else(|| panic!("no {prefix}line in {}", self.stderr));
        assert!(line.ends_with(" identical"), "{line}");
        line[prefix.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// A protein database of the subject proteome, made with `LOSAT makedb` in a directory of the
/// test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bplkpflat-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("makedb")
        .arg("-in")
        .arg(subject_file())
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
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
const SHADOW: (&str, &str) = ("LOSAT_X_BPLKPFLATSHADOW", "1");

/// The switch off, on and in its shadow mode give the same report; returns the shadow's counters.
fn same_reports(
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    outfmt: &str,
) -> HashMap<String, u64> {
    let off = blastp_fast(subject, threads, &with(base, &[OFF]), outfmt);
    assert!(!off.stdout.is_empty(), "the fixture pair has rows");
    let on = blastp_fast(subject, threads, &with(base, &[ON]), outfmt);
    assert!(
        on.stdout == off.stdout,
        "{base:?} threads {threads} outfmt {outfmt}: the report differs with the switch"
    );
    let shadow = blastp_fast(subject, threads, &with(base, &[SHADOW]), outfmt);
    assert!(
        shadow.stdout == off.stdout,
        "{base:?} threads {threads} outfmt {outfmt}: the report differs in the shadow mode"
    );
    shadow.shadow("BPLKPFLAT")
}

// The block search (every block a search with its own table) at 2 and 4 threads, with blocks of
// one batch and of 20,000 letters: every build is the pool's, compared with the `par_iter` build on
// fresh vectors, every build after a block's first takes the vectors of the table before it, and
// every context's Karlin block comes from a task of the build.
#[test]
fn block_search_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for threads in [2, 4] {
        for letters in ["1", "20000"] {
            let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", letters)];
            let counts = same_reports(&subject, threads, &base, "6");
            let what = format!("threads {threads}, blocks of {letters}: {counts:?}");
            assert!(counts["tables_compared"] > 1, "{what}");
            assert_eq!(counts["outcomes_compared"], 0, "{what}");
            assert_eq!(counts["recycled"], counts["tables_compared"] - 1, "{what}");
            assert_hits_recycled(&counts, &what);
            assert_eq!(
                counts["context_sets_compared"], counts["tables_compared"],
                "{what}"
            );
            assert!(counts["contexts_compared"] > 0, "{what}");
            assert_eq!(
                counts["contexts_from_tasks"], counts["contexts_compared"],
                "{what}"
            );
        }
    }
}

// The search of all queries (a pairwise report with `LOSAT_X_BPQBLOCKFMT=0` is not written in blocks):
// one build, by the pool,
// whose contexts' Karlin blocks are tasks of the build with `LOSAT_X_BPQBLOCK` set and are
// computed as without the switch otherwise; no earlier table.
#[test]
fn search_of_all_queries_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (letters, tasks) in [("1", true), ("0", false)] {
        let base = [
            ("LOSAT_X_BPCLUT", "1"),
            ("LOSAT_X_BPQBLOCK", letters),
            ("LOSAT_X_BPQBLOCKFMT", "0"),
        ];
        let counts = same_reports(&subject, 4, &base, "0");
        let what = format!("blocks {letters}: {counts:?}");
        assert_eq!(counts["tables_compared"], 1, "{what}");
        assert_eq!(counts["recycled"], 0, "{what}");
        assert_eq!(counts["recycled_hits"], 0, "{what}");
        assert_eq!(counts["context_sets_compared"], u64::from(tasks), "{what}");
        assert_eq!(
            counts["contexts_from_tasks"], counts["contexts_compared"],
            "{what}"
        );
    }
}

// Where the pool does not build the flat table (one thread, or `LOSAT_X_BPCLUT` off), the switch
// changes nothing and its shadow compares nothing.
#[test]
fn the_switch_does_not_act_without_the_pools_build() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (threads, clut) in [(1, true), (4, false)] {
        let mut base = vec![("LOSAT_X_BPQBLOCK", "1")];
        base.push(("LOSAT_X_BPCLUT", if clut { "1" } else { "0" }));
        let counts = same_reports(&subject, threads, &base, "6");
        assert!(
            counts.values().all(|&count| count == 0),
            "threads {threads}, BPCLUT {clut}: {counts:?}"
        );
    }
}

// `LOSAT_X_BPCLUTSHADOW` with the switch: the pool's build (one-latch passes, the vectors that the
// previous block's pool build gave back, the Karlin tasks) is compared with the serial flat build
// and the reference's table, which the search uses; the two shadows together as well.
#[test]
fn the_clut_shadow_compares_the_switchs_builds() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", "1")];
    let off = blastp_fast(&subject, 4, &with(&base, &[OFF]), "6");
    for extra in [
        &[("LOSAT_X_BPCLUTSHADOW", "1"), ON][..],
        &[("LOSAT_X_BPCLUTSHADOW", "1"), SHADOW][..],
    ] {
        let run = blastp_fast(&subject, 4, &with(&base, extra), "6");
        assert!(run.stdout == off.stdout, "{extra:?}: the report differs");
        let clut = run.shadow("BPCLUT");
        assert!(clut["parallel_tables_compared"] > 1, "{extra:?}: {clut:?}");
        assert_eq!(
            clut["parallel_tables_compared"], clut["tables_compared"],
            "{extra:?}: {clut:?}"
        );
        if extra.contains(&SHADOW) {
            let counts = run.shadow("BPLKPFLAT");
            assert_eq!(
                counts["tables_compared"], clut["parallel_tables_compared"],
                "{counts:?}"
            );
            assert_eq!(
                counts["recycled"],
                counts["tables_compared"] - 1,
                "{counts:?}"
            );
            assert_hits_recycled(&counts, "BPCLUTSHADOW");
        }
    }
}

// Split query batches inside the blocks (`CHUNK_SIZE` 1000, `OVERLAP_CHUNK_SIZE` 100): each query
// chunk of a split batch is a search of its own inside its block's search (its preliminary stage
// only), whose table the pool builds and which drops that table when it ends. A chunk takes no
// spare vectors, so the block's table keeps them and the next block reuses them, with the switch,
// its shadow and `LOSAT_X_BPCLUTSHADOW` (whose pool tables give their vectors back at once): the
// shadow compares more tables than there are blocks (the chunks' too), and `recycled` is the
// blocks after the first.
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
    let clut_shadow = ("LOSAT_X_BPCLUTSHADOW", "1");
    for extra in [
        &[ON][..],
        &[SHADOW][..],
        &[clut_shadow, ON][..],
        &[clut_shadow, SHADOW][..],
    ] {
        let run = blastp_fast(&subject, 4, &with(&base, extra), "6");
        assert!(run.stdout == off.stdout, "{extra:?}: the report differs");
        assert_eq!(run.searches(), blocks, "{extra:?}");
        let mut tables = None;
        if extra.contains(&clut_shadow) {
            let clut = run.shadow("BPCLUT");
            assert!(
                clut["parallel_tables_compared"] > blocks,
                "{extra:?}: {blocks} blocks, {clut:?}"
            );
            tables = Some(clut["parallel_tables_compared"]);
        }
        if extra.contains(&SHADOW) {
            let counts = run.shadow("BPLKPFLAT");
            let what = format!("{extra:?}: {blocks} blocks, {counts:?}");
            assert!(counts["tables_compared"] > blocks, "{what}");
            assert_eq!(
                tables.unwrap_or(counts["tables_compared"]),
                counts["tables_compared"],
                "{what}"
            );
            assert_eq!(counts["recycled"], blocks - 1, "{what}");
            assert_hits_recycled(&counts, &what);
            assert_eq!(
                counts["contexts_from_tasks"], counts["contexts_compared"],
                "{what}"
            );
        }
    }
}

// `-db`: the database search in blocks runs the same lookup region per block.
#[test]
fn database_block_search_reports_with_the_switch_are_the_reports_without_it() {
    let db = subject_database();
    let subject = [OsStr::new("-db"), db.as_os_str()];
    let base = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPQBLOCK", "20000")];
    let counts = same_reports(&subject, 4, &base, "6");
    assert!(counts["tables_compared"] > 1, "{counts:?}");
    assert_eq!(
        counts["recycled"],
        counts["tables_compared"] - 1,
        "{counts:?}"
    );
    assert_hits_recycled(&counts, "-db");
    assert_eq!(
        counts["contexts_from_tasks"], counts["contexts_compared"],
        "{counts:?}"
    );
    std::fs::remove_dir_all(db.parent().unwrap()).ok();
}

// Both switches are on by default: the report with no variable, with either `=0`, with both `=0`,
// with `LOSAT_X_DEFAULT=0` (both off unless set) and with `LOSAT_X_DEFAULT=0` and both set is one
// report, in blocks of one batch at 4 threads. The shadow of `LOSAT_X_BPLKPFLAT` acts only where the
// pool builds the flat table of `LOSAT_X_BPCLUT` and the switch is on, so with no other variable it
// compares tables exactly when both are on by default, and nothing under `LOSAT_X_DEFAULT=0`.
#[test]
fn the_default_the_switches_off_and_the_global_off_give_one_report() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let block = ("LOSAT_X_BPQBLOCK", "1");
    let global_off = ("LOSAT_X_DEFAULT", "0");
    let default = blastp_fast(&subject, 4, &[block], "6");
    assert!(!default.stdout.is_empty(), "the fixture pair has rows");
    for vars in [
        vec![block, ("LOSAT_X_BPCLUT", "0")],
        vec![block, OFF],
        vec![block, ("LOSAT_X_BPCLUT", "0"), OFF],
        vec![block, global_off],
        vec![block, global_off, ("LOSAT_X_BPCLUT", "1"), ON],
    ] {
        let run = blastp_fast(&subject, 4, &vars, "6");
        assert!(
            run.stdout == default.stdout,
            "{vars:?}: the report differs from the default's"
        );
    }
    let shadow = blastp_fast(&subject, 4, &[block, SHADOW], "6");
    assert!(shadow.stdout == default.stdout, "the shadow's report");
    let counts = shadow.shadow("BPLKPFLAT");
    assert!(
        counts["tables_compared"] > 1,
        "no variable but the shadow: the default builds the flat table on the pool: {counts:?}"
    );
    let off = blastp_fast(&subject, 4, &[block, global_off, SHADOW], "6");
    assert!(
        off.stdout == default.stdout,
        "the shadow's report, global off"
    );
    let counts = off.shadow("BPLKPFLAT");
    assert!(
        counts.values().all(|&count| count == 0),
        "LOSAT_X_DEFAULT=0 builds no flat table: {counts:?}"
    );
}
