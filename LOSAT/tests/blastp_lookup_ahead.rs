//! blastp's block search with `LOSAT_X_BPLKPAHEAD` (`src/algorithm/blastp/blast_engine/
//! x_lookup_ahead.rs`): the next query block's query frames, lookup table and contexts built during
//! the block before it, on the search's pool.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`), but every run here names it (see
//! `PREVIOUS_DEFAULT`). It acts in the block search of `LOSAT_X_BPQBLOCK` with more than one
//! thread: for blastp-fast where `LOSAT_X_BPLKPFLAT` acts on the flat table of `LOSAT_X_BPCLUT`,
//! and for exact blastp where its table is built in place (`LOSAT_X_BPLUT`, outside
//! `LOSAT_DIAGNOSTICS`); these are in the strict set too, so "off" for any of them is `=0`. These
//! tests compare the reports, standard error and exit codes with the switch off, on and in its
//! shadow mode (`LOSAT_X_BPLKPAHEADSHADOW`, which builds every block's values again inline and
//! compares them), with `-subject` and `-db`, and read the shadow's line at exit: `blocks` (blocks
//! searched where the switch acts), `prebuilt` (blocks that took values built ahead), `dropped`
//! (builds ahead that gave none), `compared`. Every run states `LOSAT_X_BPQBLOCK`, and
//! `LOSAT_X_BPCLUT` and `LOSAT_X_BPLKPFLAT` (blastp-fast) or `LOSAT_X_BPLUT` (exact blastp),
//! itself, so that the tests do not depend on the defaults of the strict set. Each test writes its
//! files to a directory of its own.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215
//! ```c++
//!     // 5. Create the lookup table
//!     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
//!         LookupTableWrap* lut =
//!             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//! ```

use std::collections::HashMap;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of the test `tag` (removed first).
fn test_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bplkpahead-{tag}-{}", std::process::id()));
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// What one run printed and how it ended.
#[derive(Debug, PartialEq)]
struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

impl Run {
    /// Standard error without the `[TIMING]` lines and the shadow lines of the switch and of
    /// `LOSAT_X_BPLKPFLATSHADOW` (whose counts include the builds ahead).
    fn stderr_compared(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| {
                !line.starts_with("[TIMING]")
                    && !line.starts_with("[X_SHADOW] BPLKPAHEAD ")
                    && !line.starts_with("[X_SHADOW] BPLKPFLAT ")
            })
            .map(|line| format!("{line}\n"))
            .collect()
    }

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
    /// `LOSAT_TIMING`): the blocks of a block search.
    fn searches(&self) -> u64 {
        self.stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] blastp_total:"))
            .count() as u64
    }
}

/// One blastp run (`task`) of the fixture pair (two related nimavirus proteomes) against `subject`
/// (`-subject FILE` or `-db NAME`) with `vars` as its only `LOSAT_*` variables that choose a code
/// path, `BATCH_SIZE` 3000 (several query batches), `LOSAT_TIMING` and the arguments `extra`
/// (`-outfmt 6` unless they name another).
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file in `dir` is named.
fn blastp(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
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
        .arg(fixture("MeenMJNV.faa"))
        .args(subject)
        .args(["-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra);
    if !extra.contains(&"-outfmt") {
        command.args(["-outfmt", "6"]);
    }
    let output = command
        .env("BATCH_SIZE", "3000")
        .env("LOSAT_TIMING", "1")
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

const ON: (&str, &str) = ("LOSAT_X_BPLKPAHEAD", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_BPLKPAHEADSHADOW", "1");
/// Where the switch acts: the block search with the pool's build of the flat table.
const ACTS: [(&str, &str); 2] = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPLKPFLAT", "1")];
/// Where the switch acts for exact blastp: the block search with the table built in place.
const EXACT_ACTS: [(&str, &str); 1] = [("LOSAT_X_BPLUT", "1")];
/// Both tasks with the variables under which the switch acts for them.
const TASKS: [(&str, &[(&str, &str)]); 2] = [("blastp-fast", &ACTS), ("blastp", &EXACT_ACTS)];
/// The composition-based redo's nested scopes (redo windows) run beside the build ahead.
// TODO after the rebase onto main: also LOSAT_X_BPSUBJRUN (its subject runs are not on this base).
const REDOPAR: (&str, &str) = ("LOSAT_X_BPREDOPAR", "1");

/// The switch off, on and in its shadow mode give the same report, standard error (but the
/// `[TIMING]` lines and the shadow's line) and exit code; returns the off run and the shadow run.
fn same_runs(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    extra: &[&str],
) -> (Run, Run) {
    let off = blastp(dir, task, subject, threads, base, extra);
    let what = format!("{task} {base:?} {extra:?} threads {threads}");
    for switch in [ON, SHADOW] {
        let run = blastp(dir, task, subject, threads, &with(base, &[switch]), extra);
        assert!(run.stdout == off.stdout, "{what} {switch:?}: the report");
        assert_eq!(
            run.stderr_compared(),
            off.stderr_compared(),
            "{what} {switch:?}: standard error"
        );
        assert_eq!(run.code, off.code, "{what} {switch:?}: exit code");
        if switch == SHADOW {
            return (off, run);
        }
    }
    unreachable!()
}

/// The shadow's counters where the switch acts on `searches` blocks: every block counted, every
/// block after the first built ahead (taken or dropped), every value taken compared, and at least
/// one taken.
fn assert_built_ahead(counts: &HashMap<String, u64>, searches: u64, what: &str) {
    assert!(searches > 1, "{what}: {searches} blocks");
    assert_eq!(counts["blocks"], searches, "{what}: {counts:?}");
    assert_eq!(
        counts["prebuilt"] + counts["dropped"],
        searches - 1,
        "{what}: {counts:?}"
    );
    assert!(counts["prebuilt"] > 0, "{what}: {counts:?}");
    assert_eq!(counts["compared"], counts["prebuilt"], "{what}: {counts:?}");
}

// The block search with blocks of one batch and of 20,000 letters: each block after the first
// takes the values built ahead during the block before it, and the shadow finds them equal to its
// own build. Both tasks: blastp-fast's flat table and exact blastp's table. One case per path (each
// search of the fixture pair costs CI seconds): four threads (every gate is `threads() > 1`; two
// threads, where the build ahead has the block's one helper: the `-db` test and the `/dev/full`
// test); the build beside the redo windows of `LOSAT_X_BPREDOPAR`: the `-db` test.
#[test]
fn block_search_reports_with_the_switch_are_the_reports_without_it() {
    let dir = test_dir("blocks");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (task, acts) in TASKS {
        for letters in ["1", "20000"] {
            let base = with(acts, &[("LOSAT_X_BPQBLOCK", letters)]);
            let (off, shadow) = same_runs(&dir, task, &subject, 4, &base, &[]);
            assert!(!off.stdout.is_empty(), "the fixture pair has rows");
            let what = format!("{task}, threads 4, blocks of {letters}");
            let counts = shadow.shadow("BPLKPAHEAD");
            assert_built_ahead(&counts, off.searches(), &what);
            // Exact blastp in blocks of one batch: a block whose table's overflow exceeds the
            // room (the largest so far plus a half) builds its table itself, after its
            // thread took the neighbour lists back.
            if task == "blastp" && letters == "1" {
                assert!(counts["dropped"] > 0, "{what}: {counts:?}");
            }
        }
    }
    // Exact blastp's pairwise report written in blocks (`LOSAT_X_BPQBLOCKFMT`): the build ahead
    // beside the block's report.
    let base = with(&EXACT_ACTS, &[("LOSAT_X_BPQBLOCK", "1")]);
    let (off, shadow) = same_runs(&dir, "blastp", &subject, 4, &base, &["-outfmt", "0"]);
    assert!(off.searches() > 1, "{} blocks", off.searches());
    assert_built_ahead(
        &shadow.shadow("BPLKPAHEAD"),
        off.searches(),
        "blastp -outfmt 0, threads 4",
    );
    std::fs::remove_dir_all(&dir).ok();
}

// `-db`: the database search in blocks runs the same block loop and lookup region. With
// `LOSAT_X_BPREDOPAR`: the build ahead beside the redo windows' nested scopes of the block's search
// (the hang of a waiting helper would show here), blastp-fast at two threads (the block's one
// helper; its build ahead waits for its own nested passes, exact blastp's is a leaf task) and exact
// blastp at four.
#[test]
fn database_block_search_reports_with_the_switch_are_the_reports_without_it() {
    let dir = test_dir("db");
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
    let subject = [OsStr::new("-db"), name.as_os_str()];
    for ((task, acts), threads) in TASKS.into_iter().zip([2, 4]) {
        let base = with(acts, &[("LOSAT_X_BPQBLOCK", "1"), REDOPAR]);
        let (off, shadow) = same_runs(&dir, task, &subject, threads, &base, &[]);
        assert!(!off.stdout.is_empty(), "the fixture pair has rows");
        let what = format!("{task} -db, threads {threads}");
        assert_built_ahead(&shadow.shadow("BPLKPAHEAD"), off.searches(), &what);
    }
    std::fs::remove_dir_all(&dir).ok();
}

// One block (blocks larger than the queries): the block plan takes the search of all queries for it
// (`x_query_block_plan`), as for `LOSAT_X_BPQBLOCK=0` and for a pairwise report with
// `LOSAT_X_BPQBLOCKFMT=0` (not written in blocks): one search, nothing built ahead, nothing counted.
// The last block of a block search builds nothing ahead either (`assert_built_ahead`: blocks - 1
// builds). The three are refusals of `x_query_block_plan`, which does not read the task, and the
// switch's state is made in the block loop only (`XLookupAhead::new`): each runs once, the tasks
// in turn.
#[test]
fn one_block_and_the_search_of_all_queries_build_nothing_ahead() {
    let dir = test_dir("one");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let [(fast, fast_acts), (exact, exact_acts)] = TASKS;
    for (task, acts, blocks, extra) in [
        (
            fast,
            fast_acts,
            &[("LOSAT_X_BPQBLOCK", "100000000")][..],
            &[][..],
        ),
        (exact, exact_acts, &[("LOSAT_X_BPQBLOCK", "0")][..], &[][..]),
        (
            fast,
            fast_acts,
            &[("LOSAT_X_BPQBLOCK", "1"), ("LOSAT_X_BPQBLOCKFMT", "0")][..],
            &["-outfmt", "0"][..],
        ),
    ] {
        let base = with(acts, blocks);
        let (off, shadow) = same_runs(&dir, task, &subject, 4, &base, extra);
        assert_eq!(off.searches(), 1, "{task} {blocks:?} {extra:?}");
        let counts = shadow.shadow("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "{task} {blocks:?} {extra:?}: {counts:?}"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Where `LOSAT_X_BPLKPFLAT` does not act on the pool's flat table (one thread, the switch =0 or
// `LOSAT_X_BPCLUT=0`, `LOSAT_X_BPCLUTSHADOW`, whose search uses the reference builder's table), the
// switch changes nothing for blastp-fast and its shadow counts nothing.
#[test]
fn the_switch_does_not_act_without_the_pools_flat_table() {
    let dir = test_dir("off");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let block = ("LOSAT_X_BPQBLOCK", "1");
    let clut = ("LOSAT_X_BPCLUT", "1");
    let flat = ("LOSAT_X_BPLKPFLAT", "1");
    // Both are in the strict set (on by default): "off" is `=0`.
    let clut_off = ("LOSAT_X_BPCLUT", "0");
    let flat_off = ("LOSAT_X_BPLKPFLAT", "0");
    let clut_shadow = ("LOSAT_X_BPCLUTSHADOW", "1");
    for (task, threads, base) in [
        ("blastp-fast", 1, vec![block, clut, flat]),
        ("blastp-fast", 4, vec![block, clut, flat_off]),
        ("blastp-fast", 4, vec![block, clut_off, flat]),
        ("blastp-fast", 4, vec![block, clut_shadow, flat]),
    ] {
        let (off, shadow) = same_runs(&dir, task, &subject, threads, &base, &[]);
        assert!(
            off.searches() > 1,
            "{task} {base:?}: {} blocks",
            off.searches()
        );
        let counts = shadow.shadow("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "{task} threads {threads} {base:?}: {counts:?}"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Exact blastp where its table is not built in place by the pool's search (one thread,
// `LOSAT_X_BPLUTSHADOW`, which builds the reference builder's table, and `LOSAT_DIAGNOSTICS`, whose
// reference builder prints its phases in order): the switch changes nothing and its shadow counts
// nothing. `LOSAT_X_BPLUT=0` is the shadow's branch (`XAaSpares::x_applies`: mode 1 only). The flat
// table's switches do not decide it.
#[test]
fn the_switch_does_not_act_for_exact_blastp_without_the_in_place_table() {
    let dir = test_dir("exact-off");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let block = ("LOSAT_X_BPQBLOCK", "1");
    let lut = ("LOSAT_X_BPLUT", "1");
    for (threads, base) in [
        (1, vec![block, lut]),
        (4, vec![block, ("LOSAT_X_BPLUTSHADOW", "1")]),
        (4, vec![block, lut, ("LOSAT_DIAGNOSTICS", "1")]),
    ] {
        let (off, shadow) = same_runs(&dir, "blastp", &subject, threads, &base, &[]);
        assert!(off.searches() > 1, "{base:?}: {} blocks", off.searches());
        let counts = shadow.shadow("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "threads {threads} {base:?}: {counts:?}"
        );
    }
    // Without the flat table's switches the switch still acts for exact blastp.
    let base = [
        block,
        lut,
        ("LOSAT_X_BPCLUT", "0"),
        ("LOSAT_X_BPLKPFLAT", "0"),
    ];
    let (off, shadow) = same_runs(&dir, "blastp", &subject, 4, &base, &[]);
    assert_built_ahead(
        &shadow.shadow("BPLKPAHEAD"),
        off.searches(),
        "exact, BPCLUT=0 BPLKPFLAT=0",
    );
    std::fs::remove_dir_all(&dir).ok();
}

// Split query batches inside the blocks (`CHUNK_SIZE` 1000, `OVERLAP_CHUNK_SIZE` 100): each query
// chunk builds its own table inside its block's search (nothing ahead) while the next block's table
// is built ahead; for blastp-fast with the shadow of `LOSAT_X_BPLKPFLAT` too, which then also
// checks the builds ahead (one-latch passes on the second spare's vectors, poisoned first), and for
// exact blastp.
#[test]
fn query_chunks_of_split_batches_build_their_own_tables() {
    let dir = test_dir("chunks");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let chunks = [
        ("LOSAT_X_BPQBLOCK", "1"),
        ("CHUNK_SIZE", "1000"),
        ("OVERLAP_CHUNK_SIZE", "100"),
    ];
    for (task, flat_shadow) in [
        ("blastp-fast", false),
        ("blastp-fast", true),
        ("blastp", false),
    ] {
        let acts: &[(&str, &str)] = if task == "blastp" { &EXACT_ACTS } else { &ACTS };
        let mut base = with(acts, &chunks);
        if flat_shadow {
            base.push(("LOSAT_X_BPLKPFLATSHADOW", "1"));
        }
        let (off, shadow) = same_runs(&dir, task, &subject, 4, &base, &[]);
        assert!(!off.stdout.is_empty(), "the fixture pair has rows");
        let what = format!("{task} split batches, BPLKPFLATSHADOW {flat_shadow}");
        let counts = shadow.shadow("BPLKPAHEAD");
        assert_built_ahead(&counts, off.searches(), &what);
        if flat_shadow {
            let flat = shadow.shadow("BPLKPFLAT");
            assert!(
                flat["tables_compared"] > off.searches() + counts["prebuilt"],
                "{what}: {flat:?}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

// A block that fails after the next block's build has started (the report written to `/dev/full`
// fails at the first block's flush): the run ends with that block's error, standard error and exit
// code as without the switch; the values built ahead are dropped and no shadow line is printed.
#[cfg(target_os = "linux")]
#[test]
fn a_block_that_fails_after_the_next_build_started_ends_the_run_alike() {
    let dir = test_dir("full");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for ((task, acts), threads) in TASKS.into_iter().flat_map(|t| [(t, 2), (t, 4)]) {
        let base = with(acts, &[("LOSAT_X_BPQBLOCK", "1")]);
        let (off, shadow) = same_runs(&dir, task, &subject, threads, &base, &["-out", "/dev/full"]);
        assert_ne!(
            off.code,
            Some(0),
            "{task} threads {threads}: {}",
            off.stderr
        );
        assert!(off.stdout.is_empty());
        assert!(
            !shadow.stderr.contains("[X_SHADOW] BPLKPAHEAD"),
            "{}",
            shadow.stderr
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}
