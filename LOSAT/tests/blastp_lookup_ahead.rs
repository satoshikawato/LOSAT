//! blastp's block search with `LOSAT_X_BPLKPAHEAD` (`src/algorithm/blastp/blast_engine/
//! x_lookup_ahead.rs`): the next query block's query frames, lookup table and contexts built during
//! the block before it, on the search's pool.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`), but every run here names it (see
//! `PREVIOUS_DEFAULT`). It acts in the block search of `LOSAT_X_BPQBLOCK` with more than one
//! thread: for blastp-fast where `LOSAT_X_BPLKPFLAT` acts on the flat table of `LOSAT_X_BPCLUT`,
//! and for exact blastp, whose table is built in place outside `LOSAT_DIAGNOSTICS`; the switches
//! of the flat table are in the strict set too, so "off" for any of them is `=0`. These
//! tests compare the reports, standard error and exit codes with the switch off and on, with
//! `-subject` and `-db`, and read the `[X_STATS] BPLKPAHEAD` line of the run with the switch on
//! (printed with `LOSAT_X_STATS`, when a counter is not zero) at exit: `blocks` (blocks searched
//! where the switch acts), `prebuilt` (blocks that took values built ahead), `dropped` (builds
//! ahead that gave none). That the values built ahead are the values the block builds itself is a
//! unit test of `x_lookup_ahead.rs`. Every run states `LOSAT_X_BPQBLOCK`, and
//! `LOSAT_X_BPCLUT` and `LOSAT_X_BPLKPFLAT` (blastp-fast),
//! itself, so that the tests do not depend on the defaults of the strict set. Each test writes its
//! files to a directory of its own.
//!
//! The queries are records of `MeenMJNV.faa`, not all 101 (one run on all 101 queries costs 2-8
//! CPU s, one on `four_queries` 0.05-0.45 s, one on `drop_queries` 0.4-0.7 s): `four_queries` for
//! every path that needs blocks, `drop_queries` for the one rule that reads the tables' sizes (the
//! room of a build ahead). The other rules depend on the letters only through counts that the four
//! queries keep: a block search needs more than one block, and a split batch at least two chunks
//! (`CHUNK_SIZE` 400).
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
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
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

/// The queries of a search (a file of records of `MeenMJNV.faa`) and its `BATCH_SIZE`.
struct Queries {
    file: PathBuf,
    batch_size: &'static str,
}

/// The records of `MeenMJNV.faa` from accession `first` to `last`, in the file `name` of `dir`.
fn records(dir: &Path, name: &str, first: &str, last: &str) -> PathBuf {
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
    let path = dir.join(name);
    std::fs::write(&path, records[position(first)..=position(last)].concat()).unwrap();
    path
}

/// The queries of most searches: four consecutive records, BDT62581.1 (698 letters, one subject),
/// BDT62582.1 (268, another), BDT62583.1 (629, no hits) and BDT62584.1 (944, five subjects), each a
/// batch of its own at `BATCH_SIZE` 250 (a batch takes records while it holds fewer letters): four
/// blocks of one batch at `LOSAT_X_BPQBLOCK=1`, two of two batches at `TWO_BATCH_BLOCKS`. No build
/// ahead lacks the room (the vectors made ready for it hold at least 2^18 entries, more than any of
/// their tables needs).
fn four_queries(dir: &Path) -> Queries {
    Queries {
        file: records(dir, "four.faa", "BDT62581.1", "BDT62584.1"),
        batch_size: "250",
    }
}

/// `LOSAT_X_BPQBLOCK` for `four_queries` in blocks of two batches (a block takes batches while it
/// holds fewer letters): {BDT62581.1, BDT62582.1} (966 letters) and {BDT62583.1, BDT62584.1}
/// (1,573), the second built ahead.
const TWO_BATCH_BLOCKS: &str = "900";

/// The queries where a build ahead lacks the room: BDT62535.1 to BDT62547.1, the fixture's third to
/// fifth query batches at `BATCH_SIZE` 3000 (3,183, 9,157 and 3,131 letters), three blocks of one
/// batch at `LOSAT_X_BPQBLOCK=1`. The second batch's table (it holds BDT62546.1, 6,458 letters) has
/// more hits (exact blastp: overflow entries) than the vectors made ready for it during the first
/// block hold (`x_lookup_ahead.rs` `x_make_room`: new ones of twice the largest table so far, at
/// least 2^18 entries), so that block builds its table itself; the third takes the values built
/// ahead on the vectors the second made ready. With all 101 queries in blocks of one batch the same
/// batch (the fourth block) is the one dropped.
fn drop_queries(dir: &Path) -> Queries {
    Queries {
        file: records(dir, "drop.faa", "BDT62535.1", "BDT62547.1"),
        batch_size: "3000",
    }
}

/// What one run printed and how it ended.
#[derive(Debug, PartialEq)]
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
    /// `LOSAT_TIMING`): the blocks of a block search.
    fn searches(&self) -> u64 {
        self.stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] blastp_total:"))
            .count() as u64
    }
}

/// One blastp run (`task`) of `queries` against `subject` (`-subject FILE` or `-db NAME`; the
/// fixture pair is two related nimavirus proteomes) with `vars` as its only `LOSAT_*` variables
/// that choose a code path, the queries' `BATCH_SIZE`, `LOSAT_TIMING` and the arguments `extra`
/// (`-outfmt 6` unless they name another).
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file in `dir` is named.
fn blastp(
    dir: &Path,
    task: &str,
    queries: &Queries,
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
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
/// Where the switch acts: the block search with the pool's build of the flat table.
const ACTS: [(&str, &str); 2] = [("LOSAT_X_BPCLUT", "1"), ("LOSAT_X_BPLKPFLAT", "1")];
/// Where the switch acts for exact blastp: the block search with the table built in place.
const EXACT_ACTS: [(&str, &str); 0] = [];
/// Both tasks with the variables under which the switch acts for them.
const TASKS: [(&str, &[(&str, &str)]); 2] = [("blastp-fast", &ACTS), ("blastp", &EXACT_ACTS)];
/// The composition-based redo's nested scopes (redo windows) run beside the build ahead.
// TODO after the rebase onto main: also LOSAT_X_BPSUBJRUN (its subject runs are not on this base).
const REDOPAR: (&str, &str) = ("LOSAT_X_BPREDOPAR", "1");

/// The switch off and on give the same report, standard error (but the `[TIMING]` and `[X_STATS]`
/// lines) and exit code; returns the off run and the run with the switch on (which also has
/// `LOSAT_X_STATS`).
fn same_runs(
    dir: &Path,
    task: &str,
    queries: &Queries,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    extra: &[&str],
) -> (Run, Run) {
    let off = blastp(dir, task, queries, subject, threads, base, extra);
    let what = format!("{task} {base:?} {extra:?} threads {threads}");
    let on = blastp(
        dir,
        task,
        queries,
        subject,
        threads,
        &with(base, &[ON, STATS]),
        extra,
    );
    assert!(on.stdout == off.stdout, "{what}: the report");
    assert_eq!(
        on.stderr_compared(),
        off.stderr_compared(),
        "{what}: standard error"
    );
    assert_eq!(on.code, off.code, "{what}: exit code");
    (off, on)
}

/// The counters of `[X_STATS] BPLKPAHEAD` where the switch acts on `searches` blocks: every block
/// counted, every block after the first built ahead (taken or dropped), and at least one taken.
fn assert_built_ahead(counts: &HashMap<String, u64>, searches: u64, what: &str) {
    assert!(searches > 1, "{what}: {searches} blocks");
    assert_eq!(counts["blocks"], searches, "{what}: {counts:?}");
    assert_eq!(
        counts["prebuilt"] + counts["dropped"],
        searches - 1,
        "{what}: {counts:?}"
    );
    assert!(counts["prebuilt"] > 0, "{what}: {counts:?}");
}

// The block search with blocks of one batch and of two: each block after the first takes the values
// built ahead during the block before it where they fit (below). Both tasks: blastp-fast's flat
// table and exact blastp's table. Blocks of one batch with `drop_queries`, whose second block's build ahead lacks the room; blocks of two batches with
// `four_queries` (blastp-fast's build ahead makes its table of the block's query batches). One
// case per path: four threads (every gate is `threads() > 1`; two threads, where the build ahead
// has the block's one helper: the `-db` test and the `/dev/full` test); the build beside the redo
// windows of `LOSAT_X_BPREDOPAR`: the `-db` test.
#[test]
fn block_search_reports_with_the_switch_are_the_reports_without_it() {
    let dir = test_dir("blocks");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let (four, drop) = (four_queries(&dir), drop_queries(&dir));
    for (task, acts) in TASKS {
        for (queries, letters, blocks, dropped) in
            [(&drop, "1", 3, 1), (&four, TWO_BATCH_BLOCKS, 2, 0)]
        {
            let base = with(acts, &[("LOSAT_X_BPQBLOCK", letters)]);
            let (off, on) = same_runs(&dir, task, queries, &subject, 4, &base, &[]);
            assert!(!off.stdout.is_empty(), "the fixture pair has rows");
            let what = format!("{task}, threads 4, blocks of {letters}");
            assert_eq!(off.searches(), blocks, "{what}");
            let counts = on.printed("BPLKPAHEAD");
            assert_built_ahead(&counts, off.searches(), &what);
            // `drop_queries`: the second block, whose table's hits (exact blastp: overflow) do not
            // fit the vectors made ready, builds its table itself (exact blastp: after its thread
            // took the neighbour lists back); the third takes the values built ahead.
            assert_eq!(counts["dropped"], dropped, "{what}: {counts:?}");
        }
    }
    // Exact blastp's pairwise report written in blocks: the build ahead beside the block's report.
    let base = with(&EXACT_ACTS, &[("LOSAT_X_BPQBLOCK", "1")]);
    let (off, on) = same_runs(&dir, "blastp", &four, &subject, 4, &base, &["-outfmt", "0"]);
    assert_eq!(off.searches(), 4, "blastp -outfmt 0");
    assert_built_ahead(
        &on.printed("BPLKPAHEAD"),
        off.searches(),
        "blastp -outfmt 0, threads 4",
    );
    std::fs::remove_dir_all(&dir).ok();
}

// `-db`: the database search in blocks runs the same block loop and lookup region. With
// `LOSAT_X_BPREDOPAR`: the build ahead beside the redo windows' nested scopes of the block's search
// (the hang of a waiting helper would show here; every block of `four_queries` redoes its 8-58
// matches on the pool, one `[TIMING] x_bpredopar` line each), blastp-fast at two threads (the
// block's one helper; its build ahead waits for its own nested passes, exact blastp's is a leaf
// task) and exact blastp at four.
#[test]
fn database_block_search_reports_with_the_switch_are_the_reports_without_it() {
    let dir = test_dir("db");
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
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
    let subject = [OsStr::new("-db"), name.as_os_str()];
    let queries = four_queries(&dir);
    for ((task, acts), threads) in TASKS.into_iter().zip([2, 4]) {
        let base = with(acts, &[("LOSAT_X_BPQBLOCK", "1"), REDOPAR]);
        let (off, on) = same_runs(&dir, task, &queries, &subject, threads, &base, &[]);
        assert!(!off.stdout.is_empty(), "the fixture pair has rows");
        let what = format!("{task} -db, threads {threads}");
        assert_eq!(off.searches(), 4, "{what}");
        let redone = off
            .stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] x_bpredopar:"))
            .count() as u64;
        assert_eq!(
            redone,
            off.searches(),
            "{what}: every block's redo on the pool"
        );
        assert_built_ahead(&on.printed("BPLKPAHEAD"), off.searches(), &what);
    }
    std::fs::remove_dir_all(&dir).ok();
}

// One block (blocks larger than the queries): the block plan takes the search of all queries for it
// (`x_query_block_plan`), as for `LOSAT_X_BPQBLOCK=0` and for the built-in block (50,000 letters,
// above the four queries' 2,539; the block switch stays on, `BATCH_SIZE` does not change the block
// size) with a pairwise report: one search, nothing built ahead, nothing counted.
// The last block of a block search builds nothing ahead either (`assert_built_ahead`: blocks - 1
// builds). The three are refusals of `x_query_block_plan`, which does not read the task, and the
// switch's state is made in the block loop only (`XLookupAhead::new`): each runs once, the tasks
// in turn.
#[test]
fn one_block_and_the_search_of_all_queries_build_nothing_ahead() {
    let dir = test_dir("one");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let queries = four_queries(&dir);
    let [(fast, fast_acts), (exact, exact_acts)] = TASKS;
    for (task, acts, blocks, extra) in [
        (
            fast,
            fast_acts,
            &[("LOSAT_X_BPQBLOCK", "100000000")][..],
            &[][..],
        ),
        (exact, exact_acts, &[("LOSAT_X_BPQBLOCK", "0")][..], &[][..]),
        (fast, fast_acts, &[][..], &["-outfmt", "0"][..]),
    ] {
        let base = with(acts, blocks);
        let (off, on) = same_runs(&dir, task, &queries, &subject, 4, &base, extra);
        assert_eq!(off.searches(), 1, "{task} {blocks:?} {extra:?}");
        let counts = on.stats("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "{task} {blocks:?} {extra:?}: {counts:?}"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Where `LOSAT_X_BPLKPFLAT` does not act on the pool's flat table (one thread, the switch =0 or
// `LOSAT_X_BPCLUT=0`), the switch changes nothing for blastp-fast and the run prints no
// `[X_STATS] BPLKPAHEAD` line.
#[test]
fn the_switch_does_not_act_without_the_pools_flat_table() {
    let dir = test_dir("off");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let queries = four_queries(&dir);
    let block = ("LOSAT_X_BPQBLOCK", "1");
    let clut = ("LOSAT_X_BPCLUT", "1");
    let flat = ("LOSAT_X_BPLKPFLAT", "1");
    // Both are in the strict set (on by default): "off" is `=0`.
    let clut_off = ("LOSAT_X_BPCLUT", "0");
    let flat_off = ("LOSAT_X_BPLKPFLAT", "0");
    for (task, threads, base) in [
        ("blastp-fast", 1, vec![block, clut, flat]),
        ("blastp-fast", 4, vec![block, clut, flat_off]),
        ("blastp-fast", 4, vec![block, clut_off, flat]),
    ] {
        let (off, on) = same_runs(&dir, task, &queries, &subject, threads, &base, &[]);
        assert_eq!(off.searches(), 4, "{task} {base:?}");
        let counts = on.stats("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "{task} threads {threads} {base:?}: {counts:?}"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Exact blastp where its table is not built in place by the pool's search (one thread, and
// `LOSAT_DIAGNOSTICS`, whose reference builder prints its phases in order): the switch changes
// nothing and the run prints no `[X_STATS] BPLKPAHEAD` line (`XAaSpares::x_applies`: the in-place
// build only). The flat table's switches do not decide it.
#[test]
fn the_switch_does_not_act_for_exact_blastp_without_the_in_place_table() {
    let dir = test_dir("exact-off");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let queries = four_queries(&dir);
    let block = ("LOSAT_X_BPQBLOCK", "1");
    for (threads, base) in [
        (1, vec![block]),
        (4, vec![block, ("LOSAT_DIAGNOSTICS", "1")]),
    ] {
        let (off, on) = same_runs(&dir, "blastp", &queries, &subject, threads, &base, &[]);
        assert_eq!(off.searches(), 4, "{base:?}");
        let counts = on.stats("BPLKPAHEAD");
        assert!(
            counts.values().all(|&count| count == 0),
            "threads {threads} {base:?}: {counts:?}"
        );
    }
    // Without the flat table's switches the switch still acts for exact blastp.
    let base = [block, ("LOSAT_X_BPCLUT", "0"), ("LOSAT_X_BPLKPFLAT", "0")];
    let (off, on) = same_runs(&dir, "blastp", &queries, &subject, 4, &base, &[]);
    assert_built_ahead(
        &on.printed("BPLKPAHEAD"),
        off.searches(),
        "exact, BPCLUT=0 BPLKPFLAT=0",
    );
    std::fs::remove_dir_all(&dir).ok();
}

// Split query batches inside the blocks (`CHUNK_SIZE` 400, `OVERLAP_CHUNK_SIZE` 100: a batch is
// split into letters / 300 chunks where that is two or more, `calculate_num_chunks`): each query
// chunk builds its own table inside its block's search (nothing ahead) while the next block's table
// is built ahead; for blastp-fast the pool's builds counted by `[X_STATS] BPLKPFLAT` then include
// the builds ahead (one-latch passes on the second spare's vectors), and for exact blastp. `four_queries` at `BATCH_SIZE` 500: three blocks of one batch, {BDT62581.1} into 2
// chunks, {BDT62582.1, BDT62583.1} (897 letters) into 2, the first holding parts of both queries,
// and {BDT62584.1} into 3.
#[test]
fn query_chunks_of_split_batches_build_their_own_tables() {
    let dir = test_dir("chunks");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let queries = Queries {
        batch_size: "500",
        ..four_queries(&dir)
    };
    let chunks = [
        ("LOSAT_X_BPQBLOCK", "1"),
        ("CHUNK_SIZE", "400"),
        ("OVERLAP_CHUNK_SIZE", "100"),
    ];
    for task in ["blastp-fast", "blastp"] {
        let acts: &[(&str, &str)] = if task == "blastp" { &EXACT_ACTS } else { &ACTS };
        let base = with(acts, &chunks);
        let (off, on) = same_runs(&dir, task, &queries, &subject, 4, &base, &[]);
        assert!(!off.stdout.is_empty(), "the fixture pair has rows");
        let what = format!("{task} split batches");
        assert_eq!(off.searches(), 3, "{what}");
        let counts = on.printed("BPLKPAHEAD");
        assert_built_ahead(&counts, off.searches(), &what);
        if task == "blastp-fast" {
            let flat = on.printed("BPLKPFLAT");
            assert!(
                flat["tables"] > off.searches() + counts["prebuilt"],
                "{what}: {flat:?}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

// A block that fails after the next block's build has started (the report written to `/dev/full`
// fails at the first block's flush; BDT62581.1, the first block of `four_queries`, has a row): the
// run ends with that block's error, standard error and exit code as without the switch; the values
// built ahead are dropped and no `[X_STATS] BPLKPAHEAD` line is printed.
#[cfg(target_os = "linux")]
#[test]
fn a_block_that_fails_after_the_next_build_started_ends_the_run_alike() {
    let dir = test_dir("full");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let queries = four_queries(&dir);
    for ((task, acts), threads) in TASKS.into_iter().flat_map(|t| [(t, 2), (t, 4)]) {
        let base = with(acts, &[("LOSAT_X_BPQBLOCK", "1")]);
        let (off, on) = same_runs(
            &dir,
            task,
            &queries,
            &subject,
            threads,
            &base,
            &["-out", "/dev/full"],
        );
        assert_ne!(
            off.code,
            Some(0),
            "{task} threads {threads}: {}",
            off.stderr
        );
        assert!(off.stdout.is_empty());
        assert!(!on.stderr.contains("[X_STATS] BPLKPAHEAD"), "{}", on.stderr);
    }
    std::fs::remove_dir_all(&dir).ok();
}
