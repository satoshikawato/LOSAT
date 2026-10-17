//! blastp's composition-based redo with `LOSAT_X_BPREDOPAR` (`src/algorithm/blastp/x_redo_par.rs`)
//! on whole searches.
//!
//! The switch is not in the strict set (`src/utils/xswitch.rs`): it is off unless its variable is
//! set. With it on, more than one thread and at least two matches, the matches of the redo run on
//! the pool in jobs of up to four consecutive matches of one query, in windows of jobs, and the
//! early termination test and the heap insertions are replayed in the order of the one-thread
//! loop. These tests compare the reports with the switch off at one thread (the one-thread loop)
//! and with the switch off and on at one, two and four threads, for blastp and blastp-fast, alone
//! and with the block search of `LOSAT_X_BPQBLOCK`; a search with `-max_target_seqs 1`, where the
//! early termination test ends matches; one query; the release of the preliminary lists
//! (`LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`) and the per-worker SEG cache (`LOSAT_X_SEGSHARE=0`);
//! small windows (`LOSAT_X_BPREDOPARWINDOW=1`) that end inside queries and leave out the matches
//! the test has ended before them, at two and eight threads; and the shadow mode
//! (`LOSAT_X_BPREDOPARSHADOW`), also with small windows, ended matches and released lists. The
//! `[TIMING] x_bpredopar:` lines show that the path ran and what it did, not only that the report
//! is equal. `LOSAT_X_BPQBLOCK` is in the strict set, on by default with blocks of 50,000 query
//! letters (two blocks of the fixture's 77,448): the tests that count one search of all queries
//! set `LOSAT_X_BPQBLOCK=0`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3447-3449,3525-3529
//! ```c
//! #pragma omp for schedule(static)
//!         for (b = 0; b < numMatches; ++b) {
//! ...
//!                 if (BlastCompo_EarlyTermination(
//!                         localMatch->best_evalue,
//!                         redoneMatches,
//!                         numQueries
//!                 )) {
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The first protein of the query fixture, alone in a file.
fn one_query() -> &'static PathBuf {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).expect("query fixture");
        let second = text[1..].find("\n>").map_or(text.len(), |at| at + 2);
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("blastp_redo_par");
        std::fs::create_dir_all(&dir).expect("work directory");
        let file = dir.join("one_query.faa");
        std::fs::write(&file, &text[..second]).expect("query file");
        file
    })
}

/// The counts of one `[TIMING] x_bpredopar:` line, or their sums.
#[derive(Clone, Copy, Debug, Default)]
struct RedoPar {
    matches: usize,
    windows: usize,
    /// Windows that end inside a query.
    split: usize,
    /// Matches the early termination test ended ...
    ended: usize,
    /// ... before their window (not computed) ...
    skipped: usize,
    /// ... and at their replay (computed).
    wasted: usize,
}

impl RedoPar {
    fn parse(line: &str) -> Self {
        let count = |name: &str| -> usize {
            line.split(' ')
                .find_map(|field| field.strip_prefix(name))
                .and_then(|value| value.parse().ok())
                .unwrap_or_else(|| panic!("{name} in {line}"))
        };
        Self {
            matches: count("matches="),
            windows: count("windows="),
            split: count("split="),
            ended: count("ended="),
            skipped: count("skipped="),
            wasted: count("wasted="),
        }
    }

    fn sum(lines: &[Self]) -> Self {
        lines.iter().fold(Self::default(), |sum, line| Self {
            matches: sum.matches + line.matches,
            windows: sum.windows + line.windows,
            split: sum.split + line.split,
            ended: sum.ended + line.ended,
            skipped: sum.skipped + line.skipped,
            wasted: sum.wasted + line.wasted,
        })
    }
}

/// What one run printed: the report, the counts of each `[TIMING] x_bpredopar:` line (one per
/// search the path ran), the lists `LOSAT_X_BPPFREE` freed (`[TIMING] x_bppfree:`, summed), and
/// standard error.
struct Run {
    stdout: Vec<u8>,
    redo_par: Vec<RedoPar>,
    freed_lists: usize,
    stderr: String,
    /// With `LOSAT_X_BPREDOLPT`: the windows computed with it (`lpt_windows=` of the
    /// `[TIMING] x_bpredopar:` lines, summed); `None` when no line has the field.
    lpt_windows: Option<usize>,
}

/// One blastp run of `query` against the fixture subjects (two related nimavirus proteomes) with
/// `vars` as its only `LOSAT_*` variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays:
/// a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI)
/// needs that file name as soon as a search has more than one thread, and it selects nothing;
/// without the caller's value a file of the tests is named. `BATCH_SIZE` is the test's own (3000:
/// several query batches of several queries). The run must exit with 0.
fn blastp(query: &Path, threads: u32, vars: &[(&str, &str)], args: &[&str]) -> Run {
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
            std::env::temp_dir().join("losat-bpredopar-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-num_threads", &threads.to_string()])
        .args(args)
        .env("BATCH_SIZE", "3000")
        .env("LOSAT_TIMING", "1")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} {args:?} threads {threads}: {stderr}"
    );
    let redo_par = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bpredopar: "))
        .map(RedoPar::parse)
        .collect();
    let freed_lists = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bppfree: freed_lists="))
        .map(|count| count.trim().parse::<usize>().expect("a count"))
        .sum();
    let lpt_windows = stderr
        .lines()
        .filter(|line| line.starts_with("[TIMING] x_bpredopar: "))
        .filter_map(|line| {
            line.split(' ')
                .find_map(|field| field.strip_prefix("lpt_windows="))
                .map(|count| count.parse::<usize>().expect("a count"))
        })
        .reduce(|sum, count| sum + count);
    Run {
        stdout: output.stdout,
        redo_par,
        freed_lists,
        stderr,
        lpt_windows,
    }
}

const SWITCH: (&str, &str) = ("LOSAT_X_BPREDOPAR", "1");
/// Windows of one job per thread.
const SMALL_WINDOWS: (&str, &str) = ("LOSAT_X_BPREDOPARWINDOW", "1");
/// One search of all queries (the block search of `LOSAT_X_BPQBLOCK` off).
const ONE_SEARCH: (&str, &str) = ("LOSAT_X_BPQBLOCK", "0");

/// The reports with the switch off and on at 1, 2 and 4 threads equal the one-thread report; the
/// path runs (one line per search, at least two matches) at 2 and 4 threads only.
fn same_rows_at_every_thread_count(query: &Path, vars: &[(&str, &str)], args: &[&str]) {
    let reference = blastp(query, 1, vars, args);
    assert!(!reference.stdout.is_empty(), "{vars:?} {args:?}: a report");
    assert!(reference.redo_par.is_empty(), "no line without the switch");
    for threads in [1, 2, 4] {
        let off = blastp(query, threads, vars, args);
        assert!(
            off.stdout == reference.stdout,
            "{vars:?} {args:?} threads {threads}: the switch-off report differs"
        );
        assert!(off.redo_par.is_empty());
        let with_switch: Vec<(&str, &str)> = vars.iter().copied().chain([SWITCH]).collect();
        let on = blastp(query, threads, &with_switch, args);
        assert!(
            on.stdout == reference.stdout,
            "{vars:?} {args:?} threads {threads}: the report differs with the switch"
        );
        if threads == 1 {
            assert!(
                on.redo_par.is_empty(),
                "one thread: the one-thread loop runs"
            );
        } else {
            assert!(
                !on.redo_par.is_empty() && RedoPar::sum(&on.redo_par).matches > 1,
                "{vars:?} {args:?} threads {threads}: the path did not run: {}",
                on.stderr
            );
        }
        let with_lpt: Vec<(&str, &str)> = with_switch.iter().copied().chain([LPT]).collect();
        let lpt = blastp(query, threads, &with_lpt, args);
        assert!(
            lpt.stdout == reference.stdout,
            "{vars:?} {args:?} threads {threads}: the report differs with the jobs largest first"
        );
        assert_lpt_ran(
            &lpt,
            threads > 1,
            &format!("{vars:?} {args:?} threads {threads}"),
        );
    }
}

/// `LOSAT_X_BPREDOLPT` computed windows in the run (at more than one thread), or the redo path
/// did not run (one thread).
fn assert_lpt_ran(run: &Run, threads: bool, what: &str) {
    if threads {
        assert!(
            run.lpt_windows.is_some_and(|windows| windows > 0),
            "{what}: no window computed with the jobs largest first: {}",
            run.stderr
        );
    } else {
        assert!(run.redo_par.is_empty(), "{what}: one thread");
    }
}

#[test]
fn match_jobs_give_the_rows_of_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    for task in ["blastp", "blastp-fast"] {
        let args = ["-task", task, "-evalue", "0.001", "-outfmt", "6"];
        same_rows_at_every_thread_count(&query, &[], &args);
    }
}

#[test]
fn match_jobs_in_query_blocks_give_the_rows_of_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    for task in ["blastp", "blastp-fast"] {
        let args = ["-task", task, "-evalue", "0.001", "-outfmt", "6"];
        same_rows_at_every_thread_count(&query, &[("LOSAT_X_BPQBLOCK", "1")], &args);
    }
}

#[test]
fn match_jobs_of_one_query_give_its_rows() {
    let args = ["-evalue", "10", "-outfmt", "6"];
    same_rows_at_every_thread_count(one_query(), &[], &args);
}

/// With `-max_target_seqs 1` the early termination test ends matches of the last query of a
/// batch (and of a query alone in its batch) once the heaps of the batch hold a subject within
/// the cutoff; the path leaves them out before their window or drops them at their replay.
#[test]
fn early_termination_ends_the_matches_it_ends_in_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    for max_target_seqs in ["1", "5"] {
        let args = [
            "-evalue",
            "10",
            "-max_target_seqs",
            max_target_seqs,
            "-outfmt",
            "6",
        ];
        let reference = blastp(&query, 1, &[], &args);
        for (threads, lpt) in [(2, false), (4, false), (2, true), (4, true)] {
            let vars = if lpt { vec![SWITCH, LPT] } else { vec![SWITCH] };
            let on = blastp(&query, threads, &vars, &args);
            assert!(
                on.stdout == reference.stdout,
                "-max_target_seqs {max_target_seqs} threads {threads} {vars:?}: the report differs"
            );
            if lpt {
                assert_lpt_ran(&on, true, &format!("-max_target_seqs {max_target_seqs}"));
            }
            let ended = RedoPar::sum(&on.redo_par).ended;
            assert!(
                !on.redo_par.is_empty(),
                "-max_target_seqs {max_target_seqs} threads {threads}: the path did not run"
            );
            if max_target_seqs == "1" {
                assert!(
                    ended > 0,
                    "threads {threads}: no match ended: {:?}",
                    on.redo_par
                );
            }
        }
    }
}

#[test]
fn released_lists_the_per_worker_seg_cache_and_pairwise_reports_are_unchanged() {
    let query = fixture("MeenMJNV.faa");
    for (vars, args) in [
        (
            vec![("LOSAT_X_BPPFREE", "1"), ("LOSAT_X_BPPCOMPACT", "1")],
            ["-evalue", "0.001", "-outfmt", "6"],
        ),
        (
            vec![("LOSAT_X_SEGSHARE", "0")],
            ["-evalue", "0.001", "-outfmt", "6"],
        ),
        (vec![], ["-evalue", "0.001", "-outfmt", "0"]),
    ] {
        let reference = blastp(&query, 1, &vars, &args);
        let with_switch: Vec<(&str, &str)> = vars.iter().copied().chain([SWITCH]).collect();
        let on = blastp(&query, 4, &with_switch, &args);
        assert!(
            on.stdout == reference.stdout,
            "{vars:?} {args:?}: the report differs with the switch"
        );
        assert!(!on.redo_par.is_empty(), "{vars:?} {args:?}: the path ran");
        let with_lpt: Vec<(&str, &str)> = with_switch.iter().copied().chain([LPT]).collect();
        let lpt = blastp(&query, 4, &with_lpt, &args);
        assert!(
            lpt.stdout == reference.stdout,
            "{vars:?} {args:?}: the report differs with the jobs largest first"
        );
        assert_eq!(
            lpt.freed_lists, reference.freed_lists,
            "{vars:?}: lists freed"
        );
        assert_lpt_ran(&lpt, true, &format!("{vars:?} {args:?}"));
    }
}

/// The input of the small-window tests: `-evalue 10` gives the last queries of some batches (all
/// of whose queries hit) more matches than a window of one job per thread holds, so with
/// `-max_target_seqs 1` the test ends matches in later windows than the one that fills the heaps.
const SMALL_WINDOW_ARGS: [&str; 6] = ["-evalue", "10", "-max_target_seqs", "1", "-outfmt", "6"];
const RELEASE: [(&str, &str); 2] = [("LOSAT_X_BPPFREE", "1"), ("LOSAT_X_BPPCOMPACT", "1")];

/// Windows of one job per thread end inside queries; the matches that the test has ended before a
/// window are not computed; the preliminary lists are released after their last match. The
/// report, the lists freed and the matches ended are those of the one-thread loop.
#[test]
fn small_windows_split_queries_and_leave_out_the_matches_the_test_ended() {
    let query = fixture("MeenMJNV.faa");
    let reference = blastp(&query, 1, &RELEASE, &SMALL_WINDOW_ARGS);
    assert!(reference.freed_lists > 0, "{}", reference.stderr);
    let vars: Vec<(&str, &str)> = RELEASE
        .iter()
        .copied()
        .chain([SWITCH, SMALL_WINDOWS])
        .collect();
    let with_lpt: Vec<(&str, &str)> = vars.iter().copied().chain([LPT]).collect();
    let mut ended = Vec::new();
    for (threads, vars) in [(2, &vars), (8, &vars), (2, &with_lpt), (8, &with_lpt)] {
        let on = blastp(&query, threads, vars, &SMALL_WINDOW_ARGS);
        assert!(
            on.stdout == reference.stdout,
            "threads {threads}: the report differs"
        );
        assert_eq!(
            on.freed_lists, reference.freed_lists,
            "threads {threads}: lists freed"
        );
        let counts = RedoPar::sum(&on.redo_par);
        // Ended before their window (not computed) and at their replay (computed and dropped).
        assert!(
            counts.windows > 1
                && counts.split > 0
                && counts.ended > 0
                && counts.skipped > 0
                && counts.wasted > 0,
            "threads {threads}: {counts:?}"
        );
        assert_eq!(counts.ended, counts.skipped + counts.wasted, "{counts:?}");
        if vars.contains(&LPT) {
            assert_lpt_ran(&on, true, &format!("small windows, threads {threads}"));
        }
        ended.push(counts.ended);
    }
    assert!(
        ended.iter().all(|&count| count == ended[0]),
        "the matches ended at 2 and 8 threads, with and without the jobs largest first: {ended:?}"
    );
}

/// The counts of the shadow mode's line at exit.
fn shadow_counts(stderr: &str) -> [u64; 5] {
    let line = stderr
        .lines()
        .find_map(|line| line.strip_prefix("[X_SHADOW] BPREDOPAR "))
        .unwrap_or_else(|| panic!("no shadow line: {stderr}"));
    let mut counts = [0u64; 5];
    for (count, name) in counts.iter_mut().zip([
        "searches=",
        "queries=",
        "matches_compared=",
        "lists_compared=",
        "diffs=",
    ]) {
        *count = line
            .split(' ')
            .find_map(|field| field.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("{name} in {line}"));
    }
    counts
}

#[test]
fn the_shadow_mode_compares_every_redone_match_and_hit_list() {
    let query = fixture("MeenMJNV.faa");
    let args = ["-evalue", "0.001", "-outfmt", "6"];
    // One search of all queries, the default blocks of 50,000 letters, and blocks of one batch.
    for (vars, one_search) in [
        (vec![ONE_SEARCH], true),
        (vec![], false),
        (vec![("LOSAT_X_BPQBLOCK", "1")], false),
    ] {
        let reference = blastp(&query, 1, &vars, &args);
        let with_shadow: Vec<(&str, &str)> = vars
            .iter()
            .copied()
            .chain([("LOSAT_X_BPREDOPARSHADOW", "1")])
            .collect();
        let shadow = blastp(&query, 4, &with_shadow, &args);
        assert!(
            shadow.stdout == reference.stdout,
            "{vars:?}: the report differs in the shadow mode"
        );
        let [searches, queries, matches, lists, diffs] = shadow_counts(&shadow.stderr);
        assert!(
            searches > 0 && queries > 0 && matches > 0 && lists > 0 && diffs == 0,
            "{vars:?}: {}",
            shadow.stderr
        );
        if one_search {
            assert_eq!(searches, 1, "one search of all queries");
        } else {
            assert!(searches > 1, "{vars:?}: one search per block");
        }
    }
}

/// The shadow mode where the replay ends matches (before their window and at their replay), the
/// windows end inside queries and the one-thread loop releases its lists (the shadow keeps them
/// for its own path).
#[test]
fn the_shadow_mode_compares_with_ended_matches_released_lists_and_small_windows() {
    let query = fixture("MeenMJNV.faa");
    let one_search: Vec<(&str, &str)> = RELEASE.iter().copied().chain([ONE_SEARCH]).collect();
    let reference = blastp(&query, 1, &one_search, &SMALL_WINDOW_ARGS);
    let vars: Vec<(&str, &str)> = one_search
        .iter()
        .copied()
        .chain([("LOSAT_X_BPREDOPARSHADOW", "1"), SMALL_WINDOWS])
        .collect();
    let shadow = blastp(&query, 4, &vars, &SMALL_WINDOW_ARGS);
    assert!(
        shadow.stdout == reference.stdout,
        "the report differs in the shadow mode"
    );
    let [searches, queries, matches, lists, diffs] = shadow_counts(&shadow.stderr);
    assert!(
        searches == 1 && queries > 0 && matches > 0 && lists > 0 && diffs == 0,
        "{}",
        shadow.stderr
    );
    // Every match that the replay redoes is compared; the others were ended.
    let counts = RedoPar::sum(&shadow.redo_par);
    assert!(
        counts.split > 0 && counts.skipped > 0,
        "the shadow's path: {counts:?}"
    );
    assert_eq!(
        counts.matches - counts.ended,
        matches as usize,
        "{counts:?}"
    );
    assert_eq!(shadow.freed_lists, reference.freed_lists, "lists freed");
}

/// `LOSAT_X_BPREDOLPT`: each window's jobs in one `rayon::scope` of one task per thread, taken
/// one at a time, largest first (`src/algorithm/blastp/x_redo_par.rs`, `x_compute_window_lpt`).
/// It acts only where `LOSAT_X_BPREDOPAR` (or its shadow mode) runs. The tests above run it beside
/// the switch; the test below runs it in the shadow mode.
const LPT: (&str, &str) = ("LOSAT_X_BPREDOLPT", "1");

/// The shadow mode checks the windows computed with `LOSAT_X_BPREDOLPT` (every match the replay
/// redoes against the one-thread computation, every hit list against the one-thread loop's), and
/// its line ends with `lpt=1`: one search or blocks of queries, and small windows that end inside
/// queries with matches ended before their window and at their replay and released lists.
#[test]
fn the_shadow_mode_checks_the_windows_computed_largest_first() {
    let query = fixture("MeenMJNV.faa");
    let args = ["-evalue", "0.001", "-outfmt", "6"];
    let shadow_lpt = [("LOSAT_X_BPREDOPARSHADOW", "1"), LPT];
    for (vars, args, small) in [
        (vec![], &args[..], false),
        (vec![("LOSAT_X_BPQBLOCK", "1")], &args[..], false),
        (RELEASE.to_vec(), &SMALL_WINDOW_ARGS[..], true),
    ] {
        let reference = blastp(&query, 1, &vars, args);
        let mut with_shadow: Vec<(&str, &str)> = vars.iter().copied().chain(shadow_lpt).collect();
        if small {
            with_shadow.push(SMALL_WINDOWS);
        }
        let shadow = blastp(&query, 4, &with_shadow, args);
        assert!(
            shadow.stdout == reference.stdout,
            "{with_shadow:?}: the report differs in the shadow mode"
        );
        let [searches, queries, matches, lists, diffs] = shadow_counts(&shadow.stderr);
        assert!(
            searches > 0 && queries > 0 && matches > 0 && lists > 0 && diffs == 0,
            "{with_shadow:?}: {}",
            shadow.stderr
        );
        assert!(
            shadow
                .stderr
                .lines()
                .any(|line| line.starts_with("[X_SHADOW] BPREDOPAR ") && line.ends_with(" lpt=1")),
            "{with_shadow:?}: the shadow line says lpt=1: {}",
            shadow.stderr
        );
        assert_lpt_ran(&shadow, true, &format!("{with_shadow:?}"));
        let counts = RedoPar::sum(&shadow.redo_par);
        assert_eq!(
            counts.matches - counts.ended,
            matches as usize,
            "{with_shadow:?}: every match the replay redoes is compared: {counts:?}"
        );
        if small {
            assert!(
                counts.split > 0 && counts.skipped > 0,
                "the shadow's path: {counts:?}"
            );
            assert_eq!(shadow.freed_lists, reference.freed_lists, "lists freed");
        }
    }
}
