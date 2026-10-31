//! blastp's composition-based redo with `LOSAT_X_BPREDOPAR` (`src/algorithm/blastp/x_redo_par.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`; on unless `=0`); the searches here
//! start from `PREVIOUS_DEFAULT`, which sets it to `0`, and `SWITCH` turns it on. With it on, more
//! than one thread and at least two matches, the matches of the redo run on
//! the pool in jobs of up to four consecutive matches of one query, in windows of jobs, and the
//! early termination test and the heap insertions are replayed in the order of the one-thread
//! loop. These tests compare the reports with the switch off at one thread (the one-thread loop)
//! and with the switch on at one and four threads, for blastp in the default blocks of
//! `LOSAT_X_BPQBLOCK` and blastp-fast in blocks of one batch (there also with the switch off at
//! four threads); a search with
//! `-max_target_seqs 1`, where the early termination test ends matches; one query; the release of
//! the preliminary lists and pairwise output;
//! small windows (`LOSAT_X_BPREDOPARWINDOW=1`) that end inside queries and leave out the matches
//! the test has ended before them, at two and eight threads, and in one search of all queries with
//! the jobs largest first at four; and the windows that expand the compact lists as they reach
//! them against the one-thread search. The
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

/// The first protein of the query fixture, alone in a file.
fn one_query() -> &'static PathBuf {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).expect("query fixture");
        let second = text[1..].find("\n>").map_or(text.len(), |at| at + 2);
        // One directory per process: nextest runs each test in a process of its own, and a
        // process that wrote these files again while another one's search read them cut its input.
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join("blastp_redo_par")
            .join(std::process::id().to_string());
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
/// search the path ran), the lists the redo freed (`[TIMING] x_bppfree:`, summed), and
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
        .envs(PREVIOUS_DEFAULT.iter().copied())
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

/// The reports with the switch on at 1 and 4 threads, and with the jobs largest first at 4, equal
/// the one-thread report without it, which is returned; the path runs (one line per search, at
/// least two matches) at 4 threads only. The path runs at more than one thread
/// (`blastp_should_parallelize_stage` in `blast_engine.rs`); more threads change only the shape
/// of its windows (jobs per window, the first cap on jobs at risk, the tasks of the jobs largest
/// first), and 2 and 8 threads run in the tests of ended matches and small windows.
fn same_rows_at_one_and_four_threads(query: &Path, vars: &[(&str, &str)], args: &[&str]) -> Run {
    let reference = blastp(query, 1, vars, args);
    assert!(!reference.stdout.is_empty(), "{vars:?} {args:?}: a report");
    assert!(reference.redo_par.is_empty(), "no line without the switch");
    let with_switch: Vec<(&str, &str)> = vars.iter().copied().chain([SWITCH]).collect();
    for threads in [1, 4] {
        let on = blastp(query, threads, &with_switch, args);
        assert!(
            on.stdout == reference.stdout,
            "{vars:?} {args:?} threads {threads}: the report differs with the switch"
        );
        if threads == 1 {
            // `LOSAT_X_BPREDOLPT` acts only where the path runs, so it is not run here.
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
    }
    let with_lpt: Vec<(&str, &str)> = with_switch.iter().copied().chain([LPT]).collect();
    let lpt = blastp(query, 4, &with_lpt, args);
    assert!(
        lpt.stdout == reference.stdout,
        "{vars:?} {args:?}: the report differs with the jobs largest first"
    );
    assert_lpt_ran(&lpt, &format!("{vars:?} {args:?}"));
    reference
}

/// `LOSAT_X_BPREDOLPT` computed windows in the run (at more than one thread).
fn assert_lpt_ran(run: &Run, what: &str) {
    assert!(
        run.lpt_windows.is_some_and(|windows| windows > 0),
        "{what}: no window computed with the jobs largest first: {}",
        run.stderr
    );
}

// The task changes the matches, not the path (`x_redo_par.rs` has no task-specific code): one
// task per block plan.
#[test]
fn match_jobs_give_the_rows_of_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    let args = ["-task", "blastp", "-evalue", "0.001", "-outfmt", "6"];
    same_rows_at_one_and_four_threads(&query, &[], &args);
}

#[test]
fn match_jobs_in_query_blocks_give_the_rows_of_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    let blocks = [("LOSAT_X_BPQBLOCK", "1")];
    let args = ["-task", "blastp-fast", "-evalue", "0.001", "-outfmt", "6"];
    let reference = same_rows_at_one_and_four_threads(&query, &blocks, &args);
    // The switch off at 4 threads: the redo the switch replaces, the query-parallel redo in the
    // blocks of several queries and the match-parallel redo in the two blocks of one query
    // (`[losat-thread-stage]` stages kappa-query-redo and kappa-match-redo). The other blastp
    // tests run it with `LOSAT_X_BPREDOPAR=0`, mostly without the one-thread report.
    let off = blastp(&query, 4, &blocks, &args);
    assert!(
        off.stdout == reference.stdout,
        "threads 4: the switch-off report differs"
    );
    assert!(off.redo_par.is_empty(), "no line without the switch");
}

#[test]
fn match_jobs_of_one_query_give_its_rows() {
    let args = ["-evalue", "10", "-outfmt", "6"];
    same_rows_at_one_and_four_threads(one_query(), &[], &args);
}

/// With `-max_target_seqs 1` the early termination test ends matches of the last query of a
/// batch (and of a query alone in its batch) once the heaps of the batch hold a subject within
/// the cutoff; the path leaves them out before their window or drops them at their replay. With
/// `-max_target_seqs 1` at 2 threads and, with the jobs largest first, at 4; with 5 at 4 threads.
#[test]
fn early_termination_ends_the_matches_it_ends_in_the_one_thread_redo() {
    let query = fixture("MeenMJNV.faa");
    for (max_target_seqs, runs) in [
        ("1", &[(2, false), (4, true)][..]),
        ("5", &[(4, false)][..]),
    ] {
        let args = [
            "-evalue",
            "10",
            "-max_target_seqs",
            max_target_seqs,
            "-outfmt",
            "6",
        ];
        let reference = blastp(&query, 1, &[], &args);
        for &(threads, lpt) in runs {
            let vars = if lpt { vec![SWITCH, LPT] } else { vec![SWITCH] };
            let on = blastp(&query, threads, &vars, &args);
            assert!(
                on.stdout == reference.stdout,
                "-max_target_seqs {max_target_seqs} threads {threads} {vars:?}: the report differs"
            );
            if lpt {
                assert_lpt_ran(&on, &format!("-max_target_seqs {max_target_seqs}"));
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

/// The released lists with the jobs largest first; pairwise output with the windows computed by
/// `par_iter` (`LOSAT_X_BPREDOLPT` off). The lists are released in the replay, which does not
/// depend on how a window was computed.
#[test]
fn released_lists_and_pairwise_reports_are_unchanged() {
    let query = fixture("MeenMJNV.faa");
    for (vars, args, lpt) in [
        (vec![], ["-evalue", "0.001", "-outfmt", "6"], true),
        (vec![], ["-evalue", "0.001", "-outfmt", "0"], false),
    ] {
        let reference = blastp(&query, 1, &vars, &args);
        let mut with_switch: Vec<(&str, &str)> = vars.iter().copied().chain([SWITCH]).collect();
        if lpt {
            with_switch.push(LPT);
        }
        let on = blastp(&query, 4, &with_switch, &args);
        assert!(
            on.stdout == reference.stdout,
            "{with_switch:?} {args:?}: the report differs with the switch"
        );
        assert!(!on.redo_par.is_empty(), "{vars:?} {args:?}: the path ran");
        assert_eq!(
            on.freed_lists, reference.freed_lists,
            "{with_switch:?}: lists freed"
        );
        if lpt {
            assert_lpt_ran(&on, &format!("{vars:?} {args:?}"));
        }
    }
}

/// The input of the small-window tests: `-evalue 10` gives the last queries of some batches (all
/// of whose queries hit) more matches than a window of one job per thread holds, so with
/// `-max_target_seqs 1` the test ends matches in later windows than the one that fills the heaps.
const SMALL_WINDOW_ARGS: [&str; 6] = ["-evalue", "10", "-max_target_seqs", "1", "-outfmt", "6"];
/// The variables of the release tests: none, the lists are released by default.
const RELEASE: [(&str, &str); 0] = [];

/// Windows of one job per thread end inside queries; the matches that the test has ended before a
/// window are not computed; the preliminary lists are released after their last match. The
/// report, the lists freed and the matches ended are those of the one-thread loop, at 2 threads
/// and, with the jobs largest first, at 8 (windows of 2 and 8 jobs).
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
    for (threads, vars) in [(2, &vars), (8, &with_lpt)] {
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
            assert_lpt_ran(&on, &format!("small windows, threads {threads}"));
        }
        ended.push(counts.ended);
    }
    assert!(
        ended.iter().all(|&count| count == ended[0]),
        "the matches ended at 2 threads and at 8 with the jobs largest first: {ended:?}"
    );
}

/// Small windows in one search of all queries (the block search of `LOSAT_X_BPQBLOCK` off), with
/// the jobs largest first: the matches the replay ends (before their window and at their replay),
/// windows that end inside queries and released lists give the report and the lists freed of the
/// one-thread loop at 4 threads. The test above runs the default blocks of queries.
#[test]
fn small_windows_in_one_search_of_all_queries_give_the_one_thread_report() {
    let query = fixture("MeenMJNV.faa");
    let one_search: Vec<(&str, &str)> = RELEASE.iter().copied().chain([ONE_SEARCH]).collect();
    let reference = blastp(&query, 1, &one_search, &SMALL_WINDOW_ARGS);
    let vars: Vec<(&str, &str)> = one_search
        .iter()
        .copied()
        .chain([SWITCH, SMALL_WINDOWS, LPT])
        .collect();
    let on = blastp(&query, 4, &vars, &SMALL_WINDOW_ARGS);
    assert!(on.stdout == reference.stdout, "the report differs");
    let counts = RedoPar::sum(&on.redo_par);
    assert!(
        on.redo_par.len() == 1 && counts.split > 0 && counts.skipped > 0,
        "one search of all queries: {counts:?}"
    );
    assert_lpt_ran(&on, "small windows in one search");
    assert_eq!(on.freed_lists, reference.freed_lists, "lists freed");
}

// The windows expand the compact preliminary lists as they reach them instead of every list
// before the redo (`src/algorithm/blastp/x_redo_par.rs`, `XExpandAhead`; the variable
// `LOSAT_X_BPREDOWINEXPAND` graduated on 2026-10-27). The one-thread loop expands each query's list
// when it reaches it, so the reports and the lists freed of a windowed search are compared with
// the one-thread search.

/// The compact lists expanded in a run (`[TIMING] x_bppcompact: expanded_lists=`,
/// summed).
fn expanded_lists(stderr: &str) -> usize {
    stderr
        .lines()
        .filter_map(|line| line.strip_prefix("[TIMING] x_bppcompact: expanded_lists="))
        .map(|count| count.trim().parse::<usize>().expect("a count"))
        .sum()
}

/// The report, the lists expanded and freed of the windows that expand the lists as they reach
/// them equal those of the one-thread search, in small windows that end inside queries and in the
/// default windows of the default blocks.
#[test]
fn windows_that_expand_the_compact_lists_as_they_reach_them_give_the_one_thread_report() {
    let query = fixture("MeenMJNV.faa");
    let default_args = ["-evalue", "0.001", "-outfmt", "6"];
    // Small windows at 8 threads (windows of 8 jobs); the default windows at 2 (128 jobs: 5 and 3
    // windows in the two blocks, against 2 and 1 at 8 threads).
    for (vars, args, small, threads) in [
        (
            vec![SWITCH, SMALL_WINDOWS, ONE_SEARCH],
            &SMALL_WINDOW_ARGS[..],
            true,
            8,
        ),
        (vec![SWITCH], &default_args[..], false, 2),
    ] {
        let what = format!("{vars:?} threads {threads}");
        let reference = blastp(&query, 1, &vars, args);
        let on = blastp(&query, threads, &vars, args);
        assert!(!on.redo_par.is_empty(), "{what}: the path ran");
        if small {
            let counts = RedoPar::sum(&on.redo_par);
            assert!(
                counts.windows > 1 && counts.split > 0 && counts.skipped > 0,
                "{what}: {counts:?}"
            );
        }
        assert!(on.stdout == reference.stdout, "{what}: the report differs");
        assert_eq!(on.freed_lists, reference.freed_lists, "{what}: lists freed");
        // Each compact list is expanded once.
        assert_eq!(
            expanded_lists(&on.stderr),
            expanded_lists(&reference.stderr),
            "{what}: lists expanded"
        );
        assert!(
            expanded_lists(&reference.stderr) > 0,
            "{what}: {}",
            reference.stderr
        );
    }
}

/// `LOSAT_X_BPREDOLPT`: each window's jobs in one `rayon::scope` of one task per thread, taken
/// one at a time, largest first (`src/algorithm/blastp/x_redo_par.rs`, `x_compute_window_lpt`).
/// It acts only where `LOSAT_X_BPREDOPAR` runs. The tests above run it beside the switch.
const LPT: (&str, &str) = ("LOSAT_X_BPREDOLPT", "1");
