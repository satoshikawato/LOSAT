//! blastp's composition-based redo with `LOSAT_X_BPREDOPAR` (`src/algorithm/blastp/x_redo_par.rs`)
//! on whole searches.
//!
//! The switch is not in the strict set (`src/utils/xswitch.rs`): it is off unless its variable is
//! set. With it on, more than one thread and at least two matches, the matches of the redo run on
//! the pool in jobs of up to four consecutive matches of one query, and the early termination test and the heap insertions are replayed in the
//! order of the one-thread loop. These tests compare the reports with the switch off at one thread
//! (the one-thread loop) and with the switch off and on at one, two and four threads, for blastp
//! and blastp-fast, alone and with the block search of `LOSAT_X_BPQBLOCK`; a search with
//! `-max_target_seqs 1`, where the early termination test ends matches; one query; the release of
//! the preliminary lists (`LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`) and the per-worker SEG cache
//! (`LOSAT_X_SEGSHARE=0`); and the shadow mode (`LOSAT_X_BPREDOPARSHADOW`). The
//! `[TIMING] x_bpredopar:` lines show that the path ran, not only that the report is equal.
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

/// What one run printed: the report, and the (matches, windows, ended) of each
/// `[TIMING] x_bpredopar` line (one per search the path ran), and standard error.
struct Run {
    stdout: Vec<u8>,
    redo_par: Vec<(usize, usize, usize)>,
    stderr: String,
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
        .map(|line| {
            let count = |name: &str| -> usize {
                line.split(' ')
                    .find_map(|field| field.strip_prefix(name))
                    .and_then(|value| value.parse().ok())
                    .unwrap_or_else(|| panic!("{name} in {line}"))
            };
            (count("matches="), count("windows="), count("ended="))
        })
        .collect();
    Run {
        stdout: output.stdout,
        redo_par,
        stderr,
    }
}

const SWITCH: (&str, &str) = ("LOSAT_X_BPREDOPAR", "1");

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
                !on.redo_par.is_empty()
                    && on
                        .redo_par
                        .iter()
                        .map(|(matches, ..)| matches)
                        .sum::<usize>()
                        > 1,
                "{vars:?} {args:?} threads {threads}: the path did not run: {}",
                on.stderr
            );
        }
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
/// the cutoff; the path computes them and the replay drops them.
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
        for threads in [2, 4] {
            let on = blastp(&query, threads, &[SWITCH], &args);
            assert!(
                on.stdout == reference.stdout,
                "-max_target_seqs {max_target_seqs} threads {threads}: the report differs"
            );
            let ended: usize = on.redo_par.iter().map(|(_, _, ended)| ended).sum();
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
    }
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
    for vars in [vec![], vec![("LOSAT_X_BPQBLOCK", "1")]] {
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
        if vars.is_empty() {
            assert_eq!(searches, 1, "one search of all queries");
        } else {
            assert!(searches > 1, "one search per block");
        }
    }
}
