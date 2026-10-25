//! The linking groups of TBLASTX as sub-slices of the sorted HSP list
//! (`src/algorithm/tblastx/x_link_slices.rs`) on whole searches.
//!
//! The linking call of the in-place path (`x_hsp_inplace.rs`) links its groups where they lie in
//! the sorted list instead of cutting them into Vecs and copying the linked groups into a new
//! list. Every test reads the `[X_STATS] LINKSLICES` line to check that the path ran, not only
//! that the report is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
//! ```c
//! if (H->prev != NULL &&
//!     ((H->hsp->context/strand_factor) !=
//!      (H->prev->hsp->context/strand_factor) ||
//!      (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
//! { /* If frame switches, then start new list. */
//!    hp_frame_number[cur_frame]--;
//!    hp_frame_number[++cur_frame]++;
//!    hp_frame_start[cur_frame] = H;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
//! ```c
//! for (frame_index=0; frame_index<num_query_frames; frame_index++)
//! {
//!    hp_start = s_LinkHSPStructReset(hp_start);
//!    hp_start->next = hp_frame_start[frame_index];
//!    hp_frame_start[frame_index]->prev = hp_start;
//!    number_of_hsps = hp_frame_number[frame_index];
//! ```
//! The expected report (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was
//! produced by NCBI BLAST+ 2.17.0.
//!
//! The fixture pair has one query and one subject: each run makes the subject's two linking
//! calls, so a run with the path on counts two linkings.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 2] = [1, 4];

/// One TBLASTX run of the fixture pair (`-outfmt 6` on standard output) with the given
/// variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing. The run must exit with 0. Returns the report and
/// the standard error.
fn tblastx(threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("tblastx")
        .arg("-query")
        .arg(fixture("short_hsp.query.fna"))
        .arg("-subject")
        .arg(fixture("short_hsp.subject.fna"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

fn ncbi_report() -> Vec<u8> {
    std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output")
}

/// The counters `names` of the line that starts with `prefix`, or `None` without such a line.
fn counters<const N: usize>(stderr: &str, prefix: &str, names: [&str; N]) -> Option<[u64; N]> {
    let line = stderr.lines().find(|line| line.starts_with(prefix))?;
    Some(names.map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("counter")
    }))
}

/// `[X_STATS] LINKSLICES`: linkings, groups, parallel_linkings, legacy_fallbacks.
fn stats(stderr: &str) -> Option<[u64; 4]> {
    counters(
        stderr,
        "[X_STATS] LINKSLICES ",
        [
            "linkings",
            "groups",
            "parallel_linkings",
            "legacy_fallbacks",
        ],
    )
}

/// `[X_STATS] HSPINPLACE`: linkings, hits, groups.
fn in_place_stats(stderr: &str) -> Option<[u64; 3]> {
    counters(
        stderr,
        "[X_STATS] HSPINPLACE ",
        ["linkings", "hits", "groups"],
    )
}

/// `[X_STATS] PAIRPAR`: units. At one thread the cost rule sends the fixture pair to the reference
/// frame loop (`units=0`), at four threads it runs its units.
fn check_path(stderr: &str, threads: u32, what: &str) {
    let [units] = counters(stderr, "[X_STATS] PAIRPAR ", ["units"])
        .unwrap_or_else(|| panic!("{what}: no [X_STATS] PAIRPAR line:\n{stderr}"));
    if threads == 1 {
        assert_eq!(units, 0, "{what}: the frame loop:\n{stderr}");
    } else {
        assert!(units > 0, "{what}: the unit path:\n{stderr}");
    }
}

/// The checks every run shares: the stats line is there, the subject's two
/// linking calls went through the path, every group of the in-place path was linked as a slice
/// (the groups of the `HSPINPLACE` line of the same run), no call went back to the cut, and the
/// pool was used only with more than one thread. Returns (linkings, groups).
fn check_stats(stderr: &str, threads: u32, what: &str) -> (u64, u64) {
    let [linkings, groups, parallel_linkings, legacy_fallbacks] =
        stats(stderr).unwrap_or_else(|| panic!("{what}: no [X_STATS] LINKSLICES line:\n{stderr}"));
    check_path(stderr, threads, what);
    assert_eq!(linkings, 2, "{what}:\n{stderr}");
    assert!(groups >= linkings, "{what}:\n{stderr}");
    assert_eq!(legacy_fallbacks, 0, "{what}:\n{stderr}");
    assert!(parallel_linkings <= linkings, "{what}:\n{stderr}");
    if threads == 1 {
        assert_eq!(parallel_linkings, 0, "{what}:\n{stderr}");
    }
    if let Some([in_place_linkings, _, in_place_groups]) = in_place_stats(stderr) {
        assert_eq!(in_place_linkings, linkings, "{what}:\n{stderr}");
        assert_eq!(in_place_groups, groups, "{what}:\n{stderr}");
    }
    (linkings, groups)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start = s_LinkHSPStructReset(hp_start);
//    hp_start->next = hp_frame_start[frame_index];
//    hp_frame_start[frame_index]->prev = hp_start;
//    number_of_hsps = hp_frame_number[frame_index];
// ```
// The default search: the report is NCBI's, both linking calls link their groups as slices, and
// the numbers are the same with one thread and four.
#[test]
fn the_default_search_gives_ncbis_report_with_the_groups_linked_as_slices() {
    let expected = ncbi_report();
    let mut seen = Vec::new();
    for threads in THREADS {
        let (report, stderr) = tblastx(threads, &[STATS]);
        assert_eq!(report, expected, "threads {threads}");
        let what = format!("threads {threads}");
        let (_, groups) = check_stats(&stderr, threads, &what);
        assert!(in_place_stats(&stderr).is_some(), "{what}:\n{stderr}");
        seen.push(groups);
    }
    assert_eq!(seen[0], seen[1], "the groups do not depend on the threads");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,483-486
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ...
// /* Sort by (reverse) position. */
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// The slices together with the other variables of the linking call: the literal port
// (LOSAT_LINKING_LEGACY: its groups are linked on copies of the slices and written back, never
// handed back to the cut) and the literal port beside the default kernel
// (LOSAT_LINK_FAST_SHADOW). Threads 1 and 4 cover the reference frame loop and the unit path. The
// report is NCBI's and the path ran in every run.
#[test]
fn the_path_with_the_other_linking_variables() {
    let expected = ncbi_report();
    let others: [(&str, &str); 2] = [
        ("LOSAT_LINKING_LEGACY", "1"),
        ("LOSAT_LINK_FAST_SHADOW", "1"),
    ];
    for threads in THREADS {
        for other in others {
            let what = format!("{other:?} threads {threads}");
            let (report, stderr) = tblastx(threads, &[STATS, other]);
            assert_eq!(report, expected, "{what}");
            check_stats(&stderr, threads, &what);
        }
    }
}

/// The `name value` pairs of a `[TIMING] x_link1:` line in milliseconds, up to the parenthesis,
/// and the line itself.
fn link_timing(stderr: &str) -> (Vec<(String, i64)>, String) {
    let line = stderr
        .lines()
        .find_map(|line| line.strip_prefix("[TIMING] x_link1: "))
        .unwrap_or_else(|| panic!("no [TIMING] x_link1 line:\n{stderr}"));
    let mut fields = line.split_whitespace();
    let mut pairs = Vec::new();
    while let (Some(name), Some(value)) = (fields.next(), fields.next()) {
        let Some(seconds) = value.strip_suffix('s') else {
            break;
        };
        let (whole, fraction) = seconds
            .split_once('.')
            .unwrap_or_else(|| panic!("no decimal point in {value}: {line}"));
        assert_eq!(fraction.len(), 3, "{line}");
        let millis = whole.parse::<i64>().expect("seconds") * 1000
            + fraction.parse::<i64>().expect("milliseconds");
        pairs.push((name.to_string(), millis));
    }
    (pairs, line.to_string())
}

/// The value of `name=<integer>` in the parenthesis that ends a linking line.
fn count(line: &str, name: &str) -> u64 {
    let key = format!("{name}=");
    let tail = line
        .split_once('(')
        .and_then(|(_, tail)| tail.split_once(key.as_str()))
        .unwrap_or_else(|| panic!("no {name}= in {line}"))
        .1;
    tail.split(|c: char| !c.is_ascii_digit())
        .next()
        .and_then(|digits| digits.parse().ok())
        .unwrap_or_else(|| panic!("no number after {name}= in {line}"))
}

// No NCBI counterpart (a timer). The `x_link1` line of LOSAT_TIMING is printed with its seven
// steps in their order, and the steps add up to its `sum` (seven values rounded to a
// millisecond against their rounded sum). No time is compared with a bound.
#[test]
fn the_timer_counts_the_steps_of_the_linking_call() {
    const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
    const STEPS: [&str; 7] = [
        "reset",
        "sort",
        "groups",
        "kernels",
        "assemble",
        "final_sorts",
        "replay",
    ];
    let expected = ncbi_report();
    for threads in THREADS {
        let what = format!("timing, threads {threads}");
        let (report, stderr) = tblastx(threads, &[TIMING]);
        assert_eq!(report, expected, "{what}");
        let (pairs, line) = link_timing(&stderr);
        let names: Vec<&str> = pairs.iter().map(|(name, _)| name.as_str()).collect();
        let expected_names: Vec<&str> = STEPS.iter().copied().chain(["sum"]).collect();
        assert_eq!(names, expected_names, "{what}: {line}");
        let (steps, sum) = pairs.split_at(STEPS.len());
        let step_sum: i64 = steps.iter().map(|(_, value)| value).sum();
        assert!(
            (step_sum - sum[0].1).abs() <= 4,
            "{what}: the steps add up to {step_sum} ms, sum is {} ms: {line}",
            sum[0].1
        );
        assert!(count(&line, "calls") >= 1, "{what}: {line}");
    }
}
