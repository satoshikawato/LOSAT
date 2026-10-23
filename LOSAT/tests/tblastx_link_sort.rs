//! The index sorts of the TBLASTX HSP list (`src/algorithm/tblastx/x_link_sort.rs`) on whole
//! searches.
//!
//! The linking call of the in-place path (`x_hsp_inplace.rs`) sorts its output order once where
//! the reference sorts twice, and an index sort of 32,768 indices or more runs as a stable
//! parallel sort when its call may use the search pool. Every test reads the `[X_STATS] LINKSORT`
//! line to check that the path ran, not only that the report is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
//! ```c
//! /* Sort by (reverse) position. */
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTbx);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
//! ```c
//! if (kTranslatedQuery) {
//!    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_RevCompareHSPsTransl);
//!    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//!          s_FwdCompareHSPsTransl);
//! ```
//! The expected report (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was
//! produced by NCBI BLAST+ 2.17.0.
//!
//! The fixture pair is small: every HSP list is far shorter than 32,768, so no sort of these
//! searches goes to the pool (`parallel_sorts=0` is asserted). The parallel sort is covered by
//! the unit tests of `x_link_sort.rs` (each comparator, lists below and above the threshold,
//! pools of 1, 2 and 8 threads) and by searches of large pairs outside this suite. What these
//! tests cover is the dispatch at the two sort sites, the single final sort and the summary
//! line, at one thread and on a pool of four.

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

/// `[X_STATS] LINKSORT`: sorts, indices, parallel_sorts, parallel_indices, composite_sorts.
fn stats(stderr: &str) -> Option<[u64; 5]> {
    counters(
        stderr,
        "[X_STATS] LINKSORT ",
        [
            "sorts",
            "indices",
            "parallel_sorts",
            "parallel_indices",
            "composite_sorts",
        ],
    )
}

/// `[X_STATS] PAIRPAR`: units, fallback_cost. At one thread the default search takes the reference
/// frame loop (the cost rule: `units=0 fallback_cost=1`), at four threads its units (`units>0`).
fn check_path(stderr: &str, threads: u32, what: &str) {
    let [units, fallback_cost] = counters(stderr, "[X_STATS] PAIRPAR ", ["units", "fallback_cost"])
        .unwrap_or_else(|| panic!("{what}: no [X_STATS] PAIRPAR line:\n{stderr}"));
    if threads == 1 {
        assert_eq!(
            (units, fallback_cost),
            (0, 1),
            "{what}: the frame loop:\n{stderr}"
        );
    } else {
        assert!(units > 0, "{what}: the unit path:\n{stderr}");
    }
}

/// The checks every run shares: the stats line is there, index sorts ran,
/// none went to the pool (the lists are short), and each of the subject's two linking calls
/// made its one final sort beside its first sort. Returns (sorts, indices).
fn check_stats(stderr: &str, what: &str) -> (u64, u64) {
    let [sorts, indices, parallel_sorts, parallel_indices, composite_sorts] =
        stats(stderr).unwrap_or_else(|| panic!("{what}: no [X_STATS] LINKSORT line:\n{stderr}"));
    assert_eq!(
        (parallel_sorts, parallel_indices),
        (0, 0),
        "{what}:\n{stderr}"
    );
    assert_eq!(composite_sorts, 2, "{what}:\n{stderr}");
    assert!(sorts >= 2 * composite_sorts, "{what}:\n{stderr}");
    assert!(indices >= 2 * sorts, "{what}:\n{stderr}");
    (sorts, indices)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// The default search: the report is NCBI's and the subject's two linking calls each sort their
// output order once. One thread runs the reference frame loop, four threads the unit path (the
// `[X_STATS] PAIRPAR` line tells which). The counters also count the index sorts of the seed
// stage's merge: with `LOSAT_X_TBXKMERGE=0` the unit path (four threads) makes the sorts of the
// reference frame loop (one thread) and the counts are equal; with that switch on (the default)
// its merge replaces some of them, so four threads make no more than one.
#[test]
fn the_default_search_gives_ncbis_report_with_one_final_sort_per_linking_call() {
    let expected = ncbi_report();
    for merge in [&[][..], &[("LOSAT_X_TBXKMERGE", "0")][..]] {
        let mut seen = Vec::new();
        for threads in THREADS {
            let what = format!("{merge:?} threads {threads}");
            let mut switches = vec![STATS];
            switches.extend_from_slice(merge);
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            check_path(&stderr, threads, &what);
            seen.push(check_stats(&stderr, &what));
        }
        if merge.is_empty() {
            assert!(
                seen[1].0 <= seen[0].0 && seen[1].1 <= seen[0].1,
                "the merge adds no sort: {seen:?}"
            );
        } else {
            assert_eq!(seen[0], seen[1], "the sorts do not depend on the threads");
        }
    }
}
