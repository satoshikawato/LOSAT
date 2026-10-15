//! The index sorts of the TBLASTX HSP list with `LOSAT_X_LINKSORT`
//! (`src/algorithm/tblastx/x_link_sort.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. With it on, the linking call of the
//! in-place path (`x_hsp_inplace.rs`) sorts its output order once where it sorts twice, and an
//! index sort of 32,768 indices or more runs as a stable parallel sort when its call may use the
//! search pool. `LOSAT_X_LINKSORTSHADOW` runs
//! every index sort both ways and stops at the first difference. Every test reads the
//! `[X_STATS] LINKSORT` or `[X_SHADOW] LINKSORT` line to check that the path ran, not only that
//! the report is equal. Both lines are printed only while the switch (or its shadow) is on.
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
//! tests cover is the dispatch at the two sort sites, the single final sort, the shadow and the
//! summary lines, at one thread and on a pool of four.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

const SWITCH: (&str, &str) = ("LOSAT_X_LINKSORT", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_LINKSORT", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const IN_PLACE: (&str, &str) = ("LOSAT_X_HSPINPLACE", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_LINKSORTSHADOW", "1");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 2] = [1, 4];

/// One TBLASTX run of the fixture pair (`-outfmt 6` on standard output) with the given
/// variables as its only `LOSAT_*` variables that choose a code path. `LOSAT_X_DEFAULT` of the
/// test's own environment is dropped as well: every test states the default it runs with.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing. The run must exit with 0. Returns the report and
/// the standard error.
fn tblastx(threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// `[X_SHADOW] LINKSORT`: sorts, indices_compared.
fn shadow(stderr: &str) -> Option<[u64; 2]> {
    counters(
        stderr,
        "[X_SHADOW] LINKSORT ",
        ["sorts", "indices_compared"],
    )
}

/// `[X_SHADOW] HSPINPLACE`: sorts, linkings.
fn in_place_shadow(stderr: &str) -> Option<[u64; 2]> {
    counters(stderr, "[X_SHADOW] HSPINPLACE ", ["sorts", "linkings"])
}

/// The checks every run with the path on shares: the stats line is there, index sorts ran,
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
// The default search (no variable for the switch): the report is NCBI's, the subject's two
// linking calls each sort their output order once, and the sorts are the same with one thread
// and four.
#[test]
fn the_switch_gives_ncbis_report_with_one_final_sort_per_linking_call() {
    let expected = ncbi_report();
    let mut seen = Vec::new();
    for threads in THREADS {
        let (report, stderr) = tblastx(threads, &[STATS]);
        assert_eq!(report, expected, "threads {threads}");
        assert_eq!(shadow(&stderr), None, "threads {threads}:\n{stderr}");
        seen.push(check_stats(&stderr, &format!("threads {threads}")));
    }
    assert_eq!(seen[0], seen[1], "the sorts do not depend on the threads");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
// ```c
// /* Sort by (reverse) position. */
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// The shadow stays opt-in (set to any value) and wins over the switch's variable, whatever it
// says. It compares every index sort with today's statements: the run ends without a
// difference, the report is NCBI's, and the shadow line counts every sort the path made (the
// stats line of the same run).
#[test]
fn the_shadow_compares_every_index_sort() {
    let expected = ncbi_report();
    for threads in THREADS {
        let (report, stderr) = tblastx(threads, &[SHADOW]);
        assert_eq!(report, expected, "threads {threads}");
        let [sorts, _] = shadow(&stderr)
            .unwrap_or_else(|| panic!("threads {threads}: no [X_SHADOW] LINKSORT line:\n{stderr}"));
        assert!(sorts > 0, "threads {threads}:\n{stderr}");
        assert_eq!(stats(&stderr), None, "threads {threads}:\n{stderr}");

        let (report, stderr) = tblastx(threads, &[SHADOW, STATS]);
        assert_eq!(report, expected, "threads {threads}");
        let what = format!("shadow, threads {threads}");
        let (sorts, indices) = check_stats(&stderr, &what);
        assert_eq!(shadow(&stderr), Some([sorts, indices]), "{what}:\n{stderr}");

        // The variable of the switch changes nothing next to its shadow.
        let (report, stderr) = tblastx(threads, &[SWITCH, SHADOW, STATS]);
        assert_eq!(report, expected, "threads {threads}");
        assert_eq!(shadow(&stderr), Some([sorts, indices]), "{what}:\n{stderr}");

        // The shadow wins over `LOSAT_X_LINKSORT=0` and runs under `LOSAT_X_DEFAULT=0` with the
        // in-place path set, where the switch alone would be off.
        let (report, stderr) = tblastx(threads, &[SWITCH_OFF, SHADOW, STATS]);
        assert_eq!(report, expected, "threads {threads}");
        assert_eq!(shadow(&stderr), Some([sorts, indices]), "{what}:\n{stderr}");
        let what = format!("shadow, LOSAT_X_DEFAULT=0 with the in-place path, threads {threads}");
        let (report, stderr) = tblastx(threads, &[DEFAULT_OFF, IN_PLACE, SHADOW, STATS]);
        assert_eq!(report, expected, "{what}");
        check_stats(&stderr, &what);
        assert!(shadow(&stderr).is_some(), "{what}:\n{stderr}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// `LOSAT_X_HSPINPLACESHADOW` runs every sort call and every linking call of the in-place path
// also through the reference bodies and compares the lists record by record: with the switch
// on, that compares this path's sorts with the reference bodies end to end. The reference
// linking bodies sort through the same dispatch, so more index sorts are counted than without
// that shadow. Both shadows end without a difference and the report is NCBI's.
#[test]
fn the_shadow_of_the_in_place_path_accepts_the_sorts() {
    let expected = ncbi_report();
    const IN_PLACE_SHADOW: (&str, &str) = ("LOSAT_X_HSPINPLACESHADOW", "1");
    for threads in THREADS {
        let (_, alone) = tblastx(threads, &[STATS]);
        let (sorts_alone, _) = check_stats(&alone, &format!("alone, threads {threads}"));
        for path in [&[][..], &[SHADOW][..]] {
            let what = format!("{path:?} with the in-place shadow, threads {threads}");
            let switches = [&[STATS, IN_PLACE_SHADOW][..], path].concat();
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            let (sorts, indices) = check_stats(&stderr, &what);
            assert!(sorts > sorts_alone, "{what}:\n{stderr}");
            let [compared_sorts, compared_linkings] = in_place_shadow(&stderr)
                .unwrap_or_else(|| panic!("{what}: no [X_SHADOW] HSPINPLACE line:\n{stderr}"));
            assert!(compared_sorts > 0, "{what}:\n{stderr}");
            assert_eq!(compared_linkings, 2, "{what}:\n{stderr}");
            if !path.is_empty() {
                assert_eq!(shadow(&stderr), Some([sorts, indices]), "{what}:\n{stderr}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
// ```c
// if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
// ```
// `LOSAT_X_PAIRPAR=0` leaves the seed stage to the reference frame loop, which appends and
// sorts the list of each subject frame through the same dispatch; the linking calls are the
// in-place path's as before. The report is NCBI's with the switch (on by default, no variable)
// and with its shadow.
#[test]
fn the_frame_loop_sorts_through_the_same_dispatch() {
    let expected = ncbi_report();
    const FRAME_LOOP: (&str, &str) = ("LOSAT_X_PAIRPAR", "0");
    for threads in THREADS {
        for path in [&[][..], &[SHADOW][..]] {
            let what = format!("{path:?} {FRAME_LOOP:?} threads {threads}");
            let switches = [&[STATS, FRAME_LOOP][..], path].concat();
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            let (sorts, indices) = check_stats(&stderr, &what);
            if !path.is_empty() {
                assert_eq!(shadow(&stderr), Some([sorts, indices]), "{what}:\n{stderr}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTransl);
//    qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_FwdCompareHSPsTransl);
// ```
// `LOSAT_X_DEFAULT=0` (with the switch set: the run that the variable gave before the switch
// joined the set) and `LOSAT_X_HSPINPLACE=0` (the switch on by default) leave every sort and
// every linking call to the reference bodies, which this switch does not touch: it changes
// nothing there, and its lines show that no index sort went through it (every counter 0). The
// report is NCBI's.
#[test]
fn the_reference_bodies_do_not_reach_the_switch() {
    let expected = ncbi_report();
    let references: [&[(&str, &str)]; 2] = [&[DEFAULT_OFF, SWITCH], &[("LOSAT_X_HSPINPLACE", "0")]];
    for reference in references {
        for threads in THREADS {
            let what = format!("{reference:?} threads {threads}");
            let switches = [&[STATS][..], reference].concat();
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            assert_eq!(stats(&stderr), Some([0; 5]), "{what}:\n{stderr}");
            assert_eq!(shadow(&stderr), None, "{what}:\n{stderr}");

            let switches = [&[SHADOW, STATS][..], reference].concat();
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            assert_eq!(stats(&stderr), Some([0; 5]), "{what}:\n{stderr}");
            assert_eq!(shadow(&stderr), Some([0; 2]), "{what}:\n{stderr}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
// ```c
// /* Sort by (reverse) position. */
// if (kTranslatedQuery) {
//    qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//          s_RevCompareHSPsTbx);
// ```
// With `LOSAT_X_LINKSORT=0`, or with `LOSAT_X_DEFAULT=0` and no variable for the switch (the
// in-place path set, so that its sort sites are reached), the statements of the in-place path
// run as before the switch existed: no LINKSORT line is printed, and the report is NCBI's.
// `LOSAT_X_LINKSORT=0` wins over the default, and `LOSAT_X_LINKSORT=1` over `LOSAT_X_DEFAULT=0`
// (the run with the in-place path set that follows).
#[test]
fn the_switch_follows_its_variable_and_the_global_one() {
    let expected = ncbi_report();
    let offs: [&[(&str, &str)]; 4] = [
        &[SWITCH_OFF],
        &[DEFAULT_OFF],
        &[DEFAULT_OFF, IN_PLACE],
        &[DEFAULT_OFF, IN_PLACE, SWITCH_OFF],
    ];
    for off in offs {
        for threads in THREADS {
            let what = format!("{off:?} threads {threads}");
            let switches = [&[STATS][..], off].concat();
            let (report, stderr) = tblastx(threads, &switches);
            assert_eq!(report, expected, "{what}");
            assert_eq!(stats(&stderr), None, "{what}:\n{stderr}");
            assert_eq!(shadow(&stderr), None, "{what}:\n{stderr}");
        }
    }
    // The run that the variable gave before the switch joined the set: the path runs.
    for threads in THREADS {
        let what = format!("LOSAT_X_DEFAULT=0 with the variable, threads {threads}");
        let (report, stderr) = tblastx(threads, &[DEFAULT_OFF, IN_PLACE, SWITCH, STATS]);
        assert_eq!(report, expected, "{what}");
        check_stats(&stderr, &what);
        assert_eq!(shadow(&stderr), None, "{what}:\n{stderr}");
    }
}
