//! The order of the TBLASTX two-hit stage (`LOSAT_X_SEEDBUCKET`,
//! `src/algorithm/tblastx/x_seed_bucket.rs`) on a whole search, against the output of NCBI
//! BLAST+ 2.17.0.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): the bucketed order runs by default
//! where TBLASTX is in its scan loop with a diagonal table of at least 2^21 cells. The searches
//! of the fixture suites do not get there (the unit path of `LOSAT_X_PAIRPAR` takes them, and
//! their tables are smaller), so these tests bring a small search there: `LOSAT_X_PAIRPAR=0` for
//! the scan loop, `LOSAT_X_SEEDBUCKET_MIN_CELLS=0` for any table size, and
//! `LOSAT_X_SEEDBUCKET_CELLS=256` with `LOSAT_X_SEEDBUCKET_BUDGET=64` for several buckets and
//! several flushes per scan. Every test reads the `[X_STATS] SEEDBUCKET` line to check which
//! order ran, not only that the output is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
//! ```c
//! for (i = 0; i < hits; ++i) {
//!     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ...
//!     diag_coord = (query_offset - subject_offset) & diag_mask;
//! ```
//! NCBI handles the hits of a scan in scan order. The bucketed order handles them one range of
//! `diag_coord` after the other and must give the same report. The expected output
//! (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was produced by NCBI
//! BLAST+ 2.17.0, which runs the loop above.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

/// The variables that bring the fixture search to the scan loop with several buckets and
/// several flushes per scan, and print the counters. None of them names the switch itself.
const SCAN_LOOP: [(&str, &str); 5] = [
    ("LOSAT_X_PAIRPAR", "0"),
    ("LOSAT_X_SEEDBUCKET_MIN_CELLS", "0"),
    ("LOSAT_X_SEEDBUCKET_CELLS", "256"),
    ("LOSAT_X_SEEDBUCKET_BUDGET", "64"),
    ("LOSAT_X_STATS", "1"),
];

/// One TBLASTX run of the fixture pair (`-outfmt 6` on standard output) with `SCAN_LOOP` and
/// the given variables as its only `LOSAT_*` variables. `LOSAT_X_DEFAULT` of the test's own
/// environment is dropped as well: every test states the default it runs with. Returns the
/// report and the standard error.
fn tblastx(switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
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
        .envs(SCAN_LOOP)
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{switches:?}: {stderr}");
    (output.stdout, stderr)
}

fn ncbi_report() -> Vec<u8> {
    std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output")
}

/// One counter of the `[X_STATS] SEEDBUCKET` line.
fn counter(stderr: &str, name: &str) -> u64 {
    let key = format!("{name}=");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] SEEDBUCKET "))
        .unwrap_or_else(|| panic!("no [X_STATS] SEEDBUCKET line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|field| field.strip_prefix(key.as_str()))
        .unwrap_or_else(|| panic!("no {name} in {line}"))
        .parse()
        .expect("counter")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
// ```c
// for (i = 0; i < hits; ++i) {
//     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//     diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// With no variable for the switch the scans run in the bucketed order, more than one flush in
// at least one of them, and the report is NCBI's.
#[test]
fn the_bucketed_order_runs_by_default_and_matches_ncbi() {
    let (report, stderr) = tblastx(&[]);
    assert_eq!(report, ncbi_report(), "default");
    let scans = counter(&stderr, "bucketed_scans");
    assert!(scans > 0, "the bucketed order did not run:\n{stderr}");
    assert!(
        counter(&stderr, "flushes") > scans,
        "no scan flushed more than once:\n{stderr}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
// ```c
// for (i = 0; i < hits; ++i) {
//     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//     diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// `LOSAT_X_SEEDBUCKET=0` and `LOSAT_X_DEFAULT=0` each leave the scans to the port of this loop
// (scan order), and `LOSAT_X_SEEDBUCKET=1` under `LOSAT_X_DEFAULT=0` runs the bucketed order as
// it did before the switch joined the set. The report is NCBI's in each case.
#[test]
fn the_switch_follows_its_variable_and_the_global_one() {
    let expected = ncbi_report();
    for off in [("LOSAT_X_SEEDBUCKET", "0"), ("LOSAT_X_DEFAULT", "0")] {
        let (report, stderr) = tblastx(&[off]);
        assert_eq!(report, expected, "{off:?}");
        assert_eq!(counter(&stderr, "bucketed_scans"), 0, "{off:?}:\n{stderr}");
    }
    let as_before = [("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_SEEDBUCKET", "1")];
    let (report, stderr) = tblastx(&as_before);
    assert_eq!(report, expected, "{as_before:?}");
    assert!(
        counter(&stderr, "bucketed_scans") > 0,
        "{as_before:?}:\n{stderr}"
    );
}
