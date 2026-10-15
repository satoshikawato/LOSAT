//! The eight Newton lanes (`LOSAT_X_NEWTONLANES8`,
//! `src/core/composition_adjustment/x_newton_lanes.rs`) on whole blastp and tblastn searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. It acts on aarch64 only, and only
//! where the Newton lanes of `LOSAT_X_NEWTONLANES` (also in the set) run; on other targets it is
//! not read and the lanes stay four wide. It chooses how many independent problems share a
//! vector, so the property under test is that the report does not depend on it: the default run
//! (eight lanes on aarch64), `LOSAT_X_NEWTONLANES8=0` (four lanes) and `LOSAT_X_DEFAULT=0` (no
//! lanes) give the same bytes, and the shadow, which compares every stored result that is used
//! with the reference port, does not change them. The `[X_STATS] NEWTONLANES8` line shows that
//! the default run takes eight-lane batches on aarch64 and none elsewhere. On x86-64 the default
//! and `LOSAT_X_NEWTONLANES8=0` run the same path, so these tests mean most on the arm64 runners
//! of platform-checks. That the bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_NEWTONLANES8", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_NEWTONLANESSHADOW", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One search: the program, the fixture whose first records are the queries, how many of them,
/// and the subject fixture.
type Case = (&'static str, &'static str, usize, &'static str);

/// blastp of five proteins of one nimavirus against the proteome of another, and tblastn of 20
/// proteins of one nimavirus against the genome of another. Both run the composition-based redo
/// (the default `-comp_based_stats 2`) on many matches, so the Newton lanes get batches of five
/// or more problems (on x86-64 with `LOSAT_X_STATS=1`: blastp 451 problems in 30 batches,
/// tblastn 637 in 29).
const CASES: [Case; 2] = [
    ("blastp", "MeenMJNV.faa", 5, "MejoMJNV.faa"),
    ("tblastn", "PemoMJNVA.faa", 20, "MelaMJNV.fasta"),
];

/// A query file of the first records of a fixture, removed when dropped.
struct Queries(PathBuf);

impl Queries {
    fn new((program, fixture_name, records, _): Case, test: &str) -> Self {
        let text = std::fs::read_to_string(fixture(fixture_name)).expect("FASTA");
        let mut kept = String::new();
        let mut seen = 0;
        for line in text.lines() {
            if line.starts_with('>') {
                seen += 1;
            }
            if seen > records {
                break;
            }
            kept.push_str(line);
            kept.push('\n');
        }
        assert!(seen >= records, "{fixture_name}: {seen} records");
        let path = std::env::temp_dir().join(format!(
            "losat_newtonlanes8_{}_{test}_{program}.faa",
            std::process::id()
        ));
        std::fs::write(&path, kept).expect("write FASTA");
        Self(path)
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// One search (`-outfmt 6` on standard output) with the given variables as its only `LOSAT_*`
/// variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing; without the caller's value
/// a file of the tests is named. `BATCH_SIZE`, `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` are the
/// tests' own. The run must exit with 0. Returns the report and the standard error.
fn search(
    (program, _, _, subject): Case,
    queries: &Queries,
    threads: u32,
    vars: &[(&str, &str)],
) -> (Vec<u8>, String) {
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
            std::env::temp_dir().join("losat-newtonlanes8-worker.log"),
        );
    }
    let output = command
        .arg(program)
        .arg("-query")
        .arg(&queries.0)
        .arg("-subject")
        .arg(fixture(subject))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{program} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The number after `key=` on the standard-error line that starts with `prefix`, or `None` if
/// there is no such line.
fn counter(stderr: &str, prefix: &str, key: &str) -> Option<u64> {
    let line = stderr.lines().find(|line| line.starts_with(prefix))?;
    let value = line
        .split_whitespace()
        .find_map(|field| field.strip_prefix(key)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no {key}= in {line}"));
    Some(value.parse().expect("a counter"))
}

// The default search (eight lanes on aarch64, four elsewhere), `LOSAT_X_NEWTONLANES8=0` (four
// lanes) and `LOSAT_X_DEFAULT=0` (no lanes) give the same report, and none prints a line on
// standard error.
#[test]
fn the_report_does_not_depend_on_the_lane_width() {
    for case in CASES {
        let queries = Queries::new(case, "width");
        for threads in [1, 4] {
            let (default, default_stderr) = search(case, &queries, threads, &[]);
            assert!(!default.is_empty(), "{case:?} threads {threads}: no hits");
            assert_eq!(default_stderr, "", "{case:?} threads {threads}");
            for off in [SWITCH_OFF, DEFAULT_OFF] {
                let (report, stderr) = search(case, &queries, threads, &[off]);
                assert_eq!(report, default, "{case:?} {off:?} threads {threads}");
                assert_eq!(stderr, "", "{case:?} {off:?} threads {threads}");
            }
        }
    }
}

// With `LOSAT_X_STATS=1` the lanes report their batches. The default run takes eight-lane
// batches on aarch64 and prints no `NEWTONLANES8` line elsewhere; `LOSAT_X_NEWTONLANES8=0` runs
// four lanes only. `LOSAT_X_DEFAULT=0` with both lane variables set to a value other than 0 is the
// run those variables gave before the switch joined the set. The report is the same in every run.
#[test]
fn the_default_takes_eight_lane_batches_on_aarch64_only() {
    const LANES: &str = "[X_STATS] NEWTONLANES ";
    const LANES8: &str = "[X_STATS] NEWTONLANES8 ";
    let before = [
        DEFAULT_OFF,
        ("LOSAT_X_NEWTONLANES", "1"),
        ("LOSAT_X_NEWTONLANES8", "1"),
        STATS,
    ];
    for case in CASES {
        let queries = Queries::new(case, "stats");
        let (expected, _) = search(case, &queries, 1, &[]);
        for (vars, eight) in [
            (&[STATS][..], cfg!(target_arch = "aarch64")),
            (&[SWITCH_OFF, STATS][..], false),
            (&before[..], cfg!(target_arch = "aarch64")),
        ] {
            let (report, stderr) = search(case, &queries, 1, vars);
            assert_eq!(report, expected, "{case:?} {vars:?}");
            let batches = counter(&stderr, LANES, "batches")
                .unwrap_or_else(|| panic!("{case:?} {vars:?}: no lanes line:\n{stderr}"));
            assert!(batches > 0, "{case:?} {vars:?}: {stderr}");
            let batches8 = counter(&stderr, LANES8, "batches8");
            if eight {
                assert!(
                    batches8.is_some_and(|n| n > 0),
                    "{case:?} {vars:?}: {stderr}"
                );
            } else {
                assert_eq!(batches8, None, "{case:?} {vars:?}: {stderr}");
            }
        }
    }
}

// The shadow (set to any value, not in the set) compares every stored result that is used with
// the reference port, eight-lane results included on aarch64, prints its summary line, and
// does not change the report.
#[test]
fn the_shadow_compares_the_lane_results_and_the_report_is_unchanged() {
    for case in CASES {
        let queries = Queries::new(case, "shadow");
        let (expected, _) = search(case, &queries, 1, &[]);
        let (report, stderr) = search(case, &queries, 1, &[SHADOW]);
        assert_eq!(report, expected, "{case:?}");
        let compared = counter(&stderr, "[X_SHADOW] NEWTONLANES ", "calls_compared")
            .unwrap_or_else(|| panic!("{case:?}: no summary line:\n{stderr}"));
        assert!(compared > 0, "{case:?}: {stderr}");
    }
}
