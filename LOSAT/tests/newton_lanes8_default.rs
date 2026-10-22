//! The eight Newton lanes (`src/core/composition_adjustment/x_newton_lanes.rs`) on whole blastp and
//! tblastn searches.
//!
//! The lanes have no variable since `LOSAT_X_NEWTONLANES` and `LOSAT_X_NEWTONLANES8` graduated
//! (2026-10-22). On aarch64 batches of five or more problems share eight lanes; on other targets
//! the lanes stay four wide. The `[X_STATS] NEWTONLANES8` line shows that a default run takes
//! eight-lane batches on aarch64 and none elsewhere, so this test means most on the arm64 runners
//! of platform-checks. That the bytes equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

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
/// variables. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
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

// With `LOSAT_X_STATS=1` the lanes report their batches. The default run takes eight-lane
// batches on aarch64 and prints no `NEWTONLANES8` line elsewhere. The report is the same with
// and without the report of the counters.
#[test]
fn the_default_takes_eight_lane_batches_on_aarch64_only() {
    const LANES: &str = "[X_STATS] NEWTONLANES ";
    const LANES8: &str = "[X_STATS] NEWTONLANES8 ";
    for case in CASES {
        let queries = Queries::new(case, "stats");
        let (expected, expected_stderr) = search(case, &queries, 1, &[]);
        assert!(!expected.is_empty(), "{case:?}: no hits");
        assert_eq!(expected_stderr, "", "{case:?}");
        let (report, stderr) = search(case, &queries, 1, &[STATS]);
        assert_eq!(report, expected, "{case:?}");
        let batches = counter(&stderr, LANES, "batches")
            .unwrap_or_else(|| panic!("{case:?}: no lanes line:\n{stderr}"));
        assert!(batches > 0, "{case:?}: {stderr}");
        let batches8 = counter(&stderr, LANES8, "batches8");
        if cfg!(target_arch = "aarch64") {
            assert!(batches8.is_some_and(|n| n > 0), "{case:?}: {stderr}");
        } else {
            assert_eq!(batches8, None, "{case:?}: {stderr}");
        }
    }
}
