//! The diagonal-hash fold of the strips paths split by bucket over the search pool
//! (`LOSAT_X_PAIRPARFOLD`, `src/algorithm/blastn/blast_engine/x_pair_seed_part.rs`) on whole
//! searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses between
//! two LOSAT code paths with the same result: the word hits of the strips paths
//! (PAIRPAR, which has no variable, for megablast and blastn; `LOSAT_X_PAIRPARDISC` for
//! dc-megablast) folded into the diagonal hash by the owner thread in scan order, or by every
//! thread of the pool, each into the buckets of its partition, with the hits and the table's
//! cell numbers put back in the reference order. So the property under test is that the report
//! does not depend on it: with the switch at 2, 4 and 8 threads (also with `LOSAT_X_DIAGLV`),
//! and (dc-megablast) with the strips shadow `LOSAT_X_PAIRPARDISCSHADOW`, which runs the
//! partitioned stage on copies and compares the hit list and the whole diagonal hash with the
//! reference fold for every chunk, the searches give the bytes of the search without it at one
//! thread. That those bytes equal NCBI's is the business of the BLASTN fixtures.
//!
//! The query (LC738874, 287 kb) makes a query block longer than 8000 letters, so the searches
//! use the diagonal hash, and the one subject (LC738868) is searched by the owner of the pool,
//! which the strips paths need.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), and
/// LOSAT_X_MBLUTONCE on 2026-10-18 (perf/bn-memory-rules), each `0`, the default before, unless the
/// search sets it: every search here runs as it was written, with the switches off that it does not
/// name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_MBLUTONCE", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_PAIRPARFOLD", "1");
const PAIRPARDISC: (&str, &str) = ("LOSAT_X_PAIRPARDISC", "1");
const DIAGLV: (&str, &str) = ("LOSAT_X_DIAGLV", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of this test process for one test (`losat-pairparfold-<test>-<pid>`), removed
/// when the test ends; a failing test keeps it for its files.
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-pairparfold-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }
}

impl std::ops::Deref for TestDir {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// One search of LC738874 (query) against LC738868 (subject), `-outfmt 6` on standard output,
/// with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays (or is set to a file of `dir`): a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name
/// as soon as a search has more than one thread, and it selects nothing. Returns the report and
/// the standard error.
fn blastn(dir: &Path, task: &str, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-task", task, "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The fields of a strips shadow's line (`prefix`, the only line of `stderr`): the chunks
/// compared, then `fold=`, `merges=` and `partitions=`.
fn shadow_fields(stderr: &str, prefix: &str, at: &str) -> (u64, String, u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(prefix))
        .unwrap_or_else(|| panic!("{at}: no shadow line:\n{stderr}"));
    assert!(
        stderr.lines().all(|l| l == line),
        "{at}: unexpected standard error line:\n{stderr}"
    );
    let field = |name: &str| -> String {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("{at}: no {name} in: {line}"))
            .to_string()
    };
    let number = |name: &str| -> u64 {
        field(name)
            .parse()
            .unwrap_or_else(|_| panic!("{at}: {name} in: {line}"))
    };
    let (merges, partitions) = (number("merges="), number("partitions="));
    let chunks = line
        .split("seed chunks equal=")
        .nth(1)
        .and_then(|rest| rest.split_whitespace().next())
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("{at}: no seed chunks in: {line}"));
    (chunks, field("fold="), merges, partitions)
}

/// For one task: the search without the switch at one thread, then at 2, 4 and 8 threads
/// without the switch, with it, with it and `LOSAT_X_DIAGLV`, and, if there is one, with it and
/// the strips shadow (`shadow`, whose line starts with `prefix`): the same bytes every time,
/// nothing on standard error without the shadow, and the shadow's line says that it compared
/// chunks folded in partitions, one per thread. `extra` goes to every run but the first.
fn check(dir: &Path, task: &str, extra: &[(&str, &str)], shadow: Option<((&str, &str), &str)>) {
    let (expected, expected_stderr) = blastn(dir, task, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [2u32, 4, 8] {
        let at = format!("{task} {extra:?} threads {threads}");
        let (report, stderr) = blastn(dir, task, threads, extra);
        assert!(
            report == expected,
            "{at}: the report differs without the switch"
        );
        assert_eq!(stderr, "", "{at}");
        for more in [&[SWITCH][..], &[SWITCH, DIAGLV][..]] {
            let mut on = extra.to_vec();
            on.extend_from_slice(more);
            let (report, stderr) = blastn(dir, task, threads, &on);
            assert!(report == expected, "{at}: the report differs with {more:?}");
            assert_eq!(stderr, "", "{at} {more:?}");
        }
        let Some((shadow, prefix)) = shadow else {
            continue;
        };
        let mut on = extra.to_vec();
        on.extend_from_slice(&[SWITCH, shadow]);
        let (report, stderr) = blastn(dir, task, threads, &on);
        assert!(
            report == expected,
            "{at}: the report differs with {shadow:?}"
        );
        let (chunks, fold, merges, partitions) = shadow_fields(&stderr, prefix, &at);
        assert!(chunks > 0, "{at}: no chunk compared:\n{stderr}");
        assert_eq!(fold, "part", "{at}");
        assert_eq!(
            merges, chunks,
            "{at}: chunks not folded in partitions:\n{stderr}"
        );
        assert_eq!(partitions, u64::from(threads), "{at}");
    }
}

#[test]
fn megablast_gives_the_same_report_with_the_partitioned_fold() {
    let dir = TestDir::new("megablast");
    check(&dir, "megablast", &[], None);
}

#[test]
fn blastn_gives_the_same_report_with_the_partitioned_fold() {
    let dir = TestDir::new("blastn");
    check(&dir, "blastn", &[], None);
}

#[test]
fn dc_megablast_gives_the_same_report_with_the_partitioned_fold() {
    let dir = TestDir::new("dc");
    check(
        &dir,
        "dc-megablast",
        &[PAIRPARDISC],
        Some((
            ("LOSAT_X_PAIRPARDISCSHADOW", "1"),
            "[X_PAIRPARDISCSHADOW] blastn:",
        )),
    );
}

// Without `LOSAT_X_PAIRPARDISC` dc-megablast has no strips path, so the switch changes nothing.
#[test]
fn dc_megablast_without_its_strips_path_is_not_taken() {
    let dir = TestDir::new("dcoff");
    let (expected, _) = blastn(&dir, "dc-megablast", 1, &[]);
    let (report, stderr) = blastn(&dir, "dc-megablast", 4, &[SWITCH]);
    assert!(report == expected, "the report differs");
    assert_eq!(stderr, "");
}

// LOSAT_TIMING: the strips paths' `[TIMING] x_pairpar_seed` / `x_pairpardisc_seed` lines end
// with the partitioned fold's fields (`partitions= merges= merge= helpers_left=`, and `runs=
// steals=` of LOSAT_X_PAIRPARHOME, on by default) only with the switch; without it they end with
// `owner_wait=`, as before the switch existed. (Names and counts only; no duration is compared.)
#[test]
fn timing_lines_have_the_fold_fields_only_with_the_switch() {
    let dir = TestDir::new("timing");
    let timing = ("LOSAT_TIMING", "1");
    for (task, extra, prefix) in [
        ("megablast", &[][..], "[TIMING] x_pairpar_seed:"),
        (
            "dc-megablast",
            &[PAIRPARDISC][..],
            "[TIMING] x_pairpardisc_seed:",
        ),
    ] {
        for on in [false, true] {
            let mut switches = extra.to_vec();
            switches.push(timing);
            if on {
                switches.push(SWITCH);
            }
            let at = format!("{task} switch {on}");
            let (_, stderr) = blastn(&dir, task, 4, &switches);
            let line = stderr
                .lines()
                .find(|line| line.starts_with(prefix))
                .unwrap_or_else(|| panic!("{at}: no {prefix} line:\n{stderr}"));
            let names: Vec<&str> = line
                .split_whitespace()
                .skip(2)
                .map(|word| word.split('=').next().unwrap_or(word))
                .collect();
            let mut expected = vec![
                "chunks",
                "strips",
                "seeds",
                "skipped",
                "fold",
                "owner_scan",
                "owner_wait",
            ];
            if on {
                expected.extend([
                    "partitions",
                    "merges",
                    "merge",
                    "helpers_left",
                    "runs",
                    "steals",
                ]);
            }
            assert_eq!(names, expected, "{at}: {line}");
            if on {
                let number = |name: &str| -> u64 {
                    line.split_whitespace()
                        .find_map(|word| word.strip_prefix(name))
                        .and_then(|value| value.parse().ok())
                        .unwrap_or_else(|| panic!("{at}: {name} in: {line}"))
                };
                assert_eq!(number("partitions="), 4, "{at}: {line}");
                assert!(number("chunks=") > 0, "{at}: {line}");
                assert_eq!(number("merges="), number("chunks="), "{at}: {line}");
                assert!(number("runs=") >= number("steals="), "{at}: {line}");
            }
        }
    }
}
