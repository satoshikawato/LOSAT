//! The stage lines of `LOSAT_TIMING` for TBLASTX (`src/algorithm/tblastx/x_search_timing.rs`)
//! on whole searches.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The timer reads a clock and prints; the
//! tests check that the report is the same with it and without it, that the lines that were
//! printed before are still printed in the same form, and that the two new lines account for
//! their intervals: the three parts of `x_search` add up to `search_total`, and the stages of
//! `x_subject` add up to its `sum`. The expected output of the fixture pair
//! (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was produced by NCBI
//! BLAST+ 2.17.0.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const THREADS: [u32; 2] = [1, 4];

/// The lines `LOSAT_TIMING` printed for a TBLASTX batch before the stage lines were added, in
/// their order.
const OLD_LINES: [&str; 10] = [
    "read_queries",
    "build_lookup",
    "read_subjects",
    "scan_subject",
    "ungapped_extend",
    "reevaluate",
    "sum_stats_linking",
    "identity_calc",
    "search_total",
    "total",
];

/// The stages of the `x_subject` line, in their order.
const STAGES: [&str; 15] = [
    "enter",
    "translate",
    "prepare",
    "units",
    "units_release",
    "merge",
    "frame_loop",
    "seed_release",
    "link1",
    "reap",
    "reevaluate",
    "link2",
    "convert",
    "emit",
    "end",
];

/// One TBLASTX run of the fixture query against `subject` (`-outfmt 6` on standard output)
/// with the given variables as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays:
/// a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects
/// nothing. Returns the report and the standard error.
fn tblastx(subject: &Path, threads: u32, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        .arg(subject)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{variables:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

fn ncbi_report() -> Vec<u8> {
    std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output")
}

/// The text after `[TIMING] <name>: ` of the only line with that name.
fn timing_line<'a>(stderr: &'a str, name: &str) -> &'a str {
    let prefix = format!("[TIMING] {name}: ");
    let mut lines = stderr
        .lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()));
    let line = lines
        .next()
        .unwrap_or_else(|| panic!("no [TIMING] {name} line:\n{stderr}"));
    assert!(
        lines.next().is_none(),
        "two [TIMING] {name} lines:\n{stderr}"
    );
    line
}

/// `<seconds>s` at the start of `text`, in milliseconds.
fn millis(text: &str) -> i64 {
    let number = text
        .split('s')
        .next()
        .unwrap_or_else(|| panic!("no seconds in {text}"));
    let (whole, fraction) = number
        .split_once('.')
        .unwrap_or_else(|| panic!("no decimal point in {text}"));
    assert_eq!(fraction.len(), 3, "{text}");
    whole.parse::<i64>().expect("seconds") * 1000 + fraction.parse::<i64>().expect("milliseconds")
}

/// The `name value` pairs of a stage line (`setup 0.001s subjects 0.002s ...`), values in
/// milliseconds, up to the first field that is not such a pair.
fn pairs(line: &str) -> Vec<(&str, i64)> {
    let mut fields = line.split_whitespace();
    let mut found = Vec::new();
    while let (Some(name), Some(value)) = (fields.next(), fields.next()) {
        if !value.ends_with('s') || !value.contains('.') {
            break;
        }
        found.push((name, millis(value)));
    }
    found
}

/// Checks the two stage lines of one run against the old lines of the same run.
fn check_stage_lines(stderr: &str, subjects: u64, what: &str) {
    let search_total = millis(timing_line(stderr, "search_total"));
    let search = pairs(timing_line(stderr, "x_search"));
    let names: Vec<&str> = search.iter().map(|(name, _)| *name).collect();
    assert_eq!(names, ["setup", "subjects", "collect"], "{what}:\n{stderr}");
    let parts: i64 = search.iter().map(|(_, value)| value).sum();
    // Three values rounded to a millisecond against their rounded sum.
    assert!(
        (parts - search_total).abs() <= 2,
        "{what}: x_search adds up to {parts} ms, search_total is {search_total} ms:\n{stderr}"
    );

    let subject_line = timing_line(stderr, "x_subject");
    let subject = pairs(subject_line);
    let names: Vec<&str> = subject.iter().map(|(name, _)| *name).collect();
    let expected: Vec<&str> = STAGES.iter().copied().chain(["sum"]).collect();
    assert_eq!(names, expected, "{what}:\n{stderr}");
    let (stages, sum) = subject.split_at(STAGES.len());
    let stage_sum: i64 = stages.iter().map(|(_, value)| value).sum();
    // Fifteen values rounded to a millisecond against their rounded sum.
    assert!(
        (stage_sum - sum[0].1).abs() <= 8,
        "{what}: the stages add up to {stage_sum} ms, sum is {} ms:\n{stderr}",
        sum[0].1
    );
    assert!(
        subject_line.ends_with(&format!(" (subjects={subjects})")),
        "{what}:\n{stderr}"
    );
}

// No NCBI counterpart (a timer). With LOSAT_TIMING the report is NCBI's, the lines that were
// there before are printed once each in their order and form, and the two stage lines follow
// them and account for their intervals. The default search runs on the unit path, so the
// reference frame loop takes no measurable part of the subject.
#[test]
fn the_stage_lines_follow_the_old_lines_and_add_up() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    for threads in THREADS {
        let what = format!("threads {threads}");
        let (report, stderr) = tblastx(&subject, threads, &[TIMING]);
        assert_eq!(report, expected, "{what}");
        let timing: Vec<&str> = stderr
            .lines()
            .filter_map(|line| line.strip_prefix("[TIMING] "))
            .map(|line| line.split(':').next().unwrap_or(line))
            .collect();
        let old: Vec<&str> = timing.iter().copied().take(OLD_LINES.len()).collect();
        assert_eq!(old, OLD_LINES, "{what}:\n{stderr}");
        assert_eq!(
            &timing[OLD_LINES.len()..OLD_LINES.len() + 2],
            ["x_search", "x_subject"],
            "{what}:\n{stderr}"
        );
        for name in OLD_LINES {
            let value = timing_line(&stderr, name);
            assert!(millis(value) >= 0, "{what}: {name}: {value}");
        }
        for name in [
            "scan_subject",
            "ungapped_extend",
            "reevaluate",
            "sum_stats_linking",
        ] {
            assert!(
                timing_line(&stderr, name).contains("s (calls="),
                "{what}: {name}:\n{stderr}"
            );
        }
        check_stage_lines(&stderr, 1, &what);
    }
}

// No NCBI counterpart (a timer). The paths a subject can take all reach the stage lines: the
// reference frame loop (`LOSAT_X_PAIRPAR=0`), the previous paths (`LOSAT_X_DEFAULT=0`), and
// several subject records with one too short to search (counted as a subject, all of it in the
// last stage). The report does not depend on the timer.
#[test]
fn every_path_of_a_subject_is_timed_and_the_report_is_the_same() {
    let single = fixture("short_hsp.subject.fna");
    let first: String = std::fs::read_to_string(&single)
        .expect("FASTA")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let reversed: String = first.chars().rev().collect();
    let several = std::env::temp_dir().join(format!(
        "losat_search_timing_{}_several.fna",
        std::process::id()
    ));
    std::fs::write(
        &several,
        format!(">first\n{first}\n>short\nACGTACGTAC\n>third\n{reversed}\n"),
    )
    .expect("write FASTA");
    for (subject, subjects) in [(&single, 1), (&several, 3)] {
        for path in [
            &[][..],
            &[("LOSAT_X_PAIRPAR", "0")][..],
            &[("LOSAT_X_DEFAULT", "0")][..],
        ] {
            for threads in THREADS {
                let what = format!("{subjects} subject(s) {path:?} threads {threads}");
                let (expected, plain_stderr) = tblastx(subject, threads, path);
                assert!(!expected.is_empty(), "{what}: no hits");
                assert!(
                    !plain_stderr.contains("[TIMING]"),
                    "{what}: a [TIMING] line without LOSAT_TIMING:\n{plain_stderr}"
                );
                let variables = [&[TIMING][..], path].concat();
                let (report, stderr) = tblastx(subject, threads, &variables);
                assert_eq!(report, expected, "{what}");
                check_stage_lines(&stderr, subjects, &what);
            }
        }
    }
    let _ = std::fs::remove_file(&several);
}
