//! The stage lines of `LOSAT_TIMING` for TBLASTX (`src/algorithm/tblastx/x_search_timing.rs`,
//! `x_link_timing.rs`) on whole searches.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The timer reads a clock and prints; the
//! tests check that the report is the same with it and without it, that the lines that were
//! printed before are still printed in the same form, and that the new lines account for their
//! intervals: the three parts of `x_search` add up to `search_total`, the stages of `x_subject`
//! add up to its `sum`, and the stages of `x_link1` and `x_link2` add up to their `sum`, which
//! lies inside `link1` and `link2` of `x_subject`. The expected output of the fixture pair
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
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
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

/// The stages of the `x_link1` and `x_link2` lines, in their order.
const LINK_STAGES: [&str; 7] = [
    "reset",
    "sort",
    "groups",
    "kernels",
    "assemble",
    "final_sorts",
    "replay",
];

/// One TBLASTX run of the fixture query against `subject` (`-outfmt 6` on standard output)
/// with the given variables as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays:
/// a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects
/// nothing. Returns the report and the standard error.
fn tblastx(subject: &Path, threads: u32, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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

/// Checks the two linking lines of one run against its `x_subject` line. `timed` says whether
/// the run's linking calls take the path the timer is in (the in-place path): then the first
/// call is counted for at least one subject and at most `subjects`, and the second for no more
/// subjects than the first; otherwise both lines are all zero.
fn check_link_lines(stderr: &str, subjects: u64, timed: bool, what: &str) {
    let subject = pairs(timing_line(stderr, "x_subject"));
    let mut calls = Vec::new();
    for (name, stage) in [("x_link1", "link1"), ("x_link2", "link2")] {
        let line = timing_line(stderr, name);
        let fields = pairs(line);
        let names: Vec<&str> = fields.iter().map(|(field, _)| *field).collect();
        let expected: Vec<&str> = LINK_STAGES.iter().copied().chain(["sum"]).collect();
        assert_eq!(names, expected, "{what}: {name}:\n{stderr}");
        let (stages, sum) = fields.split_at(LINK_STAGES.len());
        let stage_sum: i64 = stages.iter().map(|(_, value)| value).sum();
        // Seven values rounded to a millisecond against their rounded sum.
        assert!(
            (stage_sum - sum[0].1).abs() <= 4,
            "{what}: {name}: the stages add up to {stage_sum} ms, sum is {} ms:\n{stderr}",
            sum[0].1
        );
        // The timed part of a linking call lies inside the call the subject's clock times.
        let whole = subject
            .iter()
            .find(|(field, _)| *field == stage)
            .unwrap_or_else(|| panic!("{what}: no {stage} in x_subject:\n{stderr}"))
            .1;
        assert!(
            sum[0].1 <= whole + 1,
            "{what}: {name}: sum {} ms, {stage} {whole} ms:\n{stderr}",
            sum[0].1
        );
        let tail = line.split_once('(').map(|(_, tail)| tail).unwrap_or("");
        let fields: Vec<&str> = tail.split_whitespace().collect();
        assert!(
            fields.len() == 8
                && fields[0].starts_with("calls=")
                && fields[1].starts_with("hsps=")
                && fields[2].starts_with("groups=")
                && fields[3].starts_with("largest_group=")
                && fields[4] == "kernel_sum"
                && fields[6] == "kernel_max"
                && fields[7].ends_with("s)"),
            "{what}: {name}:\n{stderr}"
        );
        let (kernel_sum, kernel_max) = (millis(fields[5]), millis(fields[7]));
        assert!(kernel_max <= kernel_sum + 1, "{what}: {name}:\n{stderr}");
        let (n, hsps, groups, largest) = (
            count(line, "calls"),
            count(line, "hsps"),
            count(line, "groups"),
            count(line, "largest_group"),
        );
        // Every timed call has a list of at least one HSP in at least one group.
        assert!(
            n <= groups && groups <= hsps && n <= largest && largest <= hsps,
            "{what}: {name}:\n{stderr}"
        );
        if n == 0 {
            assert_eq!((sum[0].1, stage_sum, hsps), (0, 0, 0), "{what}: {name}");
        }
        calls.push(n);
    }
    if timed {
        assert!(
            1 <= calls[0] && calls[0] <= subjects && calls[1] <= calls[0],
            "{what}: calls {calls:?}:\n{stderr}"
        );
    } else {
        assert_eq!(calls, [0, 0], "{what}:\n{stderr}");
    }
}

// No NCBI counterpart (a timer). With LOSAT_TIMING the report is NCBI's, the lines that were
// there before are printed once each in their order and form, and the stage lines follow them
// and account for their intervals. The default search runs on the unit path, so the
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
            &timing[OLD_LINES.len()..OLD_LINES.len() + 4],
            ["x_search", "x_subject", "x_link1", "x_link2"],
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
        check_link_lines(&stderr, 1, true, &what);
    }
}

// No NCBI counterpart (a timer). The paths a subject can take all reach the stage lines: the
// default (the reference frame loop at one thread, the unit path at four), `LOSAT_X_PAIRPARGATE=0`
// (the unit path at one thread too), the debugging output `LOSAT_DEBUG_CUTOFFS` (the one path whose
// linking calls take the reference body, which is not timed), and several subject records with one
// too short to search (counted as a subject, all of it in the last stage). The linking lines count the calls of the in-place linking
// path, whichever seed stage ran before it, and are all zero with `LOSAT_DEBUG_CUTOFFS`. The report
// does not depend on the timer. At one thread the cost rule sends the subject to the reference
// frame loop (`[X_STATS] PAIRPAR units=0`; not with `LOSAT_X_PAIRPARGATE=0`, which turns the cost
// rule off). The `LOSAT_DEBUG_CUTOFFS` arm is the only search of this file whose linking calls
// take the reference body, so the report of each other arm is also compared with the default
// path's report of the same subject and thread count (the debug lines are on stderr).
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
    let mut default_reports = std::collections::HashMap::new();
    for (subject, subjects) in [(&single, 1), (&several, 3)] {
        for path in [
            &[][..],
            &[("LOSAT_X_PAIRPARGATE", "0")][..],
            &[("LOSAT_DEBUG_CUTOFFS", "1")][..],
        ] {
            for threads in THREADS {
                let what = format!("{subjects} subject(s) {path:?} threads {threads}");
                let (expected, plain_stderr) =
                    tblastx(subject, threads, &[path, &[STATS][..]].concat());
                assert!(!expected.is_empty(), "{what}: no hits");
                if threads == 1 && !path.contains(&("LOSAT_X_PAIRPARGATE", "0")) {
                    let line = plain_stderr
                        .lines()
                        .find(|line| line.starts_with("[X_STATS] PAIRPAR "))
                        .unwrap_or_else(|| panic!("{what}: no PAIRPAR line:\n{plain_stderr}"));
                    assert!(
                        line.split_whitespace().any(|field| field == "units=0"),
                        "{what}: the frame loop:\n{plain_stderr}"
                    );
                }
                if path.is_empty() {
                    default_reports.insert((subjects, threads), expected.clone());
                } else {
                    assert_eq!(
                        expected,
                        default_reports[&(subjects, threads)],
                        "{what}: the report of the default path"
                    );
                }
                assert!(
                    !plain_stderr.contains("[TIMING]"),
                    "{what}: a [TIMING] line without LOSAT_TIMING:\n{plain_stderr}"
                );
                let variables = [&[TIMING][..], path].concat();
                let (report, stderr) = tblastx(subject, threads, &variables);
                assert_eq!(report, expected, "{what}");
                check_stage_lines(&stderr, subjects, &what);
                let timed = !path.contains(&("LOSAT_DEBUG_CUTOFFS", "1"));
                check_link_lines(&stderr, subjects, timed, &what);
            }
        }
    }
    let _ = std::fs::remove_file(&several);
}
