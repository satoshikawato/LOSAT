//! The level bounds of the diagonal hash (`LOSAT_X_DIAGLV`,
//! `src/algorithm/blastn/blast_engine/x_diag_lv.rs`) on whole searches.
//!
//! The switch is in the strict set since 2026-10-17; these searches run with it `0` unless they
//! set it (`PREVIOUS_DEFAULT`). It chooses between
//! two LOSAT code paths with the same result: a seed's retrieve of its diagonal walks the chain
//! of the diagonal's bucket, or returns "no cell" without the walk when every cell of the
//! bucket lies below the seed's level bound. So the property under test is that the report does
//! not depend on it: with the switch the searches give the bytes of the search without it at one
//! thread, and the reference fold ran with the bounds (`LOSAT_X_REACH`). That the skipped walks
//! change no decision is checked by the unit tests of `x_diag_lv.rs`. That those bytes equal
//! NCBI's is the business of the BLASTN fixtures.
//!
//! The query (LC738874, 287 kb) makes a query block longer than 8000 letters, so the searches
//! use the diagonal hash. The folds reached: at one thread the reference fold
//! (`process_offset_pairs`: its first retrieve for megablast, blastn and dc-megablast, whose
//! tables have a word extension or a template; its second for blastn-short, whose table has
//! none); at four threads with one subject the strips folds (PAIRPAR, which has no
//! variable, for megablast and blastn; `LOSAT_X_PAIRPARDISC` for dc-megablast when it is set), and
//! with many subjects the reference fold on every worker, whose hash offset moves on from subject
//! to subject. The searches that set `LOSAT_X_BNDCONE` take the dc-megablast strips stage at one
//! thread too.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), on
/// 2026-10-18 `LOSAT_X_MBLUTONCE` (perf/bn-memory-rules), on 2026-10-22 `LOSAT_X_DIAGLV2` and
/// `LOSAT_X_FOLDLEAN` (perf/dc-fold-lean), and on 2026-10-24 `LOSAT_X_BNDCONE`
/// (perf/dc-one-thread), each `0`, the default before, unless the search sets it: every search
/// here runs as it was written, with the switches off that it does not name (the default runs in
/// `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BNDCONE", "0"),
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_DIAGLV2", "0"),
    ("LOSAT_X_FOLDLEAN", "0"),
    ("LOSAT_X_MBLUTONCE", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_DIAGLV", "1");
const PAIRPARDISC: (&str, &str) = ("LOSAT_X_PAIRPARDISC", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of this test process for one test (`losat-diaglv-<test>-<pid>`), removed when
/// the test ends; a failing test keeps it for its files.
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-diaglv-{test}-{}", std::process::id()));
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

/// The letters of a FASTA file's records, one string per record.
fn records(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    let mut out: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            out.push(String::new());
        } else if let Some(last) = out.last_mut() {
            last.push_str(line.trim());
        }
    }
    out
}

/// One `blastn` search, `-outfmt 6` on standard output, with the given variables as its only
/// `LOSAT_*` variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays (or is set to a
/// file of `dir`): a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn blastn(
    dir: &Path,
    query: &Path,
    subject: &Path,
    task: &str,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    blastn_with(dir, query, subject, task, threads, switches, &[])
}

/// `blastn` with `args` after the others (`-lcase_masking`).
fn blastn_with(
    dir: &Path,
    query: &Path,
    subject: &Path,
    task: &str,
    threads: u32,
    switches: &[(&str, &str)],
    args: &[&str],
) -> (Vec<u8>, String) {
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
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-task", task, "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .args(args)
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {switches:?} {args:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For one task, query and subject: the search without the switch at one thread, then with
/// the switch at 1 and 4 threads (the same bytes, nothing on standard error; at one thread the
/// reference fold ran with the bounds), and the search without the switch at 4 threads; `extra`
/// adds variables to every run but the first (the strips fold of dc-megablast).
fn check(dir: &Path, query: &Path, subject: &Path, task: &str, extra: &[(&str, &str)]) -> Vec<u8> {
    let (expected, expected_stderr) = blastn(dir, query, subject, task, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [1, 4] {
        let at = format!("{task} {extra:?} threads {threads}");
        let reach = reach_support::Reach::new("diaglv");
        let mut on = extra.to_vec();
        on.push(SWITCH);
        on.push(reach.var());
        let (report, stderr) = blastn(dir, query, subject, task, threads, &on);
        assert!(
            report == expected,
            "{at}: the report differs with {SWITCH:?}"
        );
        assert_eq!(stderr, "", "{at}");
        if threads == 1 {
            assert!(
                reach.sides("LOSAT_X_DIAGLV").0,
                "{at}: the reference fold did not run with the bounds"
            );
        }
    }
    let (report, stderr) = blastn(dir, query, subject, task, 4, extra);
    assert!(
        report == expected,
        "{task} {extra:?}: the report differs at 4 threads"
    );
    assert_eq!(stderr, "", "{task}");
    expected
}

#[test]
fn one_subject_gives_the_same_report_with_the_level_bounds() {
    let dir = TestDir::new("one");
    let query = fasta("LC738874.fasta");
    let subject = fasta("LC738868.fasta");
    for task in ["megablast", "dc-megablast", "blastn"] {
        check(&dir, &query, &subject, task, &[]);
    }
    // the strips fold of dc-megablast
    check(&dir, &query, &subject, "dc-megablast", &[PAIRPARDISC]);
}

// The memo of the latest insert (`LOSAT_X_DIAGLV2`, with the level bounds) and the lean
// instances of the dc-megablast strips folds (`LOSAT_X_FOLDLEAN`): dc-megablast at one thread
// (the reference fold) and at four (the strips fold of `LOSAT_X_PAIRPARDISC`) gives the bytes of
// the search with neither, each switch alone and both, at four threads also with the fold split
// over partitions (`LOSAT_X_PAIRPARFOLD`); the memo exists in the strips stage at four threads
// (`LOSAT_X_REACH`) and not at one (the reference fold does not read it) unless
// `LOSAT_X_BNDCONE` runs the strips stage there too. megablast keeps no memo (a one-hit table).
#[test]
fn the_memo_and_the_lean_folds_give_the_same_report() {
    let dir = TestDir::new("memo");
    let query = fasta("LC738874.fasta");
    let subject = fasta("LC738868.fasta");
    let task = "dc-megablast";
    let (expected, expected_stderr) = blastn(&dir, &query, &subject, task, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    const MEMO: (&str, &str) = ("LOSAT_X_DIAGLV2", "1");
    const LEAN: (&str, &str) = ("LOSAT_X_FOLDLEAN", "1");
    const FOLD: (&str, &str) = ("LOSAT_X_PAIRPARFOLD", "1");
    for threads in [1, 4] {
        let folds: &[bool] = if threads > 1 {
            &[false, true]
        } else {
            &[false]
        };
        for &fold in folds {
            for arm in [vec![SWITCH, MEMO], vec![LEAN], vec![SWITCH, MEMO, LEAN]] {
                let mut on = arm.clone();
                on.push(PAIRPARDISC);
                if fold {
                    on.push(FOLD);
                }
                let at = format!("{task} {on:?} threads {threads}");
                let (report, stderr) = blastn(&dir, &query, &subject, task, threads, &on);
                assert!(report == expected, "{at}: the report differs");
                assert_eq!(stderr, "", "{at}");
            }
        }
        let reach = reach_support::Reach::new("diaglv-memo");
        let on = [SWITCH, MEMO, LEAN, PAIRPARDISC, reach.var()];
        let at = format!("{task} {on:?} threads {threads}");
        let (report, stderr) = blastn(&dir, &query, &subject, task, threads, &on);
        assert!(report == expected, "{at}: the report differs");
        assert_eq!(stderr, "", "{at}");
        assert_eq!(
            reach.sides("LOSAT_X_DIAGLV2").0,
            threads > 1,
            "{at}: the memo of the strips stage"
        );
    }
    let reach = reach_support::Reach::new("diaglv-memo-one");
    let on = [
        SWITCH,
        MEMO,
        LEAN,
        PAIRPARDISC,
        ("LOSAT_X_BNDCONE", "1"),
        reach.var(),
    ];
    let at = format!("{task} {on:?} threads 1");
    let (report, stderr) = blastn(&dir, &query, &subject, task, 1, &on);
    assert!(report == expected, "{at}: the report differs");
    assert_eq!(stderr, "", "{at}");
    assert!(
        reach.sides("LOSAT_X_DIAGLV2").0,
        "{at}: the strips stage at one thread kept no memo"
    );
    let reach = reach_support::Reach::new("diaglv-memo-megablast");
    let on = [SWITCH, MEMO, reach.var()];
    let at = format!("megablast {on:?}");
    let (_, stderr) = blastn(&dir, &query, &subject, "megablast", 1, &on);
    assert_eq!(stderr, "", "{at}");
    assert!(
        !reach.sides("LOSAT_X_DIAGLV2").0,
        "{at}: a one-hit table kept a memo"
    );
}

// blastn-short: a lookup table of the word size (7) without the two-stage word extension, so
// the reference fold reads the diagonal hash through its second retrieve. The query is the
// first 9,000 letters of LC738874 (both strands: a block of 18,001 letters).
#[test]
fn blastn_short_gives_the_same_report_with_the_level_bounds() {
    let dir = TestDir::new("short");
    let letters = records(&fasta("LC738874.fasta")).remove(0);
    let query = dir.join("query.fna");
    std::fs::write(&query, format!(">q9000\n{}\n", &letters[..9_000])).expect("write the query");
    check(&dir, &query, &fasta("LC738868.fasta"), "blastn-short", &[]);
}

// Many subjects: 61 pieces of LC738868 (5,000 letters every 5,000, and an 8-letter piece
// shorter than the word, whose search only moves the hash offset on), searched one after the
// other at one thread and by the workers at four, so that every worker's diagonal hash keeps
// its cells and bounds from subject to subject with a growing offset. With `LOSAT_X_BNDCONE` at
// one thread every subject takes the dc-megablast strips stage in turn on that one table, and
// the memo answers; with `-lcase_masking` and every
// third piece in lowercase those subjects take the reference fold between the stages, whose
// memo is made again at each stage's start (`[TIMING] x_pairpardisc_seed:` counts the chunks the
// stage took).
#[test]
fn many_subjects_give_the_same_report_with_the_level_bounds() {
    let dir = TestDir::new("many");
    let letters = records(&fasta("LC738868.fasta")).remove(0);
    let subjects = dir.join("subjects.fna");
    let mut text = String::new();
    for index in 0..60 {
        let start = index * 5_000;
        let end = (start + 5_000).min(letters.len());
        text.push_str(&format!(">s{index}\n{}\n", &letters[start..end]));
        if index == 20 {
            text.push_str(&format!(">short\n{}\n", &letters[100..108]));
        }
    }
    std::fs::write(&subjects, text).expect("write the subjects");
    let query = fasta("LC738874.fasta");
    let mut dc = Vec::new();
    for task in ["megablast", "dc-megablast"] {
        dc = check(&dir, &query, &subjects, task, &[]);
    }
    const ONE: [(&str, &str); 4] = [
        ("LOSAT_X_DIAGLV2", "1"),
        ("LOSAT_X_FOLDLEAN", "1"),
        PAIRPARDISC,
        ("LOSAT_X_BNDCONE", "1"),
    ];
    let task = "dc-megablast";
    let reach = reach_support::Reach::new("diaglv-many");
    let mut on = ONE.to_vec();
    on.push(SWITCH);
    on.push(reach.var());
    let at = format!("{task} {on:?} threads 1");
    let (report, stderr) = blastn(&dir, &query, &subjects, task, 1, &on);
    assert!(report == dc, "{at}: the report differs");
    assert_eq!(stderr, "", "{at}");
    assert!(
        reach.sides("LOSAT_X_DIAGLV2").0,
        "{at}: the strips stage kept no memo"
    );

    let masked = dir.join("masked.fna");
    let mut text = String::new();
    for (index, record) in std::fs::read_to_string(&subjects)
        .expect("read the subjects")
        .split('>')
        .filter(|record| !record.is_empty())
        .enumerate()
    {
        let (name, letters) = record.split_once('\n').expect("a record");
        let letters = letters.trim_end();
        let letters = if index % 3 == 0 {
            letters.to_ascii_lowercase()
        } else {
            letters.to_string()
        };
        text.push_str(&format!(">{name}\n{letters}\n"));
    }
    std::fs::write(&masked, text).expect("write the masked subjects");
    let lcase = ["-lcase_masking"];
    let (expected, expected_stderr) = blastn_with(&dir, &query, &masked, task, 1, &[], &lcase);
    assert!(!expected.is_empty(), "{task} {lcase:?}: no hits");
    assert_eq!(expected_stderr, "", "{task} {lcase:?}");
    let reach = reach_support::Reach::new("diaglv-masked");
    let mut on = ONE.to_vec();
    on.push(SWITCH);
    on.push(("LOSAT_TIMING", "1"));
    on.push(reach.var());
    let at = format!("{task} {on:?} {lcase:?} threads 1");
    let (report, stderr) = blastn_with(&dir, &query, &masked, task, 1, &on, &lcase);
    assert!(report == expected, "{at}: the report differs");
    assert!(
        reach.sides("LOSAT_X_DIAGLV2").0,
        "{at}: the strips stage kept no memo"
    );
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] x_pairpardisc_seed: chunks="))
        .unwrap_or_else(|| panic!("{at}: no stage line:\n{stderr}"));
    let chunks: u64 = line["[TIMING] x_pairpardisc_seed: chunks=".len()..]
        .split_whitespace()
        .next()
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("{at}: {line}"));
    // 61 records, every third in lowercase (21, among them the 8-letter one): the other 40 take
    // the stage, one chunk each
    assert_eq!(chunks, 40, "{at}: {line}");
}
