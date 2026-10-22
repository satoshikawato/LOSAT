//! The level bounds of the diagonal hash (`LOSAT_X_DIAGLV`,
//! `src/algorithm/blastn/blast_engine/x_diag_lv.rs`) on whole searches.
//!
//! The switch is in the strict set since 2026-10-17; these searches run with it `0` unless they
//! set it (`PREVIOUS_DEFAULT`). It chooses between
//! two LOSAT code paths with the same result: a seed's retrieve of its diagonal walks the chain
//! of the diagonal's bucket, or returns "no cell" without the walk when every cell of the
//! bucket lies below the seed's level bound. So the property under test is that the report does
//! not depend on it: with the switch, and with its shadow (`LOSAT_X_DIAGLVSHADOW`, which walks
//! every skipped bucket, checks each cell against the bound and checks that the cell the walk
//! returns gives the seed's decisions of no cell), the searches give the bytes of the search
//! without it at one thread, and the shadow skips walks. That those bytes equal
//! NCBI's is the business of the BLASTN fixtures.
//!
//! The query (LC738874, 287 kb) makes a query block longer than 8000 letters, so the searches
//! use the diagonal hash. The folds reached: at one thread the reference fold
//! (`process_offset_pairs`: its first retrieve for megablast, blastn and dc-megablast, whose
//! tables have a word extension or a template; its second for blastn-short, whose table has
//! none); at four threads with one subject the strips folds (`LOSAT_X_PAIRPAR`, on by default,
//! for megablast and blastn; `LOSAT_X_PAIRPARDISC` for dc-megablast when it is set), and with
//! many subjects the reference fold on every worker, whose hash offset moves on from subject to
//! subject.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), on
/// 2026-10-18 the one whose shadow line `LOSAT_X_PAIRPARSHADOW` adds (`LOSAT_X_MBLUTONCE`,
/// perf/bn-memory-rules), and on 2026-10-22 `LOSAT_X_DIAGLV2` and `LOSAT_X_FOLDLEAN`
/// (perf/dc-fold-lean), each `0`, the default before, unless the search sets it: every search
/// here runs as it was written, with the switches off that it does not name (the default runs in
/// `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
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
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_DIAGLV", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_DIAGLVSHADOW", "1");
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

/// The counters of the shadow's line (`[X_SHADOW] DIAGLV retrieves= skipped= verified=
/// decisions= cells= memo= memofound= identical`), which must be the only line of `stderr`
/// unless `others` lines start with one of the other prefixes: the retrieves and the skipped
/// ones, after
/// checking that every skipped walk passed both checks (its bucket's levels and its
/// decisions).
fn shadow_counts(stderr: &str, others: &[&str], at: &str) -> (u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] DIAGLV "))
        .unwrap_or_else(|| panic!("{at}: no shadow line:\n{stderr}"));
    assert!(line.ends_with(" identical"), "{at}: {line}");
    for other in stderr.lines().filter(|l| *l != line) {
        assert!(
            others.iter().any(|prefix| other.starts_with(prefix)),
            "{at}: unexpected standard error line: {other}"
        );
    }
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("{at}: no {name} in: {line}"))
    };
    let (retrieves, skipped) = (field("retrieves="), field("skipped="));
    assert_eq!(field("verified="), skipped, "{at}: {line}");
    assert_eq!(field("decisions="), skipped, "{at}: {line}");
    assert!(field("cells=") <= skipped, "{at}: {line}");
    (retrieves, skipped)
}

/// For one task, query and subject: the search without the switch at one thread, then with
/// the switch at 1 and 4 threads (the same bytes, nothing on standard error), the shadow at 1
/// and 4 threads (the same bytes, walks skipped and every skipped walk checked), and the
/// search without the switch at 4 threads; `extra` adds variables to every run but the first
/// (the strips fold of dc-megablast).
fn check(dir: &Path, query: &Path, subject: &Path, task: &str, extra: &[(&str, &str)]) -> Vec<u8> {
    let (expected, expected_stderr) = blastn(dir, query, subject, task, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [1, 4] {
        let at = format!("{task} {extra:?} threads {threads}");
        let mut on = extra.to_vec();
        on.push(SWITCH);
        let (report, stderr) = blastn(dir, query, subject, task, threads, &on);
        assert!(
            report == expected,
            "{at}: the report differs with {SWITCH:?}"
        );
        assert_eq!(stderr, "", "{at}");
        let mut shadow = extra.to_vec();
        shadow.push(SHADOW);
        let (report, stderr) = blastn(dir, query, subject, task, threads, &shadow);
        assert!(
            report == expected,
            "{at}: the report differs with {SHADOW:?}"
        );
        let (retrieves, skipped) = shadow_counts(&stderr, &[], &at);
        assert!(skipped > 0 && skipped <= retrieves, "{at}: {stderr}");
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

// With the strips folds at 4 threads, their own shadows compare the hit list and the whole
// diagonal hash of every chunk, the level bounds included, with the reference fold, both with
// the switch.
#[test]
fn the_strips_shadows_compare_the_tables_with_the_level_bounds() {
    let dir = TestDir::new("strips");
    let query = fasta("LC738874.fasta");
    let subject = fasta("LC738868.fasta");
    // (dc-megablast with the memo, LOSAT_X_DIAGLV2, and the lean folds, LOSAT_X_FOLDLEAN: the
    // strips shadow also compares the memo of every chunk's table)
    for (task, strips_shadow, prefix, more) in [
        (
            "megablast",
            ("LOSAT_X_PAIRPARSHADOW", "1"),
            "[X_PAIRPARSHADOW] blastn:",
            &[][..],
        ),
        (
            "dc-megablast",
            ("LOSAT_X_PAIRPARDISCSHADOW", "1"),
            "[X_PAIRPARDISCSHADOW] blastn:",
            &[("LOSAT_X_DIAGLV2", "1"), ("LOSAT_X_FOLDLEAN", "1")][..],
        ),
    ] {
        let (expected, _) = blastn(&dir, &query, &subject, task, 1, &[]);
        let mut on = vec![SHADOW, strips_shadow];
        on.extend_from_slice(more);
        let (report, stderr) = blastn(&dir, &query, &subject, task, 4, &on);
        assert!(report == expected, "{task}: the report differs");
        let at = format!("{task} with {strips_shadow:?}");
        // (`LOSAT_X_PAIRPARSHADOW` also names the shadow of TBLASTX's own `PAIRPAR`, whose
        // line says it compared nothing)
        let (_, skipped) = shadow_counts(&stderr, &[prefix, "[X_SHADOW] PAIRPAR "], &at);
        assert!(skipped > 0, "{at}: {stderr}");
        let line = stderr
            .lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("{at}: no strips shadow line:\n{stderr}"));
        assert!(
            !line.contains("seed chunks equal=0 "),
            "{at}: no chunk compared: {line}"
        );
    }
}

/// The `memo=` and `memofound=` counters of the shadow's line (`LOSAT_X_DIAGLV2`).
fn memo_counts(stderr: &str, at: &str) -> (u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] DIAGLV "))
        .unwrap_or_else(|| panic!("{at}: no shadow line:\n{stderr}"));
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("{at}: no {name} in: {line}"))
    };
    (field("memo="), field("memofound="))
}

// The memo of the latest insert (`LOSAT_X_DIAGLV2`, with the level bounds) and the lean
// instances of the dc-megablast strips folds (`LOSAT_X_FOLDLEAN`): dc-megablast at one thread
// (the reference fold) and at four (the strips fold of `LOSAT_X_PAIRPARDISC`) gives the bytes of
// the search with neither, each switch alone and both, at four threads also with the fold split
// over partitions (`LOSAT_X_PAIRPARFOLD`); with the shadow, which checks every answer of the memo,
// the memo answers, also for the latest insert's diagonal. megablast keeps no memo (a one-hit
// table).
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
        let on = [SHADOW, MEMO, LEAN, PAIRPARDISC];
        let at = format!("{task} {on:?} threads {threads}");
        let (report, stderr) = blastn(&dir, &query, &subject, task, threads, &on);
        assert!(report == expected, "{at}: the report differs");
        shadow_counts(&stderr, &[], &at);
        let (memo, found) = memo_counts(&stderr, &at);
        assert!(memo > 0 && found > 0 && found <= memo, "{at}: {stderr}");
    }
    let on = [SHADOW, MEMO];
    let (_, stderr) = blastn(&dir, &query, &subject, "megablast", 1, &on);
    let at = format!("megablast {on:?}");
    shadow_counts(&stderr, &[], &at);
    assert_eq!(memo_counts(&stderr, &at), (0, 0), "{at}: {stderr}");
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
// its cells and bounds from subject to subject with a growing offset.
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
    for task in ["megablast", "dc-megablast"] {
        check(&dir, &query, &subjects, task, &[]);
    }
}
