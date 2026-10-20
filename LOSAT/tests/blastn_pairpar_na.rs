//! The one-hit seed stage of blastn's tables that are not two-stage on the search pool
//! (`LOSAT_X_PAIRPARNA`, `src/algorithm/blastn/blast_engine/x_pair_seed_na.rs`) on whole
//! searches.
//!
//! The switch is in the strict set (on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not
//! set). It chooses between two LOSAT code paths with the same result: the subject scan cut
//! into strips that the pool threads scan, and the diagonal-hash fold in the scan order (by the
//! owner, or split by diagonal-hash bucket with `LOSAT_X_PAIRPARFOLD`), or the scan and fold of
//! the reference code. So the property under test is that the report does not depend on it:
//! `-task blastn-short` (word 7) and word sizes 6 and 8 with the switch at 1, 2 and 4 threads
//! give the bytes of the search at one thread (where the stage does not act), and so do the
//! switch at `0` and `LOSAT_X_DEFAULT=0` at 4 threads, and the shadow run (`LOSAT_X_PAIRPARNASHADOW`,
//! which computes the stage on copies and compares the hit list and the whole diagonal hash
//! with the reference code for every chunk) gives them too and compares at least one chunk.
//! That those bytes equal NCBI's is the business of the blastn-short fixtures
//! (`tests/blastn_regression_fixtures.py`).
//!
//! The input is one subject (the stage takes a search of one subject by the owner of a pool of
//! two threads or more) and 820 short queries cut from it, 19,673 letters (the stage takes the
//! diagonal hash only, a query block longer than 8000 letters): 24-mers on both strands,
//! queries of 18 to 30 letters, and queries with a low-complexity half that DUST masks (so
//! `s_TypeOfWord` checks the masks).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const SWITCH: (&str, &str) = ("LOSAT_X_PAIRPARNA", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_PAIRPARNASHADOW", "1");
const OFF: (&str, &str) = ("LOSAT_X_PAIRPARNA", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// The letters of the subject fixture LC738868 (306,008 bases, no ambiguity).
fn subject_letters() -> Vec<u8> {
    let text = std::fs::read(fixture("LC738868.fasta")).expect("subject fixture");
    text.split(|&b| b == b'\n')
        .skip(1)
        .flatten()
        .copied()
        .filter(u8::is_ascii_alphabetic)
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// A query file written once per test binary: `count` queries cut from the subject at evenly
/// spaced offsets, of the lengths `length(i)`, every third one reverse-complemented, and every
/// fortieth with its first half replaced by a dinucleotide repeat (DUST masks it).
fn query_file(name: &str, count: usize, length: impl Fn(usize) -> usize) -> PathBuf {
    let subject = subject_letters();
    let step = (subject.len() - 40) / count;
    let mut text = String::new();
    for i in 0..count {
        let len = length(i);
        let mut seq = subject[i * step..i * step + len].to_vec();
        if i % 3 == 1 {
            seq = reverse_complement(&seq);
        }
        if i % 40 == 7 {
            for (k, b) in seq.iter_mut().take(len / 2).enumerate() {
                *b = if k % 2 == 0 { b'A' } else { b'T' };
            }
        }
        text.push_str(&format!(
            ">{name}_{i}\n{}\n",
            String::from_utf8(seq).unwrap()
        ));
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("blastn_pairpar_na");
    std::fs::create_dir_all(&dir).expect("work directory");
    let file = dir.join(format!("{name}.fa"));
    std::fs::write(&file, text).expect("query file");
    file
}

/// 600 24-mers, then 200 queries of 18 to 30 letters, then 20 more 24-mers (19,673 letters).
fn queries() -> &'static PathBuf {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| {
        query_file("short", 820, |i| {
            if (600..800).contains(&i) {
                18 + (i * 7) % 13
            } else {
                24
            }
        })
    })
}

/// 100 24-mers (2,400 letters: a query block under 8000 letters takes the diagonal array).
fn few_queries() -> &'static PathBuf {
    static FILE: OnceLock<PathBuf> = OnceLock::new();
    FILE.get_or_init(|| query_file("few", 100, |_| 24))
}

/// One search of `query` against LC738868, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built
/// with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that
/// file name as soon as a search has more than one thread, and it selects nothing. Returns the
/// report and the standard error.
fn blastn(
    query: &Path,
    task_args: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fixture("LC738868.fasta"))
        .args(task_args)
        .args(["-evalue", "0.001", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task_args:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The `chunks=` count of the stage's `LOSAT_TIMING` line (`[TIMING] x_pairparna_seed:
/// chunks=N ...`) in `stderr`.
fn stage_chunks(stderr: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] x_pairparna_seed: chunks="))
        .unwrap_or_else(|| panic!("no x_pairparna_seed line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|word| word.strip_prefix("chunks="))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no chunks= in: {line}"))
}

/// For one word size: the search without a variable at one thread, then with the switch at 1,
/// 2 and 4 threads, without a variable, with the owner's fold (`LOSAT_X_PAIRPARFOLD=0`), with
/// the switch at `0` and with `LOSAT_X_DEFAULT=0` at 4 threads (the same bytes, nothing on
/// standard error), then the shadow at 4 threads and at 1 thread (the same bytes, and its
/// summary line with at least one chunk compared).
fn check(task_args: &[&str]) {
    let query = queries();
    let (expected, expected_stderr) = blastn(query, task_args, 1, &[]);
    assert!(!expected.is_empty(), "{task_args:?}: no hits");
    assert_eq!(expected_stderr, "", "{task_args:?}");
    let owner_fold = [SWITCH, ("LOSAT_X_PAIRPARFOLD", "0")];
    for (threads, switches) in [
        (1, &[SWITCH][..]),
        (2, &[SWITCH][..]),
        (4, &[SWITCH][..]),
        (4, &[][..]),
        (4, &owner_fold[..]),
        (4, &[OFF][..]),
        (4, &[DEFAULT_OFF][..]),
    ] {
        let (report, stderr) = blastn(query, task_args, threads, switches);
        assert!(
            report == expected,
            "{task_args:?} threads {threads}: the report differs with {switches:?}"
        );
        assert_eq!(stderr, "", "{task_args:?} threads {threads} {switches:?}");
    }
    for threads in [4, 1] {
        let (report, stderr) = blastn(query, task_args, threads, &[SHADOW]);
        assert!(
            report == expected,
            "{task_args:?} threads {threads}: the report differs with {SHADOW:?}"
        );
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[X_PAIRPARNASHADOW] blastn: seed chunks equal="))
            .unwrap_or_else(|| panic!("{task_args:?}: no summary line:\n{stderr}"));
        assert!(
            !line.starts_with("[X_PAIRPARNASHADOW] blastn: seed chunks equal=0 "),
            "{task_args:?} threads {threads}: {line}"
        );
        assert_eq!(stderr.lines().count(), 1, "{task_args:?}: {stderr}");
    }
}

#[test]
fn blastn_short_gives_the_same_report() {
    check(&["-task", "blastn-short"]);
}

// Word sizes 6 and 8 with `-task blastn-short`: a direct table whose lookup width is the word
// size (`coordination.rs` `choose_na_lookup_table`: width 8 for word 8 from 8,500 table
// entries on, width 6 for word 6).
#[test]
fn word_sizes_6_and_8_give_the_same_report() {
    check(&["-task", "blastn-short", "-word_size", "6"]);
    check(&["-task", "blastn-short", "-word_size", "8"]);
}

// The switch-on runs of `check` compare bytes only, which a stage that silently stopped
// running would pass too: with LOSAT_TIMING the stage prints how many chunks it took. It runs
// by default (no variable) and with the switch set; with the switch at `0` or
// `LOSAT_X_DEFAULT=0` it does not (no line). It takes none at one thread (no pool) and none for
// a query block under 8000 letters (the diagonal array), and the report is the same there.
#[test]
fn the_stage_runs_only_in_its_scope() {
    let task_args = ["-task", "blastn-short"];
    let timing = ("LOSAT_TIMING", "1");
    let (expected, _) = blastn(queries(), &task_args, 1, &[]);
    for switches in [&[timing][..], &[SWITCH, timing][..]] {
        let (report, stderr) = blastn(queries(), &task_args, 4, switches);
        assert!(report == expected, "the report differs with {switches:?}");
        assert!(
            stage_chunks(&stderr) > 0,
            "{switches:?}: the stage took no chunk:\n{stderr}"
        );
    }
    for switches in [[OFF, timing], [DEFAULT_OFF, timing]] {
        let (report, stderr) = blastn(queries(), &task_args, 4, &switches);
        assert!(report == expected, "the report differs with {switches:?}");
        assert!(
            !stderr.contains("x_pairparna_seed"),
            "{switches:?}: the stage ran:\n{stderr}"
        );
    }
    let (report, stderr) = blastn(queries(), &task_args, 1, &[SWITCH, timing]);
    assert!(report == expected, "the report differs at one thread");
    assert_eq!(stage_chunks(&stderr), 0, "{stderr}");
    let (expected, _) = blastn(few_queries(), &task_args, 1, &[]);
    assert!(!expected.is_empty(), "few queries: no hits");
    let (report, stderr) = blastn(few_queries(), &task_args, 4, &[SWITCH, timing]);
    assert!(report == expected, "few queries: the report differs");
    assert_eq!(stage_chunks(&stderr), 0, "{stderr}");
}
