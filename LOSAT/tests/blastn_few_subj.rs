//! A `-subject` search with fewer 5 Mb work units than threads, and a subject batch of fewer chunks
//! than threads, searched one at a time on the whole pool (`LOSAT_X_BNFEWSUBJ`,
//! `src/algorithm/blastn/blast_engine/x_few_subj.rs`), on whole searches.
//!
//! The switch is in the strict set (on unless it is `0`). It chooses between two LOSAT shapes of the
//! search with the same result: the subjects (or the chunks of a batch) one at a time, each with the
//! subject scan in strips (`LOSAT_X_PAIRPAR`, `LOSAT_X_PAIRPARDISC`), or the subjects (the chunks) as
//! tasks of the pool without strips. So the property under test is that the report does not depend
//! on it or on the thread count: the search with the switch at 2, 3 and 4 threads gives the bytes
//! of the search at one thread and of the search without it, and the shadows of the strips
//! (`LOSAT_X_PAIRPARSHADOW`, `LOSAT_X_PAIRPARDISCSHADOW`, which compute the stage on copies and
//! compare the hit list and the whole diagonal hash with the reference code for every chunk) give
//! them too and compare the chunks. With `LOSAT_TIMING` the strips report how many chunks they
//! took, which shows which shape ran: the rule needs fewer units than threads (two subjects of one
//! unit each: 3 threads or more), and no lowercase mask with `-lcase_masking`. That those bytes
//! equal NCBI's is the business of the fixtures (`tests/blastn_regression_fixtures.py`).
//!
//! The query is LC738874 (287 kb, longer than 8000 letters: the strips take the diagonal hash
//! only). The subjects are LC738868 and LC738869 (two records), and one record of 5.2 Mb (LC738868,
//! pseudo-random bases, LC738869, pseudo-random bases: two chunks of `MAX_DBSEQ_LEN`), written per
//! call.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_BNFEWSUBJ", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const SHADOWS: &[(&str, &str)] = &[
    ("LOSAT_X_PAIRPARSHADOW", "1"),
    ("LOSAT_X_PAIRPARDISCSHADOW", "1"),
];

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

/// A FASTA file of `records` (name, letters) in a directory of this process, with a name of its
/// own per call (tests of this binary run in parallel).
fn subject_file(records: &[(&str, String)]) -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!("losat-blastn-few-subj-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!(
        "subject-{}.fasta",
        CALLS.fetch_add(1, Ordering::Relaxed)
    ));
    let mut text = String::new();
    for (name, letters) in records {
        text.push('>');
        text.push_str(name);
        text.push('\n');
        for line in letters.as_bytes().chunks(80) {
            text.push_str(std::str::from_utf8(line).unwrap());
            text.push('\n');
        }
    }
    std::fs::write(&path, text).unwrap();
    path
}

/// LC738868 and LC738869, two records of one 5 Mb unit each.
fn two_subjects() -> PathBuf {
    subject_file(&[
        ("LC738868", letters("LC738868.fasta")),
        ("LC738869", letters("LC738869.fasta")),
    ])
}

/// One record of 5,200,000 letters (two chunks): LC738868, pseudo-random bases, LC738869,
/// pseudo-random bases.
fn long_subject() -> PathBuf {
    let (first, second) = (letters("LC738868.fasta"), letters("LC738869.fasta"));
    let filler = 5_200_000 - first.len() - second.len();
    let mut state: u32 = 12345;
    let mut random = (0..filler).map(|_| {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345) & 0x7fff_ffff;
        b"ACGT"[((state >> 16) & 3) as usize] as char
    });
    let mut letters = first;
    letters.extend(random.by_ref().take(filler / 2));
    letters.push_str(&second);
    letters.extend(random);
    assert_eq!(letters.len(), 5_200_000);
    subject_file(&[("long", letters)])
}

/// One search of LC738874 against `subject`, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects nothing.
/// Returns the report and the standard error.
fn blastn(
    subject: &Path,
    task_args: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    blastn_query(
        &fasta("LC738874.fasta"),
        subject,
        task_args,
        threads,
        switches,
    )
}

/// `blastn` with another query file.
fn blastn_query(
    query: &Path,
    subject: &Path,
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
        .arg(subject)
        .args(task_args)
        .args(["-outfmt", "6"])
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

/// The chunks the strips took (`chunks=` of `[TIMING] x_pairpar_seed:` plus that of
/// `[TIMING] x_pairpardisc_seed:`) in `stderr`.
fn strip_chunks(stderr: &str) -> u64 {
    [
        "[TIMING] x_pairpar_seed: chunks=",
        "[TIMING] x_pairpardisc_seed: chunks=",
    ]
    .iter()
    .map(|prefix| {
        let line = stderr
            .lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("no {prefix} line:\n{stderr}"));
        line[prefix.len()..]
            .split_whitespace()
            .next()
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("no count in: {line}"))
    })
    .sum()
}

/// The chunks the strips' shadows compared (`seed chunks equal=N` of the two summary lines) and the
/// merges of their partitioned folds (`merges=M`), after checking that every other line of `stderr`
/// is a summary line of another shadow that the same variables turn on, without a difference.
fn shadow_chunks(stderr: &str) -> (u64, u64) {
    let field = |line: &str, name: &str| -> u64 {
        line.split(name)
            .nth(1)
            .and_then(|rest| rest.split_whitespace().next())
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    let (mut chunks, mut merges, mut summaries) = (0, 0, 0);
    for line in stderr.lines() {
        if line.starts_with("[X_PAIRPARSHADOW] blastn: ")
            || line.starts_with("[X_PAIRPARDISCSHADOW] blastn: ")
        {
            chunks += field(line, "seed chunks equal=");
            merges += field(line, " merges=");
            summaries += 1;
        } else {
            assert!(
                line.starts_with("[X_SHADOW] ") && line.ends_with("identical)")
                    || line.starts_with("[X_SHADOW] ") && line.ends_with(" identical"),
                "{stderr}"
            );
        }
    }
    assert_eq!(summaries, 2, "{stderr}");
    (chunks, merges)
}

/// For one task and one subject file: the one-thread search; the switch's searches at `threads`
/// with the chunks the strips took (`expect_strips`: more than zero, else zero); the search
/// without the switch at the largest count (no strips); the shadows at the largest count.
fn check(subject: &Path, task: &str, threads_and_strips: &[(u32, bool)]) {
    let task_args = ["-task", task];
    let (expected, expected_stderr) = blastn(subject, &task_args, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for &(threads, expect_strips) in threads_and_strips {
        let (report, stderr) = blastn(subject, &task_args, threads, &[TIMING]);
        assert!(
            report == expected,
            "{task} threads {threads}: the report differs from one thread's"
        );
        let chunks = strip_chunks(&stderr);
        assert_eq!(
            chunks > 0,
            expect_strips,
            "{task} threads {threads}: strips took {chunks} chunks"
        );
    }
    let most = threads_and_strips.iter().map(|&(t, _)| t).max().unwrap();
    let (report, stderr) = blastn(subject, &task_args, most, &[OFF, TIMING]);
    assert!(
        report == expected,
        "{task} threads {most}: the report differs with {OFF:?}"
    );
    assert_eq!(strip_chunks(&stderr), 0, "{task} threads {most} {OFF:?}");
    let (report, stderr) = blastn(subject, &task_args, most, SHADOWS);
    assert!(
        report == expected,
        "{task} threads {most}: the report differs with the shadows"
    );
    let (chunks, merges) = shadow_chunks(&stderr);
    assert!(chunks > 0, "{task}: {stderr}");
    // LOSAT_X_PAIRPARFOLD (strict set) splits each chunk's fold over the pool: every chunk the
    // shadows compared went through a partitioned fold, the later ones from the diagonal hash that
    // the chunks (and subjects) before them filled.
    assert_eq!(merges, chunks, "{task}: {stderr}");
}

// Two subjects of one unit each: the subjects one at a time from 3 threads on (fewer units than
// threads); at 2 threads (as many units as threads) they stay tasks of the pool.
#[test]
fn two_subjects_dc_megablast() {
    check(
        &two_subjects(),
        "dc-megablast",
        &[(2, false), (3, true), (4, true)],
    );
}

#[test]
fn two_subjects_megablast() {
    check(
        &two_subjects(),
        "megablast",
        &[(2, false), (3, true), (4, true)],
    );
}

#[test]
fn two_subjects_blastn() {
    check(&two_subjects(), "blastn", &[(2, false), (4, true)]);
}

// One subject of two chunks: the chunks one at a time from 3 threads on (fewer chunks than
// threads); at 2 threads they stay tasks of the pool.
#[test]
fn a_subject_of_two_chunks_dc_megablast() {
    check(&long_subject(), "dc-megablast", &[(2, false), (4, true)]);
}

#[test]
fn a_subject_of_two_chunks_megablast() {
    check(&long_subject(), "megablast", &[(3, true)]);
}

// With `-lcase_masking` the rule keeps the tasks of the pool (a masked chunk has no strips): no
// strips with or without lowercase letters, and the one-thread report.
#[test]
fn lowercase_masking_keeps_the_subject_tasks() {
    let first = letters("LC738868.fasta");
    let mut masked = first[..1000].to_ascii_lowercase();
    masked.push_str(&first[1000..]);
    let subject = subject_file(&[
        ("LC738868", masked),
        ("LC738869", letters("LC738869.fasta")),
    ]);
    let task_args = ["-task", "megablast", "-lcase_masking"];
    let (expected, _) = blastn(&subject, &task_args, 1, &[]);
    assert!(!expected.is_empty());
    let (report, stderr) = blastn(&subject, &task_args, 4, &[TIMING]);
    assert!(report == expected, "the report differs at 4 threads");
    assert_eq!(strip_chunks(&stderr), 0, "{stderr}");
}

// Both shapes in one search: the 5.2 Mb record (two chunks) and LC738869, three units: at 4 threads
// the subjects one at a time and the long subject's chunks one at a time, every chunk with strips;
// each later chunk's partitioned fold starts from the diagonal hash that the chunks before it
// filled (one scratch for the whole loop, as at one thread).
#[test]
fn a_long_subject_and_another_take_both_shapes() {
    let (first, second) = (letters("LC738868.fasta"), letters("LC738869.fasta"));
    let long = std::fs::read_to_string(long_subject())
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>();
    let subject = subject_file(&[("long", long), ("LC738869", second), ("LC738868", first)]);
    for task in ["dc-megablast", "megablast"] {
        check(&subject, task, &[(4, true)]);
        let (_, stderr) = blastn(&subject, &["-task", task], 4, SHADOWS);
        let (chunks, merges) = shadow_chunks(&stderr);
        assert_eq!((chunks, merges), (4, 4), "{task}: {stderr}");
    }
}

// The per-query gapped stages (LOSAT_X_QGROUPPAR, strict set) run on the chunks after the first of
// a subject searched by its owner (before, only on one-chunk batches): 27 queries of 10 kb cut from
// LC738868 (3,000 letters on, as `blastn_qgroup_par.rs` cuts them) against a record of 4,990,000
// pseudo-random bases and LC738868 (two chunks; the hits in the second), blastn at 4 threads, with
// the stages' shadow comparing them to the owner loop.
#[test]
fn per_query_stages_on_a_later_chunk_match_their_shadow() {
    let first = letters("LC738868.fasta");
    let pieces: Vec<(String, String)> = (0..27)
        .map(|i| {
            let start = 3000 + i * 10_000;
            (format!("p{i}"), first[start..start + 10_000].to_string())
        })
        .collect();
    let pieces: Vec<(&str, String)> = pieces
        .iter()
        .map(|(name, letters)| (name.as_str(), letters.clone()))
        .collect();
    let query = subject_file(&pieces);
    let mut state: u32 = 12345;
    let mut letters_of_tail: String = (0..4_990_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345) & 0x7fff_ffff;
            b"ACGT"[((state >> 16) & 3) as usize] as char
        })
        .collect();
    letters_of_tail.push_str(&first);
    let subject = subject_file(&[("tail", letters_of_tail)]);
    let task_args = ["-task", "blastn"];
    let (expected, _) = blastn_query(&query, &subject, &task_args, 1, &[]);
    assert!(!expected.is_empty());
    let (report, stderr) = blastn_query(
        &query,
        &subject,
        &task_args,
        4,
        &[("LOSAT_X_QGROUPPARSHADOW", "1"), TIMING],
    );
    assert!(report == expected, "the report differs with the shadow");
    assert!(strip_chunks(&stderr) >= 2, "{stderr}");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] QGROUPPAR prelim="))
        .unwrap_or_else(|| panic!("no QGROUPPAR line:\n{stderr}"));
    assert!(line.ends_with(" identical"), "{line}");
    let prelim: u64 = line["[X_SHADOW] QGROUPPAR prelim=".len()..]
        .split_whitespace()
        .next()
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("{line}"));
    assert!(prelim >= 2, "{line}");
}

// A query split into query chunks (`CHUNK_SIZE`, as NCBI reads it): the batch's own call has no
// lookup table (its traceback only) and keeps the subjects as tasks of the pool; each query chunk's
// call takes the one-subject shape with strips. The report is the one-thread one.
#[test]
fn query_chunks_take_the_shape_and_the_batch_keeps_its_traceback() {
    let subject = two_subjects();
    let task_args = ["-task", "megablast"];
    let split = ("CHUNK_SIZE", "100000");
    let (expected, _) = blastn(&subject, &task_args, 1, &[split]);
    assert!(!expected.is_empty());
    let (report, stderr) = blastn(&subject, &task_args, 4, &[split, TIMING]);
    assert!(report == expected, "the report differs at 4 threads");
    assert!(strip_chunks(&stderr) >= 4, "{stderr}");
    let (report, stderr) = blastn(&subject, &task_args, 4, &[split, OFF, TIMING]);
    assert!(report == expected, "the report differs with {OFF:?}");
    assert_eq!(strip_chunks(&stderr), 0, "{stderr}");
}
