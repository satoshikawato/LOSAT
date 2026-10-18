//! NCBI's query batches searched side by side ahead of blastn's batch loop
//! (`LOSAT_X_QBATCHPAR`, `src/algorithm/blastn/blast_engine/x_qbatch_par.rs`) on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses between
//! two LOSAT code paths with the same result: the loop searching every batch itself, or slots
//! searching the planned batches ahead of it, whose results the loop takes only under the key it
//! computes itself. So the property under test is that nothing the run writes depends on it:
//! standard output (or the `-out` file), standard error (byte for byte without `LOSAT_X_STATS`;
//! with it, without the `[X_STATS]` lines it adds) and the exit status equal those of the same
//! search without the switch, for every task, outfmt 0/6/7, fixed and adaptive batch sizes (with
//! the mixer's alternation of tiny and full batches, and a batch split into query chunks),
//! `-max_target_seqs 1`, warnings of later batches (the reader's, and an invalid query's),
//! `-query_loc` skipping records, `LOSAT_TIMING` blocks, an error at a batch with slots searching
//! past it (each kind of error of a batch), corrupted predictions (`LOSAT_X_QBATCHPARMISS`, two
//! slots), dc-megablast's `LOSAT_X_PAIRPARDISC` inside a slot's pool, and the shadow
//! (`LOSAT_X_QBATCHPARSHADOW`), which also compares every batch taken from a slot with the loop's
//! own search of it. With `LOSAT_X_STATS` the summary line (`[X_STATS] QBATCHPAR ...`) shows that
//! the slots' path ran (`precomputed=` batches taken from a slot), or why it did not, and which
//! layout the slot model judged (`units=`, `rule_k=`, `rule_m=`, `predicted=`). That the bytes
//! equal NCBI's is the business of the BLASTN fixtures.
//!
//! The slot model (`choose_slots`) keeps the loop on many of the small megablast and
//! dc-megablast inputs here at four threads, and where it chooses slots for them it does so by a
//! few percent (predicted 0.93-0.96 of the loop's span against its 0.97 gate, depending on the
//! parity of the batch count), so the tests whose property needs the slots' path on those inputs
//! set two slots (`SLOTS2`) or four; the `-task blastn` tests and
//! `the_rule_chooses_the_slots_on_whole_runs` (which also runs megablast and dc-megablast where
//! the model itself starts slots, at eight threads) run on the model's own choice.
//!
//! The queries are 32 pieces of LC738874 and the subject LC738868 (one subject, so its 5 Mb work
//! units are fewer than the threads). `-task blastn` runs on pieces of LC738884 with one query
//! per batch, and on both sets of pieces with several queries per batch (the adaptive batch size
//! and fixed ones).

#![cfg(feature = "parallel")]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_QBATCHPARSHADOW", "1");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
/// Two slots (`LOSAT_X_QBATCHPARSLOTS`, tuning), whatever the slot model would choose.
const SLOTS2: (&str, &str) = ("LOSAT_X_QBATCHPARSLOTS", "2");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of a one-record FASTA fixture.
fn genome(name: &str) -> Vec<u8> {
    let text = std::fs::read(fixture(name)).expect("read fixture");
    text.split(|&byte| byte == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A file name in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_qbatchpar_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A FASTA file in the temporary directory, removed when dropped.
struct TempFasta(PathBuf);

impl TempFasta {
    fn new(tag: &str, records: &[(String, Vec<u8>)]) -> Self {
        let path = temp_path(&format!("{tag}.fa"));
        let mut text = Vec::new();
        for (title, letters) in records {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            for line in letters.chunks(70) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        }
        std::fs::write(&path, text).expect("write temporary FASTA");
        TempFasta(path)
    }

    /// The file with `text` appended after the records.
    fn with_tail(tag: &str, records: &[(String, Vec<u8>)], tail: &str) -> Self {
        let file = TempFasta::new(tag, records);
        let mut text = std::fs::read(&file.0).expect("read temporary FASTA");
        text.extend_from_slice(tail.as_bytes());
        std::fs::write(&file.0, text).expect("write temporary FASTA");
        file
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// 32 queries cut from LC738874, 1,500 to 3,900 letters each (about 86 kb).
fn queries() -> Vec<(String, Vec<u8>)> {
    pieces("LC738874.fasta")
}

/// 32 pieces of a genome fixture, as `queries`.
fn pieces(name: &str) -> Vec<(String, Vec<u8>)> {
    let letters = genome(name);
    (0..32)
        .map(|index| {
            let start = index * 8_000;
            let length = 1_500 + (index * 373) % 2_400;
            (
                format!("q{index} piece {start}"),
                letters[start..start + length].to_vec(),
            )
        })
        .collect()
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

impl Run {
    /// Standard error without the counter lines of `LOSAT_X_STATS` and the shadow's line.
    fn stderr_without_summary(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| {
                !line.starts_with("[X_STATS] ") && !line.starts_with("[X_SHADOW] QBATCHPAR")
            })
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// The fields of the summary line `prefix ...` (`name=value`).
    fn summary(&self, prefix: &str) -> HashMap<String, String> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("no {prefix} line:\n{}", self.stderr));
        line.split_whitespace()
            .filter_map(|word| word.split_once('='))
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    }

    fn stat(&self, name: &str) -> u64 {
        let fields = self.summary("[X_STATS] QBATCHPAR");
        fields
            .get(name)
            .unwrap_or_else(|| panic!("no {name} in the summary:\n{}", self.stderr))
            .parse()
            .expect("a number")
    }
}

/// `LOSAT_X_QBATCHPARMEM` as the slot model's default budget was before 2026-10-17 for the task
/// and the thread count of `args`: 96 MiB (blastn, blastn-short) or 28 MiB (megablast,
/// dc-megablast) per thread beyond the first. Since then a query file that NCBI's blastn would not
/// split by queries (all the files here) gets no default budget (`x_qbatch_par.rs`,
/// `Profile::of_search`); these tests check the slots, so their runs keep the earlier budget,
/// unless a test sets the variable itself (`env` comes after it).
fn earlier_budget(args: &[&str]) -> String {
    let value = |flag: &str| {
        args.iter()
            .position(|arg| *arg == flag)
            .and_then(|index| args.get(index + 1))
            .copied()
    };
    let threads: u64 = value("-num_threads")
        .and_then(|threads| threads.parse().ok())
        .unwrap_or(1);
    let per_thread = match value("-task").unwrap_or("megablast") {
        "megablast" | "dc-megablast" => 28,
        _ => 96,
    };
    (per_thread * threads.saturating_sub(1)).to_string()
}

/// One `LOSAT blastn` run with `args`, with `env` as its only variables that choose a code path
/// or NCBI's batches (and the earlier default budget, `earlier_budget`). `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn blastn(args: &[&str], env: &[(&str, &str)]) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .args(args)
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .env("LOSAT_X_QBATCHPARMEM", earlier_budget(args))
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// The run with `switches` (and `LOSAT_X_STATS`, for its counters) gives what the run without
/// them gives (standard output, standard error without the counter and shadow lines, exit
/// status); returns it.
fn same_as_without(args: &[&str], env: &[(&str, &str)], switches: &[(&str, &str)]) -> Run {
    let reference = blastn(args, env);
    let mut with = env.to_vec();
    with.extend_from_slice(switches);
    with.push(STATS);
    let run = blastn(args, &with);
    let what = format!("{args:?} {env:?} {switches:?}");
    assert_eq!(
        run.code, reference.code,
        "{what}: exit status\n{}",
        run.stderr
    );
    assert!(
        run.stdout == reference.stdout,
        "{what}: standard output differs ({} and {} bytes)",
        run.stdout.len(),
        reference.stdout.len()
    );
    assert_eq!(
        run.stderr_without_summary(),
        reference.stderr,
        "{what}: standard error"
    );
    // No batch a slot searched panicked (kept or dropped), and the shadow found no difference.
    if run.stderr.contains("[X_STATS] QBATCHPAR") {
        assert_eq!(run.stat("panicked"), 0, "{what}: {}", run.stderr);
    }
    if run.stderr.contains("[X_SHADOW] QBATCHPAR") {
        let shadow = run.summary("[X_SHADOW] QBATCHPAR");
        assert_eq!(
            (shadow["diffs"].as_str(), shadow["panicked"].as_str()),
            ("0", "0"),
            "{what}"
        );
    }
    run
}

fn path(file: &TempFasta) -> &str {
    file.0.to_str().expect("UTF-8 temporary path")
}

// Megablast and dc-megablast with each outfmt, with the adaptive batch size (its first batch is
// 5,000 letters, so two batches here) and fixed ones; with `BATCH_SIZE=1` every query is a batch
// and the slots take all but the first.
#[test]
fn every_outfmt_and_batch_size_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &queries());
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    for task in ["megablast", "dc-megablast"] {
        for outfmt in ["0", "6", "7"] {
            for batch in [None, Some("1"), Some("7000")] {
                if outfmt != "6" && batch == Some("7000") {
                    continue;
                }
                let args = [
                    "-task",
                    task,
                    "-query",
                    path(&query),
                    "-subject",
                    subject,
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    "4",
                ];
                let env: Vec<(&str, &str)> =
                    batch.map(|size| ("BATCH_SIZE", size)).into_iter().collect();
                let switches: &[(&str, &str)] = if batch == Some("1") {
                    &[SWITCH, SLOTS2]
                } else {
                    &[SWITCH]
                };
                let run = same_as_without(&args, &env, switches);
                if batch == Some("1") {
                    assert_eq!(run.stat("batches"), 32, "{task} {outfmt}");
                    assert!(
                        run.stat("precomputed") > 0,
                        "{task} {outfmt}: {}",
                        run.stderr
                    );
                }
            }
        }
    }
}

// The mixer's alternation of tiny and full batches (dc-megablast: a batch with few extensions
// per letter makes the next size overflow to the smallest, 100 letters, one query, and a small
// batch is followed by a full one). The queries are pieces of E. coli, unrelated to the
// subject, so their batches have few extensions. `CHUNK_SIZE` makes the full batch 20,000
// letters (the mixer's maximum is the chunk size less 1000) without splitting it into query
// chunks.
#[test]
fn the_mixers_tiny_and_full_batches_are_taken_from_the_slots() {
    let query = TempFasta::new("q", &pieces("MG1655.fna"));
    let subject = fixture("LC738868.fasta");
    let args = [
        "-task",
        "dc-megablast",
        "-query",
        path(&query),
        "-subject",
        subject.to_str().expect("UTF-8 path"),
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let env = [("CHUNK_SIZE", "21000")];
    let run = same_as_without(&args, &env, &[SWITCH, ("LOSAT_X_QBATCHPARSLOTS", "2")]);
    assert!(run.stat("batches") >= 6, "{}", run.stderr);
    assert!(run.stat("precomputed") > 0, "{}", run.stderr);
    let run = same_as_without(&args, &env, &[SHADOW, SLOTS2]);
    let shadow = run.summary("[X_SHADOW] QBATCHPAR");
    assert_ne!(shadow["compared"], "0", "{}", run.stderr);
    assert_eq!(shadow["diffs"], "0");
}

// `-max_target_seqs 1` against three subjects (three work units, below the four threads).
#[test]
fn one_target_sequence_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &queries());
    let subjects: Vec<(String, Vec<u8>)> = ["LC738868.fasta", "LC738869.fasta", "LC738870.fasta"]
        .iter()
        .map(|name| (name.to_string(), genome(name)[..60_000].to_vec()))
        .collect();
    let subject = TempFasta::new("s", &subjects);
    for task in ["megablast", "dc-megablast"] {
        for outfmt in ["6", "7"] {
            let args = [
                "-task",
                task,
                "-query",
                path(&query),
                "-subject",
                path(&subject),
                "-max_target_seqs",
                "1",
                "-outfmt",
                outfmt,
                "-num_threads",
                "4",
            ];
            let run = same_as_without(&args, &[("BATCH_SIZE", "3000")], &[SWITCH, SLOTS2]);
            assert!(
                run.stat("precomputed") > 0,
                "{task} {outfmt}: {}",
                run.stderr
            );
        }
    }
}

/// The searches of `query` against `subject` with each outfmt (outfmt 0 also into an `-out`
/// file, whose outfmt 0 prolog is written with the reports), with and without the switch, the
/// corrupted predictions and the shadow; they must fail as without the switch.
fn errors_as_without_the_switch(query: &TempFasta, subject: &str, extra: &[&str], batch: &str) {
    let out = temp_path("out.txt");
    let out = out.to_str().expect("UTF-8 path").to_string();
    for outfmt in ["0", "6", "7"] {
        for to_file in [false, true] {
            if to_file && outfmt != "0" {
                continue;
            }
            let mut args = vec![
                "-query",
                path(query),
                "-subject",
                subject,
                "-outfmt",
                outfmt,
            ];
            args.extend_from_slice(extra);
            args.extend_from_slice(&["-num_threads", "4"]);
            if to_file {
                args.extend_from_slice(&["-out", &out]);
            }
            let env = [("BATCH_SIZE", batch)];
            let reference = blastn(&args, &env);
            assert_ne!(reference.code, Some(0), "{args:?}: the search must fail");
            let reference_file = to_file.then(|| std::fs::read(&out).expect("read -out"));
            for switches in [
                vec![SWITCH],
                vec![
                    SWITCH,
                    ("LOSAT_X_QBATCHPARSLOTS", "2"),
                    ("LOSAT_X_QBATCHPARMISS", "2"),
                ],
                vec![SHADOW, SLOTS2],
            ] {
                let mut with = env.to_vec();
                with.extend_from_slice(&switches);
                let run = blastn(&args, &with);
                let what = format!("{args:?} {switches:?}");
                assert_eq!(run.code, reference.code, "{what}\n{}", run.stderr);
                assert!(run.stdout == reference.stdout, "{what}: standard output");
                // Without LOSAT_X_STATS the switch writes nothing of its own: the same bytes
                // (the shadow writes its line).
                if switches.contains(&SHADOW) {
                    assert_eq!(run.stderr_without_summary(), reference.stderr, "{what}");
                } else {
                    assert_eq!(run.stderr, reference.stderr, "{what}");
                }
                if let Some(expected) = &reference_file {
                    assert!(
                        std::fs::read(&out).expect("read -out") == *expected,
                        "{what}: the -out file differs"
                    );
                }
            }
        }
    }
    let _ = std::fs::remove_file(&out);
}

// An error at a batch: the loop meets it at the same batch, after the reports of the batches
// before, as without the switch. The reader error after the last record, a last batch of
// queries without letters, and a last batch of records that `-query_loc` skips (the batches a
// batch size reaches end there: skipped records and records without letters add no letters).
#[test]
fn an_error_at_the_last_batch_comes_as_without_the_switch() {
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    let records = queries();
    let reader_error = TempFasta::with_tail("bad", &records, ">bad\n@@@@\n");
    errors_as_without_the_switch(&reader_error, subject, &[], "3000");
    let mut empty = records.clone();
    empty.push(("e1 empty".to_string(), Vec::new()));
    empty.push(("e2 empty".to_string(), Vec::new()));
    let empty = TempFasta::new("empty", &empty);
    errors_as_without_the_switch(&empty, subject, &[], "1");
    let mut short = records.clone();
    short.push(("s1 short".to_string(), records[0].1[..400].to_vec()));
    short.push(("s2 short".to_string(), records[1].1[..300].to_vec()));
    let short = TempFasta::new("short", &short);
    errors_as_without_the_switch(&short, subject, &["-query_loc", "1001-1400"], "1");
}

// An error inside a batch's search, with slots searching the batches after it: subjects with
// fewer letters than subjects stop the first batch with a valid query (NCBI's "average subject
// length" error), here the second batch, after the report of the first (a query of N only,
// which is invalid and not searched). The batches the slots searched past it are dropped.
#[test]
fn an_error_inside_a_batch_with_slots_past_it_comes_as_without_the_switch() {
    let mut records = vec![("n40 only N".to_string(), vec![b'N'; 40])];
    records.extend(queries());
    let query = TempFasta::new("q", &records);
    let subject = TempFasta::new(
        "s",
        &[
            ("s1 one".to_string(), b"A".to_vec()),
            ("s2 none".to_string(), Vec::new()),
        ],
    );
    errors_as_without_the_switch(&query, path(&subject), &[], "1");
    let args = [
        "-query",
        path(&query),
        "-subject",
        path(&subject),
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SWITCH, SLOTS2]);
    assert_eq!(run.stat("precomputed"), 1, "{}", run.stderr);
    assert_eq!(run.stat("batches"), 1, "{}", run.stderr);
}

// Two slots with corrupted predictions: the slots' batches under keys the loop never computes
// are discarded, the speculation stops (two misses within eight batches, or the waste bound),
// and the reports do not change.
#[test]
fn mispredicted_batches_are_discarded() {
    let query = TempFasta::new("q", &queries());
    let subject = fixture("LC738868.fasta");
    let args = [
        "-query",
        path(&query),
        "-subject",
        subject.to_str().expect("UTF-8 path"),
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    // With two slots a plan holds four batches, so three predictions; every second (or fifth)
    // corrupted makes the plan's next key wrong at every second (fifth) batch.
    for every in ["2", "5"] {
        let run = same_as_without(
            &args,
            &[("BATCH_SIZE", "1")],
            &[
                SWITCH,
                ("LOSAT_X_QBATCHPARSLOTS", "2"),
                ("LOSAT_X_QBATCHPARMISS", every),
            ],
        );
        assert!(run.stat("mispredicted") > 0, "{every}: {}", run.stderr);
        // Two misses within eight batches, or the waste bound if the discards reach it first.
        let off = run.summary("[X_STATS] QBATCHPAR")["off"].clone();
        assert!(off == "guard" || off == "waste", "{}", run.stderr);
        let run = same_as_without(
            &args,
            &[("BATCH_SIZE", "1")],
            &[SHADOW, SLOTS2, ("LOSAT_X_QBATCHPARMISS", every)],
        );
        let shadow = run.summary("[X_SHADOW] QBATCHPAR");
        assert_eq!(shadow["diffs"], "0");
    }
}

// The shadow compares every batch it takes from a slot with the loop's own search of it.
#[test]
fn the_shadow_compares_every_batch_taken_from_a_slot() {
    let query = TempFasta::new("q", &queries());
    let subject = fixture("LC738868.fasta");
    for task in ["megablast", "dc-megablast"] {
        let args = [
            "-task",
            task,
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ];
        // Two slots of two threads, and four of one thread (a plain thread with a one-thread
        // pool: PAIRPARDISC and DISCLUTPAR take their serial branch there).
        for slots in [SLOTS2, ("LOSAT_X_QBATCHPARSLOTS", "4")] {
            let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SHADOW, slots]);
            let shadow = run.summary("[X_SHADOW] QBATCHPAR");
            let what = format!("{task} {slots:?}: {}", run.stderr);
            assert_eq!(shadow["diffs"], "0", "{what}");
            assert_eq!(shadow["compared"], run.stat("precomputed").to_string());
            assert_eq!(shadow["compared"], "31", "{what}");
            let k: u64 = slots.1.parse().expect("a slot count");
            assert_eq!((run.stat("k"), run.stat("m")), (k, 4 / k), "{what}");
        }
    }
}

/// The queries with warnings in later batches: the reader's (invalid residues, record 12) and an
/// invalid query's (N only, record 8; not searched), as NCBI writes them before the report of
/// the query.
fn queries_with_warnings() -> TempFasta {
    let mut records = queries();
    records[8] = ("q8 only N".to_string(), vec![b'N'; 60]);
    let mut text = Vec::new();
    for (index, (title, letters)) in records.iter().enumerate() {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for (line, chunk) in letters.chunks(70).enumerate() {
            text.extend_from_slice(chunk);
            if index == 12 && line == 3 {
                text.extend_from_slice(b"123");
            }
            text.push(b'\n');
        }
    }
    let file = TempFasta::new("warn", &[]);
    std::fs::write(&file.0, text).expect("write temporary FASTA");
    file
}

/// `text` with the seconds of `LOSAT_TIMING` lines (`0.123s`) replaced by `#s`, and without the
/// end-of-run lines of other switches (`[TIMING] x_...`: their work counters count what each
/// thread did, so they differ when slots search the batches).
fn without_seconds(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("[TIMING] x_"))
        .map(|line| {
            if !line.starts_with("[TIMING]") {
                return format!("{line}\n");
            }
            let bytes = line.as_bytes();
            let mut masked = String::new();
            let mut at = 0;
            while at < bytes.len() {
                let digits = bytes[at..]
                    .iter()
                    .take_while(|b| b.is_ascii_digit())
                    .count();
                if digits > 0 && bytes.get(at + digits) == Some(&b'.') {
                    let fraction = bytes[at + digits + 1..]
                        .iter()
                        .take_while(|b| b.is_ascii_digit())
                        .count();
                    if fraction > 0 && bytes.get(at + digits + 1 + fraction) == Some(&b's') {
                        masked.push_str("#s");
                        at += digits + 2 + fraction;
                        continue;
                    }
                }
                masked.push(bytes[at] as char);
                at += 1;
            }
            masked.push('\n');
            masked
        })
        .collect()
}

// `LOSAT_TIMING`: each batch's block comes when the loop takes the batch, in batch order and in
// its place among the warnings, and no block of a discarded batch: standard error is that of
// the search without the switch, but for the seconds (and the end-of-run counters of other
// switches).
#[test]
fn timing_blocks_come_in_batch_order() {
    let query = queries_with_warnings();
    let subject = fixture("LC738868.fasta");
    let args = [
        "-query",
        path(&query),
        "-subject",
        subject.to_str().expect("UTF-8 path"),
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let env = [("BATCH_SIZE", "3000"), ("LOSAT_TIMING", "1")];
    let reference = blastn(&args, &env);
    assert!(
        reference.stderr.matches("[TIMING] total:").count() > 4,
        "{}",
        reference.stderr
    );
    assert!(
        reference.stderr.contains("FASTA-Reader"),
        "{}",
        reference.stderr
    );
    for switches in [
        vec![SWITCH, ("LOSAT_X_QBATCHPARSLOTS", "4")],
        vec![
            SWITCH,
            ("LOSAT_X_QBATCHPARSLOTS", "2"),
            ("LOSAT_X_QBATCHPARMISS", "2"),
        ],
    ] {
        let mut with = env.to_vec();
        with.extend_from_slice(&switches);
        let run = blastn(&args, &with);
        assert!(run.stdout == reference.stdout, "{switches:?}");
        assert_eq!(
            without_seconds(&run.stderr),
            without_seconds(&reference.stderr),
            "{switches:?}"
        );
        with.push(STATS);
        let run = blastn(&args, &with);
        assert!(run.stat("precomputed") > 0, "{}", run.stderr);
    }
}

// Without LOSAT_X_STATS the switch writes nothing of its own: standard error is the same bytes
// as without it, with warnings of later batches, at one thread (where the path is off) and four.
#[test]
fn the_switch_writes_nothing_of_its_own() {
    let query = queries_with_warnings();
    let subject = fixture("LC738868.fasta");
    for threads in ["1", "4"] {
        let args = [
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-outfmt",
            "7",
            "-num_threads",
            threads,
        ];
        let env = [("BATCH_SIZE", "1")];
        let reference = blastn(&args, &env);
        assert!(
            reference.stderr.contains("FASTA-Reader"),
            "{}",
            reference.stderr
        );
        assert!(
            reference.stderr.contains("q8 only N"),
            "{}",
            reference.stderr
        );
        let run = blastn(&args, &[env[0], SWITCH, SLOTS2]);
        assert_eq!(run.code, reference.code);
        assert!(run.stdout == reference.stdout, "{threads}");
        assert_eq!(run.stderr, reference.stderr, "{threads}");
    }
}

// Warnings of records in later batches (the reader's, an invalid query's) come where they come
// without the switch, with the batches taken from slots.
#[test]
fn warnings_of_later_batches_come_as_without_the_switch() {
    let query = queries_with_warnings();
    let subject = fixture("LC738868.fasta");
    for outfmt in ["0", "6"] {
        let args = [
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-outfmt",
            outfmt,
            "-num_threads",
            "4",
        ];
        let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SWITCH, SLOTS2]);
        assert!(run.stat("precomputed") > 0, "{}", run.stderr);
        let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SHADOW, SLOTS2]);
        assert_eq!(run.summary("[X_SHADOW] QBATCHPAR")["diffs"], "0");
    }
}

// `-task blastn`, one query per batch (pieces of LC738884), with the shadow. Batches of several
// queries are in the next test.
#[test]
fn blastn_task_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &pieces("LC738884.fasta"));
    let subject = fixture("LC738868.fasta");
    for outfmt in ["6", "7"] {
        let args = [
            "-task",
            "blastn",
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-outfmt",
            outfmt,
            "-num_threads",
            "4",
        ];
        let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SWITCH]);
        assert_eq!(run.stat("precomputed"), 31, "{}", run.stderr);
        let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SHADOW]);
        assert_eq!(run.summary("[X_SHADOW] QBATCHPAR")["compared"], "31");
        assert_eq!(run.summary("[X_SHADOW] QBATCHPAR")["diffs"], "0");
    }
}

// `-task blastn` with several queries per batch (pieces of LC738874 and of LC738884, outfmt 0, 6
// and 7). With the adaptive batch size the first batch is 5,000 letters (two queries) and the
// second holds the other 30: two batches, so the slots have nothing to search ahead (the path is
// off, `off=single`). With `BATCH_SIZE` 5,000 and 10,000 there are 14 and 8 batches of two to
// four queries, taken from the slots, and the shadow compares each with the loop's own search.
#[test]
fn blastn_task_with_several_queries_per_batch_gives_the_reports_without_the_switch() {
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    for name in ["LC738874.fasta", "LC738884.fasta"] {
        let query = TempFasta::new("q", &pieces(name));
        for outfmt in ["6", "7", "0"] {
            let args = [
                "-task",
                "blastn",
                "-query",
                path(&query),
                "-subject",
                subject,
                "-outfmt",
                outfmt,
                "-num_threads",
                "4",
            ];
            let run = same_as_without(&args, &[], &[SWITCH]);
            assert_eq!(run.stat("batches"), 2, "{name} {outfmt}: {}", run.stderr);
            assert_eq!(run.stat("precomputed"), 0, "{name} {outfmt}");
            for (size, batches) in [("5000", 14), ("10000", 8)] {
                let env = [("BATCH_SIZE", size)];
                let run = same_as_without(&args, &env, &[SWITCH]);
                assert_eq!(run.stat("batches"), batches, "{name} {outfmt} {size}");
                assert_eq!(
                    run.stat("precomputed"),
                    batches - 1,
                    "{name} {outfmt} {size}"
                );
                let run = same_as_without(&args, &env, &[SHADOW]);
                let shadow = run.summary("[X_SHADOW] QBATCHPAR");
                assert_eq!(
                    shadow["compared"],
                    (batches - 1).to_string(),
                    "{name} {outfmt} {size}: {}",
                    run.stderr
                );
                assert_eq!(shadow["diffs"], "0");
            }
        }
    }
}

// The slot model on whole runs (no `LOSAT_X_QBATCHPARSLOTS`). It takes slots for `-task blastn`
// with batches of 10,000 letters at four threads (eight batches, six units after the first:
// 2 x 2 predicted at 0.88 of the loop) and for megablast and dc-megablast with one query per
// batch at eight threads (31 units: 2 x 4 predicted at 0.68); it keeps the loop for megablast
// with batches of 20,000 letters at four threads (three units: the best, 3 x 1, predicted at
// 1.19) and for dc-megablast's mixer with one 50 kb query among the 2-4 kb ones
// (`CHUNK_SIZE=21000`: one unit after the first batch). The reports are those without the
// switch either way (and the shadow's, where slots run), and the summary line states the
// units, the layout judged and its prediction.
#[test]
fn the_rule_chooses_the_slots_on_whole_runs() {
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    let blastn_queries = TempFasta::new("q", &pieces("LC738884.fasta"));
    let plain = TempFasta::new("q", &queries());
    let mut records = queries();
    let long = genome("LC738874.fasta")[100_000..150_000].to_vec();
    records.insert(10, ("long 50 kb".to_string(), long));
    let with_long = TempFasta::new("q", &records);
    // (task, query, variables, threads, expected off, units, judged layout)
    let cases = [
        (
            "blastn",
            &blastn_queries,
            vec![("BATCH_SIZE", "10000")],
            "4",
            "none",
            6,
            Some((2, 2)),
        ),
        (
            "megablast",
            &plain,
            vec![("BATCH_SIZE", "1")],
            "8",
            "none",
            31,
            Some((2, 4)),
        ),
        (
            "dc-megablast",
            &plain,
            vec![("BATCH_SIZE", "1")],
            "8",
            "none",
            31,
            Some((2, 4)),
        ),
        (
            "megablast",
            &plain,
            vec![("BATCH_SIZE", "20000")],
            "4",
            "model",
            3,
            Some((3, 1)),
        ),
        (
            "dc-megablast",
            &with_long,
            vec![("CHUNK_SIZE", "21000")],
            "4",
            "single",
            1,
            None,
        ),
    ];
    for (task, query, env, threads, off, units, layout) in cases {
        let args = [
            "-task",
            task,
            "-query",
            path(query),
            "-subject",
            subject,
            "-outfmt",
            "6",
            "-num_threads",
            threads,
        ];
        let run = same_as_without(&args, &env, &[SWITCH]);
        let fields = run.summary("[X_STATS] QBATCHPAR");
        let what = format!("{task} {env:?} {threads} threads: {}", run.stderr);
        assert_eq!(fields["off"], off, "{what}");
        assert_eq!(run.stat("units"), units, "{what}");
        let judged = (run.stat("rule_k") as usize, run.stat("rule_m") as usize);
        match layout {
            None => {
                assert_eq!(judged, (0, 0), "{what}");
                assert_eq!(fields["predicted"], "-", "{what}");
            }
            Some(layout) => {
                assert_eq!(judged, layout, "{what}");
                let predicted: f64 = fields["predicted"].parse().expect("a ratio");
                assert_eq!(predicted < 0.97, off == "none", "{what}");
            }
        }
        if off == "none" {
            let (k, m) = layout.expect("a layout");
            assert_eq!(
                (run.stat("k") as usize, run.stat("m") as usize),
                (k, m),
                "{what}"
            );
            assert_eq!(run.stat("precomputed"), run.stat("batches") - 1, "{what}");
            // The shadow on the model's own layout: every batch from a slot compared.
            let run = same_as_without(&args, &env, &[SHADOW]);
            let shadow = run.summary("[X_SHADOW] QBATCHPAR");
            assert_eq!(shadow["diffs"], "0", "{what}");
            assert_eq!(
                shadow["compared"],
                (run.stat("batches") - 1).to_string(),
                "{what}"
            );
        } else {
            assert_eq!((run.stat("k"), run.stat("precomputed")), (0, 0), "{what}");
        }
    }
}

// A batch that NCBI splits into query chunks, searched in a slot: a 50,000-letter record in the
// middle of the input with `CHUNK_SIZE=21000` (a batch of more than 41,800 letters is split).
// With `BATCH_SIZE=1` the record is a batch of its own; with the mixer (dc-megablast, whose
// full batch is 20,000 letters with this `CHUNK_SIZE`) the batch that holds it is split too, and
// the planner predicts it without extensions.
#[test]
fn a_batch_split_into_query_chunks_in_a_slot() {
    let mut records = queries();
    let long = genome("LC738874.fasta")[100_000..150_000].to_vec();
    records.insert(10, ("long 50 kb".to_string(), long));
    let query = TempFasta::new("q", &records);
    let subject = fixture("LC738868.fasta");
    for (task, batch) in [("megablast", Some("1")), ("dc-megablast", None)] {
        let args = [
            "-task",
            task,
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ];
        let mut env = vec![("CHUNK_SIZE", "21000")];
        env.extend(batch.map(|size| ("BATCH_SIZE", size)));
        let run = same_as_without(&args, &env, &[SWITCH, SLOTS2]);
        assert!(run.stat("precomputed") > 0, "{task}: {}", run.stderr);
        let run = same_as_without(&args, &env, &[SHADOW, SLOTS2]);
        assert_eq!(run.summary("[X_SHADOW] QBATCHPAR")["diffs"], "0");
    }
}

// `-query_loc` that skips records in the middle of the input (they join the next batch), and a
// successful search into an `-out` file.
#[test]
fn skipped_records_and_an_out_file_give_the_reports_without_the_switch() {
    let mut records = queries();
    for index in [5, 6, 13] {
        records[index].1.truncate(600);
    }
    let query = TempFasta::new("q", &records);
    let subject = fixture("LC738868.fasta");
    let out = temp_path("out.txt");
    let out = out.to_str().expect("UTF-8 path").to_string();
    for outfmt in ["0", "6"] {
        let args = [
            "-query",
            path(&query),
            "-subject",
            subject.to_str().expect("UTF-8 path"),
            "-query_loc",
            "1001-1400",
            "-outfmt",
            outfmt,
            "-num_threads",
            "4",
            "-out",
            &out,
        ];
        let env = [("BATCH_SIZE", "1")];
        let reference = blastn(&args, &env);
        assert_eq!(reference.code, Some(0), "{}", reference.stderr);
        let expected = std::fs::read(&out).expect("read -out");
        assert!(!expected.is_empty());
        let run = blastn(&args, &[env[0], SWITCH, SLOTS2]);
        assert_eq!(run.code, reference.code);
        assert_eq!(run.stderr, reference.stderr);
        assert!(
            std::fs::read(&out).expect("read -out") == expected,
            "{outfmt}: -out"
        );
        let run = blastn(&args, &[env[0], SWITCH, SLOTS2, STATS]);
        assert!(run.stat("precomputed") > 0, "{}", run.stderr);
        assert!(
            std::fs::read(&out).expect("read -out") == expected,
            "{outfmt}: -out"
        );
    }
    let _ = std::fs::remove_file(&out);
}

// dc-megablast's `LOSAT_X_PAIRPARDISC` inside a slot's pool of two threads: batches of 10,000
// letters or more (the stage takes a query block over 8,000 letters), two slots of two threads
// at four threads. The stage's shadow compares the chunks of every batch, the loop's first
// and the slots'.
#[test]
fn pairpardisc_runs_inside_a_slots_pool() {
    let query = TempFasta::new("q", &queries());
    let subject = fixture("LC738868.fasta");
    let args = [
        "-task",
        "dc-megablast",
        "-query",
        path(&query),
        "-subject",
        subject.to_str().expect("UTF-8 path"),
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let env = [("BATCH_SIZE", "10000")];
    let switches = [
        SWITCH,
        ("LOSAT_X_QBATCHPARSLOTS", "2"),
        ("LOSAT_X_PAIRPARDISC", "1"),
        ("LOSAT_X_PAIRPARDISCSHADOW", "1"),
    ];
    let reference = blastn(&args, &env);
    let mut with = env.to_vec();
    with.extend_from_slice(&switches);
    with.push(STATS);
    let run = blastn(&args, &with);
    assert_eq!(run.code, reference.code, "{}", run.stderr);
    assert!(run.stdout == reference.stdout);
    let batches = run.stat("batches");
    assert_eq!((run.stat("k"), run.stat("m")), (2, 2), "{}", run.stderr);
    assert!(run.stat("precomputed") > 0, "{}", run.stderr);
    let line = run
        .stderr
        .lines()
        .find(|line| line.starts_with("[X_PAIRPARDISCSHADOW] blastn: seed chunks equal="))
        .unwrap_or_else(|| panic!("no PAIRPARDISC shadow line:\n{}", run.stderr));
    let chunks: u64 = line
        .trim_start_matches("[X_PAIRPARDISCSHADOW] blastn: seed chunks equal=")
        .split_whitespace()
        .next()
        .and_then(|value| value.parse().ok())
        .expect("a chunk count");
    // The loop searches the first batch; every other chunk was a slot's (the last batch may be
    // below the stage's 8,000 letters).
    assert!(
        chunks > 1 && chunks + 1 >= batches,
        "{line} for {batches} batches"
    );
}

// Where the loop keeps every batch: one thread, as many subject work units as threads, and the
// thread diagnostics (written while a batch is searched).
#[test]
fn the_loop_keeps_the_batches_where_the_path_is_off() {
    let query = TempFasta::new("q", &queries());
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    let base = ["-query", path(&query), "-outfmt", "6", "-subject", subject];
    // (threads, variables, why the path is off)
    let cases = [
        ("1", None, "threads"),
        ("4", Some(("LOSAT_WASI_THREADS_DEBUG", "1")), "debug"),
        (
            "4",
            Some(("LOSAT_TRACE_BLASTN_SUBJECT", "no_such_subject")),
            "debug",
        ),
    ];
    for (threads, variable, off) in cases {
        let mut args = base.to_vec();
        args.extend_from_slice(&["-num_threads", threads]);
        let mut env: Vec<(&str, &str)> = variable.into_iter().collect();
        env.push(("BATCH_SIZE", "1"));
        let run = same_as_without(&args, &env, &[SWITCH]);
        assert_eq!(run.stat("precomputed"), 0);
        assert_eq!(
            run.summary("[X_STATS] QBATCHPAR")["off"],
            off,
            "{}",
            run.stderr
        );
    }
    let subjects: Vec<(String, Vec<u8>)> = (0..4)
        .map(|index| {
            (
                format!("s{index}"),
                genome("LC738868.fasta")[index * 50_000..][..40_000].to_vec(),
            )
        })
        .collect();
    let subjects = TempFasta::new("s4", &subjects);
    let args = [
        "-query",
        path(&query),
        "-outfmt",
        "6",
        "-subject",
        path(&subjects),
        "-num_threads",
        "4",
    ];
    let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SWITCH]);
    assert_eq!(run.summary("[X_STATS] QBATCHPAR")["off"], "subjects");
    // With LOSAT_X_BNSTORE the records keep no letters (the packed subjects hold them): the work
    // units are counted from the stored letters, so the subjects still keep the loop (before
    // 2026-10-17 the records' empty letters counted 0 units here).
    let stored = [("BATCH_SIZE", "1"), ("LOSAT_X_BNSTORE", "1")];
    let run = same_as_without(&args, &stored, &[SWITCH]);
    assert_eq!(
        run.summary("[X_STATS] QBATCHPAR")["off"],
        "subjects",
        "{}",
        run.stderr
    );
    assert_eq!(
        run.summary("[X_STATS] BNSTORE")["subjects"],
        "4",
        "{}",
        run.stderr
    );
}

// Debug output (written while a batch is searched) and a BLAST database (`-db`, made here by
// `LOSAT makedb`) keep the loop: the reports are those without the switch. (`-verbose` is not a
// command-line option of LOSAT's blastn.)
#[test]
fn debug_output_and_a_database_keep_the_loop() {
    let records: Vec<(String, Vec<u8>)> = queries().into_iter().take(4).collect();
    let query = TempFasta::new("q", &records);
    let subject = fixture("LC738868.fasta");
    let subject = subject.to_str().expect("UTF-8 path");
    for variable in [
        ("LOSAT_DEBUG_BLASTN", "1"),
        ("LOSAT_DEBUG_COORDS", "1"),
        ("BLEMIR_DEBUG", "1"),
    ] {
        let args = [
            "-query",
            path(&query),
            "-subject",
            subject,
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ];
        let env = [variable, ("BATCH_SIZE", "1")];
        let reference = blastn(&args, &env);
        let run = blastn(&args, &[env[0], env[1], SWITCH, STATS]);
        assert_eq!(run.code, reference.code, "{variable:?}");
        assert!(run.stdout == reference.stdout, "{variable:?}");
        assert_eq!(
            run.summary("[X_STATS] QBATCHPAR")["off"],
            "debug",
            "{variable:?}"
        );
    }
    let db = temp_path("db");
    std::fs::create_dir_all(&db).expect("make the database directory");
    let name = db.join("subject");
    let name = name.to_str().expect("UTF-8 path");
    let made = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .args(["makedb", "-in", subject, "-dbtype", "nucl", "-out", name])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        made.status.success(),
        "{}",
        String::from_utf8_lossy(&made.stderr)
    );
    let args = [
        "-query",
        path(&query),
        "-db",
        name,
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let run = same_as_without(&args, &[("BATCH_SIZE", "1")], &[SWITCH]);
    assert_eq!(run.summary("[X_STATS] QBATCHPAR")["off"], "db");
    let _ = std::fs::remove_dir_all(&db);
}
