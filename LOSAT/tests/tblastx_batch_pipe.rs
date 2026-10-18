//! The batch pipe of TBLASTX (`LOSAT_X_TBXPIPE`, `src/algorithm/tblastx/x_batch_pipe.rs`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_TBXPIPE` is
//! `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it on and two or more threads, the next query
//! batch starts when the batch before has ended its seed stage, and the batch loop takes the
//! results in its order. `[X_STATS] TBXPIPE` counts the planned and taken batches and how each was
//! started; `LOSAT_X_TBXPIPESHADOW=1` searches every taken batch again on the loop and compares the
//! two results; `LOSAT_X_TBXPIPEHSPS=0` admits no early start.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-195
//! ```c
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!                 results = lcl_blast.Run();
//! ```
//! The property under test is that the reports, the messages and the exit status do not depend on
//! the switch or on the thread count: the searches are compared with the same LOSAT searches with
//! the switch off. Each query record below is longer than TBLASTX's batch size (10,002 letters),
//! so each is a batch of its own.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const PIPE_OFF: (&str, &str) = ("LOSAT_X_TBXPIPE", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBXPIPESHADOW", "1");
const NO_EARLY: (&str, &str) = ("LOSAT_X_TBXPIPEHSPS", "0");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 5] = [1, 2, 4, 8, 16];

fn letters(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    std::fs::read(path)
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(String, &[u8])]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-tbxpipe-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let path = dir.join(file);
    let mut text = Vec::new();
    for (title, seq) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).expect("write FASTA");
    path
}

/// Four query batches of 10,500 letters each from MelaMJNV, then the records `tail`.
fn queries(file: &str, tail: &[(String, &[u8])]) -> PathBuf {
    let mela = letters("MelaMJNV.fasta");
    let mut records: Vec<(String, &[u8])> = (0..4)
        .map(|i| (format!("mela_{i}"), &mela[i * 60_000..i * 60_000 + 10_500]))
        .collect();
    records.extend(tail.iter().cloned());
    write_fasta(file, &records)
}

/// Two subjects of 45,000 letters each from PemoMJNVA, written once (the tests run in parallel).
fn subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let pemo = letters("PemoMJNVA.fasta");
        write_fasta(
            "subjects.fa",
            &[
                ("pemo_a".to_string(), &pemo[..45_000]),
                ("pemo_b".to_string(), &pemo[150_000..195_000]),
            ],
        )
    })
}

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path
/// (`LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread, and it selects nothing). `LOSAT_X_PAIRPARGATE=0`:
/// the batches are short, and the gate would send them to the reference frame loop; the pipe acts
/// on both, and the unit path is the one long batches take.
fn tblastx(query: &Path, threads: u32, outfmt: &str, switches: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env("LOSAT_X_PAIRPARGATE", "0");
    command.envs(switches.iter().copied());
    command
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects())
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The counters `names` of the `[X_STATS] TBXPIPE` line.
fn counters<const N: usize>(output: &Output, names: [&str; N]) -> [u64; N] {
    let text = stderr(output);
    let found = text
        .lines()
        .find(|line| line.starts_with("[X_STATS] TBXPIPE "))
        .unwrap_or_else(|| panic!("no TBXPIPE line:\n{text}"));
    names.map(|name| {
        let key = format!("{name}=");
        found
            .split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {found}"))
            .parse()
            .expect("a count")
    })
}

#[test]
fn the_report_does_not_depend_on_the_pipe_at_any_thread_count() {
    let query = queries("four.fa", &[]);
    let off = tblastx(&query, 1, "6", &[PIPE_OFF]);
    assert!(off.status.success(), "{}", stderr(&off));
    assert!(!off.stdout.is_empty(), "the pairs have hits");
    for threads in THREADS {
        let on = tblastx(&query, threads, "6", &[STATS]);
        assert!(on.status.success(), "threads {threads}: {}", stderr(&on));
        assert_eq!(on.stdout, off.stdout, "threads {threads}");
        let shadow = tblastx(&query, threads, "6", &[SHADOW, STATS]);
        assert!(
            shadow.status.success(),
            "threads {threads}: {}",
            stderr(&shadow)
        );
        assert_eq!(shadow.stdout, off.stdout, "threads {threads}, shadow");
        if threads == 1 {
            // One thread: the loop searches every batch itself, as before.
            assert!(!stderr(&on).contains("TBXPIPE"), "{}", stderr(&on));
            continue;
        }
        let [planned, taken, early, loop_searched, key_mismatch] =
            counters(&on, ["planned", "taken", "early", "loop", "key_mismatch"]);
        assert_eq!(
            (planned, taken, key_mismatch),
            (4, 4, 0),
            "threads {threads}"
        );
        assert_eq!(early + loop_searched, 4, "threads {threads}");
        assert!(early >= 1, "threads {threads}: no batch started early");
        assert!(
            stderr(&shadow).contains("[X_SHADOW] TBXPIPE batches=4 compared=4 diffs=0"),
            "threads {threads}: {}",
            stderr(&shadow)
        );
    }
}

#[test]
fn the_pairwise_report_does_not_depend_on_the_pipe() {
    let query = queries("four_pairwise.fa", &[]);
    let off = tblastx(&query, 4, "0", &[PIPE_OFF]);
    let on = tblastx(&query, 4, "0", &[]);
    assert!(
        off.status.success() && on.status.success(),
        "{}",
        stderr(&on)
    );
    assert_eq!(on.stdout, off.stdout);
}

#[test]
fn without_room_for_held_hsps_the_loop_starts_every_batch() {
    let query = queries("four_no_early.fa", &[]);
    let off = tblastx(&query, 4, "6", &[PIPE_OFF]);
    let held = tblastx(&query, 4, "6", &[NO_EARLY, STATS]);
    assert!(held.status.success(), "{}", stderr(&held));
    assert_eq!(held.stdout, off.stdout);
    assert_eq!(
        counters(&held, ["planned", "taken", "early", "loop"]),
        [4, 4, 0, 4]
    );
}

#[test]
fn the_previous_default_keeps_the_loop() {
    let query = queries("four_previous.fa", &[]);
    let off = tblastx(&query, 4, "6", &[PIPE_OFF]);
    let previous = tblastx(&query, 4, "6", &[PREVIOUS, STATS]);
    assert!(previous.status.success(), "{}", stderr(&previous));
    assert_eq!(previous.stdout, off.stdout);
    assert!(
        !stderr(&previous).contains("TBXPIPE"),
        "{}",
        stderr(&previous)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652 (a query without
// letters gets a message; a batch whose queries all lack a valid context then stops the search
// after the reports of the batches before, `BLAST engine error:`, exit 3).
#[test]
fn a_batch_that_fails_its_checks_stops_the_search_where_the_loop_stopped_it() {
    let query = queries("four_then_empty.fa", &[("empty".to_string(), &[][..])]);
    let off = tblastx(&query, 8, "6", &[PIPE_OFF]);
    let on = tblastx(&query, 8, "6", &[STATS]);
    assert_eq!(on.status.code(), off.status.code());
    assert_eq!(off.status.code(), Some(3), "{}", stderr(&off));
    assert_eq!(on.stdout, off.stdout);
    let without_counters: String = stderr(&on)
        .lines()
        .filter(|line| !line.starts_with("[X_STATS]"))
        .map(|line| format!("{line}\n"))
        .collect();
    assert_eq!(without_counters, stderr(&off));
    let [planned, taken] = counters(&on, ["planned", "taken"]);
    assert_eq!((planned, taken), (4, 4));
}
