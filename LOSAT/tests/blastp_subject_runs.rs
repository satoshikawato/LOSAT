//! blastp's parallel `-subject` loop in runs of subjects with `LOSAT_X_BPSUBJRUN`
//! (`src/algorithm/blastp/blast_engine/x_subject_runs.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, `LOSAT_X_BPSUBJRUN=0`
//! turns it off, and `LOSAT_X_DEFAULT=0` leaves it off unless it is set. It acts only where the
//! subjects of a `-subject` search are searched on the pool (more than one thread and more than
//! one subject). With it on, the subjects are cut into runs of
//! consecutive subjects, each run searched with one diagonal table that is filled once, at the
//! one-thread order's offset of the run's first subject. These tests compare the reports with the
//! switch off (`LOSAT_X_BPSUBJRUN=0`) and on (blastp and blastp-fast, one hit and two hits, the block search of
//! `LOSAT_X_BPQBLOCK`, `-db`) and run its shadow mode (`LOSAT_X_BPSUBJRUNSHADOW`), which checks
//! every later subject of a run's offset and table and compares every subject's result with the
//! result from a cleared table; its line at exit shows that the path ran. A difference stops the
//! run.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-172
//! ```c
//!         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//!             ewp->diag_table->offset = ewp->diag_table->window;
//!             s_BlastDiagClear(ewp->diag_table);
//!         } else {
//!             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//!         }
//! ```

use std::collections::HashMap;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The subjects of the fixture's subject proteome (`MejoMJNV.faa`).
const SUBJECTS: u64 = 124;

/// What one run printed: the report and standard error.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
}

impl Run {
    /// The counters of the `[X_SHADOW] BPSUBJRUN` line, which must end with `identical`.
    fn shadow(&self) -> HashMap<String, u64> {
        let prefix = "[X_SHADOW] BPSUBJRUN ";
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("no {prefix}line in {}", self.stderr));
        assert!(line.ends_with(" identical"), "{line}");
        line[prefix.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
    }
}

/// One blastp run of the fixture pair (two related nimavirus proteomes) with `-task task`,
/// against `subject` (`-subject FILE` or `-db NAME`), with `vars` as its only `LOSAT_*` variables
/// that choose a code path, `BATCH_SIZE` 3000 (several query batches), `extra` arguments and
/// `-outfmt outfmt`. `LOSAT_X_BPQBLOCK` is 0 (one search of all queries) unless `vars` sets it, so
/// that the counts do not depend on the strict set's default. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing; without the caller's value
/// a file of the tests is named. The run must exit with 0.
fn blastp(
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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
            std::env::temp_dir().join("losat-bpsubjrun-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .args(["-task", task])
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
        .args(subject)
        .args(["-evalue", "0.001", "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .env("BATCH_SIZE", "3000")
        .env("LOSAT_X_BPQBLOCK", "0")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {vars:?} {extra:?} threads {threads}: {stderr}"
    );
    Run {
        stdout: output.stdout,
        stderr,
    }
}

fn subject_file() -> PathBuf {
    fixture("MejoMJNV.faa")
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

const ON: (&str, &str) = ("LOSAT_X_BPSUBJRUN", "1");
const OFF: (&str, &str) = ("LOSAT_X_BPSUBJRUN", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BPSUBJRUNSHADOW", "1");

/// The switch off, on and in its shadow mode give the same report; returns the shadow's counters.
fn same_reports(
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
) -> HashMap<String, u64> {
    let what = format!("{task} {base:?} {extra:?} threads {threads} outfmt {outfmt}");
    let off = blastp(task, subject, threads, &with(base, &[OFF]), extra, outfmt);
    assert!(!off.stdout.is_empty(), "{what}: the fixture pair has rows");
    let on = blastp(task, subject, threads, &with(base, &[ON]), extra, outfmt);
    assert!(
        on.stdout == off.stdout,
        "{what}: the report differs with the switch"
    );
    let shadow = blastp(
        task,
        subject,
        threads,
        &with(base, &[SHADOW]),
        extra,
        outfmt,
    );
    assert!(
        shadow.stdout == off.stdout,
        "{what}: the report differs in the shadow mode"
    );
    shadow.shadow()
}

/// The counters of a shadow run where the switch acted in `searches` searches of `subjects`
/// subjects each (at most one run per subject, at least one per search): every subject compared,
/// every subject but the first of each run searched without a fill and with its table checked.
fn assert_counts(counts: &HashMap<String, u64>, searches: u64, subjects: u64, what: &str) {
    assert_eq!(counts["searches"], searches, "{what}: {counts:?}");
    assert_eq!(
        counts["subjects"],
        subjects * searches,
        "{what}: {counts:?}"
    );
    assert_eq!(
        counts["subjects_compared"], counts["subjects"],
        "{what}: {counts:?}"
    );
    assert!(
        (searches..counts["subjects"]).contains(&counts["runs"]),
        "{what}: {counts:?}"
    );
    assert_eq!(
        counts["fills_saved"],
        counts["subjects"] - counts["runs"],
        "{what}: {counts:?}"
    );
    assert_eq!(
        counts["tables_checked"], counts["fills_saved"],
        "{what}: {counts:?}"
    );
}

/// `assert_counts` for the fixture's subjects, with no INT4_MAX / 4 reset (with the default window
/// the fixture's offsets stay far below it).
fn assert_engaged(counts: &HashMap<String, u64>, searches: u64, what: &str) {
    assert_counts(counts, searches, SUBJECTS, what);
    assert_eq!(counts["resets_in_runs"], 0, "{what}: {counts:?}");
}

// blastp and blastp-fast at 2 and 4 threads, the reference loops (the two-hit loop of
// `blastp_for_each_offset_pair` and the reference compressed scan; `LOSAT_X_DEFAULT=0`, so that they
// stay the reference loops where the strict set turns `LOSAT_X_BPSEED` on) and the fused loop of
// `LOSAT_X_BPSEED` (with the flat table of `LOSAT_X_BPCLUT` for blastp-fast): one search of all
// queries, run by the switch.
#[test]
fn search_of_all_queries_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (task, base) in [
        ("blastp", &[("LOSAT_X_DEFAULT", "0")][..]),
        ("blastp", &[("LOSAT_X_BPSEED", "1")][..]),
        ("blastp-fast", &[("LOSAT_X_DEFAULT", "0")][..]),
        (
            "blastp-fast",
            &[("LOSAT_X_BPSEED", "1"), ("LOSAT_X_BPCLUT", "1")][..],
        ),
    ] {
        for threads in [2, 4] {
            let counts = same_reports(task, &subject, threads, base, &[], "6");
            assert_engaged(&counts, 1, &format!("{task} {base:?} threads {threads}"));
        }
    }
}

// The switch is on by default: the report with no variable, with `LOSAT_X_BPSUBJRUN=0`, with
// `LOSAT_X_DEFAULT=0` (off unless set) and with `LOSAT_X_DEFAULT=0` and the switch set are one
// report, at 4 threads for blastp and blastp-fast.
#[test]
fn the_default_the_switch_off_and_the_global_off_give_one_report() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for task in ["blastp", "blastp-fast"] {
        let default = blastp(task, &subject, 4, &[], &[], "6");
        assert!(
            !default.stdout.is_empty(),
            "{task}: the fixture pair has rows"
        );
        for vars in [
            &[OFF][..],
            &[("LOSAT_X_DEFAULT", "0")][..],
            &[("LOSAT_X_DEFAULT", "0"), ON][..],
        ] {
            let run = blastp(task, &subject, 4, vars, &[], "6");
            assert!(
                run.stdout == default.stdout,
                "{task} {vars:?}: the report differs from the default's"
            );
        }
    }
}

// One hit (`-window_size 0`): an entry left by an earlier subject of a run reads as a cleared one
// there too (`diff >= 0`).
#[test]
fn one_hit_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for base in [&[][..], &[("LOSAT_X_BPSEED", "1")][..]] {
        let counts = same_reports("blastp", &subject, 3, base, &["-window_size", "0"], "6");
        assert_engaged(&counts, 1, &format!("one hit {base:?}"));
    }
}

// The block search (`LOSAT_X_BPQBLOCK`, every block a search of all subjects) and the pairwise
// report, which is not written in blocks.
#[test]
fn block_search_and_pairwise_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [("LOSAT_X_BPQBLOCK", "1"), ("LOSAT_X_BPSEED", "1")];
    let counts = same_reports("blastp", &subject, 4, &base, &[], "6");
    assert!(counts["searches"] > 1, "{counts:?}");
    assert_engaged(&counts, counts["searches"], "block search");
    let counts = same_reports("blastp-fast", &subject, 2, &base, &[], "0");
    assert_engaged(&counts, 1, "pairwise report");
}

// `-window_size 10000000`: each subject moves the diagonal offset by its length and 10,000,000, so
// the offset passes INT4_MAX / 4 twice in the fixture's subjects (after the 54th and the 108th
// subject; at 2 threads both inside runs). `finish_subject` (`Blast_ExtendWordExit`) puts it back
// to the window and clears the table inside the run, and the later subjects of the run start from
// there as in the one-thread loop. Each run allocates diagonal tables of 64 MiB (two per thread in
// the shadow mode).
#[test]
fn the_int4_max_quarter_reset_inside_a_run_gives_the_reports_without_the_switch() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let counts = same_reports(
        "blastp",
        &subject,
        2,
        &[],
        &["-window_size", "10000000"],
        "6",
    );
    assert_counts(&counts, 1, SUBJECTS, "window 10000000");
    assert!(counts["resets_in_runs"] > 0, "{counts:?}");
}

/// The fixture's subject proteome with a record without letters before its 61st record, in a file
/// of the test; returns its path.
fn subject_file_with_an_empty_record() -> PathBuf {
    let text = std::fs::read_to_string(subject_file()).unwrap();
    let mut out = String::new();
    let mut records = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            if records == 60 {
                out.push_str(">empty_record\n");
            }
            records += 1;
        }
        out.push_str(line);
        out.push('\n');
    }
    assert_eq!(records, SUBJECTS as usize);
    let path =
        std::env::temp_dir().join(format!("losat-bpsubjrun-empty-{}.faa", std::process::id()));
    std::fs::write(&path, out).unwrap();
    path
}

// A subject without letters (NCBI's warning, not searched, the diagonal offset not moved) inside a
// run: at 2 threads the record at index 60 falls in the run of indices 55..67 (the cut depends only
// on the lengths and the thread count).
#[test]
fn a_subject_without_letters_inside_a_run_gives_the_reports_without_the_switch() {
    let path = subject_file_with_an_empty_record();
    let subject = [OsStr::new("-subject"), path.as_os_str()];
    let counts = same_reports("blastp", &subject, 2, &[], &[], "6");
    assert_counts(&counts, 1, SUBJECTS + 1, "empty record");
    assert_eq!(counts["resets_in_runs"], 0, "{counts:?}");
    let _ = std::fs::remove_file(&path);
}

/// A protein database of the subject proteome, made with `LOSAT makedb` in a directory of the
/// test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpsubjrun-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("makedb")
        .arg("-in")
        .arg(subject_file())
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

// Where the subjects are not searched on the pool (one thread; a database, whose blocks are
// searched with one table each already), the switch changes nothing and its shadow compares
// nothing.
#[test]
fn the_switch_does_not_act_at_one_thread_or_on_a_database() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let counts = same_reports("blastp", &subject, 1, &[], &[], "6");
    assert!(counts.values().all(|&count| count == 0), "{counts:?}");
    let database = subject_database();
    let db = [OsStr::new("-db"), database.as_os_str()];
    let counts = same_reports("blastp", &db, 4, &[], &[], "6");
    assert!(counts.values().all(|&count| count == 0), "{counts:?}");
    let _ = std::fs::remove_dir_all(database.parent().unwrap());
}
