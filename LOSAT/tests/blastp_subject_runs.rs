//! blastp's parallel `-subject` loop in runs of subjects with `LOSAT_X_BPSUBJRUN`
//! (`src/algorithm/blastp/blast_engine/x_subject_runs.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on by default, `LOSAT_X_BPSUBJRUN=0`
//! turns it off. It acts only where the
//! subjects of a `-subject` search are searched on the pool (more than one thread and more than
//! one subject). With it on, the subjects are cut into runs of
//! consecutive subjects, each run searched with one diagonal table that is filled once, at the
//! one-thread order's offset of the run's first subject. These tests compare the reports with the
//! switch off (`LOSAT_X_BPSUBJRUN=0`) and on (the default; blastp and blastp-fast, one hit and two
//! hits, the block search of `LOSAT_X_BPQBLOCK`); the reach rows of the run with the switch on
//! (`LOSAT_X_REACH`, `tests/reach_support`) show that the path ran (and, at one thread and on
//! `-db`, that it did not). Each kernel and case runs at one thread count: the runs are cut by the
//! subjects' lengths and the thread count only, and the unit tests of `x_subject_runs.rs` cut them
//! at several thread counts and check each later subject's offset and table.
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

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The subjects of the fixture's subject proteome (`MejoMJNV.faa`).
const SUBJECTS: u64 = 124;

/// What one run printed: the report.
struct Run {
    stdout: Vec<u8>,
}

/// One blastp run of the fixture pair (two related nimavirus proteomes) with `-task task`,
/// against `subject` (`-subject FILE` or `-db NAME`), with `vars` as its only `LOSAT_*` variables
/// that choose a code path, `BATCH_SIZE` 3000 (several query batches), `extra` arguments and
/// `-outfmt outfmt`. `LOSAT_X_BPQBLOCK` is 0 (one search of all queries) unless `vars` sets it, so
/// that the searches do not depend on the strict set's default. `LOSAT_BLASTX_WORKER_LOG` stays: a
/// binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI)
/// needs that file name as soon as a search has more than one thread, and it selects nothing;
/// without the caller's value a file of the tests is named. The run must exit with 0.
fn blastp(
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
    }
}

fn subject_file() -> PathBuf {
    fixture("MejoMJNV.faa")
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

const OFF: (&str, &str) = ("LOSAT_X_BPSUBJRUN", "0");

/// The switch off and on (the default: the switch is in the strict set) give the same report;
/// returns the reach rows of the run with the switch on (`LOSAT_X_BPSUBJRUN`).
fn same_reports(
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    extra: &[&str],
    outfmt: &str,
) -> reach_support::Reach {
    let what = format!("{task} {base:?} {extra:?} threads {threads} outfmt {outfmt}");
    let off = blastp(task, subject, threads, &with(base, &[OFF]), extra, outfmt);
    assert!(!off.stdout.is_empty(), "{what}: the fixture pair has rows");
    let reach = reach_support::Reach::new("subject_runs");
    let on = blastp(
        task,
        subject,
        threads,
        &with(base, &[reach.var()]),
        extra,
        outfmt,
    );
    assert!(
        on.stdout == off.stdout,
        "{what}: the report differs with the switch"
    );
    reach
}

/// The runs of the switch were searched (its new side, not the one-thread loop's).
fn assert_engaged(reach: &reach_support::Reach, what: &str) {
    let sides = reach.sides("LOSAT_X_BPSUBJRUN");
    assert!(sides.0, "{what}: the runs were not searched: {sides:?}");
}

// blastp and blastp-fast, the reference loops (the two-hit loop of
// `blastp_for_each_offset_pair` and the reference compressed scan; `LOSAT_X_BPSEED=0`, and
// `LOSAT_X_BPCLUT=0` for blastp-fast) and the fused loop of `LOSAT_X_BPSEED` (the default; with
// the flat table of `LOSAT_X_BPCLUT` for blastp-fast): one search of all queries, run by the
// switch, each kernel at 2 or 4 threads.
#[test]
fn search_of_all_queries_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (task, kernel, base, threads) in [
        ("blastp", "reference", &[("LOSAT_X_BPSEED", "0")][..], 2),
        ("blastp", "fused", &[][..], 4),
        (
            "blastp-fast",
            "reference",
            &[("LOSAT_X_BPSEED", "0"), ("LOSAT_X_BPCLUT", "0")][..],
            4,
        ),
        ("blastp-fast", "fused", &[][..], 2),
    ] {
        let reach = same_reports(task, &subject, threads, base, &[], "6");
        assert_engaged(&reach, &format!("{task} {kernel} threads {threads}"));
    }
}

// One hit (`-window_size 0`, the one-hit loop of `blastp_for_each_offset_pair`; `LOSAT_X_BPSEED`
// does not act at window 0): an entry left by an earlier subject of a run reads as a cleared one
// there too (`diff >= 0`).
#[test]
fn one_hit_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let reach = same_reports("blastp", &subject, 3, &[], &["-window_size", "0"], "6");
    assert_engaged(&reach, "one hit");
}

// The block search (`LOSAT_X_BPQBLOCK`, every block a search of all subjects). A report that is
// not written in blocks is one search of all queries, as in the tests above.
#[test]
fn block_search_reports_with_the_switch_are_the_reports_without_it() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let base = [("LOSAT_X_BPQBLOCK", "1")];
    let reach = same_reports("blastp", &subject, 4, &base, &[], "6");
    assert!(
        reach.sides("LOSAT_X_BPQBLOCK").0,
        "the block search did not run"
    );
    assert_engaged(&reach, "block search");
}

// `-window_size 10000000`: each subject moves the diagonal offset by its length and 10,000,000, so
// the offset passes INT4_MAX / 4 twice in the fixture's subjects (after the 54th and the 108th
// subject; at 2 threads both inside runs). `finish_subject` (`Blast_ExtendWordExit`) puts it back
// to the window and clears the table inside the run, and the later subjects of the run start from
// there as in the one-thread loop. Each run allocates a diagonal table of 64 MiB.
#[test]
fn the_int4_max_quarter_reset_inside_a_run_gives_the_reports_without_the_switch() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let reach = same_reports(
        "blastp",
        &subject,
        2,
        &[],
        &["-window_size", "10000000"],
        "6",
    );
    assert_engaged(&reach, "window 10000000");
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
    let reach = same_reports("blastp", &subject, 2, &[], &[], "6");
    assert_engaged(&reach, "empty record");
    let _ = std::fs::remove_file(&path);
}

/// A protein database of the subject proteome, made with `losat makedb` in a directory of the
/// test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpsubjrun-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(subject_file())
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

// Where the subjects are not searched on the pool (one thread; a database, whose blocks are
// searched with one table each already), the switch does not act: the run with the switch on
// takes no new side (`LOSAT_X_REACH`). The one-thread `-subject` loop marks the reference side
// (the run went through the switch's decision); the `-db` search marks none.
#[test]
fn the_switch_does_not_act_at_one_thread_or_on_a_database() {
    let subject = subject_file();
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    let database = subject_database();
    let db = [OsStr::new("-db"), database.as_os_str()];
    for (subject, threads) in [(&subject, 1), (&db, 4)] {
        let reach = reach_support::Reach::new("subject_runs_none");
        let run = blastp("blastp", subject, threads, &[reach.var()], &[], "6");
        assert!(
            !run.stdout.is_empty(),
            "threads {threads}: the fixture pair has rows"
        );
        let sides = reach.sides("LOSAT_X_BPSUBJRUN");
        let want = if threads == 1 {
            (false, true)
        } else {
            (false, false)
        };
        assert_eq!(sides, want, "threads {threads}: the runs were searched?");
    }
    let _ = std::fs::remove_dir_all(database.parent().unwrap());
}
