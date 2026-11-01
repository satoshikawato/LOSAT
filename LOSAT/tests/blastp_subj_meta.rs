//! blastp's block search (`LOSAT_X_BPQBLOCK`) with the subjects' tabular IDs and titles made once
//! for all blocks (`src/algorithm/blastp/blast_engine/x_subject_meta.rs`) on whole searches.
//!
//! The stored values belong to the `-subject` block search (no switch of their own since
//! 2026-10-31). These tests compare the report, standard error and exit code of that block search
//! with the search of all queries (`LOSAT_X_BPQBLOCK=0`); the reach rows of the block run
//! (`LOSAT_X_REACH`, `tests/reach_support`) show that it took blocks. A subject whose ID cannot be
//! written is refused before any block runs, as in the search of all queries.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
//! ```c++
//!         } catch (const CException&) {
//!             list<CRef<CSeq_id> > subject_ids;
//!             CRef<CSeq_id> id(new CSeq_id());
//!             id->Assign(align.GetSeq_id(1));
//!             subject_ids.push_back(id);
//!             SetSubjectId(subject_ids);
//!             bioseqs_found = false;
//!         }
//! ```

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), each
/// `0`, the default before, unless the search sets it: every search here runs as it was written,
/// with the switches off that it does not name (the default runs in `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of one test (test-unique: tests of this file run side by side).
fn test_dir(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpsubjmeta-{test}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// What one run printed, and its exit code.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

/// One blastp run of the fixture queries (`MeenMJNV.faa`) against `subject` (`-subject FILE` or
/// `-db NAME`) with `-task task`, with `vars` as its only `LOSAT_*` variables that choose a code
/// path, `BATCH_SIZE` 3000 (several query batches, so several blocks) and `-outfmt outfmt`.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing; without the caller's value a file of the test's directory is
/// named.
fn blastp(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    outfmt: &str,
) -> Run {
    blastp_with(dir, task, subject, threads, vars, outfmt, &[])
}

/// `blastp` with the arguments `extra` after the others (e.g. `-max_target_seqs 3`).
fn blastp_with(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    vars: &[(&str, &str)],
    outfmt: &str,
    extra: &[&str],
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
        command.env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
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
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// `base` with `extra` added.
fn with<'a>(base: &[(&'a str, &'a str)], extra: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
    base.iter().chain(extra).copied().collect()
}

/// The block search (`base` names `LOSAT_X_BPQBLOCK`) and the search of all queries
/// (`LOSAT_X_BPQBLOCK=0`) give the same report, standard error and exit code; returns the run of
/// the search of all queries and the reach rows of the block run.
fn same_runs(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    outfmt: &str,
) -> (Run, reach_support::Reach) {
    same_runs_with(dir, task, subject, threads, base, outfmt, &[])
}

/// `same_runs` with the arguments `extra` in every run.
fn same_runs_with(
    dir: &Path,
    task: &str,
    subject: &[&OsStr],
    threads: u32,
    base: &[(&str, &str)],
    outfmt: &str,
    extra: &[&str],
) -> (Run, reach_support::Reach) {
    let what = format!("{task} {base:?} threads {threads} outfmt {outfmt} {extra:?}");
    let all = blastp_with(
        dir,
        task,
        subject,
        threads,
        &with(base, &[("LOSAT_X_BPQBLOCK", "0")]),
        outfmt,
        extra,
    );
    let reach = reach_support::Reach::new("subj_meta");
    let blocks = blastp_with(
        dir,
        task,
        subject,
        threads,
        &with(base, &[reach.var()]),
        outfmt,
        extra,
    );
    assert!(blocks.stdout == all.stdout, "{what}: the report differs");
    assert_eq!(blocks.stderr, all.stderr, "{what}: standard error differs");
    assert_eq!(blocks.code, all.code, "{what}: the exit code differs");
    (all, reach)
}

// The block search in blocks of one query batch and of about 20,000 letters: every block's output
// stage takes the stored IDs and titles, and the rows and standard error are those of the search
// of all queries. The IDs and titles are made once per block search (`XBlockSubjects::new`) and read at the
// start of each block's search and by its output, the same code at any thread count and in both
// tasks (the task chooses only the lookup table, which neither reads nor writes these values):
// blastp at 4 threads in both block sizes, blastp-fast at 1 thread in blocks of 20,000.
#[test]
fn block_search_reports_are_the_reports_of_the_search_of_all_queries() {
    let dir = test_dir("blocks");
    let subject = fixture("MejoMJNV.faa");
    let subject = [OsStr::new("-subject"), subject.as_os_str()];
    for (task, threads, sizes) in [
        ("blastp", 4, &["1", "20000"][..]),
        ("blastp-fast", 1, &["20000"][..]),
    ] {
        for &letters in sizes {
            let base = [("LOSAT_X_BPQBLOCK", letters)];
            let (all, reach) = same_runs(&dir, task, &subject, threads, &base, "6");
            let what = format!("{task} threads {threads} blocks of {letters}");
            assert_eq!(all.code, Some(0), "{what}: {}", all.stderr);
            assert!(!all.stdout.is_empty(), "{what}: the fixture pair has rows");
            assert!(
                reach.sides("LOSAT_X_BPQBLOCK").0,
                "{what}: the block search did not run"
            );
        }
    }
    // With the parallel subject loop in runs (`LOSAT_X_BPSUBJRUN`), as the workstation runs it.
    let base = [("LOSAT_X_BPQBLOCK", "1"), ("LOSAT_X_BPSUBJRUN", "1")];
    let (_, reach) = same_runs(&dir, "blastp-fast", &subject, 4, &base, "6");
    assert!(reach.sides("LOSAT_X_BPQBLOCK").0);
    let _ = std::fs::remove_dir_all(&dir);
}

/// The fixture's subject proteome with the title of the subject that the first row of the report
/// names replaced by `Subject_\x81 ...` (a `Subject_` first word whose byte 0x81 is in no encoding
/// NCBI guesses: no tabular ID can be written for it), in `dir`; returns its path.
fn subject_file_with_an_undecodable_reported_title(dir: &Path, first_row: &[u8]) -> PathBuf {
    let sseqid = first_row.split(|&byte| byte == b'\t').nth(1).unwrap();
    let text = std::fs::read(fixture("MejoMJNV.faa")).unwrap();
    let mut out = Vec::new();
    let mut replaced = 0;
    for line in text.split_inclusive(|&byte| byte == b'\n') {
        let header = line.strip_prefix(b">");
        if header.is_some_and(|header| header.starts_with(sseqid) && header[sseqid.len()] == b' ') {
            out.extend_from_slice(b">Subject_\x81 undecodable title\n");
            replaced += 1;
        } else {
            out.extend_from_slice(line);
        }
    }
    assert_eq!(replaced, 1, "{}", String::from_utf8_lossy(sseqid));
    let path = dir.join("undecodable.faa");
    std::fs::write(&path, out).unwrap();
    path
}

// A reported subject whose tabular ID cannot be written: the run is refused with the same message,
// exit code and rows (none: LOSAT rejects it before any format is written) in the block search and
// in the search of all queries. The block search is not taken for such a subject set
// (`x_query_block_plan`), so no block runs and no stored value is made. The check of the stored
// values with such records is tested in `x_subject_meta.rs`. One search setting: the plan refuses
// such a subject set whatever the block size, before any block runs, at any thread count.
#[test]
fn a_reported_subject_without_a_tabular_id_is_refused_alike_in_blocks() {
    let dir = test_dir("undecodable");
    let subject = fixture("MejoMJNV.faa");
    let reference = blastp(
        &dir,
        "blastp",
        &[OsStr::new("-subject"), subject.as_os_str()],
        1,
        &[],
        "6",
    );
    assert_eq!(reference.code, Some(0), "{}", reference.stderr);
    let first_row = reference
        .stdout
        .split(|&byte| byte == b'\n')
        .next()
        .unwrap();
    let path = subject_file_with_an_undecodable_reported_title(&dir, first_row);
    let subject = [OsStr::new("-subject"), path.as_os_str()];
    let base = [("LOSAT_X_BPQBLOCK", "1")];
    let (all, _) = same_runs(&dir, "blastp", &subject, 4, &base, "6");
    assert_ne!(all.code, Some(0));
    assert!(
        all.stderr
            .contains("has a title whose first word starts with 'Subject_'"),
        "{}",
        all.stderr
    );
    assert!(all.stdout.is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

/// The fixture's subject proteome with the titles of the subjects that the report `reference`
/// names first replaced by `titles` (a title is the header's text after `>`, so a title without a
/// space replaces the ID too), in `dir`; returns its path.
fn subject_file_with_retitled_reported_subjects(
    dir: &Path,
    reference: &[u8],
    titles: &[&[u8]],
) -> PathBuf {
    let mut reported: Vec<&[u8]> = Vec::new();
    for row in reference.split(|&byte| byte == b'\n') {
        if let Some(sseqid) = row.split(|&byte| byte == b'\t').nth(1) {
            if !reported.contains(&sseqid) {
                reported.push(sseqid);
            }
        }
    }
    assert!(
        reported.len() >= titles.len(),
        "{} reported",
        reported.len()
    );
    let text = std::fs::read(fixture("MejoMJNV.faa")).unwrap();
    let mut out = Vec::new();
    let mut replaced = 0;
    for line in text.split_inclusive(|&byte| byte == b'\n') {
        let first_word = line.strip_prefix(b">").map(|header| {
            let end = header
                .iter()
                .position(|&byte| byte == b' ' || byte == b'\n')
                .unwrap_or(header.len());
            &header[..end]
        });
        match first_word
            .and_then(|word| reported.iter().take(titles.len()).position(|r| *r == word))
        {
            Some(index) => {
                out.extend_from_slice(b">");
                out.extend_from_slice(titles[index]);
                out.push(b'\n');
                replaced += 1;
            }
            None => out.extend_from_slice(line),
        }
    }
    assert_eq!(replaced, titles.len());
    let path = dir.join("retitled.faa");
    std::fs::write(&path, out).unwrap();
    path
}

// Subjects whose tabular ID is not the title's first word: an empty title (`unnamed`), a
// whitespace-only title, `Subject_` first words (NCBI generates their definition line: `Subject_99`,
// a cp1252 / latin-1 byte, `Subject_ `, `Subject_3`) and `lcl|Subject_7`. With the stored IDs the
// block search writes the same rows as the search of all queries, in the columns that show the ID, accession and title, with
// `-max_target_seqs 3`, in blastp-fast (the IDs are made the same way in both tasks; its report,
// as blastp's, shows all seven retitled subjects).
#[test]
fn generated_subject_ids_are_written_alike_in_blocks() {
    let dir = test_dir("generated-ids");
    let original = fixture("MejoMJNV.faa");
    let reference = blastp(
        &dir,
        "blastp",
        &[OsStr::new("-subject"), original.as_os_str()],
        1,
        &[],
        "6",
    );
    assert_eq!(reference.code, Some(0), "{}", reference.stderr);
    let titles: [&[u8]; 7] = [
        b"",
        b"Subject_99 some desc",
        b"Subject_\xe9 latin1 desc",
        b"lcl|Subject_7 foo",
        b"   ",
        b"Subject_ x",
        b"Subject_3",
    ];
    let path = subject_file_with_retitled_reported_subjects(&dir, &reference.stdout, &titles);
    let subject = [OsStr::new("-subject"), path.as_os_str()];
    let outfmt = "6 qseqid sseqid sacc saccver stitle slen evalue bitscore";
    let base = [("LOSAT_X_BPQBLOCK", "1")];
    let task = "blastp-fast";
    let (all, reach) = same_runs_with(
        &dir,
        task,
        &subject,
        4,
        &base,
        outfmt,
        &["-max_target_seqs", "3"],
    );
    assert_eq!(all.code, Some(0), "{task}: {}", all.stderr);
    assert!(
        reach.sides("LOSAT_X_BPQBLOCK").0,
        "{task}: the block search did not run"
    );
    // The generated IDs are in the report, so the stored IDs were written.
    for generated in [&b"unnamed"[..], b"Subject_99"] {
        assert!(
            all.stdout.windows(generated.len()).any(|w| w == generated),
            "{task}: {} is not in the report",
            String::from_utf8_lossy(generated)
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}
