//! TBLASTN's `-db` search with each database sequence read where
//! it is searched (`tblastn/x_db_stream.rs`) on a whole run of the binary.
//!
//! The stream graduated on 2026-10-27 (the switch was `LOSAT_X_TBNDBSTREAM`, with a shadow mode
//! that compared every read with a copy of the database; the copy and the shadow are deleted, and
//! the frozen NCBI `-db` fixtures hold the reports). The properties under test are LOSAT's own: the
//! reports of the search at four threads are those of the search at one thread, with
//! composition-based statistics (the redo reads its subjects) and without (`-comp_based_stats 0`:
//! the traceback reads them), in one query batch and in several (`BATCH_SIZE`: the plans kept
//! between batches are not read again), on a database of several volumes with ambiguity codes and
//! a subject long enough to be split into units across the threads; the stages that let a subject
//! go after its last match (`LOSAT_X_TBNSTAGEDROP`) give the reports of those that hold it to their
//! end; and a database that cannot be read fails with the same error.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path, then `env`.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env_remove("BATCH_SIZE");
    command.envs(env.iter().copied());
    command
}

/// A temporary directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-tblastn-db-stream-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The database of `fasta_file` made by `losat makeblastdb` in `dir` under `name`, in volumes of at
/// most `max_file_sz` bytes.
fn makedb(dir: &Path, fasta_file: &Path, name: &str, max_file_sz: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat(&[])
        .arg("makeblastdb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "nucl", "-out", &db, "-max_file_sz", max_file_sz])
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

/// The subjects: the slices of `tslice` (ambiguity codes, an N block) and the MelaMJNV genome
/// (287 kb, a subject of its own volume, split into units at four threads); the queries: the
/// seven of `tslice_q.faa` and the first 20 proteins of PemoMJNVA, so that several queries align
/// with the same subjects.
fn inputs(dir: &Path) -> (String, PathBuf) {
    let mut subjects = std::fs::read(fasta("db_search/tslice.fna")).unwrap();
    subjects.extend_from_slice(&std::fs::read(fasta("fasta/MelaMJNV.fasta")).unwrap());
    let subjects_file = dir.join("subjects.fna");
    std::fs::write(&subjects_file, subjects).unwrap();
    let db = makedb(dir, &subjects_file, "subjects", "20000");
    let mut queries = std::fs::read_to_string(fasta("db_search/tslice_q.faa")).unwrap();
    let proteins = std::fs::read_to_string(fasta("fasta/PemoMJNVA.faa")).unwrap();
    for record in proteins.split_inclusive('>').skip(1).take(20) {
        queries.push('>');
        queries.push_str(record.trim_end_matches('>'));
    }
    let query = dir.join("queries.faa");
    std::fs::write(&query, queries).unwrap();
    (db, query)
}

/// One search of `query` on `db` with `args` at `threads` threads; returns the exit status, the
/// report and the standard error.
fn tblastn(
    query: &Path,
    db: &str,
    args: &[&str],
    threads: u32,
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let output = losat(env)
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .args(["-db", db])
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

#[test]
fn the_reports_do_not_depend_on_the_threads_or_the_query_batches() {
    let dir = temp_dir("reports");
    let (db, query) = inputs(&dir);
    for args in [
        &["-outfmt", "6"][..],
        &["-outfmt", "0"][..],
        &["-outfmt", "6", "-comp_based_stats", "0"][..],
    ] {
        let reference = tblastn(&query, &db, args, 1, &[]);
        assert_eq!(reference.0, Some(0), "{args:?}: {}", reference.2);
        assert!(
            reference.1.split(|&byte| byte == b'\n').count() > 20,
            "{args:?}: too few lines"
        );
        for batch in [None, Some("2000")] {
            let batch_env: Vec<(&str, &str)> =
                batch.map(|b| ("BATCH_SIZE", b)).into_iter().collect();
            for threads in [1, 4] {
                let run = tblastn(&query, &db, args, threads, &batch_env);
                assert!(
                    run == reference,
                    "{args:?} batch {batch:?} {threads} threads"
                );
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `LOSAT_X_TBNSTAGEDROP`: the redo and the traceback of `-comp_based_stats 0` let a subject go
/// after its last match; the reports are those of the stages that hold every subject until they
/// end, at one and four threads. The 27 queries align with the same subjects, so most subjects
/// have several matches.
#[test]
fn a_stage_lets_a_subject_go_after_its_last_match_with_the_same_reports() {
    let dir = temp_dir("stagedrop");
    let (db, query) = inputs(&dir);
    for args in [
        &["-outfmt", "6"][..],
        &["-outfmt", "6", "-comp_based_stats", "0"][..],
    ] {
        for threads in [1, 4] {
            let held = tblastn(&query, &db, args, threads, &[("LOSAT_X_TBNSTAGEDROP", "0")]);
            assert_eq!(held.0, Some(0), "{}", held.2);
            let dropped = tblastn(&query, &db, args, threads, &[]);
            assert!(dropped == held, "{args:?} {threads} threads");
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `LOSAT_X_TBNPLANCUT`: a thread's block ends with its chunk of subjects; the reports are those
/// of the block of 256 KiB kept to the next read, in one query batch and in several, with
/// composition-based statistics and without.
#[test]
fn a_block_that_ends_with_its_chunk_gives_the_same_reports() {
    let dir = temp_dir("plan-cut");
    let (db, query) = inputs(&dir);
    for threads in [2, 4, 8] {
        for args in [
            &["-outfmt", "6"][..],
            &["-outfmt", "6", "-comp_based_stats", "0"],
        ] {
            for batch in [&[][..], &[("BATCH_SIZE", "2000")]] {
                let run = |off: &[(&str, &str)]| {
                    let env: Vec<_> = batch.iter().chain(off).copied().collect();
                    tblastn(&query, &db, args, threads, &env)
                };
                let (cut, kept) = (run(&[]), run(&[("LOSAT_X_TBNPLANCUT", "0")]));
                assert_eq!(cut.0, Some(0), "{}", cut.2);
                assert!(!cut.1.is_empty());
                assert_eq!(cut.1, kept.1, "threads {threads} {args:?} {batch:?}");
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `LOSAT_X_TBNLOOPFREE`: a thread frees the buffers it kept between subjects when its loop over
/// the chunks ends, and makes new ones for the next query batch; the reports are those of the
/// buffers kept to the next batch, with composition-based statistics and without: one query batch,
/// batches searched one after another on one pool (2 threads with 3 batches; 4 and 8 threads with
/// 6 batches and no pools of their own), and batches on pools of their own (8 threads).
#[test]
fn buffers_freed_at_the_end_of_a_threads_loop_give_the_same_reports() {
    let dir = temp_dir("loop-free");
    let (db, query) = inputs(&dir);
    let own_pools_off = [("BATCH_SIZE", "2000"), ("LOSAT_X_TBNOWNPOOLS", "0")];
    let cases: [(u32, &[(&str, &str)]); 7] = [
        (2, &[]),
        (4, &[]),
        (8, &[]),
        (2, &[("BATCH_SIZE", "5000")]),
        (4, &own_pools_off),
        (8, &own_pools_off),
        (8, &[("BATCH_SIZE", "2000")]),
    ];
    for (threads, batch) in cases {
        for args in [
            &["-outfmt", "6"][..],
            &["-outfmt", "6", "-comp_based_stats", "0"],
        ] {
            let run = |off: &[(&str, &str)]| {
                let env: Vec<_> = batch.iter().chain(off).copied().collect();
                tblastn(&query, &db, args, threads, &env)
            };
            let (freed, kept) = (run(&[]), run(&[("LOSAT_X_TBNLOOPFREE", "0")]));
            assert_eq!(freed.0, Some(0), "{}", freed.2);
            assert!(!freed.1.is_empty());
            assert_eq!(freed.1, kept.1, "threads {threads} {args:?} {batch:?}");
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one (as
// `tblastn_db_search.rs` `the_sequence_file_is_read_by_the_first_batch_with_a_valid_context`).
#[test]
fn an_unreadable_database_fails_the_same_batch_with_the_same_error() {
    let dir = temp_dir("unreadable");
    let (db, query) = inputs(&dir);
    // The last volume's sequence file cut short: the index names bytes the file does not hold
    // (LOSAT's refusal, at the first batch with a valid context).
    let last = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.to_string_lossy().ends_with(".nsq"))
        .max()
        .unwrap();
    let bytes = std::fs::read(&last).unwrap();
    std::fs::write(&last, &bytes[..bytes.len() / 2]).unwrap();
    let all_x = format!(">xx all X\n{}\n", "X".repeat(32));
    let mut both = all_x.into_bytes();
    both.extend_from_slice(&std::fs::read(&query).unwrap());
    let query = dir.join("xx.faa");
    std::fs::write(&query, both).unwrap();
    let env = [("BATCH_SIZE", "1")];
    let one = tblastn(&query, &db, &[], 1, &env);
    assert_ne!(one.0, Some(0), "{}", one.2);
    assert!(one.2.contains("passes the end of the file"), "{}", one.2);
    for threads in [2, 4] {
        assert!(
            tblastn(&query, &db, &[], threads, &env) == one,
            "{threads} threads"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// A sequence longer than a block of the check is a block of its own, which the check reads without
// its 2-bit data (`common::database::x_check_nucleotide_stores`, `SeqDb::x_nucleotide_check`): only
// the last 2-bit byte and the ambiguity data.
#[test]
fn a_damaged_sequence_of_its_own_block_fails_or_reports_the_same_at_every_thread_count() {
    let dir = temp_dir("damaged");
    let (_, query) = inputs(&dir);
    // The MelaMJNV genome with an ambiguity code every 1000 bases (runs of one, the old format),
    // the last sequence and alone in its volume.
    let mut letters: Vec<u8> = std::fs::read_to_string(fasta("fasta/MelaMJNV.fasta"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect();
    let runs = (500..letters.len()).step_by(1000).count();
    for at in (500..letters.len()).step_by(1000) {
        letters[at] = b"NRYKMSWBDHV"[at / 1000 % 11];
    }
    let mut subjects = std::fs::read(fasta("db_search/tslice.fna")).unwrap();
    subjects.extend_from_slice(b">mela_amb MelaMJNV with ambiguity codes\n");
    for line in letters.chunks(80) {
        subjects.extend_from_slice(line);
        subjects.push(b'\n');
    }
    let subjects_file = dir.join("subjects_amb.fna");
    std::fs::write(&subjects_file, subjects).unwrap();
    let db = makedb(&dir, &subjects_file, "subjects_amb", "20000");
    let last = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            let name = path.to_string_lossy();
            name.contains("subjects_amb.") && name.ends_with(".nsq")
        })
        .max()
        .unwrap();
    let original = std::fs::read(&last).unwrap();
    // The sentinel byte, the 2-bit data, then the ambiguity data: a count word and a word per run.
    let ambiguity = 1 + letters.len() / 4 + 1;
    assert_eq!(original.len(), ambiguity + 4 + 4 * runs, "one sequence");
    let count =
        |bytes: &[u8]| u32::from_be_bytes(bytes[ambiguity..ambiguity + 4].try_into().unwrap());
    assert_eq!(count(&original) as usize, runs);
    let mut long_count = original.clone();
    long_count[ambiguity..ambiguity + 4].copy_from_slice(&(count(&original) + 100).to_be_bytes());
    let mut remainder = original.clone();
    remainder[ambiguity - 1] ^= 1;
    let all_x = format!(">xx all X\n{}\n", "X".repeat(32));
    let mut both = all_x.into_bytes();
    both.extend_from_slice(&std::fs::read(&query).unwrap());
    let query = dir.join("xx.faa");
    std::fs::write(&query, both).unwrap();
    let env = [("BATCH_SIZE", "1")];
    for (what, bytes, error) in [
        (
            "count past the end",
            long_count,
            Some("the ambiguity data of database sequence 5 cannot be read"),
        ),
        ("remainder", remainder, None),
        (
            "cut in the ambiguity data",
            original[..ambiguity + 8].to_vec(),
            Some("passes the end of the file"),
        ),
        (
            "cut in the 2-bit data",
            original[..ambiguity / 2].to_vec(),
            Some("passes the end of the file"),
        ),
    ] {
        std::fs::write(&last, bytes).unwrap();
        let copy = tblastn(&query, &db, &[], 1, &env);
        match error {
            Some(error) => {
                assert_ne!(copy.0, Some(0), "{what}: {}", copy.2);
                assert!(copy.2.contains(error), "{what}: {}", copy.2);
            }
            None => {
                assert_eq!(copy.0, Some(0), "{what}: {}", copy.2);
                assert!(!copy.1.is_empty(), "{what}");
            }
        }
        for threads in [2, 4] {
            assert!(
                tblastn(&query, &db, &[], threads, &env) == copy,
                "{what}: {threads} threads"
            );
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
