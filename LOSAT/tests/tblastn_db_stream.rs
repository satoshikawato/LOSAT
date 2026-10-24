//! EXPERIMENT (LOSAT_X_TBNDBSTREAM): TBLASTN's `-db` search with each database sequence read where
//! it is searched (`tblastn/x_db_stream.rs`) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the reports are those of the search on the copy of
//! the whole database (`LOSAT_X_TBNDBSTREAM=0`) at one and four threads, with composition-based
//! statistics (the redo reads its subjects) and without (`-comp_based_stats 0`: the traceback reads
//! them), in one query batch and in several (`BATCH_SIZE`: the plans kept between batches are not
//! read again), on a database of several volumes with ambiguity codes and a subject long enough to
//! be split into units across the threads; the shadow mode (`LOSAT_X_TBNDBSTREAMSHADOW=1`) finds
//! every read identical to the copy, with no more reads at four threads than at one; and a
//! database that cannot be read fails with the same error.

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

/// The database of `fasta_file` made by `losat makedb` in `dir` under `name`, in volumes of at
/// most `max_file_sz` bytes.
fn makedb(dir: &Path, fasta_file: &Path, name: &str, max_file_sz: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat(&[])
        .arg("makedb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "nucl", "-out", &db, "-max_file_sz", max_file_sz])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
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
fn the_stream_gives_the_reports_of_the_copy() {
    let dir = temp_dir("reports");
    let (db, query) = inputs(&dir);
    for args in [
        &["-outfmt", "6"][..],
        &["-outfmt", "0"][..],
        &["-outfmt", "6", "-comp_based_stats", "0"][..],
    ] {
        for batch in [None, Some("2000")] {
            let batch_env: Vec<(&str, &str)> =
                batch.map(|b| ("BATCH_SIZE", b)).into_iter().collect();
            let mut copy_env = batch_env.clone();
            copy_env.push(("LOSAT_X_TBNDBSTREAM", "0"));
            for threads in [1, 4] {
                let copy = tblastn(&query, &db, args, threads, &copy_env);
                assert_eq!(copy.0, Some(0), "{args:?} {batch:?} {threads}: {}", copy.2);
                assert!(
                    copy.1.split(|&byte| byte == b'\n').count() > 20,
                    "{args:?}: too few lines"
                );
                let stream = tblastn(&query, &db, args, threads, &batch_env);
                assert!(stream == copy, "{args:?} batch {batch:?} {threads} threads");
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_shadow_mode_finds_every_read_identical() {
    let dir = temp_dir("shadow");
    let (db, query) = inputs(&dir);
    for (args, batch) in [
        (&["-outfmt", "6"][..], None),
        (
            &["-outfmt", "6", "-comp_based_stats", "0"][..],
            Some("2000"),
        ),
    ] {
        let batch_env: Vec<(&str, &str)> = batch.map(|b| ("BATCH_SIZE", b)).into_iter().collect();
        let mut copy_env = batch_env.clone();
        copy_env.push(("LOSAT_X_TBNDBSTREAM", "0"));
        let mut shadow_env = batch_env.clone();
        shadow_env.push(("LOSAT_X_TBNDBSTREAMSHADOW", "1"));
        let mut reads_at = Vec::new();
        for threads in [1, 4] {
            let copy = tblastn(&query, &db, args, threads, &copy_env);
            let shadow = tblastn(&query, &db, args, threads, &shadow_env);
            assert_eq!(shadow.0, Some(0), "{}", shadow.2);
            assert_eq!(shadow.1, copy.1, "{args:?} {threads} threads");
            let reads =
                shadow_reads(&shadow.2).unwrap_or_else(|| panic!("no shadow line: {}", shadow.2));
            assert!(reads > 0, "{}", shadow.2);
            reads_at.push(reads);
            // Without the shadow variable, nothing is printed.
            assert!(!tblastn(&query, &db, args, threads, &batch_env)
                .2
                .contains("TBNDBSTREAM"));
        }
        // Threads that ask for a sequence at once share one read (`XDbStream::read_shared`).
        assert!(
            reads_at[1] <= reads_at[0],
            "{args:?}: reads at 1 / 4 threads {reads_at:?}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// The reads of `[X_SHADOW] TBNDBSTREAM sequences=6 reads=N bases=B identical`, the last line of
/// `stderr`.
fn shadow_reads(stderr: &str) -> Option<u64> {
    let line = stderr.lines().last()?;
    let rest = line.strip_prefix("[X_SHADOW] TBNDBSTREAM sequences=6 reads=")?;
    let (reads, rest) = rest.split_once(" bases=")?;
    let (bases, rest) = rest.split_once(' ')?;
    bases.parse::<u64>().ok()?;
    (rest == "identical").then_some(())?;
    reads.parse().ok()
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
    let copy = tblastn(&query, &db, &[], 1, &[env[0], ("LOSAT_X_TBNDBSTREAM", "0")]);
    assert_ne!(copy.0, Some(0), "{}", copy.2);
    assert!(copy.2.contains("passes the end of the file"), "{}", copy.2);
    for threads in [1, 4] {
        assert!(
            tblastn(&query, &db, &[], threads, &env) == copy,
            "{threads} threads"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
