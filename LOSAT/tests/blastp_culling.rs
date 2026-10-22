//! BLASTP `-culling_limit` (`src/algorithm/blastp/blast_engine/culling.rs`) on whole searches:
//! the option is read as NCBI reads it (non-negative, hexadecimal too, 0 is no culling), the
//! search with it culls (the preliminary writer and the pipe of `hspfilter_culling.c`), and the
//! report does not depend on the thread count, on the subject source (`-subject` or `-db`) or on
//! the representation of the preliminary lists (`LOSAT_X_DEFAULT=0`). That the bytes equal NCBI's
//! is the business of the parity probe (`O3b`); the unit tests of `culling.rs` pin the key, the
//! length of the tree, the shape of the final hit list and the order of the pipe.
//!
//! NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c (the 300 subjects of
//! `e2e_many_subject.faa` hold related proteins; NCBI 2.17.0 keeps 300 / 1 / 2 / 5 rows of
//! them with no limit / 1 / 2 / 5, outfmt 6).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta/outfmt0")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastp_culling_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The binary with `env` as its only variables that choose a code path.
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: String,
    stderr: String,
    code: Option<i32>,
}

fn blastp(args: &[&str], env: &[(&str, &str)]) -> Run {
    let output = losat(env)
        .arg("blastp")
        .args(args)
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// The subject database of the fixture, removed when dropped.
struct Db {
    dir: PathBuf,
    name: PathBuf,
}

impl Db {
    fn new() -> Self {
        let dir = temp_path("db");
        std::fs::create_dir_all(&dir).unwrap();
        let name = dir.join("many");
        let made = losat(&[])
            .arg("makedb")
            .arg("-in")
            .arg(fixture("e2e_many_subject.faa"))
            .args(["-dbtype", "prot", "-out"])
            .arg(&name)
            .output()
            .expect("run LOSAT makedb");
        assert!(
            made.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&made.stderr)
        );
        Db { dir, name }
    }
}

impl Drop for Db {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn subject_search(extra: &[&str], env: &[(&str, &str)]) -> Run {
    let query = fixture("e2e_many_query.faa");
    let subject = fixture("e2e_many_subject.faa");
    let mut args = vec![
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    args.extend_from_slice(extra);
    blastp(&args, env)
}

fn db_search(db: &Db, extra: &[&str], env: &[(&str, &str)]) -> Run {
    let query = fixture("e2e_many_query.faa");
    let mut args = vec![
        "-query",
        query.to_str().unwrap(),
        "-db",
        db.name.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    args.extend_from_slice(extra);
    blastp(&args, env)
}

fn rows(run: &Run) -> usize {
    assert_eq!(run.code, Some(0), "{}", run.stderr);
    run.stdout.lines().count()
}

// NCBI c++/src/algo/blast/blastinput/blast_args.cpp:3292-3333: `-culling_limit` is an integer
// with the constraint `>= 0`; the parser reads hexadecimal values (blastn's fixture
// `culling.hex_0x2`) and a negative or non-numeric value is an argument error.
#[test]
fn the_option_is_read_as_ncbi_reads_it() {
    for (value, expected) in [("0", 300), ("1", 1), ("2", 2), ("0x2", 2)] {
        let run = subject_search(&["-culling_limit", value], &[]);
        assert_eq!(rows(&run), expected, "{value}");
    }
    for bad in ["-1", "1.5", "abc", "2147483648"] {
        let run = subject_search(&["-culling_limit", bad], &[]);
        assert_eq!(run.code, Some(2), "{bad}: {}", run.stderr);
        assert!(
            !run.stderr.contains("is not supported by"),
            "{bad}: {}",
            run.stderr
        );
    }
}

// NCBI `SetCullingLimit(0)` is ignored: the report is the one without the option.
#[test]
fn a_limit_of_zero_is_no_culling() {
    let none = subject_search(&[], &[]);
    let zero = subject_search(&["-culling_limit", "0"], &[]);
    assert_eq!(none.stdout, zero.stdout);
    assert_eq!(rows(&none), 300);
}

// hspfilter_culling.c: with a limit of 1 or 2 the related subjects are dominated by the best
// HSPs of the query's tree, in the preliminary stage and again after the traceback.
#[test]
fn the_limit_culls_the_related_subjects() {
    assert_eq!(rows(&subject_search(&["-culling_limit", "1"], &[])), 1);
    assert_eq!(rows(&subject_search(&["-culling_limit", "2"], &[])), 2);
    assert_eq!(rows(&subject_search(&["-culling_limit", "5"], &[])), 5);
}

// The report does not depend on the thread count, on the representation of the preliminary
// lists, or on the subject source.
#[test]
fn the_report_is_the_same_at_every_thread_count_and_source() {
    let db = Db::new();
    for limit in ["1", "2", "5"] {
        let reference = subject_search(&["-culling_limit", limit], &[]);
        assert_eq!(rows(&reference), limit.parse::<usize>().unwrap());
        for threads in ["1", "4"] {
            let args = ["-culling_limit", limit, "-num_threads", threads];
            for env in [
                &[][..],
                &[("LOSAT_X_DEFAULT", "0")][..],
                &[("LOSAT_X_BPPCOMPACT", "shadow")][..],
            ] {
                let subject = subject_search(&args, env);
                assert_eq!(subject.stdout, reference.stdout, "{limit} -subject {threads} {env:?}");
                let database = db_search(&db, &args, env);
                assert_eq!(database.stdout, reference.stdout, "{limit} -db {threads} {env:?}");
            }
        }
    }
}

// A query split into chunks has a writer for each chunk (CHUNK_SIZE as NCBI's environment
// variable): the query of the fixture is shorter than the chunk, so every chunk size runs the
// same search; this pins that a split search with the option runs and reports.
#[test]
fn a_split_query_search_runs_with_the_option() {
    let split = subject_search(
        &["-culling_limit", "1"],
        &[("CHUNK_SIZE", "100"), ("OVERLAP_CHUNK_SIZE", "20")],
    );
    assert!(rows(&split) >= 1);
}
