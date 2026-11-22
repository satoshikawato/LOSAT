//! `-import_search_strategy FILE` of blastn, blastp, blastx, tblastn and tblastx (line RG, item
//! F2-2) against the bytes that NCBI BLAST+ 2.17.0 `-num_threads 1` writes: the standard output,
//! the standard error and the exit status of each row of
//! `tests/fixtures/strategy/import/import_cases.tsv` (its header says what the columns mean; the
//! strategy files are NCBI's own exports, `ss/`, and each output is stored once, `bytes/`). The rows are the RG-F2 probes of these programs
//! whose LOSAT result equals NCBI's: imports with `-db` and `-subject` strategies, the command
//! line's overrides (`-query`, `-db`, `-subject`, the formatting options), the warnings of the
//! options a strategy ignores, the strategies' query ranges and filters, the archive (outfmt 11)
//! of an imported run, the reader's errors, another program's strategy with options of the
//! command line (`C.*.from.*`), the parameters that the setters take (`S1.*`: best hits, database
//! masks, the repeat filter, the MegaBLAST index), the command line's index and database mask
//! over a strategy's (`B.*`), and the formatting options that a strategy's hit list size above
//! 500 leaves unread in the pairwise formats (`Y.*`).
//!
//! A strategy's subjects are searched once whatever `-num_threads` says (Owner 2026-11-20 F2-Q1
//! "Q1=C"; NCBI prints every alignment once per thread or aborts): the runs at 2, 4 and 8 threads
//! write the frozen `-num_threads 1` bytes.

use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixtures() -> PathBuf {
    engine().join("tests/fixtures/strategy/import")
}

struct Case {
    name: String,
    program: String,
    args: Vec<String>,
    rc: i32,
    out: String,
    err: String,
}

fn cases() -> Vec<Case> {
    let table =
        std::fs::read_to_string(fixtures().join("import_cases.tsv")).expect("import_cases.tsv");
    table
        .lines()
        .filter(|line| !(line.starts_with('#') || line.starts_with("name\t") || line.is_empty()))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 6, "import_cases.tsv row: {line}");
            Case {
                name: cells[0].to_string(),
                program: cells[1].to_string(),
                args: cells[2]
                    .split(' ')
                    .map(|arg| arg.replace("%20", " "))
                    .collect(),
                rc: cells[3].parse().expect("rc"),
                out: cells[4].to_string(),
                err: cells[5].to_string(),
            }
        })
        .collect()
}

/// NCBI's bytes of the file `name` of `bytes/` (`-`: none) in LOSAT's words.
fn frozen(name: &str) -> Vec<u8> {
    if name == "-" {
        return Vec::new();
    }
    crate::outfmt11_blastx::expected_losat_bytes(
        std::fs::read(fixtures().join("bytes").join(name)).expect("frozen bytes"),
    )
}

/// One run of the case from `LOSAT/` without LOSAT's switches and batch variables; standard
/// input is empty.
fn run(case: &Case, threads: Option<&str>) -> Output {
    let ss = fixtures().join("ss");
    let ss = ss.to_str().expect("UTF-8 path");
    let mut command = Command::new(BIN);
    command
        .current_dir(engine())
        .arg(&case.program)
        .args(case.args.iter().map(|arg| arg.replace("{ss}", ss)))
        .stdin(Stdio::null());
    if let Some(threads) = threads {
        command.args(["-num_threads", threads]);
    }
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE" | "BLASTDB"
            )
        {
            command.env_remove(name);
        }
    }
    command.output().expect("run losat")
}

fn check(case: &Case, threads: Option<&str>) {
    let output = run(case, threads);
    let label = format!("{} threads {threads:?}", case.name);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(
        output.status.code(),
        Some(case.rc),
        "{label}: exit status; stderr: {stderr}"
    );
    let want_out = frozen(&case.out);
    assert!(
        output.stdout == want_out,
        "{label}: stdout differs ({} bytes, NCBI wrote {}): {}",
        output.stdout.len(),
        want_out.len(),
        crate::outfmt11_blastx::first_difference(&output.stdout, &want_out)
    );
    let stderr = stderr.replace(fixtures().join("ss").to_str().unwrap(), "{SS}");
    assert_eq!(
        stderr,
        String::from_utf8_lossy(&frozen(&case.err)),
        "{label}: stderr"
    );
}

#[test]
fn imports_equal_ncbis() {
    let cases = cases();
    assert!(cases.len() >= 100, "import_cases.tsv rows");
    let failed: Vec<&str> = cases
        .iter()
        .filter(|case| std::panic::catch_unwind(|| check(case, None)).is_err())
        .map(|case| case.name.as_str())
        .collect();
    assert!(failed.is_empty(), "rows that differ: {failed:?}");
}

/// F2-Q1 (Owner 2026-11-20 "Q1=C"): with the strategy's subjects, `-num_threads` 2, 4 and 8
/// write NCBI's `-num_threads 1` bytes in each program (NCBI prints each alignment once per
/// thread, or aborts), also where blastn runs a TBLASTX strategy.
#[test]
fn a_strategys_subjects_are_searched_once_at_any_thread_count() {
    let cases = cases();
    let names = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
        .map(|program| format!("I.{program}.sub.alone"));
    for name in names
        .iter()
        .map(String::as_str)
        .chain(["C.blastn.from.tblastx.sub"])
    {
        let case = cases
            .iter()
            .find(|case| case.name == name)
            .unwrap_or_else(|| panic!("row {name}"));
        for threads in ["2", "4", "8"] {
            check(case, Some(threads));
        }
    }
}
