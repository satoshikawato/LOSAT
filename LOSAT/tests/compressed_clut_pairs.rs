//! The neighbouring words of the compressed lookup tables take the cells of their two innermost
//! positions from a row per pair of query letters (`LOSAT_X_CLUTPAIRS`,
//! `src/algorithm/tblastx/lookup/x_clut_direct.rs` `XPairRow`, `XWords::neighbors_paired`) on whole
//! searches of the three programs that build those tables.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053
//! The switch is in the strict set (on unless it is `0`). A query offset gets the same cells as
//! from the loops, in another order within the offset, which changes no table (the lib tests
//! `paired_cells_are_the_loops_cells` and `pair_tables_hold_every_matrix` compare the cells and the
//! rows with the loops'), so the property under test is that the report and the messages do not
//! depend on it: blastp `-task blastp-fast`, blastx `-task blastx-fast` and tblastn `-task
//! tblastn-fast` (word size 5), blastx and tblastn with `-word_size 6`, at 1 and 2 threads (the
//! serial builders and the pool's). That reports equal NCBI's is the business of the programs'
//! fixtures.
//!
//! With `LOSAT_X_STATS` the `[X_STATS] CLUTPAIRS tables= rows=` line counts the pair tables and the
//! rows the searches made, so the default arm is shown to read the rows, and the `=0` arm prints no
//! such line. `LOSAT_X_STATS` changes neither the report nor the other messages. A build without
//! the `parallel` feature runs the one-thread cases only.
//!
//! The proteins are the first 30 of AP027131 (queries) and the first 150 of the related AP027133
//! (subjects); the nucleotide sequences are the first 30,000 bases of AP027131 (blastx's query)
//! and the first 60,000 of AP027133 (tblastn's subject).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_CLUTPAIRS", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_compressed_clut_pairs_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Records `from..to` (0-based) of a FASTA file, as text.
fn records(name: &str, from: usize, to: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut out = String::new();
    let mut record = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            record += 1;
        }
        if record > from && record <= to {
            out.push_str(line);
            out.push('\n');
        }
    }
    assert!(record >= to, "{name} has fewer than {to} records");
    out
}

/// The first `bases` letters of a one-record FASTA file, as a record of 60 letters per line.
fn prefix(name: &str, bases: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut lines = text.lines();
    let title = lines.next().expect("a title line");
    let letters: String = lines
        .flat_map(|line| line.trim().chars())
        .take(bases)
        .collect();
    assert_eq!(letters.len(), bases, "{name} is shorter than {bases}");
    let mut out = format!("{title}\n");
    for chunk in letters.as_bytes().chunks(60) {
        out.push_str(std::str::from_utf8(chunk).unwrap());
        out.push('\n');
    }
    out
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    proteins: PathBuf,
    subjects: PathBuf,
    genome: PathBuf,
    subject_genome: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let proteins = dir.join("proteins.faa");
        let subjects = dir.join("subjects.faa");
        let genome = dir.join("genome.fna");
        let subject_genome = dir.join("subject_genome.fna");
        std::fs::write(&proteins, records("AP027131.faa", 0, 30)).unwrap();
        std::fs::write(&subjects, records("AP027133.faa", 0, 150)).unwrap();
        std::fs::write(&genome, prefix("AP027131.fasta", 30_000)).unwrap();
        std::fs::write(&subject_genome, prefix("AP027133.fasta", 60_000)).unwrap();
        Inputs {
            dir,
            proteins,
            subjects,
            genome,
            subject_genome,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The binary with `env` as its only variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// names `worker_log`, in a directory the caller removes: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn losat(env: &[(&str, &str)], worker_log: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
        .env("LOSAT_BLASTX_WORKER_LOG", worker_log)
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

impl Run {
    /// Standard error without the `[X_STATS]` lines.
    fn messages(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[X_STATS]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// The numbers after `tables=` and `rows=` in the `[X_STATS] CLUTPAIRS` line, or `None`
    /// without that line.
    fn pairs(&self) -> Option<(u64, u64)> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with("[X_STATS] CLUTPAIRS "))?;
        let number = |name: &str| -> u64 {
            let key = format!(" {name}=");
            let at = line
                .find(&key)
                .unwrap_or_else(|| panic!("no {name} in {line}"))
                + key.len();
            line[at..]
                .split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        Some((number("tables"), number("rows")))
    }
}

/// `program` with `args` and `env` (the worker log is in a directory of its own, removed after
/// the run).
fn run(program: &str, args: &[&str], env: &[(&str, &str)]) -> Run {
    let worker_log = temp_path("worker");
    std::fs::create_dir_all(&worker_log).unwrap();
    let output = losat(env, &worker_log.join("worker.log"))
        .arg(program)
        .args(args)
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_dir_all(&worker_log);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// Default and `=0` give the same report and messages for the compressed tables of blastp,
/// blastx and tblastn, at 1 and 2 threads; the default reads pair rows, `=0` none.
#[test]
fn compressed_table_reports_do_not_depend_on_the_pair_rows() {
    let inputs = Inputs::new();
    let path = |p: &PathBuf| p.to_str().unwrap().to_owned();
    let (proteins, subjects) = (path(&inputs.proteins), path(&inputs.subjects));
    let (genome, subject_genome) = (path(&inputs.genome), path(&inputs.subject_genome));
    let cases: Vec<(&str, Vec<&str>)> = vec![
        (
            "blastp",
            vec![
                "-task",
                "blastp-fast",
                "-query",
                &proteins,
                "-subject",
                &subjects,
            ],
        ),
        (
            "blastx",
            vec![
                "-task",
                "blastx-fast",
                "-query",
                &genome,
                "-subject",
                &subjects,
            ],
        ),
        (
            "blastx",
            vec!["-word_size", "6", "-query", &genome, "-subject", &subjects],
        ),
        (
            "tblastn",
            vec![
                "-task",
                "tblastn-fast",
                "-query",
                &proteins,
                "-subject",
                &subject_genome,
            ],
        ),
        (
            "tblastn",
            vec![
                "-word_size",
                "6",
                "-query",
                &proteins,
                "-subject",
                &subject_genome,
            ],
        ),
    ];
    for (program, case_args) in &cases {
        for threads in ["1", "2"] {
            // A build without the `parallel` feature refuses two threads.
            if threads == "2" && !cfg!(feature = "parallel") {
                continue;
            }
            let case = format!("{program} {} n{threads}", case_args[..2].join(" "));
            let mut args = case_args.clone();
            args.extend(["-outfmt", "6", "-evalue", "10", "-num_threads", threads]);
            let on = run(program, &args, &[]);
            assert_eq!(on.code, Some(0), "{case}: {}", on.stderr);
            assert!(!on.stdout.is_empty(), "{case}: no hits");
            assert_eq!(
                on.pairs(),
                None,
                "{case}: a CLUTPAIRS line without LOSAT_X_STATS"
            );
            let off = run(program, &args, &[OFF]);
            assert_eq!(off.code, Some(0), "{case} off: {}", off.stderr);
            assert!(off.stdout == on.stdout, "{case} off: report differs");
            assert_eq!(off.stderr, on.stderr, "{case} off: messages");
            let counted = run(program, &args, &[STATS]);
            assert_eq!(counted.code, Some(0), "{case} stats: {}", counted.stderr);
            assert!(counted.stdout == on.stdout, "{case} stats: report differs");
            assert_eq!(counted.messages(), on.stderr, "{case} stats: messages");
            let (tables, rows) = counted
                .pairs()
                .unwrap_or_else(|| panic!("{case}: no CLUTPAIRS line: {}", counted.stderr));
            assert!(
                tables >= 1 && rows >= 1,
                "{case}: tables {tables} rows {rows}"
            );
            let off_counted = run(program, &args, &[OFF, STATS]);
            assert!(
                off_counted.stdout == on.stdout,
                "{case} off stats: report differs"
            );
            assert_eq!(off_counted.pairs(), None, "{case} off: pair rows read");
        }
    }
}
