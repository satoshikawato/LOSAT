//! Two TBLASTX switches of the strict set on whole runs of the binary: `LOSAT_X_TBXNEXT` (each
//! thread of the pool takes the next subject when it is free,
//! `src/algorithm/tblastx/x_subject_next.rs`) and `LOSAT_X_TBXFRAMES` (a subject's preliminary
//! frames made in windows and freed after the seed stage, the frames of its letters made only when
//! HSPs reach the re-evaluation, `Subjects::x_preliminary_frames` / `x_letter_frames`). Both are on
//! by default (`src/utils/xswitch.rs`); `LOSAT_X_<NAME>=0` turns one off and `LOSAT_X_DEFAULT=0`
//! runs the previous paths.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ```
//! Neither switch may change a report: the tests compare the reports by default (both on), with
//! each one off, and with both set on the previous paths, at one and several threads, against the
//! one-thread search on the previous paths (`LOSAT_X_DEFAULT=0`), on a set laid out as
//! a genome database is (two long records first, then many short ones), with ambiguity codes in
//! the long records (N blocks, single IUPAC codes, one of them in a run across a translation
//! window), through `-subject` and through a database of the same records (`LOSAT makedb`).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path, then `switches`.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    for (name, value) in switches {
        command.env(name, value);
    }
    command
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-tblastx-subject-next-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn letters(name: &str) -> Vec<u8> {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

fn record(defline: &str, letters: &[u8]) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The first 6 kb of MjeNMV; the MelaMJNV and PemoMJNVA genomes with ambiguity codes, then 60
/// pieces of 1,500 bases of both (every fourth with one N).
fn inputs(dir: &Path) -> (PathBuf, PathBuf) {
    let query = dir.join("q.fna");
    std::fs::write(&query, record("mje6k", &letters("MjeNMV.fasta")[..6000])).unwrap();
    let mut mela = letters("MelaMJNV.fasta");
    let pemo = letters("PemoMJNVA.fasta");
    let iupac = b"NRYKMSWBDHV";
    for (n, at) in (1_000..mela.len()).step_by(997).enumerate() {
        mela[at] = iupac[n % iupac.len()];
    }
    // N blocks, one across the border of the first translation window (196,608 bases).
    for (start, len) in [(20_000, 300), (196_500, 250), (120_000, 40)] {
        if start + len <= mela.len() {
            mela[start..start + len].fill(b'N');
        }
    }
    let mut set =
        record("mela MelaMJNV with ambiguity codes", &mela) + &record("pemo PemoMJNVA", &pemo);
    for piece in 0..60 {
        let source = if piece % 2 == 0 { &mela } else { &pemo };
        let start = (piece * 7_919) % (source.len() - 1_500);
        let mut letters = source[start..start + 1_500].to_vec();
        if piece % 4 == 1 {
            letters[700] = b'N';
        }
        set += &record(&format!("piece{piece}"), &letters);
    }
    let set_file = dir.join("set.fna");
    std::fs::write(&set_file, set).unwrap();
    (query, set_file)
}

fn makedb(dir: &Path, fasta_file: &Path) -> String {
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat(&[])
        .arg("makedb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

fn tblastx(
    query: &Path,
    subjects: &[&str],
    outfmt: &str,
    threads: u32,
    switches: &[(&str, &str)],
) -> Vec<u8> {
    let output = losat(switches)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .args(subjects)
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{switches:?} {subjects:?} {outfmt} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

const NEXT: (&str, &str) = ("LOSAT_X_TBXNEXT", "1");
const FRAMES: (&str, &str) = ("LOSAT_X_TBXFRAMES", "1");
const NEXT_OFF: (&str, &str) = ("LOSAT_X_TBXNEXT", "0");
const FRAMES_OFF: (&str, &str) = ("LOSAT_X_TBXFRAMES", "0");
/// The previous paths: every switch of the strict set off unless it is set.
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");

fn same_reports(subjects: &[&str], query: &Path) {
    for outfmt in ["6", "0"] {
        let reference = tblastx(query, subjects, outfmt, 1, &[PREVIOUS]);
        assert!(
            reference.iter().filter(|&&b| b == b'\n').count() > 20,
            "the search finds hits"
        );
        let runs: &[(u32, &[(&str, &str)])] = if outfmt == "6" {
            &[
                (1, &[]),
                (2, &[]),
                (4, &[]),
                (8, &[]),
                (4, &[NEXT_OFF]),
                (8, &[FRAMES_OFF]),
                (1, &[FRAMES_OFF]),
                (2, &[NEXT_OFF, FRAMES_OFF]),
                (2, &[PREVIOUS, NEXT]),
                (1, &[PREVIOUS, FRAMES]),
                (3, &[PREVIOUS, NEXT, FRAMES]),
            ]
        } else {
            &[(1, &[]), (4, &[]), (4, &[PREVIOUS, NEXT, FRAMES])]
        };
        for &(threads, switches) in runs {
            let report = tblastx(query, subjects, outfmt, threads, switches);
            assert!(
                report == reference,
                "{switches:?} threads {threads} outfmt {outfmt} {subjects:?}: the report differs"
            );
        }
    }
}

#[test]
fn subject_records_give_the_same_reports_with_the_switches() {
    let dir = temp_dir("subject");
    let (query, set) = inputs(&dir);
    same_reports(&["-subject", set.to_str().unwrap()], &query);
}

#[test]
fn a_database_gives_the_same_reports_with_the_switches() {
    let dir = temp_dir("db");
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set);
    same_reports(&["-db", &db], &query);
}
