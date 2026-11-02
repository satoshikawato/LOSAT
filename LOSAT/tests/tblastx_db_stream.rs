//! TBLASTX's `-db` search with each database sequence read where it is searched and reported
//! (`common::x_db_stream`, `tblastx::subjects::Subjects::Stream`) on a whole run of the binary.
//!
//! The stream replaced a copy of the whole database made before the first query batch; the unit
//! tests compare its frames and its errors with those of the copy, and the frozen NCBI `-db`
//! fixtures hold the reports. The properties under test are LOSAT's own: the reports at four
//! threads are those at one thread, in one query batch and in batches of one query (`BATCH_SIZE`:
//! each batch reads the sequences it reports), in every report form that reads the letters, with
//! each thread taking the next subject (`LOSAT_X_TBXNEXT`), with the board of the pairs
//! (`LOSAT_X_TBXBOARD`) and with the batches searched while the one before is reported
//! (`LOSAT_X_TBXPIPE`) and without; a custom tabular report that reads no letters gives the
//! columns of one that reads them; and a damaged database fails the same batch with the same error
//! at every thread count. The database has several volumes, a sequence with N blocks and single
//! ambiguity codes in a block of its own, and short ones in blocks of several.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
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
        "losat-tblastx-db-stream-{tag}-{}",
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

/// The MelaMJNV genome with an ambiguity code every 1000 bases (runs of one, the old format), and
/// the number of its runs.
fn mela_with_codes() -> (Vec<u8>, usize) {
    let mut mela = letters("MelaMJNV.fasta");
    let at: Vec<usize> = (500..mela.len()).step_by(1000).collect();
    for &at in &at {
        mela[at] = b"NRYKMSWBDHV"[at / 1000 % 11];
    }
    (mela, at.len())
}

/// The subjects, in volumes of at most 20,000 bytes: the PemoMJNVA genome with N blocks (a block of
/// its own), 40 pieces of 1,500 bases of both genomes (every fourth with one N; blocks of several),
/// and last, alone in its volume, MelaMJNV with its codes (`mela_with_codes`). The queries: 60 N
/// (no valid context), the first 6 kb of MjeNMV and 3 kb of each genome.
fn inputs(dir: &Path) -> (String, PathBuf) {
    let (mela, _) = mela_with_codes();
    let mut pemo = letters("PemoMJNVA.fasta");
    for (start, length) in [(20_000, 300), (60_000, 40)] {
        pemo[start..start + length].fill(b'N');
    }
    let mut set = record("pemo PemoMJNVA with N blocks", &pemo);
    for piece in 0..40 {
        let source = if piece % 2 == 0 { &mela } else { &pemo };
        let start = (piece * 7_919) % (source.len() - 1_500);
        let mut letters = source[start..start + 1_500].to_vec();
        if piece % 4 == 1 {
            letters[700] = b'N';
        }
        set += &record(&format!("piece{piece}"), &letters);
    }
    set += &record("mela_amb MelaMJNV with ambiguity codes", &mela);
    let set_file = dir.join("set.fna");
    std::fs::write(&set_file, set).unwrap();
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat(&[])
        .arg("makeblastdb")
        .arg("-in")
        .arg(&set_file)
        .args(["-dbtype", "nucl", "-out", &db, "-max_file_sz", "20000"])
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let queries = record("nn all N", &[b'N'; 60])
        + &record("mje6k", &letters("MjeNMV.fasta")[..6000])
        + &record("mela3k", &mela[100_000..103_000])
        + &record("pemo3k", &letters("PemoMJNVA.fasta")[30_000..33_000]);
    let query = dir.join("q.fna");
    std::fs::write(&query, queries).unwrap();
    (db, query)
}

/// One search of `query` on `db` with `outfmt` at `threads` threads; returns the exit status, the
/// report and the standard error.
fn tblastx(
    query: &Path,
    db: &str,
    outfmt: &str,
    threads: u32,
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let output = losat(env)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .args(["-db", db, "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// A custom tabular report of every field that reads the letters of the subjects.
const LETTERS: &str =
    "6 qseqid sseqid pident nident mismatch positive ppos qstart qend sstart send btop qseq sseq evalue bitscore";

// The query batches change the rows of TBLASTX (as they did with the copy of the database), so
// each search is compared with the search at one thread in the same batches.
#[test]
fn the_reports_do_not_depend_on_the_threads_or_the_switches() {
    let dir = temp_dir("reports");
    let (db, query) = inputs(&dir);
    let batches = [("BATCH_SIZE", "1")];
    for outfmt in [LETTERS, "0", "7", "5"] {
        for env in [&[][..], &batches[..]] {
            let reference = tblastx(&query, &db, outfmt, 1, env);
            assert_eq!(reference.0, Some(0), "{outfmt}: {}", reference.2);
            assert!(reference.1.len() > 10_000, "{outfmt}: a short report");
            assert!(
                tblastx(&query, &db, outfmt, 4, env) == reference,
                "{outfmt}: 4 threads {env:?}"
            );
        }
    }
    let reference = tblastx(&query, &db, LETTERS, 1, &batches);
    for switch in ["LOSAT_X_TBXNEXT", "LOSAT_X_TBXBOARD", "LOSAT_X_TBXPIPE"] {
        let env = [("BATCH_SIZE", "1"), (switch, "0")];
        assert!(
            tblastx(&query, &db, LETTERS, 4, &env) == reference,
            "{switch}=0"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A custom tabular report without a field that reads the letters of the subjects does not read
/// the sequences for the report (`tabular::reads_letters`), with comments and as CSV too; its
/// columns are those of the report that reads them.
#[test]
fn a_report_without_letters_gives_the_columns_of_one_with_them() {
    let dir = temp_dir("columns");
    let (db, query) = inputs(&dir);
    let env = [("BATCH_SIZE", "1")];
    let with = tblastx(&query, &db, LETTERS, 4, &env);
    assert_eq!(with.0, Some(0), "{}", with.2);
    let columns = |report: &[u8], delimiter: char, picked: &[usize]| -> Vec<String> {
        String::from_utf8_lossy(report)
            .lines()
            .filter(|row| !row.starts_with('#'))
            .map(|row| {
                let fields: Vec<&str> = row.split(delimiter).collect();
                picked
                    .iter()
                    .map(|&i| fields[i])
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect()
    };
    let rows = columns(&with.1, '\t', &[0, 1, 14, 15]);
    assert!(rows.len() > 100, "{} rows", rows.len());
    for (outfmt, delimiter) in [
        ("6 qseqid sseqid evalue bitscore score", '\t'),
        ("7 qseqid sseqid evalue bitscore", '\t'),
        ("10 qseqid sseqid evalue bitscore", ','),
    ] {
        let without = tblastx(&query, &db, outfmt, 4, &env);
        assert_eq!(without.0, Some(0), "{outfmt}: {}", without.2);
        assert_eq!(
            columns(&without.1, delimiter, &[0, 1, 2, 3]),
            rows,
            "{outfmt}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one (as
// `tblastx_db_search.rs` `the_sequence_file_is_read_by_the_first_batch_with_a_valid_context`).
// MelaMJNV, a block of its own, is checked without its 2-bit data
// (`common::database::x_check_nucleotide_stores`, `SeqDb::x_nucleotide_check`): only the last
// 2-bit byte and the ambiguity data. The reports of the batches before the failed one are written
// with NCBI's errors (an engine error) and not with LOSAT's refusals, whether the reports are
// written batch by batch (outfmt 7) or after the loop (outfmt 5).
#[test]
fn a_damaged_database_fails_the_same_batch_with_the_same_error() {
    let dir = temp_dir("damaged");
    let (db, query) = inputs(&dir);
    let (mela, runs) = mela_with_codes();
    let last = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.to_string_lossy().ends_with(".nsq"))
        .max()
        .unwrap();
    let original = std::fs::read(&last).unwrap();
    // The sentinel byte, the 2-bit data, then the ambiguity data: a count word and a word per run.
    let ambiguity = 1 + mela.len() / 4 + 1;
    assert_eq!(original.len(), ambiguity + 4 + 4 * runs, "one sequence");
    let count =
        |bytes: &[u8]| u32::from_be_bytes(bytes[ambiguity..ambiguity + 4].try_into().unwrap());
    let mut long_count = original.clone();
    long_count[ambiguity..ambiguity + 4].copy_from_slice(&(count(&original) + 100).to_be_bytes());
    let mut remainder = original.clone();
    remainder[ambiguity - 1] ^= 1;
    let env = [("BATCH_SIZE", "1")];
    for (what, bytes, error) in [
        (
            "count past the end",
            Some(long_count),
            Some("the ambiguity data of database sequence 41 cannot be read"),
        ),
        ("remainder", Some(remainder), None),
        (
            "cut in the ambiguity data",
            Some(original[..ambiguity + 8].to_vec()),
            Some("passes the end of the file"),
        ),
        (
            "cut in the 2-bit data",
            Some(original[..ambiguity / 2].to_vec()),
            Some("passes the end of the file"),
        ),
        ("missing", None, Some("Database memory map file error")),
    ] {
        match &bytes {
            Some(bytes) => std::fs::write(&last, bytes).unwrap(),
            None => std::fs::remove_file(&last).unwrap(),
        }
        for outfmt in ["7", "5"] {
            let what = format!("{what} outfmt {outfmt}");
            let one = tblastx(&query, &db, outfmt, 1, &env);
            match error {
                Some(error) => {
                    assert_ne!(one.0, Some(0), "{what}: {}", one.2);
                    assert!(one.2.contains(error), "{what}: {}", one.2);
                    // NCBI's error comes after the first batch's report, LOSAT's refusal stops
                    // the run before it.
                    let engine = one.0 == Some(3);
                    assert_eq!(engine, error == "Database memory map file error", "{what}");
                    let report = String::from_utf8_lossy(&one.1);
                    assert_eq!(report.contains("nn all N"), engine, "{what}: {report}");
                    assert!(!report.contains("mje6k"), "{what}: {report}");
                }
                None => {
                    assert_eq!(one.0, Some(0), "{what}: {}", one.2);
                    assert!(!one.1.is_empty(), "{what}");
                }
            }
            for (threads, switch) in [
                (2, None),
                (4, None),
                (4, Some("LOSAT_X_TBXPIPE")),
                (4, Some("LOSAT_X_TBXBOARD")),
            ] {
                let env: Vec<_> = env
                    .iter()
                    .copied()
                    .chain(switch.map(|s| (s, "0")))
                    .collect();
                assert!(
                    tblastx(&query, &db, outfmt, threads, &env) == one,
                    "{what}: {threads} threads {switch:?}"
                );
            }
        }
        std::fs::write(&last, &original).unwrap();
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
