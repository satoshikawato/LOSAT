//! `LOSAT_X_BNEXTFIRST` (`pair_extend::extend_offset_pairs` extends a megablast table's lookup
//! hit to its word before it looks up the hit's context, where the lookup width is below the word
//! size): the report of each lookup width and stride, with `-subject` and `-db`, equals the
//! report with `=0` at one thread (where the side is the one the table's widths give:
//! `LOSAT_X_REACH`) and the default's at four threads. Which width a search gets is NCBI's choice
//! from the word size and the query's words (`BlastChooseNaLookupTable`, blast_nalookup.c:46-180,
//! and the stride `word_length - lut_word_length + 1`, :1334). That those bytes equal NCBI's is the
//! business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_BNEXTFIRST";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-extfirst-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// The letters of the first record of a FASTA file.
fn letters(path: &Path) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

fn database(dir: &TestDir, subject: &Path) -> PathBuf {
    let db = dir.path("subject_db");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(subject)
        .args(["-dbtype", "nucl", "-out"])
        .arg(&db)
        .output()
        .expect("run makeblastdb");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

fn search(
    query: &Path,
    subjects: (&str, &Path),
    options: &[&str],
    threads: u32,
    vars: &[(&str, &str)],
) -> Vec<u8> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg(subjects.0)
        .arg(subjects.1)
        .args(options)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{options:?} {vars:?} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// Each (query letters of Sakai from 1 Mb, options, the width and stride NCBI takes for it,
/// whether the width is below the word size).
const CASES: &[(usize, &[&str], &str, bool)] = &[
    (
        20_000,
        &["-task", "blastn", "-word_size", "10"],
        "width 10 of 10",
        false,
    ),
    (
        6_000,
        &["-task", "blastn", "-word_size", "10"],
        "width 9 of 10, stride 2",
        true,
    ),
    (
        50_000,
        &["-task", "blastn"],
        "width 10 of 11, stride 2",
        true,
    ),
    (
        15_000,
        &["-task", "blastn", "-word_size", "12"],
        "width 10 of 12, stride 3",
        true,
    ),
    (200_000, &["-task", "blastn"], "width 11 of 11", false),
    (
        100_000,
        &["-task", "megablast"],
        "width 11 of 28, stride 18",
        true,
    ),
    (
        50_000,
        &["-task", "megablast", "-word_size", "13"],
        "width 11 of 13, stride 3",
        true,
    ),
    (
        200_000,
        &["-task", "megablast"],
        "width 12 of 28, stride 17",
        true,
    ),
    (
        200_000,
        &[
            "-task",
            "megablast",
            "-word_size",
            "32",
            "-ungapped",
            "-perc_identity",
            "95",
        ],
        "width 12 of 32, stride 21 (mlst's search)",
        true,
    ),
];

#[test]
fn the_word_extension_first_gives_the_same_report() {
    let dir = TestDir::new("each");
    let sakai = letters(&fasta("Sakai.fna"));
    let genome = letters(&fasta("MG1655.fna"));
    // Several records, so that `-db` scans many subjects with one table.
    let subject = dir.path("mg1655_500k.fna");
    let records: String = (0..5)
        .map(|i| {
            let from = 1_000_000 + i * 100_000;
            format!(">mg{i}\n{}\n", &genome[from..from + 100_000])
        })
        .collect();
    std::fs::write(&subject, records).expect("write the subjects");
    let db = database(&dir, &subject);
    for (n, &(length, options, what, extends)) in CASES.iter().enumerate() {
        let query = dir.path(&format!("sakai_{n}.fna"));
        std::fs::write(
            &query,
            format!(">sakai_{n}\n{}\n", &sakai[1_000_000..1_000_000 + length]),
        )
        .expect("write the query");
        for subjects in [("-subject", subject.as_path()), ("-db", db.as_path())] {
            let what = format!("{what} {}", subjects.0);
            let expected = search(&query, subjects, options, 1, &[(SWITCH, "0")]);
            assert!(!expected.is_empty(), "{what}: no hits");
            let reach = reach_support::Reach::new("extfirst");
            let report = search(&query, subjects, options, 1, &[reach.var()]);
            assert!(report == expected, "{what}: default differs");
            assert_eq!(reach.sides(SWITCH), (extends, !extends), "{what}: the side");
            // The pool's stages at more threads.
            let report = search(&query, subjects, options, 4, &[]);
            assert!(report == expected, "{what} threads 4: default differs");
        }
    }
}
