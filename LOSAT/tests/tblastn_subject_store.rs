//! TBLASTN subjects held packed (`LOSAT_X_TBNSTORE`, `src/algorithm/tblastn/x_subject_store.rs`,
//! `src/sequence/packed_subject.rs`), on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
//! ```c
//!     virtual SBlastSequence GetCompressedPlusStrand() {
//!         SBlastSequence retval(size());
//!         string ncbi4na = kEmptyStr;
//!         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
//!         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
//!         return retval;
//! ```
//! NCBI keeps the subjects in its scope and makes each one's search forms from it. With
//! `LOSAT_X_TBNSTORE` set (to any value; the switch is not in the strict set) the reader's records
//! keep their deflines and the letters are held packed (2 bits per base, the runs of ambiguous
//! letters, the lower-case intervals), made again when a subject is searched, redone or reported.
//! The switch must change nothing a user sees, so the tests compare the report and the messages
//! with those of a run without the variable, and read the `[X_STATS] TBNSTORE` line to check that
//! the subjects were held packed.
//!
//! The subjects are pieces of `tests/fasta/MelaMJNV.fasta` with what the packed form keeps apart
//! from the bases: lower-case runs (soft masking with `-lcase_masking`), every ambiguous IUPAC
//! letter and `-`, runs of `N`, one subject too short for a word, two of a few codons and one
//! without letters (its warning comes from the length the store keeps). The queries are proteins
//! of `tests/fasta/PemoMJNVA.faa`: ten of them (one query batch), or all 115 (four batches).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const PIECE: usize = 6_000;

/// The records of a FASTA file as (title line, letters).
fn records(path: &Path) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).expect("fixture");
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        match line.strip_prefix('>') {
            Some(title) => out.push((title.to_string(), String::new())),
            None => out
                .last_mut()
                .expect("title line first")
                .1
                .push_str(line.trim()),
        }
    }
    out
}

/// The files of the tests, written once: ten proteins, all proteins, the subjects.
fn inputs() -> &'static (PathBuf, PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_subject_store");
        std::fs::create_dir_all(&dir).expect("work directory");
        let proteins = records(&fasta.join("PemoMJNVA.faa"));
        let mut ten = String::new();
        let mut all = String::new();
        for (index, (title, letters)) in proteins.iter().enumerate() {
            let record = format!(">{title}\n{letters}\n");
            if index < 10 {
                ten.push_str(&record);
            }
            all.push_str(&record);
        }
        let genome = records(&fasta.join("MelaMJNV.fasta"))[0]
            .1
            .clone()
            .into_bytes();
        assert!(genome.len() >= 4 * 60_000 + PIECE, "genome too short");
        let piece = |index: usize| genome[index * 60_000..index * 60_000 + PIECE].to_vec();
        let lowercase = |letters: &mut [u8], from: usize, to: usize| {
            letters[from..to].make_ascii_lowercase();
        };
        let plain = piece(0);
        let mut lower = piece(1);
        lowercase(&mut lower, 0, 40);
        lowercase(&mut lower, 1_500, 1_800);
        lowercase(&mut lower, 4_200, 4_230);
        lowercase(&mut lower, PIECE - 25, PIECE);
        let mut ambiguous = piece(2);
        for (count, position) in (11..PIECE).step_by(97).enumerate() {
            ambiguous[position] = b"NRYKMSWBDHV-nrykmswbdhv"[count % 23];
        }
        lowercase(&mut ambiguous, 2_000, 2_300);
        let mut runs = piece(3);
        runs[3_000..3_040].fill(b'N');
        runs[5_000..5_007].fill(b'n');
        lowercase(&mut runs, 900, 1_200);
        let mut subjects = String::new();
        let mut long = String::new();
        for (title, letters) in [
            ("plain", plain.as_slice()),
            ("lower", lower.as_slice()),
            ("no_letters", b"".as_slice()),
            ("ambiguous", ambiguous.as_slice()),
            ("runs", runs.as_slice()),
            ("too_short_for_a_word", b"ACGTN".as_slice()),
            ("one_word", b"ATGNCGTRAAC".as_slice()),
            ("few_codons", b"ATGAAACGTYTGGCNnnAAATAG".as_slice()),
        ] {
            let record = format!(
                ">{title}\n{}\n",
                String::from_utf8(letters.to_vec()).expect("ASCII letters")
            );
            if letters.len() == PIECE {
                long.push_str(&record);
            }
            subjects.push_str(&record);
        }
        let files = (
            dir.join("ten.faa"),
            dir.join("all.faa"),
            dir.join("subjects.fna"),
        );
        std::fs::write(&files.0, ten).expect("query file");
        std::fs::write(&files.1, all).expect("query file");
        std::fs::write(&files.2, subjects).expect("subject file");
        std::fs::write(long_subjects(&files.2), long).expect("subject file");
        files
    })
}

/// The subject file of the four pieces alone (a subject range must lie in every subject).
fn long_subjects(subjects: &Path) -> PathBuf {
    subjects.with_file_name("long.fna")
}

/// One TBLASTN run (report on standard output) with `options` and with `variables` and
/// `LOSAT_X_STATS` as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build
/// with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so
/// the caller's value is kept, or a file of the tests is named. Returns the report, the standard
/// error without the `[X_STATS]` lines, and the `[X_STATS] TBNSTORE` line.
fn tblastn(
    query: &Path,
    subjects: &Path,
    options: &[&str],
    variables: &[(&str, &str)],
) -> (Vec<u8>, String, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects)
        .args(options)
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{options:?} {variables:?}: {stderr}"
    );
    let messages: String = stderr
        .lines()
        .filter(|line| !line.starts_with("[X_STATS]"))
        .map(|line| format!("{line}\n"))
        .collect();
    let stats: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNSTORE"))
        .collect();
    assert_eq!(stats.len(), 1, "one TBNSTORE line:\n{stderr}");
    (output.stdout, messages, stats[0].to_string())
}

/// The number after `field=` in a `[X_STATS]` line.
fn counter(line: &str, field: &str) -> u64 {
    line.split_whitespace()
        .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no {field} in {line}"))
        .parse()
        .expect("a number")
}

/// The subjects named in the second column of a tabular report.
fn hit_subjects(report: &[u8]) -> std::collections::BTreeSet<String> {
    String::from_utf8_lossy(report)
        .lines()
        .filter(|row| !row.starts_with('#'))
        .filter_map(|row| row.split('\t').nth(1).map(str::to_string))
        .collect()
}

const STORE: (&str, &str) = ("LOSAT_X_TBNSTORE", "1");

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
// ```c
//         catch(CBlastException & e ) {
//         	// Skip bad subject sequence
//         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
//         		seqblk_vec->push_back(subj);
//         ...
//         		warning += "Subject sequence contains no data";
// ```
// Every search below gives one report and the same messages (the warning of the subject without
// letters among them) without the variable and with it: tabular, commented tabular and pairwise
// reports, the default options and `-lcase_masking` at one and three threads, another genetic
// code, composition mode 0, and four query batches. With the variable the eight subjects are held
// packed and their letters are made again (at least once per searched subject); without it none.
#[test]
fn the_report_and_the_messages_are_the_same_with_the_subjects_held_packed() {
    let (ten, all, subjects) = inputs();
    let searches: [(&Path, &[&str]); 8] = [
        (ten, &["-outfmt", "6", "-num_threads", "1"]),
        (ten, &["-outfmt", "6", "-num_threads", "3"]),
        (
            ten,
            &["-outfmt", "7", "-num_threads", "1", "-lcase_masking"],
        ),
        (
            ten,
            &["-outfmt", "0", "-num_threads", "3", "-lcase_masking"],
        ),
        (
            ten,
            &["-outfmt", "0", "-num_threads", "1", "-db_gencode", "4"],
        ),
        (
            ten,
            &[
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-comp_based_stats",
                "0",
            ],
        ),
        (
            ten,
            &[
                "-outfmt",
                "6",
                "-num_threads",
                "3",
                "-comp_based_stats",
                "0",
            ],
        ),
        (
            all,
            &[
                "-outfmt",
                "6",
                "-num_threads",
                "3",
                "-lcase_masking",
                "-db_gencode",
                "11",
            ],
        ),
    ];
    for (query, options) in searches {
        let (expected, messages, stats) = tblastn(query, subjects, options, &[]);
        assert_eq!(counter(&stats, "subjects"), 0, "{options:?}");
        assert_eq!(counter(&stats, "letters_made"), 0, "{options:?}");
        assert!(
            messages.contains("no_letters: Subject sequence contains no data"),
            "{options:?}: {messages}"
        );
        // the tabular reports (`-outfmt` is the first option of every search)
        if options[1] != "0" {
            let hits = hit_subjects(&expected);
            for subject in ["plain", "lower", "ambiguous", "runs"] {
                assert!(hits.contains(subject), "{options:?}: no hit on {subject}");
            }
        }

        let (report, store_messages, stats) = tblastn(query, subjects, options, &[STORE]);
        assert!(
            report == expected,
            "{options:?}: the packed subjects change the report"
        );
        assert_eq!(store_messages, messages, "{options:?}");
        assert_eq!(counter(&stats, "subjects"), 8, "{options:?}");
        assert!(counter(&stats, "stored_bytes") > 0, "{options:?}");
        // the preliminary stage reads the resolved bases unless it reads the case; the later
        // stages of the subjects with hits read the letters
        if options.contains(&"-lcase_masking") {
            assert!(counter(&stats, "letters_made") >= 7, "{options:?}: {stats}");
            assert_eq!(counter(&stats, "resolved_made"), 0, "{options:?}");
        } else {
            assert!(
                counter(&stats, "resolved_made") >= 7,
                "{options:?}: {stats}"
            );
            assert!(counter(&stats, "letters_made") >= 4, "{options:?}: {stats}");
        }
    }
}

// No NCBI counterpart: where the switch acts. A subject range (`-subject_loc`) keeps the
// letters (the records are cut before the search), so the variable changes nothing there and no
// subject is held packed. With `LOSAT_X_DEFAULT=0` (the reference loops) the subjects are still
// held packed and the report is the same.
#[test]
fn a_subject_range_keeps_the_letters_and_the_reference_loops_read_the_packed_subjects() {
    let (ten, _, subjects) = inputs();
    let long = long_subjects(subjects);
    let ranged = ["-outfmt", "6", "-subject_loc", "100-5000"];
    let (expected, messages, _) = tblastn(ten, &long, &ranged, &[]);
    let (report, store_messages, stats) = tblastn(ten, &long, &ranged, &[STORE]);
    assert!(report == expected);
    assert_eq!(store_messages, messages);
    assert_eq!(counter(&stats, "subjects"), 0);
    assert_eq!(counter(&stats, "letters_made"), 0);
    assert_eq!(counter(&stats, "resolved_made"), 0);

    let options = ["-outfmt", "6", "-num_threads", "1"];
    let (expected, messages, _) = tblastn(ten, subjects, &options, &[]);
    let previous = ("LOSAT_X_DEFAULT", "0");
    let (report, store_messages, stats) = tblastn(ten, subjects, &options, &[previous, STORE]);
    assert!(report == expected);
    assert_eq!(store_messages, messages);
    assert_eq!(counter(&stats, "subjects"), 8);
    assert!(counter(&stats, "letters_made") + counter(&stats, "resolved_made") >= 7);
}
