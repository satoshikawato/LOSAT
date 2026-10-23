//! BLASTN subjects held packed (`LOSAT_X_BNSTORE`, `src/algorithm/blastn/x_subject_store.rs`,
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
//! `LOSAT_X_BNSTORE` on (the strict set since 2026-10-18: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set) the reader's records
//! keep their deflines and the letters are held packed (2 bits per base, the runs of ambiguous
//! letters, the lower-case intervals): the preliminary stage reads the packed bases, the traceback
//! makes the blastna letters of a subject, the pairwise report the letters of each row. The switch
//! must change nothing a user sees, so the tests compare the report and the messages with those of
//! a run with `LOSAT_X_BNSTORE=0`, and read the `[X_STATS] BNSTORE` line to check that the
//! subjects were held packed.
//!
//! The subjects are pieces of `tests/fasta/MelaMJNV.fasta` with what the packed form keeps apart
//! from the bases: lower-case runs (soft masking with `-lcase_masking`, lower-case rows in the
//! pairwise report), every ambiguous IUPAC letter (in the aligned regions, so the traceback scores
//! them and the report shows them), runs of `N`, one subject too short for a word and one without
//! letters (its warning comes from the length the packed form keeps). The queries are parts of the
//! pieces with substitutions and gaps, on both strands.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), and this
/// file's own on 2026-10-18 (perf/bn-memory-rules), each `0`, the default before, unless the
/// search sets it: every search here runs as it was written, with the switches off that it does
/// not name (the default runs in `default_bundle.rs` and
/// `the_subjects_are_held_packed_by_default`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BNSTORE", "0"),
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const PIECE: usize = 6_000;

/// The letters of the first record of a FASTA file.
fn first_record(path: &Path) -> Vec<u8> {
    let text = std::fs::read_to_string(path).expect("fixture");
    text.lines()
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

/// The complement of a letter of `ACGT` (the queries are made from the genome's bases).
fn complement(letter: u8) -> u8 {
    match letter.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// The files of the tests, written once: the queries and the subjects.
fn inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        // One directory per process: nextest runs each test in a process of its own, and a
        // process that wrote these files again while another one's search read them cut its input.
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join("blastn_subject_store")
            .join(std::process::id().to_string());
        std::fs::create_dir_all(&dir).expect("work directory");
        let genome = first_record(&fasta.join("MelaMJNV.fasta"));
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
            ambiguous[position] = b"NRYKMSWBDHVnrykmswbdhv"[count % 22];
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
            ("too_short_for_a_word", b"ACGTNACGTR".as_slice()),
        ] {
            let record = format!(
                ">{title} piece of MelaMJNV\n{}\n",
                String::from_utf8(letters.to_vec()).expect("ASCII letters")
            );
            if letters.len() == PIECE {
                long.push_str(&record);
            }
            subjects.push_str(&record);
        }
        // The queries: 2,000 bases from inside each piece with a substitution every 53 bases,
        // a deletion of 3 and an insertion of 2; every other one reverse-complemented.
        let mut queries = String::new();
        for index in 0..4 {
            let start = index * 60_000 + 1_000 + 300 * index;
            let mut letters: Vec<u8> = genome[start..start + 2_000]
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect();
            for position in (17..letters.len()).step_by(53) {
                letters[position] = complement(letters[position]);
            }
            letters.drain(700..703);
            letters.splice(1_400..1_400, *b"GA");
            if index % 2 == 1 {
                letters = letters
                    .iter()
                    .rev()
                    .map(|&letter| complement(letter))
                    .collect();
            }
            queries.push_str(&format!(
                ">query_{index}\n{}\n",
                String::from_utf8(letters).expect("ASCII letters")
            ));
        }
        let files = (dir.join("queries.fna"), dir.join("subjects.fna"));
        std::fs::write(&files.0, queries).expect("query file");
        std::fs::write(&files.1, subjects).expect("subject file");
        std::fs::write(long_subjects(&files.1), long).expect("subject file");
        files
    })
}

/// The subject file of the four pieces alone (a subject range must lie in every subject).
fn long_subjects(subjects: &Path) -> PathBuf {
    subjects.with_file_name("long.fna")
}

/// One BLASTN run (report on standard output) with `options` and with `variables` and
/// `LOSAT_X_STATS` as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build
/// with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so
/// the caller's value is kept, or a file of the tests is named. Returns the report, the standard
/// error without the `[X_STATS]` lines, and the `[X_STATS] BNSTORE` line.
fn blastn(
    query: &Path,
    subjects: &Path,
    options: &[&str],
    variables: &[(&str, &str)],
) -> (Vec<u8>, String, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects)
        .args(options)
        .env("LOSAT_X_STATS", "1")
        .envs(PREVIOUS_DEFAULT.iter().copied())
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
        .filter(|line| line.starts_with("[X_STATS] BNSTORE"))
        .collect();
    assert_eq!(stats.len(), 1, "one BNSTORE line:\n{stderr}");
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

const STORE: (&str, &str) = ("LOSAT_X_BNSTORE", "1");

/// The options of one search, and its `LOSAT_*` variables besides `LOSAT_X_BNSTORE`.
type Search<'a> = (&'a [&'a str], &'a [(&'a str, &'a str)]);

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
// letters among them) with the variable `0` (the default before 2026-10-18) and `1`: tabular, commented tabular and pairwise
// reports, megablast, blastn and dc-megablast, `-lcase_masking`, one and three threads, and three
// threads with every subject searched before the hand-over (`LOSAT_X_BNSTREAM=0`; it is on by
// default). With `1` the six
// subjects are held packed, the traceback makes the blastna letters of each subject with hits and
// the pairwise report the letters of each row; with `0` none.
#[test]
fn the_report_and_the_messages_are_the_same_with_the_subjects_held_packed() {
    let (queries, subjects) = inputs();
    let searches: [Search; 9] = [
        (&["-outfmt", "6", "-num_threads", "1"], &[]),
        (&["-outfmt", "6", "-num_threads", "3"], &[]),
        (
            &["-outfmt", "7", "-num_threads", "1", "-lcase_masking"],
            &[],
        ),
        (
            &["-outfmt", "0", "-num_threads", "3", "-lcase_masking"],
            &[],
        ),
        (&["-outfmt", "0", "-num_threads", "1"], &[]),
        (
            &["-outfmt", "0", "-num_threads", "1", "-task", "blastn"],
            &[],
        ),
        (
            &[
                "-outfmt",
                "6",
                "-num_threads",
                "3",
                "-task",
                "dc-megablast",
                "-lcase_masking",
            ],
            &[],
        ),
        (
            &["-outfmt", "0", "-num_threads", "3", "-task", "blastn"],
            &[("LOSAT_X_BNSTREAM", "0")],
        ),
        (&["-outfmt", "6", "-num_threads", "1", "-dust", "no"], &[]),
    ];
    for (options, variables) in searches {
        let (expected, messages, stats) = blastn(queries, subjects, options, variables);
        assert_eq!(counter(&stats, "subjects"), 0, "{options:?}");
        assert_eq!(counter(&stats, "blastna_made"), 0, "{options:?}");
        assert!(
            messages.contains("no_letters piece of MelaMJNV: Subject sequence contains no data"),
            "{options:?}: {messages}"
        );
        // the tabular reports (`-outfmt` is the first option of every search)
        if options[1] != "0" {
            let hits = hit_subjects(&expected);
            for subject in ["plain", "lower", "ambiguous", "runs"] {
                assert!(hits.contains(subject), "{options:?}: no hit on {subject}");
            }
        }

        let with_store: Vec<(&str, &str)> = variables.iter().copied().chain([STORE]).collect();
        let (report, store_messages, stats) = blastn(queries, subjects, options, &with_store);
        assert!(
            report == expected,
            "{options:?}: the packed subjects change the report"
        );
        assert_eq!(store_messages, messages, "{options:?}");
        assert_eq!(counter(&stats, "subjects"), 6, "{options:?}");
        assert!(counter(&stats, "stored_bytes") > 0, "{options:?}");
        // the traceback of the four subjects with hits, and the rows of the pairwise report
        assert!(counter(&stats, "blastna_made") >= 4, "{options:?}: {stats}");
        if options[1] == "0" {
            assert!(counter(&stats, "rows_made") >= 4, "{options:?}: {stats}");
        } else {
            assert_eq!(counter(&stats, "rows_made"), 0, "{options:?}: {stats}");
        }
    }
}

// No NCBI counterpart: the switch is on by default (an empty value, as unset): the six subjects
// are held packed; with `LOSAT_X_DEFAULT=0` and no value they are not; the report and the
// messages are the same.
#[test]
fn the_subjects_are_held_packed_by_default() {
    let (queries, subjects) = inputs();
    let options = ["-outfmt", "6", "-num_threads", "3"];
    let default = ("LOSAT_X_BNSTORE", "");
    let (report, messages, stats) = blastn(queries, subjects, &options, &[default]);
    assert_eq!(counter(&stats, "subjects"), 6, "{stats}");
    let (old_report, old_messages, old_stats) = blastn(
        queries,
        subjects,
        &options,
        &[default, ("LOSAT_X_DEFAULT", "0")],
    );
    assert_eq!(counter(&old_stats, "subjects"), 0, "{old_stats}");
    assert!(report == old_report, "the default changes the report");
    assert_eq!(messages, old_messages);
}

// No NCBI counterpart: where the switch acts. A subject range (`-subject_loc`) keeps the
// letters (the records are cut before the search), so the variable changes nothing there and no
// subject is held packed. With `LOSAT_X_DEFAULT=0` (the reference loops) the subjects are still
// held packed and the report is the same.
#[test]
fn a_subject_range_keeps_the_letters_and_the_reference_loops_read_the_packed_subjects() {
    let (queries, subjects) = inputs();
    let long = long_subjects(subjects);
    let ranged = ["-outfmt", "0", "-subject_loc", "100-5000", "-lcase_masking"];
    let (expected, messages, _) = blastn(queries, &long, &ranged, &[]);
    let (report, store_messages, stats) = blastn(queries, &long, &ranged, &[STORE]);
    assert!(report == expected);
    assert_eq!(store_messages, messages);
    assert_eq!(counter(&stats, "subjects"), 0);
    assert_eq!(counter(&stats, "blastna_made"), 0);
    assert_eq!(counter(&stats, "rows_made"), 0);

    let options = ["-outfmt", "0", "-num_threads", "1"];
    let (expected, messages, _) = blastn(queries, subjects, &options, &[]);
    let previous = ("LOSAT_X_DEFAULT", "0");
    let (report, store_messages, stats) = blastn(queries, subjects, &options, &[previous, STORE]);
    assert!(report == expected);
    assert_eq!(store_messages, messages);
    assert_eq!(counter(&stats, "subjects"), 6);
    assert!(counter(&stats, "blastna_made") >= 4);
    assert!(counter(&stats, "rows_made") >= 4);
}
