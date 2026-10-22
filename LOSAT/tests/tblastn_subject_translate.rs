//! The six translated frames of a TBLASTN subject read from tables (`LOSAT_X_TBNXLATE`,
//! `src/algorithm/tblastn/x_subject_translate.rs`), on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1071-1074,1080-1090
//! ```c
//!    } else {
//!       translation_table = s_BlastGetTranslationTable(genetic_code, FALSE);
//!       translation_table_rc = s_BlastGetTranslationTable(genetic_code, TRUE);
//!    }
//! ...
//!    for (context = 0; context < NUM_FRAMES; ++context) {
//!       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//!       if (encoding == eBlastEncodingNcbi2na) {
//!          if (frame > 0) {
//!             length =
//!                BLAST_TranslateCompressedSequence(translation_table,
//!                   nucl_length, nucl_seq, frame, translation_buffer+offset);
//!          } else {
//!             length =
//!                BLAST_TranslateCompressedSequence(translation_table_rc,
//!                   nucl_length, nucl_seq, frame, translation_buffer+offset);
//! ```
//! NCBI translates the compressed subject with one table read per codon. LOSAT's reference path
//! reaches the same residues through letters; with `LOSAT_X_TBNXLATE` (a switch of the strict set:
//! on by default, off with `LOSAT_X_TBNXLATE=0`, and off with `LOSAT_X_DEFAULT=0` unless its
//! variable is set) the frames are read from the tables, and with `LOSAT_X_TBNXLATESHADOW` both
//! translations are made and compared for every subject plan of the unit path (`LOSAT_X_TBNPAR`,
//! on by default). The switch must change nothing in the report, so the tests compare the report
//! with the one of a run with `LOSAT_X_TBNXLATE=0`, and read the `[X_STATS] TBNXLATE` and
//! `[X_SHADOW] TBNXLATE` lines to check which translation ran.
//!
//! The subjects are pieces of `tests/fasta/MelaMJNV.fasta` with what the two translations
//! treat differently: lowercase runs (soft masking with `-lcase_masking`), every ambiguous
//! IUPAC letter and runs of `N` (one draw each from the generator seeded with the subject's
//! length), one subject too short for a word (it gets no plan) and two of a few codons. The
//! queries are proteins of `tests/fasta/PemoMJNVA.faa`: ten of them (one query batch), or all
//! 115 (four batches, whose plans go through the cells that keep them for the later batches,
//! `LOSAT_X_TBNSSIDE`).

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
        // One directory per process: nextest runs each test in a process of its own, and a
        // process that wrote these files again while another one's search read them cut its input.
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join("tblastn_subject_translate")
            .join(std::process::id().to_string());
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
            ambiguous[position] = b"NRYKMSWBDHVnrykmswbdhv"[count % 22];
        }
        lowercase(&mut ambiguous, 2_000, 2_300);
        let mut runs = piece(3);
        runs[3_000..3_040].fill(b'N');
        runs[5_000..5_007].fill(b'n');
        lowercase(&mut runs, 900, 1_200);
        let mut subjects = String::new();
        for (title, letters) in [
            ("plain", plain.as_slice()),
            ("lower", lower.as_slice()),
            ("ambiguous", ambiguous.as_slice()),
            ("runs", runs.as_slice()),
            ("too_short_for_a_word", b"ACGTN".as_slice()),
            ("one_word", b"ATGNCGTRAAC".as_slice()),
            ("few_codons", b"ATGAAACGTYTGGCNnnAAATAG".as_slice()),
        ] {
            subjects.push_str(&format!(
                ">{title}\n{}\n",
                String::from_utf8(letters.to_vec()).expect("ASCII letters")
            ));
        }
        let files = (
            dir.join("ten.faa"),
            dir.join("all.faa"),
            dir.join("subjects.fna"),
        );
        std::fs::write(&files.0, ten).expect("query file");
        std::fs::write(&files.1, all).expect("query file");
        std::fs::write(&files.2, subjects).expect("subject file");
        files
    })
}

/// One TBLASTN run (`-outfmt 6` on standard output) with `options` and with `variables` and
/// `LOSAT_X_STATS` as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build
/// with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so
/// the caller's value is kept, or a file of the tests is named. Returns the report and the
/// standard error.
fn tblastn(query: &Path, options: &[&str], variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(&inputs().2)
        .args(["-outfmt", "6"])
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
    (output.stdout, stderr)
}

/// The number after `field=` in the one line of the standard error that starts with `line`.
fn counter(stderr: &str, line: &str, field: &str) -> u64 {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|text| text.starts_with(line))
        .collect();
    assert_eq!(lines.len(), 1, "one `{line}` line:\n{stderr}");
    lines[0]
        .split_whitespace()
        .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no {field} in {}", lines[0]))
        .parse()
        .expect("a number")
}

/// The subjects named in the second column of a report.
fn hit_subjects(report: &[u8]) -> std::collections::BTreeSet<String> {
    String::from_utf8_lossy(report)
        .lines()
        .filter_map(|row| row.split('\t').nth(1).map(str::to_string))
        .collect()
}

const TABLES: (&str, &str) = ("LOSAT_X_TBNXLATE", "1");
const REFERENCE: (&str, &str) = ("LOSAT_X_TBNXLATE", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBNXLATESHADOW", "1");
/// The subjects keep their letters: with the packed subjects of the default (`LOSAT_X_TBNSTORE`,
/// `tests/tblastn_subject_store.rs`) the frames are translated from the packed bases whatever
/// `LOSAT_X_TBNXLATE` says.
const LETTERS: (&str, &str) = ("LOSAT_X_TBNSTORE", "0");

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
// Every search below gives one report with `LOSAT_X_TBNXLATE=0`, with the frames from the tables
// (no variable: the default) and in the shadow mode: the default options and `-lcase_masking` at
// one and three threads, another genetic code, and four query batches. The run with
// `LOSAT_X_TBNXLATE=0` makes no plan from the tables; the default run makes one for each of the
// six subjects that are long enough for a word (more if a batch makes its own); the shadow mode
// compares them. Every run keeps the subjects' letters (`LETTERS`).
#[test]
fn the_report_is_the_same_with_the_frames_from_the_tables_and_in_the_shadow_mode() {
    let (ten, all, _) = inputs();
    let searches: [(&Path, &[&str]); 6] = [
        (ten, &["-num_threads", "1"]),
        (ten, &["-num_threads", "3"]),
        (ten, &["-num_threads", "1", "-lcase_masking"]),
        (ten, &["-num_threads", "3", "-lcase_masking"]),
        (ten, &["-num_threads", "1", "-db_gencode", "4"]),
        (
            all,
            &["-num_threads", "3", "-lcase_masking", "-db_gencode", "11"],
        ),
    ];
    for (query, options) in searches {
        let (expected, stderr) = tblastn(query, options, &[REFERENCE, LETTERS]);
        assert_eq!(
            counter(&stderr, "[X_STATS] TBNXLATE", "table_plans"),
            0,
            "{options:?}"
        );
        let hits = hit_subjects(&expected);
        for subject in ["plain", "lower", "ambiguous", "runs"] {
            assert!(hits.contains(subject), "{options:?}: no hit on {subject}");
        }

        let (report, stderr) = tblastn(query, options, &[LETTERS]);
        assert!(
            report == expected,
            "{options:?}: the tables change the report"
        );
        let table_plans = counter(&stderr, "[X_STATS] TBNXLATE", "table_plans");
        assert!(table_plans >= 6, "{options:?}: {table_plans} plans");
        assert!(!stderr.contains("[X_SHADOW] TBNXLATE"));

        let (report, stderr) = tblastn(query, options, &[SHADOW, LETTERS]);
        assert!(
            report == expected,
            "{options:?}: the shadow mode changes the report"
        );
        let compared = counter(&stderr, "[X_SHADOW] TBNXLATE", "plans_compared");
        assert!(compared >= 6, "{options:?}: {compared} plans compared");
    }
}

// No NCBI counterpart: where the switch acts. The frames are read from the tables in the
// subject plans of the unit path (`LOSAT_X_TBNPAR`, on by default). With `LOSAT_X_DEFAULT=0`
// the reference loop searches each subject and translates it itself, so the variable changes
// nothing there: no plan is made, the shadow mode says that it compared nothing, and the report
// is the same. With `LOSAT_X_TBNPAR` set as well the unit path runs again, and the switch is off
// unless its variable is set, as before it joined the set.
#[test]
fn the_switch_acts_on_the_subject_plans_of_the_unit_path_only() {
    let (ten, _, _) = inputs();
    let options = ["-num_threads", "1"];
    let (expected, _) = tblastn(ten, &options, &[]);
    let previous = ("LOSAT_X_DEFAULT", "0");
    for variables in [&[previous][..], &[previous, TABLES], &[previous, SHADOW]] {
        let (report, stderr) = tblastn(ten, &options, variables);
        assert!(report == expected, "{variables:?}");
        assert_eq!(
            counter(&stderr, "[X_STATS] TBNXLATE", "table_plans"),
            0,
            "{variables:?}"
        );
        if variables.contains(&SHADOW) {
            assert_eq!(counter(&stderr, "[X_SHADOW] TBNXLATE", "plans_compared"), 0);
            assert!(stderr.contains("nothing compared"), "{stderr}");
        }
    }
    let unit_path = ("LOSAT_X_TBNPAR", "1");
    let (report, stderr) = tblastn(ten, &options, &[previous, unit_path]);
    assert!(report == expected);
    assert_eq!(counter(&stderr, "[X_STATS] TBNXLATE", "table_plans"), 0);
    let (report, stderr) = tblastn(ten, &options, &[previous, unit_path, TABLES]);
    assert!(report == expected);
    assert!(counter(&stderr, "[X_STATS] TBNXLATE", "table_plans") >= 6);
}
