//! `losat makeblastdb -mask_data`, `-mask_id` and `-mask_desc` through the command line: the
//! words go through the parser and `makeblastdb::run_with` (the entry point of the command and of
//! the library's `run_makeblastdb`), and the exit status, the log on the standard output and the
//! text on the standard error stream are compared with what NCBI BLAST+ 2.17.0's makeblastdb wrote
//! for the same words.
//!
//! * `the_command_line_matches_makeblastdb_on_the_cases_of_the_writer_fixtures`: the 36 cases of
//!   `tests/seqdb/mask_writer` (`expected/<case>.{rc,stderr,stdout}`), and the database each one
//!   writes (files, algorithms and ranges as NCBI's `blastdbcmd` reads them back).
//! * `the_command_line_matches_makeblastdb_on_every_mask_parse_file`: the 128 runs of
//!   `tests/seqdb/mask_parse` (`expected.tsv`, `expected_stdout.tsv`).
//! * `the_command_line_reaches_the_sequences_that_ncbi_reaches_with_parse_seqids`: a few of the
//!   cases of `tests/seqdb/mask_ids`, with `-parse_seqids`, in the default format version 5.
//! * `the_command_line_reads_the_options_as_makeblastdb_does`: a few answers of `x_ProcessMaskData`
//!   to its options (the rest are in `tests_e2e`).
//! * `errors_of_makeblastdb_come_before_the_refusal_of_masks_ordered_by_address`: the 5 cases of
//!   `tests/seqdb/mask_order` in format versions 5 and 4 (`expected/<case>.v<N>.{rc,stderr,stdout,files}`).
//!
//! The mask and input files are named by absolute paths (a test cannot change the directory of the
//! process); the directory is cut off the messages, which then are those of a run in the directory
//! that holds the files, as the expected texts were made.

use std::collections::BTreeMap;
use std::path::Path;

use super::tests::column_date;
use super::tests_e2e::{
    assert_database, assert_reads_back, assert_same_files, classify, option_root, ranges_text,
    read, read_expected, stored_title, unescape, volume_bases, writer_cases, Class,
    DESTRUCTOR_MESSAGE,
};
use super::tests_ids::{masked_rows, rows, Row};
use crate::seqdb::{SeqDb, SeqType};
use crate::writedb::makeblastdb::run_with;
use crate::writedb::tests::{
    files_named, log_with_markers, parse_makedb_words, replace_bytes, seqdb_fixtures, stored_date,
    time_of, TempDir, LINUX,
};

/// What a run of the command left on its two streams, with the directory of the input files cut
/// off and the log's time, path and seconds masked.
struct Output {
    status: i32,
    stdout: String,
    stderr: String,
}

/// Runs `losat makeblastdb` with `words` (after the program name) at the clock time `date`.
/// `dirs` are the directories of the input files, whose names are cut off the streams.
fn run_words(words: &[String], date: &str, dirs: &[&Path]) -> Output {
    let args = parse_makedb_words(words);
    let now = time_of(date);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let status = run_with(
        &args,
        None,
        Some(&now),
        LINUX,
        &mut std::io::empty(),
        &mut out,
        &mut err,
    );
    let cut = |bytes: Vec<u8>| {
        let mut bytes = bytes;
        for dir in dirs {
            let prefix = format!("{}/", dir.display());
            bytes = replace_bytes(&bytes, prefix.as_bytes(), b"");
        }
        bytes
    };
    Output {
        status,
        stdout: String::from_utf8_lossy(&log_with_markers(&cut(out))).into_owned(),
        stderr: String::from_utf8_lossy(&cut(err)).into_owned(),
    }
}

/// `-in` of one or several files: quoted, as a path may hold a space.
fn in_words(files: &[std::path::PathBuf]) -> String {
    files
        .iter()
        .map(|file| format!("\"{}\"", file.display()))
        .collect::<Vec<_>>()
        .join(" ")
}

/// `-mask_data` of the files named in `list`, in the directory `root`.
fn mask_words(root: &Path, list: &str) -> String {
    list.split(',')
        .map(|name| format!("{}/{name}", root.display()))
        .collect::<Vec<_>>()
        .join(",")
}

/// The ranges of every sequence of the database as NCBI's `blastdbcmd` read them from the
/// fixture (`expected/<case>.masks.tsv`).
fn assert_ranges_read_back(base: &str, protein: bool, expected_dir: &Path, name: &str) {
    let seqtype = if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let db = SeqDb::open(base, seqtype).unwrap();
    let rows = read(expected_dir, &format!("{name}.masks.tsv"));
    assert!(!rows.is_empty(), "{name}");
    for line in rows.lines() {
        let mut fields = line.split('\t');
        let algo: i32 = fields.next().unwrap().parse().unwrap();
        let oid: u32 = fields.next().unwrap().parse().unwrap();
        let expected = fields.next().unwrap();
        assert_eq!(
            ranges_text(
                &db.mask_ranges_at(db.locate_oid(oid).unwrap(), algo)
                    .unwrap()
            ),
            expected,
            "{name}: algorithm {algo}, sequence {oid}"
        );
    }
}

#[test]
fn the_command_line_matches_makeblastdb_on_the_cases_of_the_writer_fixtures() {
    let dir = seqdb_fixtures().join("mask_writer");
    let seqdb = seqdb_fixtures();
    let expected_dir = dir.join("expected");
    let cases = writer_cases();
    assert_eq!(cases.len(), 36);
    let (mut built, mut failed) = (0, 0);
    for case in &cases {
        let name = &case.name;
        let rc: i32 = read(&expected_dir, &format!("{name}.rc"))
            .trim()
            .parse()
            .unwrap();
        let stderr = read(&expected_dir, &format!("{name}.stderr"));
        let stdout = read(&expected_dir, &format!("{name}.stdout"));
        let tmp = TempDir::new(name);
        let base = tmp.base(name);

        let mut words: Vec<String> = [
            "-in",
            &in_words(&case.inputs),
            "-dbtype",
            if case.protein { "prot" } else { "nucl" },
            "-out",
            &base,
            "-blastdb_version",
            "4",
            "-title",
            &case.title,
        ]
        .map(str::to_string)
        .to_vec();
        if let Some(list) = &case.mask_data {
            words.extend(["-mask_data".to_string(), mask_words(&case.root, list)]);
        }
        if let Some(ids) = &case.mask_id {
            words.extend(["-mask_id".to_string(), ids.clone()]);
        }
        if let Some(descs) = &case.mask_desc {
            words.extend(["-mask_desc".to_string(), descs.clone()]);
        }
        if case.max_file_size != 3_000_000_000 {
            words.extend([
                "-max_file_sz".to_string(),
                format!("{}B", case.max_file_size),
            ]);
        }

        // The date of the new database is the time of the run: the one NCBI's files hold.
        let letter = if case.protein { 'p' } else { 'n' };
        let seqtype = if case.protein {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        let fixture_bases = volume_bases(&dir.join("db"), name, case.protein);
        let first_index = dir
            .join("db")
            .join(format!("{}.{letter}in", fixture_bases[0]));
        let date = if rc == 0 {
            stored_date(&first_index, seqtype)
        } else {
            "Oct 20, 2026  1:00 PM".to_string()
        };

        let run = run_words(&words, &date, &[&dir, &seqdb]);
        assert_eq!(run.status, rc, "{name}: exit status\n{}", run.stderr);
        assert_eq!(run.stderr, stderr, "{name}: standard error");
        assert_eq!(run.stdout, stdout, "{name}: standard output");
        if rc == 0 {
            built += 1;
            // makeblastdb also writes the metadata file, which the fixtures do not keep.
            for extension in ["njs", "pjs"] {
                let _ = std::fs::remove_file(tmp.0.join(format!("{name}.{extension}")));
            }
            assert_eq!(
                stored_title(
                    &tmp.0.join(format!("{}.{letter}in", fixture_bases[0])),
                    seqtype
                ),
                case.title,
                "{name}: the title"
            );
            // The name of an algorithm made from a program other than dust, seg and
            // windowmasker is the name of its mask file as given, which is a path here: the
            // table of algorithms and the column files differ from NCBI's by that name alone.
            // Those cases compare the ranges only (the library test `tests_e2e` has the files).
            let named_by_path =
                std::fs::read_to_string(expected_dir.join(format!("{name}.algorithms.txt")))
                    .is_ok_and(|table| table.contains(".asnt"));
            if named_by_path {
                assert_eq!(name, "ids_mix");
                assert_ranges_read_back(&base, case.protein, &expected_dir, name);
                continue;
            }
            if expected_dir.join(format!("{name}.masks.tsv")).exists() {
                assert_reads_back(&base, case.protein, &expected_dir, name);
            }
            // The files, byte for byte. The index date is the run's clock, as in NCBI's files; the
            // column files store the time of their creation to the second, which the command line
            // cannot be given, so that date is set to the one in NCBI's file first.
            for base in &fixture_bases {
                let written = tmp.0.join(format!("{base}.{letter}aa"));
                if let (Some(want), Some(have)) = (
                    column_date(&dir.join("db").join(format!("{base}.{letter}aa"))),
                    column_date(&written),
                ) {
                    let bytes = std::fs::read(&written).unwrap();
                    std::fs::write(
                        &written,
                        replace_bytes(&bytes, have.as_bytes(), want.as_bytes()),
                    )
                    .unwrap();
                }
            }
            assert_same_files(&tmp.0, name, &dir.join("db"), name, case.version);
        } else {
            failed += 1;
            assert_eq!(
                files_named(&tmp.0, &format!("{name}.")),
                Vec::<String>::new(),
                "{name}: files left behind"
            );
        }
    }
    assert_eq!((built, failed), (25, 11), "built and failed cases");
}

#[test]
fn the_command_line_matches_makeblastdb_on_every_mask_parse_file() {
    let seqdb = seqdb_fixtures();
    let parse_dir = seqdb.join("mask_parse");
    let expected = read_expected(&parse_dir);
    let logs: BTreeMap<String, String> = read(&parse_dir, "expected_stdout.tsv")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let (file, log) = line.split_once('\t').unwrap();
            (file.to_string(), unescape(log))
        })
        .collect();
    assert_eq!(logs.len(), expected.len());
    let mut classes: BTreeMap<String, usize> = BTreeMap::new();
    let mut compared = 0;
    for (file, ncbi) in &expected {
        let class = classify(file, ncbi);
        *classes.entry(format!("{class:?}")).or_default() += 1;
        let protein = ncbi.dbtype == "prot";
        let fasta = if parse_dir.join(&ncbi.fasta).exists() {
            parse_dir.join(&ncbi.fasta)
        } else {
            seqdb.join(&ncbi.fasta)
        };
        let tmp = TempDir::new("mask-parse-cli");
        let base = tmp.base("db");
        // The default format version (5).
        let words: Vec<String> = [
            "-in",
            &in_words(&[fasta]),
            "-dbtype",
            &ncbi.dbtype,
            "-out",
            &base,
            "-title",
            &ncbi.fasta,
            "-mask_data",
            &mask_words(&parse_dir, file),
        ]
        .map(str::to_string)
        .to_vec();
        let run = run_words(&words, "Oct 20, 2026  1:00 PM", &[&parse_dir]);
        let context = format!("{file}: {}", run.stderr);
        // A failed build of format version 5 may leave the three tax-id files, as makeblastdb
        // 2.17.0 does (`db.not`, `db.ntf`, `db.nto` after "Masked data offsets out of bounds"):
        // nothing else is left.
        let no_files = || {
            let left = files_named(&tmp.0, "db.");
            assert!(
                left.iter()
                    .all(|f| ["db.not", "db.ntf", "db.nto"].contains(&f.as_str())),
                "{file}: files left behind: {left:?}"
            );
        };
        match class {
            Class::UnknownEncoding | Class::NoDatabase | Class::Equal => {
                assert_eq!(run.status, ncbi.rc, "{context}");
                assert_eq!(run.stderr, ncbi.stderr, "{file}: standard error");
                assert_eq!(run.stdout, logs[file], "{file}: standard output");
                if run.status == 0 && class == Class::Equal {
                    compared += assert_database(&base, protein, file, ncbi);
                } else {
                    no_files();
                }
            }
            Class::Serial | Class::OutOfMemory => {
                assert_eq!(run.status, 255, "{context}");
                assert!(
                    run.stderr
                        .starts_with(&format!("Error: mask data in '{file}': ")),
                    "{context}"
                );
                assert_eq!(
                    run.stderr.ends_with(DESTRUCTOR_MESSAGE),
                    ncbi.stderr.ends_with(DESTRUCTOR_MESSAGE),
                    "{file}: the message of the destructor"
                );
                assert_eq!(run.stdout, logs[file], "{file}: standard output");
                no_files();
            }
            Class::Unsupported => {
                assert_eq!(run.status, 1, "{context}");
                assert!(
                    run.stderr
                        .starts_with("losat makeblastdb: mask data that LOSAT does not read: "),
                    "{context}"
                );
                no_files();
            }
        }
    }
    assert!(compared > 300, "{compared} rows compared");
    assert_eq!(classes["Equal"], 109);
    assert_eq!(classes.values().sum::<usize>(), 162);
}

#[test]
fn a_sample_of_the_mask_ids_cases_runs_through_the_command_line() {
    let dir = seqdb_fixtures().join("mask_ids");
    let case_dir = dir.join("case");
    // The first case of each fate, and the masks of the real tools (`A_*`, `B_*`, `N_*`).
    let all = rows();
    let mut chosen: Vec<&Row> = Vec::new();
    for fate in ["refuse", "NOMATCH", "match"] {
        let found = all.iter().find(|row| {
            !row.loss.starts_with("env:")
                && if fate == "match" {
                    row.loss == "match" && row.expected != "NOMATCH"
                } else if fate == "NOMATCH" {
                    row.loss == "match" && row.expected == "NOMATCH"
                } else {
                    row.loss == fate
                }
        });
        chosen.extend(found);
    }
    chosen.extend(all.iter().filter(|row| {
        ["A_dust_parsed", "B_prot_seg_parsed", "N_dust_plain"].contains(&row.name.as_str())
    }));
    assert_eq!(chosen.len(), 6);
    for row in chosen {
        let tmp = TempDir::new("mask-ids-cli");
        let base = tmp.base("db");
        let mut words: Vec<String> = [
            "-in",
            &in_words(&[dir.join("fasta").join(&row.fasta)]),
            "-dbtype",
            if row.protein { "prot" } else { "nucl" },
            "-out",
            &base,
            "-title",
            &row.fasta,
            "-mask_data",
            &mask_words(&case_dir, &row.files),
        ]
        .map(str::to_string)
        .to_vec();
        if row.parse {
            words.push("-parse_seqids".to_string());
        }
        let run = run_words(&words, "Oct 20, 2026  1:00 PM", &[&case_dir]);
        let context = format!("{}: {}", row.name, run.stderr);
        if row.loss == "refuse" {
            // LOSAT's refusal alone: status 1, no destructor message of NCBI's after it.
            assert_eq!(run.status, 1, "{context}");
            assert!(run.stderr.starts_with("losat makeblastdb: "), "{context}");
            assert!(run.stderr.contains("memory addresses"), "{context}");
            assert_eq!(run.stderr.lines().count(), 1, "{context}");
        } else if row.expected == "NOMATCH" {
            assert_eq!(run.status, 0, "{context}");
            assert_eq!(run.stderr, DESTRUCTOR_MESSAGE, "{context}");
        } else {
            assert_eq!(run.status, 0, "{context}");
            assert_eq!(run.stderr, "", "{context}");
            assert_eq!(
                masked_rows(&base, row.protein, &row.algorithms),
                row.expected,
                "{context}"
            );
        }
    }
}

/// The runs of `tests/seqdb/mask_order` (review 2 of M4a-1, probes t1, t1b, t3): masks of
/// `gi|5` and `gb|AY000005.1|` (one packed number) that NCBI writes in the order of the addresses
/// of their handles, and an error of makeblastdb that does not depend on that order later in the
/// run (a range of the same sequence or of a later algorithm past the end, a later sequence, a
/// later input file without sequences). NCBI's order never reaches a file there, so LOSAT gives
/// makeblastdb's exit status, messages, log and files. Where nothing else fails (`tie`), LOSAT
/// refuses at the close, after the log of the whole build, and leaves the files of NCBI's failed
/// build of the same input (`oob_later_sequence`).
#[test]
fn errors_of_makeblastdb_come_before_the_refusal_of_masks_ordered_by_address() {
    let dir = seqdb_fixtures().join("mask_order");
    let expected_dir = dir.join("expected");
    let (mut compared, mut refused) = (0, 0);
    for line in read(&dir, "cases.tsv").lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (name, inputs, masks) = (fields[0], fields[1], fields[2]);
        for version in ["5", "4"] {
            let stem = format!("{name}.v{version}");
            let expected = |extension: &str| read(&expected_dir, &format!("{stem}.{extension}"));
            let tmp = TempDir::new("mask-order");
            let in_files: Vec<_> = inputs.split(' ').map(|file| dir.join(file)).collect();
            let words: Vec<String> = [
                "-in",
                &in_words(&in_files),
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-title",
                "in.fna",
                "-mask_data",
                &mask_words(&dir, masks),
                "-out",
                &tmp.base("db"),
                "-blastdb_version",
                version,
            ]
            .map(str::to_string)
            .to_vec();
            let run = run_words(&words, "Oct 20, 2026  1:00 PM", &[&dir]);
            let left: String = files_named(&tmp.0, "db.")
                .iter()
                .map(|file| format!("{file}\n"))
                .collect();
            let context = format!("{stem}: {}", run.stderr);
            assert_eq!(run.stdout, expected("stdout"), "{stem}: standard output");
            if name == "tie" {
                // NCBI builds it, in the order of the addresses.
                assert_eq!(expected("rc"), "0\n");
                assert_eq!(run.status, 1, "{context}");
                assert_eq!(
                    run.stderr,
                    "losat makeblastdb: the masks of the Seq-ids gi|5 and gb|AY000005.1| of the \
                     sequence gi|5 are written in an order that NCBI takes from the memory \
                     addresses of their handles (equal packed numbers); LOSAT does not write this \
                     database\n",
                    "{stem}"
                );
                assert_eq!(
                    left,
                    read(
                        &expected_dir,
                        &format!("oob_later_sequence.v{version}.files")
                    ),
                    "{stem}: files left"
                );
                refused += 1;
            } else {
                assert_eq!(format!("{}\n", run.status), expected("rc"), "{context}");
                assert_eq!(run.stderr, expected("stderr"), "{stem}: standard error");
                assert_eq!(left, expected("files"), "{stem}: files left");
                compared += 1;
            }
        }
    }
    assert_eq!((compared, refused), (8, 2));
}

#[test]
fn the_command_line_reads_the_options_as_makeblastdb_does() {
    let tmp = option_root();
    let fasta = seqdb_fixtures().join("mask_nucl.fna");
    let root = tmp.0.clone();
    let missing = |name: &str| {
        format!("Error: [makeblastdb] Ignoring mask file '{name}' as it does not exist.\n")
    };
    // (-mask_data, -mask_id, -mask_desc, status, standard error, the "Mask file:" lines of the log):
    // the answers of makeblastdb 2.17.0 recorded in `tests_e2e` (`opt.py` of the M4 brief C).
    #[allow(clippy::type_complexity)]
    let cases: Vec<(&str, Option<&str>, Option<&str>, i32, String, Vec<&str>)> = vec![
        (
            "",
            None,
            None,
            1,
            "BLAST options error: mask_data option found, but no files were specified.\n".into(),
            vec![],
        ),
        (
            "nonexist.asnt",
            None,
            None,
            0,
            missing("nonexist.asnt"),
            vec![],
        ),
        (
            "nonexist.asnt,dust.asnt",
            None,
            None,
            0,
            missing("nonexist.asnt"),
            vec!["dust.asnt"],
        ),
        (
            "dust.asnt",
            Some("5"),
            None,
            1,
            "BLAST options error: mask_id can not be numeric.\n".into(),
            vec![],
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a"),
            None,
            1,
            "BLAST options error: the size of mask_id does not match that of mask_data.\n".into(),
            vec![],
        ),
        (
            "dust.asnt",
            Some("a"),
            Some("x,y"),
            1,
            "BLAST options error: the size of mask_desc does not match that of mask_data.\n".into(),
            vec![],
        ),
        (
            "empty.txt",
            None,
            None,
            1,
            "BLAST options error: Unknown encoding for mask data in 'empty.txt'\n".into(),
            vec![],
        ),
        (
            "dust.asnt",
            Some("alpha"),
            Some("first"),
            0,
            String::new(),
            vec!["dust.asnt"],
        ),
    ];
    for (mask_data, mask_id, mask_desc, status, stderr, files) in cases {
        let tmp_out = TempDir::new("mask-options-cli");
        let mut words: Vec<String> = [
            "-in",
            &in_words(std::slice::from_ref(&fasta)),
            "-dbtype",
            "nucl",
            "-out",
            &tmp_out.base("db"),
            "-title",
            "mask_nucl.fna",
        ]
        .map(str::to_string)
        .to_vec();
        let list = if mask_data.is_empty() {
            String::new()
        } else {
            mask_words(&root, mask_data)
        };
        words.extend(["-mask_data".to_string(), list]);
        if let Some(ids) = mask_id {
            words.extend(["-mask_id".to_string(), ids.to_string()]);
        }
        if let Some(descs) = mask_desc {
            words.extend(["-mask_desc".to_string(), descs.to_string()]);
        }
        let run = run_words(&words, "Oct 20, 2026  1:00 PM", &[&root]);
        let context = format!("{mask_data:?} {mask_id:?} {mask_desc:?}");
        assert_eq!(run.status, status, "{context}: {}", run.stderr);
        assert_eq!(run.stderr, stderr, "{context}: standard error");
        let logged: Vec<&str> = run
            .stdout
            .lines()
            .filter_map(|line| line.strip_prefix("Mask file: "))
            .collect();
        assert_eq!(logged, files, "{context}: the mask files in the log");
    }
}
