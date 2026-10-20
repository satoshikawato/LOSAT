//! End-to-end tests of `-mask_data`: makeblastdb's flow from the mask files to the database,
//! through the library (`process_mask_data`, `BuildDatabase::{add_fasta, end_build}`,
//! `post_destruction_diagnostics`), in the order the command line will call them (see the
//! module documentation of `process`). The tests run no NCBI program and read the files that
//! NCBI's tools and the fixture scripts wrote.
//!
//! * `the_writer_rebuilds_makeblastdb_s_masked_databases_byte_for_byte`: the databases of
//!   `LOSAT/tests/seqdb/db` that `make_mask_fixtures.py` built in format version 4, from the
//!   mask files NCBI's `dustmasker` / `segmasker` wrote (`mask_dust.asnb`, `mask_seg.asnb`, binary)
//!   and the hand-written `mask_wm*.asnt`; every file of the database is compared.
//! * `the_writer_matches_makeblastdb_on_the_cases_of_the_writer_fixtures`: the 36 cases of
//!   `make_mask_writer_fixtures.py` (`tests/seqdb/mask_writer`): exit status, standard error, the
//!   files, and what a reader finds in the database.
//! * `every_mask_parse_file_runs_end_to_end`: every file of `tests/seqdb/mask_parse/expected.tsv`
//!   (`make_mask_parse_fixtures.py`) with the FASTA it was run with.
//! * `process_mask_data_follows_makeblastdb_on_its_options`: the option handling of
//!   `x_ProcessMaskData` (lists, empty tokens, numeric ids, files that do not exist), against the
//!   answers of makeblastdb 2.17.0.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use super::process::{process_mask_data, MaskDataError, MaskDataOptions};
use super::tests::{bioseq, column_date};
use crate::seqdb::{MaskRanges, SeqDb, SeqType};
use crate::writedb::build_db::BuildDatabase;
use crate::writedb::fasta_source::Warn;
use crate::writedb::makeblastdb::MakedbError;
use crate::writedb::mask_info::MaskFileError;
use crate::writedb::tests::{files_named, seqdb_fixtures, stored_date, time_of, TempDir};
use crate::writedb::{WriteDb, WriteDbError};

/// What makeblastdb is given: the FASTA files, and the options `x_ProcessMaskData` reads.
struct Case {
    name: String,
    protein: bool,
    /// The FASTA files in the order of `-in`.
    inputs: Vec<PathBuf>,
    /// `-mask_data`, `-mask_id`, `-mask_desc`: the option values as the command line has them.
    mask_data: Option<String>,
    mask_id: Option<String>,
    mask_desc: Option<String>,
    /// The directory the mask file names are looked up in.
    root: PathBuf,
    max_file_size: u64,
    title: String,
}

impl Case {
    fn options(&self) -> MaskDataOptions<'_> {
        MaskDataOptions {
            mask_data: self.mask_data.as_deref(),
            mask_id: self.mask_id.as_deref(),
            mask_desc: self.mask_desc.as_deref(),
            gi_mask_name: None,
            root: Some(&self.root),
        }
    }
}

/// The result of a run: makeblastdb's exit status, the text on its standard error stream and
/// the log (its standard output).
struct Run {
    status: i32,
    stderr: Vec<u8>,
    log: String,
}

impl Run {
    fn stderr_text(&self) -> String {
        String::from_utf8_lossy(&self.stderr).into_owned()
    }

    /// The names after "Mask file: " in the log.
    fn mask_files(&self) -> Vec<&str> {
        self.log
            .lines()
            .filter_map(|line| line.strip_prefix("Mask file: "))
            .collect()
    }
}

/// How a step of `x_BuildDatabase` ends it.
enum Stop {
    /// `CATCH_ALL` prints the error and sets the status.
    Error(MakedbError),
    /// The exception that `CATCH_ALL` ignores: status 0, nothing printed.
    Silent,
}

/// `x_BuildDatabase` from the construction of the `CBuildDatabase` on, as `LOSAT makedb` will
/// call it with `-mask_data`: the maximum file size, the mask files, the input, the end of the
/// build (the log head is written by `BuildDatabase::new`).
fn build(db: &mut BuildDatabase<'_>, case: &Case, warn: Warn<'_>) -> Result<(), Stop> {
    db.set_max_file_size(case.max_file_size);
    process_mask_data(db, &case.options(), &mut *warn)
        .map_err(|error| error.into_makedb_error().map_or(Stop::Silent, Stop::Error))?;
    for input in &case.inputs {
        let file = std::fs::File::open(input).unwrap();
        db.add_fasta(file, &mut *warn)
            .map_err(|error| Stop::Error(error.into()))?;
    }
    db.end_build(false, &mut *warn)
        .map_err(|error| Stop::Error(error.into()))?;
    Ok(())
}

/// A whole run: `build`, then what `Run` of makeblastdb does with its outcome (`CATCH_ALL`), then
/// the destructor of the `CBuildDatabase`, in that order. `column_dates` are the dates of the
/// column files of the volumes to reproduce (`None`: the clock).
fn run_case(case: &Case, out_base: &str, index_date: &str, column_dates: &[Option<String>]) -> Run {
    let mut log = Vec::new();
    let mut messages = Vec::new();
    let status = {
        let mut warn = |line: &[u8]| {
            messages.extend_from_slice(line);
            messages.push(b'\n');
        };
        let warn_ref: Warn<'_> = &mut warn;
        let mut db = BuildDatabase::new(
            out_base,
            &case.title,
            case.protein,
            &mut log,
            &time_of(index_date),
        )
        .unwrap();
        for (i, date) in column_dates.iter().enumerate() {
            if let Some(date) = date {
                db.set_column_date(i, date.as_bytes());
            }
        }
        let status = match build(&mut db, case, &mut *warn_ref) {
            Ok(()) | Err(Stop::Silent) => 0,
            Err(Stop::Error(error)) => {
                // CATCH_ALL
                warn_ref(error.to_string().as_bytes());
                error.exit_code()
            }
        };
        db.post_destruction_diagnostics(&mut *warn_ref);
        status
    };
    Run {
        status,
        stderr: messages,
        log: String::from_utf8_lossy(&log).into_owned(),
    }
}

const DESTRUCTOR_MESSAGE: &str = "Error: [makeblastdb] No sequences matched any of the masks \
     provided.\nPlease ensure that the -parse_seqids option is used in the\nfiltering program as \
     well as makeblastdb.\n";

const DB_EXTENSIONS: [&str; 14] = [
    "nin", "nsq", "nhr", "naa", "nab", "nac", "nal", "pin", "psq", "phr", "paa", "pab", "pac",
    "pal",
];

fn is_database_file(name: &str) -> bool {
    name.rsplit_once('.')
        .is_some_and(|(_, ext)| DB_EXTENSIONS.contains(&ext))
}

/// The names of the volumes of a database in `dir`: those the alias file
/// lists, or the database itself.
fn volume_bases(dir: &Path, name: &str, protein: bool) -> Vec<String> {
    let alias = dir.join(format!("{name}.{}al", if protein { 'p' } else { 'n' }));
    match std::fs::read_to_string(alias) {
        Ok(text) => text
            .lines()
            .find_map(|line| line.strip_prefix("DBLIST "))
            .expect("an alias file lists its volumes")
            .split(' ')
            .map(str::to_string)
            .collect(),
        Err(_) => vec![name.to_string()],
    }
}

fn stored_title(index_file: &Path, seqtype: SeqType) -> String {
    let bytes = std::fs::read(index_file).unwrap();
    let index = crate::seqdb::index::IndexFile::parse(&bytes, seqtype, index_file).unwrap();
    String::from_utf8(index.title).unwrap()
}

/// Every file `out_dir/name.*` equals the fixture's file of the same
/// extension and the fixture has no other file.
fn assert_same_files(out_dir: &Path, name: &str, fixture_dir: &Path, fixture: &str) {
    let suffix = |file: &str, stem: &str| file[stem.len()..].to_string();
    let written = files_named(out_dir, &format!("{name}."));
    let expected: Vec<String> = files_named(fixture_dir, &format!("{fixture}."))
        .iter()
        .filter(|file| is_database_file(file))
        .map(|file| suffix(file, fixture))
        .collect();
    let actual: Vec<String> = written.iter().map(|file| suffix(file, name)).collect();
    assert_eq!(actual, expected, "the files of {fixture}");
    for file in &written {
        let rest = suffix(file, name);
        let expected = std::fs::read(fixture_dir.join(format!("{fixture}{rest}"))).unwrap();
        let actual = std::fs::read(out_dir.join(file)).unwrap();
        assert!(
            actual == expected,
            "{fixture}{rest}: {} bytes written, {} expected; first difference at byte {:?}",
            actual.len(),
            expected.len(),
            actual.iter().zip(&expected).position(|(a, b)| a != b)
        );
    }
}

fn ranges_text(result: &MaskRanges) -> String {
    match result {
        MaskRanges::Ranges(ranges) if ranges.is_empty() => "none".to_string(),
        MaskRanges::Ranges(ranges) => ranges.iter().map(|(a, b)| format!("{a}-{b};")).collect(),
        other => panic!("{other:?}"),
    }
}

fn seqtype_of(protein: bool) -> SeqType {
    if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    }
}

/// The database `out_base` reads back as NCBI's blastdbcmd reads the
/// fixture: the table of algorithms and the ranges of every sequence.
fn assert_reads_back(out_base: &str, protein: bool, expected_dir: &Path, fixture: &str) {
    let db = SeqDb::open_with_masks(out_base, seqtype_of(protein)).unwrap();
    let algorithms_file = expected_dir.join(format!("{fixture}.algorithms.txt"));
    let table = db.mask_algorithms().unwrap();
    if algorithms_file.exists() {
        let expected = std::fs::read_to_string(&algorithms_file).unwrap();
        assert_eq!(table.descriptions().unwrap(), expected, "{fixture}");
    } else {
        assert!(table.ids().is_empty(), "{fixture}");
    }
    let rows = std::fs::read_to_string(expected_dir.join(format!("{fixture}.masks.tsv")))
        .unwrap_or_default();
    for line in rows.lines() {
        let mut fields = line.split('\t');
        let algo: i32 = fields.next().unwrap().parse().unwrap();
        let oid: u32 = fields.next().unwrap().parse().unwrap();
        let expected = fields.next().unwrap();
        let result = db.mask_ranges(oid, algo).unwrap();
        assert_eq!(
            ranges_text(&result),
            expected,
            "{fixture}: algorithm {algo}, sequence {oid}"
        );
    }
}

/// Rebuilds the NCBI database `fixture` of `fixture_dir` from `case` and
/// compares every file. The index date, the title and the dates of the
/// column files are the ones in NCBI's files.
fn assert_rebuilds(case: &Case, fixture_dir: &Path, fixture: &str, expected_stderr: &[u8]) {
    let letter = if case.protein { 'p' } else { 'n' };
    let seqtype = seqtype_of(case.protein);
    let bases = volume_bases(fixture_dir, fixture, case.protein);
    let first_index = fixture_dir.join(format!("{}.{letter}in", bases[0]));
    let index_date = stored_date(&first_index, seqtype);
    let column_dates: Vec<Option<String>> = bases
        .iter()
        .map(|base| column_date(&fixture_dir.join(format!("{base}.{letter}aa"))))
        .collect();
    let tmp = TempDir::new(&case.name);
    let run = run_case(case, &tmp.base(fixture), &index_date, &column_dates);
    assert_eq!(run.status, 0, "{fixture}: {}", run.stderr_text());
    assert_eq!(
        run.stderr_text(),
        String::from_utf8_lossy(expected_stderr),
        "{fixture}: standard error"
    );
    // Every mask file of these cases holds an object: the log names it.
    let names: Vec<&str> = case
        .mask_data
        .as_deref()
        .map_or(Vec::new(), |list| list.split(',').collect());
    assert_eq!(
        run.mask_files(),
        names,
        "{fixture}: the mask files in the log"
    );
    assert_same_files(&tmp.0, fixture, fixture_dir, fixture);
    assert_eq!(
        stored_title(&first_index, seqtype),
        case.title,
        "{fixture}: the title is the one of the fixture"
    );
    let expected_dir = fixture_dir.parent().unwrap().join("expected");
    if expected_dir.join(format!("{fixture}.masks.tsv")).exists() {
        assert_reads_back(&tmp.base(fixture), case.protein, &expected_dir, fixture);
    }
}

fn read(dir: &Path, name: &str) -> String {
    std::fs::read_to_string(dir.join(name)).unwrap_or_else(|error| panic!("{name}: {error}"))
}

/// The databases of `make_mask_fixtures.py` in format version 4. The mask
/// files are the ones that script gave makeblastdb (written by NCBI's
/// dustmasker and segmasker in binary, and hand-written text for windowmasker),
/// named here by their files in `tests/seqdb`.
fn fixture_case(name: &str) -> Case {
    let seqdb = seqdb_fixtures();
    let (protein, masks, ids, descs, max_file_size): (
        bool,
        &[&str],
        Option<&str>,
        Option<&str>,
        u64,
    ) = match name {
        "mask_nucl_v4" => (
            false,
            &["mask_dust.asnb", "mask_wm.asnt"],
            None,
            None,
            3_000_000_000,
        ),
        "mask_nucl_mv" => (false, &["mask_dust.asnb", "mask_wm.asnt"], None, None, 400),
        "mask_nucl_ids" => (
            false,
            &["mask_dust.asnb", "mask_wm.asnt"],
            Some("dustX,wmY"),
            Some("dust: from dustmasker,hand written"),
            3_000_000_000,
        ),
        "mask_nucl_z" => (
            false,
            &["mask_wm_opts.asnt"],
            Some("wmZ"),
            None,
            3_000_000_000,
        ),
        "mask_nucl_y" => (
            false,
            &["mask_wm_opts.asnt"],
            Some("wmY"),
            None,
            3_000_000_000,
        ),
        "mask_nucl_opts" => (
            false,
            &["mask_wm_opts.asnt", "mask_wm.asnt"],
            None,
            None,
            3_000_000_000,
        ),
        "mask_nucl_colon" => (false, &["mask_wm_colon.asnt"], None, None, 3_000_000_000),
        "mask_nucl_colon_id" => (
            false,
            &["mask_wm_colon.asnt"],
            Some("c:id"),
            None,
            3_000_000_000,
        ),
        "mask_prot_v4" => (true, &["mask_seg.asnb"], None, None, 3_000_000_000),
        "mask_nucl_plain" => (false, &[], None, None, 3_000_000_000),
        other => panic!("no fixture {other}"),
    };
    // makeblastdb's default title is the name of the input file; the two
    // fixtures with a title of their own state it in their index file.
    let fasta = if protein {
        "mask_prot.faa"
    } else {
        "mask_nucl.fna"
    };
    let title = match name {
        "mask_nucl_v4" => "LOSAT mask fixture, nucleotide, version 4",
        "mask_nucl_mv" => "masked nucleotide in volumes",
        _ => fasta,
    };
    Case {
        name: name.to_string(),
        protein,
        inputs: vec![seqdb.join(fasta)],
        mask_data: (!masks.is_empty()).then(|| masks.join(",")),
        mask_id: ids.map(str::to_string),
        mask_desc: descs.map(str::to_string),
        root: seqdb,
        max_file_size,
        title: title.to_string(),
    }
}

/// Every database that `make_mask_fixtures.py` built in format version 4,
/// from its FASTA file and the mask files NCBI's tools wrote (binary for dustmasker and
/// segmasker).
#[test]
fn the_writer_rebuilds_makeblastdb_s_masked_databases_byte_for_byte() {
    let fixtures = seqdb_fixtures().join("db");
    for name in [
        "mask_nucl_v4",
        "mask_nucl_mv",
        "mask_nucl_ids",
        "mask_nucl_z",
        "mask_nucl_y",
        "mask_nucl_opts",
        "mask_nucl_colon",
        "mask_nucl_colon_id",
        "mask_prot_v4",
        "mask_nucl_plain",
    ] {
        assert_rebuilds(&fixture_case(name), &fixtures, name, b"");
    }
}

/// The cases of `make_mask_writer_fixtures.py`.
fn writer_cases() -> Vec<Case> {
    let seqdb = seqdb_fixtures();
    let dir = seqdb.join("mask_writer");
    let table = read(&dir, "cases.tsv");
    let resolve = |file: &str| {
        if dir.join(file).exists() {
            dir.join(file)
        } else {
            seqdb.join(file)
        }
    };
    let option = |field: &str| (field != "-").then(|| field.to_string());
    table
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 8, "{line}");
            Case {
                name: f[0].to_string(),
                protein: f[1] == "prot",
                inputs: f[2].split(' ').map(resolve).collect(),
                mask_data: option(f[3]),
                mask_id: option(f[4]),
                mask_desc: option(f[5]),
                root: dir.clone(),
                max_file_size: f[6]
                    .strip_suffix('B')
                    .map_or(3_000_000_000, |n| n.parse().unwrap()),
                title: f[7].to_string(),
            }
        })
        .collect()
}

#[test]
fn the_writer_matches_makeblastdb_on_the_cases_of_the_writer_fixtures() {
    let dir = seqdb_fixtures().join("mask_writer");
    let cases = writer_cases();
    assert_eq!(cases.len(), 36);
    let mut built = 0;
    let mut failed = 0;
    for case in &cases {
        let name = &case.name;
        let expected_status = read(&dir.join("expected"), &format!("{name}.rc"));
        let expected_stderr =
            std::fs::read(dir.join("expected").join(format!("{name}.stderr"))).unwrap();
        if expected_status.trim() == "0" {
            built += 1;
            assert_rebuilds(case, &dir.join("db"), name, &expected_stderr);
        } else {
            failed += 1;
            assert_eq!(expected_status.trim(), "1", "{name}");
            assert!(files_named(&dir.join("db"), &format!("{name}.")).is_empty());
            let tmp = TempDir::new(name);
            let run = run_case(case, &tmp.base(name), "Oct 20, 2026  1:00 PM", &[]);
            assert_eq!(run.status, 1, "{name}: the run should fail");
            assert_eq!(
                run.stderr_text(),
                String::from_utf8_lossy(&expected_stderr),
                "{name}: standard error"
            );
            // The files of a failed run are removed.
            assert_eq!(
                files_named(&tmp.0, &format!("{name}.")),
                Vec::<String>::new(),
                "{name}: files left behind"
            );
        }
    }
    assert_eq!((built, failed), (25, 11), "built and failed cases");
}

// ---------------------------------------------------------------------
// tests/seqdb/mask_parse/expected.tsv
// ---------------------------------------------------------------------

/// What makeblastdb did with a file (`expected.tsv`).
struct Expected {
    dbtype: String,
    fasta: String,
    rc: i32,
    stderr: String,
    no_db: bool,
    /// `MASK` rows: algorithm id -> oid -> ranges.
    masks: BTreeMap<i32, BTreeMap<u32, String>>,
}

fn unescape(field: &str) -> String {
    if field == "\\e" {
        return String::new();
    }
    let bytes = field.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'\\' {
            out.push(bytes[i]);
            i += 1;
            continue;
        }
        match bytes[i + 1] {
            b'\\' => out.push(b'\\'),
            b't' => out.push(b'\t'),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b's' => out.push(b' '),
            b'e' => {}
            b'x' => {
                out.push(u8::from_str_radix(&field[i + 2..i + 4], 16).unwrap());
                i += 2;
            }
            other => panic!("unknown escape \\{}", other as char),
        }
        i += 2;
    }
    // The standard error of NCBI's programs is bytes; the tests compare it with ASCII.
    String::from_utf8_lossy(&out).into_owned()
}

fn read_expected(parse_dir: &Path) -> BTreeMap<String, Expected> {
    let text = read(parse_dir, "expected.tsv");
    let mut runs: BTreeMap<String, Expected> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "RUN" => {
                runs.insert(
                    f[1].to_string(),
                    Expected {
                        dbtype: f[2].to_string(),
                        fasta: f[3].to_string(),
                        rc: f[4].parse().unwrap(),
                        stderr: unescape(f[5]),
                        no_db: false,
                        masks: BTreeMap::new(),
                    },
                );
            }
            "NODB" => runs.get_mut(f[1]).unwrap().no_db = true,
            "MASK" => {
                runs.get_mut(f[1])
                    .unwrap()
                    .masks
                    .entry(f[2].parse().unwrap())
                    .or_default()
                    .insert(f[3].parse().unwrap(), unescape(f[4]));
            }
            other => panic!("unknown row {other}"),
        }
    }
    runs
}

/// Files that makeblastdb reads and LOSAT does not (`MaskFileError::Unsupported`): LOSAT's
/// outcome is its own refusal (status 1, "LOSAT makedb: mask data that LOSAT does not read: ..."),
/// whatever NCBI did with them. A new file that LOSAT refuses must be added here.
const LOSAT_UNSUPPORTED: [&str; 10] = [
    "h_unsupported_whole.asnt",
    "h_unsupported_pnt.asnt",
    "h_unsupported_mix.asnt",
    "h_unsupported_null.asnt",
    "h_unsupported_feat.asnt",
    "b_small_whole.asnb",
    "h_different_other_ids.asnt",
    // NCBI replaces a byte that is not printable in a `VisibleString` with '#' and logs it, and
    // warns about a UTF8String tag for a VisibleString; those lines are its only standard error.
    // LOSAT does not write them, so it refuses the string.
    "h_string_nonprint.asnt",
    "h_utf8_in_options.asnt",
    "b_utf8_options.asnb",
];

/// Files with which NCBI ends in `std::bad_alloc` (status 4, "BLAST ran out of memory"): a string
/// whose length cannot be allocated. Whether an allocation of a given size fails depends on the
/// host, so LOSAT does not tell these from the other damaged files and ends with status 255.
const NCBI_OUT_OF_MEMORY: [&str; 1] = ["b_length_huge.asnb"];

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
enum Class {
    /// "Unknown encoding for mask data in ..." (status 1): text and status equal NCBI's.
    UnknownEncoding,
    /// NCBI's serial reader throws (status 255): LOSAT's status, its own text.
    Serial,
    /// `std::bad_alloc` in NCBI (status 4): LOSAT's status 255, its own text.
    OutOfMemory,
    /// A missing continuation (`more TRUE` and the end of the file): status 0, no database.
    NoDatabase,
    /// LOSAT refuses a file that NCBI reads.
    Unsupported,
    /// Status, standard error and database equal NCBI's.
    Equal,
}

fn classify(file: &str, run: &Expected) -> Class {
    if run.rc == 1 && run.stderr.contains("Unknown encoding for mask data in '") {
        Class::UnknownEncoding
    } else if run.rc == 255 && run.stderr.starts_with("Error: NCBI C++ Exception:") {
        Class::Serial
    } else if run.rc == 4 {
        assert!(
            NCBI_OUT_OF_MEMORY.contains(&file),
            "{file}: a new file with NCBI's status 4"
        );
        Class::OutOfMemory
    } else if run.rc == 0 && run.no_db {
        Class::NoDatabase
    } else if LOSAT_UNSUPPORTED.contains(&file) {
        Class::Unsupported
    } else {
        assert!(
            !run.stderr.contains("NCBI C++ Exception")
                && !run.stderr.contains("(803.")
                && matches!(run.rc, 0 | 1),
            "{file}: a new kind of outcome: status {}, {}",
            run.rc,
            run.stderr
        );
        Class::Equal
    }
}

/// Every file of `tests/seqdb/mask_parse` goes through the flow of makeblastdb with
/// `-mask_data <file>` (the FASTA it was run with, no `-parse_seqids`).
#[test]
fn every_mask_parse_file_runs_end_to_end() {
    let seqdb = seqdb_fixtures();
    let parse_dir = seqdb.join("mask_parse");
    let expected = read_expected(&parse_dir);
    let mut tally: BTreeMap<String, usize> = BTreeMap::new();
    // The rows of blastdbcmd compared with the databases LOSAT wrote, and those NCBI listed
    // for databases that were written.
    let (mut compared, mut listed) = (0, 0);
    for (file, ncbi) in &expected {
        let class = classify(file, ncbi);
        *tally.entry(format!("{class:?}")).or_default() += 1;
        let protein = ncbi.dbtype == "prot";
        let fasta = if parse_dir.join(&ncbi.fasta).exists() {
            parse_dir.join(&ncbi.fasta)
        } else {
            seqdb.join(&ncbi.fasta)
        };
        let case = Case {
            name: file.clone(),
            protein,
            inputs: vec![fasta],
            mask_data: Some(file.clone()),
            mask_id: None,
            mask_desc: None,
            root: parse_dir.clone(),
            max_file_size: 3_000_000_000,
            title: ncbi.fasta.clone(),
        };
        let tmp = TempDir::new("mask-parse");
        let run = run_case(&case, &tmp.base("db"), "Oct 20, 2026  1:00 PM", &[]);
        let stderr = run.stderr_text();
        let no_files = || {
            assert_eq!(
                files_named(&tmp.0, "db."),
                Vec::<String>::new(),
                "{file}: files left behind"
            );
        };
        // The text of the destructor of the CBuildDatabase ends the standard error of both
        // programs when mask data were read before the failure.
        let ends_alike = || {
            assert_eq!(
                stderr.ends_with(DESTRUCTOR_MESSAGE),
                ncbi.stderr.ends_with(DESTRUCTOR_MESSAGE),
                "{file}: the destructor message\nLOSAT: {stderr}\nNCBI: {}",
                ncbi.stderr
            );
        };
        match class {
            Class::UnknownEncoding | Class::NoDatabase => {
                assert_eq!(run.status, ncbi.rc, "{file}: {stderr}");
                assert_eq!(stderr, ncbi.stderr, "{file}: standard error");
                no_files();
            }
            Class::Serial | Class::OutOfMemory => {
                assert_eq!(run.status, 255, "{file}: {stderr}");
                assert!(
                    stderr.starts_with(&format!("Error: mask data in '{file}': ")),
                    "{file}: {stderr}"
                );
                if class == Class::OutOfMemory {
                    assert_eq!(ncbi.rc, 4, "{file}");
                    assert_eq!(ncbi.stderr, "BLAST ran out of memory\n", "{file}");
                } else {
                    assert_eq!(ncbi.rc, 255, "{file}");
                }
                ends_alike();
                no_files();
            }
            Class::Unsupported => {
                assert_eq!(run.status, 1, "{file}: {stderr}");
                assert!(
                    stderr.starts_with("LOSAT makedb: mask data that LOSAT does not read: "),
                    "{file}: {stderr}"
                );
                no_files();
            }
            Class::Equal => {
                assert_eq!(run.status, ncbi.rc, "{file}: {stderr}");
                assert_eq!(stderr, ncbi.stderr, "{file}: standard error");
                if run.status == 0 {
                    compared += assert_database(&tmp.base("db"), protein, file, ncbi);
                    listed += ncbi.masks.values().map(BTreeMap::len).sum::<usize>();
                } else {
                    no_files();
                }
            }
        }
    }
    assert_eq!(compared, listed);
    assert!(compared > 300, "{compared} rows compared");
    // The classes and their sizes: a change in the fixtures or in the flow shows here.
    let tally: Vec<(String, usize)> = tally.into_iter().collect();
    assert_eq!(
        tally,
        [
            ("Equal".to_string(), 68),
            ("NoDatabase".to_string(), 3),
            ("OutOfMemory".to_string(), 1),
            ("Serial".to_string(), 39),
            ("UnknownEncoding".to_string(), 7),
            ("Unsupported".to_string(), 10),
        ]
    );
}

/// The database that LOSAT wrote reads back as blastdbcmd read NCBI's: the algorithms it lists
/// and the ranges of every sequence.
/// Returns the number of rows compared.
fn assert_database(out_base: &str, protein: bool, file: &str, ncbi: &Expected) -> usize {
    let db = SeqDb::open_with_masks(out_base, seqtype_of(protein))
        .unwrap_or_else(|error| panic!("{file}: {error}"));
    if ncbi.masks.is_empty() {
        // No sequence of the database has a range (the masks name other sequences, or there
        // are none): blastdbcmd lists no algorithm, and the reader finds none.
        let listed = db.mask_algorithms().map_or(Vec::new(), |table| table.ids());
        assert_eq!(
            listed,
            Vec::<i32>::new(),
            "{file}: algorithms without ranges"
        );
        return 0;
    }
    let table = db.mask_algorithms().unwrap();
    let listed: BTreeSet<i32> = table.ids().into_iter().collect();
    let expected: BTreeSet<i32> = ncbi.masks.keys().copied().collect();
    assert_eq!(listed, expected, "{file}: the algorithms");
    let mut compared = 0;
    for (algo, rows) in &ncbi.masks {
        for (oid, text) in rows {
            let result = db.mask_ranges(*oid, *algo).unwrap();
            assert_eq!(
                ranges_text(&result),
                *text,
                "{file}: algorithm {algo}, sequence {oid}"
            );
            compared += 1;
        }
    }
    compared
}

// ---------------------------------------------------------------------
// The options of x_ProcessMaskData
// ---------------------------------------------------------------------

/// The directory the option cases look their mask files up in: `dust.asnt` and `wm.asnt`
/// of `tests/seqdb/mask_writer`, a directory, a file that is not an ASN.1 file, an empty file
/// and (on Unix) a link to nothing.
fn option_root() -> TempDir {
    let tmp = TempDir::new("mask-options");
    let writer = seqdb_fixtures().join("mask_writer");
    for file in ["dust.asnt", "wm.asnt"] {
        std::fs::copy(writer.join(file), tmp.0.join(file)).unwrap();
    }
    std::fs::create_dir(tmp.0.join("adir")).unwrap();
    std::fs::write(tmp.0.join("empty.txt"), b"").unwrap();
    #[cfg(unix)]
    std::os::unix::fs::symlink("nowhere.asnt", tmp.0.join("dangling.asnt")).unwrap();
    tmp
}

fn option_case(
    root: &Path,
    mask_data: &str,
    mask_id: Option<&str>,
    mask_desc: Option<&str>,
) -> Case {
    Case {
        name: "options".to_string(),
        protein: false,
        inputs: vec![seqdb_fixtures().join("mask_nucl.fna")],
        mask_data: Some(mask_data.to_string()),
        mask_id: mask_id.map(str::to_string),
        mask_desc: mask_desc.map(str::to_string),
        root: root.to_path_buf(),
        max_file_size: 3_000_000_000,
        title: "mask_nucl.fna".to_string(),
    }
}

const MISSING_EMPTY: &str = "Error: [makeblastdb] Ignoring mask file '' as it does not exist.\n";

fn missing(name: &str) -> String {
    format!("Error: [makeblastdb] Ignoring mask file '{name}' as it does not exist.\n")
}

/// The answers of makeblastdb 2.17.0 to `-in mask_nucl.fna -dbtype nucl -mask_data A [-mask_id B
/// [-mask_desc C]]` (status, standard error, the "Mask file:" lines of the log, whether the
/// database has a mask column), found with the oracle in a directory with the same files
/// (`opt.py` of the M4 brief C: the cases `empty_value` ... `missing_all`).
#[test]
fn process_mask_data_follows_makeblastdb_on_its_options() {
    let tmp = option_root();
    let no_stderr = String::new;
    // (mask_data, mask_id, mask_desc, status, standard error, mask files, column)
    #[allow(clippy::type_complexity)]
    let mut cases: Vec<(
        &str,
        Option<&str>,
        Option<&str>,
        i32,
        String,
        Vec<&str>,
        bool,
    )> = vec![
        // The list of files.
        (
            "",
            None,
            None,
            1,
            "BLAST options error: mask_data option found, but no files were specified.\n".into(),
            vec![],
            false,
        ),
        (
            " \t ",
            None,
            None,
            1,
            "BLAST options error: mask_data option found, but no files were specified.\n".into(),
            vec![],
            false,
        ),
        (",", None, None, 0, MISSING_EMPTY.repeat(2), vec![], false),
        (
            "dust.asnt,",
            None,
            None,
            0,
            MISSING_EMPTY.into(),
            vec!["dust.asnt"],
            true,
        ),
        (
            "nonexist.asnt",
            None,
            None,
            0,
            missing("nonexist.asnt"),
            vec![],
            false,
        ),
        (
            "nonexist.asnt,dust.asnt",
            None,
            None,
            0,
            missing("nonexist.asnt"),
            vec!["dust.asnt"],
            true,
        ),
        (
            "nonexist.asnt,nonexist2.asnt",
            None,
            None,
            0,
            missing("nonexist.asnt") + &missing("nonexist2.asnt"),
            vec![],
            false,
        ),
        (
            " dust.asnt ",
            None,
            None,
            0,
            no_stderr(),
            vec!["dust.asnt"],
            true,
        ),
        // Only the ends of the whole option are trimmed.
        (
            "dust.asnt, wm.asnt",
            None,
            None,
            0,
            missing(" wm.asnt"),
            vec!["dust.asnt"],
            true,
        ),
        // Only a regular file exists: a directory and (the oracle) a device, a link to nothing.
        ("adir", None, None, 0, missing("adir"), vec![], false),
        (
            "/dev/null",
            None,
            None,
            0,
            missing("/dev/null"),
            vec![],
            false,
        ),
        (
            "empty.txt",
            None,
            None,
            1,
            "BLAST options error: Unknown encoding for mask data in 'empty.txt'\n".into(),
            vec![],
            false,
        ),
        // -mask_id
        (
            "dust.asnt",
            Some("5"),
            None,
            1,
            "BLAST options error: mask_id can not be numeric.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,-0"),
            None,
            1,
            "BLAST options error: mask_id can not be numeric.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,2147483647"),
            None,
            1,
            "BLAST options error: mask_id can not be numeric.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,2147483648"),
            None,
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a, 5"),
            None,
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,+"),
            None,
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a5,5a"),
            None,
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a"),
            None,
            1,
            "BLAST options error: the size of mask_id does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b,"),
            None,
            1,
            "BLAST options error: the size of mask_id does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt",
            Some(""),
            None,
            1,
            "BLAST options error: the size of mask_id does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            None,
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "nonexist.asnt,dust.asnt",
            Some("a,b"),
            None,
            0,
            missing("nonexist.asnt"),
            vec!["dust.asnt"],
            true,
        ),
        // -mask_desc
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            Some("x"),
            1,
            "BLAST options error: the size of mask_desc does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            Some("x,y,"),
            1,
            "BLAST options error: the size of mask_desc does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            Some(""),
            1,
            "BLAST options error: the size of mask_desc does not match that of mask_data.\n".into(),
            vec![],
            false,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            Some(","),
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
        (
            "dust.asnt,wm.asnt",
            Some("a,b"),
            Some("d1, d2"),
            0,
            no_stderr(),
            vec!["dust.asnt", "wm.asnt"],
            true,
        ),
    ];
    #[cfg(unix)]
    cases.push((
        "dangling.asnt",
        None,
        None,
        0,
        missing("dangling.asnt"),
        vec![],
        false,
    ));
    for (mask_data, mask_id, mask_desc, status, stderr, mask_files, column) in cases {
        let label = format!("{mask_data:?} {mask_id:?} {mask_desc:?}");
        let case = option_case(&tmp.0, mask_data, mask_id, mask_desc);
        let out = TempDir::new("mask-options-out");
        let run = run_case(&case, &out.base("db"), "Oct 20, 2026  1:00 PM", &[]);
        assert_eq!(run.status, status, "{label}: {}", run.stderr_text());
        assert_eq!(run.stderr_text(), stderr, "{label}: standard error");
        assert_eq!(
            run.mask_files(),
            mask_files,
            "{label}: the mask files in the log"
        );
        let written = files_named(&out.0, "db.");
        if status == 0 {
            assert!(
                written.iter().any(|f| f == "db.nin"),
                "{label}: {written:?}"
            );
            assert_eq!(
                written.iter().any(|f| f == "db.naa"),
                column,
                "{label}: the mask column in {written:?}"
            );
        } else {
            assert_eq!(written, Vec::<String>::new(), "{label}: files left behind");
        }
    }
}

/// The oracle (file with mode 000 that makeblastdb cannot open): "Unknown encoding", like a
/// file that is not an ASN.1 file. Skipped where the process can open it anyway (root).
#[cfg(unix)]
#[test]
fn a_mask_file_that_cannot_be_opened_is_an_unknown_encoding() {
    use std::os::unix::fs::PermissionsExt;
    let tmp = option_root();
    let file = tmp.0.join("noperm.asnt");
    std::fs::copy(tmp.0.join("dust.asnt"), &file).unwrap();
    std::fs::set_permissions(&file, std::fs::Permissions::from_mode(0)).unwrap();
    if std::fs::File::open(&file).is_ok() {
        return;
    }
    let case = option_case(&tmp.0, "noperm.asnt", None, None);
    let out = TempDir::new("mask-options-out");
    let run = run_case(&case, &out.base("db"), "Oct 20, 2026  1:00 PM", &[]);
    assert_eq!(run.status, 1);
    assert_eq!(
        run.stderr_text(),
        "BLAST options error: Unknown encoding for mask data in 'noperm.asnt'\n"
    );
    assert_eq!(files_named(&out.0, "db."), Vec::<String>::new());
}

/// `-gi_mask_name` (reachable only with `-gi_mask`, which LOSAT does not offer): the list is
/// split and checked like the others, and the name replaces the file name in the registration.
#[test]
fn gi_mask_names_are_split_and_replace_the_file_name() {
    let tmp = option_root();
    let writer = seqdb_fixtures().join("mask_writer");
    std::fs::copy(writer.join("oth_d.asnt"), tmp.0.join("oth_d.asnt")).unwrap();
    let mut log = Vec::new();
    let out = TempDir::new("mask-options-out");
    let mut db = BuildDatabase::new(
        &out.base("db"),
        "t",
        false,
        &mut log,
        &time_of("Oct 20, 2026  1:00 PM"),
    )
    .unwrap();
    let mut warn = |_: &[u8]| {};
    let mut options = MaskDataOptions {
        mask_data: Some("dust.asnt,wm.asnt"),
        gi_mask_name: Some("x"),
        root: Some(&tmp.0),
        ..MaskDataOptions::default()
    };
    let error = process_mask_data(&mut db, &options, &mut warn).unwrap_err();
    assert_eq!(
        error.to_string(),
        "BLAST options error: the size of gi_mask_name does not match that of mask_data."
    );
    assert_eq!(error.into_makedb_error().unwrap().exit_code(), 1);
    // `other` is the program whose name is stored: oth_d.asnt, registered as "gi name".
    options.mask_data = Some("oth_d.asnt");
    options.gi_mask_name = Some("gi name");
    process_mask_data(&mut db, &options, &mut warn).unwrap();
    // The same file under its own name is another algorithm, the same name again a duplicate.
    options.gi_mask_name = None;
    process_mask_data(&mut db, &options, &mut warn).unwrap();
    options.gi_mask_name = Some("gi name");
    let error = process_mask_data(&mut db, &options, &mut warn).unwrap_err();
    assert_eq!(
        error.to_string(),
        "BLAST Database creation error: Duplicate masking algorithm found."
    );
}

// ---------------------------------------------------------------------
// The mapping of the outcomes
// ---------------------------------------------------------------------

/// `MaskDataError::into_makedb_error`: the one place that turns an outcome of
/// `x_ProcessMaskData` into makeblastdb's text and status.
#[test]
fn outcomes_become_makeblastdb_s_status_and_text() {
    let io = || std::io::Error::from(std::io::ErrorKind::PermissionDenied);
    let path = || PathBuf::from("m.asnb");
    let describe = |error: MaskDataError| {
        error
            .into_makedb_error()
            .map(|e| (e.exit_code(), e.to_string()))
    };
    assert_eq!(
        describe(MaskDataError::Ncbi(WriteDbError::Options("o".into()))),
        Some((1, "BLAST options error: o".to_string()))
    );
    assert_eq!(
        describe(MaskDataError::Ncbi(WriteDbError::Creation("c".into()))),
        Some((1, "BLAST Database creation error: c".to_string()))
    );
    assert_eq!(
        describe(MaskDataError::File(MaskFileError::Malformed {
            path: path(),
            message: "bad".into()
        })),
        Some((255, "Error: mask data in 'm.asnb': bad".to_string()))
    );
    assert_eq!(
        describe(MaskDataError::File(MaskFileError::Unsupported(
            "`whole`".into()
        ))),
        Some((
            1,
            "LOSAT makedb: mask data that LOSAT does not read: `whole`".to_string()
        ))
    );
    assert_eq!(
        describe(MaskDataError::File(MaskFileError::EofInContinuation {
            path: path()
        })),
        None
    );
    // The reader's own NCBI-worded error is the same as one raised here.
    let from_reader = MaskDataError::from(MaskFileError::Ncbi(WriteDbError::Options("u".into())));
    assert!(matches!(from_reader, MaskDataError::Ncbi(_)));
}

/// A database built without mask data lists and writes no column files.
#[test]
fn a_database_without_mask_data_has_no_column_files() {
    let tmp = TempDir::new("mask-plain");
    let mut db = WriteDb::new(&tmp.base("db"), true, b"t", b"Oct 12, 2026  3:22 PM", 100);
    for i in 1..=5 {
        db.add_sequence(bioseq(i, 40)).unwrap();
    }
    db.close().unwrap();
    for file in db.list_files() {
        assert!(
            !file.ends_with("aa") && !file.ends_with("ab") && !file.ends_with("ac"),
            "{file}"
        );
    }
    assert!(files_named(&tmp.0, "db.")
        .iter()
        .all(|f| is_database_file(f)));
    assert!(files_named(&tmp.0, "db.").iter().all(|f| {
        let ext = f.rsplit('.').next().unwrap();
        ["pin", "psq", "phr", "pal"].contains(&ext)
    }));
}

/// The message of the destructor of `CBuildDatabase`: posted only when mask data
/// were read and no sequence matched.
#[test]
fn the_message_of_a_run_without_matching_masks_is_ncbi_s() {
    let mut log = Vec::new();
    let tmp = TempDir::new("mask-destructor");
    let mut db = BuildDatabase::new(
        &tmp.base("db"),
        "t",
        true,
        &mut log,
        &time_of("Oct 12, 2026  3:22 PM"),
    )
    .unwrap();
    let mut messages: Vec<u8> = Vec::new();
    let mut warn = |line: &[u8]| {
        messages.extend_from_slice(line);
        messages.push(b'\n');
    };
    // No mask data source: no message.
    db.post_destruction_diagnostics(&mut warn);
    assert!(!db.found_matching_masks());
    db.set_mask_data_source(super::MaskedRangeSet::new());
    db.post_destruction_diagnostics(&mut warn);
    assert_eq!(String::from_utf8_lossy(&messages), DESTRUCTOR_MESSAGE);
}
