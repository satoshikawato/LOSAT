//! Tests of the database writer against databases written by makeblastdb
//! 2.17.0: the fixtures of `LOSAT/tests/seqdb` and the runs recorded in
//! `tests/makeblastdb_app_cases.rs`. The tests run no NCBI program.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use super::bioseq::{Bioseq, SeqData};
use super::build_db::BuildDatabase;
use super::makeblastdb::{run_with, MakedbArgs};
use super::platform::Platform;
use super::time::LocalTime;
use super::writer::{DbVersion, Indices, WriteDb};
use crate::seqdb::SeqType;

mod makeblastdb_app_cases;
use makeblastdb_app_cases::APP_CASES;
mod makeblastdb_v5_cases;
use makeblastdb_v5_cases::{Expected, V5Case, V5_CASES};

/// The platform of every expected file here: makeblastdb 2.17.0 on Linux x86-64. The tests
/// write its files on whatever platform they run on (`Platform`).
const LINUX: Platform = Platform::LINUX_X86_64;

fn seqdb_fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb")
}

/// A directory of this test run, removed when the value is dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.subsec_nanos());
        let dir = std::env::temp_dir().join(format!(
            "losat-writedb-{label}-{}-{}-{nanos}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn base(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The creation date stored in an index file of makeblastdb.
fn stored_date(index_file: &Path, seqtype: SeqType) -> String {
    let bytes = std::fs::read(index_file).unwrap();
    let index = crate::seqdb::index::IndexFile::parse(&bytes, seqtype, index_file).unwrap();
    let date: Vec<u8> = index.date.iter().copied().take_while(|&b| b != 0).collect();
    String::from_utf8(date).unwrap()
}

/// The local time that `LocalTime::index_date` formats as `date`
/// ("Oct 12, 2026  3:22 PM").
fn time_of(date: &str) -> LocalTime {
    let months = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let words: Vec<&str> = date.split_whitespace().collect();
    let [month, day, year, clock, half] = words.as_slice() else {
        panic!("unexpected date {date:?}");
    };
    let (hour, minute) = clock.split_once(':').unwrap();
    let hour12: u32 = hour.parse().unwrap();
    let time = LocalTime {
        year: year.parse().unwrap(),
        month: months.iter().position(|m| m == month).unwrap() as u32 + 1,
        day: day.trim_end_matches(',').parse().unwrap(),
        hour: hour12 % 12 + if *half == "PM" { 12 } else { 0 },
        minute: minute.parse().unwrap(),
        second: 0,
    };
    assert_eq!(time.index_date(), date);
    time
}

fn files_named(dir: &Path, prefix: &str) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with(prefix))
        .collect();
    names.sort();
    names
}

/// The files of a database that the fixtures of `LOSAT/tests/seqdb` keep: the
/// index, sequence, header and alias files and the LMDB file of version 5
/// (`make_fixtures.py`, `KEEP`).
fn kept_by_fixtures(file: &str) -> bool {
    ["in", "sq", "hr", "al", "db"]
        .iter()
        .any(|ext| file.ends_with(&format!(".p{ext}")) || file.ends_with(&format!(".n{ext}")))
}

/// Every file `out_dir/name.*` that the fixtures keep equals the fixture's
/// file of the same extension, and the fixture has no other such file.
fn assert_same_database(out_dir: &Path, name: &str, fixture_dir: &Path, fixture: &str) {
    let written: Vec<String> = files_named(out_dir, &format!("{name}."))
        .into_iter()
        .filter(|file| kept_by_fixtures(file))
        .collect();
    assert!(!written.is_empty(), "no file {name}.* was written");
    let suffix = |file: &str, stem: &str| file[stem.len()..].to_string();
    for file in &written {
        let rest = suffix(file, name);
        let expected = std::fs::read(fixture_dir.join(format!("{fixture}{rest}")))
            .unwrap_or_else(|error| panic!("{fixture}{rest}: {error}"));
        let actual = std::fs::read(out_dir.join(file)).unwrap();
        assert!(
            actual == expected,
            "{fixture}{rest}: {} bytes written, {} expected; first difference at byte {:?}",
            actual.len(),
            expected.len(),
            actual.iter().zip(&expected).position(|(a, b)| a != b)
        );
    }
    let expected: Vec<String> = files_named(fixture_dir, &format!("{fixture}."))
        .iter()
        .filter(|file| kept_by_fixtures(file))
        .map(|file| suffix(file, fixture))
        .collect();
    let actual: Vec<String> = written.iter().map(|file| suffix(file, name)).collect();
    assert_eq!(actual, expected, "the files of {fixture}");
}

/// Records of a plain FASTA file (a definition line, then sequence lines).
fn plain_fasta(path: &Path) -> Vec<(Vec<u8>, Vec<u8>)> {
    let text = std::fs::read(path).unwrap();
    let mut records: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for line in text.split(|&b| b == b'\n') {
        if let Some(title) = line.strip_prefix(b">") {
            records.push((title.to_vec(), Vec::new()));
        } else if let Some(last) = records.last_mut() {
            last.1.extend_from_slice(line);
        }
    }
    records
}

fn protein_bioseqs(path: &Path) -> Vec<Bioseq> {
    plain_fasta(path)
        .into_iter()
        .enumerate()
        .map(|(i, (title, residues))| Bioseq {
            local_id: (i + 1).to_string(),
            title,
            is_aa: true,
            length: residues.len() as u64,
            // The reader hands over capital letters.
            seq_data: Some(SeqData::Ncbieaa(residues.to_ascii_uppercase())),
            parsed: None,
        })
        .collect()
}

/// The writer alone, from sequences built by hand: the single-volume and the
/// six-volume protein databases of `LOSAT/tests/seqdb` (version 4) and the
/// single-volume one of version 5 with its LMDB file.
#[test]
fn the_writer_reproduces_makeblastdb_s_protein_volumes() {
    let fixtures = seqdb_fixtures();
    let cases = [
        (
            "prot_v4",
            "prot_v4.pin",
            "LOSAT seqdb fixture, protein, version 4",
            3_000_000_000u64,
            DbVersion::V4,
        ),
        (
            "prot_mv",
            "prot_mv.00.pin",
            "protein in volumes",
            300,
            DbVersion::V4,
        ),
        (
            "prot_v5",
            "prot_v5.pin",
            "prot.faa",
            3_000_000_000u64,
            DbVersion::V5,
        ),
    ];
    for (fixture, index_file, title, max_file_size, version) in cases {
        let date = stored_date(&fixtures.join("db").join(index_file), SeqType::Protein);
        let tmp = TempDir::new(fixture);
        let mut db = WriteDb::new(
            &tmp.base(fixture),
            true,
            title.as_bytes(),
            date.as_bytes(),
            max_file_size,
            version,
            LINUX,
            Indices::default(),
        );
        for bioseq in protein_bioseqs(&fixtures.join("prot.faa")) {
            db.extract_bioseq_deflines(&bioseq);
            db.add_sequence(bioseq).unwrap();
        }
        db.close().unwrap();
        assert!(db.take_diagnostics().is_empty());
        // The alias file of several volumes names them, so the database is
        // built under the fixture's name.
        assert_same_database(&tmp.0, fixture, &fixtures.join("db"), fixture);
    }
}

/// The platform reaches the files where NCBI's source asks it (`Platform`): with Windows' values
/// (a 4-byte `long`, map sizes of 500000 bytes, LMDB files as long as the map) the date in the
/// index file has 4 NUL bytes more or fewer, and the two LMDB files are 500000 bytes long; the
/// sequence and header files are Linux's. What makeblastdb.exe writes was not compared.
#[test]
fn the_platform_reaches_the_index_header_and_the_lmdb_files() {
    use super::platform::{ncbi_platform, NcbiOs, PageSize};
    let fixtures = seqdb_fixtures();
    let windows = ncbi_platform(NcbiOs::Windows, PageSize::Fixed(4096), 4, u64::MAX);
    let mut built = Vec::new();
    for platform in [LINUX, windows] {
        let tmp = TempDir::new("platform");
        let mut db = WriteDb::new(
            &tmp.base("db"),
            true,
            b"prot.faa",
            b"Oct 12, 2026  9:00 PM",
            3_000_000_000,
            DbVersion::V5,
            platform,
            Indices::default(),
        );
        for bioseq in protein_bioseqs(&fixtures.join("prot.faa")) {
            db.extract_bioseq_deflines(&bioseq);
            db.add_sequence(bioseq).unwrap();
        }
        db.close().unwrap();
        let read = |ext: &str| std::fs::read(tmp.0.join(format!("db.{ext}"))).unwrap();
        built.push([
            read("pin"),
            read("psq"),
            read("phr"),
            read("pdb"),
            read("ptf"),
        ]);
    }
    let [linux, windows] = [&built[0], &built[1]];
    assert_eq!(linux[0].len().abs_diff(windows[0].len()), 4, ".pin");
    assert_eq!(
        (&linux[1], &linux[2]),
        (&windows[1], &windows[2]),
        ".psq, .phr"
    );
    assert_eq!((linux[3].len(), linux[4].len()), (20480, 16384));
    assert_eq!((windows[3].len(), windows[4].len()), (500_000, 500_000));
}

/// One build through `BuildDatabase` with the FASTA source: the result, the
/// log, and the messages for the standard error stream.
fn build(
    fasta: &Path,
    is_protein: bool,
    title: &str,
    max_file_size: u64,
    date: &str,
    out_base: &str,
) -> (Result<(), super::WriteDbError>, String, Vec<u8>) {
    build_version(
        fasta,
        is_protein,
        title,
        max_file_size,
        date,
        out_base,
        DbVersion::V4,
    )
}

fn build_version(
    fasta: &Path,
    is_protein: bool,
    title: &str,
    max_file_size: u64,
    date: &str,
    out_base: &str,
    version: DbVersion,
) -> (Result<(), super::WriteDbError>, String, Vec<u8>) {
    let mut log = Vec::new();
    let mut messages = Vec::new();
    let result = (|| {
        let mut db = BuildDatabase::new(
            out_base,
            title,
            is_protein,
            &mut log,
            &time_of(date),
            version,
            LINUX,
            Indices::default(),
        )?;
        db.set_max_file_size(max_file_size);
        let mut warn = |line: &[u8]| {
            messages.extend_from_slice(line);
            messages.push(b'\n');
        };
        db.add_fasta(std::fs::File::open(fasta).unwrap(), &mut warn)?;
        db.end_build(false, &mut warn).map(|_| ())
    })();
    (result, String::from_utf8_lossy(&log).into_owned(), messages)
}

/// The databases of `LOSAT/tests/seqdb` that makeblastdb wrote without parsed
/// Seq-ids, in format versions 4 and 5, from their FASTA files.
#[test]
fn fasta_files_become_makeblastdb_s_databases() {
    let fixtures = seqdb_fixtures();
    let cases = [
        (
            "prot_v4",
            "prot.faa",
            true,
            "prot_v4.pin",
            "LOSAT seqdb fixture, protein, version 4",
            3_000_000_000u64,
            DbVersion::V4,
        ),
        (
            "nucl_v4",
            "nucl.fna",
            false,
            "nucl_v4.nin",
            "LOSAT seqdb fixture, nucleotide, version 4",
            3_000_000_000,
            DbVersion::V4,
        ),
        (
            "prot_mv",
            "prot.faa",
            true,
            "prot_mv.00.pin",
            "protein in volumes",
            300,
            DbVersion::V4,
        ),
        (
            "prot_v5",
            "prot.faa",
            true,
            "prot_v5.pin",
            "prot.faa",
            3_000_000_000,
            DbVersion::V5,
        ),
        (
            "nucl_v5",
            "nucl.fna",
            false,
            "nucl_v5.nin",
            "nucl.fna",
            3_000_000_000,
            DbVersion::V5,
        ),
        (
            "nucl_mv",
            "nucl.fna",
            false,
            "nucl_mv.00.nin",
            "nucleotide in volumes",
            300,
            DbVersion::V5,
        ),
    ];
    for (fixture, fasta, is_protein, index_file, title, max_file_size, version) in cases {
        let seqtype = if is_protein {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        let date = stored_date(&fixtures.join("db").join(index_file), seqtype);
        let tmp = TempDir::new(fixture);
        let (result, log, messages) = build_version(
            &fixtures.join(fasta),
            is_protein,
            title,
            max_file_size,
            &date,
            &tmp.base(fixture),
            version,
        );
        result.unwrap();
        assert!(
            messages.is_empty(),
            "{}",
            String::from_utf8_lossy(&messages)
        );
        // The alias file of several volumes names them, so the database is
        // built under the fixture's name.
        assert_same_database(&tmp.0, fixture, &fixtures.join("db"), fixture);
        let mol_type = if is_protein { "Protein" } else { "Nucleotide" };
        assert!(
            log.starts_with("\n\nBuilding a new DB, current time: "),
            "{log:?}"
        );
        assert!(log.contains(&format!(
            "\nNew DB title:  {title}\nSequence type: {mol_type}\nKeep MBits: T\n"
        )));
        assert!(log.contains("\nAdding sequences from FASTA; added "));
        assert!(log.ends_with(" seconds.\n\n\n"), "{log:?}");
    }
}

/// A database of the same name is removed before the new one is written:
/// the volumes and the alias file of a six-volume database do not survive a
/// single-volume build.
#[test]
fn an_existing_database_of_the_same_name_is_removed() {
    let fixtures = seqdb_fixtures();
    let date = stored_date(&fixtures.join("db/prot_v4.pin"), SeqType::Protein);
    let tmp = TempDir::new("replace");
    let base = tmp.base("db");
    let fasta = fixtures.join("prot.faa");
    let (result, _, _) = build(&fasta, true, "protein in volumes", 300, &date, &base);
    result.unwrap();
    assert_eq!(files_named(&tmp.0, "db.").len(), 6 * 3 + 1);
    // A file makeblastdb would also remove, and one it would not.
    std::fs::write(tmp.0.join("db.00.pjs"), b"{}").unwrap();
    std::fs::write(tmp.0.join("db.txt"), b"keep").unwrap();
    let (result, log, _) = build(
        &fasta,
        true,
        "LOSAT seqdb fixture, protein, version 4",
        3_000_000_000,
        &date,
        &base,
    );
    result.unwrap();
    assert!(
        log.contains("\nDeleted existing Protein BLAST database named "),
        "{log:?}"
    );
    assert_eq!(
        files_named(&tmp.0, "db."),
        ["db.phr", "db.pin", "db.psq", "db.txt"]
    );
    std::fs::remove_file(tmp.0.join("db.txt")).unwrap();
    assert_same_database(&tmp.0, "db", &fixtures.join("db"), "prot_v4");
}

/// The log with the markers of `makeblastdb_app_cases.rs` in place of the
/// build time, the output path and the seconds of the timed line.
fn log_with_markers(log: &[u8]) -> Vec<u8> {
    let lines: Vec<Vec<u8>> = log
        .split(|&b| b == b'\n')
        .map(|line| {
            let text = String::from_utf8_lossy(line);
            if text.starts_with("Building a new DB, current time: ") {
                b"Building a new DB, current time: <T>".to_vec()
            } else if text.starts_with("New DB name:   ") {
                b"New DB name:   <OUT>".to_vec()
            } else if let (true, Some(at)) = (
                text.starts_with("Adding sequences from FASTA; added "),
                text.find(" sequences in "),
            ) {
                format!("{} sequences in <S> seconds.", &text[..at]).into_bytes()
            } else {
                line.to_vec()
            }
        })
        .collect();
    lines.join(&b'\n')
}

fn replace_bytes(text: &[u8], from: &[u8], to: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.windows(from.len()).position(|w| w == from) {
        out.extend_from_slice(&rest[..at]);
        out.extend_from_slice(to);
        rest = &rest[at + from.len()..];
    }
    out.extend_from_slice(rest);
    out
}

/// `losat makeblastdb` from its arguments to its exit status on the inputs of
/// `makeblastdb_app_cases.rs`: the files, the log, the messages and the exit
/// status are makeblastdb's.
#[test]
fn makedb_runs_as_makeblastdb_on_edge_inputs() {
    for case in APP_CASES {
        let tmp = TempDir::new(case.name);
        let input = tmp.base("in.fa");
        std::fs::write(&input, case.input).unwrap();
        let out_dir = tmp.0.join("out");
        std::fs::create_dir(&out_dir).unwrap();
        let base = out_dir.join("db").to_string_lossy().into_owned();
        let args = MakedbArgs {
            input: input.clone(),
            input_type: "fasta".to_string(),
            dbtype: if case.protein { "prot" } else { "nucl" }.to_string(),
            title: Some(case.title.to_string()),
            out: Some(base.clone()),
            blastdb_version: 4,
            max_file_sz: case.max_file_sz.unwrap_or("3GB").to_string(),
            logfile: None,
            parse_seqids: false,
            hash_index: false,
            taxid: None,
            taxid_map: None,
        };
        // A run that writes no database stores no date; any time will do.
        let now = time_of(if case.date.is_empty() {
            "Oct 12, 2026  9:00 PM"
        } else {
            case.date
        });
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = run_with(
            &args,
            Some(&now),
            LINUX,
            &mut std::io::empty(),
            &mut out,
            &mut err,
        );
        let name = case.name;
        let show = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
        let err = replace_bytes(&err, input.as_bytes(), b"<IN>");
        let err = replace_bytes(&err, base.as_bytes(), b"<OUT>");
        assert!(
            err == case.stderr,
            "{name}: standard error\n{}\nexpected\n{}",
            show(&err),
            show(case.stderr)
        );
        let log = log_with_markers(&out);
        assert!(
            log == case.stdout,
            "{name}: standard output\n{}\nexpected\n{}",
            show(&log),
            show(case.stdout)
        );
        assert_eq!(status, case.status, "{name}: exit status");
        // The recorded files are those of the volumes; makeblastdb also
        // writes the metadata file after every build that made a volume.
        let mut written = files_named(&out_dir, "db.");
        let metadata = written
            .iter()
            .position(|file| file == "db.pjs" || file == "db.njs");
        assert_eq!(
            metadata.is_some(),
            case.status == 0 && !case.files.is_empty(),
            "{name}: metadata file"
        );
        if let Some(at) = metadata {
            written.remove(at);
        }
        let expected: Vec<String> = case
            .files
            .iter()
            .map(|(rest, _)| format!("db{rest}"))
            .collect();
        assert_eq!(written, expected, "{name}: files");
        for (rest, bytes) in case.files {
            let actual = std::fs::read(out_dir.join(format!("db{rest}"))).unwrap();
            assert!(
                actual == *bytes,
                "{name}: db{rest}: {} bytes written, {} expected; first difference at byte {:?}",
                actual.len(),
                bytes.len(),
                actual.iter().zip(*bytes).position(|(a, b)| a != b)
            );
        }
    }
}

/// One run of `losat makeblastdb` for a protein database titled "t": the exit
/// status, the log and the messages of the standard error stream.
fn makedb(input: &str, out: &str, max_file_sz: &str, stdin: &[u8]) -> (i32, String, String) {
    let args = MakedbArgs {
        input: input.to_string(),
        input_type: "fasta".to_string(),
        dbtype: "prot".to_string(),
        title: Some("t".to_string()),
        out: Some(out.to_string()),
        blastdb_version: 4,
        max_file_sz: max_file_sz.to_string(),
        logfile: None,
        parse_seqids: false,
        hash_index: false,
        taxid: None,
        taxid_map: None,
    };
    let now = time_of("Oct 12, 2026  9:00 PM");
    let (mut log, mut messages) = (Vec::new(), Vec::new());
    let mut stdin = stdin;
    let status = run_with(
        &args,
        Some(&now),
        LINUX,
        &mut stdin,
        &mut log,
        &mut messages,
    );
    (
        status,
        String::from_utf8_lossy(&log).into_owned(),
        String::from_utf8_lossy(&messages).into_owned(),
    )
}

fn prot_fasta() -> String {
    seqdb_fixtures()
        .join("prot.faa")
        .to_string_lossy()
        .into_owned()
}

/// The name of the database as the log gives it ("New DB name:   ...").
fn logged_db_name(log: &str) -> &str {
    log.lines()
        .find_map(|line| line.strip_prefix("New DB name:   "))
        .unwrap_or_else(|| panic!("no database name in {log:?}"))
}

/// makeblastdb splits the output name at every space and refuses a name of
/// more than one part, an empty part included (2.17.0: `-out "db "`,
/// `-out " db"`, `-out "a  b"` and `-out "my data/db"` all end like this).
#[test]
fn an_output_name_with_a_space_is_refused() {
    let fasta = prot_fasta();
    let tmp = TempDir::new("space");
    for out in [
        format!("{} ", tmp.base("db")),
        format!(" {}", tmp.base("db")),
        format!("{}  b", tmp.base("a")),
        format!("{} data/db", tmp.base("my")),
        " ".to_string(),
    ] {
        let (status, log, messages) = makedb(&fasta, &out, "3GB", b"");
        assert_eq!(status, 1, "{out:?}");
        assert_eq!(
            messages, "BLAST options error: Please provide a database name using -out\n",
            "{out:?}"
        );
        assert_eq!(log, "", "{out:?}");
        assert!(
            std::fs::read_dir(&tmp.0).unwrap().next().is_none(),
            "{out:?}"
        );
    }
}

/// The output path is one database name for LOSAT, also with a space or a
/// double quote in a directory's name: the database that is there is found
/// and removed, and the new one opens. (makeblastdb 2.17.0 run in a directory
/// `my data` leaves the old volumes and the old alias file beside the new
/// volume and ends with exit status 2.)
#[test]
fn a_database_in_a_directory_with_a_space_is_replaced() {
    let fixtures = seqdb_fixtures();
    let date = stored_date(&fixtures.join("db/prot_v4.pin"), SeqType::Protein);
    let fasta = fixtures.join("prot.faa");
    let mut names = vec!["my data"];
    if cfg!(unix) {
        names.push("q\"d");
    }
    for name in names {
        let tmp = TempDir::new("spacedir");
        let dir = tmp.0.join(name);
        std::fs::create_dir(&dir).unwrap();
        let base = dir.join("db").to_string_lossy().into_owned();
        let (result, _, _) = build(&fasta, true, "protein in volumes", 300, &date, &base);
        result.unwrap();
        assert_eq!(files_named(&dir, "db.").len(), 6 * 3 + 1, "{name}");
        // In NCBI's syntax of a list of names the path is not one name.
        assert!(crate::seqdb::SeqDb::open(&base, SeqType::Protein).is_err());
        let old = crate::seqdb::SeqDb::open_one_name(&base, SeqType::Protein).unwrap();
        assert_eq!(old.volumes().len(), 6, "{name}");

        let (result, log, _) = build(&fasta, true, "t", 3_000_000_000, &date, &base);
        result.unwrap();
        assert!(
            log.contains("\nDeleted existing Protein BLAST database named "),
            "{name}: {log:?}"
        );
        assert_eq!(
            files_named(&dir, "db."),
            ["db.phr", "db.pin", "db.psq"],
            "{name}"
        );
        let new = crate::seqdb::SeqDb::open_one_name(&base, SeqType::Protein).unwrap();
        assert_eq!((new.volumes().len(), new.num_seqs()), (1, 10), "{name}");
    }
}

/// The `files` of a metadata file.
fn metadata_files(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap();
    let list = &text[text.find("\"files\": [").unwrap()..];
    list.lines()
        .skip(1)
        .take_while(|line| !line.trim_start().starts_with(']'))
        .map(|line| {
            line.trim()
                .trim_end_matches(',')
                .trim_matches('"')
                .to_string()
        })
        .collect()
}

/// makeblastdb writes `NAME.pjs` after every build. The removal of the old
/// database takes `NAME.pjs` only when the old database has one volume, so
/// the file of a database of several volumes is still there when the new one
/// is listed: makeblastdb 2.17.0 (`-max_file_sz 300`, then a single volume)
/// lists `db.pjs` among the new database's files. A JSON file of the other
/// molecule type is not this database's and stays.
#[test]
fn the_json_file_lists_the_files_of_the_database() {
    let fasta = prot_fasta();
    for second_size in ["3GB", "400"] {
        let tmp = TempDir::new("json");
        let base = tmp.base("db");
        let (status, _, messages) = makedb(&fasta, &base, "300", b"");
        assert_eq!((status, messages.as_str()), (0, ""));
        assert_eq!(files_named(&tmp.0, "db.").len(), 6 * 3 + 2);
        let mut volumes: Vec<String> = files_named(&tmp.0, "db.");
        volumes.retain(|file| file != "db.pjs");
        assert_eq!(metadata_files(&tmp.0.join("db.pjs")), volumes);
        std::fs::write(tmp.0.join("db.njs"), b"{}").unwrap();

        let (status, _, messages) = makedb(&fasta, &base, second_size, b"");
        assert_eq!((status, messages.as_str()), (0, ""), "{second_size}");
        let files = files_named(&tmp.0, "db.");
        assert!(files.iter().any(|file| file == "db.njs"), "{files:?}");
        if second_size == "3GB" {
            assert_eq!(files, ["db.njs", "db.phr", "db.pin", "db.pjs", "db.psq"]);
            assert_eq!(
                metadata_files(&tmp.0.join("db.pjs")),
                ["db.phr", "db.pin", "db.pjs", "db.psq"]
            );
        } else {
            let mut volumes = files.clone();
            volumes.retain(|file| file != "db.pjs" && file != "db.njs");
            assert_eq!(metadata_files(&tmp.0.join("db.pjs")), volumes);
        }
    }
}

/// After the build makeblastdb opens the database it wrote. An alias file
/// whose volumes are missing cannot be resolved, so it is not removed before
/// the build and hides the new volume after it. makeblastdb 2.17.0 with a
/// `db.pal` of "DBLIST db.00 db.01": the three files are written, the alias
/// file and the old JSON file stay, the message below, exit status 2.
#[test]
fn a_new_database_that_cannot_be_opened_is_an_error() {
    let fasta = prot_fasta();
    let tmp = TempDir::new("stale-alias");
    let base = tmp.base("db");
    std::fs::write(tmp.0.join("db.pal"), b"TITLE x\nDBLIST db.00 db.01\n").unwrap();
    std::fs::write(tmp.0.join("db.pjs"), b"{}").unwrap();
    let (status, log, messages) = makedb(&fasta, &base, "3GB", b"");
    assert_eq!(status, 2);
    assert_eq!(
        messages,
        format!(
            "BLAST Database error: Could not find volume or alias file (db.00) referenced in \
             alias file ({}).\n",
            logged_db_name(&log)
        )
    );
    assert!(!log.contains("Deleted existing"), "{log:?}");
    assert!(log.ends_with(" seconds.\n\n\n"), "{log:?}");
    assert_eq!(
        files_named(&tmp.0, "db."),
        ["db.pal", "db.phr", "db.pin", "db.pjs", "db.psq"]
    );
}

/// A file of a volume that cannot be created ends the build (LOSAT's own
/// error: makeblastdb does not test its output files), and the files of that
/// volume that were created before it are removed.
#[test]
fn a_volume_that_cannot_be_created_leaves_no_file() {
    let fasta = prot_fasta();
    let tmp = TempDir::new("planted");
    let base = tmp.base("db");
    std::fs::create_dir(tmp.0.join("db.00.phr")).unwrap();
    let (status, log, messages) = makedb(&fasta, &base, "3GB", b"");
    assert_eq!(status, 255);
    assert!(
        messages.starts_with("Error: ") && messages.contains("db.00.phr: "),
        "{messages:?}"
    );
    assert_eq!(messages.lines().count(), 1, "{messages:?}");
    assert!(
        log.ends_with("\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n"),
        "{log:?}"
    );
    assert_eq!(files_named(&tmp.0, "db."), ["db.00.phr"]);
    assert!(tmp.0.join("db.00.phr").is_dir());
}

/// `-in -` reads the standard input: the same database as from the file.
#[test]
fn the_standard_input_gives_the_database_of_the_file() {
    let fasta = prot_fasta();
    let bytes = std::fs::read(&fasta).unwrap();
    for max_file_sz in ["3GB", "300"] {
        let tmp = TempDir::new("stdin");
        let (from_file, from_stdin) = (tmp.0.join("file"), tmp.0.join("stdin"));
        for dir in [&from_file, &from_stdin] {
            std::fs::create_dir(dir).unwrap();
        }
        let base = |dir: &Path| dir.join("db").to_string_lossy().into_owned();
        let (status, log, messages) = makedb(&fasta, &base(&from_file), max_file_sz, b"");
        assert_eq!((status, messages.as_str()), (0, ""));
        // The input of the file run is not read from the standard input.
        assert!(
            log.contains("\nAdding sequences from FASTA; added 10 sequences in "),
            "{log:?}"
        );
        let (status, log, messages) = makedb("-", &base(&from_stdin), max_file_sz, &bytes);
        assert_eq!((status, messages.as_str()), (0, ""));
        assert!(
            log.contains("\nAdding sequences from FASTA; added 10 sequences in "),
            "{log:?}"
        );
        let files = files_named(&from_file, "db.");
        assert_eq!(files.len(), if max_file_sz == "300" { 20 } else { 4 });
        assert_eq!(files_named(&from_stdin, "db."), files);
        for file in &files {
            assert!(
                std::fs::read(from_file.join(file)).unwrap()
                    == std::fs::read(from_stdin.join(file)).unwrap(),
                "{max_file_sz}: {file}"
            );
        }
    }
}

/// An error of the input after volumes were written removes them all. The
/// texts are makeblastdb's (2.17.0, `-max_file_sz 300`, the ten sequences of
/// `prot.faa` followed by ">bad" and "MK?V@L."): the log ends with the two
/// blank lines of the end of the build, without "No volumes were created.".
#[test]
fn a_late_error_removes_the_volumes_already_written() {
    let mut input = std::fs::read(prot_fasta()).unwrap();
    // Without the last record the same input makes six volumes.
    let tmp = TempDir::new("late");
    let fasta = tmp.base("in.fa");
    std::fs::write(&fasta, &input).unwrap();
    let (status, _, _) = makedb(&fasta, &tmp.base("db"), "300", b"");
    assert_eq!(status, 0);
    assert_eq!(files_named(&tmp.0, "db.").len(), 6 * 3 + 2);
    for file in files_named(&tmp.0, "db.") {
        std::fs::remove_file(tmp.0.join(file)).unwrap();
    }

    input.extend_from_slice(b">bad\nMK?V@L.\n");
    std::fs::write(&fasta, &input).unwrap();
    let (status, log, messages) = makedb(&fasta, &tmp.base("db"), "300", b"");
    assert_eq!(status, 1);
    assert_eq!(
        messages,
        "BLAST Database creation error: Near line 30, there's a line that doesn't look like \
         plausible data, but it's not marked as defline or comment.\n"
    );
    assert_eq!(
        String::from_utf8(log_with_markers(log.as_bytes())).unwrap(),
        "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <OUT>\nNew DB title:  t\n\
         Sequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\n\n\n"
    );
    assert_eq!(files_named(&tmp.0, "db."), [] as [&str; 0]);
}

/// The errors makeblastdb reports through its last `catch` clauses end with
/// exit status 255 and "Error: ": a `-max_file_sz` that is no size (LOSAT
/// tests it before it touches a file), and a directory that cannot be used.
#[test]
fn other_errors_end_with_status_255() {
    let fasta = prot_fasta();
    let tmp = TempDir::new("other");
    for value in [
        "abc",
        "",
        "-1",
        "1e3",
        " 1KB",
        "1Ki",
        "18446744073709551616",
    ] {
        let (status, log, messages) = makedb(&fasta, &tmp.base("db"), value, b"");
        assert_eq!(status, 255, "{value:?}");
        assert_eq!(
            messages,
            format!("Error: -max_file_sz: Cannot convert string '{value}' to Uint8\n")
        );
        assert_eq!(log, "", "{value:?}");
    }
    // A value that is a size, but too large, is an error of the options.
    let (status, _, messages) = makedb(&fasta, &tmp.base("db"), "4GiB", b"");
    assert_eq!(
        (status, messages.as_str()),
        (1, "BLAST options error: max_file_sz must be < 4 GiB\n")
    );
    // The parent of the database is a file.
    std::fs::write(tmp.0.join("f"), b"").unwrap();
    let out = tmp.0.join("f").join("db").to_string_lossy().into_owned();
    let (status, log, messages) = makedb(&fasta, &out, "3GB", b"");
    assert_eq!(
        (status, log.as_str(), messages.as_str()),
        (255, "", "Error: Failed to create directory 'f'\n")
    );
    assert_eq!(
        files_named(&tmp.0, ""),
        ["f"],
        "nothing but the planted file"
    );
}

/// The directory of the output name is made as makeblastdb makes it
/// (2.17.0): `-out sub/` is the database `sub` and makes no directory;
/// missing directories are made; and for `-out nodir/../db` makeblastdb makes
/// the normalized path (nothing: it exists) and then fails on the path as it
/// was given, which does not resolve, before it writes a line of the log.
#[test]
fn directories_of_the_output_name() {
    let fasta = prot_fasta();
    let tmp = TempDir::new("dirs");
    let sep = std::path::MAIN_SEPARATOR;

    let (status, _, messages) = makedb(&fasta, &format!("{}{sep}", tmp.base("sub")), "3GB", b"");
    assert_eq!((status, messages.as_str()), (0, ""));
    assert_eq!(
        files_named(&tmp.0, "sub"),
        ["sub.phr", "sub.pin", "sub.pjs", "sub.psq"]
    );

    let nested = tmp.0.join("a").join("b");
    let out = nested.join("db").to_string_lossy().into_owned();
    let (status, _, messages) = makedb(&fasta, &out, "3GB", b"");
    assert_eq!((status, messages.as_str()), (0, ""));
    assert_eq!(
        files_named(&nested, "db."),
        ["db.phr", "db.pin", "db.pjs", "db.psq"]
    );
    #[cfg(unix)]
    {
        // The directories get the mode rwxr-xr-x whatever the umask.
        use std::os::unix::fs::PermissionsExt;
        for dir in [tmp.0.join("a"), nested.clone()] {
            let mode = std::fs::metadata(&dir).unwrap().permissions().mode();
            assert_eq!(mode & 0o7777, 0o755, "{}", dir.display());
        }
    }

    // On Windows the system resolves `nodir\..` without looking at `nodir`.
    #[cfg(unix)]
    {
        for (given, name, made, not_made) in [
            ("nodir/../db", "..", None, "nodir"),
            ("x/b/../c/db", "c", Some("x/c"), "x/b"),
        ] {
            let out = format!("{}/{given}", tmp.0.display());
            let (status, log, messages) = makedb(&fasta, &out, "3GB", b"");
            assert_eq!(
                (status, log.as_str(), messages),
                (
                    255,
                    "",
                    format!("Error: You do not have write permissions on '{name}'\n")
                ),
                "{given}"
            );
            if let Some(made) = made {
                assert!(tmp.0.join(made).is_dir(), "{given}");
                assert!(files_named(&tmp.0.join(made), "").is_empty(), "{given}");
            }
            assert!(!tmp.0.join(not_made).exists(), "{given}");
        }
        assert!(files_named(&tmp.0, "db.").is_empty());
    }
}

/// More than 100 volumes are renamed with leading zeros when the database
/// is closed. The names are makeblastdb's (2.17.0, 1,001 volumes): the
/// volumes below 100 get four digits, those from 100 to 999 keep their three
/// digits behind two zeros, and those from 1,000 keep their name.
#[test]
fn more_than_a_hundred_volumes_get_makeblastdb_s_names() {
    let tmp = TempDir::new("volumes");
    let mut db = WriteDb::new(
        &tmp.base("db"),
        true,
        b"t",
        b"Oct 12, 2026  9:00 PM",
        1,
        DbVersion::V4,
        LINUX,
        Indices::default(),
    );
    for i in 0..1001 {
        let bioseq = Bioseq {
            local_id: (i + 1).to_string(),
            title: format!("p{i}").into_bytes(),
            is_aa: true,
            length: 2,
            seq_data: Some(SeqData::Ncbieaa(b"MK".to_vec())),
            parsed: None,
        };
        db.extract_bioseq_deflines(&bioseq);
        db.add_sequence(bioseq).unwrap();
    }
    db.close().unwrap();
    let files = files_named(&tmp.0, "db.");
    assert_eq!(files.len(), 3 * 1001 + 1);
    for name in [
        "db.0000.pin",
        "db.0099.psq",
        "db.00100.phr",
        "db.00999.pin",
        "db.1000.psq",
        "db.pal",
    ] {
        assert!(files.iter().any(|file| file == name), "{name}");
    }
    let alias = std::fs::read_to_string(tmp.0.join("db.pal")).unwrap();
    let list = alias
        .lines()
        .find_map(|line| line.strip_prefix("DBLIST "))
        .unwrap();
    let volumes: Vec<&str> = list.split(' ').collect();
    assert_eq!(volumes.len(), 1001);
    assert_eq!(volumes[..2], ["db.0000", "db.0001"]);
    assert_eq!(volumes[99..101], ["db.0099", "db.00100"]);
    assert_eq!(volumes[999..], ["db.00999", "db.1000"]);
}

/// The parsed `losat makeblastdb` arguments of a makeblastdb command line (the
/// words after the program's name).
fn parse_makedb_words(words: &[String]) -> MakedbArgs {
    use crate::cli::{Cli, Commands};
    let mut argv = vec!["losat".to_string(), "makeblastdb".to_string()];
    argv.extend(words.iter().cloned());
    match crate::cli::try_parse_from::<Cli, _, _>(argv) {
        Ok(cli) => match cli.command {
            Commands::Makedb(args) => args,
            _ => unreachable!("a makeblastdb command line"),
        },
        Err(error) => panic!("{words:?}: {error}"),
    }
}

/// FNV-1a, 64 bits.
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// A stream or log of a run with the markers of `makeblastdb_v5_cases.rs`: the output
/// directory `dir` as `<D>`, the build time as `<T>`, the seconds of the timed line as `<S>`.
/// On Windows the directory's `\` are written as `/`, as in the recorded Linux runs.
fn v5_markers(bytes: &[u8], dir: &str) -> Vec<u8> {
    let mut text = replace_bytes(bytes, dir.as_bytes(), b"<D>");
    if cfg!(windows) {
        text = replace_bytes(&text, b"\\", b"/");
        text = replace_bytes(&text, dir.replace('\\', "/").as_bytes(), b"<D>");
    }
    let lines: Vec<Vec<u8>> = text
        .split(|&b| b == b'\n')
        .map(|line| {
            let line_text = String::from_utf8_lossy(line);
            const BUILDING: &str = "Building a new DB, current time: ";
            if line_text.starts_with(BUILDING) {
                format!("{BUILDING}<T>").into_bytes()
            } else if let (Some(at), true) = (
                line_text.find(" sequences in "),
                line_text.ends_with(" seconds."),
            ) {
                format!("{} sequences in <S> seconds.", &line_text[..at]).into_bytes()
            } else {
                line.to_vec()
            }
        })
        .collect();
    lines.join(&b'\n')
}

/// Replays a case of `makeblastdb_v5_cases.rs` through the parser and `run_with` in a directory
/// of its own: every run's exit status and streams, the log file and every file are
/// makeblastdb's (the date of the database is the recorded one).
fn run_v5_case(case: &V5Case) {
    let name = case.name;
    let tmp = TempDir::new(&format!("v5-{name}"));
    let dir = tmp.0.to_string_lossy().into_owned();
    let now = time_of(if case.date.is_empty() {
        "Oct 12, 2026  9:00 PM"
    } else {
        case.date
    });
    let show = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    for (k, run) in case.runs.iter().enumerate() {
        let words: Vec<String> = run
            .words
            .iter()
            .map(|word| match word.replace("{D}", &dir) {
                // makeblastdb's null device on Windows is `NUL`; `/dev/null` names it on Unix.
                null if cfg!(windows) && null == "/dev/null" => "NUL".to_string(),
                word => word,
            })
            .collect();
        let args = parse_makedb_words(&words);
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = run_with(
            &args,
            Some(&now),
            LINUX,
            &mut std::io::empty(),
            &mut out,
            &mut err,
        );
        let (out, err) = (v5_markers(&out, &dir), v5_markers(&err, &dir));
        assert!(
            err == run.stderr,
            "{name} run {k}: standard error\n{}\nexpected\n{}",
            show(&err),
            show(run.stderr)
        );
        assert!(
            out == run.stdout,
            "{name} run {k}: standard output\n{}\nexpected\n{}",
            show(&out),
            show(run.stdout)
        );
        assert_eq!(status, run.status, "{name} run {k}: exit status");
    }
    let log = std::fs::read(tmp.0.join("mk.log"))
        .ok()
        .map(|log| v5_markers(&log, &dir));
    assert!(
        log.as_deref() == case.log,
        "{name}: log file\n{:?}\nexpected\n{:?}",
        log.as_deref().map(show),
        case.log.map(show)
    );
    let mut written: Vec<String> = std::fs::read_dir(&tmp.0)
        .unwrap()
        .map(|entry| entry.unwrap())
        .filter(|entry| entry.file_type().unwrap().is_file())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|file| file != "mk.log")
        .collect();
    written.sort();
    let expected: Vec<&str> = case.files.iter().map(|(file, _)| *file).collect();
    assert_eq!(written, expected, "{name}: files");
    for (file, content) in case.files {
        let actual = std::fs::read(tmp.0.join(file)).unwrap();
        match content {
            Expected::Bytes(bytes) => assert!(
                actual == *bytes,
                "{name}: {file}: {} bytes written, {} expected; first difference at byte {:?}",
                actual.len(),
                bytes.len(),
                actual.iter().zip(*bytes).position(|(a, b)| a != b)
            ),
            Expected::Digest(len, digest) => assert_eq!(
                (actual.len(), fnv1a(&actual)),
                (*len, *digest),
                "{name}: {file}: size and FNV-1a"
            ),
        }
    }
}

fn run_v5_cases(prefix: &str) {
    let cases: Vec<&V5Case> = V5_CASES
        .iter()
        .filter(|case| case.name.starts_with(prefix))
        .collect();
    assert!(!cases.is_empty(), "no case {prefix}*");
    for case in cases {
        run_v5_case(case);
    }
}

/// `-logfile` as makeblastdb 2.17.0 has it: a new file and two runs appending to one (the log
/// and the messages, nothing on the standard streams), `/dev/null` (`NUL` on Windows), `-` (the
/// standard streams, as without `-logfile`), a warning and an error that go to the file, the
/// empty name (exit status 1, nothing shown) and a file in a missing directory.
#[test]
fn logfile_runs_as_makeblastdb() {
    run_v5_cases("L_");
}

/// Several files in one `-in` as makeblastdb 2.17.0 has them: the default title is the whole
/// value, a quoted name, a record without data in the second file, a missing and an empty
/// second file.
#[test]
fn several_input_files_run_as_makeblastdb() {
    run_v5_cases("M_");
}

/// The metadata file of version 5 (every field and the list of files, the LMDB and lookup
/// files included) for protein and nucleotide databases of one volume and of several.
#[test]
fn the_metadata_files_of_version_5_are_makeblastdb_s() {
    run_v5_cases("J_");
}

/// One of several `-in` files that exists and cannot be opened: LOSAT stops the build with
/// exit status 255 and removes what the earlier files added, as for a late error (makeblastdb
/// reads nothing from it and goes on; PD-LOSAT-DATABASE-SEARCH, DECISIONS 20 of the makedb-v5
/// line).
#[cfg(unix)]
#[test]
fn an_input_file_that_cannot_be_opened_stops_the_build() {
    use std::os::unix::fs::PermissionsExt;
    let tmp = TempDir::new("unreadable");
    let locked = tmp.base("locked.fa");
    std::fs::write(&locked, b">l\nMKV\n").unwrap();
    std::fs::set_permissions(&locked, std::fs::Permissions::from_mode(0o000)).unwrap();
    if std::fs::File::open(&locked).is_ok() {
        // A user who reads every file (root): there is no file that cannot be opened.
        return;
    }
    let base = tmp.base("db");
    for version in ["4", "5"] {
        let words: Vec<String> = [
            "-in",
            &format!("tests/writedb_v5/a.fa {locked}"),
            "-dbtype",
            "prot",
            "-out",
            &base,
            "-title",
            "t",
            "-blastdb_version",
            version,
        ]
        .iter()
        .map(|word| word.to_string())
        .collect();
        let args = parse_makedb_words(&words);
        let now = time_of("Oct 12, 2026  9:00 PM");
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = run_with(
            &args,
            Some(&now),
            LINUX,
            &mut std::io::empty(),
            &mut out,
            &mut err,
        );
        let err = String::from_utf8_lossy(&err).into_owned();
        let out = String::from_utf8_lossy(&out).into_owned();
        assert_eq!(status, 255, "v{version}: {err}");
        assert!(
            err.starts_with(&format!(
                "Warning: [makeblastdb] \nInput file {locked} does NOT appear to be FASTA \
                 (processing anyway).\n"
            )),
            "v{version}: {err}"
        );
        assert!(
            err.ends_with("\n") && err.contains(&format!("\nError: {locked}: ")),
            "v{version}: {err}"
        );
        // The first file was read into the database before the second failed.
        assert!(
            out.contains("Adding sequences from FASTA; added 2 sequences in "),
            "v{version}: {out}"
        );
        // The volumes and the LMDB file are removed. Of version 5 the tax-id lookup files and
        // the `.ptf` stay, as after makeblastdb's own late errors (`ListFiles` names the `.pdb`
        // only; the late-error cases E_late and E_late1 of the makedb-v5 S2 gate).
        let left: Vec<String> = files_named(&tmp.0, "db.");
        let expected: &[&str] = if version == "4" {
            &[]
        } else {
            &["db.pot", "db.ptf", "db.pto"]
        };
        assert_eq!(left, expected, "v{version}");
    }
}

/// The driver of the makedb-v5 line's S2 gate
/// (`makedb-v5-20261029/scripts/s2_gate.py`), not a test of its own: for each
/// line of the manifest named by `LOSAT_S2_GATE_MANIFEST` (a directory, the
/// build date makeblastdb stored, then one JSON array of words per run) it runs
/// `losat makeblastdb` through `run_with` in that directory's terms, with that date,
/// and writes `so<k>.txt`, `se<k>.txt` and `rc<k>.txt` there. The script
/// compares every file with makeblastdb's.
#[test]
#[ignore = "driver of an external gate (LOSAT_S2_GATE_MANIFEST)"]
fn s2_gate() {
    let manifest = std::env::var("LOSAT_S2_GATE_MANIFEST").expect("LOSAT_S2_GATE_MANIFEST");
    for line in std::fs::read_to_string(manifest).unwrap().lines() {
        let mut fields = line.split('\t');
        let dir = PathBuf::from(fields.next().unwrap());
        let date = fields.next().unwrap();
        let now = time_of(if date.is_empty() {
            "Oct 12, 2026  9:00 PM"
        } else {
            date
        });
        for (k, run) in fields.enumerate() {
            let words: Vec<String> = parse_json_words(run);
            let args = parse_makedb_words(&words);
            let (mut out, mut err) = (Vec::new(), Vec::new());
            let status = run_with(
                &args,
                Some(&now),
                LINUX,
                &mut std::io::empty(),
                &mut out,
                &mut err,
            );
            std::fs::write(dir.join(format!("so{k}.txt")), out).unwrap();
            std::fs::write(dir.join(format!("se{k}.txt")), err).unwrap();
            std::fs::write(dir.join(format!("rc{k}.txt")), status.to_string()).unwrap();
        }
    }
}

/// The driver of the S3a gate of the makedb-v5 line
/// (`makedb-v5-20261029/scripts/s3a_gate.py`): as `s2_gate`, for the manifest named by
/// `LOSAT_S3A_GATE_MANIFEST`, but the word `-parse_seqids`, which the command line
/// refuses, sets `MakedbArgs::parse_seqids` instead (the library-only path of S3a).
#[test]
#[ignore = "driver of an external gate (LOSAT_S3A_GATE_MANIFEST)"]
fn s3a_gate() {
    let manifest = std::env::var("LOSAT_S3A_GATE_MANIFEST").expect("LOSAT_S3A_GATE_MANIFEST");
    for line in std::fs::read_to_string(manifest).unwrap().lines() {
        let mut fields = line.split('\t');
        let dir = PathBuf::from(fields.next().unwrap());
        let date = fields.next().unwrap();
        let now = time_of(if date.is_empty() {
            "Oct 12, 2026  9:00 PM"
        } else {
            date
        });
        for (k, run) in fields.enumerate() {
            let words: Vec<String> = parse_json_words(run);
            let parse_seqids = words.iter().any(|word| word == "-parse_seqids");
            let words: Vec<String> = words
                .into_iter()
                .filter(|word| word != "-parse_seqids")
                .collect();
            let mut args = parse_makedb_words(&words);
            args.parse_seqids = parse_seqids;
            let (mut out, mut err) = (Vec::new(), Vec::new());
            let status = run_with(
                &args,
                Some(&now),
                LINUX,
                &mut std::io::empty(),
                &mut out,
                &mut err,
            );
            std::fs::write(dir.join(format!("so{k}.txt")), out).unwrap();
            std::fs::write(dir.join(format!("se{k}.txt")), err).unwrap();
            std::fs::write(dir.join(format!("rc{k}.txt")), status.to_string()).unwrap();
        }
    }
}

/// The driver of the S3b gate of the makedb-v5 line
/// (`makedb-v5-20261029/scripts/s3b_gate.py`): as `s3a_gate`, for the manifest named by
/// `LOSAT_S3B_GATE_MANIFEST`. The words go through LOSAT's command line parser (a parser
/// error is written as `main` writes it, with its exit code) except the two that the
/// command line refuses until S3c: `-parse_seqids` sets `MakedbArgs::parse_seqids` and
/// `-taxid_map <file>` sets `MakedbArgs::taxid_map`.
#[test]
#[ignore = "driver of an external gate (LOSAT_S3B_GATE_MANIFEST)"]
fn s3b_gate() {
    use crate::cli::{Cli, Commands};
    let manifest = std::env::var("LOSAT_S3B_GATE_MANIFEST").expect("LOSAT_S3B_GATE_MANIFEST");
    for line in std::fs::read_to_string(manifest).unwrap().lines() {
        let mut fields = line.split('\t');
        let dir = PathBuf::from(fields.next().unwrap());
        let date = fields.next().unwrap();
        let now = time_of(if date.is_empty() {
            "Oct 12, 2026  9:00 PM"
        } else {
            date
        });
        for (k, run) in fields.enumerate() {
            let words: Vec<String> = parse_json_words(run);
            let mut parse_seqids = false;
            let mut taxid_map = None;
            let mut rest = Vec::new();
            let mut it = words.into_iter();
            while let Some(word) = it.next() {
                match word.as_str() {
                    "-parse_seqids" => parse_seqids = true,
                    "-taxid_map" => taxid_map = it.next(),
                    _ => rest.push(word),
                }
            }
            let mut argv = vec!["losat".to_string(), "makeblastdb".to_string()];
            argv.extend(rest);
            let (mut out, mut err) = (Vec::new(), Vec::new());
            let status = match crate::cli::try_parse_from::<Cli, _, _>(argv) {
                Ok(cli) => {
                    let Commands::Makedb(mut args) = cli.command else {
                        unreachable!("a makeblastdb command line");
                    };
                    args.parse_seqids = parse_seqids;
                    args.taxid_map = taxid_map;
                    run_with(
                        &args,
                        Some(&now),
                        LINUX,
                        &mut std::io::empty(),
                        &mut out,
                        &mut err,
                    )
                }
                Err(error) => {
                    let message = crate::cli::render_message(&error);
                    if error.use_stderr() {
                        err.extend_from_slice(message.as_bytes());
                    } else {
                        out.extend_from_slice(message.as_bytes());
                    }
                    error.exit_code()
                }
            };
            std::fs::write(dir.join(format!("so{k}.txt")), out).unwrap();
            std::fs::write(dir.join(format!("se{k}.txt")), err).unwrap();
            std::fs::write(dir.join(format!("rc{k}.txt")), status.to_string()).unwrap();
        }
    }
}

/// A JSON array of strings as the gate script writes it (`json.dumps` with
/// `ensure_ascii`: the escapes `\"`, `\\`, `\uXXXX` and the one-letter ones).
fn parse_json_words(text: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut chars = text
        .trim()
        .trim_start_matches('[')
        .trim_end_matches(']')
        .chars();
    let mut current: Option<String> = None;
    while let Some(c) = chars.next() {
        match (&mut current, c) {
            (None, '"') => current = Some(String::new()),
            (None, _) => {}
            (Some(word), '"') => {
                words.push(std::mem::take(word));
                current = None;
            }
            (Some(word), '\\') => match chars.next().unwrap() {
                'n' => word.push('\n'),
                't' => word.push('\t'),
                'r' => word.push('\r'),
                'u' => {
                    let hex: String = chars.by_ref().take(4).collect();
                    word.push(char::from_u32(u32::from_str_radix(&hex, 16).unwrap()).unwrap());
                }
                other => word.push(other),
            },
            (Some(word), c) => word.push(c),
        }
    }
    words
}
