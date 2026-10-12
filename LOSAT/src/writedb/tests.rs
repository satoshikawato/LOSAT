//! Tests of the database writer against databases written by makeblastdb
//! 2.17.0 (`LOSAT/tests/seqdb`, `LOSAT/tests/writedb`). The tests read the
//! committed files and run no NCBI program.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use super::bioseq::{Bioseq, SeqData};
use super::build_db::BuildDatabase;
use super::time::LocalTime;
use super::writer::WriteDb;
use crate::seqdb::SeqType;

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

/// Every file `out_dir/name.*` equals the fixture's file of the same
/// extension, and the fixture has no other file LOSAT writes.
fn assert_same_database(out_dir: &Path, name: &str, fixture_dir: &Path, fixture: &str) {
    let written = files_named(out_dir, &format!("{name}."));
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
    let losat_writes = |file: &String| {
        ["in", "sq", "hr", "al"]
            .iter()
            .any(|ext| file.ends_with(&format!(".p{ext}")) || file.ends_with(&format!(".n{ext}")))
    };
    let expected: Vec<String> = files_named(fixture_dir, &format!("{fixture}."))
        .iter()
        .filter(|file| losat_writes(file))
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
        })
        .collect()
}

/// The writer alone, from sequences built by hand: the single-volume and the
/// six-volume protein databases of `LOSAT/tests/seqdb`.
#[test]
fn the_writer_reproduces_makeblastdb_s_protein_volumes() {
    let fixtures = seqdb_fixtures();
    let cases = [
        (
            "prot_v4",
            "prot_v4.pin",
            "LOSAT seqdb fixture, protein, version 4",
            3_000_000_000u64,
        ),
        ("prot_mv", "prot_mv.00.pin", "protein in volumes", 300),
    ];
    for (fixture, index_file, title, max_file_size) in cases {
        let date = stored_date(&fixtures.join("db").join(index_file), SeqType::Protein);
        let tmp = TempDir::new(fixture);
        let mut db = WriteDb::new(
            &tmp.base(fixture),
            true,
            title.as_bytes(),
            date.as_bytes(),
            max_file_size,
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
    let mut log = Vec::new();
    let mut messages = Vec::new();
    let result = (|| {
        let mut db = BuildDatabase::new(out_base, title, is_protein, &mut log, &time_of(date))?;
        db.set_max_file_size(max_file_size);
        let mut warn = |line: &[u8]| {
            messages.extend_from_slice(line);
            messages.push(b'\n');
        };
        let added = db.add_fasta(std::fs::File::open(fasta).unwrap(), &mut warn);
        messages.extend_from_slice(&db.take_diagnostics());
        added?;
        let ended = db.end_build(false);
        messages.extend_from_slice(&db.take_diagnostics());
        ended.map(|_| ())
    })();
    (result, String::from_utf8_lossy(&log).into_owned(), messages)
}

/// The databases of `LOSAT/tests/seqdb` that makeblastdb wrote in format
/// version 4, from their FASTA files.
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
        ),
        (
            "nucl_v4",
            "nucl.fna",
            false,
            "nucl_v4.nin",
            "LOSAT seqdb fixture, nucleotide, version 4",
            3_000_000_000,
        ),
        (
            "prot_mv",
            "prot.faa",
            true,
            "prot_mv.00.pin",
            "protein in volumes",
            300,
        ),
    ];
    for (fixture, fasta, is_protein, index_file, title, max_file_size) in cases {
        let seqtype = if is_protein {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        let date = stored_date(&fixtures.join("db").join(index_file), seqtype);
        let tmp = TempDir::new(fixture);
        let (result, log, messages) = build(
            &fixtures.join(fasta),
            is_protein,
            title,
            max_file_size,
            &date,
            &tmp.base(fixture),
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
