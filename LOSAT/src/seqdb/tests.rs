//! Tests of the database reader.
//!
//! The fixtures under `LOSAT/tests/seqdb` were written by NCBI BLAST+ 2.17.0
//! (`make_fixtures.py`): the databases by makeblastdb, `expected/*.tsv` by
//! blastdbcmd. The tests compare what LOSAT reads with the FASTA files the
//! databases were built from and with what NCBI reads back.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use super::alias::{test_read_values, test_split_quoted};
use super::header::decode_defline_set;
use super::index::IndexFile;
use super::*;
use crate::core::blast_encoding::{encode_iupac_to_blastna, iupacna_to_ncbi4na};

/// NCBIstdaa: the residue at position `i` has the code `i`
/// (c++/src/algo/blast/core/blast_encoding.c:115-118, `NCBISTDAA_TO_AMINOACID`).
const NCBISTDAA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb")
}

fn db_name(name: &str) -> String {
    fixture_dir()
        .join("db")
        .join(name)
        .to_str()
        .expect("fixture path is UTF-8")
        .to_string()
}

struct Record {
    title: Vec<u8>,
    seq: Vec<u8>,
}

fn read_fasta(path: &Path) -> Vec<Record> {
    let mut records: Vec<Record> = Vec::new();
    let file = std::io::BufReader::new(std::fs::File::open(path).expect("open FASTA"));
    for line in file.split(b'\n') {
        let mut line = line.expect("read FASTA");
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if let Some(title) = line.strip_prefix(b">") {
            records.push(Record {
                title: title.to_vec(),
                seq: Vec::new(),
            });
        } else if let Some(record) = records.last_mut() {
            record
                .seq
                .extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    records
}

/// One line of `expected/NAME.entries.tsv`: blastdbcmd `%o %l %i %s %t`.
struct Entry {
    ordinal: u32,
    length: usize,
    id: Vec<u8>,
    seq: Vec<u8>,
    title: Vec<u8>,
}

fn expected_entries(name: &str) -> Vec<Entry> {
    let path = fixture_dir()
        .join("expected")
        .join(format!("{name}.entries.tsv"));
    let bytes = std::fs::read(&path).expect("read expected entries");
    bytes
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            let cols: Vec<&[u8]> = line.splitn(5, |&b| b == b'\t').collect();
            assert_eq!(cols.len(), 5, "{name}: five columns");
            let number = |col: &[u8]| std::str::from_utf8(col).unwrap().parse::<u64>().unwrap();
            Entry {
                ordinal: number(cols[0]) as u32,
                length: number(cols[1]) as usize,
                id: cols[2].to_vec(),
                seq: cols[3].to_vec(),
                title: cols[4].to_vec(),
            }
        })
        .collect()
}

fn expected_info(name: &str) -> HashMap<String, String> {
    let path = fixture_dir()
        .join("expected")
        .join(format!("{name}.info.tsv"));
    std::fs::read_to_string(&path)
        .expect("read expected info")
        .lines()
        .map(|line| {
            let (key, value) = line.split_once('\t').expect("key and value");
            (key.to_string(), value.to_string())
        })
        .collect()
}

fn protein_codes(letters: &[u8]) -> Vec<u8> {
    letters
        .iter()
        .map(|letter| {
            let upper = letter.to_ascii_uppercase();
            NCBISTDAA
                .iter()
                .position(|&aa| aa == upper)
                .unwrap_or_else(|| panic!("letter {:?} is not in NCBIstdaa", upper as char))
                as u8
        })
        .collect()
}

/// Checks sequence `i` of `block` against the letters of the source.
fn check_sequence(block: &SequenceBlock, i: usize, letters: &[u8], what: &str) {
    match block.seqtype() {
        SeqType::Protein => {
            let codes = protein_codes(letters);
            assert_eq!(block.protein(i), codes.as_slice(), "{what}: residues");
            let with_sentinels = block.protein_with_sentinels(i);
            assert_eq!(with_sentinels.len(), codes.len() + 2, "{what}: sentinels");
            assert_eq!(with_sentinels[0], 0, "{what}: sentinel before");
            assert_eq!(with_sentinels[codes.len() + 1], 0, "{what}: sentinel after");
            assert_eq!(&with_sentinels[1..=codes.len()], codes.as_slice());
        }
        SeqType::Nucleotide => {
            let upper = letters.to_ascii_uppercase();
            assert_eq!(block.nucleotide_len(i), upper.len(), "{what}: length");
            // The sequence with its ambiguity codes, as the `-subject` path encodes it.
            let blastna = block.nucleotide_blastna(i).expect("BLASTNA");
            assert_eq!(blastna.len(), upper.len() + 2, "{what}: BLASTNA length");
            assert_eq!(blastna[0], BLASTNA_SENTINEL);
            assert_eq!(blastna[upper.len() + 1], BLASTNA_SENTINEL);
            assert_eq!(
                &blastna[1..=upper.len()],
                encode_iupac_to_blastna(&upper).as_slice(),
                "{what}: BLASTNA"
            );
            // The 2-bit data: the base itself, or for an ambiguity code one of
            // the bases it stands for (makeblastdb picks one at random).
            let packed = block.nucleotide_packed(i);
            assert_eq!(packed.len(), upper.len() / 4 + 1, "{what}: packed bytes");
            assert_eq!(usize::from(packed[packed.len() - 1] & 3), upper.len() % 4);
            let mut ambiguous = 0;
            for (position, &letter) in upper.iter().enumerate() {
                let code = (packed[position >> 2] >> (6 - 2 * (position & 3))) & 3;
                let allowed = iupacna_to_ncbi4na(letter).expect("IUPAC letter");
                assert!(
                    allowed & (1 << code) != 0,
                    "{what}: base {position} ({}) stored as {code}",
                    letter as char
                );
                if !matches!(letter, b'A' | b'C' | b'G' | b'T') {
                    ambiguous += 1;
                }
            }
            let covered: u32 = block
                .ambiguity_runs(i)
                .expect("ambiguity runs")
                .iter()
                .map(|run| run.length)
                .sum();
            assert_eq!(covered as usize, ambiguous, "{what}: ambiguous bases");
            // A sequence without ambiguous bases normally has no ambiguity
            // data; makeblastdb gives some of them a count of zero (see
            // `both_ambiguity_formats_are_read`).
            if ambiguous > 0 {
                assert!(
                    !block.ambiguity_bytes(i).is_empty(),
                    "{what}: ambiguity data"
                );
            }
        }
    }
}

/// Reads the whole database in blocks of at most `max_bytes` and checks every
/// sequence against `entries`.
fn check_database(db: &SeqDb, entries: &[Entry], max_bytes: u64, name: &str) {
    assert_eq!(db.num_seqs() as usize, entries.len(), "{name}: sequences");
    let mut next = 0u32;
    let mut blocks = 0;
    while let Some(range) = db.next_block(next, max_bytes) {
        assert_eq!(range.first_ordinal, next, "{name}: blocks are contiguous");
        assert!(!range.is_empty(), "{name}: a block holds a sequence");
        let volume = &db.volumes()[range.volume];
        assert!(
            range.end <= volume.num_oids(),
            "{name}: block inside its volume"
        );
        if range.len() > 1 {
            assert!(
                volume.block_bytes(range.first, range.end) <= max_bytes,
                "{name}: block within the budget"
            );
        }
        let block = db.read_block(range).expect("read block");
        assert_eq!(block.len(), range.len() as usize);
        assert_eq!(block.first_ordinal(), next);
        for i in 0..block.len() {
            let entry = &entries[next as usize + i];
            let what = format!("{name} ordinal {}", entry.ordinal);
            assert_eq!(entry.ordinal as usize, next as usize + i);
            assert_eq!(entry.seq.len(), entry.length, "{what}: NCBI length");
            check_sequence(&block, i, &entry.seq, &what);
        }
        next = range.end_ordinal();
        blocks += 1;
    }
    assert_eq!(next, db.num_seqs(), "{name}: every sequence read");
    assert!(db.next_block(next, max_bytes).is_none());
    if max_bytes == 0 {
        assert_eq!(blocks, entries.len(), "{name}: one sequence per block");
    }
}

fn check_deflines_and_single_reads(db: &SeqDb, entries: &[Entry], name: &str) {
    for entry in entries {
        let what = format!("{name} ordinal {}", entry.ordinal);
        let deflines = db.deflines(entry.ordinal).expect("deflines");
        assert_eq!(deflines.len(), 1, "{what}: one definition line");
        assert_eq!(deflines[0].title, entry.title, "{what}: title");
        assert_eq!(deflines[0].bl_ord_id, entry.ordinal, "{what}: BL_ORD_ID");
        assert_eq!(
            entry.id,
            format!("gnl|BL_ORD_ID|{}", entry.ordinal).into_bytes(),
            "{what}: NCBI id"
        );
        let block = db.read_sequence(entry.ordinal).expect("read sequence");
        assert_eq!(block.len(), 1);
        check_sequence(&block, 0, &entry.seq, &what);
    }
    assert!(matches!(
        db.deflines(db.num_seqs()),
        Err(SeqDbError::Range(_))
    ));
    assert!(matches!(
        db.read_sequence(db.num_seqs()),
        Err(SeqDbError::Range(_))
    ));
}

fn check_info(db: &SeqDb, name: &str) {
    let info = expected_info(name);
    assert_eq!(db.title(), info["title"].as_bytes(), "{name}: title");
    assert_eq!(
        db.num_seqs().to_string(),
        info["sequences"],
        "{name}: sequences"
    );
    assert_eq!(
        db.total_length().to_string(),
        info["letters"],
        "{name}: letters"
    );
    assert_eq!(
        db.max_length().to_string(),
        info["longest"],
        "{name}: longest"
    );
    assert_eq!(db.date(), info["date"], "{name}: date");
    for volume in db.volumes() {
        assert_eq!(volume.index().version.to_string(), info["version"]);
    }
}

fn check_against_ncbi(name: &str, seqtype: SeqType) -> (SeqDb, Vec<Entry>) {
    let db = SeqDb::open(&db_name(name), seqtype).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(db.seqtype(), seqtype);
    let entries = expected_entries(name);
    check_info(&db, name);
    for max_bytes in [0, 64, 400, 1 << 20] {
        check_database(&db, &entries, max_bytes, name);
    }
    check_deflines_and_single_reads(&db, &entries, name);
    (db, entries)
}

/// The sequences and titles NCBI reads back are those of the FASTA, in order.
fn check_entries_against_fasta(entries: &[Entry], fasta: &str, name: &str) {
    let records = read_fasta(&fixture_dir().join(fasta));
    assert_eq!(entries.len(), records.len(), "{name}: FASTA records");
    for (entry, record) in entries.iter().zip(&records) {
        assert_eq!(entry.title, record.title, "{name}: FASTA title");
        assert_eq!(
            entry.seq,
            record.seq.to_ascii_uppercase(),
            "{name}: FASTA sequence"
        );
    }
}

fn volume_names(db: &SeqDb) -> Vec<String> {
    db.volumes()
        .iter()
        .map(|v| v.base().file_name().unwrap().to_str().unwrap().to_string())
        .collect()
}

#[test]
fn single_volume_databases_match_ncbi_and_the_fasta() {
    for (name, seqtype, fasta) in [
        ("prot_v4", SeqType::Protein, "prot.faa"),
        ("prot_v5", SeqType::Protein, "prot.faa"),
        ("nucl_v4", SeqType::Nucleotide, "nucl.fna"),
        ("nucl_v5", SeqType::Nucleotide, "nucl.fna"),
    ] {
        let (db, entries) = check_against_ncbi(name, seqtype);
        check_entries_against_fasta(&entries, fasta, name);
        assert_eq!(db.volumes().len(), 1);
        let index = db.volumes()[0].index();
        let version5 = name.ends_with("v5");
        assert_eq!(index.version, if version5 { 5 } else { 4 });
        assert_eq!(index.volume_number.is_some(), version5);
        assert_eq!(index.lmdb_file.is_some(), version5);
    }
}

#[test]
fn version_4_and_5_store_the_same_sequence_and_header_bytes() {
    let dir = fixture_dir().join("db");
    for (a, b) in [
        ("prot_v4.psq", "prot_v5.psq"),
        ("prot_v4.phr", "prot_v5.phr"),
        ("nucl_v4.nhr", "nucl_v5.nhr"),
    ] {
        assert_eq!(
            std::fs::read(dir.join(a)).unwrap(),
            std::fs::read(dir.join(b)).unwrap(),
            "{a} = {b}"
        );
    }
}

#[test]
fn multi_volume_databases_match_ncbi_and_the_fasta() {
    for (name, seqtype, fasta, version) in [
        ("prot_mv", SeqType::Protein, "prot.faa", 4),
        ("nucl_mv", SeqType::Nucleotide, "nucl.fna", 5),
    ] {
        let (db, entries) = check_against_ncbi(name, seqtype);
        check_entries_against_fasta(&entries, fasta, name);
        assert!(db.volumes().len() > 2, "{name}: several volumes");
        let info = expected_info(name);
        assert_eq!(
            volume_names(&db).join(" "),
            info["volumes"],
            "{name}: volume order"
        );
        for volume in db.volumes() {
            assert_eq!(volume.index().version, version);
        }
        // A block never crosses a volume, whatever the budget.
        let whole = db.next_block(0, u64::MAX).unwrap();
        assert_eq!(whole.end, db.volumes()[0].num_oids());
    }
}

#[test]
fn alias_files_give_ncbis_volume_order_title_and_date() {
    // No TITLE: the titles below are joined with "; ". The volumes are opened
    // in the order of their names, not in DBLIST order (prot_v4 is listed
    // first and read last).
    let (db, entries) = check_against_ncbi("alias_nested", SeqType::Protein);
    assert_eq!(
        volume_names(&db),
        [
            "prot_mv.00",
            "prot_mv.01",
            "prot_mv.02",
            "prot_mv.03",
            "prot_mv.04",
            "prot_mv.05",
            "prot_old",
            "prot_v4"
        ]
    );
    assert_eq!(entries.len(), 30);
    assert_eq!(
        db.title(),
        b"protein in volumes; LOSAT seqdb fixture, protein, version 4; LOSAT seqdb fixture, protein, version 4"
    );

    // TITLE given twice (the later line wins), tabs, CR LF, an unknown key.
    // prot_old carries an earlier date than prot_v4: the later date is reported.
    let (db, _) = check_against_ncbi("alias_titled", SeqType::Protein);
    assert_eq!(db.title(), b"the second TITLE line wins");
    assert_eq!(volume_names(&db), ["prot_old", "prot_v4"]);
    let old = SeqDb::open(&db_name("prot_old"), SeqType::Protein).unwrap();
    assert_eq!(old.date(), "Jan 5, 2020  11:07 PM");
    assert_ne!(db.date(), old.date());
    check_against_ncbi("prot_old", SeqType::Protein);
}

#[test]
fn both_ambiguity_formats_are_read() {
    let db = SeqDb::open(&db_name("nucl_v4"), SeqType::Nucleotide).unwrap();
    let block = db.read_block(db.next_block(0, u64::MAX).unwrap()).unwrap();
    assert_eq!(block.len(), 12);
    // n05: runs of 5, 1 and 1 fit the old format (four bytes per run).
    let old = block.ambiguity_bytes(4);
    assert_eq!(old[0] & 0x80, 0, "old format");
    assert_eq!(
        block.ambiguity_runs(4).unwrap(),
        [
            AmbiguityRun {
                offset: 30,
                length: 5,
                ncbi4na: 15
            },
            AmbiguityRun {
                offset: 55,
                length: 1,
                ncbi4na: 5
            },
            AmbiguityRun {
                offset: 59,
                length: 1,
                ncbi4na: 10
            },
        ]
    );
    // n06: a run of 40 takes the new format (eight bytes per run).
    let new = block.ambiguity_bytes(5);
    assert_eq!(new[0] & 0x80, 0x80, "new format");
    let runs = block.ambiguity_runs(5).unwrap();
    assert_eq!(runs.len() * 8 + 4, new.len());
    assert!(runs.contains(&AmbiguityRun {
        offset: 37,
        length: 40,
        ncbi4na: 15
    }));
    // n01: no ambiguity data.
    assert!(block.ambiguity_bytes(0).is_empty());
    assert!(block.ambiguity_runs(0).unwrap().is_empty());
    // n07, the one base "G": makeblastdb 2.17.0 writes ambiguity data that
    // hold the count zero.
    assert_eq!(block.ambiguity_bytes(6), [0, 0, 0, 0]);
    assert!(block.ambiguity_runs(6).unwrap().is_empty());
    assert_eq!(block.nucleotide_blastna(6).unwrap(), [15, 2, 15]);
}

fn open_error(name: &str, seqtype: SeqType) -> SeqDbError {
    match SeqDb::open(&db_name(name), seqtype) {
        Ok(_) => panic!("{name} opened"),
        Err(error) => error,
    }
}

#[test]
fn databases_that_need_unread_data_are_rejected() {
    for name in ["prot_ids_v4", "prot_ids_v5"] {
        let error = open_error(name, SeqType::Protein);
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("parsed Seq-ids")),
            "{name}: {error}"
        );
    }
    for (name, key) in [("alias_nseq", "NSEQ"), ("alias_gilist", "GILIST")] {
        let error = open_error(name, SeqType::Protein);
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains(key)),
            "{name}: {error}"
        );
    }
    // NCBI counts prot_v4 twice here (30 sequences reported, 20 read).
    let error = open_error("alias_twice", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("more than once")),
        "{error}"
    );
}

#[test]
fn invalid_databases_give_ncbis_errors() {
    let error = open_error("alias_mixed", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::Format { message, .. }
            if message == "DB list contains both Version 4 and Version 5 dbs"),
        "{error}"
    );
    let error = open_error("alias_loop_a", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::Format { message, .. }
            if message == "Illegal configuration: DB alias files are mutually recursive."),
        "{error}"
    );
    let error = open_error("alias_empty", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::Format { message, .. }
            if message.starts_with("No database names were found in alias file [")),
        "{error}"
    );
    let error = open_error("alias_missing", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::NotFound(message)
            if message.starts_with("Could not find volume or alias file (no_such_volume) referenced in alias file (")),
        "{error}"
    );
    let error = open_error("no_such_database", SeqType::Protein);
    assert!(
        matches!(&error, SeqDbError::NotFound(message)
            if message.starts_with("No alias or index file found for protein database [")),
        "{error}"
    );
    // A protein database asked for as nucleotide has no .nin file.
    let error = open_error("prot_v4", SeqType::Nucleotide);
    assert!(
        matches!(&error, SeqDbError::NotFound(message)
            if message.starts_with("No alias or index file found for nucleotide database [")),
        "{error}"
    );
}

#[test]
fn index_file_errors() {
    let path = fixture_dir().join("db/prot_v4.pin");
    let bytes = std::fs::read(&path).unwrap();
    let index = IndexFile::parse(&bytes, SeqType::Protein, &path).unwrap();
    assert_eq!(index.num_oids, 10);
    assert_eq!(index.volume_length, 801);
    assert_eq!(index.max_length, 300);
    assert_eq!(index.title, b"LOSAT seqdb fixture, protein, version 4");
    assert_eq!(index.sequence_range(0), (1, 63));
    assert_eq!(index.ambiguity_range(0), None);

    let message = |result: Result<IndexFile, SeqDbError>| match result {
        Err(SeqDbError::Format { message, .. }) => message,
        other => panic!("expected a format error, got {other:?}"),
    };
    assert_eq!(
        message(IndexFile::parse(&bytes, SeqType::Nucleotide, &path)),
        "Error: requested sequence type does not match DB."
    );
    let mut version3 = bytes.clone();
    version3[3] = 3;
    assert_eq!(
        message(IndexFile::parse(&version3, SeqType::Protein, &path)),
        "Error: Not a valid version 4 or 5 database."
    );
    for cut in [0, 3, 20, bytes.len() - 1] {
        assert!(
            message(IndexFile::parse(&bytes[..cut], SeqType::Protein, &path))
                .starts_with("index file ends inside a field")
        );
    }
    // Offsets that go backwards.
    let mut backwards = bytes.clone();
    let last = backwards.len() - 4;
    backwards[last..].copy_from_slice(&0u32.to_be_bytes());
    assert!(
        message(IndexFile::parse(&backwards, SeqType::Protein, &path))
            .contains("not in ascending order")
    );

    let path = fixture_dir().join("db/nucl_v5.nin");
    let bytes = std::fs::read(&path).unwrap();
    let index = IndexFile::parse(&bytes, SeqType::Nucleotide, &path).unwrap();
    assert_eq!(index.version, 5);
    assert_eq!(index.lmdb_file.as_deref(), Some(b"nucl_v5.ndb".as_slice()));
    assert_eq!(index.sequence_range(0), (1, 12));
    assert_eq!(index.ambiguity_range(0), Some((12, 12)));
}

fn temp_dir(tag: &str) -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "losat-seqdb-{tag}-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn copy_database(name: &str, exts: &[&str], dir: &Path) {
    for ext in exts {
        let file = format!("{name}.{ext}");
        std::fs::copy(fixture_dir().join("db").join(&file), dir.join(&file)).unwrap();
    }
}

#[test]
fn column_files_and_missing_files() {
    let dir = temp_dir("files");
    copy_database("prot_v4", &["pin", "psq", "phr"], &dir);
    let name = dir.join("prot_v4").to_str().unwrap().to_string();
    assert!(SeqDb::open(&name, SeqType::Protein).is_ok());

    // A column file (mask data): rejected.
    std::fs::write(dir.join("prot_v4.paa"), b"column").unwrap();
    let error = SeqDb::open(&name, SeqType::Protein).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("column data")),
        "{error}"
    );
    std::fs::remove_file(dir.join("prot_v4.paa")).unwrap();

    // No sequence file: the database opens (its headers are read) and the
    // first block read gives NCBI's message.
    std::fs::remove_file(dir.join("prot_v4.psq")).unwrap();
    {
        let db = SeqDb::open(&name, SeqType::Protein).unwrap();
        let error = db.read_sequence(0).unwrap_err();
        assert!(
            matches!(&error, SeqDbError::NotFound(message)
                if message.starts_with("Error: File (") && message.ends_with("prot_v4.psq) not found.")),
            "{error}"
        );
    }

    // A sequence file cut short: a format error, not a read past the end.
    copy_database("prot_v4", &["psq"], &dir);
    let whole = std::fs::read(dir.join("prot_v4.psq")).unwrap();
    std::fs::write(dir.join("prot_v4.psq"), &whole[..whole.len() - 5]).unwrap();
    {
        let db = SeqDb::open(&name, SeqType::Protein).unwrap();
        assert!(db.read_sequence(0).is_ok());
        assert!(matches!(
            db.read_sequence(9),
            Err(SeqDbError::Format { message, .. }) if message.contains("passes the end of the file")
        ));
    }
    // The databases are closed before their files are removed (Windows).
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn an_alias_that_names_itself_reads_the_volume_beside_it() {
    // DBLIST naming the alias file's own name means the volume of that name
    // (seqdbalias.cpp:774-803), as in NCBI's own "nt.nal".
    let dir = temp_dir("self");
    copy_database("prot_v4", &["pin", "psq", "phr"], &dir);
    std::fs::write(
        dir.join("prot_v4.pal"),
        b"TITLE through the alias\nDBLIST prot_v4\n",
    )
    .unwrap();
    {
        let db = SeqDb::open(dir.join("prot_v4").to_str().unwrap(), SeqType::Protein).unwrap();
        assert_eq!(db.title(), b"through the alias");
        assert_eq!(db.num_seqs(), 10);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn alias_lines_are_read_as_ncbi_reads_them() {
    let values = test_read_values(
        b"#comment\n  # indented comment\nTITLE  a\ttitle \t\n\n\tignored line\nDBLIST \"a b\" c\r\nKEY\nTITLE again\n",
    );
    let pairs: Vec<(&[u8], &[u8])> = values
        .iter()
        .map(|(k, v)| (k.as_slice(), v.as_slice()))
        .collect();
    assert_eq!(
        pairs,
        [
            (b"TITLE".as_slice(), b"again".as_slice()),
            (b"DBLIST".as_slice(), b"\"a b\" c".as_slice()),
            (b"KEY".as_slice(), b"".as_slice()),
        ]
    );
    assert_eq!(
        test_split_quoted(b"one \"two words\"  three\"four\" five"),
        [b"one".as_slice(), b"two words", b"three", b"four", b"five"]
    );
    assert!(test_split_quoted(b"  ").is_empty());
    assert_eq!(
        test_split_quoted(b"\"open quote"),
        [b"open quote".as_slice()]
    );
}

/// BER bytes of one `Blast-def-line-set` with one definition line.
fn defline_record(title: &[u8], seqid: &[u8]) -> Vec<u8> {
    let mut record = vec![0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a];
    if title.len() < 128 {
        record.push(title.len() as u8);
    } else {
        record.extend_from_slice(&[0x82, (title.len() >> 8) as u8, title.len() as u8]);
    }
    record.extend_from_slice(title);
    record.extend_from_slice(&[0, 0, 0xa1, 0x80, 0x30, 0x80]);
    record.extend_from_slice(seqid);
    record.extend_from_slice(&[0, 0, 0, 0, 0xa2, 0x80, 0x02, 0x01, 0x00, 0, 0, 0, 0, 0, 0]);
    record
}

fn bl_ord_id(integer: &[u8]) -> Vec<u8> {
    let mut seqid = vec![0xaa, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x09];
    seqid.extend_from_slice(b"BL_ORD_ID");
    seqid.extend_from_slice(&[0, 0, 0xa1, 0x80, 0xa0, 0x80, 0x02, integer.len() as u8]);
    seqid.extend_from_slice(integer);
    seqid.extend_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0]);
    seqid
}

#[test]
fn header_records_are_decoded() {
    let path = Path::new("test.phr");
    let decode = |record: &[u8]| decode_defline_set(record, path);

    // The first record of the fixture, byte for byte.
    let fixture = std::fs::read(fixture_dir().join("db/prot_v4.phr")).unwrap();
    let first = defline_record(b"p01 first protein [Test organism]", &bl_ord_id(&[0]));
    assert_eq!(&fixture[..first.len()], first.as_slice());
    let deflines = decode(&first).unwrap();
    assert_eq!(deflines.len(), 1);
    assert_eq!(deflines[0].title, b"p01 first protein [Test organism]");
    assert_eq!(deflines[0].volume_ordinal, 0);

    // Ordinal ids of one to four content bytes; a long title (two length bytes).
    for (integer, value) in [
        (vec![0x7f], 127u32),
        (vec![0x00, 0x80], 128),
        (vec![0x01, 0x00, 0x00], 65_536),
        (vec![0x7f, 0xff, 0xff, 0xff], i32::MAX as u32),
    ] {
        let title = vec![b't'; 300];
        let deflines = decode(&defline_record(&title, &bl_ord_id(&integer))).unwrap();
        assert_eq!(deflines[0].volume_ordinal, value);
        assert_eq!(deflines[0].title, title);
    }
    assert!(decode(&[]).unwrap().is_empty());

    // Definite lengths for the constructed values (NCBI reads both forms).
    // The db is "BL", not "BL_ORD_ID": a parsed Seq-id.
    let definite = [
        0x30, 0x1c, 0x30, 0x1a, 0xa0, 0x03, 0x1a, 0x01, b'x', 0xa1, 0x13, 0x30, 0x11, 0xaa, 0x0f,
        0x30, 0x0d, 0xa0, 0x04, 0x1a, 0x02, b'B', b'L', 0xa1, 0x05, 0xa0, 0x03, 0x02, 0x01, 0x07,
    ];
    assert!(matches!(
        decode(&definite),
        Err(SeqDbError::Unsupported { .. })
    ));

    // A local id (makeblastdb -parse_seqids): unsupported, not a format error.
    let local = [
        0xa0, 0x80, 0xa1, 0x80, 0x1a, 0x03, b'p', b'0', b'1', 0, 0, 0, 0,
    ];
    assert!(matches!(
        decode(&defline_record(b"title", &local)),
        Err(SeqDbError::Unsupported { .. })
    ));
    // Two ids in one definition line: unsupported.
    let mut two = bl_ord_id(&[1]);
    two.extend_from_slice(&bl_ord_id(&[2]));
    assert!(matches!(
        decode(&defline_record(b"title", &two)),
        Err(SeqDbError::Unsupported { .. })
    ));
    // A negative ordinal id, a cut record, bytes after the record.
    assert!(matches!(
        decode(&defline_record(b"title", &bl_ord_id(&[0xff]))),
        Err(SeqDbError::Format { .. })
    ));
    for cut in 1..first.len() {
        assert!(
            matches!(decode(&first[..cut]), Err(SeqDbError::Format { .. })),
            "record cut at {cut}"
        );
    }
    let mut trailing = first.clone();
    trailing.push(0);
    assert!(matches!(decode(&trailing), Err(SeqDbError::Format { .. })));
}

/// `expected/date_pairs.tsv`: the date NCBI's blastdbcmd reports for a database
/// of two volumes with the given date strings, or `ERROR` where it fails.
#[test]
fn the_database_date_is_the_one_ncbi_reports() {
    let path = fixture_dir().join("expected/date_pairs.tsv");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut errors = 0;
    for line in text.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 3, "{line:?}");
        let got = date::database_date(&[cols[0], cols[1]]);
        match cols[2] {
            "ERROR" => {
                assert!(got.is_err(), "{line:?}: NCBI fails, LOSAT gives {got:?}");
                errors += 1;
            }
            reported => assert_eq!(got.as_deref(), Ok(reported), "{line:?}"),
        }
    }
    assert_eq!(text.lines().count(), 29);
    assert_eq!(errors, 14);
}

/// Compares whole databases outside the repository with their FASTA files.
///
/// `LOSAT_SEQDB_FULL_CHECK` lists `prot|nucl,DATABASE,FASTA` triples separated
/// by `;`. Every sequence is read in blocks of 16 MiB and compared with the
/// FASTA record of the same position; every title is compared too.
#[test]
#[ignore = "needs LOSAT_SEQDB_FULL_CHECK (databases outside the repository)"]
fn full_databases_match_their_fasta() {
    let spec = std::env::var("LOSAT_SEQDB_FULL_CHECK").expect("LOSAT_SEQDB_FULL_CHECK");
    for triple in spec.split(';').filter(|s| !s.is_empty()) {
        let parts: Vec<&str> = triple.split(',').collect();
        assert_eq!(parts.len(), 3, "prot|nucl,DATABASE,FASTA");
        let seqtype = match parts[0] {
            "prot" => SeqType::Protein,
            "nucl" => SeqType::Nucleotide,
            other => panic!("unknown sequence type {other}"),
        };
        let db = SeqDb::open(parts[1], seqtype).unwrap_or_else(|e| panic!("{}: {e}", parts[1]));
        let records = read_fasta(Path::new(parts[2]));
        assert_eq!(
            db.num_seqs() as usize,
            records.len(),
            "{}: sequences",
            parts[1]
        );
        let mut next = 0u32;
        let mut blocks = 0u32;
        let mut letters = 0u64;
        let mut longest = 0usize;
        let mut ambiguous_sequences = 0u32;
        while let Some(range) = db.next_block(next, 16 << 20) {
            let block = db.read_block(range).expect("read block");
            for i in 0..block.len() {
                let ordinal = next + i as u32;
                let record = &records[ordinal as usize];
                let what = format!("{} ordinal {ordinal}", parts[1]);
                check_sequence(&block, i, &record.seq, &what);
                letters += record.seq.len() as u64;
                longest = longest.max(record.seq.len());
                if seqtype == SeqType::Nucleotide && !block.ambiguity_bytes(i).is_empty() {
                    ambiguous_sequences += 1;
                }
                let deflines = db.deflines(ordinal).expect("deflines");
                assert_eq!(deflines.len(), 1, "{what}: one definition line");
                assert_eq!(deflines[0].title, record.title, "{what}: title");
                assert_eq!(deflines[0].bl_ord_id, ordinal, "{what}: BL_ORD_ID");
            }
            next = range.end_ordinal();
            blocks += 1;
        }
        assert_eq!(next, db.num_seqs());
        assert_eq!(db.total_length(), letters, "{}: letters", parts[1]);
        assert_eq!(db.max_length() as usize, longest, "{}: longest", parts[1]);
        println!(
            "{}: {} sequences, {} letters, longest {}, {} volumes, {} blocks, {} sequences with ambiguity data, title {:?}, date {:?}: all equal to {}",
            parts[1],
            db.num_seqs(),
            letters,
            longest,
            db.volumes().len(),
            blocks,
            ambiguous_sequences,
            String::from_utf8_lossy(db.title()),
            db.date(),
            parts[2]
        );
    }
}
