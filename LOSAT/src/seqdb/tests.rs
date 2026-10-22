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
use super::header::{decode_defline_set, decode_defline_taxids};
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
            assert_eq!(
                block.nucleotide_len(i).unwrap(),
                upper.len(),
                "{what}: length"
            );
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
                volume.block_bytes(range.first, range.end).unwrap() <= max_bytes,
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

/// The length and the letters of ranges of a nucleotide sequence, read without the rest of
/// the sequence, are those NCBI reads back (the whole sequence, ambiguity codes included).
fn check_nucleotide_ranges(db: &SeqDb, entries: &[Entry], name: &str) {
    for entry in entries {
        let what = format!("{name} ordinal {}", entry.ordinal);
        assert_eq!(
            db.nucleotide_length(entry.ordinal).unwrap(),
            entry.seq.len(),
            "{what}: length"
        );
        let n = entry.seq.len();
        let mut ranges = vec![(0, n), (0, 0), (n, n)];
        for (from, to) in [
            (0, 1),
            (1, 5),
            (3, 4),
            (2, 9),
            (n / 2, n),
            (n.saturating_sub(7), n),
        ] {
            if from <= to && to <= n {
                ranges.push((from, to));
            }
        }
        for step in [1usize, 3, 5, 61] {
            for from in (0..n).step_by(step.max(1)) {
                ranges.push((from, (from + step + 2).min(n)));
            }
        }
        for (from, to) in ranges {
            assert_eq!(
                db.nucleotide_iupac_range(entry.ordinal, from, to).unwrap(),
                entry.seq[from..to],
                "{what}: letters {from}..{to}"
            );
        }
        assert!(db.nucleotide_iupac_range(entry.ordinal, 0, n + 1).is_err());
        assert!(db.nucleotide_iupac_range(entry.ordinal, 2, 1).is_err());
    }
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
    if seqtype == SeqType::Nucleotide {
        check_nucleotide_ranges(&db, &entries, name);
    }
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

    // Volumes with different proteins. NCBI compares volume names without a
    // final ".pdb" (and ".pal", ".pin", ...): "ord.pdb" comes before
    // "ord-extra", although the DBLIST and plain byte order say otherwise.
    // The title still follows the DBLIST. One name carries a directory.
    let (db, entries) = check_against_ncbi("alias_order", SeqType::Protein);
    assert_eq!(volume_names(&db), ["ord.pdb", "ord-extra"]);
    assert_eq!(db.title(), b"next three proteins; first four proteins");
    let titles: Vec<&[u8]> = entries.iter().map(|e| &e.title[..3]).collect();
    assert_eq!(
        titles,
        [b"p01", b"p02", b"p03", b"p04", b"p05", b"p06", b"p07"]
    );

    // A volume without sequences, first in the order.
    let (db, _) = check_against_ncbi("alias_emptyvol", SeqType::Protein);
    assert_eq!(volume_names(&db), ["empty", "prot_v4"]);
    assert_eq!(db.volumes()[0].num_oids(), 0);
    assert_eq!(db.next_block(0, 0).unwrap().volume, 1);
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
    for (name, key) in [("alias_gilist", "GILIST")] {
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
fn alias_nseq_and_length_that_equal_the_sums_are_read() {
    // blastdb_aliastool writes NSEQ 20 and LENGTH 1602, the sums over prot_v4
    // and prot_old. NCBI's walk takes the values instead of the sums of
    // that node, so the database is the two volumes.
    let (db, entries) = check_against_ncbi("alias_nseq", SeqType::Protein);
    assert_eq!(db.num_seqs(), 20);
    assert_eq!(db.total_length(), 1602);
    assert_eq!(db.title(), b"alias with NSEQ and LENGTH");
    assert_eq!(volume_names(&db), ["prot_old", "prot_v4"]);
    let mut volume_entries = expected_entries("prot_old");
    volume_entries.extend(expected_entries("prot_v4"));
    assert_eq!(entries.len(), volume_entries.len());
    for (entry, volume_entry) in entries.iter().zip(&volume_entries) {
        // The id of a record is its ordinal id, which differs.
        assert_eq!(entry.length, volume_entry.length);
        assert_eq!(entry.seq, volume_entry.seq);
        assert_eq!(entry.title, volume_entry.title);
    }

    // A nested alias file whose inner and outer nodes both have right values
    // (20 and 1602 inside; 30 and 2403 over the inner node and prot_mv).
    let (db, entries) = check_against_ncbi("alias_sums_nested", SeqType::Protein);
    assert_eq!(db.num_seqs(), 30);
    assert_eq!(db.total_length(), 2403);
    assert_eq!(entries.len(), 30);
    assert_eq!(db.volumes().len(), 8);
}

#[test]
fn alias_nseq_and_length_that_differ_from_the_sums_are_rejected() {
    for (name, key, shown) in [
        ("alias_nseq_high", "NSEQ", "(21)"),
        ("alias_length_low", "LENGTH", "(1601)"),
        ("alias_nseq_plus", "NSEQ", "(+20)"),
        ("alias_nseq_exp", "NSEQ", "(2e1)"),
        ("alias_length_space", "LENGTH", "(1 602)"),
        // The inner node is right and the outer one has LENGTH 2402 (sum 2403).
        ("alias_sums_nested_bad", "LENGTH", "(2402)"),
        // The inner alias file is wrong, although the outer node's values would
        // make NCBI skip it.
        ("alias_sums_inner_bad", "NSEQ", "(21)"),
    ] {
        let error = open_error(name, SeqType::Protein);
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. }
                if message.contains(&format!("alias file key {key} {shown}"))),
            "{name}: {error}"
        );
        let SeqDbError::Unsupported { path, .. } = &error else {
            unreachable!()
        };
        assert!(
            path.extension().is_some_and(|ext| ext == "pal"),
            "{name}: the error names an alias file: {}",
            path.display()
        );
    }
    // The wrong value is named together with the sum.
    let error = open_error("alias_nseq_high", SeqType::Protein);
    assert!(
        error.to_string().contains(
            "alias file key NSEQ (21) differs from the number of sequences of its databases (20)"
        ),
        "{error}"
    );
    let error = open_error("alias_length_low", SeqType::Protein);
    assert!(
        error.to_string().contains(
            "alias file key LENGTH (1601) differs from the total length of its databases (1602)"
        ),
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
    assert_eq!(index.sequence_range(0).unwrap(), (1, 63));
    assert_eq!(index.ambiguity_range(0).unwrap(), None);

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
    assert_eq!(index.sequence_range(0).unwrap(), (1, 12));
    assert_eq!(index.ambiguity_range(0).unwrap(), Some((12, 12)));
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

/// The taxid of a definition line is skipped by `decode_defline_set` and read apart by
/// `decode_defline_taxids` (the XML2 and JSON formats write it).
#[test]
fn the_taxid_of_a_definition_line_is_read_apart() {
    let path = Path::new("test.phr");
    let plain = defline_record(b"t1 title", &bl_ord_id(&[0x05]));
    // `taxid [2] INTEGER` 0, then 9606 (0x2586), then a negative value read as it is
    let mut tagged = Vec::new();
    for (taxid, bytes) in [(0, vec![0x00]), (9606, vec![0x25, 0x86]), (-1, vec![0xff])] {
        let mut record = plain.clone();
        let at = record
            .windows(5)
            .position(|window| window == [0xa2, 0x80, 0x02, 0x01, 0x00])
            .unwrap();
        let mut member = vec![0xa2, 0x80, 0x02, bytes.len() as u8];
        member.extend_from_slice(&bytes);
        record.splice(at..at + 5, member);
        tagged.push((taxid, record));
    }
    for (taxid, record) in &tagged {
        assert_eq!(decode_defline_taxids(record, path).unwrap(), [*taxid]);
        // the title and ordinal id are read the same, and the default decode leaves the taxid 0
        let deflines = decode_defline_set(record, path).unwrap();
        assert_eq!(deflines[0].title, b"t1 title");
        assert_eq!(deflines[0].volume_ordinal, 5);
        assert_eq!(deflines[0].taxid, 0);
    }
    // a definition line without the member
    let mut record = defline_record(b"t1", &bl_ord_id(&[0x05]));
    let at = record
        .windows(5)
        .position(|window| window == [0xa2, 0x80, 0x02, 0x01, 0x00])
        .unwrap();
    record.drain(at..at + 7);
    assert_eq!(decode_defline_taxids(&record, path).unwrap(), [0]);
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

#[test]
fn every_filtering_or_statistics_alias_key_is_rejected() {
    let dir = temp_dir("keys");
    copy_database("prot_v4", &["pin", "psq", "phr"], &dir);
    let name = dir.join("keyed").to_str().unwrap().to_string();
    for key in [
        "FIRST_OID",
        "LAST_OID",
        "GILIST",
        "TILIST",
        "SEQIDLIST",
        "TAXIDLIST",
        "OIDLIST",
        "MEMB_BIT",
        "MASKLIST",
        "VDBLIST",
        "NSEQ",
        "LENGTH",
        "STATS_NSEQ",
        "STATS_TOTLEN",
        "MIN_SEQ_LENGTH",
        "OID_MASK_TYPE",
    ] {
        std::fs::write(
            dir.join("keyed.pal"),
            format!("TITLE t\nDBLIST prot_v4\n{key} 1\n"),
        )
        .unwrap();
        let error = SeqDb::open(&name, SeqType::Protein).unwrap_err();
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains(key)),
            "{key}: {error}"
        );
    }
    // Keys NCBI does not apply are ignored: an unknown key, and MAX_SEQ_LENGTH.
    std::fs::write(
        dir.join("keyed.pal"),
        b"TITLE t\nDBLIST prot_v4\nMAX_SEQ_LENGTH 5\nSOME_KEY x\n",
    )
    .unwrap();
    assert_eq!(
        SeqDb::open(&name, SeqType::Protein).unwrap().max_length(),
        300
    );

    // NSEQ and LENGTH that equal the sums are read; a repeated key takes the
    // later value (seqdbalias.cpp:446-449), here the right one.
    std::fs::write(
        dir.join("keyed.pal"),
        b"TITLE t\nDBLIST prot_v4\nNSEQ 99\nLENGTH 801\nNSEQ 10\n",
    )
    .unwrap();
    let db = SeqDb::open(&name, SeqType::Protein).unwrap();
    assert_eq!((db.num_seqs(), db.total_length()), (10, 801));
    std::fs::write(
        dir.join("keyed.pal"),
        b"TITLE t\nDBLIST prot_v4\nNSEQ 10\nNSEQ 99\n",
    )
    .unwrap();
    let error = SeqDb::open(&name, SeqType::Protein).unwrap_err();
    assert!(error.to_string().contains("NSEQ (99)"), "{error}");

    // A group alias file in the directory: rejected for an alias file and for
    // a plain volume alike.
    std::fs::write(dir.join("index.alx"), b"ALIAS_FILE keyed.pal\n").unwrap();
    for database in ["keyed", "prot_v4"] {
        let error =
            SeqDb::open(dir.join(database).to_str().unwrap(), SeqType::Protein).unwrap_err();
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("group alias file")),
            "{database}: {error}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_name_is_a_list_as_for_ncbi() {
    let prot_v4 = db_name("prot_v4");
    // Double quotes keep one name together and are removed.
    let quoted = format!("\"{prot_v4}\"");
    assert_eq!(
        SeqDb::open(&quoted, SeqType::Protein).unwrap().num_seqs(),
        10
    );
    // Two names: several databases in one search are not supported.
    let two = format!("{prot_v4} {}", db_name("prot_old"));
    let error = SeqDb::open(&two, SeqType::Protein).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("several databases")),
        "{error}"
    );
    let error = SeqDb::open("  ", SeqType::Protein).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::NotFound(message)
            if message == "No database names were passed to CSeqDB::CSeqDB()."),
        "{error}"
    );
}

/// One database name is a path taken as it is: a space or a double quote in
/// it does not split it (the entry point of `LOSAT makedb`, whose output name
/// is one name).
#[test]
fn one_name_is_not_a_list() {
    let mut names = vec!["my data"];
    if cfg!(unix) {
        names.push("q\"d");
    }
    for name in names {
        let dir = temp_dir("one-name").join(name);
        std::fs::create_dir_all(&dir).unwrap();
        copy_database("prot_v4", &["pin", "phr", "psq"], &dir);
        let path = dir.join("prot_v4").to_str().unwrap().to_string();
        assert!(SeqDb::open(&path, SeqType::Protein).is_err(), "{name}");
        let db = SeqDb::open_one_name(&path, SeqType::Protein).unwrap();
        assert_eq!(db.num_seqs(), 10, "{name}");
        let (volumes, alias_files) =
            super::alias::find_volume_paths(&path, SeqType::Protein).unwrap();
        assert_eq!((volumes, alias_files), (vec![path], Vec::new()), "{name}");
        std::fs::remove_dir_all(dir.parent().unwrap()).unwrap();
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:831-841
// ```c++
//         // Always check local copy first unless full path is given
//         if (!m_SkipLocal[i]) {
//             // GCC 3.0.4 requires the extra parentheses below.
//             CSeqDB_BasePath local_base((CSeqDB_DirName(CDir::GetCwd())),
//                                    CSeqDB_BaseName(m_DBList[i].FindBaseName()));
//             CSeqDB_Path new_local_vol_path(local_base, prot_nucl, 'i', 'n' );
//             if (m_Atlas.DoesFileExist(new_local_vol_path)) {
//                 bp = CSeqDB_BasePath(new_local_vol_path.FindBasePath());
//                 found = true;
//             }
//         }
// ```
/// A volume that an alias file names is taken from the current directory
/// when a volume of that name is there (NCBI's rule, for reading). The
/// removal of a database before a build must not follow it: it lists the
/// volumes beside the alias file. The current directory of the test process
/// is not changed; the resolver takes the directory as an argument.
#[test]
fn volumes_to_remove_are_those_beside_the_alias_file() {
    let root = temp_dir("local-copy");
    let (sub, here) = (root.join("sub"), root.join("here"));
    for dir in [&sub, &here] {
        std::fs::create_dir(dir).unwrap();
    }
    std::fs::write(sub.join("db.pal"), b"TITLE t\nDBLIST db.00 db.01\n").unwrap();
    // The resolver looks for the index files only.
    for file in [
        sub.join("db.00.pin"),
        sub.join("db.01.pin"),
        here.join("db.00.pin"),
    ] {
        std::fs::write(file, b"").unwrap();
    }
    let path = |dir: &Path, name: &str| dir.join(name).to_str().unwrap().to_string();
    let name = path(&sub, "db");

    let read = super::alias::resolve_one_name(
        &name,
        SeqType::Protein,
        super::alias::LocalCopy::Dir(here.to_str().unwrap()),
    )
    .unwrap();
    assert_eq!(read.volumes, [path(&here, "db.00"), path(&sub, "db.01")]);

    let (volumes, alias_files) = super::alias::find_volume_paths(&name, SeqType::Protein).unwrap();
    assert_eq!(volumes, [path(&sub, "db.00"), path(&sub, "db.01")]);
    assert_eq!(alias_files, [path(&sub, "db.pal")]);

    // A volume that is missing beside the alias file is not looked for
    // elsewhere: the database cannot be listed, and nothing is removed.
    std::fs::remove_file(sub.join("db.00.pin")).unwrap();
    let error = super::alias::find_volume_paths(&name, SeqType::Protein).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::NotFound(message) if message
            == &format!("Could not find volume or alias file (db.00) referenced in alias file ({name}).")),
        "{error}"
    );
    std::fs::remove_dir_all(&root).unwrap();
}

#[test]
fn a_version_5_database_needs_its_lmdb_file() {
    let dir = temp_dir("lmdb");
    copy_database("nucl_v5", &["nin", "nsq", "nhr"], &dir);
    let name = dir.join("nucl_v5").to_str().unwrap().to_string();
    let error = SeqDb::open(&name, SeqType::Nucleotide).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::NotFound(message)
            if message.starts_with("File ") && message.contains("nucl_v5.ndb not found. If you renamed")),
        "{error}"
    );
    copy_database("nucl_v5", &["ndb"], &dir);
    assert_eq!(
        SeqDb::open(&name, SeqType::Nucleotide).unwrap().num_seqs(),
        12
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn ambiguity_data_that_do_not_fit_the_sequence() {
    // n05 of the fixture (74 bases) has the ambiguity words
    // 00000003 f400001e 50000037 a000003b at offset 64 of the sequence file:
    // three runs in the old format (N x 5 at 30, R at 55, Y at 59).
    let dir = temp_dir("ambiguity");
    copy_database("nucl_v4", &["nin", "nsq", "nhr"], &dir);
    let name = dir.join("nucl_v4").to_str().unwrap().to_string();
    let path = dir.join("nucl_v4.nsq");
    let original = std::fs::read(&path).unwrap();
    assert_eq!(
        original[64..80],
        [0, 0, 0, 3, 0xf4, 0, 0, 0x1e, 0x50, 0, 0, 0x37, 0xa0, 0, 0, 0x3b]
    );
    let read = |bytes: &[u8]| {
        std::fs::write(&path, bytes).unwrap();
        let db = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
        db.read_sequence(4).unwrap()
    };
    let intact = read(&original).nucleotide_blastna(0).unwrap();
    // BLASTNA: N = 14, R = 4, Y = 5; the stored 2-bit bases are 0..=3.
    assert_eq!(intact[1 + 30..1 + 35], [14; 5]);
    assert_eq!((intact[1 + 55], intact[1 + 59]), (4, 5));

    // A run that starts past the end of the sequence ends the loop: the run
    // after it is not applied either (seqdbvol.cpp:837-838, `break`).
    // Runs now: N x 5 at 30, R at 74 (the sequence has 74 bases), Y at 59.
    let mut past_end = original.clone();
    past_end[72..76].copy_from_slice(&[0x50, 0, 0, 74]);
    let clipped = read(&past_end).nucleotide_blastna(0).unwrap();
    assert_eq!(clipped[1 + 30..1 + 35], [14; 5]);
    assert!(clipped[1 + 55] < 4, "R is not restored");
    assert!(
        clipped[1 + 59] < 4,
        "the run after the break is not applied"
    );

    // A run that crosses the end is applied up to the end.
    let mut crossing = original.clone();
    crossing[76..80].copy_from_slice(&[0xf4, 0, 0, 72]);
    let crossed = read(&crossing).nucleotide_blastna(0).unwrap();
    assert_eq!(crossed.len(), 74 + 2);
    assert_eq!(crossed[1 + 72..1 + 74], [14, 14]);
    assert_eq!(crossed[75], BLASTNA_SENTINEL);

    // A count that passes the end of the ambiguity data.
    let mut long_count = original.clone();
    long_count[64..68].copy_from_slice(&[0, 0, 0, 9]);
    let block = read(&long_count);
    assert!(matches!(
        block.nucleotide_blastna(0),
        Err(SeqDbError::Format { message, .. }) if message.contains("shorter than their count")
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn header_members_lengths_and_bytes() {
    let path = Path::new("test.phr");
    let decode = |record: &[u8]| decode_defline_set(record, path);

    // Definite lengths throughout, and no title (it is optional).
    let seqid = [
        0xaa, 0x18, 0x30, 0x16, 0xa0, 0x0b, 0x1a, 0x09, b'B', b'L', b'_', b'O', b'R', b'D', b'_',
        b'I', b'D', 0xa1, 0x07, 0xa0, 0x05, 0x02, 0x03, 0x01, 0x00, 0x01,
    ];
    let mut record = vec![0x30, 0x20, 0x30, 0x1e, 0xa1, 0x1c, 0x30, 0x1a];
    record.extend_from_slice(&seqid);
    let deflines = decode(&record).unwrap();
    assert_eq!(deflines.len(), 1);
    assert!(deflines[0].title.is_empty());
    assert_eq!(deflines[0].volume_ordinal, 65_537);

    // taxid, memberships, links and other-info are skipped.
    let mut members = vec![0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x01, b't', 0, 0];
    members.extend_from_slice(&[0xa1, 0x80, 0x30, 0x80]);
    members.extend_from_slice(&bl_ord_id(&[3]));
    members.extend_from_slice(&[0, 0, 0, 0]);
    members.extend_from_slice(&[0xa2, 0x80, 0x02, 0x02, 0x25, 0x8a, 0, 0]);
    for member in [0xa3, 0xa4, 0xa5] {
        members.extend_from_slice(&[member, 0x80, 0x30, 0x80, 0x02, 0x01, 0x01, 0x02, 0x01, 0x02]);
        members.extend_from_slice(&[0, 0, 0, 0]);
    }
    members.extend_from_slice(&[0, 0, 0, 0]);
    let deflines = decode(&members).unwrap();
    assert_eq!(deflines[0].title, b"t");
    assert_eq!(deflines[0].volume_ordinal, 3);

    // Members out of order, and values nested too deeply in a skipped member.
    let mut out_of_order = vec![0x30, 0x80, 0x30, 0x80, 0xa2, 0x80, 0x02, 0x01, 0x00, 0, 0];
    out_of_order.extend_from_slice(&[0xa0, 0x80, 0x1a, 0x01, b't', 0, 0]);
    assert!(matches!(
        decode(&out_of_order),
        Err(SeqDbError::Format { .. })
    ));
    let mut deep = vec![0x30, 0x80, 0x30, 0x80, 0xa2, 0x80];
    deep.extend(std::iter::repeat_n([0x30u8, 0x80], 100_000).flatten());
    assert!(matches!(
        decode(&deep),
        Err(SeqDbError::Format { message, .. }) if message.contains("nested too deeply")
    ));

    // Bytes outside ' ' ..= '~' in a title are read as '#', as NCBI reads them.
    let deflines = decode(&defline_record(
        b"p10\tlast \xc2\xb5 protein\x7f",
        &bl_ord_id(&[0]),
    ))
    .unwrap();
    assert_eq!(deflines[0].title, b"p10#last ## protein#");
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

// ---------------------------------------------------------------------------
// The database search path: what NCBI BLAST+ 2.17.0 does with the name given
// to `-db`. `expected/search_path/*.tsv` were written from runs of NCBI's
// blastp and blastdbcmd (script `db3_path_cases.py` of the campaign's task
// folder; see `tests/seqdb/README.md`). The tests build the same directory
// tree and never change the current directory or the environment of the test
// process: the search path is a value. They run on Unix: the directory
// names, the symbolic links and the file modes of the tables are Unix's.

/// A field of the tables: backslash escapes, `\s` for a space at either end
/// and `\e` for an empty field.
#[cfg(unix)]
fn unescape(field: &str) -> String {
    let mut out = String::with_capacity(field.len());
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some('s') => out.push(' '),
            Some('e') => {}
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

/// The rows of `expected/search_path/NAME` (`header`: the first line names the
/// columns and is left out).
#[cfg(unix)]
fn search_path_table(name: &str, header: bool) -> Vec<Vec<String>> {
    let path = fixture_dir().join("expected/search_path").join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
    text.split('\n')
        .skip(usize::from(header))
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(unescape).collect())
        .collect()
}

/// Builds the directory tree of `expected/search_path/layout.tsv` under a new
/// temporary directory and returns its path without symbolic links, as the
/// operating system reports a current directory.
#[cfg(unix)]
fn search_path_sandbox(tag: &str, layout: &str) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let root = temp_dir(tag).canonicalize().unwrap();
    let root_str = root.to_str().expect("temporary directory is UTF-8");
    let source = fixture_dir().join("db");
    for op in search_path_table(layout, false) {
        let op: Vec<&str> = op.iter().map(String::as_str).collect();
        match op.as_slice() {
            ["dir", path] => std::fs::create_dir_all(root.join(path)).unwrap(),
            ["copy", from, to] => {
                std::fs::copy(source.join(from), root.join(to)).unwrap();
            }
            ["text", to, content] => {
                std::fs::write(root.join(to), content.replace("{R}", root_str)).unwrap();
            }
            ["symlink", target, to] => {
                std::os::unix::fs::symlink(target.replace("{R}", root_str), root.join(to)).unwrap();
            }
            ["mode", path, bits] => {
                let mode = u32::from_str_radix(bits, 8).unwrap();
                std::fs::set_permissions(root.join(path), std::fs::Permissions::from_mode(mode))
                    .unwrap();
            }
            ["patch", from, to, offset, byte] => {
                let mut bytes = std::fs::read(source.join(from)).unwrap();
                bytes[offset.parse::<usize>().unwrap()] = u8::from_str_radix(byte, 16).unwrap();
                std::fs::write(root.join(to), bytes).unwrap();
            }
            ["trunc", from, to, length] => {
                let bytes = std::fs::read(source.join(from)).unwrap();
                let length: i64 = length.parse().unwrap();
                let keep = if length >= 0 {
                    length as usize
                } else {
                    bytes.len() - length.unsigned_abs() as usize
                };
                std::fs::write(root.join(to), &bytes[..keep]).unwrap();
            }
            other => panic!("unknown layout line {other:?}"),
        }
    }
    root
}

/// The search path of a table row: `<unset>` and `<none>` stand for a
/// variable that is not set and a registry without the entry.
#[cfg(unix)]
fn row_search_path(cwd: &Path, blastdb: &str, registry: &str) -> SearchPath {
    SearchPath::new(
        cwd,
        (blastdb != "<unset>").then(|| std::ffi::OsStr::new(blastdb)),
        (registry != "<none>").then_some(registry),
    )
    .unwrap()
}

/// Compares what LOSAT does with one row with what NCBI did.
/// `row`: NCBI's outcome, title, sequences, volumes, error, exit status.
#[cfg(unix)]
fn check_search_path_row(
    id: &str,
    name: &str,
    search: &SearchPath,
    row: &[String],
    losat: &str,
    place: &dyn Fn(&str) -> String,
) {
    let (ncbi, title, sequences, volumes, error, status) =
        (&row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
    // The tables hold the value of blastp's `-db`, which the program trims.
    let name = db_argument_name(name);
    let opened = SeqDb::open_in(name, SeqType::Protein, search);

    // A volume without its sequence file: NCBI's search fails when it reads
    // the first sequence, with the search stage's own text and exit status.
    if let Some(text) = losat.strip_prefix("search:") {
        let db = opened.unwrap_or_else(|e| panic!("{id}: LOSAT does not open it: {e:?}"));
        let read = db.read_sequence(0).expect_err("the first read fails");
        assert_eq!(read.ncbi_search_message(), Some(text), "{id}: {read:?}");
        assert_eq!(*error, format!("BLAST engine error: {text}"), "{id}");
        assert_eq!(status, "3", "{id}: exit status of an engine error");
        return;
    }
    // A volume without its header file: NCBI opens the file at the first read of
    // a header, which a search makes after its search (`SeqDb::open_in`); the
    // error is the reader's.
    if losat == "header" {
        let db = opened.unwrap_or_else(|e| panic!("{id}: LOSAT does not open it: {e:?}"));
        db.read_sequence(0)
            .unwrap_or_else(|e| panic!("{id}: the sequence file: {e:?}"));
        let read = db.deflines(0).expect_err("the first header read fails");
        assert_eq!(read.ncbi_message(), Some(place(error)), "{id}: {read:?}");
        assert_eq!(status, "2", "{id}: exit status of a database error");
        return;
    }
    // Databases LOSAT refuses on purpose, whatever NCBI does with them.
    if let Some((kind, text)) = losat.split_once(':') {
        let error = match opened {
            Ok(_) => panic!("{id}: opened, expected {losat}"),
            Err(error) => error,
        };
        let matches = match (kind, &error) {
            ("unsupported", SeqDbError::Unsupported { message, .. })
            | ("format", SeqDbError::Format { message, .. }) => message.contains(text),
            _ => false,
        };
        assert!(matches, "{id}: expected {losat}, got {error:?}");
        assert_eq!(error.ncbi_message(), None, "{id}: a refusal of LOSAT");
        return;
    }

    match ncbi.as_str() {
        "ok" => {
            let db = opened.unwrap_or_else(|e| panic!("{id}: NCBI opens it, LOSAT: {e:?}"));
            assert_eq!(String::from_utf8_lossy(db.title()), *title, "{id}: title");
            assert_eq!(db.num_seqs().to_string(), *sequences, "{id}: sequences");
            let listed = super::alias::test_info_volumes(name, SeqType::Protein, search).unwrap();
            let expected: Vec<String> = volumes.split(" | ").map(place).collect();
            assert_eq!(listed, expected, "{id}: volumes");
        }
        "error" => {
            assert_eq!(status, "2", "{id}: exit status of a database error");
            let error_text = match opened {
                Ok(_) => panic!("{id}: LOSAT opens it, NCBI: {error}"),
                Err(e) => e
                    .ncbi_message()
                    .unwrap_or_else(|| panic!("{id}: not an NCBI error: {e:?}")),
            };
            assert_eq!(error_text, place(error), "{id}: error text");
        }
        other => panic!("{id}: outcome {other} needs a `losat` column"),
    }
}

/// `expected/search_path/cases.tsv`: one row per `-db` name with a current
/// directory, `BLASTDB` and the registry value.
#[cfg(unix)]
#[test]
fn a_name_is_found_where_ncbi_finds_it() {
    let root = search_path_sandbox("search-path", "layout.tsv");
    let root_str = root.to_str().unwrap().to_string();
    let place = |value: &str| value.replace("{R}", &root_str);
    let rows = search_path_table("cases.tsv", true);
    for row in &rows {
        assert_eq!(row.len(), 13, "{row:?}");
        let (id, cwd, blastdb, registry, name) = (&row[0], &row[1], &row[2], &row[3], &row[4]);
        let search = row_search_path(&root.join(cwd), &place(blastdb), &place(registry));
        check_search_path_row(id, &place(name), &search, &row[5..11], &row[12], &place);
    }
    assert!(rows.len() > 100, "{} cases", rows.len());
    std::fs::remove_dir_all(&root).unwrap();
}

/// `expected/search_path/relative.tsv`: relative directories in `BLASTDB` and
/// in the registry value, run from the crate root, where Cargo runs the tests
/// of the crate. The fixtures of `tests/seqdb/db` are read in place.
#[cfg(unix)]
#[test]
fn relative_directories_are_used_as_they_are_written() {
    let cwd = std::env::current_dir().unwrap();
    let crate_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    if cwd.canonicalize().ok() != crate_root.canonicalize().ok() {
        eprintln!("skipped: the tests do not run in the crate root");
        return;
    }
    let cwd_str = cwd.to_str().unwrap().to_string();
    let place = |value: &str| value.replace("{P}", &cwd_str);
    let rows = search_path_table("relative.tsv", true);
    for row in &rows {
        assert_eq!(row.len(), 11, "{row:?}");
        let (id, blastdb, registry, name) = (&row[0], &row[1], &row[2], &row[3]);
        let search = row_search_path(&cwd, blastdb, registry);
        check_search_path_row(id, name, &search, &row[4..10], "", &place);
    }
    assert!(rows.len() >= 15, "{} cases", rows.len());
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5932-5946 and
// c++/src/app/blast/blast_app_util.hpp:185-203 (quoted at `DbFile::open` and
// at `SeqDbError::ncbi_message`).
/// `expected/search_path/open_cases.tsv`: files without read permission and a
/// directory that cannot be searched. LOSAT opens a volume's files at other
/// moments than NCBI (NCBI maps the sequence and header files when the
/// search reads them), so the text compared is that of the first error LOSAT
/// meets when it opens the database and reads its first sequence and title.
/// The rows with a limit of open files are not run here (the limit belongs
/// to the process); `ncbis_errors_are_told_from_losats` covers their text.
#[cfg(unix)]
#[test]
fn files_that_cannot_be_opened_give_ncbis_errors() {
    use std::os::unix::fs::PermissionsExt;
    let root = search_path_sandbox("search-path-open", "open_layout.tsv");
    let root_str = root.to_str().unwrap().to_string();
    let place = |value: &str| value.replace("{R}", &root_str);
    // The administrator reads files of any mode: nothing to compare then.
    let readable_anyway = std::fs::File::open(root.join("d/nopin.pin")).is_ok();
    let rows = search_path_table("open_cases.tsv", true);
    let mut compared = 0;
    for row in &rows {
        assert_eq!(row.len(), 12, "{row:?}");
        let (id, cwd, blastdb, name, limit) = (&row[0], &row[1], &row[2], &row[3], &row[4]);
        let (blastp, blastp_error, losat) = (&row[5], &row[6], &row[11]);
        if !limit.is_empty() || readable_anyway {
            continue;
        }
        let search = row_search_path(&root.join(cwd), &place(blastdb), "<none>");
        let first_error = SeqDb::open_in(&place(name), SeqType::Protein, &search).and_then(|db| {
            db.read_sequence(0)?;
            db.deflines(0).map(|_| db.num_seqs())
        });
        if let Some(text) = losat.strip_prefix("search:") {
            // The sequence file: NCBI's search stage reports the failure in
            // its own words and with exit status 3.
            let error = first_error.expect_err("the first read fails");
            assert_eq!(error.ncbi_search_message(), Some(text), "{id}: {error:?}");
            assert_eq!(*blastp_error, format!("BLAST engine error: {text}"), "{id}");
            assert_eq!(row[7], "3", "{id}: exit status of an engine error");
        } else if let Some((_, text)) = losat.split_once(':') {
            assert!(
                matches!(&first_error, Err(SeqDbError::Unsupported { message, .. }) if message.contains(text)),
                "{id}: expected {losat}, got {first_error:?}"
            );
        } else if blastp == "ok" {
            assert!(first_error.is_ok(), "{id}: {first_error:?}");
        } else {
            assert_eq!(
                blastp, "error",
                "{id}: an outcome that needs a `losat` column"
            );
            let message = first_error
                .err()
                .and_then(|error| error.ncbi_message())
                .unwrap_or_else(|| panic!("{id}: no NCBI error"));
            assert_eq!(message, place(blastp_error), "{id}");
            assert_eq!(row[7], "2", "{id}: exit status of a database error");
        }
        compared += 1;
    }
    assert!(readable_anyway || compared >= 8, "{compared} cases");
    // Directories are opened again before the tree is removed.
    let closed = root.join("d/closed");
    std::fs::set_permissions(&closed, std::fs::Permissions::from_mode(0o755)).unwrap();
    std::fs::remove_dir_all(&root).unwrap();
}

/// The string of the search path, and the errors of a search path that
/// cannot be used.
#[test]
fn the_search_path_string_and_its_errors() {
    let cwd = std::env::current_dir().unwrap();
    let cwd_str = cwd.to_str().unwrap();
    let search = SearchPath::new(&cwd, None, None).unwrap();
    // The directory is the one the operating system reports; its spelling
    // is compared on Unix only (the Windows rules of NCBI's path
    // normalisation are ported from the source and have no oracle run).
    if cfg!(unix) {
        assert_eq!(search.as_str(), format!("{cwd_str}::"));
    }
    // The path of the process starts with the same directory, whatever
    // `BLASTDB` holds.
    assert!(SearchPath::from_environment(None)
        .unwrap()
        .as_str()
        .starts_with(cwd_str));
    // A registry entry adds its value and a separator, also when it is empty.
    let with_entry = SearchPath::new(&cwd, None, Some("")).unwrap();
    let separator = if cfg!(windows) { ';' } else { ':' };
    assert_eq!(
        with_entry.as_str(),
        format!("{}{separator}", search.as_str())
    );

    // `CSeqDB::CSeqDB` (seqdb.cpp:155-165) and `x_ExpandAliases` on an empty
    // list (seqdbalias.cpp:740-754).
    let message = |name: &str| {
        SeqDb::open_in(name, SeqType::Protein, &search)
            .unwrap_err()
            .ncbi_message()
    };
    assert_eq!(message("").as_deref(), Some("Database name is required."));
    assert_eq!(
        message("  ").as_deref(),
        Some("No database names were passed to CSeqDB::CSeqDB().")
    );
    // A program trims the value of `-db` before it opens the database
    // (blast_args.cpp:2472-2476).
    assert_eq!(db_argument_name(" \t db one\r\n\x0b\x0c "), "db one");
    assert_eq!(db_argument_name(" \t "), "");
    assert_eq!(db_argument_name("\u{a0}db"), "\u{a0}db");
    let error = SeqDb::open_in("no_such_database", SeqType::Nucleotide, &search).unwrap_err();
    assert_eq!(
        error.ncbi_message(),
        Some(format!(
            "No alias or index file found for nucleotide database [no_such_database] in search \
             path [{}]. Please verify the spelling of the BLAST database and its molecule type.",
            search.as_str()
        ))
    );
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// (the texts of `eTooManyOpenFiles` and `eMemoryMappingFailure`, quoted at
// `SeqDbError::ncbi_message`).
/// Which errors are NCBI's own, and their texts.
#[test]
fn ncbis_errors_are_told_from_losats() {
    let path = PathBuf::from("/data/db.pin");
    let format = |message: &str| SeqDbError::Format {
        path: path.clone(),
        message: message.to_string(),
    };
    for message in [
        "Error: Not a valid version 4 or 5 database.",
        "Error: requested sequence type does not match DB.",
        "DB list contains both Version 4 and Version 5 dbs",
        "Illegal configuration: DB alias files are mutually recursive.",
        "No database names were found in alias file [/data/alias].",
    ] {
        assert_eq!(format(message).ncbi_message().as_deref(), Some(message));
    }
    assert_eq!(
        format("index file ends inside a field").ncbi_message(),
        None
    );
    assert_eq!(
        SeqDbError::NotFound("Error: File (/data/db.psq) not found.".to_string()).ncbi_message(),
        Some("Error: File (/data/db.psq) not found.".to_string())
    );
    let unsupported = SeqDbError::Unsupported {
        path: path.clone(),
        message: "x".to_string(),
    };
    assert_eq!(unsupported.ncbi_message(), None);
    assert_eq!(unsupported.ncbi_search_message(), None);
    assert_eq!(SeqDbError::Range("x".to_string()).ncbi_message(), None);
    // Inside the search NCBI replaces the reader's text (local_blast.cpp:240-252).
    assert_eq!(
        SeqDbError::NotFound("Error: File (/data/db.psq) not found.".to_string())
            .ncbi_search_message(),
        Some("Database memory map file error")
    );
    let io = |code: i32| std::io::Error::from_raw_os_error(code);
    assert_eq!(
        SeqDbError::Io {
            path: path.clone(),
            source: io(5)
        }
        .ncbi_message(),
        None
    );
    // An error without an operating system code is LOSAT's.
    assert_eq!(
        SeqDbError::Open {
            path: path.clone(),
            source: std::io::Error::other("x")
        }
        .ncbi_message(),
        None
    );
    #[cfg(unix)]
    {
        // EACCES (13), EMFILE (24) and ENFILE (23) on Linux and macOS.
        assert_eq!(
            SeqDbError::Open {
                path: path.clone(),
                source: io(13)
            }
            .ncbi_message()
            .as_deref(),
            Some(
                "Database memory map file error - CMemoryFile: Cannot memory map file \
                 '/data/db.pin':Permission denied"
            )
        );
        assert_eq!(
            SeqDbError::Open {
                path: path.clone(),
                source: io(13)
            }
            .ncbi_search_message(),
            Some("Database memory map file error")
        );
        for code in [23, 24] {
            let error = SeqDbError::Open {
                path: path.clone(),
                source: io(code),
            };
            let text = "Too many open files, please raise the open file limit";
            assert_eq!(error.ncbi_message().as_deref(), Some(text));
            assert_eq!(error.ncbi_search_message(), Some(text));
        }
    }
}

/// Opens a database as a search program would and prints what a program
/// would report, for checks by hand against NCBI that need a property of the
/// process (its current directory, `BLASTDB`, its limit of open files).
///
/// `LOSAT_SEQDB_OPEN_PROBE` is `prot|nucl,NAME`; the registry value comes
/// from `LOSAT_SEQDB_OPEN_PROBE_REGISTRY` when it is set. Every sequence and
/// title is read, so that all the files of all the volumes are opened.
#[test]
#[ignore = "needs LOSAT_SEQDB_OPEN_PROBE (a check by hand)"]
fn open_probe() {
    let spec = std::env::var("LOSAT_SEQDB_OPEN_PROBE").expect("LOSAT_SEQDB_OPEN_PROBE");
    let (kind, name) = spec.split_once(',').expect("prot|nucl,NAME");
    let seqtype = match kind {
        "prot" => SeqType::Protein,
        "nucl" => SeqType::Nucleotide,
        other => panic!("unknown sequence type {other}"),
    };
    let registry = std::env::var("LOSAT_SEQDB_OPEN_PROBE_REGISTRY").ok();
    let outcome = SearchPath::from_environment(registry.as_deref())
        .and_then(|search| SeqDb::open_in(db_argument_name(name), seqtype, &search))
        .and_then(|db| {
            for ordinal in 0..db.num_seqs() {
                db.read_sequence(ordinal)?;
                db.deflines(ordinal)?;
            }
            Ok(db)
        });
    match outcome {
        Ok(db) => println!(
            "probe: ok: {} sequences in {} volume(s), title {:?}",
            db.num_seqs(),
            db.volumes().len(),
            String::from_utf8_lossy(db.title())
        ),
        Err(error) => match error.ncbi_message() {
            Some(text) => println!("probe: BLAST Database error: {text}"),
            None => println!("probe: LOSAT: {error}"),
        },
    }
}
