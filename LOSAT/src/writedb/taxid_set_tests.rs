//! Tests of `-taxid` and `-taxid_map` (`taxid_set.rs`): the keys and the matches makeblastdb
//! 2.17.0 gave in makedb-v5 survey B fixture t01_taxmap (blastdbcmd `%T` of its database) and
//! the S3b extra fixtures x_taxmap_* (gate s3b/gate).

use super::*;
use crate::seqid::model::{Dbtag, ObjectId, TextseqId};

fn key(acc: &str) -> String {
    let mut warnings = Vec::new();
    let key = accession_to_key(acc.as_bytes(), &mut |w| warnings.push(w)).unwrap();
    assert!(warnings.is_empty(), "{acc}: {warnings:?}");
    String::from_utf8(key).unwrap()
}

fn map(text: &str) -> Result<TaxIdSet, TaxIdMapError> {
    let mut set = TaxIdSet::default();
    set.set_mapping_from_file(&mut text.as_bytes(), &mut |_| {})?;
    Ok(set)
}

fn defline(ids: Vec<ParsedSeqId>) -> BlastDefline {
    BlastDefline {
        ids,
        ..BlastDefline::default()
    }
}

fn refseq(acc: &str, version: i32) -> ParsedSeqId {
    ParsedSeqId::Other(TextseqId {
        accession: Some(acc.as_bytes().to_vec()),
        version: Some(version),
        ..TextseqId::default()
    })
}

/// `AccessionToKey`: a number is a GI, an accession loses its version and is upper-cased, a
/// general id keeps its tag only, a pdb id is `mol_chain`, a word that is no Seq-id is retried
/// as a local id.
#[test]
fn map_keys_are_makeblastdb_s() {
    for (acc, expected) in [
        ("123456", "123456"),
        ("0", ""),
        ("NP_000549.7", "NP_000549"),
        ("np_000551", "NP_000551"),
        ("sp|P69906|HBB_HUMAN", "P69906"),
        ("HBA_HUMAN", "HBA_HUMAN"),
        ("myseq1", "myseq1"),
        ("gnl|mydb|tag2", "tag2"),
        ("1ABC_A", "1ABC_A"),
        ("XP_011520000.1", "XP_011520000"),
        ("UniRef50_A0A5A9P0L4", "UniRef50_A0A5A9P0L4"),
        ("lcl|777", "777"),
        ("gi|42", "42"),
        ("a|b|c", ""),
    ] {
        assert_eq!(key(acc), expected, "{acc}");
    }
}

/// A key of more digits than an Int8 holds ends makeblastdb (x_taxmap_gibig, exit 255).
#[test]
fn a_gi_key_past_int8_is_an_exception() {
    let error = accession_to_key(b"99999999999999999999", &mut |_| {}).unwrap_err();
    assert!(
        error.ends_with(
            "line 852: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot \
             convert string '99999999999999999999' to Int8, overflow (m_Pos = 18)\n"
        ),
        "{error}"
    );
}

/// The taxid of a line: `StringToInt` with leading spaces allowed; its errors are NCBI's
/// (t03_taxmap_badnum, x_taxmap_trail, x_taxmap_minus, x_taxmap_overflow).
#[test]
fn taxids_convert_as_ncbi_s_string_to_int() {
    assert_eq!(string_to_taxid(b"9606"), Ok(9606));
    assert_eq!(string_to_taxid(b"+9"), Ok(9));
    assert_eq!(string_to_taxid(b"-8"), Ok(-8));
    assert_eq!(string_to_taxid(b"\x0c19"), Ok(19));
    assert_eq!(string_to_taxid(b"12\x0034"), Ok(12));
    for (text, line, function, message) in [
        (
            &b"abc"[..],
            861,
            "StringToInt8",
            "Cannot convert string 'abc' to Int8 (m_Pos = 0)",
        ),
        (
            b"12ab",
            870,
            "StringToInt8",
            "Cannot convert string '12ab' to Int8 (m_Pos = 2)",
        ),
        (
            b"-",
            861,
            "StringToInt8",
            "Cannot convert string '-' to Int8 (m_Pos = 1)",
        ),
        (
            b"3000000000",
            640,
            "StringToInt",
            "Cannot convert string '3000000000' to int, overflow (m_Pos = 0)",
        ),
    ] {
        let error = string_to_taxid(text).unwrap_err();
        assert!(
            error.ends_with(&format!(
                "line {line}: Error: (CStringException::eConvert) ncbi::NStr::{function}() - \
                 {message}\n"
            )),
            "{error}"
        );
        assert!(error.starts_with("NCBI C++ Exception:\n    T0 \""));
    }
}

/// The map file: lines of exactly two tokens split at each blank or tab, the last duplicate
/// key wins, other lines are ignored, a line of white space only is refused (makeblastdb
/// crashes there, x_taxmap_blank).
#[test]
fn map_lines_follow_ncbi_s_split() {
    let mut set =
        map("NP_000001 5\nNP_000001 6\nNP_000002  7\nNP_000003 8 9\nlonely\n\nNP_000004\t10\r\n")
            .unwrap();
    let mut deflines = DeflineSet {
        deflines: vec![
            defline(vec![refseq("NP_000001", 1)]),
            defline(vec![refseq("NP_000002", 1)]),
            defline(vec![refseq("NP_000003", 1)]),
            defline(vec![refseq("NP_000004", 2)]),
        ],
    };
    set.fix_tax_id(&mut deflines);
    let taxids: Vec<Option<i32>> = deflines.deflines.iter().map(|d| d.taxid).collect();
    assert_eq!(taxids, [Some(6), Some(0), Some(0), Some(10)]);
    assert!(set.has_ever_fixed_id());
    assert_eq!(
        map("NP_000001 5\n \nNP_000002 6\n").unwrap_err(),
        TaxIdMapError::BlankLine(2)
    );
}

/// `x_SelectBestTaxid`: the first id of the definition line whose key (or its part before
/// the first '.') is in the map; a general id matches by its tag whatever its database; no
/// match gives 0 and leaves `HasEverFixedId` false (t02_taxmap_nomatch); `-taxid` wins.
#[test]
fn deflines_take_the_taxid_of_their_first_matching_id() {
    let mut set = map("tag2 1009\nseq 1012\n777 1022\n").unwrap();
    assert!(!set.has_ever_fixed_id());
    let general = ParsedSeqId::General(Dbtag {
        db: b"otherdb".to_vec(),
        tag: ObjectId::Str(b"tag2".to_vec()),
    });
    let mut deflines = DeflineSet {
        deflines: vec![
            defline(vec![refseq("NP_9", 1)]),
            defline(vec![refseq("NP_8", 1), general]),
            defline(vec![ParsedSeqId::Local(ObjectId::Str(b"seq.5".to_vec()))]),
            defline(vec![ParsedSeqId::Gi(777)]),
        ],
    };
    set.fix_tax_id(&mut deflines);
    let taxids: Vec<Option<i32>> = deflines.deflines.iter().map(|d| d.taxid).collect();
    assert_eq!(taxids, [Some(0), Some(1009), Some(1012), Some(1022)]);
    assert!(set.has_ever_fixed_id());

    let mut none = map("unknown_id 5555\n").unwrap();
    let mut deflines = DeflineSet {
        deflines: vec![defline(vec![refseq("NP_9", 1)])],
    };
    none.fix_tax_id(&mut deflines);
    assert_eq!(deflines.deflines[0].taxid, Some(0));
    assert!(!none.has_ever_fixed_id());

    let mut global = TaxIdSet::with_global(9606);
    global.fix_tax_id(&mut deflines);
    assert_eq!(deflines.deflines[0].taxid, Some(9606));
    // An empty map matches "everything" (`m_TaxIdMap.empty()`).
    assert!(map("").unwrap().has_ever_fixed_id());
}

/// A map file whose reading fails after `data` (a directory, EIO): `read` gives the bytes, then
/// the error.
struct FailingReader {
    data: &'static [u8],
}

impl std::io::Read for FailingReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.data.is_empty() {
            return Err(std::io::Error::other("read error"));
        }
        let n = buf.len().min(self.data.len());
        buf[..n].copy_from_slice(&self.data[..n]);
        self.data = &self.data[n..];
        Ok(n)
    }
}

/// A read error ends the reading (`while (f && !f.eof())`, taxid_set.cpp:47): the entries read
/// so far stay, the line read up to the error is taken, no error is reported; an error at the
/// start is an empty map (makeblastdb 2.17.0 with a directory as `-taxid_map`: exit status 0,
/// no message; review B S2).
#[test]
fn a_read_error_keeps_the_map_read_so_far() {
    let read = |data: &'static [u8]| {
        let mut set = TaxIdSet::default();
        let mut reader = std::io::BufReader::with_capacity(4, FailingReader { data });
        set.set_mapping_from_file(&mut reader, &mut |_| {}).unwrap();
        set
    };
    let set = read(b"");
    assert!(set.map.is_empty());
    assert!(set.has_ever_fixed_id());
    let set = read(b"NP_000001.1 9606\nlowerid 10090\nP12345 5");
    assert_eq!(set.map.get(b"NP_000001"), Some(9606));
    assert_eq!(set.map.get(b"lowerid"), Some(10090));
    assert_eq!(set.map.get(b"P12345"), Some(5));
    assert!(!set.has_ever_fixed_id());
    // A directory, as std::fs opens it on Linux.
    #[cfg(target_os = "linux")]
    {
        let dir = std::env::temp_dir();
        let file = std::fs::File::open(&dir).unwrap();
        let mut set = TaxIdSet::default();
        set.set_mapping_from_file(&mut std::io::BufReader::new(file), &mut |_| {})
            .unwrap();
        assert!(set.map.is_empty());
    }
}

/// The map in the arena (review B S5): the last line of a key wins, a key is found by its bytes
/// only, and a UniRef50-like map takes less than 50 bytes per entry with its table.
#[test]
fn the_map_keeps_the_last_taxid_of_a_key_in_little_memory() {
    let set = map("A 1\nB 2\nA 3\n").unwrap();
    assert_eq!(set.map.get(b"A"), Some(3));
    assert_eq!(set.map.get(b"B"), Some(2));
    assert_eq!(set.map.get(b"a"), None);
    assert_eq!(set.map.len(), 2);
    let mut text = String::new();
    let n = 200_000u32;
    for i in 0..n {
        text.push_str(&format!("UniRef50_A{i:09} {}\n", i % 5000));
    }
    let set = map(&text).unwrap();
    assert_eq!(set.map.len(), n as usize);
    for i in (0..n).step_by(997) {
        assert_eq!(
            set.map.get(format!("UniRef50_A{i:09}").as_bytes()),
            Some((i % 5000) as i32)
        );
    }
    // 19-byte keys: the bytes, one length byte, the 4-byte value and the slots.
    let per_entry = set.map.heap_bytes() as f64 / f64::from(n);
    assert!(per_entry < 50.0, "{per_entry} bytes per entry");
}
