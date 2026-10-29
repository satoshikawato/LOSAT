//! Tests of `CWriteDB_LMDB::x_InsertEntry` (`WriteDbLmdb::insert_entries`) and of the compact
//! `m_list` (`KeyValuePair`): the elements each kind of Seq-id gives, their order, and the
//! memory per key. The expected keys are those of makeblastdb 2.17.0's `.pdb` and `.pos` files
//! for the same ids (makedb-v5 S3c extra fixture x_keys_all, `-parse_seqids -blastdb_version 5`;
//! the byte-for-byte checks of whole databases are in `writedb::tests::parse_seqids_v5`).

use super::*;
use crate::seqid::accguide::parse_flags;
use crate::seqid::parse::seq_id_from_text;
use crate::writedb::platform::LmdbPlatform;

fn id(text: &str) -> ParsedSeqId {
    seq_id_from_text(text.as_bytes(), parse_flags::DEFAULT, &mut |_| {})
        .unwrap_or_else(|e| panic!("{text}: {e:?}"))
}

/// A writer in a directory of its own (removed when dropped).
struct Writer {
    lmdb: Option<WriteDbLmdb>,
    dir: PathBuf,
}

impl Writer {
    fn new(label: &str) -> Self {
        use std::sync::atomic::{AtomicU32, Ordering};
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-lmdb-keys-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let lmdb = WriteDbLmdb::new(
            &dir.join("db.pdb"),
            700_000_000,
            &LmdbPlatform::LINUX_X86_64,
            &LmdbEnvSettings::default(),
        )
        .unwrap();
        Writer {
            lmdb: Some(lmdb),
            dir,
        }
    }

    fn lmdb(&mut self) -> &mut WriteDbLmdb {
        self.lmdb.as_mut().unwrap()
    }

    /// The elements of `m_list` in insertion order: (id, OID, `saveToOidList`).
    fn list(&mut self) -> Vec<(String, u32, bool)> {
        let lmdb = self.lmdb();
        let mut out = Vec::new();
        let mut key = Vec::new();
        for kv in &lmdb.list {
            kv.id_into(&lmdb.ids, &mut key);
            out.push((
                String::from_utf8(key.clone()).unwrap(),
                kv.oid(),
                kv.save_to_oid_list(),
            ));
        }
        out
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        drop(self.lmdb.take());
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn entries(list: &[(&str, bool)], oid: u32) -> Vec<(String, u32, bool)> {
    list.iter()
        .map(|(id, save)| (id.to_string(), oid, *save))
        .collect()
}

// The keys of each kind of Seq-id (writedb_lmdb.cpp:120-196, quoted at `insert_seq_id`), as
// makeblastdb's acc2oid keys and `.pos` ids of x_keys_all show them: `saveToOidList` is true
// for the ids that `.pos` lists.
#[test]
fn each_kind_of_seq_id_gives_makeblastdb_s_keys() {
    let cases: &[(&str, &[(&str, bool)])] = &[
        ("pir||S12345", &[("pir||S12345", true)]),
        ("prf||1234567A", &[("prf||1234567A", true)]),
        (
            "pdb|1abc|A",
            &[("1abc", false), ("1ABC", false), ("1abc_A", true)],
        ),
        ("pdb|2XYZ|C", &[("2XYZ", false), ("2XYZ_C", true)]),
        (
            "pdb|3def|",
            &[("3def", false), ("3DEF", false), ("3def", true)],
        ),
        (
            "pdb|4ghi|AA",
            &[("4ghi", false), ("4GHI", false), ("4ghi_AA", true)],
        ),
        ("ref|NP_000551|", &[("NP_000551", true)]),
        (
            "ref|NP_000552.3|",
            &[("NP_000552", false), ("NP_000552.3", true)],
        ),
        (
            "gb|aab12345.1|",
            &[("AAB12345", false), ("AAB12345.1", true)],
        ),
        ("sp|P69905|HBA_HUMAN", &[("P69905", true)]),
        ("tr|A0A000|A0A000_YEAST", &[("A0A000", true)]),
        ("lcl|lowerid", &[("lowerid", true), ("LOWERID", false)]),
        ("lcl|UPPERID", &[("UPPERID", true)]),
        (
            "lcl|MiXeD_1.2",
            &[("MiXeD_1.2", true), ("MIXED_1.2", false)],
        ),
        ("lcl|12345", &[("12345", true)]),
        (
            "gnl|mydb|tag1",
            &[("mydb:tag1", true), ("MYDB:TAG1", false)],
        ),
        ("gnl|MYDB|TAG2", &[("MYDB:TAG2", true)]),
        ("gnl|Db|77", &[("Db:77", true), ("DB:77", false)]),
        ("pat|US|1234567|1", &[("US1234567_1", true)]),
        ("pgp|EP|0238993|7", &[("EP0238993_7", true)]),
        ("gi|555", &[]),
    ];
    for (oid, (text, expected)) in cases.iter().enumerate() {
        let mut w = Writer::new("kinds");
        let added = w.lmdb().insert_entries(&[id(text)], oid as u32).unwrap();
        assert_eq!(added, expected.len(), "{text}");
        assert_eq!(w.list(), entries(expected, oid as u32), "{text}");
    }
}

// `InsertEntries` takes the ids of every definition line of a sequence in order; a GI adds
// nothing and two ids of one sequence may give the same key (an equal (id, OID) pair, which
// the transaction skips: `put_batch`).
#[test]
fn the_ids_of_a_sequence_go_in_order() {
    let mut w = Writer::new("order");
    let ids = [
        id("gi|555"),
        id("ref|NP_000901.1|"),
        id("ref|NP_000901.2|"),
        id("lcl|dupkey"),
    ];
    assert_eq!(w.lmdb().insert_entries(&ids, 31).unwrap(), 6);
    assert_eq!(
        w.list(),
        entries(
            &[
                ("NP_000901", false),
                ("NP_000901.1", true),
                ("NP_000901", false),
                ("NP_000901.2", true),
                ("dupkey", true),
                ("DUPKEY", false),
            ],
            31
        )
    );
    // The bytes of each Seq-id are stored once.
    assert_eq!(w.lmdb().ids, b"NP_000901.1NP_000901.2dupkey".to_vec());
}

// SKeyValuePair::cmp_key (writedb_lmdb.hpp:116-127): the strings bytewise (upper case before
// lower case, a prefix first), then the OIDs by their bytes from the least significant one
// (OID 256 = 00 01 00 00 sorts before OID 1 = 01 00 00 00).
#[test]
fn the_list_sorts_as_std_string_then_oid_bytes() {
    let mut w = Writer::new("sort");
    let lmdb = w.lmdb();
    lmdb.insert_entries(&[id("lcl|abc")], 1).unwrap();
    lmdb.insert_entries(&[id("lcl|ABC")], 2).unwrap();
    lmdb.insert_entries(&[id("lcl|ab")], 3).unwrap();
    lmdb.insert_entries(&[id("lcl|abc")], 256).unwrap();
    lmdb.insert_entries(&[id("lcl|aBc")], 4).unwrap();
    let ids = &lmdb.ids;
    lmdb.list.sort_unstable_by(|a, b| kv_cmp(ids, *a, *b));
    let sorted: Vec<(String, u32)> = w.list().into_iter().map(|(i, o, _)| (i, o)).collect();
    let expect = [
        ("AB", 3),
        ("ABC", 256),
        ("ABC", 1),
        ("ABC", 2),
        ("ABC", 4),
        ("aBc", 4),
        ("ab", 3),
        ("abc", 256),
        ("abc", 1),
    ];
    assert_eq!(
        sorted,
        expect
            .iter()
            .map(|(i, o)| (i.to_string(), *o))
            .collect::<Vec<_>>()
    );
}

// Memory per key: 12 bytes per element of `m_list` (NCBI: a 40-byte element and the string's
// heap block); the bytes of a Seq-id once. A UniRef50-like local id (two keys) costs its 19
// bytes and 24 bytes of elements; a RefSeq accession with a version 14 + 24.
#[test]
fn keys_take_twelve_bytes_and_the_id_once() {
    assert_eq!(std::mem::size_of::<KeyValuePair>(), 12);
    let mut w = Writer::new("memory");
    let n = 10_000u32;
    for oid in 0..n {
        let text = format!("lcl|UniRef50_A0A{oid:07}");
        w.lmdb().insert_entries(&[id(&text)], oid).unwrap();
    }
    let lmdb = w.lmdb();
    assert_eq!(lmdb.list.len(), 2 * n as usize);
    assert_eq!(lmdb.ids.len(), 19 * n as usize);
    let per_sequence = (lmdb.list.len() * 12 + lmdb.ids.len()) as f64 / n as f64;
    assert!(per_sequence <= 43.0, "{per_sequence}");
}

// An id longer than the 22 bits of a length is refused (makeblastdb would fail on it later:
// an LMDB key holds at most 511 bytes).
#[test]
fn an_id_past_the_word_s_length_is_refused() {
    assert!(KeyValuePair::new(0, KV_MAX_LEN as usize, 0, false, true).is_ok());
    assert!(KeyValuePair::new(0, KV_MAX_LEN as usize + 1, 0, false, true).is_err());
    assert!(KeyValuePair::new(KV_MAX_OFF as usize, 1, 0, false, true).is_ok());
    assert!(KeyValuePair::new(KV_MAX_OFF as usize + 1, 1, 0, false, true).is_err());
    let kv = KeyValuePair::new(123_456_789_012, 4_000_000, 2_147_483_647, true, false).unwrap();
    assert_eq!(
        (
            kv.off(),
            kv.len(),
            kv.oid(),
            kv.upper(),
            kv.save_to_oid_list()
        ),
        (123_456_789_012, 4_000_000, 2_147_483_647, true, false)
    );
}
