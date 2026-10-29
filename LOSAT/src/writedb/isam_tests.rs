//! Tests of the ISAM writers (`isam.rs`) against the files makeblastdb 2.17.0 wrote for the
//! same ids (makedb-v5 survey B fixtures p14_large, p22_oid_order, p23_nonascii, p24_gi_range,
//! `-parse_seqids -blastdb_version 4`; the S3b probe of NCBI's writedb unit-test ids).

use super::*;
use crate::seqid::model::{ObjectId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId};

/// A directory of this test, removed when dropped.
struct Dir(std::path::PathBuf);

impl Dir {
    fn new(label: &str) -> Self {
        use std::sync::atomic::{AtomicU32, Ordering};
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-isam-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn base(&self) -> String {
        self.0.join("db").to_string_lossy().into_owned()
    }

    fn read(&self, name: &str) -> Vec<u8> {
        std::fs::read(self.0.join(name)).unwrap()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

fn sp(acc: &str, name: &str) -> ParsedSeqId {
    ParsedSeqId::Swissprot(TextseqId {
        name: Some(name.as_bytes().to_vec()),
        accession: Some(acc.as_bytes().to_vec()),
        release: Some(b"reviewed".to_vec()),
        version: None,
    })
}

fn refseq(acc: &str, version: i32) -> ParsedSeqId {
    ParsedSeqId::Other(TextseqId {
        accession: Some(acc.as_bytes().to_vec()),
        version: Some(version),
        ..TextseqId::default()
    })
}

fn local(text: &[u8]) -> ParsedSeqId {
    ParsedSeqId::Local(ObjectId::Str(text.to_vec()))
}

fn records(list: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    for record in list {
        out.extend_from_slice(record);
        out.push(b'\n');
    }
    out
}

fn acc_isam(dir: &Dir) -> StringIsam {
    StringIsam::new(
        IsamType::Acc,
        &dir.base(),
        true,
        0,
        3_000_000_000,
        StringIsamScratch::default(),
    )
}

/// p22_oid_order: the OID is compared as decimal text (`same_name\x0210` before
/// `same_name\x022`), and one sample page.
#[test]
fn records_sort_with_the_oid_as_decimal_text() {
    let dir = Dir::new("p22");
    let mut isam = acc_isam(&dir);
    for oid in 0..12u32 {
        isam.add_ids(oid, &[sp(&format!("Q{oid:05}"), "SAME_NAME")])
            .unwrap();
    }
    isam.close().unwrap();
    assert_eq!(
        dir.read("db.00.psi"),
        hex("0000000100000002000001000000001800000001000000400000100000000000000000000000000000000100000000340000003d713030303030023000")
    );
    let mut expected: Vec<Vec<u8>> = (0..12)
        .map(|oid| format!("q{oid:05}\x02{oid}").into_bytes())
        .collect();
    for oid in [0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9] {
        expected.push(format!("same_name\x02{oid}").into_bytes());
    }
    let list: Vec<&[u8]> = expected.iter().map(Vec::as_slice).collect();
    assert_eq!(dir.read("db.00.psd"), records(&list));
}

/// p23_nonascii: a byte of 0x80 or more in the first 6 bytes sorts before ASCII (signed
/// `char`), keys are lower-cased, an accession gives `acc` and `acc.ver`.
#[test]
fn a_high_byte_in_the_prefix_sorts_first() {
    let dir = Dir::new("p23");
    let mut isam = acc_isam(&dir);
    isam.add_ids(0, &[local(b"zzz")]).unwrap();
    isam.add_ids(1, &[local(b"\xe9abc")]).unwrap();
    isam.add_ids(2, &[local(b"aaa")]).unwrap();
    isam.add_ids(3, &[refseq("NP_000001", 1)]).unwrap();
    isam.close().unwrap();
    assert_eq!(
        dir.read("db.00.psi"),
        hex("00000001000000020000004c000000080000000100000040000010000000000000000000000000000000004c000000340000003be9616263023100")
    );
    assert_eq!(
        dir.read("db.00.psd"),
        records(&[
            b"\xe9abc\x021",
            b"aaa\x022",
            b"lcl|\xe9abc\x021",
            b"lcl|aaa\x022",
            b"lcl|zzz\x020",
            b"np_000001\x023",
            b"np_000001.1\x023",
            b"zzz\x020",
        ])
    );
}

/// p24_gi_range: GIs sort as Int8 and are written as Int4 (2^31 negative); `.pog` holds the
/// first GI of each OID.
#[test]
fn gis_sort_as_int8_and_are_written_as_int4() {
    let dir = Dir::new("p24");
    let mut gi = NumericIsam::new(IsamType::Gi, &dir.base(), true, 0, 3_000_000_000);
    let mut pog = GiIndex::new(&dir.base(), true, 0, 3_000_000_000);
    let oids: [Vec<ParsedSeqId>; 4] = [
        vec![ParsedSeqId::Gi(0), refseq("NP_000010", 1)],
        vec![ParsedSeqId::Gi(2_147_483_648), refseq("NP_000011", 1)],
        vec![ParsedSeqId::Gi(4_294_967_295)],
        vec![ParsedSeqId::Gi(99)],
    ];
    for (oid, ids) in oids.iter().enumerate() {
        gi.add_gis(oid as u32, ids);
        pog.add_gi(ids.iter().find_map(ParsedSeqId::gi).unwrap_or(-1));
    }
    gi.close().unwrap();
    pog.close().unwrap();
    assert_eq!(
        dir.read("db.00.pni"),
        hex("0000000100000000000000200000000400000001000001000000000000000000000000000000000000000000ffffffff00000000")
    );
    assert_eq!(
        dir.read("db.00.pnd"),
        hex("000000000000000000000063000000038000000000000001ffffffff00000002")
    );
    assert_eq!(
        dir.read("db.00.pog"),
        hex("00000001000000000000000400000004000000000000000000000000000000000000000080000000ffffffff00000063")
    );
}

/// p14_large: 672 string records (11 sample pages, a pdb mol shared by 12 OIDs, one name on
/// 12 OIDs) and 300 GIs (2 pages); the files are makeblastdb's (FNV-1a of each).
#[test]
fn several_sample_pages_are_makeblastdb_s() {
    let dir = Dir::new("p14");
    let base = dir.base();
    let mut acc = acc_isam(&dir);
    let mut gi = NumericIsam::new(IsamType::Gi, &base, true, 0, 3_000_000_000);
    let mut pog = GiIndex::new(&base, true, 0, 3_000_000_000);
    let mut oid = 0u32;
    let mut add = |ids: Vec<ParsedSeqId>| {
        acc.add_ids(oid, &ids).unwrap();
        gi.add_gis(oid, &ids);
        pog.add_gi(ids.iter().find_map(ParsedSeqId::gi).unwrap_or(-1));
        oid += 1;
    };
    for i in 0..300 {
        add(vec![
            ParsedSeqId::Gi(100_000 + i * 7),
            refseq(&format!("NP_{i:06}"), 1),
        ]);
    }
    for chain in b"ABCDEFGHIJKL" {
        add(vec![ParsedSeqId::Pdb(PdbSeqId {
            mol: b"9BIG".to_vec(),
            chain: Some(i32::from(*chain)),
            rel: None,
            chain_id: Some(vec![*chain]),
        })]);
    }
    for i in 0..12 {
        add(vec![sp(&format!("Q{i:05}"), "SHARED_NAME")]);
    }
    acc.close().unwrap();
    gi.close().unwrap();
    pog.close().unwrap();
    for (name, len, hash) in [
        ("db.00.psi", 279, 0x9044_1e0e_aadf_62a8u64),
        ("db.00.psd", 9656, 0x1e09_ce87_d44b_f267),
        ("db.00.pni", 60, 0x7f1e_cb8d_1197_1b3b),
        ("db.00.pnd", 2400, 0xda20_0286_98d7_04c0),
        ("db.00.pog", 1328, 0x2710_8af8_747b_c3c6),
    ] {
        let bytes = dir.read(name);
        assert_eq!((bytes.len(), fnv1a(&bytes)), (len, hash), "{name}");
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1500-1541
// ```c++
// BOOST_AUTO_TEST_CASE(PDBIdLowerCase)
// {
//     ...
//         I1("pdb|3E3Q|b"), T1("Lower case chain b");
//     ...
//         CRef<CSeqDB> rd ...
//         rd.AccessionToOids("3e3q_b", oids);
//
//         BOOST_REQUIRE(oids.size() == 1);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1318-1359
// ```c++
// BOOST_AUTO_TEST_CASE(UsPatId)
// {
//
//     CRef<CSeq_id> seqid(new CSeq_id("pat|us|123|456"));
//     ...
//     CSeqDB seqdb("uspatid", CSeqDB::eProtein);
//     int oid(-1);
//     bool found = seqdb.SeqidToOid(*seqid, oid);
//
//     BOOST_REQUIRE_EQUAL(found, true);
//     BOOST_REQUIRE_EQUAL(oid,   0);
// ```
// NCBI's cases read the database back with CSeqDB (the ISAM lookups the files allow). Here
// the files are compared with makeblastdb's for the same two ids (s3b/probe/ncbiut, where
// `blastdbcmd -entry '3e3q_b,pat|us|123|456'` finds OIDs 0 and 1); the ledger of NCBI unit
// cases does not cover the seqdb_writer tests, so this is no ledger port.
#[test]
fn ncbi_s_pdb_lower_case_chain_and_us_patent_keys() {
    let dir = Dir::new("ncbiut");
    let mut isam = acc_isam(&dir);
    isam.add_ids(
        0,
        &[ParsedSeqId::Pdb(PdbSeqId {
            mol: b"3E3Q".to_vec(),
            chain: Some(i32::from(b'b')),
            rel: None,
            chain_id: Some(b"b".to_vec()),
        })],
    )
    .unwrap();
    isam.add_ids(
        1,
        &[ParsedSeqId::Patent(PatentSeqId {
            seqid: 456,
            cit: crate::seqid::model::IdPat {
                country: b"us".to_vec(),
                id: PatentNumber::Number(b"123".to_vec()),
                doc_type: None,
            },
        })],
    )
    .unwrap();
    isam.close().unwrap();
    assert_eq!(
        dir.read("db.00.psd"),
        records(&[
            b"3e3q\x020",
            b"3e3q b\x020",
            b"3e3q|b\x020",
            b"pat|us|123|456\x021",
            b"pdb|3e3q|b\x020",
        ])
    );
    assert_eq!(
        dir.read("db.00.psi"),
        hex("0000000100000002000000370000000500000001000000400000100000000000000000000000000000000037000000340000003b33653371023000")
    );
}

/// A key of 256 bytes is written; one of 257 bytes, which makeblastdb's `char buf[256]` cannot
/// hold (it stops with a buffer overflow), is refused.
#[test]
fn keys_longer_than_256_bytes_are_refused() {
    let dir = Dir::new("long");
    let mut isam = acc_isam(&dir);
    isam.add_string_data(0, &[b'k'; 256]).unwrap();
    let error = isam.add_string_data(1, &[b'k'; 257]).unwrap_err();
    assert!(error.to_string().contains("257 bytes long"), "{error}");
}

/// `CanFit`: the data size and `(num + 1) * bytes_per_elem` must stay under the limit.
#[test]
fn can_fit_is_strict() {
    assert!(can_fit(0, 1, 1024, 2049));
    assert!(!can_fit(0, 1, 1024, 2048));
    assert!(can_fit(100, 0, 8, 109));
    assert!(!can_fit(100, 0, 8, 108));
}

/// Keys of one OID are added once (`m_OidStringData`), after lower-casing.
#[test]
fn a_key_is_added_once_per_oid() {
    let dir = Dir::new("dedup");
    let mut isam = acc_isam(&dir);
    isam.add_string_data(0, b"ABC").unwrap();
    isam.add_string_data(0, b"abc").unwrap();
    isam.add_string_data(1, b"abc").unwrap();
    for i in 0..40u32 {
        isam.add_string_data(2, format!("K{}", i % 37).as_bytes())
            .unwrap();
    }
    isam.close().unwrap();
    let psd = dir.read("db.00.psd");
    assert_eq!(psd.iter().filter(|&&b| b == b'\n').count(), 2 + 37);
    assert!(psd.starts_with(b"abc\x020\nabc\x021\n"));
}

/// The memory of the string index for UniRef50-like local ids (`lcl|UniRef50_<10>` and the
/// bare id: two records per sequence): entries of 16 bytes and the key bytes past the first 8.
#[test]
fn string_index_memory_per_sequence() {
    let dir = Dir::new("mem");
    let mut isam = acc_isam(&dir);
    let n = 100_000u32;
    for oid in 0..n {
        isam.add_ids(oid, &[local(format!("UniRef50_A{oid:09}").as_bytes())])
            .unwrap();
    }
    let per_sequence = isam.record_bytes() as f64 / f64::from(n);
    // 2 x 16 (entries) + (23 - 8 + 1) + (19 - 8 + 1) = 60; the chunks add at most 1 MiB each
    // to it (one entry chunk and one tail chunk being filled).
    assert_eq!(per_sequence, 60.0);
}
