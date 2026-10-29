//! Tests of the LMDB files of version 5 against makeblastdb 2.17.0: the inputs are made in the
//! test code (the FASTA fixtures of the makedb-v5 survey A, whose ids and tax ids follow simple
//! rules), the expected files are makeblastdb's (`LOSAT/tests/writedb_lmdb/`, or their FNV-1a
//! digests for the larger ones). The tests run no NCBI program.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use super::{
    lmdb_file_name, max_entry_per_txn, string_to_int, string_to_uint, string_to_uint8, txn_batches,
    LmdbEnvSettings, WriteDbLmdb, WriteDbTaxId, DEFAULT_LMDB_MAP_SIZE, DEFAULT_TAXID_MAP_SIZE,
};

/// A directory of this test run, removed when the value is dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-lmdb-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The `m_list` entries of one OID: (id, `saveToOidList`).
type Entries = Vec<(Vec<u8>, bool)>;

/// The expected size and FNV-1a digest of a file (extension, size, digest).
type Digest = (&'static str, usize, u64);

/// The logical content of a database's LMDB files, as `CWriteDB_Impl` hands it over.
#[derive(Default, Debug)]
struct Content {
    /// 'p' (protein: `.pdb`) or 'n' (nucleotide: `.ndb`).
    dbtype: char,
    /// `BLASTDB_LMDB_MAP_SIZE` (0: the default map sizes).
    map_size: u64,
    max_txn: Option<String>,
    volumes: Vec<(String, u32)>,
    /// Per OID: the `m_list` entries (id, `saveToOidList`) and the tax ids.
    oids: Vec<(Entries, Vec<i32>)>,
}

/// Write the LMDB and lookup files of `c` as `<dir>/db.*`, in the order of
/// `CWriteDB_Impl::x_Publish` and `CWriteDB_Impl::Close`.
fn write(c: &Content, dir: &Path) -> Result<(), crate::writedb::WriteDbError> {
    let t = c.dbtype;
    let pdb = dir.join(format!("db.{t}db"));
    let ptf = lmdb_file_name(&pdb, "tf");
    let settings = LmdbEnvSettings {
        max_lmdb_txn_entry: c.max_txn.clone(),
        ..Default::default()
    };
    let (m1, m2) = if c.map_size > 0 {
        (c.map_size, c.map_size)
    } else {
        (DEFAULT_LMDB_MAP_SIZE, DEFAULT_TAXID_MAP_SIZE)
    };
    let mut lmdb = WriteDbLmdb::new(&pdb, m1, 4096, &settings)?;
    let mut tax = WriteDbTaxId::new(&ptf, m2, 4096, &settings)?;
    for (oid, (entries, taxids)) in c.oids.iter().enumerate() {
        for (id, save) in entries {
            lmdb.insert_entry(id, oid as u32, *save);
        }
        tax.insert_entries(taxids, oid as u32);
    }
    let names: Vec<String> = c.volumes.iter().map(|v| v.0.clone()).collect();
    let counts: Vec<u32> = c.volumes.iter().map(|v| v.1).collect();
    lmdb.insert_volumes_info(&names, &counts)?;
    lmdb.close()?;
    tax.close()
}

/// The first offset where two files differ.
fn first_difference(a: &[u8], b: &[u8]) -> Option<usize> {
    if a == b {
        return None;
    }
    Some(
        a.iter()
            .zip(b)
            .position(|(x, y)| x != y)
            .unwrap_or(a.len().min(b.len())),
    )
}

/// FNV-1a, 64 bits.
fn fnv1a(b: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &x in b {
        h ^= x as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

// ---- the inputs of survey A (gen_fixtures.py), made here ---------------------------------

const SEVEN: [i32; 7] = [9606, 10090, 7227, 562, 4932, 3702, 16_777_216];

/// A text Seq-id with a version: `x_InsertEntry` adds the id without the version and the
/// id with it (only the latter goes into `.pos`).
fn text_id(acc: &str, version: u32) -> Entries {
    vec![
        (acc.as_bytes().to_vec(), false),
        (format!("{acc}.{version}").into_bytes(), true),
    ]
}

/// `f4_ids{n}.fa`: `>XP_%09d.1` with `i*7+13`; tax ids from `-taxid_map` (`distinct`:
/// `1 + 37 i`, `seven`: SEVEN[i % 7]) or 0.
fn f4(n: usize, taxid: Option<&str>) -> Content {
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), n as u32)],
        ..Default::default()
    };
    for i in 0..n {
        let t = match taxid {
            None => vec![0],
            Some("distinct") => vec![1 + i as i32 * 37],
            Some(_) => vec![SEVEN[i % 7]],
        };
        c.oids
            .push((text_id(&format!("XP_{:09}", i * 7 + 13), 1), t));
    }
    c
}

/// `f5_dup{nv}.fa` (`gen_fixtures.py` main: 20 fillers follow) or the `extras` version
/// (500, 501: no fillers): `ref|NP_000001.v|` for v = 1..nv.
fn f5(nv: u32, fillers: bool) -> Content {
    let mut c = Content {
        dbtype: 'p',
        ..Default::default()
    };
    for v in 1..=nv {
        c.oids.push((text_id("NP_000001", v), vec![0]));
    }
    if fillers {
        for i in 0..20 {
            c.oids.push((text_id(&format!("XP_{i:06}"), 1), vec![0]));
        }
    }
    c.volumes = vec![("db".to_string(), c.oids.len() as u32)];
    c
}

/// `f5_pdb.fa`: `pdb|<1 + m % 9>X<m:02>|<chain>` for 40 molecules with 1-8 chains; the keys
/// are the molecule (not in `.pos`) and `<mol>_<chain>`.
fn f5pdb() -> Content {
    let mut c = Content {
        dbtype: 'p',
        ..Default::default()
    };
    for m in 0..40 {
        for ch in "ABCDEFGH".chars().take(1 + m % 8) {
            let mol = format!("{}X{:02}", 1 + m % 9, m);
            c.oids.push((
                vec![
                    (mol.clone().into_bytes(), false),
                    (format!("{mol}_{ch}").into_bytes(), true),
                ],
                vec![0],
            ));
        }
    }
    c.volumes = vec![("db".to_string(), c.oids.len() as u32)];
    c
}

/// The files of the expected database `name` in `LOSAT/tests/writedb_lmdb/`.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/writedb_lmdb")
        .join(name)
}

/// Write `c` and compare every file with the fixture files `<name>.<ext>`, or with the
/// (size, FNV-1a) digests given.
fn check(name: &str, c: &Content, digests: &[Digest]) {
    let dir = TempDir::new(name);
    write(c, &dir.0).unwrap();
    let t = c.dbtype;
    let exts = ["db", "tf", "os", "ot", "to"];
    let mut compared = 0;
    for ext in exts {
        let ext = format!("{t}{ext}");
        let ours = std::fs::read(dir.0.join(format!("db.{ext}"))).ok();
        if let Some(&(_, size, h)) = digests.iter().find(|d| d.0 == ext) {
            let ours = ours.unwrap_or_else(|| panic!("{name}: no .{ext} written"));
            assert_eq!((ours.len(), fnv1a(&ours)), (size, h), "{name}: .{ext}");
            compared += 1;
            continue;
        }
        let expected = std::fs::read(fixture(&format!("{name}.{ext}"))).ok();
        match (ours, expected) {
            (Some(o), Some(e)) => {
                assert_eq!(
                    first_difference(&o, &e),
                    None,
                    "{name}: .{ext} ({} vs {} bytes)",
                    o.len(),
                    e.len()
                );
                compared += 1;
            }
            (None, None) => {}
            (o, e) => panic!(
                "{name}: .{ext} written {} expected {}",
                o.is_some(),
                e.is_some()
            ),
        }
    }
    assert!(compared >= 2, "{name}: nothing compared");
    assert!(!dir.0.join("db.pdb-lock").exists() && !dir.0.join("db.ptf-lock").exists());
}

#[test]
fn txn_batches_follow_unsigned_int_arithmetic() {
    assert_eq!(txn_batches(5, 2).unwrap(), vec![(0, 2), (2, 4), (4, 5)]);
    assert_eq!(txn_batches(5, 40_000).unwrap(), vec![(0, 5)]);
    // MAX_LMDB_TXN_ENTRY=-1 is 4294967295 as `unsigned int`: one batch.
    assert_eq!(txn_batches(5, u32::MAX).unwrap(), vec![(0, 5)]);
    assert!(txn_batches(0, 7).unwrap().is_empty());
    assert!(txn_batches(3, 0).is_err());
    // 3000000000 + 3000000000 wraps in `unsigned int`.
    assert!(txn_batches(4_000_000_000, 3_000_000_000).is_err());
    assert_eq!(
        txn_batches(3_000_000_000, 2_000_000_000).unwrap(),
        vec![(0, 2_000_000_000), (2_000_000_000, 3_000_000_000)]
    );
}

/// `f9_multi.fa` with `f9_multi.map` (`t9`) or `f9_neg.map` (`t9neg`): three OIDs with
/// several deflines each.
fn f9(taxids: [Vec<i32>; 3]) -> Content {
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), 3)],
        ..Default::default()
    };
    let ids: [&[u32]; 3] = [&[1, 2], &[3], &[4, 5, 6]];
    for (k, t) in taxids.into_iter().enumerate() {
        let mut e = Vec::new();
        for &n in ids[k] {
            e.extend(text_id(&format!("NP_{n:06}"), 1));
        }
        c.oids.push((e, t));
    }
    c
}

/// The `.ptf` of a database without tax ids: one key, 0 (the same 16384 bytes for any size).
const PTF_TAXID0: Digest = ("ptf", 16384, 0xc29f655d5a3815f8);

fn with_max(mut c: Content, max: &str) -> Content {
    c.max_txn = Some(max.to_string());
    c
}

// The expected files are those of makeblastdb 2.17.0 (`-blastdb_version 5`, page size 4096)
// in the makedb-v5 survey A (`survey/A/run_fx.sh`).

#[test]
fn one_volume_without_parse_seqids_has_no_acc2oid() {
    // f1np: 3 proteins, no -parse_seqids: volinfo and volname only; .ptf holds tax id 0.
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), 3)],
        ..Default::default()
    };
    c.oids = (0..3).map(|_| (Vec::new(), vec![0])).collect();
    check(
        "f1np",
        &c,
        &[
            ("pdb", 20480, 0x4ecbae656a52eaf5),
            ("pot", 44, 0x4f24a59f1f248446),
            ("ptf", 16384, 0xc29f655d5a3815f8),
            ("pto", 16, 0xb623a479d2984da5),
        ],
    );
}

#[test]
fn small_acc2oid_reuses_the_main_leaf_buffer() {
    // f1p: txn 2 pops the buffers of txn 1 (pages 4, 3, 2) for pages 5, 6 and 7: page 7 (the
    // FREE_DBI leaf) still shows the volname/volinfo nodes of txn 1 in its gap.
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), 3)],
        ..Default::default()
    };
    for i in 1..=3 {
        c.oids.push((text_id(&format!("XP_{i:06}"), 1), vec![0]));
    }
    check("f1p", &c, &[PTF_TAXID0]);
}

#[test]
fn map_size_from_blastdb_lmdb_map_size() {
    // f1pmap: BLASTDB_LMDB_MAP_SIZE=123456789 changes only the mapsize of the meta pages.
    let mut c = Content {
        dbtype: 'p',
        map_size: 123_456_789,
        volumes: vec![("db".to_string(), 3)],
        ..Default::default()
    };
    for i in 1..=3 {
        c.oids.push((text_id(&format!("XP_{i:06}"), 1), vec![0]));
    }
    let dir = TempDir::new("f1pmap");
    write(&c, &dir.0).unwrap();
    let ours = std::fs::read(dir.0.join("db.pdb")).unwrap();
    let mut expected = std::fs::read(fixture("f1p.pdb")).unwrap();
    for off in [32, 4096 + 32] {
        expected[off..off + 8].copy_from_slice(&123_456_789u64.to_le_bytes());
    }
    assert_eq!(first_difference(&ours, &expected), None);
}

#[test]
fn fifty_proteins_and_a_nucleotide_database() {
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), 50)],
        ..Default::default()
    };
    for i in 0..50 {
        c.oids.push((text_id(&format!("XP_{i:06}"), 1), vec![0]));
    }
    check(
        "f2p",
        &c,
        &[
            ("pdb", 32768, 0xe432ad80dd5f9835),
            ("pos", 1008, 0x386fbd3b423ea42b),
            ("pot", 608, 0xf2964a90dcd0b8c4),
            ("ptf", 16384, 0xc29f655d5a3815f8),
            ("pto", 204, 0x1a5f70c525627df6),
        ],
    );
    // f7n: `>NM_%06d.1` nucleotides: the same layout in .ndb/.ntf/.nos/.not/.nto.
    let mut n = Content {
        dbtype: 'n',
        volumes: vec![("db".to_string(), 5)],
        ..Default::default()
    };
    for i in 0..5 {
        n.oids.push((text_id(&format!("NM_{i:06}"), 1), vec![0]));
    }
    check(
        "f7n",
        &n,
        &[
            ("ndb", 32768, 0x6315ddf59a910882),
            ("nos", 108, 0xe38e089651c04430),
            ("not", 68, 0x91b20d91d89a2711),
            ("ntf", 16384, 0xc29f655d5a3815f8),
            ("nto", 24, 0x1d0b69bb4de52ca4),
        ],
    );
}

#[test]
fn six_hundred_sixty_seven_volumes_split_the_integer_key_trees() {
    // f3v600: -max_file_sz 1kB: 667 volumes db.000..db.666 (3 OIDs each, 2 in the last);
    // volname and volinfo split into branch pages, their leaves interleaved.
    let mut c = Content {
        dbtype: 'p',
        ..Default::default()
    };
    for v in 0..667 {
        c.volumes
            .push((format!("db.{v:03}"), if v == 666 { 2 } else { 3 }));
    }
    c.oids = (0..2000).map(|_| (Vec::new(), vec![0])).collect();
    check(
        "f3v600",
        &c,
        &[
            ("pdb", 49152, 0xfcbd73704a2fe1e4),
            ("pot", 24008, 0x20d29f559fb1c705),
            ("ptf", 16384, 0xc29f655d5a3815f8),
            ("pto", 8004, 0x40f8b7d37fbe21d8),
        ],
    );
}

#[test]
fn many_small_transactions_reuse_freed_pages() {
    // f4p200t7: MAX_LMDB_TXN_ENTRY=7: 400 entries in 58 acc2oid transactions; pages freed by
    // transaction k are taken again from transaction k+2 on.
    check("f4p200t7", &with_max(f4(200, None), "7"), &[PTF_TAXID0]);
    // f4p2000t100: 40 transactions of 100 entries over a two-level tree.
    check(
        "f4p2000t100",
        &with_max(f4(2000, None), "100"),
        &[
            ("pdb", 172032, 0x2d91ad38cec695ff),
            ("pos", 46008, 0x00f982cc16af69e1),
            ("pot", 24008, 0x20d29f559fb1c705),
            ("ptf", 16384, 0xc29f655d5a3815f8),
            ("pto", 8004, 0x40f8b7d37fbe21d8),
        ],
    );
}

#[test]
fn dups_in_sub_pages_and_sub_databases() {
    // f5d2, f5d3, f5d100: one key with 2, 3 and 100 OIDs in a LEAF2 sub-page (built in
    // me_pbuf, with its stale bytes); f5d500: a sub-page of 2016 bytes; f5d501: the 501st dup
    // converts it to a sub-database; f5d3000: a sub-database of two levels.
    let digests: [(&str, Content, [Digest; 5]); 6] = [
        (
            "f5d2",
            f5(2, true),
            [
                ("pdb", 32768, 0xa34aea5de452c6c3),
                ("pos", 448, 0x219d11a25cad6245),
                ("pot", 272, 0xaef304cbdaa11d44),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 92, 0xf82fa0167cd877b2),
            ],
        ),
        (
            "f5d3",
            f5(3, true),
            [
                ("pdb", 32768, 0xb8a67c16d7d75586),
                ("pos", 468, 0xfffb46d82f911ff7),
                ("pot", 284, 0x90d693198943b912),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 96, 0x2412ad6cf2b613e5),
            ],
        ),
        (
            "f5d100",
            f5(100, true),
            [
                ("pdb", 32768, 0x82b872320b4f841a),
                ("pos", 2500, 0xd9abbcc9d5be323c),
                ("pot", 1448, 0x51885a02f9d686c5),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 484, 0x26a175b28489066d),
            ],
        ),
        (
            "f5d500",
            f5(500, false),
            [
                ("pdb", 49152, 0xc3af0366c09ba83b),
                ("pos", 10900, 0x11e40b9405052f5d),
                ("pot", 6008, 0xc92402e70d5655b1),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 2004, 0x49e07fc119691022),
            ],
        ),
        (
            "f5d501",
            f5(501, false),
            [
                ("pdb", 53248, 0x35f7e70140366e91),
                ("pos", 10922, 0x21cc01507b296c4f),
                ("pot", 6020, 0xc31b544181dd95e2),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 2008, 0x7bc92a165d7b3938),
            ],
        ),
        (
            "f5d3000",
            f5(3000, true),
            [
                ("pdb", 139264, 0x3db324e78f4ddec0),
                ("pos", 68301, 0x41cdee639879e1fa),
                ("pot", 36248, 0x7f81d41fce089a39),
                ("ptf", 16384, 0xc29f655d5a3815f8),
                ("pto", 12084, 0x67deca3b737b9d1c),
            ],
        ),
    ];
    for (name, c, d) in &digests {
        check(name, c, d);
    }
    // f5d600t50: MAX_LMDB_TXN_ENTRY=50: the sub-database grows over 26 transactions.
    check("f5d600t50", &with_max(f5(600, true), "50"), &[PTF_TAXID0]);
}

#[test]
fn pdb_chains_share_the_molecule_key() {
    // f5pdbt3: MAX_LMDB_TXN_ENTRY=3: 35 sub-pages of up to 8 chains, 120 transactions.
    check("f5pdbt3", &with_max(f5pdb(), "3"), &[PTF_TAXID0]);
}

#[test]
fn several_tax_ids_per_oid() {
    // t9: OID 0 has 9606 and 10090, OID 2 has 1, 256 and 65536 (keys ordered by their
    // little-endian bytes: 65536 < 256 < 1 < 10090 < 9606).
    check(
        "t9",
        &f9([vec![9606, 10090], vec![9606], vec![1, 256, 65536]]),
        &[],
    );
    // t9neg: -5 and 2147483647.
    check(
        "t9neg",
        &f9([vec![-5, 9606], vec![2_147_483_647], vec![0]]),
        &[
            ("pdb", 32768, 0xc32a7e7e461fe8d9),
            ("pos", 104, 0x0ea2cbbfdd8d242f),
            ("pot", 48, 0x406f905876b9aaaa),
            ("ptf", 16384, 0x63d24d406311734f),
            ("pto", 32, 0x222c83ca5a120676),
        ],
    );
}

#[test]
fn tax_id_map_inserts_at_random_positions() {
    // tm2000d: 2000 distinct tax ids: mid-page splits; tm2000s: 7 tax ids; tm2000dt500: 4
    // transactions that copy most leaves and keep FREE_DBI records.
    let pdb = ("pdb", 143360, 0xdc911154718804dd);
    let pos = ("pos", 46008, 0x00f982cc16af69e1);
    check(
        "tm2000d",
        &f4(2000, Some("distinct")),
        &[
            pdb,
            pos,
            ("pot", 24008, 0x69a2718cdd1ea20c),
            ("ptf", 86016, 0x0f96d8f6323cb904),
            ("pto", 16000, 0x57cab55e94999605),
        ],
    );
    check(
        "tm2000s",
        &f4(2000, Some("seven")),
        &[
            pdb,
            pos,
            ("pot", 24008, 0x6ba5e131244b0036),
            ("ptf", 16384, 0x3e8d884ada3cabcd),
            ("pto", 8028, 0x6195a1f0fc2ab934),
        ],
    );
    check(
        "tm2000dt500",
        &with_max(f4(2000, Some("distinct")), "500"),
        &[
            ("pdb", 172032, 0x84e12dcc333c3829),
            pos,
            ("pot", 24008, 0x69a2718cdd1ea20c),
            ("ptf", 204800, 0x0916c668f9dcb97d),
            ("pto", 16000, 0x57cab55e94999605),
        ],
    );
}

#[test]
fn large_trees_and_overflow_freelist_records() {
    // f4p20000: acc2oid of three levels; tm60000d: 60000 tax ids in 2 transactions (a FREE_DBI
    // record in an overflow page); tm60000dt5ka: 12 transactions with 1- and 2-page overflow
    // records.
    check(
        "f4p20000",
        &f4(20000, None),
        &[
            ("pdb", 1134592, 0xcef0307021d6607d),
            ("pos", 460008, 0x39772ae74f05b8bc),
            ("pot", 240008, 0x8a7bc0204955b545),
            ("ptf", 16384, 0xc29f655d5a3815f8),
            ("pto", 80004, 0x493bfad7fa42b87b),
        ],
    );
    let pdb = ("pdb", 3366912, 0xcc9d838f76e6520d);
    let pos = ("pos", 1380008, 0x234329188ee903fc);
    let pot = ("pot", 720008, 0xa29c5c985be670f3);
    let pto = ("pto", 480000, 0xce39285b78ec1725);
    check(
        "tm60000d",
        &f4(60000, Some("distinct")),
        &[pdb, pos, pot, ("ptf", 3219456, 0x76aa2bb7adcf893c), pto],
    );
    check(
        "tm60000dt5ka",
        &with_max(f4(60000, Some("distinct")), "5000"),
        &[
            ("pdb", 3366912, 0x7f6aba40e0fdf646),
            pos,
            pot,
            ("ptf", 6320128, 0x365209bc3c8fd635),
            pto,
        ],
    );
}

#[test]
fn environment_values_are_parsed_as_ncbi_parses_them() {
    // NStr::StringToInt / StringToUInt / StringToUInt8 with default flags (makeblastdb
    // 2.17.0: MAX_LMDB_TXN_ENTRY="+3" is 3, " 3" and "abc" stop it, "3000000000" overflows
    // int, BLASTDB_LMDB_MAP_SIZE="-5" stops it).
    assert_eq!(string_to_int("+3"), Ok(3));
    assert_eq!(string_to_int("-1"), Ok(-1));
    assert!(string_to_int(" 3").is_err());
    assert!(string_to_int("3 ").is_err());
    assert!(string_to_int("abc").is_err());
    assert!(string_to_int("").is_err());
    assert!(string_to_int("3000000000")
        .unwrap_err()
        .ends_with("to int, overflow"));
    assert_eq!(string_to_uint("4294967295"), Ok(u32::MAX));
    assert!(string_to_uint("4294967296")
        .unwrap_err()
        .ends_with("to unsigned int, overflow"));
    assert!(string_to_uint8("-5").is_err());
    assert_eq!(string_to_uint8("+40960"), Ok(40960));
    assert_eq!(string_to_uint8("18446744073709551615"), Ok(u64::MAX));
    assert!(string_to_uint8("18446744073709551616").is_err());
    // m_MaxEntryPerTxn is unsigned: -1 is 4294967295.
    let s = |v: &str| LmdbEnvSettings {
        max_lmdb_txn_entry: Some(v.to_string()),
        ..Default::default()
    };
    assert_eq!(
        max_entry_per_txn(&LmdbEnvSettings::default()).ok(),
        Some(40_000)
    );
    assert_eq!(max_entry_per_txn(&s("-1")).ok(), Some(u32::MAX));
    assert!(max_entry_per_txn(&s("7x")).is_err());
}

#[test]
fn lookup_file_names_replace_the_last_two_letters() {
    let p = Path::new("/d/db.00.pdb");
    assert_eq!(lmdb_file_name(p, "os"), Path::new("/d/db.00.pos"));
    assert_eq!(
        lmdb_file_name(&lmdb_file_name(p, "tf"), "to"),
        Path::new("/d/db.00.pto")
    );
    assert_eq!(lmdb_file_name(Path::new("x.ndb"), "ot"), Path::new("x.not"));
}

#[test]
fn an_oid_without_ids_in_the_middle_stops_the_pos_file() {
    // f8gimid: a GI-only OID between two OIDs with accessions: makeblastdb aborts with
    // "Input id list not in ascending oid order" (a GI-only OID at the end is fine: f8gitail).
    let mut c = Content {
        dbtype: 'p',
        volumes: vec![("db".to_string(), 3)],
        ..Default::default()
    };
    c.oids.push((text_id("NP_000001", 1), vec![0]));
    c.oids.push((Vec::new(), vec![0]));
    c.oids.push((text_id("NP_000002", 1), vec![0]));
    let dir = TempDir::new("f8gimid");
    let e = write(&c, &dir.0).unwrap_err();
    assert_eq!(
        e.to_string(),
        "Error: Input id list not in ascending oid order"
    );
}

#[test]
fn zero_entries_per_transaction_is_refused() {
    // MAX_LMDB_TXN_ENTRY=0: makeblastdb never leaves its transaction loop.
    let mut c = f4(3, None);
    c.max_txn = Some("0".to_string());
    let dir = TempDir::new("max0");
    assert!(write(&c, &dir.0).is_err());
}

// ---- the gate of the makedb-v5 line (scripts/s1_gate.py) -----------------------------------

/// Parse a manifest of `s1_gate.py`: the content and the oracle files (extension, path).
fn parse_manifest(text: &str) -> (String, Content, Vec<(String, PathBuf)>) {
    let hex = |s: &str| -> Vec<u8> {
        (0..s.len() / 2)
            .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap())
            .collect()
    };
    let mut name = String::new();
    let mut c = Content::default();
    let mut oracle = Vec::new();
    for line in text.lines() {
        let (tag, rest) = line.split_once(' ').unwrap_or((line, ""));
        match tag {
            "CASE" => name = rest.to_string(),
            "DBTYPE" => c.dbtype = rest.chars().next().unwrap(),
            "MAPSIZE" => c.map_size = rest.parse().unwrap(),
            "MAXTXN" => c.max_txn = (rest != "-").then(|| rest.to_string()),
            "ORACLE" => {
                let (ext, path) = rest.split_once(' ').unwrap();
                oracle.push((ext.to_string(), PathBuf::from(path)));
            }
            "VOL" => {
                let (n, k) = rest.split_once(' ').unwrap();
                c.volumes
                    .push((String::from_utf8(hex(n)).unwrap(), k.parse().unwrap()));
            }
            "OID" => c.oids.push((Vec::new(), Vec::new())),
            "E" => {
                let (s, id) = rest.split_once(' ').unwrap();
                c.oids.last_mut().unwrap().0.push((hex(id), s == "1"));
            }
            "T" => {
                c.oids.last_mut().unwrap().1 = rest
                    .split_whitespace()
                    .map(|t| t.parse().unwrap())
                    .collect();
            }
            _ => {}
        }
    }
    (name, c, oracle)
}

/// The gate: `LOSAT_S1_GATE_MANIFEST` names the manifests of `s1_gate.py`; one line per
/// fixture and file goes to `LOSAT_S1_GATE_RESULT` (and to the standard output).
#[test]
#[ignore]
fn s1_gate() {
    let Some(mdir) = std::env::var_os("LOSAT_S1_GATE_MANIFEST") else {
        return;
    };
    let mut files: Vec<PathBuf> = std::fs::read_dir(&mdir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "txt"))
        .collect();
    files.sort();
    let mut lines = Vec::new();
    for f in files {
        let (name, c, oracle) = parse_manifest(&std::fs::read_to_string(&f).unwrap());
        let dir = TempDir::new("gate");
        let start = std::time::Instant::now();
        if let Err(e) = write(&c, &dir.0) {
            lines.push(format!("{name} - error: {e}"));
            continue;
        }
        let secs = start.elapsed().as_secs_f64();
        for (ext, path) in oracle {
            let expected = std::fs::read(&path).unwrap();
            let line = match std::fs::read(dir.0.join(format!("db.{ext}"))) {
                Err(_) => format!("{name} .{ext} missing"),
                Ok(ours) => match first_difference(&ours, &expected) {
                    None => format!("{name} .{ext} equal"),
                    Some(off) => format!(
                        "{name} .{ext} differs at {off} (page {}, offset {}; {} vs {} bytes)",
                        off / 4096,
                        off % 4096,
                        ours.len(),
                        expected.len()
                    ),
                },
            };
            println!("{line}");
            lines.push(line);
        }
        println!("{name}: {secs:.3} s");
    }
    if let Some(r) = std::env::var_os("LOSAT_S1_GATE_RESULT") {
        std::fs::write(r, lines.join("\n") + "\n").unwrap();
    }
    let bad = lines.iter().filter(|l| !l.ends_with(" equal")).count();
    assert_eq!(bad, 0, "{bad} files differ");
}
