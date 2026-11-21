use super::*;
use std::fs::File;

/// Fixture databases (made by NCBI makeblastdb 2.17.0, see `tests/seqdb/README.md`).
const FX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb");
/// Values written by `dump.py` (a Python reference reader) for the fixtures.
const EXP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb/expected");

/// Path of an `expected/*.txt` file: the repository fixtures, else the directory named by
/// `LOSAT_TEST_18S_EXPECTED` (the dumps of the NCBI-distributed 18S database are too big to commit).
fn exp_path(name: &str) -> String {
    let repo = format!("{EXP}/{name}.txt");
    if std::path::Path::new(&repo).exists() {
        return repo;
    }
    let dir = std::env::var("LOSAT_TEST_18S_EXPECTED").unwrap_or_default();
    format!("{dir}/{name}.txt")
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

/// One database of an `expected/*.txt` file written by `gen/expect.py` from `dump.py`.
struct ExpDb {
    name: String,
    flags: u16,
    depth: u16,
    entries: u64,
    keys: Vec<(Vec<u8>, Vec<Vec<u8>>)>,
}

fn load_expected(name: &str) -> Vec<ExpDb> {
    let text = std::fs::read_to_string(exp_path(name)).unwrap();
    let mut dbs: Vec<ExpDb> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "META" => {}
            "DB" => dbs.push(ExpDb {
                name: f[1].to_string(),
                flags: f[2].parse().unwrap(),
                depth: f[3].parse().unwrap(),
                entries: f[4].parse().unwrap(),
                keys: Vec::new(),
            }),
            "K" => {
                let vals = f[2].split(',').map(unhex).collect();
                dbs.last_mut().unwrap().keys.push((unhex(f[1]), vals));
            }
            _ => panic!("bad expected line {line}"),
        }
    }
    dbs
}

/// Every named database: record fields, the in-order walk, `first_value`, `values`.
fn check_against_dump<R: ReadAt + ?Sized>(env: &Env<R>, exp_name: &str) {
    for ex in load_expected(exp_name) {
        let db = env
            .open_db(ex.name.as_bytes())
            .unwrap()
            .unwrap_or_else(|| panic!("{exp_name}: {} missing", ex.name));
        assert_eq!(
            (db.flags, db.depth, db.entries),
            (ex.flags, ex.depth, ex.entries),
            "{exp_name} {}",
            ex.name
        );
        let mut walked: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        env.for_each(&db, |k, v| {
            walked.push((k.to_vec(), v.to_vec()));
            Ok(())
        })
        .unwrap();
        let want: Vec<(Vec<u8>, Vec<u8>)> = ex
            .keys
            .iter()
            .flat_map(|(k, vs)| vs.iter().map(move |v| (k.clone(), v.clone())))
            .collect();
        assert_eq!(
            walked.len(),
            want.len(),
            "{exp_name} {} walk length",
            ex.name
        );
        assert!(walked == want, "{exp_name} {} walk differs", ex.name);
        assert_eq!(walked.len() as u64, db.entries);
        for (k, vs) in &ex.keys {
            assert_eq!(
                env.first_value(&db, k).unwrap().as_ref(),
                Some(&vs[0]),
                "{exp_name} first of {}",
                hex(k)
            );
            assert_eq!(
                &env.values(&db, k).unwrap(),
                vs,
                "{exp_name} values of {}",
                hex(k)
            );
        }
    }
}

fn check_file(name: &str) {
    let f = File::open(format!("{FX}/{name}")).unwrap();
    let env = Env::open(&f).unwrap();
    check_against_dump(&env, name);
}

#[test]
fn fixtures_match_dump_py() {
    for n in [
        "p5.pdb",
        "p5.ptf",
        "n5.ndb",
        "n5.ntf",
        "mv.pdb",
        "mv.ptf",
        "dup.pdb",
        "dup.ptf",
        "thr500.pdb",
        "thr501.pdb",
        "mvbig.pdb",
    ] {
        check_file(n);
    }
}

#[test]
fn in_memory_sources_agree_with_files() {
    let bytes = std::fs::read(format!("{FX}/p5.pdb")).unwrap();
    let env = Env::open(&bytes).unwrap();
    check_against_dump(&env, "p5.pdb");
    let env2 = Env::open(&bytes[..]).unwrap();
    check_against_dump(&env2, "p5.pdb");
    assert_eq!((env.page_size(), env.txnid()), (4096, 2));
}

#[test]
fn p5_literals() {
    let bytes = std::fs::read(format!("{FX}/p5.pdb")).unwrap();
    let env = Env::open(&bytes).unwrap();
    let acc = env.open_db(b"acc2oid").unwrap().unwrap();
    assert_eq!(
        (acc.flags, acc.depth, acc.entries, acc.root),
        (0x14, 1, 16, 6)
    );
    // a key with two OIDs lives in an inline sub-page; both values come back
    assert_eq!(
        env.values(&acc, b"AAA00001").unwrap(),
        vec![2u32.to_le_bytes().to_vec(), 3u32.to_le_bytes().to_vec()]
    );
    assert_eq!(
        env.first_value(&acc, b"AAA00001").unwrap(),
        Some(2u32.to_le_bytes().to_vec())
    );
    assert_eq!(
        env.first_value(&acc, b"AAA00001.2").unwrap(),
        Some(3u32.to_le_bytes().to_vec())
    );
    // keys are case sensitive; absent keys are None / empty
    assert_eq!(env.first_value(&acc, b"aaa00001").unwrap(), None);
    assert_eq!(
        env.first_value(&acc, b"MYPROT9").unwrap(),
        Some(4u32.to_le_bytes().to_vec())
    );
    assert_eq!(
        env.first_value(&acc, b"myprot9").unwrap(),
        Some(4u32.to_le_bytes().to_vec())
    );
    assert!(env.values(&acc, b"nokey").unwrap().is_empty());
    assert_eq!(env.first_value(&acc, b"AAA0000").unwrap(), None);
    assert_eq!(env.first_value(&acc, b"zzzz").unwrap(), None);
    // an empty key is MDB_BAD_VALSIZE in NCBI (mdb_cursor_set)
    assert!(env
        .first_value(&acc, b"")
        .unwrap_err()
        .contains("MDB_BAD_VALSIZE"));
    // unknown database name
    assert_eq!(env.open_db(b"taxid2offset").unwrap(), None);
    // integer keys: volume 0
    let vn = env.open_db(b"volname").unwrap().unwrap();
    assert_eq!(vn.flags, 0x08);
    assert_eq!(
        env.first_value(&vn, &0u32.to_le_bytes()).unwrap(),
        Some(b"p5".to_vec())
    );
    let vi = env.open_db(b"volinfo").unwrap().unwrap();
    assert_eq!(
        env.first_value(&vi, &0u32.to_le_bytes()).unwrap(),
        Some(7i32.to_le_bytes().to_vec())
    );
    assert_eq!(env.first_value(&vi, &1u32.to_le_bytes()).unwrap(), None);
}

#[test]
fn p5_taxid2offset_order_is_memcmp_of_le_bytes() {
    // p5.ptf: txnid 0 (page 0, empty) and 1 (page 1); page 1 is picked.
    let f = File::open(format!("{FX}/p5.ptf")).unwrap();
    let env = Env::open(&f).unwrap();
    assert_eq!(env.txnid(), 1);
    let db = env.open_db(b"taxid2offset").unwrap().unwrap();
    let mut keys = Vec::new();
    env.for_each(&db, |k, v| {
        keys.push((
            i32::from_le_bytes(k.try_into().unwrap()),
            u64::from_le_bytes(v.try_into().unwrap()),
        ));
        Ok(())
    })
    .unwrap();
    assert_eq!(
        keys.iter().map(|k| k.0).collect::<Vec<_>>(),
        vec![0, 7227, 10090, 9606]
    );
}

#[test]
fn dup_subdatabase_keys_work_unlike_ncbi() {
    let f = File::open(format!("{FX}/dup.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    let acc = env.open_db(b"acc2oid").unwrap().unwrap();
    // 700 OIDs under one key (sub-database), in memcmp order of the LE bytes: 0, 256, 512, 1, ...
    let vals = env.values(&acc, b"DUP00001").unwrap();
    assert_eq!(vals.len(), 700);
    let oids: Vec<u32> = vals
        .iter()
        .map(|v| u32::from_le_bytes(v[..].try_into().unwrap()))
        .collect();
    assert_eq!(&oids[..7], &[0, 256, 512, 1, 257, 513, 2]);
    assert_eq!(
        env.first_value(&acc, b"DUP00001").unwrap(),
        Some(0u32.to_le_bytes().to_vec())
    );
    let mut sorted = oids.clone();
    sorted.sort();
    assert_eq!(sorted, (0..700).collect::<Vec<u32>>());
    assert_eq!(
        env.first_value(&acc, b"DUP00001.5").unwrap(),
        Some(4u32.to_le_bytes().to_vec())
    );
}

#[test]
fn threshold_500_subpage_501_subdb() {
    for (name, n) in [("thr500.pdb", 500usize), ("thr501.pdb", 501)] {
        let f = File::open(format!("{FX}/{name}")).unwrap();
        let env = Env::open(&f).unwrap();
        let acc = env.open_db(b"acc2oid").unwrap().unwrap();
        let vals = env.values(&acc, b"DUP00001").unwrap();
        assert_eq!(vals.len(), n, "{name}");
        let mut seen = std::collections::BTreeSet::new();
        for v in &vals {
            seen.insert(u32::from_le_bytes(v[..].try_into().unwrap()));
        }
        assert_eq!(seen.len(), n);
    }
}

/// `first_key_at_least` (`MDB_SET_RANGE`) against the dump of a depth-two `acc2oid` (2,400
/// keys over several leaves): each key finds itself, a key just after it finds the next key
/// (the first of the next leaf after a leaf's last), and past the last key nothing.
#[test]
fn first_key_at_least_follows_the_key_order_across_leaves() {
    let f = File::open(format!("{FX}/mvbig.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    let acc = env.open_db(b"acc2oid").unwrap().unwrap();
    assert_eq!(acc.depth, 2);
    let exp = load_expected("mvbig.pdb");
    let keys: Vec<&Vec<u8>> = exp
        .iter()
        .find(|db| db.name == "acc2oid")
        .unwrap()
        .keys
        .iter()
        .map(|(key, _)| key)
        .collect();
    assert_eq!(keys.len(), 2400);
    assert_eq!(
        env.first_key_at_least(&acc, b"\x01").unwrap().as_ref(),
        Some(keys[0])
    );
    for (at, key) in keys.iter().enumerate() {
        assert_eq!(
            env.first_key_at_least(&acc, key).unwrap().as_ref(),
            Some(*key)
        );
        let after = [key.as_slice(), b"\0"].concat();
        assert_eq!(
            env.first_key_at_least(&acc, &after).unwrap().as_ref(),
            keys.get(at + 1).copied(),
            "after {}",
            String::from_utf8_lossy(key)
        );
        assert!(env.has_key_starting_with(&acc, &key[..3]).unwrap());
    }
    assert!(!env.has_key_starting_with(&acc, b"BIG00000.10").unwrap());
    assert!(!env.has_key_starting_with(&acc, b"big").unwrap());
    assert_eq!(env.first_key_at_least(&acc, b"\xff").unwrap(), None);
    assert!(env.first_key_at_least(&acc, b"").is_err());
}

#[test]
fn mvbig_integer_keys_depth_two() {
    let f = File::open(format!("{FX}/mvbig.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    let vn = env.open_db(b"volname").unwrap().unwrap();
    assert_eq!((vn.depth, vn.entries), (2, 400));
    // numeric order (memcmp order would put 256 before 1)
    let mut keys = Vec::new();
    env.for_each(&vn, |k, _| {
        keys.push(u32::from_le_bytes(k.try_into().unwrap()));
        Ok(())
    })
    .unwrap();
    assert_eq!(keys, (0..400).collect::<Vec<u32>>());
    for i in [0u32, 1, 255, 256, 257, 399] {
        let name = env.first_value(&vn, &i.to_le_bytes()).unwrap().unwrap();
        assert_eq!(name, format!("mvbig.{i:03}").into_bytes(), "volume {i}");
    }
    assert_eq!(env.first_value(&vn, &400u32.to_le_bytes()).unwrap(), None);
}

/// Needs the NCBI-distributed `18S_fungal_sequences` database: `LOSAT_TEST_18S_DB` is its path
/// without the extension, `LOSAT_TEST_18S_EXPECTED` the directory of its `dump.py` outputs
/// (`18S_fungal_sequences.ndb.txt`, `.ntf.txt`).
#[test]
#[ignore = "needs LOSAT_TEST_18S_DB and LOSAT_TEST_18S_EXPECTED"]
fn ncbi_distributed_18s_matches_dump_py() {
    let prefix = std::env::var("LOSAT_TEST_18S_DB").expect("LOSAT_TEST_18S_DB");
    let ndb = format!("{prefix}.ndb");
    let ntf = format!("{prefix}.ntf");
    let f = File::open(&ndb).unwrap();
    let env = Env::open(&f).unwrap();
    check_against_dump(&env, "18S_fungal_sequences.ndb");
    let f = File::open(&ntf).unwrap();
    let env = Env::open(&f).unwrap();
    check_against_dump(&env, "18S_fungal_sequences.ntf");
}

// ---- synthetic files: overflow values, meta choice, errors ----

fn put16(b: &mut [u8], o: usize, v: u16) {
    b[o..o + 2].copy_from_slice(&v.to_le_bytes());
}

fn meta_page(psize: usize, pgno: u64, txnid: u64, main_root: u64, last_pg: u64) -> Vec<u8> {
    let mut p = vec![0u8; psize];
    p[0..8].copy_from_slice(&pgno.to_le_bytes());
    put16(&mut p, 10, P_META);
    p[16..20].copy_from_slice(&MDB_MAGIC.to_le_bytes());
    p[20..24].copy_from_slice(&1u32.to_le_bytes());
    // FREE_DBI: md_pad = page size, flags 0x4008, empty
    p[40..44].copy_from_slice(&(psize as u32).to_le_bytes());
    put16(&mut p, 44, 0x4008);
    p[80..88].copy_from_slice(&u64::MAX.to_le_bytes());
    // MAIN_DBI
    put16(&mut p, 88 + 6, 1);
    p[88 + 40..88 + 48].copy_from_slice(&main_root.to_le_bytes());
    p[136..144].copy_from_slice(&last_pg.to_le_bytes());
    p[144..152].copy_from_slice(&txnid.to_le_bytes());
    p
}

/// Leaf page with `(flags, key, data, dsize)` nodes (keys given in order).
fn leaf_page(psize: usize, pgno: u64, nodes: &[(u16, &[u8], Vec<u8>, usize)]) -> Vec<u8> {
    let mut p = vec![0u8; psize];
    p[0..8].copy_from_slice(&pgno.to_le_bytes());
    put16(&mut p, 10, P_LEAF);
    let mut upper = psize;
    for (i, (fl, k, d, dsize)) in nodes.iter().enumerate() {
        let sz = (8 + k.len() + d.len() + 1) & !1;
        upper -= sz;
        put16(&mut p, 16 + 2 * i, upper as u16);
        put16(&mut p, upper, (*dsize & 0xffff) as u16);
        put16(&mut p, upper + 2, (*dsize >> 16) as u16);
        put16(&mut p, upper + 4, *fl);
        put16(&mut p, upper + 6, k.len() as u16);
        p[upper + 8..upper + 8 + k.len()].copy_from_slice(k);
        p[upper + 8 + k.len()..upper + 8 + k.len() + d.len()].copy_from_slice(d);
    }
    put16(&mut p, 12, (16 + 2 * nodes.len()) as u16);
    put16(&mut p, 14, upper as u16);
    p
}

/// Main database with a plain value, an overflow value (2 pages) and a named-DB record.
fn synthetic(txn0: u64, txn1: u64) -> (Vec<u8>, Vec<u8>) {
    let psize = 4096usize;
    let big: Vec<u8> = (0..5000u32).map(|i| (i % 251) as u8).collect();
    let mut dbrec = vec![0u8; 48];
    dbrec[4..6].copy_from_slice(&0x14u16.to_le_bytes());
    dbrec[40..48].copy_from_slice(&u64::MAX.to_le_bytes());
    let nodes: Vec<(u16, &[u8], Vec<u8>, usize)> = vec![
        (F_BIGDATA, b"big", 3u64.to_le_bytes().to_vec(), big.len()),
        (F_SUBDATA, b"emptydb", dbrec, 48),
        (0, b"small", b"xyz".to_vec(), 3),
    ];
    let mut file = Vec::new();
    file.extend(meta_page(psize, 0, txn0, 2, 4));
    file.extend(meta_page(psize, 1, txn1, 2, 4));
    file.extend(leaf_page(psize, 2, &nodes));
    let mut ov = vec![0u8; 2 * psize];
    ov[0..8].copy_from_slice(&3u64.to_le_bytes());
    put16(&mut ov, 10, P_OVERFLOW);
    ov[12..16].copy_from_slice(&2u32.to_le_bytes());
    ov[16..16 + big.len()].copy_from_slice(&big);
    file.extend(ov);
    (file, big)
}

#[test]
fn overflow_values_and_named_records() {
    let (file, big) = synthetic(1, 2);
    let env = Env::open(&file).unwrap();
    let main = Db {
        pad: 0,
        flags: 0,
        depth: 1,
        branch_pages: 0,
        leaf_pages: 1,
        overflow_pages: 0,
        entries: 3,
        root: 2,
    };
    assert_eq!(env.first_value(&main, b"big").unwrap(), Some(big.clone()));
    assert_eq!(env.values(&main, b"big").unwrap(), vec![big]);
    assert_eq!(
        env.first_value(&main, b"small").unwrap(),
        Some(b"xyz".to_vec())
    );
    let db = env.open_db(b"emptydb").unwrap().unwrap();
    assert_eq!(db.flags, 0x14);
    assert!(db.is_empty_tree());
    assert_eq!(env.first_value(&db, b"k").unwrap(), None);
    env.for_each(&db, |_, _| Err("must not be called".into()))
        .unwrap();
    // not a database record
    assert!(env
        .open_db(b"small")
        .unwrap_err()
        .contains("MDB_INCOMPATIBLE"));
    // corrupt overflow: claims more bytes than its pages hold
    let mut bad = file.clone();
    bad[3 * 4096 + 12..3 * 4096 + 16].copy_from_slice(&1u32.to_le_bytes());
    let env = Env::open(&bad).unwrap();
    assert!(env.first_value(&main, b"big").is_err());
}

#[test]
fn meta_pick_larger_txnid_tie_meta0() {
    // meta 0 has the larger txnid: still root 2 either way; tell the metas apart by txnid()
    let (f, _) = synthetic(5, 3);
    assert_eq!(Env::open(&f).unwrap().txnid(), 5);
    let (f, _) = synthetic(3, 5);
    assert_eq!(Env::open(&f).unwrap().txnid(), 5);
    let (f, _) = synthetic(4, 4);
    assert_eq!(Env::open(&f).unwrap().txnid(), 4);
}

#[test]
fn open_errors() {
    let (good, _) = synthetic(1, 2);
    assert!(Env::open(&Vec::<u8>::new())
        .unwrap_err()
        .contains("too small"));
    assert!(Env::open(&good[..100]).unwrap_err().contains("too small"));
    let mut bad = good.clone();
    bad[16] ^= 1;
    assert!(Env::open(&bad).unwrap_err().contains("magic"));
    let mut bad = good.clone();
    bad[20] = 2;
    assert!(Env::open(&bad)
        .unwrap_err()
        .contains("MDB_VERSION_MISMATCH"));
    let mut bad = good.clone();
    bad[10] = 0;
    assert!(Env::open(&bad).unwrap_err().contains("not a meta page"));
    let mut bad = good.clone();
    bad[40..44].copy_from_slice(&3000u32.to_le_bytes());
    assert!(Env::open(&bad).unwrap_err().contains("page size"));
    // second meta page cut off
    assert!(Env::open(&good[..4096 + 100]).is_err());
    // a text file
    assert!(Env::open(&b"hello world, this is not an LMDB file at all ".repeat(10)).is_err());
}

/// Cheap deterministic generator for the corruption tests.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Open, list and look up in a possibly corrupt file; any outcome but a panic or a hang is fine.
fn exercise(bytes: &[u8]) {
    let Ok(env) = Env::open(bytes) else { return };
    for name in [&b"acc2oid"[..], b"volname", b"volinfo", b"taxid2offset"] {
        let Ok(Some(db)) = env.open_db(name) else {
            continue;
        };
        let _ = env.for_each(&db, |_, _| Ok(()));
        let _ = env.first_value(&db, b"DUP00001");
        let _ = env.values(&db, b"AAA00001");
        let _ = env.values(&db, &7u32.to_le_bytes());
        let _ = env.first_value(&db, &9606i32.to_le_bytes());
    }
}

#[test]
fn corrupt_files_give_errors_not_panics() {
    for name in [
        "p5.pdb",
        "n5.ndb",
        "mv.pdb",
        "dup.pdb",
        "p5.ptf",
        "mvbig.pdb",
    ] {
        let orig = std::fs::read(format!("{FX}/{name}")).unwrap();
        let mut rng = Lcg(0x1234_5678 ^ orig.len() as u64);
        for round in 0..400 {
            let mut b = orig.clone();
            for _ in 0..1 + round % 4 {
                let at = (rng.next() as usize) % b.len();
                // aim half of the flips at the first bytes of pages (headers)
                let at = if rng.next().is_multiple_of(2) {
                    at & !4095 | (rng.next() as usize % 160)
                } else {
                    at
                };
                let at = at.min(b.len() - 1);
                b[at] = rng.next() as u8;
            }
            exercise(&b);
        }
        // truncation at many lengths
        for cut in (0..orig.len()).step_by(517) {
            exercise(&orig[..cut]);
        }
    }
}

#[test]
fn cyclic_tree_is_an_error_not_a_hang() {
    // mvbig volname tree has a branch root: point every child of a branch page at itself
    let mut b = std::fs::read(format!("{FX}/mvbig.pdb")).unwrap();
    let env = Env::open(&b).unwrap();
    let vn = env.open_db(b"volname").unwrap().unwrap();
    let root = vn.root as usize;
    let page = &mut b[root * 4096..(root + 1) * 4096];
    let n = (u16::from_le_bytes([page[12], page[13]]) as usize - 16) / 2;
    for i in 0..n {
        let o = u16::from_le_bytes([page[16 + 2 * i], page[17 + 2 * i]]) as usize;
        page[o..o + 2].copy_from_slice(&(root as u16).to_le_bytes());
        page[o + 2..o + 6].copy_from_slice(&[0, 0, 0, 0]);
    }
    let env = Env::open(&b).unwrap();
    assert!(env.for_each(&vn, |_, _| Ok(())).is_err());
    assert!(env.first_value(&vn, &300u32.to_le_bytes()).is_err());
}
