//! Tests of the LMDB port: the failures makeblastdb reaches and some invariants of the page
//! code. The byte-for-byte tests against makeblastdb's files are in
//! `LOSAT/src/writedb/lmdb_tests.rs`.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};

use super::page::*;
use super::*;

/// A file path of this test run, removed (with its directory) when the value is dropped.
struct TempFile(PathBuf);

impl TempFile {
    fn new(label: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-mdb-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir.join("db.pdb"))
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(self.0.parent().unwrap());
    }
}

#[test]
fn a_new_environment_is_two_meta_pages() {
    let f = TempFile::new("meta");
    let env = Env::create(&f.0, 700_000_000, 4096).unwrap();
    env.close().unwrap();
    let b = std::fs::read(&f.0).unwrap();
    assert_eq!(b.len(), 8192);
    for (i, p) in b.chunks(4096).enumerate() {
        assert_eq!(pgno(p), i as u64);
        assert_eq!(flags(p), P_META);
        assert_eq!(get_u32(p, 16), MDB_MAGIC);
        assert_eq!(get_u64(p, 32), 700_000_000);
        // FREE_DBI: md_pad = psize, md_flags = MDB_NOSUBDIR | MDB_INTEGERKEY.
        assert_eq!(get_u32(p, 40), 4096);
        assert_eq!(get_u16(p, 44), 0x4008);
        assert_eq!(get_u64(p, 40 + 40), P_INVALID);
        assert_eq!(get_u64(p, 88 + 40), P_INVALID);
        assert_eq!(get_u64(p, 136), 1);
        assert_eq!(get_u64(p, 144), 0);
    }
}

// makeblastdb stops with "mdb_put: MDB_BAD_VALSIZE: Unsupported size of key/DB name/data, or
// wrong DUPFIXED size" on a key of 512 bytes or more (survey A, f6d509-f6d600).
#[test]
fn keys_over_511_bytes_fail_as_in_makeblastdb() {
    let f = TempFile::new("valsize");
    let mut env = Env::create(&f.0, 700_000_000, 4096).unwrap();
    env.txn_begin().unwrap();
    let dbi = env
        .dbi_open("acc2oid", MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED)
        .unwrap();
    assert!(env
        .put(dbi, &[b'k'; 511], &0u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap());
    let e = env
        .put(dbi, &[b'k'; 512], &0u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap_err();
    assert_eq!(e.code(), Some(MDB_BAD_VALSIZE));
    assert_eq!(
        e.to_string(),
        "mdb_put: MDB_BAD_VALSIZE: Unsupported size of key/DB name/data, or wrong DUPFIXED size"
    );
    assert!(env
        .put(dbi, b"", &0u32.to_le_bytes(), MDB_APPENDDUP)
        .is_err());
}

// `mdb_page_alloc` fails with MDB_MAP_FULL when `pgno + num >= me_maxpg` (mdb.c:2290-2294).
#[test]
fn a_full_map_fails_where_lmdb_fails() {
    let f = TempFile::new("mapfull");
    // 5 pages: the metas (0, 1), the main leaf (2), volinfo's leaf (3); page 4 is refused.
    let mut env = Env::create(&f.0, 5 * 4096, 4096).unwrap();
    env.txn_begin().unwrap();
    let volinfo = env
        .dbi_open("volinfo", MDB_CREATE | MDB_INTEGERKEY)
        .unwrap();
    let volname = env
        .dbi_open("volname", MDB_CREATE | MDB_INTEGERKEY)
        .unwrap();
    assert!(env
        .put(volinfo, &0u32.to_le_bytes(), &1i32.to_le_bytes(), 0)
        .unwrap());
    let e = env.put(volname, &0u32.to_le_bytes(), b"db", 0).unwrap_err();
    assert_eq!(e.code(), Some(MDB_MAP_FULL));
    env.txn_abort();
    // A larger map (`mdb_env_set_mapsize` between transactions) lets the same puts through.
    env.set_mapsize(6 * 4096).unwrap();
    env.txn_begin().unwrap();
    let volinfo = env
        .dbi_open("volinfo", MDB_CREATE | MDB_INTEGERKEY)
        .unwrap();
    let volname = env
        .dbi_open("volname", MDB_CREATE | MDB_INTEGERKEY)
        .unwrap();
    assert!(env
        .put(volinfo, &0u32.to_le_bytes(), &1i32.to_le_bytes(), 0)
        .unwrap());
    assert!(env.put(volname, &0u32.to_le_bytes(), b"db", 0).unwrap());
    // A put of an existing key of a database without dups replaces its data.
    assert!(env.put(volname, &0u32.to_le_bytes(), b"dx", 0).unwrap());
    env.txn_commit().unwrap();
    env.close().unwrap();
    assert_eq!(std::fs::metadata(&f.0).unwrap().len(), 5 * 4096);
}

// A dup that exists already: MDB_APPENDDUP returns MDB_KEYEXIST (mdb.c:6745-6746), which
// lmdb::dbi_put hands back as `false` (and CWriteDB_LMDB turns into "acc2oid error for id").
#[test]
fn an_existing_dup_is_keyexist() {
    let f = TempFile::new("keyexist");
    let mut env = Env::create(&f.0, 700_000_000, 4096).unwrap();
    env.txn_begin().unwrap();
    let dbi = env
        .dbi_open("acc2oid", MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED)
        .unwrap();
    assert!(env
        .put(dbi, b"XP_1", &5u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap());
    assert!(!env
        .put(dbi, b"XP_1", &5u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap());
    assert!(env
        .put(dbi, b"XP_1", &6u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap());
    assert!(!env
        .put(dbi, b"XP_1", &6u32.to_le_bytes(), MDB_APPENDDUP)
        .unwrap());
}

#[test]
fn integer_keys_order_as_unsigned_numbers() {
    let c = Cmp::Int;
    assert!(c.cmp(&1u32.to_le_bytes(), &256u32.to_le_bytes()) < 0);
    assert!(c.cmp(&65_536u32.to_le_bytes(), &65_535u32.to_le_bytes()) > 0);
    assert_eq!(c.cmp(&7u64.to_le_bytes(), &7u64.to_le_bytes()), 0);
    assert!(c.cmp(&(u64::MAX).to_le_bytes(), &1u64.to_le_bytes()) > 0);
    let m = Cmp::Memn;
    assert!(m.cmp(b"AB", b"ABC") < 0);
    assert!(m.cmp(b"AC", b"ABC") > 0);
    assert!(m.cmp(&[0xFF], b"A") > 0);
}

// mdb_page_copy copies [0, align8(lower)) and [floor8(upper), psize): the destination keeps
// its bytes in between (the stale gap of a copied page).
#[test]
fn page_copy_keeps_the_middle_of_the_destination() {
    let psize = 4096;
    let mut src = vec![0xAAu8; psize];
    set_flags(&mut src, P_LEAF);
    set_lower(&mut src, 22);
    set_upper(&mut src, 4001);
    let mut dst = vec![0x55u8; psize];
    Env::page_copy(&mut dst, &src, psize);
    assert!(dst[..24].iter().skip(16).all(|&b| b == 0xAA));
    assert!(dst[24..4000].iter().all(|&b| b == 0x55));
    assert!(dst[4000..].iter().all(|&b| b == 0xAA));
    // A LEAF2 page (or a gap under 8 bytes) is copied up to `psize - (gap & ~7)`.
    set_flags(&mut src, P_LEAF | P_LEAF2);
    let mut dst = vec![0x55u8; psize];
    Env::page_copy(&mut dst, &src, psize);
    let unused = (4001 - 22) & !7;
    assert!(dst[..psize - unused].iter().skip(16).all(|&b| b == 0xAA));
    assert!(dst[psize - unused..].iter().all(|&b| b == 0x55));
}

#[test]
fn page_size_is_the_os_page_size_up_to_32768() {
    let p = os_page_size();
    assert!(p.is_power_of_two() && (4096..=MAX_PAGESIZE).contains(&p));
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    assert_eq!(p, 4096);
}
