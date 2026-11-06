//! Damaged copies of the files that map ids and taxids to OIDs, whose counts drove a reader to
//! allocate or loop far past what the file holds: the volume tables of the LMDB file of version
//! 5 (`tax_files::volumes`), the string ISAM index of version 4 (`StringIsam::load_data`) and the
//! taxonomy SQLite file (`TaxonomySqlite::open_named`). Each case was a process abort (a failed
//! allocation) or a loop that did not end before the fix; now each is an error.

use std::path::{Path, PathBuf};

use super::isam::{NegativeOids, SiOid, StringIsam};
use super::lmdb_read::Env;
use super::tax_files;
use super::tax_sqlite::TaxonomySqlite;

fn lmdb_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/lmdb")
}

fn ids_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids")
}

fn taxsqlite_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxsqlite")
}

/// A directory of its own for each test (memory `parallel-test-shared-file-race`), removed at
/// the end.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-damaged-ids-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Sets `md_entries` (offset 32 of the 48-byte `MDB_db` record) of the named database `name`
/// at every match of its node (`flags, key size, key, record`) in `file`, live or stale.
fn set_entries(file: &mut [u8], name: &[u8], entries: u64) -> usize {
    let mut node = vec![2u8, 0];
    node.extend_from_slice(&(name.len() as u16).to_le_bytes());
    node.extend_from_slice(name);
    let mut found = 0;
    let mut at = 0;
    while let Some(i) = file[at..]
        .windows(node.len())
        .position(|w| w == node.as_slice())
    {
        let record = at + i + node.len();
        file[record + 32..record + 40].copy_from_slice(&entries.to_le_bytes());
        found += 1;
        at = record;
    }
    found
}

/// `md_entries` of both volume tables set to a count that the file cannot hold. Before the fix the
/// count of 2^28 ran 84 s up to 4.2 GB, then aborted on an allocation of 8 GiB (under an 8 GB
/// address-space limit); 2^32 - 1 entries would ask for about 128 GiB.
#[test]
fn volume_tables_with_more_entries_than_the_file_holds_are_errors() {
    for name in ["p5.pdb", "mv.pdb"] {
        let mut file = std::fs::read(lmdb_dir().join(name)).unwrap();
        let intact = Env::open(&file).unwrap();
        let volumes = tax_files::volumes(&intact).unwrap();
        assert!(!volumes.is_empty());
        for entries in [1u64 << 28, u64::from(u32::MAX)] {
            assert!(set_entries(&mut file, b"volname", entries) > 0, "{name}");
            assert!(set_entries(&mut file, b"volinfo", entries) > 0, "{name}");
            let env = Env::open(&file).unwrap();
            let error = tax_files::volumes(&env).unwrap_err();
            assert_eq!(
                error,
                format!("Volinfo error : {entries} volumes"),
                "{name}"
            );
            assert!(!tax_files::is_ncbi_volume_error(&error), "{name}: {error}");
        }
    }
}

/// `num_terms` and `page_size` of a string ISAM index rewritten to numbers that pass the header
/// checks (`(samples - 1) * page < terms <= samples * page`) but that the data file cannot hold.
/// Before the fix the first page read asked for 11.4 GB at once (an abort under an 8 GB
/// address-space limit).
#[test]
fn string_isam_pages_larger_than_the_data_file_are_errors() {
    let dir = TempDir::new("isam");
    let index = std::fs::read(ids_dir().join("prot_v4.psi")).unwrap();
    let samples = u32::from_be_bytes(index[16..20].try_into().unwrap());
    assert_eq!(samples, 3, "the fixture's sample count");
    let page: u32 = 715_000_000;
    let terms: u32 = 2 * page + 1;
    let mut damaged = index.clone();
    damaged[12..16].copy_from_slice(&terms.to_be_bytes());
    damaged[20..24].copy_from_slice(&page.to_be_bytes());
    let (idx, dat) = (dir.0.join("d.psi"), dir.0.join("d.psd"));
    std::fs::write(&idx, &damaged).unwrap();
    std::fs::copy(ids_dir().join("prot_v4.psd"), &dat).unwrap();
    let isam = StringIsam::open(&idx, &dat).unwrap();
    // The batch functions read a whole page of `page_size` records (`load_data`).
    let mut ids = [
        SiOid::new(b"xp_000011.2".to_vec()),
        SiOid::new(b"0008123a".to_vec()),
    ];
    let error = isam.ids_to_oids(0, &mut ids).unwrap_err().to_string();
    assert!(error.contains("string ISAM data page"), "{error}");
    let mut negative = NegativeOids::new();
    let mut list = vec![b"xp_000011.2".to_vec()];
    let error = isam
        .negative_scan(0, 10, &mut list, &mut negative)
        .unwrap_err()
        .to_string();
    assert!(error.contains("string ISAM data page"), "{error}");
    // The intact index gives the OID of the same key.
    let intact = StringIsam::open(
        &ids_dir().join("prot_v4.psi"),
        &ids_dir().join("prot_v4.psd"),
    )
    .unwrap();
    let mut ids = [SiOid::new(b"xp_000011.2".to_vec())];
    intact.ids_to_oids(0, &mut ids).unwrap();
    assert_eq!(ids[0].oid, Some(11));
}

/// A taxonomy SQLite file whose `sqlite_master` tree is a chain of interior pages, each pointing
/// at the next page from every cell and from its right-most pointer (a tree of 58 * 72^5 leaves
/// on 7 pages), and whose header claims 2^32 - 1 pages.
fn branching_sqlite() -> Vec<u8> {
    let mut b = std::fs::read(taxsqlite_dir().join("p512.sqlite3")).unwrap();
    let ps = 512usize;
    assert_eq!(u16::from_be_bytes([b[16], b[17]]) as usize, ps);
    b[28..32].copy_from_slice(&u32::MAX.to_be_bytes());
    let counter: [u8; 4] = b[24..28].try_into().unwrap();
    b[92..96].copy_from_slice(&counter);
    for pg in 1..=7usize {
        let base = (pg - 1) * ps;
        let h = if pg == 1 { 100 } else { 0 };
        b[base + h..base + ps].fill(0);
        let page = &mut b[base..base + ps];
        if pg == 7 {
            page[h] = 0x0d;
            page[h + 5..h + 7].copy_from_slice(&(ps as u16).to_be_bytes());
            continue;
        }
        page[h] = 0x05;
        let k = (ps - h - 12) / 7;
        let mut content = ps;
        for i in 0..k {
            content -= 5;
            page[content..content + 4].copy_from_slice(&((pg + 1) as u32).to_be_bytes());
            page[content + 4] = 1;
            let ptr = h + 12 + 2 * i;
            page[ptr..ptr + 2].copy_from_slice(&(content as u16).to_be_bytes());
        }
        page[h + 3..h + 5].copy_from_slice(&(k as u16).to_be_bytes());
        page[h + 5..h + 7].copy_from_slice(&(content as u16).to_be_bytes());
        page[h + 8..h + 12].copy_from_slice(&((pg + 1) as u32).to_be_bytes());
    }
    b
}

/// Before the fix the walk of `sqlite_master` did not end: the tree is not cyclic, so the depth
/// limit is never reached, and the page budget (the header's page count plus 64) was 2^32 + 63.
/// SQLite refuses the file (`database disk image is malformed`).
#[test]
fn sqlite_header_page_count_past_the_file_is_an_error() {
    let file = branching_sqlite();
    let error = TaxonomySqlite::open(&file).err().expect("must not open");
    assert!(error.contains("malformed"), "{error}");
    // The header's count alone, as SQLite 3.51.1 reads it: 248 pages on a file of 247 pages is
    // malformed; on 247 pages and 100 bytes (248 pages rounded up) the file opens.
    let mut file = std::fs::read(taxsqlite_dir().join("p512.sqlite3")).unwrap();
    assert_eq!(file.len(), 247 * 512);
    let counter: [u8; 4] = file[24..28].try_into().unwrap();
    file[92..96].copy_from_slice(&counter);
    file[28..32].copy_from_slice(&247u32.to_be_bytes());
    assert!(TaxonomySqlite::open(&file).is_ok());
    file[28..32].copy_from_slice(&248u32.to_be_bytes());
    let error = TaxonomySqlite::open(&file).err().expect("248 pages on 247");
    assert!(error.contains("malformed"), "{error}");
    file.extend_from_slice(&[0; 100]);
    assert!(TaxonomySqlite::open(&file).is_ok());
}
