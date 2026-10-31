//! The names of the named databases of an LMDB file (the `.pdb`/`.ndb` file of a version 5
//! BLAST database): what `mdb_dbi_open` would find. NCBI's reader opens the database
//! `acc2oid` to look an accession up, and a database built without `-parse_seqids` has no
//! such database (`DB contains no accession info.`). This reads the main database's B-tree
//! alone, page by page, with no other part of the file: a reader of the whole file (keys,
//! values, duplicates) is line DI-8's `seqdb::lmdb_read`, which line BDC's BDC2b uses.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1036-1048
// ```c
// typedef struct MDB_meta {
//     uint32_t mm_magic;
//     uint32_t mm_version;
//     void    *mm_address;        /**< address for fixed mapping */
//     size_t   mm_mapsize;        /**< size of mmap region */
//     MDB_db   mm_dbs[CORE_DBS];  /**< first is free space, 2nd is main db */
//     pgno_t   mm_last_pg;
//     volatile txnid_t mm_txnid;  /**< txnid that committed this page */
// } MDB_meta;
// ```
// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1005-1014
// ```c
// typedef struct MDB_db {
//     uint32_t md_pad;      /**< also ksize for LEAF2 pages */
//     uint16_t md_flags;    /**< @ref mdb_dbi_open */
//     uint16_t md_depth;    /**< depth of this tree */
//     pgno_t   md_branch_pages;
//     pgno_t   md_leaf_pages;
//     pgno_t   md_overflow_pages;
//     size_t   md_entries;
//     pgno_t   md_root;     /**< the root page of this tree */
// } MDB_db;
// ```
/// `PAGEHDRSZ`: the page header (`mp_pgno`, `mp_pad`, `mp_flags`, `mp_lower`, `mp_upper`).
const PAGE_HEADER: usize = 16;
/// `offsetof(MDB_page, mp_flags)`.
const PAGE_FLAGS: usize = 10;
/// `offsetof(MDB_page, mp_lower)`.
const PAGE_LOWER: usize = 12;
/// The meta structure starts after the page header; `mm_dbs[0].md_pad` (the page size,
/// `mm_psize`) at 24 within it, `mm_dbs[MAIN_DBI].md_root` at 72 + 40, `mm_txnid` at 128.
const META_MAGIC: usize = PAGE_HEADER;
const META_PSIZE: usize = PAGE_HEADER + 24;
const META_MAIN_ROOT: usize = PAGE_HEADER + 72 + 40;
const META_TXNID: usize = PAGE_HEADER + 128;
const META_END: usize = PAGE_HEADER + 136;
const MDB_MAGIC: u32 = 0xBEEF_C0DE;
const MDB_DATA_VERSION: u32 = 1;
/// `P_INVALID`: no root page.
const P_INVALID: u64 = u64::MAX;
const P_BRANCH: u16 = 0x01;
const P_LEAF: u16 = 0x02;
const P_META: u16 = 0x08;
/// `NODESIZE`: `mn_lo`, `mn_hi`, `mn_flags`, `mn_ksize` before the key.
const NODE_SIZE: usize = 8;
/// `F_SUBDATA`: the node's data is a named database (`MDB_db`).
const F_SUBDATA: u16 = 0x02;
/// The depth NCBI's LMDB bounds a cursor at (`CURSOR_STACK`).
const MAX_DEPTH: usize = 32;

fn u16_at(bytes: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([bytes[at], bytes[at + 1]])
}

fn u32_at(bytes: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
}

fn u64_at(bytes: &[u8], at: usize) -> u64 {
    let mut word = [0u8; 8];
    word.copy_from_slice(&bytes[at..at + 8]);
    u64::from_le_bytes(word)
}

struct Pages {
    file: File,
    psize: usize,
    len: u64,
}

impl Pages {
    fn read(&mut self, pgno: u64) -> Result<Vec<u8>, String> {
        let start = pgno
            .checked_mul(self.psize as u64)
            .filter(|&start| start + self.psize as u64 <= self.len)
            .ok_or_else(|| format!("page {pgno} is outside the file"))?;
        self.file
            .seek(SeekFrom::Start(start))
            .map_err(|error| error.to_string())?;
        let mut page = vec![0u8; self.psize];
        self.file
            .read_exact(&mut page)
            .map_err(|error| error.to_string())?;
        if u64_at(&page, 0) != pgno {
            return Err(format!("page {pgno} does not carry its own number"));
        }
        Ok(page)
    }
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4164-4186
// ```c
// static MDB_meta *
// mdb_env_pick_meta(const MDB_env *env)
// {
// 	MDB_meta *const *metas = env->me_metas;
// 	return metas[ (metas[0]->mm_txnid < metas[1]->mm_txnid) ^
// 		(env->me_flags & MDB_PREVSNAPSHOT) ];
// }
// ```
/// The names of the named databases of the LMDB file `path`, in key order: the keys of the
/// main database whose node holds a database (`F_SUBDATA`), from the meta page of the newer
/// transaction. `Err`: a file that is not an LMDB data file or a tree LOSAT cannot follow.
pub(super) fn sub_database_names(path: &Path) -> Result<Vec<Vec<u8>>, String> {
    let mut file = File::open(path).map_err(|error| error.to_string())?;
    let len = file.metadata().map_err(|error| error.to_string())?.len();
    let mut head = vec![0u8; META_END];
    file.read_exact(&mut head)
        .map_err(|_| "the file is shorter than a meta page".to_string())?;
    if u32_at(&head, META_MAGIC) != MDB_MAGIC {
        return Err("not an LMDB data file (bad magic number)".to_string());
    }
    if u32_at(&head, META_MAGIC + 4) != MDB_DATA_VERSION {
        return Err("an LMDB data file of another version".to_string());
    }
    let psize = u32_at(&head, META_PSIZE) as usize;
    if psize < META_END || psize > 0x8000 || !psize.is_power_of_two() {
        return Err(format!("an LMDB page size of {psize} bytes"));
    }
    let mut pages = Pages { file, psize, len };
    let meta0 = pages.read(0)?;
    let meta1 = pages.read(1)?;
    let meta = if u64_at(&meta0, META_TXNID) < u64_at(&meta1, META_TXNID) {
        &meta1
    } else {
        &meta0
    };
    if u16_at(meta, PAGE_FLAGS) & P_META == 0 {
        return Err("the meta page is not marked as one".to_string());
    }
    let root = u64_at(meta, META_MAIN_ROOT);
    let mut names = Vec::new();
    if root != P_INVALID {
        walk(&mut pages, root, 0, &mut names)?;
    }
    Ok(names)
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1240-1262
// ```c
// 	/** Address of node \b i in page \b p */
// #define NODEPTR(p, i)	 ((MDB_node *)((char *)(p) + (p)->mp_ptrs[i] + PAGEBASE))
// 	/** Address of the key for the node */
// #define NODEKEY(node)	 (void *)((node)->mn_data)
// 	/** Get the page number pointed to by a branch node */
// #define NODEPGNO(node) \
// 	((node)->mn_lo | ((pgno_t) (node)->mn_hi << 16) | \
// 	 (PGNO_TOPBITS ? ((pgno_t) (node)->mn_flags << 32) : 0))
// ```
/// Collects the keys of the named databases below page `pgno`.
fn walk(
    pages: &mut Pages,
    pgno: u64,
    depth: usize,
    names: &mut Vec<Vec<u8>>,
) -> Result<(), String> {
    if depth > MAX_DEPTH {
        return Err("a tree deeper than LMDB allows".to_string());
    }
    let page = pages.read(pgno)?;
    let flags = u16_at(&page, PAGE_FLAGS);
    let lower = usize::from(u16_at(&page, PAGE_LOWER));
    if lower < PAGE_HEADER || lower > pages.psize {
        return Err(format!("page {pgno} has a bad lower bound"));
    }
    let numkeys = (lower - PAGE_HEADER) / 2;
    let psize = pages.psize;
    let node_at = |i: usize| -> Result<usize, String> {
        let at = usize::from(u16_at(&page, PAGE_HEADER + 2 * i));
        if at + NODE_SIZE > psize {
            return Err(format!("page {pgno} points outside itself"));
        }
        Ok(at)
    };
    if flags & P_BRANCH != 0 {
        for i in 0..numkeys {
            let node = node_at(i)?;
            let child = u64::from(u16_at(&page, node))
                | (u64::from(u16_at(&page, node + 2)) << 16)
                | (u64::from(u16_at(&page, node + 4)) << 32);
            walk(pages, child, depth + 1, names)?;
        }
    } else if flags & P_LEAF != 0 {
        for i in 0..numkeys {
            let node = node_at(i)?;
            let node_flags = u16_at(&page, node + 4);
            let ksize = usize::from(u16_at(&page, node + 6));
            if node + NODE_SIZE + ksize > psize {
                return Err(format!("page {pgno} holds a key past its end"));
            }
            if node_flags & F_SUBDATA != 0 {
                names.push(page[node + NODE_SIZE..node + NODE_SIZE + ksize].to_vec());
            }
        }
    } else {
        return Err(format!("page {pgno} is neither a branch nor a leaf"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/blastdbcmd/db")
            .join(name)
    }

    #[test]
    fn a_database_without_parsed_ids_has_no_acc2oid() {
        let names = sub_database_names(&fixture("n5n.ndb")).unwrap();
        assert_eq!(names, [b"volinfo".to_vec(), b"volname".to_vec()]);
        let names = sub_database_names(&fixture("p5.pdb")).unwrap();
        assert_eq!(
            names,
            [
                b"acc2oid".to_vec(),
                b"volinfo".to_vec(),
                b"volname".to_vec()
            ]
        );
    }

    #[test]
    fn other_files_are_refused() {
        assert!(sub_database_names(&fixture("p5.pin")).is_err());
        assert!(sub_database_names(&fixture("nosuch.pdb")).is_err());
    }
}
