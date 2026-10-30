//! Read-only reader of LMDB 0.9.24 files (64-bit, little-endian) over positioned reads.
//!
//! Scope: what NCBI's BLAST database reader does with the `.ndb/.pdb` and `.ntf/.ptf` files
//! (`MDB_NOSUBDIR|MDB_NOLOCK|MDB_RDONLY`): pick the newest meta page, find a named database in
//! the main database, look a key up (`MDB_SET`, first duplicate or all duplicates) and walk a
//! database in tree order (`MDB_FIRST`/`MDB_NEXT`). No transaction, free list, lock file or
//! writer. Pages are read one at a time through [`ReadAt`] (a small page cache keeps the last
//! few hundred); the file is never held in memory. Every offset is bounds checked, a corrupt
//! file gives `Err` and never a panic or an endless walk.
//!
//! Duplicate shapes handled (`MDB_DUPSORT`, with or without `MDB_DUPFIXED`): a plain node (one
//! value), an inline sub-page (`F_DUPDATA`) and a sub-database (`F_DUPDATA|F_SUBDATA`).
//! `F_BIGDATA` (overflow page) values are read for plain nodes.
//!
//! Difference from NCBI: NCBI fails on keys that hold a sub-database (more than about 500
//! duplicates) because its `MDB_GET_CURRENT` runs on an uninitialised sub-cursor
//! (`mdb_xcursor_init1` sets only `C_SUB` for `F_SUBDATA` nodes). [`Env::first_value`] and
//! [`Env::values`] work for those keys.

use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

/// Positioned reads over a database file ([`super::file::DbFile`]) or a buffer.
pub trait ReadAt {
    /// Total length in bytes.
    fn len(&self) -> u64;
    /// True when the source has no bytes.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Fill `buf` from `offset`; `Err` when the range is not inside the source.
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String>;
}

impl ReadAt for [u8] {
    fn len(&self) -> u64 {
        <[u8]>::len(self) as u64
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        let end = offset
            .checked_add(buf.len() as u64)
            .filter(|&e| e <= <[u8]>::len(self) as u64)
            .ok_or_else(|| {
                format!(
                    "read of {} bytes at offset {} is beyond the end of the file ({} bytes)",
                    buf.len(),
                    offset,
                    <[u8]>::len(self)
                )
            })?;
        buf.copy_from_slice(&self[offset as usize..end as usize]);
        Ok(())
    }
}

impl ReadAt for Vec<u8> {
    fn len(&self) -> u64 {
        Vec::len(self) as u64
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        self.as_slice().read_exact_at(buf, offset)
    }
}

impl ReadAt for super::file::DbFile {
    fn len(&self) -> u64 {
        super::file::DbFile::len(self)
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        super::file::DbFile::read_exact_at(self, buf, offset).map_err(|e| e.to_string())
    }
}

/// Plain files, for the tests (the library reads through [`super::file::DbFile`]).
#[cfg(test)]
impl ReadAt for std::fs::File {
    fn len(&self) -> u64 {
        self.metadata().map(|m| m.len()).unwrap_or(0)
    }
    #[cfg(unix)]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        std::os::unix::fs::FileExt::read_exact_at(self, buf, offset)
            .map_err(|e| format!("read of {} bytes at offset {}: {e}", buf.len(), offset))
    }
    #[cfg(not(any(unix, windows)))]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        use std::io::{Read, Seek, SeekFrom};
        let mut f = self;
        f.seek(SeekFrom::Start(offset))
            .and_then(|_| f.read_exact(buf))
            .map_err(|e| format!("read of {} bytes at offset {}: {e}", buf.len(), offset))
    }
    #[cfg(windows)]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
        let mut done = 0usize;
        while done < buf.len() {
            let n = std::os::windows::fs::FileExt::seek_read(
                self,
                &mut buf[done..],
                offset + done as u64,
            )
            .map_err(|e| format!("read of {} bytes at offset {}: {e}", buf.len(), offset))?;
            if n == 0 {
                return Err(format!(
                    "read of {} bytes at offset {}: unexpected end of file",
                    buf.len(),
                    offset
                ));
            }
            done += n;
        }
        Ok(())
    }
}

/// `MDB_DUPSORT` database flag.
pub const MDB_DUPSORT: u16 = 0x04;
/// `MDB_INTEGERKEY` database flag.
pub const MDB_INTEGERKEY: u16 = 0x08;
/// `MDB_DUPFIXED` database flag.
pub const MDB_DUPFIXED: u16 = 0x10;

const MDB_MAGIC: u32 = 0xBEEF_C0DE;
const MDB_DATA_VERSION: u32 = 1;
const P_INVALID: u64 = u64::MAX;
const P_BRANCH: u16 = 0x01;
const P_LEAF: u16 = 0x02;
const P_OVERFLOW: u16 = 0x04;
const P_META: u16 = 0x08;
const P_LEAF2: u16 = 0x20;
const P_SUBP: u16 = 0x40;
const F_BIGDATA: u16 = 0x01;
const F_SUBDATA: u16 = 0x02;
const F_DUPDATA: u16 = 0x04;
const PAGEHDRSZ: usize = 16;
const NODESZ: usize = 8;
const DB_RECORD_SZ: usize = 48;
const META_BYTES: usize = 152;
const MAX_DEPTH: usize = 32;
const CACHE_PAGES: usize = 256;

fn bad_valsize() -> String {
    "MDB_BAD_VALSIZE: Unsupported size of key/DB name/data, or wrong DUPFIXED size".to_string()
}

fn slice(b: &[u8], off: usize, n: usize) -> Result<&[u8], String> {
    off.checked_add(n)
        .and_then(|end| b.get(off..end))
        .ok_or_else(|| {
            format!(
                "LMDB: {n} bytes at offset {off} are outside a block of {} bytes",
                b.len()
            )
        })
}

fn rd_u16(b: &[u8], off: usize) -> Result<u16, String> {
    Ok(u16::from_le_bytes(
        slice(b, off, 2)?
            .try_into()
            .map_err(|_| "LMDB: short read")?,
    ))
}

fn rd_u32(b: &[u8], off: usize) -> Result<u32, String> {
    Ok(u32::from_le_bytes(
        slice(b, off, 4)?
            .try_into()
            .map_err(|_| "LMDB: short read")?,
    ))
}

fn rd_u64(b: &[u8], off: usize) -> Result<u64, String> {
    Ok(u64::from_le_bytes(
        slice(b, off, 8)?
            .try_into()
            .map_err(|_| "LMDB: short read")?,
    ))
}

/// The `MDB_db` record of a database (48 bytes in the meta page, or as the data of a named
/// database / sub-database node).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Db {
    /// `md_pad`: page size in `FREE_DBI`, key size of LEAF2 pages in sub-databases.
    pub pad: u32,
    /// `md_flags`: the stored database flags (`MDB_DUPSORT`, `MDB_INTEGERKEY`, ...).
    pub flags: u16,
    /// `md_depth`: depth of the tree.
    pub depth: u16,
    /// `md_branch_pages`.
    pub branch_pages: u64,
    /// `md_leaf_pages`.
    pub leaf_pages: u64,
    /// `md_overflow_pages`.
    pub overflow_pages: u64,
    /// `md_entries`: number of data items (every duplicate counts).
    pub entries: u64,
    /// `md_root`: root page, `u64::MAX` for an empty tree.
    pub root: u64,
}

impl Db {
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1005-1014
    // ```c++
    // typedef struct MDB_db {
    //     uint32_t    md_pad;        /**< also ksize for LEAF2 pages */
    //     uint16_t    md_flags;    /**< @ref mdb_dbi_open */
    //     uint16_t    md_depth;    /**< depth of this tree */
    //     pgno_t        md_branch_pages;    /**< number of internal pages */
    //     pgno_t        md_leaf_pages;        /**< number of leaf pages */
    //     pgno_t        md_overflow_pages;    /**< number of overflow pages */
    //     size_t        md_entries;        /**< number of data items */
    //     pgno_t        md_root;        /**< the root page of this tree */
    // } MDB_db;
    // ```
    fn parse(raw: &[u8]) -> Result<Db, String> {
        if raw.len() < DB_RECORD_SZ {
            return Err(format!(
                "LMDB: database record of {} bytes, expected {DB_RECORD_SZ}",
                raw.len()
            ));
        }
        Ok(Db {
            pad: rd_u32(raw, 0)?,
            flags: rd_u16(raw, 4)?,
            depth: rd_u16(raw, 6)?,
            branch_pages: rd_u64(raw, 8)?,
            leaf_pages: rd_u64(raw, 16)?,
            overflow_pages: rd_u64(raw, 24)?,
            entries: rd_u64(raw, 32)?,
            root: rd_u64(raw, 40)?,
        })
    }

    /// True when the tree has no page (`md_root == P_INVALID`).
    pub fn is_empty_tree(&self) -> bool {
        self.root == P_INVALID
    }
}

struct Meta {
    free: Db,
    main: Db,
    last_pg: u64,
    txnid: u64,
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1036-1054
// ```c++
// typedef struct MDB_meta {
//         /** Stamp identifying this as an LMDB file. It must be set
//          *    to #MDB_MAGIC. */
//     uint32_t    mm_magic;
//         /** Version number of this file. Must be set to #MDB_DATA_VERSION. */
//     uint32_t    mm_version;
//     void        *mm_address;        /**< address for fixed mapping */
//     size_t        mm_mapsize;            /**< size of mmap region */
//     MDB_db        mm_dbs[CORE_DBS];    /**< first is free space, 2nd is main db */
//     /** The size of pages used in this DB */
// #define    mm_psize    mm_dbs[FREE_DBI].md_pad
//     /** Any persistent environment flags. @ref mdb_env */
// #define    mm_flags    mm_dbs[FREE_DBI].md_flags
//     /** Last used page in the datafile.
//      *    Actually the file may be shorter if the freeDB lists the final pages.
//      */
//     pgno_t        mm_last_pg;
//     volatile txnid_t    mm_txnid;    /**< txnid that committed this page */
// } MDB_meta;
// ```
/// Read and validate the meta page that starts at `off` (header 16 bytes, then the meta).
fn read_meta<R: ReadAt + ?Sized>(src: &R, off: u64) -> Result<Meta, String> {
    let end = off
        .checked_add(META_BYTES as u64)
        .filter(|&e| e <= src.len());
    if end.is_none() {
        return Err(format!(
            "LMDB: file of {} bytes is too small for a meta page at offset {off}",
            src.len()
        ));
    }
    let mut b = [0u8; META_BYTES];
    src.read_exact_at(&mut b, off)?;
    // The checks of mdb_env_read_header: P_META flag, magic, version.
    if rd_u16(&b, 10)? & P_META == 0 {
        return Err(format!(
            "LMDB: page at offset {off} is not a meta page (MDB_INVALID)"
        ));
    }
    if rd_u32(&b, 16)? != MDB_MAGIC {
        return Err("LMDB: meta page has invalid magic (MDB_INVALID)".to_string());
    }
    let version = rd_u32(&b, 20)?;
    if version != MDB_DATA_VERSION {
        return Err(format!(
            "LMDB: database is version {version}, expected version {MDB_DATA_VERSION} (MDB_VERSION_MISMATCH)"
        ));
    }
    Ok(Meta {
        free: Db::parse(slice(&b, 40, DB_RECORD_SZ)?)?,
        main: Db::parse(slice(&b, 88, DB_RECORD_SZ)?)?,
        last_pg: rd_u64(&b, 136)?,
        txnid: rd_u64(&b, 144)?,
    })
}

/// A page (or an inline sub-page) with its parsed header.
struct Page {
    buf: Rc<Vec<u8>>,
    flags: u16,
    pad: u16,
    n: usize,
}

/// One node of a branch or leaf page (offsets into the page buffer).
struct Node {
    flags: u16,
    koff: usize,
    ksize: usize,
    doff: usize,
    dsize: usize,
    child: u64,
}

impl Page {
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:846-852
    // ```c++
    // #define METADATA(p)     ((void *)((char *)(p) + PAGEHDRSZ))
    //
    //     /** ITS#7713, change PAGEBASE to handle 65536 byte pages */
    // #define    PAGEBASE    ((MDB_DEVEL) ? PAGEHDRSZ : 0)
    //
    //     /** Number of nodes on a page */
    // #define NUMKEYS(p)     (((p)->mp_lower - (PAGEHDRSZ-PAGEBASE)) >> 1)
    // ```
    fn new(buf: Rc<Vec<u8>>) -> Result<Page, String> {
        let pad = rd_u16(&buf, 8)?;
        let flags = rd_u16(&buf, 10)?;
        let lower = rd_u16(&buf, 12)? as usize;
        if lower < PAGEHDRSZ || lower > buf.len() {
            return Err(format!(
                "LMDB: page header lower={lower} is outside a page of {} bytes",
                buf.len()
            ));
        }
        Ok(Page {
            buf,
            flags,
            pad,
            n: (lower - PAGEHDRSZ) >> 1,
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:898-960
    // ```c++
    // typedef struct MDB_node {
    //     /** part of data size or pgno
    //      *    @{ */
    // #if BYTE_ORDER == LITTLE_ENDIAN
    //     unsigned short    mn_lo, mn_hi;
    // #else
    //     unsigned short    mn_hi, mn_lo;
    // #endif
    //     /** @} */
    // /** @defgroup mdb_node Node Flags
    //  *    @ingroup internal
    //  *    Flags for node headers.
    //  *    @{
    //  */
    // #define F_BIGDATA     0x01            /**< data put on overflow page */
    // #define F_SUBDATA     0x02            /**< data is a sub-database */
    // #define F_DUPDATA     0x04            /**< data has duplicates */
    //
    // /** valid flags for #mdb_node_add() */
    // #define    NODE_ADD_FLAGS    (F_DUPDATA|F_SUBDATA|MDB_RESERVE|MDB_APPEND)
    //
    // /** @} */
    //     unsigned short    mn_flags;        /**< @ref mdb_node */
    //     unsigned short    mn_ksize;        /**< key size */
    //     char        mn_data[1];            /**< key and data are appended here */
    // } MDB_node;
    //
    //     /** Size of the node header, excluding dynamic data at the end */
    // #define NODESIZE     offsetof(MDB_node, mn_data)
    //
    //     /** Bit position of top word in page number, for shifting mn_flags */
    // #define PGNO_TOPWORD ((pgno_t)-1 > 0xffffffffu ? 32 : 0)
    //
    //     /** Size of a node in a branch page with a given key.
    //      *    This is just the node header plus the key, there is no data.
    //      */
    // #define INDXSIZE(k)     (NODESIZE + ((k) == NULL ? 0 : (k)->mv_size))
    //
    //     /** Size of a node in a leaf page with a given key and data.
    //      *    This is node header plus key plus data size.
    //      */
    // #define LEAFSIZE(k, d)     (NODESIZE + (k)->mv_size + (d)->mv_size)
    //
    //     /** Address of node \b i in page \b p */
    // #define NODEPTR(p, i)     ((MDB_node *)((char *)(p) + (p)->mp_ptrs[i] + PAGEBASE))
    //
    //     /** Address of the key for the node */
    // #define NODEKEY(node)     (void *)((node)->mn_data)
    //
    //     /** Address of the data for a node */
    // #define NODEDATA(node)     (void *)((char *)(node)->mn_data + (node)->mn_ksize)
    //
    //     /** Get the page number pointed to by a branch node */
    // #define NODEPGNO(node) \
    //     ((node)->mn_lo | ((pgno_t) (node)->mn_hi << 16) | \
    //      (PGNO_TOPWORD ? ((pgno_t) (node)->mn_flags << PGNO_TOPWORD) : 0))
    //     /** Set the page number in a branch node */
    // #define SETPGNO(node,pgno)    do { \
    //     (node)->mn_lo = (pgno) & 0xffff; (node)->mn_hi = (pgno) >> 16; \
    //     if (PGNO_TOPWORD) (node)->mn_flags = (pgno) >> PGNO_TOPWORD; } while(0)
    //
    //     /** Get the size of the data in a leaf node */
    // #define NODEDSZ(node)     ((node)->mn_lo | ((unsigned)(node)->mn_hi << 16))
    // ```
    fn node(&self, i: usize) -> Result<Node, String> {
        if self.flags & P_LEAF2 != 0 || self.flags & (P_BRANCH | P_LEAF) == 0 {
            return Err(format!(
                "LMDB: page with flags {:#x} has no nodes",
                self.flags
            ));
        }
        if i >= self.n {
            return Err(format!(
                "LMDB: node {i} is outside a page of {} nodes",
                self.n
            ));
        }
        let o = rd_u16(&self.buf, PAGEHDRSZ + 2 * i)? as usize;
        let lo = rd_u16(&self.buf, o)? as u64;
        let hi = rd_u16(&self.buf, o + 2)? as u64;
        let flags = rd_u16(&self.buf, o + 4)?;
        let ksize = rd_u16(&self.buf, o + 6)? as usize;
        let koff = o + NODESZ;
        slice(&self.buf, koff, ksize)?;
        if self.flags & P_BRANCH != 0 {
            // The flags field is the top word of the child page number.
            return Ok(Node {
                flags,
                koff,
                ksize,
                doff: 0,
                dsize: 0,
                child: lo | (hi << 16) | ((flags as u64) << 32),
            });
        }
        let dsize = (lo | (hi << 16)) as usize;
        let doff = koff + ksize;
        // For F_BIGDATA the data of the node is the 8-byte number of the overflow page.
        slice(
            &self.buf,
            doff,
            if flags & F_BIGDATA != 0 { 8 } else { dsize },
        )?;
        Ok(Node {
            flags,
            koff,
            ksize,
            doff,
            dsize,
            child: 0,
        })
    }

    fn key(&self, n: &Node) -> &[u8] {
        self.buf.get(n.koff..n.koff + n.ksize).unwrap_or(&[])
    }

    /// Node data as stored in the page (the overflow page number for F_BIGDATA).
    fn data(&self, n: &Node) -> &[u8] {
        let len = if n.flags & F_BIGDATA != 0 { 8 } else { n.dsize };
        self.buf.get(n.doff..n.doff + len).unwrap_or(&[])
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:995-995
    // ```c++
    // #define LEAF2KEY(p, i, ks)    ((char *)(p) + PAGEHDRSZ + ((i)*(ks)))
    // ```
    fn leaf2_key(&self, i: usize) -> Result<&[u8], String> {
        let ks = self.pad as usize;
        if ks == 0 {
            return Err("LMDB: LEAF2 page with key size 0".to_string());
        }
        if i >= self.n {
            return Err(format!(
                "LMDB: LEAF2 key {i} is outside a page of {} keys",
                self.n
            ));
        }
        slice(&self.buf, PAGEHDRSZ + i * ks, ks)
    }

    /// Key of entry `i` of a leaf page of a sub-database or sub-page: LEAF2 keys are the
    /// fixed-size values, in other leaves the value of a duplicate is the node key.
    fn dup_item(&self, i: usize) -> Result<&[u8], String> {
        if self.flags & P_LEAF2 != 0 {
            self.leaf2_key(i)
        } else {
            let n = self.node(i)?;
            Ok(self.key(&n))
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5276-5365
    // ```c++
    // static MDB_node *
    // mdb_node_search(MDB_cursor *mc, MDB_val *key, int *exactp)
    // {
    //     unsigned int     i = 0, nkeys;
    //     int         low, high;
    //     int         rc = 0;
    //     MDB_page *mp = mc->mc_pg[mc->mc_top];
    //     MDB_node    *node = NULL;
    //     MDB_val     nodekey;
    //     MDB_cmp_func *cmp;
    //     DKBUF;
    //
    //     nkeys = NUMKEYS(mp);
    //
    //     DPRINTF(("searching %u keys in %s %spage %"Z"u",
    //         nkeys, IS_LEAF(mp) ? "leaf" : "branch", IS_SUBP(mp) ? "sub-" : "",
    //         mdb_dbg_pgno(mp)));
    //
    //     low = IS_LEAF(mp) ? 0 : 1;
    //     high = nkeys - 1;
    //     cmp = mc->mc_dbx->md_cmp;
    //
    //     /* Branch pages have no data, so if using integer keys,
    //      * alignment is guaranteed. Use faster mdb_cmp_int.
    //      */
    //     if (cmp == mdb_cmp_cint && IS_BRANCH(mp)) {
    //         if (NODEPTR(mp, 1)->mn_ksize == sizeof(size_t))
    //             cmp = mdb_cmp_long;
    //         else
    //             cmp = mdb_cmp_int;
    //     }
    //
    //     if (IS_LEAF2(mp)) {
    //         nodekey.mv_size = mc->mc_db->md_pad;
    //         node = NODEPTR(mp, 0);    /* fake */
    //         while (low <= high) {
    //             i = (low + high) >> 1;
    //             nodekey.mv_data = LEAF2KEY(mp, i, nodekey.mv_size);
    //             rc = cmp(key, &nodekey);
    //             DPRINTF(("found leaf index %u [%s], rc = %i",
    //                 i, DKEY(&nodekey), rc));
    //             if (rc == 0)
    //                 break;
    //             if (rc > 0)
    //                 low = i + 1;
    //             else
    //                 high = i - 1;
    //         }
    //     } else {
    //         while (low <= high) {
    //             i = (low + high) >> 1;
    //
    //             node = NODEPTR(mp, i);
    //             nodekey.mv_size = NODEKSZ(node);
    //             nodekey.mv_data = NODEKEY(node);
    //
    //             rc = cmp(key, &nodekey);
    // #if MDB_DEBUG
    //             if (IS_LEAF(mp))
    //                 DPRINTF(("found leaf index %u [%s], rc = %i",
    //                     i, DKEY(&nodekey), rc));
    //             else
    //                 DPRINTF(("found branch index %u [%s -> %"Z"u], rc = %i",
    //                     i, DKEY(&nodekey), NODEPGNO(node), rc));
    // #endif
    //             if (rc == 0)
    //                 break;
    //             if (rc > 0)
    //                 low = i + 1;
    //             else
    //                 high = i - 1;
    //         }
    //     }
    //
    //     if (rc > 0) {    /* Found entry is less than the key. */
    //         i++;    /* Skip to get the smallest entry larger than key. */
    //         if (!IS_LEAF2(mp))
    //             node = NODEPTR(mp, i);
    //     }
    //     if (exactp)
    //         *exactp = (rc == 0 && nkeys > 0);
    //     /* store the key index */
    //     mc->mc_ki[mc->mc_top] = i;
    //     if (i >= nkeys)
    //         /* There is no entry larger or equal to the key. */
    //         return NULL;
    //
    //     /* nodeptr is fake for LEAF2 */
    //     return node;
    // }
    // ```
    /// `mdb_node_search`: index of the first entry >= `key` and whether it is equal.
    fn node_search(&self, key: &[u8], db_flags: u16) -> Result<(usize, bool), String> {
        let leaf2 = self.flags & P_LEAF2 != 0;
        let mut low: isize = if self.flags & P_LEAF != 0 { 0 } else { 1 };
        let mut high: isize = self.n as isize - 1;
        let mut i = 0usize;
        let mut rc = Ordering::Equal;
        while low <= high {
            i = ((low + high) >> 1) as usize;
            rc = if leaf2 {
                cmp_keys(db_flags, key, self.leaf2_key(i)?)
            } else {
                let node = self.node(i)?;
                cmp_keys(db_flags, key, self.key(&node))
            };
            match rc {
                Ordering::Equal => break,
                Ordering::Greater => low = i as isize + 1,
                Ordering::Less => high = i as isize - 1,
            }
        }
        if rc == Ordering::Greater {
            // Found entry is less than the key: skip to the smallest entry larger than key.
            i += 1;
        }
        Ok((i, rc == Ordering::Equal && self.n > 0))
    }
}

/// Key order by the stored database flags: `MDB_INTEGERKEY` = unsigned little-endian integer
/// (`mdb_cmp_cint`), otherwise `memcmp`, shorter key first (`mdb_cmp_memn`). The flags of the
/// open call do not matter: `mdb_dbi_open` takes them from the stored record.
fn cmp_keys(db_flags: u16, a: &[u8], b: &[u8]) -> Ordering {
    if db_flags & MDB_INTEGERKEY != 0 {
        cmp_cint(a, b)
    } else {
        cmp_memn(a, b)
    }
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5197-5217
// ```c++
// static int
// mdb_cmp_cint(const MDB_val *a, const MDB_val *b)
// {
// #if BYTE_ORDER == LITTLE_ENDIAN
//     unsigned short *u, *c;
//     int x;
//
//     u = (unsigned short *) ((char *) a->mv_data + a->mv_size);
//     c = (unsigned short *) ((char *) b->mv_data + a->mv_size);
//     do {
//         x = *--u - *--c;
//     } while(!x && u > (unsigned short *)a->mv_data);
//     return x;
// #else
//     unsigned short *u, *c, *end;
//     int x;
//
//     end = (unsigned short *) ((char *) a->mv_data + a->mv_size);
//     u = (unsigned short *)a->mv_data;
//     c = (unsigned short *)b->mv_data;
//     do {
// ```
fn cmp_cint(a: &[u8], b: &[u8]) -> Ordering {
    // Keys of one database have one size (4 bytes here); NCBI reads `a->mv_size` bytes of both.
    if a.len() != b.len() {
        return a.len().cmp(&b.len());
    }
    a.iter().rev().cmp(b.iter().rev())
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5225-5241
// ```c++
// static int
// mdb_cmp_memn(const MDB_val *a, const MDB_val *b)
// {
//     int diff;
//     ssize_t len_diff;
//     unsigned int len;
//
//     len = a->mv_size;
//     len_diff = (ssize_t) a->mv_size - (ssize_t) b->mv_size;
//     if (len_diff > 0) {
//         len = b->mv_size;
//         len_diff = 1;
//     }
//
//     diff = memcmp(a->mv_data, b->mv_data, len);
//     return diff ? diff : len_diff<0 ? -1 : len_diff;
// }
// ```
fn cmp_memn(a: &[u8], b: &[u8]) -> Ordering {
    a.cmp(b)
}

/// How the values of one key are stored.
enum Dups {
    /// A plain node: one value.
    One(Vec<u8>),
    /// An inline sub-page (`F_DUPDATA`).
    SubPage(Page),
    /// A sub-database (`F_DUPDATA|F_SUBDATA`).
    SubDb(Db),
}

/// Read-only LMDB environment over a [`ReadAt`] source.
pub struct Env<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    psize: usize,
    npages: u64,
    main: Db,
    txnid: u64,
    cache: RefCell<HashMap<u64, Rc<Vec<u8>>>>,
}

impl<R: ReadAt + ?Sized> std::fmt::Debug for Env<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Env")
            .field("psize", &self.psize)
            .field("npages", &self.npages)
            .field("txnid", &self.txnid)
            .field("main", &self.main)
            .finish()
    }
}

impl<'r, R: ReadAt + ?Sized> Env<'r, R> {
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3682-3738
    // ```c++
    // static int ESECT
    // mdb_env_read_header(MDB_env *env, MDB_meta *meta)
    // {
    //     MDB_metabuf    pbuf;
    //     MDB_page    *p;
    //     MDB_meta    *m;
    //     int            i, rc, off;
    //     enum { Size = sizeof(pbuf) };
    //
    //     /* We don't know the page size yet, so use a minimum value.
    //      * Read both meta pages so we can use the latest one.
    //      */
    //
    //     for (i=off=0; i<NUM_METAS; i++, off += meta->mm_psize) {
    // #ifdef _WIN32
    //         DWORD len;
    //         OVERLAPPED ov;
    //         memset(&ov, 0, sizeof(ov));
    //         ov.Offset = off;
    //         rc = ReadFile(env->me_fd, &pbuf, Size, &len, &ov) ? (int)len : -1;
    //         if (rc == -1 && ErrCode() == ERROR_HANDLE_EOF)
    //             rc = 0;
    // #else
    //         rc = pread(env->me_fd, &pbuf, Size, off);
    // #endif
    //         if (rc != Size) {
    //             if (rc == 0 && off == 0)
    //                 return ENOENT;
    //             rc = rc < 0 ? (int) ErrCode() : MDB_INVALID;
    //             DPRINTF(("read: %s", mdb_strerror(rc)));
    //             return rc;
    //         }
    //
    //         p = (MDB_page *)&pbuf;
    //
    //         if (!F_ISSET(p->mp_flags, P_META)) {
    //             DPRINTF(("page %"Z"u not a meta page", p->mp_pgno));
    //             return MDB_INVALID;
    //         }
    //
    //         m = METADATA(p);
    //         if (m->mm_magic != MDB_MAGIC) {
    //             DPUTS("meta has invalid magic");
    //             return MDB_INVALID;
    //         }
    //
    //         if (m->mm_version != MDB_DATA_VERSION) {
    //             DPRINTF(("database is version %u, expected version %u",
    //                 m->mm_version, MDB_DATA_VERSION));
    //             return MDB_VERSION_MISMATCH;
    //         }
    //
    //         if (off == 0 || m->mm_txnid > meta->mm_txnid)
    //             *meta = *m;
    //     }
    //     return 0;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3943-3949
    // ```c++
    // static MDB_meta *
    // mdb_env_pick_meta(const MDB_env *env)
    // {
    //     MDB_meta *const *metas = env->me_metas;
    //     return metas[ metas[0]->mm_txnid < metas[1]->mm_txnid ];
    // }
    //
    // ```
    /// Open an environment: validate both meta pages (magic, version), take the page size from
    /// `FREE_DBI.md_pad`, use the meta page with the larger transaction id (tie: meta 0).
    pub fn open(src: &'r R) -> Result<Self, String> {
        let m0 = read_meta(src, 0)?;
        let psize = m0.free.pad as usize;
        if !(512..=32768).contains(&psize) || !psize.is_power_of_two() {
            return Err(format!("LMDB: invalid page size {psize} in the meta page"));
        }
        let m1 = read_meta(src, psize as u64)?;
        let meta = if m0.txnid < m1.txnid { m1 } else { m0 };
        // mdb_page_get refuses pages >= mt_next_pgno = mm_last_pg + 1.
        let npages = meta.last_pg.saturating_add(1).min(src.len() / psize as u64);
        let _ = meta.free;
        Ok(Env {
            src,
            psize,
            npages,
            main: meta.main,
            txnid: meta.txnid,
            cache: RefCell::new(HashMap::new()),
        })
    }

    /// Page size of the file.
    pub fn page_size(&self) -> usize {
        self.psize
    }

    /// Transaction id of the meta page in use.
    pub fn txnid(&self) -> u64 {
        self.txnid
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5465-5479
    // ```c++
    //     if (pgno < txn->mt_next_pgno) {
    //         level = 0;
    //         p = (MDB_page *)(env->me_map + env->me_psize * pgno);
    //     } else {
    //         DPRINTF(("page %"Z"u not found", pgno));
    //         txn->mt_flags |= MDB_TXN_ERROR;
    //         return MDB_PAGE_NOTFOUND;
    //     }
    //
    // done:
    //     *ret = p;
    //     if (lvl)
    //         *lvl = level;
    //     return MDB_SUCCESS;
    // }
    // ```
    fn page(&self, pgno: u64) -> Result<Page, String> {
        if pgno >= self.npages {
            return Err(format!(
                "LMDB: page {pgno} not found (the database has {} pages) (MDB_PAGE_NOTFOUND)",
                self.npages
            ));
        }
        if let Some(buf) = self.cache.borrow().get(&pgno) {
            return Page::new(Rc::clone(buf));
        }
        let mut buf = vec![0u8; self.psize];
        self.src.read_exact_at(&mut buf, pgno * self.psize as u64)?;
        let buf = Rc::new(buf);
        let mut cache = self.cache.borrow_mut();
        if cache.len() >= CACHE_PAGES {
            cache.clear();
        }
        cache.insert(pgno, Rc::clone(&buf));
        Page::new(buf)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5484-5562
    // ```c++
    // static int
    // mdb_page_search_root(MDB_cursor *mc, MDB_val *key, int flags)
    // {
    //     MDB_page    *mp = mc->mc_pg[mc->mc_top];
    //     int rc;
    //     DKBUF;
    //
    //     while (IS_BRANCH(mp)) {
    //         MDB_node    *node;
    //         indx_t        i;
    //
    //         DPRINTF(("branch page %"Z"u has %u keys", mp->mp_pgno, NUMKEYS(mp)));
    //         /* Don't assert on branch pages in the FreeDB. We can get here
    //          * while in the process of rebalancing a FreeDB branch page; we must
    //          * let that proceed. ITS#8336
    //          */
    //         mdb_cassert(mc, !mc->mc_dbi || NUMKEYS(mp) > 1);
    //         DPRINTF(("found index 0 to page %"Z"u", NODEPGNO(NODEPTR(mp, 0))));
    //
    //         if (flags & (MDB_PS_FIRST|MDB_PS_LAST)) {
    //             i = 0;
    //             if (flags & MDB_PS_LAST) {
    //                 i = NUMKEYS(mp) - 1;
    //                 /* if already init'd, see if we're already in right place */
    //                 if (mc->mc_flags & C_INITIALIZED) {
    //                     if (mc->mc_ki[mc->mc_top] == i) {
    //                         mc->mc_top = mc->mc_snum++;
    //                         mp = mc->mc_pg[mc->mc_top];
    //                         goto ready;
    //                     }
    //                 }
    //             }
    //         } else {
    //             int     exact;
    //             node = mdb_node_search(mc, key, &exact);
    //             if (node == NULL)
    //                 i = NUMKEYS(mp) - 1;
    //             else {
    //                 i = mc->mc_ki[mc->mc_top];
    //                 if (!exact) {
    //                     mdb_cassert(mc, i > 0);
    //                     i--;
    //                 }
    //             }
    //             DPRINTF(("following index %u for key [%s]", i, DKEY(key)));
    //         }
    //
    //         mdb_cassert(mc, i < NUMKEYS(mp));
    //         node = NODEPTR(mp, i);
    //
    //         if ((rc = mdb_page_get(mc, NODEPGNO(node), &mp, NULL)) != 0)
    //             return rc;
    //
    //         mc->mc_ki[mc->mc_top] = i;
    //         if ((rc = mdb_cursor_push(mc, mp)))
    //             return rc;
    //
    // ready:
    //         if (flags & MDB_PS_MODIFY) {
    //             if ((rc = mdb_page_touch(mc)) != 0)
    //                 return rc;
    //             mp = mc->mc_pg[mc->mc_top];
    //         }
    //     }
    //
    //     if (!IS_LEAF(mp)) {
    //         DPRINTF(("internal error, index points to a %02X page!?",
    //             mp->mp_flags));
    //         mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    //         return MDB_CORRUPTED;
    //     }
    //
    //     DPRINTF(("found leaf page %"Z"u for key [%s]", mp->mp_pgno,
    //         key ? DKEY(key) : "null"));
    //     mc->mc_flags |= C_INITIALIZED;
    //     mc->mc_flags &= ~C_EOF;
    //
    //     return MDB_SUCCESS;
    // }
    // ```
    /// Leaf page whose key range covers `key` (`None` for an empty tree).
    fn leaf_for(&self, db: &Db, key: &[u8]) -> Result<Option<Page>, String> {
        if db.root == P_INVALID {
            return Ok(None);
        }
        let mut page = self.page(db.root)?;
        for _ in 0..=MAX_DEPTH {
            if page.flags & P_BRANCH != 0 {
                if page.n == 0 {
                    return Err("LMDB: empty branch page (MDB_CORRUPTED)".to_string());
                }
                let (mut i, exact) = page.node_search(key, db.flags)?;
                if i >= page.n {
                    i = page.n - 1;
                } else if !exact {
                    if i == 0 {
                        return Err(
                            "LMDB: branch search ended before the first node (MDB_CORRUPTED)"
                                .to_string(),
                        );
                    }
                    i -= 1;
                }
                let child = page.node(i)?.child;
                page = self.page(child)?;
            } else if page.flags & P_LEAF != 0 && page.flags & P_LEAF2 == 0 {
                return Ok(Some(page));
            } else {
                return Err(format!(
                    "LMDB: index points to a page with flags {:#x}, expected branch or leaf (MDB_CORRUPTED)",
                    page.flags
                ));
            }
        }
        Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)".to_string())
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6032-6047
    // ```c++
    // static int
    // mdb_cursor_set(MDB_cursor *mc, MDB_val *key, MDB_val *data,
    //     MDB_cursor_op op, int *exactp)
    // {
    //     int         rc;
    //     MDB_page    *mp;
    //     MDB_node    *leaf = NULL;
    //     DKBUF;
    //
    //     if (key->mv_size == 0)
    //         return MDB_BAD_VALSIZE;
    //
    //     if (mc->mc_xcursor)
    //         mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    //
    //     /* See if we're already on the right page */
    // ```
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6148-6153
    // ```c++
    // set2:
    //     leaf = mdb_node_search(mc, key, exactp);
    //     if (exactp != NULL && !*exactp) {
    //         /* MDB_SET specified and not an exact match. */
    //         return MDB_NOTFOUND;
    //     }
    // ```
    /// `mdb_cursor_set(MDB_SET)`: leaf page and node of an exact key match.
    fn find(&self, db: &Db, key: &[u8]) -> Result<Option<(Page, Node)>, String> {
        if key.is_empty() {
            return Err(bad_valsize());
        }
        let Some(leaf) = self.leaf_for(db, key)? else {
            return Ok(None);
        };
        let (i, exact) = leaf.node_search(key, db.flags)?;
        if !exact {
            return Ok(None);
        }
        let node = leaf.node(i)?;
        Ok(Some((leaf, node)))
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:9784-9799
    // ```c++
    //     /* Cannot mix named databases with some mainDB flags */
    //     if (txn->mt_dbs[MAIN_DBI].md_flags & (MDB_DUPSORT|MDB_INTEGERKEY))
    //         return (flags & MDB_CREATE) ? MDB_INCOMPATIBLE : MDB_NOTFOUND;
    //
    //     /* Find the DB info */
    //     dbflag = DB_NEW|DB_VALID|DB_USRVALID;
    //     exact = 0;
    //     key.mv_size = len;
    //     key.mv_data = (void *)name;
    //     mdb_cursor_init(&mc, txn, MAIN_DBI, NULL);
    //     rc = mdb_cursor_set(&mc, &key, &data, MDB_SET, &exact);
    //     if (rc == MDB_SUCCESS) {
    //         /* make sure this is actually a DB */
    //         MDB_node *node = NODEPTR(mc.mc_pg[mc.mc_top], mc.mc_ki[mc.mc_top]);
    //         if ((node->mn_flags & (F_DUPDATA|F_SUBDATA)) != F_SUBDATA)
    //             return MDB_INCOMPATIBLE;
    // ```
    /// Look a named database up in the main database (`mdb_dbi_open` of an existing database).
    /// `Ok(None)` when there is no such name (NCBI: `MDB_NOTFOUND`).
    pub fn open_db(&self, name: &[u8]) -> Result<Option<Db>, String> {
        // Cannot mix named databases with some main database flags.
        if self.main.flags & (MDB_DUPSORT | MDB_INTEGERKEY) != 0 {
            return Ok(None);
        }
        let Some((page, node)) = self.find(&self.main, name)? else {
            return Ok(None);
        };
        // Make sure this is actually a database.
        if node.flags & (F_DUPDATA | F_SUBDATA) != F_SUBDATA {
            return Err(
                "MDB_INCOMPATIBLE: Operation and DB incompatible, or DB type changed".to_string(),
            );
        }
        Ok(Some(Db::parse(page.data(&node))?))
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5755-5779
    // ```c++
    // static int
    // mdb_node_read(MDB_cursor *mc, MDB_node *leaf, MDB_val *data)
    // {
    //     MDB_page    *omp;        /* overflow page */
    //     pgno_t         pgno;
    //     int rc;
    //
    //     if (!F_ISSET(leaf->mn_flags, F_BIGDATA)) {
    //         data->mv_size = NODEDSZ(leaf);
    //         data->mv_data = NODEDATA(leaf);
    //         return MDB_SUCCESS;
    //     }
    //
    //     /* Read overflow data.
    //      */
    //     data->mv_size = NODEDSZ(leaf);
    //     memcpy(&pgno, NODEDATA(leaf), sizeof(pgno));
    //     if ((rc = mdb_page_get(mc, pgno, &omp, NULL)) != 0) {
    //         DPRINTF(("read overflow page %"Z"u failed", pgno));
    //         return rc;
    //     }
    //     data->mv_data = METADATA(omp);
    //
    //     return MDB_SUCCESS;
    // }
    // ```
    /// Data of a node: inline, or read from the overflow pages for `F_BIGDATA`.
    fn node_value(&self, page: &Page, node: &Node) -> Result<Vec<u8>, String> {
        if node.flags & F_BIGDATA == 0 {
            return Ok(page.data(node).to_vec());
        }
        let pgno = rd_u64(page.data(node), 0)?;
        if pgno >= self.npages {
            return Err(format!(
                "LMDB: overflow page {pgno} not found (MDB_PAGE_NOTFOUND)"
            ));
        }
        let mut hdr = [0u8; PAGEHDRSZ];
        let base = pgno * self.psize as u64;
        self.src.read_exact_at(&mut hdr, base)?;
        if rd_u16(&hdr, 10)? & P_OVERFLOW == 0 {
            return Err(format!(
                "LMDB: page {pgno} is not an overflow page (MDB_CORRUPTED)"
            ));
        }
        // mp_pages shares the place of mp_lower/mp_upper in the page header.
        let pages = rd_u32(&hdr, 12)? as u64;
        let room = pages
            .saturating_mul(self.psize as u64)
            .saturating_sub(PAGEHDRSZ as u64);
        if pgno.saturating_add(pages) > self.npages || node.dsize as u64 > room {
            return Err(format!(
                "LMDB: overflow value of {} bytes does not fit its {pages} pages (MDB_CORRUPTED)",
                node.dsize
            ));
        }
        let mut v = vec![0u8; node.dsize];
        self.src.read_exact_at(&mut v, base + PAGEHDRSZ as u64)?;
        Ok(v)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7531-7571
    // ```c++
    // static void
    // mdb_xcursor_init1(MDB_cursor *mc, MDB_node *node)
    // {
    //     MDB_xcursor *mx = mc->mc_xcursor;
    //
    //     if (node->mn_flags & F_SUBDATA) {
    //         memcpy(&mx->mx_db, NODEDATA(node), sizeof(MDB_db));
    //         mx->mx_cursor.mc_pg[0] = 0;
    //         mx->mx_cursor.mc_snum = 0;
    //         mx->mx_cursor.mc_top = 0;
    //         mx->mx_cursor.mc_flags = C_SUB;
    //     } else {
    //         MDB_page *fp = NODEDATA(node);
    //         mx->mx_db.md_pad = 0;
    //         mx->mx_db.md_flags = 0;
    //         mx->mx_db.md_depth = 1;
    //         mx->mx_db.md_branch_pages = 0;
    //         mx->mx_db.md_leaf_pages = 1;
    //         mx->mx_db.md_overflow_pages = 0;
    //         mx->mx_db.md_entries = NUMKEYS(fp);
    //         COPY_PGNO(mx->mx_db.md_root, fp->mp_pgno);
    //         mx->mx_cursor.mc_snum = 1;
    //         mx->mx_cursor.mc_top = 0;
    //         mx->mx_cursor.mc_flags = C_INITIALIZED|C_SUB;
    //         mx->mx_cursor.mc_pg[0] = fp;
    //         mx->mx_cursor.mc_ki[0] = 0;
    //         if (mc->mc_db->md_flags & MDB_DUPFIXED) {
    //             mx->mx_db.md_flags = MDB_DUPFIXED;
    //             mx->mx_db.md_pad = fp->mp_pad;
    //             if (mc->mc_db->md_flags & MDB_INTEGERDUP)
    //                 mx->mx_db.md_flags |= MDB_INTEGERKEY;
    //         }
    //     }
    //     DPRINTF(("Sub-db -%u root page %"Z"u", mx->mx_cursor.mc_dbi,
    //         mx->mx_db.md_root));
    //     mx->mx_dbflag = DB_VALID|DB_USRVALID|DB_DUPDATA;
    // #if UINT_MAX < SIZE_MAX
    //     if (mx->mx_dbx.md_cmp == mdb_cmp_int && mx->mx_db.md_pad == sizeof(size_t))
    //         mx->mx_dbx.md_cmp = mdb_cmp_clong;
    // #endif
    // }
    // ```
    /// Shape of the values of one leaf node of `db`.
    fn node_dups(&self, db: &Db, page: &Page, node: &Node) -> Result<Dups, String> {
        if node.flags & F_DUPDATA == 0 {
            return Ok(Dups::One(self.node_value(page, node)?));
        }
        if db.flags & MDB_DUPSORT == 0 || node.flags & F_BIGDATA != 0 {
            return Err("LMDB: duplicate node in a database without MDB_DUPSORT, or with overflow data (MDB_CORRUPTED)".to_string());
        }
        if node.flags & F_SUBDATA != 0 {
            return Ok(Dups::SubDb(Db::parse(page.data(node))?));
        }
        let sub = Page::new(Rc::new(page.data(node).to_vec()))?;
        if sub.flags & P_SUBP == 0 || sub.flags & P_LEAF == 0 {
            return Err(format!("LMDB: duplicate sub-page has flags {:#x}, expected a leaf sub-page (MDB_CORRUPTED)", sub.flags));
        }
        Ok(Dups::SubPage(sub))
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6260-6269
    // ```c++
    //     if (data) {
    //         if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    //             mdb_xcursor_init1(mc, leaf);
    //             rc = mdb_cursor_first(&mc->mc_xcursor->mx_cursor, data, NULL);
    //             if (rc)
    //                 return rc;
    //         } else {
    //             if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    //                 return rc;
    //         }
    // ```
    /// First value in the order of the duplicates (`mdb_cursor_first` of the sub-cursor).
    fn first_dup(&self, dups: &Dups) -> Result<Vec<u8>, String> {
        match dups {
            Dups::One(v) => Ok(v.clone()),
            Dups::SubPage(p) => {
                if p.n == 0 {
                    return Err("LMDB: empty duplicate sub-page (MDB_CORRUPTED)".to_string());
                }
                Ok(p.dup_item(0)?.to_vec())
            }
            Dups::SubDb(sub) => {
                if sub.root == P_INVALID {
                    return Err("LMDB: empty duplicate sub-database (MDB_CORRUPTED)".to_string());
                }
                let mut page = self.page(sub.root)?;
                for _ in 0..=MAX_DEPTH {
                    if page.flags & P_BRANCH != 0 {
                        let child = page.node(0)?.child;
                        page = self.page(child)?;
                    } else if page.flags & P_LEAF != 0 {
                        if page.n == 0 {
                            return Err(
                                "LMDB: empty leaf in a duplicate sub-database (MDB_CORRUPTED)"
                                    .to_string(),
                            );
                        }
                        return Ok(page.dup_item(0)?.to_vec());
                    } else {
                        return Err(format!(
                            "LMDB: sub-database page with flags {:#x} (MDB_CORRUPTED)",
                            page.flags
                        ));
                    }
                }
                Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)".to_string())
            }
        }
    }

    /// Call `f` for every value of `dups`, in LMDB duplicate order (the sub-cursor walk of
    /// `values` and `for_each`, cited there).
    fn each_dup(
        &self,
        dups: &Dups,
        f: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    ) -> Result<(), String> {
        match dups {
            Dups::One(v) => f(v),
            Dups::SubPage(p) => {
                for i in 0..p.n {
                    f(p.dup_item(i)?)?;
                }
                Ok(())
            }
            Dups::SubDb(sub) => {
                if sub.root == P_INVALID {
                    return Ok(());
                }
                self.walk_leaves(sub.root, &mut |leaf: &Page| {
                    for i in 0..leaf.n {
                        f(leaf.dup_item(i)?)?;
                    }
                    Ok(())
                })
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5810-5857
    // ```c++
    // static int
    // mdb_cursor_sibling(MDB_cursor *mc, int move_right)
    // {
    //     int         rc;
    //     MDB_node    *indx;
    //     MDB_page    *mp;
    //
    //     if (mc->mc_snum < 2) {
    //         return MDB_NOTFOUND;        /* root has no siblings */
    //     }
    //
    //     mdb_cursor_pop(mc);
    //     DPRINTF(("parent page is page %"Z"u, index %u",
    //         mc->mc_pg[mc->mc_top]->mp_pgno, mc->mc_ki[mc->mc_top]));
    //
    //     if (move_right ? (mc->mc_ki[mc->mc_top] + 1u >= NUMKEYS(mc->mc_pg[mc->mc_top]))
    //                : (mc->mc_ki[mc->mc_top] == 0)) {
    //         DPRINTF(("no more keys left, moving to %s sibling",
    //             move_right ? "right" : "left"));
    //         if ((rc = mdb_cursor_sibling(mc, move_right)) != MDB_SUCCESS) {
    //             /* undo cursor_pop before returning */
    //             mc->mc_top++;
    //             mc->mc_snum++;
    //             return rc;
    //         }
    //     } else {
    //         if (move_right)
    //             mc->mc_ki[mc->mc_top]++;
    //         else
    //             mc->mc_ki[mc->mc_top]--;
    //         DPRINTF(("just moving to %s index key %u",
    //             move_right ? "right" : "left", mc->mc_ki[mc->mc_top]));
    //     }
    //     mdb_cassert(mc, IS_BRANCH(mc->mc_pg[mc->mc_top]));
    //
    //     indx = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
    //     if ((rc = mdb_page_get(mc, NODEPGNO(indx), &mp, NULL)) != 0) {
    //         /* mc will be inconsistent if caller does mc_snum++ as above */
    //         mc->mc_flags &= ~(C_INITIALIZED|C_EOF);
    //         return rc;
    //     }
    //
    //     mdb_cursor_push(mc, mp);
    //     if (!move_right)
    //         mc->mc_ki[mc->mc_top] = NUMKEYS(mp)-1;
    //
    //     return MDB_SUCCESS;
    // }
    // ```
    /// Leaf pages of the tree under `root` in key order (the page sequence of
    /// `mdb_cursor_first` + `mdb_cursor_sibling`). A tree that visits more pages than the
    /// database has is a cycle.
    fn walk_leaves(
        &self,
        root: u64,
        f: &mut dyn FnMut(&Page) -> Result<(), String>,
    ) -> Result<(), String> {
        let mut budget = self.npages;
        self.walk_rec(root, 0, &mut budget, f)
    }

    fn walk_rec(
        &self,
        pgno: u64,
        depth: usize,
        budget: &mut u64,
        f: &mut dyn FnMut(&Page) -> Result<(), String>,
    ) -> Result<(), String> {
        if depth > MAX_DEPTH {
            return Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)".to_string());
        }
        if *budget == 0 {
            return Err(
                "LMDB: tree visits more pages than the database has (MDB_CORRUPTED)".to_string(),
            );
        }
        *budget -= 1;
        let page = self.page(pgno)?;
        if page.flags & P_BRANCH != 0 {
            for i in 0..page.n {
                let child = page.node(i)?.child;
                self.walk_rec(child, depth + 1, budget, f)?;
            }
            Ok(())
        } else if page.flags & P_LEAF != 0 {
            f(&page)
        } else {
            Err(format!(
                "LMDB: page {pgno} has flags {:#x}, expected branch or leaf (MDB_CORRUPTED)",
                page.flags
            ))
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6333-6359
    // ```c++
    //     case MDB_GET_CURRENT:
    //         if (!(mc->mc_flags & C_INITIALIZED)) {
    //             rc = EINVAL;
    //         } else {
    //             MDB_page *mp = mc->mc_pg[mc->mc_top];
    //             int nkeys = NUMKEYS(mp);
    //             if (!nkeys || mc->mc_ki[mc->mc_top] >= nkeys) {
    //                 mc->mc_ki[mc->mc_top] = nkeys;
    //                 rc = MDB_NOTFOUND;
    //                 break;
    //             }
    //             rc = MDB_SUCCESS;
    //             if (IS_LEAF2(mp)) {
    //                 key->mv_size = mc->mc_db->md_pad;
    //                 key->mv_data = LEAF2KEY(mp, mc->mc_ki[mc->mc_top], key->mv_size);
    //             } else {
    //                 MDB_node *leaf = NODEPTR(mp, mc->mc_ki[mc->mc_top]);
    //                 MDB_GET_KEY(leaf, key);
    //                 if (data) {
    //                     if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    //                         rc = mdb_cursor_get(&mc->mc_xcursor->mx_cursor, data, NULL, MDB_GET_CURRENT);
    //                     } else {
    //                         rc = mdb_node_read(mc, leaf, data);
    //                     }
    //                 }
    //             }
    //         }
    // ```
    /// `MDB_SET` + `MDB_GET_CURRENT`: the first duplicate of `key` (the only value of a plain
    /// key). Works for sub-page and sub-database keys. `Ok(None)`: no such key.
    pub fn first_value(&self, db: &Db, key: &[u8]) -> Result<Option<Vec<u8>>, String> {
        let Some((page, node)) = self.find(db, key)? else {
            return Ok(None);
        };
        let dups = self.node_dups(db, &page, &node)?;
        Ok(Some(self.first_dup(&dups)?))
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5860-5896
    // ```c++
    // static int
    // mdb_cursor_next(MDB_cursor *mc, MDB_val *key, MDB_val *data, MDB_cursor_op op)
    // {
    //     MDB_page    *mp;
    //     MDB_node    *leaf;
    //     int rc;
    //
    //     if ((mc->mc_flags & C_DEL && op == MDB_NEXT_DUP))
    //         return MDB_NOTFOUND;
    //
    //     if (!(mc->mc_flags & C_INITIALIZED))
    //         return mdb_cursor_first(mc, key, data);
    //
    //     mp = mc->mc_pg[mc->mc_top];
    //
    //     if (mc->mc_flags & C_EOF) {
    //         if (mc->mc_ki[mc->mc_top] >= NUMKEYS(mp)-1)
    //             return MDB_NOTFOUND;
    //         mc->mc_flags ^= C_EOF;
    //     }
    //
    //     if (mc->mc_db->md_flags & MDB_DUPSORT) {
    //         leaf = NODEPTR(mp, mc->mc_ki[mc->mc_top]);
    //         if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    //             if (op == MDB_NEXT || op == MDB_NEXT_DUP) {
    //                 rc = mdb_cursor_next(&mc->mc_xcursor->mx_cursor, data, NULL, MDB_NEXT);
    //                 if (op != MDB_NEXT || rc != MDB_NOTFOUND) {
    //                     if (rc == MDB_SUCCESS)
    //                         MDB_GET_KEY(leaf, key);
    //                     return rc;
    //                 }
    //             }
    //         } else {
    //             mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    //             if (op == MDB_NEXT_DUP)
    //                 return MDB_NOTFOUND;
    //         }
    // ```
    /// `MDB_SET` + `MDB_GET_CURRENT` + `MDB_NEXT_DUP` until the end: every value of `key` in
    /// LMDB order (for `MDB_DUPFIXED` without `MDB_INTEGERDUP`: `memcmp` of the value bytes).
    /// Empty when the key is absent.
    pub fn values(&self, db: &Db, key: &[u8]) -> Result<Vec<Vec<u8>>, String> {
        let Some((page, node)) = self.find(db, key)? else {
            return Ok(Vec::new());
        };
        let dups = self.node_dups(db, &page, &node)?;
        let mut out = Vec::new();
        self.each_dup(&dups, &mut |v| {
            out.push(v.to_vec());
            Ok(())
        })?;
        Ok(out)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6232-6273
    // ```c++
    // static int
    // mdb_cursor_first(MDB_cursor *mc, MDB_val *key, MDB_val *data)
    // {
    //     int         rc;
    //     MDB_node    *leaf;
    //
    //     if (mc->mc_xcursor)
    //         mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    //
    //     if (!(mc->mc_flags & C_INITIALIZED) || mc->mc_top) {
    //         rc = mdb_page_search(mc, NULL, MDB_PS_FIRST);
    //         if (rc != MDB_SUCCESS)
    //             return rc;
    //     }
    //     mdb_cassert(mc, IS_LEAF(mc->mc_pg[mc->mc_top]));
    //
    //     leaf = NODEPTR(mc->mc_pg[mc->mc_top], 0);
    //     mc->mc_flags |= C_INITIALIZED;
    //     mc->mc_flags &= ~C_EOF;
    //
    //     mc->mc_ki[mc->mc_top] = 0;
    //
    //     if (IS_LEAF2(mc->mc_pg[mc->mc_top])) {
    //         key->mv_size = mc->mc_db->md_pad;
    //         key->mv_data = LEAF2KEY(mc->mc_pg[mc->mc_top], 0, key->mv_size);
    //         return MDB_SUCCESS;
    //     }
    //
    //     if (data) {
    //         if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    //             mdb_xcursor_init1(mc, leaf);
    //             rc = mdb_cursor_first(&mc->mc_xcursor->mx_cursor, data, NULL);
    //             if (rc)
    //                 return rc;
    //         } else {
    //             if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    //                 return rc;
    //         }
    //     }
    //     MDB_GET_KEY(leaf, key);
    //     return MDB_SUCCESS;
    // }
    // ```
    /// `MDB_FIRST` then `MDB_NEXT` until the end: every `(key, value)` pair in tree order (a key
    /// with several duplicates appears once per duplicate). Used for `GetDBTaxIds`.
    pub fn for_each(
        &self,
        db: &Db,
        mut f: impl FnMut(&[u8], &[u8]) -> Result<(), String>,
    ) -> Result<(), String> {
        if db.root == P_INVALID {
            return Ok(());
        }
        self.walk_leaves(db.root, &mut |leaf: &Page| {
            if leaf.flags & P_LEAF2 != 0 {
                return Err("LMDB: LEAF2 page in a main tree (MDB_CORRUPTED)".to_string());
            }
            for i in 0..leaf.n {
                let node = leaf.node(i)?;
                let key = leaf.key(&node);
                let dups = self.node_dups(db, leaf, &node)?;
                self.each_dup(&dups, &mut |v| f(key, v))?;
            }
            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
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
                    let at = if rng.next() % 2 == 0 {
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
}
