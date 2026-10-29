//! Byte layout of LMDB pages, nodes and records (LMDB 0.9.24 as built into NCBI BLAST+,
//! 64-bit little-endian, `MDB_DEVEL` off). Accessors over a page buffer: a page starts at
//! offset 0 of the slice (a sub-page inside a node is a slice that starts at the node data).

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:805-848
// ```c
// typedef struct MDB_page {
//     union { pgno_t p_pgno; struct MDB_page *p_next; } mp_p;
//     uint16_t mp_pad;            /**< key size if this is a LEAF2 page */
// #define P_BRANCH 0x01 ... #define P_SUBP 0x40
//     uint16_t mp_flags;
//     union { struct { indx_t pb_lower; indx_t pb_upper; } pb; uint32_t pb_pages; } mp_pb;
//     indx_t mp_ptrs[1];          /**< dynamic size */
// } MDB_page;
// #define PAGEHDRSZ ((unsigned) offsetof(MDB_page, mp_ptrs))
// #define NUMKEYS(p) (((p)->mp_lower - (PAGEHDRSZ-PAGEBASE)) >> 1)
// #define SIZELEFT(p) (indx_t)((p)->mp_upper - (p)->mp_lower)
// ```

/// `PAGEHDRSZ`.
pub(super) const PAGEHDRSZ: usize = 16;
/// `NODESIZE` (`offsetof(MDB_node, mn_data)`).
pub(super) const NODESIZE: usize = 8;

pub(super) const P_BRANCH: u16 = 0x01;
pub(super) const P_LEAF: u16 = 0x02;
pub(super) const P_OVERFLOW: u16 = 0x04;
pub(super) const P_META: u16 = 0x08;
pub(super) const P_DIRTY: u16 = 0x10;
pub(super) const P_LEAF2: u16 = 0x20;
pub(super) const P_SUBP: u16 = 0x40;

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:906-914
// ```c
// #define F_BIGDATA	 0x01			/**< data put on overflow page */
// #define F_SUBDATA	 0x02			/**< data is a sub-database */
// #define F_DUPDATA	 0x04			/**< data has duplicates */
// ```
pub(super) const F_BIGDATA: u16 = 0x01;
pub(super) const F_SUBDATA: u16 = 0x02;
pub(super) const F_DUPDATA: u16 = 0x04;

#[inline]
pub(super) fn get_u16(b: &[u8], off: usize) -> u16 {
    u16::from_le_bytes([b[off], b[off + 1]])
}

#[inline]
pub(super) fn put_u16(b: &mut [u8], off: usize, v: u16) {
    b[off..off + 2].copy_from_slice(&v.to_le_bytes());
}

#[inline]
pub(super) fn get_u32(b: &[u8], off: usize) -> u32 {
    u32::from_le_bytes(b[off..off + 4].try_into().unwrap())
}

#[inline]
pub(super) fn put_u32(b: &mut [u8], off: usize, v: u32) {
    b[off..off + 4].copy_from_slice(&v.to_le_bytes());
}

#[inline]
pub(super) fn get_u64(b: &[u8], off: usize) -> u64 {
    u64::from_le_bytes(b[off..off + 8].try_into().unwrap())
}

#[inline]
pub(super) fn put_u64(b: &mut [u8], off: usize, v: u64) {
    b[off..off + 8].copy_from_slice(&v.to_le_bytes());
}

// Page header fields.
#[inline]
pub(super) fn pgno(p: &[u8]) -> u64 {
    get_u64(p, 0)
}
#[inline]
pub(super) fn set_pgno(p: &mut [u8], v: u64) {
    put_u64(p, 0, v)
}
#[inline]
pub(super) fn pad(p: &[u8]) -> u16 {
    get_u16(p, 8)
}
#[inline]
pub(super) fn set_pad(p: &mut [u8], v: u16) {
    put_u16(p, 8, v)
}
#[inline]
pub(super) fn flags(p: &[u8]) -> u16 {
    get_u16(p, 10)
}
#[inline]
pub(super) fn set_flags(p: &mut [u8], v: u16) {
    put_u16(p, 10, v)
}
#[inline]
pub(super) fn lower(p: &[u8]) -> u16 {
    get_u16(p, 12)
}
#[inline]
pub(super) fn set_lower(p: &mut [u8], v: u16) {
    put_u16(p, 12, v)
}
#[inline]
pub(super) fn upper(p: &[u8]) -> u16 {
    get_u16(p, 14)
}
#[inline]
pub(super) fn set_upper(p: &mut [u8], v: u16) {
    put_u16(p, 14, v)
}
/// `mp_pages` of an overflow page (the union with `mp_lower`/`mp_upper`).
#[inline]
pub(super) fn pages(p: &[u8]) -> u32 {
    get_u32(p, 12)
}
#[inline]
pub(super) fn set_pages(p: &mut [u8], v: u32) {
    put_u32(p, 12, v)
}
/// `NUMKEYS(p)`.
#[inline]
pub(super) fn numkeys(p: &[u8]) -> usize {
    (lower(p) as usize).wrapping_sub(PAGEHDRSZ) >> 1
}
/// `SIZELEFT(p)`: `(indx_t)(mp_upper - mp_lower)`.
#[inline]
pub(super) fn sizeleft(p: &[u8]) -> u16 {
    upper(p).wrapping_sub(lower(p))
}
#[inline]
pub(super) fn is_leaf(p: &[u8]) -> bool {
    flags(p) & P_LEAF == P_LEAF
}
#[inline]
pub(super) fn is_leaf2(p: &[u8]) -> bool {
    flags(p) & P_LEAF2 == P_LEAF2
}
#[inline]
pub(super) fn is_branch(p: &[u8]) -> bool {
    flags(p) & P_BRANCH == P_BRANCH
}
#[inline]
pub(super) fn is_overflow(p: &[u8]) -> bool {
    flags(p) & P_OVERFLOW == P_OVERFLOW
}
#[inline]
pub(super) fn is_subp(p: &[u8]) -> bool {
    flags(p) & P_SUBP == P_SUBP
}
/// `mp_ptrs[i]`.
#[inline]
pub(super) fn ptr(p: &[u8], i: usize) -> u16 {
    get_u16(p, PAGEHDRSZ + 2 * i)
}
#[inline]
pub(super) fn set_ptr(p: &mut [u8], i: usize, v: u16) {
    put_u16(p, PAGEHDRSZ + 2 * i, v)
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:898-990
// ```c
// typedef struct MDB_node {
//     unsigned short mn_lo, mn_hi;    /**< part of data size or pgno */
//     unsigned short mn_flags;
//     unsigned short mn_ksize;        /**< key size */
//     char mn_data[1];                /**< key and data are appended here */
// } MDB_node;
// #define NODEPTR(p, i)	 ((MDB_node *)((char *)(p) + (p)->mp_ptrs[i] + PAGEBASE))
// #define NODEKEY(node)	 (void *)((node)->mn_data)
// #define NODEDATA(node)	 (void *)((char *)(node)->mn_data + (node)->mn_ksize)
// #define NODEPGNO(node) ((node)->mn_lo | ((pgno_t) (node)->mn_hi << 16) | \
//	 (PGNO_TOPWORD ? ((pgno_t) (node)->mn_flags << PGNO_TOPWORD) : 0))
// #define NODEDSZ(node)	 ((node)->mn_lo | ((unsigned)(node)->mn_hi << 16))
// #define LEAF2KEY(p, i, ks)	((char *)(p) + PAGEHDRSZ + ((i)*(ks)))
// ```
// Node accessors take the offset of the node in the page (`mp_ptrs[i]`).

/// Offset of node `i` (`NODEPTR(p, i)`).
#[inline]
pub(super) fn nodeptr(p: &[u8], i: usize) -> usize {
    ptr(p, i) as usize
}
#[inline]
pub(super) fn node_flags(p: &[u8], n: usize) -> u16 {
    get_u16(p, n + 4)
}
#[inline]
pub(super) fn set_node_flags(p: &mut [u8], n: usize, v: u16) {
    put_u16(p, n + 4, v)
}
#[inline]
pub(super) fn node_ksize(p: &[u8], n: usize) -> usize {
    get_u16(p, n + 6) as usize
}
/// Offset of the key of the node at `n` (`NODEKEY`).
#[inline]
pub(super) fn node_key(n: usize) -> usize {
    n + NODESIZE
}
/// Offset of the data of the node at `n` (`NODEDATA`).
#[inline]
pub(super) fn node_data(p: &[u8], n: usize) -> usize {
    n + NODESIZE + node_ksize(p, n)
}
/// `NODEDSZ`.
#[inline]
pub(super) fn node_dsz(p: &[u8], n: usize) -> usize {
    get_u16(p, n) as usize | ((get_u16(p, n + 2) as usize) << 16)
}
/// `SETDSZ`.
#[inline]
pub(super) fn set_node_dsz(p: &mut [u8], n: usize, size: usize) {
    put_u16(p, n, (size & 0xffff) as u16);
    put_u16(p, n + 2, (size >> 16) as u16);
}
/// `NODEPGNO`.
#[inline]
pub(super) fn node_pgno(p: &[u8], n: usize) -> u64 {
    get_u16(p, n) as u64 | ((get_u16(p, n + 2) as u64) << 16) | ((get_u16(p, n + 4) as u64) << 32)
}
/// `SETPGNO`.
#[inline]
pub(super) fn set_node_pgno(p: &mut [u8], n: usize, pgno: u64) {
    put_u16(p, n, (pgno & 0xffff) as u16);
    put_u16(p, n + 2, ((pgno >> 16) & 0xffff) as u16);
    put_u16(p, n + 4, (pgno >> 32) as u16);
}
/// Offset of LEAF2 key `i` (`LEAF2KEY`).
#[inline]
pub(super) fn leaf2key(i: usize, ks: usize) -> usize {
    PAGEHDRSZ + i * ks
}

/// `EVEN(n)`: `((n) + 1U) & -2`.
#[inline]
pub(super) fn even(n: usize) -> usize {
    (n + 1) & !1
}

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
/// `MDB_db`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct MdbDb {
    pub pad: u32,
    pub flags: u16,
    pub depth: u16,
    pub branch_pages: u64,
    pub leaf_pages: u64,
    pub overflow_pages: u64,
    pub entries: u64,
    pub root: u64,
}

/// `sizeof(MDB_db)`.
pub(super) const MDB_DB_SIZE: usize = 48;

impl MdbDb {
    pub(super) fn to_bytes(self) -> [u8; MDB_DB_SIZE] {
        let mut b = [0u8; MDB_DB_SIZE];
        put_u32(&mut b, 0, self.pad);
        put_u16(&mut b, 4, self.flags);
        put_u16(&mut b, 6, self.depth);
        put_u64(&mut b, 8, self.branch_pages);
        put_u64(&mut b, 16, self.leaf_pages);
        put_u64(&mut b, 24, self.overflow_pages);
        put_u64(&mut b, 32, self.entries);
        put_u64(&mut b, 40, self.root);
        b
    }

    pub(super) fn from_bytes(b: &[u8]) -> Self {
        MdbDb {
            pad: get_u32(b, 0),
            flags: get_u16(b, 4),
            depth: get_u16(b, 6),
            branch_pages: get_u64(b, 8),
            leaf_pages: get_u64(b, 16),
            overflow_pages: get_u64(b, 24),
            entries: get_u64(b, 32),
            root: get_u64(b, 40),
        }
    }
}

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
/// `MDB_meta`, kept in memory as it is in the file (the map in NCBI's LMDB).
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Meta {
    pub magic: u32,
    pub version: u32,
    pub address: u64,
    pub mapsize: u64,
    pub dbs: [MdbDb; 2],
    pub last_pg: u64,
    pub txnid: u64,
}

/// `sizeof(MDB_meta)`.
pub(super) const META_SIZE: usize = 136;
/// `offsetof(MDB_meta, mm_mapsize)`.
pub(super) const META_MAPSIZE_OFF: usize = 16;

impl Meta {
    pub(super) fn to_bytes(self) -> [u8; META_SIZE] {
        let mut b = [0u8; META_SIZE];
        put_u32(&mut b, 0, self.magic);
        put_u32(&mut b, 4, self.version);
        put_u64(&mut b, 8, self.address);
        put_u64(&mut b, 16, self.mapsize);
        b[24..72].copy_from_slice(&self.dbs[0].to_bytes());
        b[72..120].copy_from_slice(&self.dbs[1].to_bytes());
        put_u64(&mut b, 120, self.last_pg);
        put_u64(&mut b, 128, self.txnid);
        b
    }
}
