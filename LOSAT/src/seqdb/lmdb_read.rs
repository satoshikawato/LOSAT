// Derived from LMDB 0.9.24 (mdb.c, lmdb.h) as bundled in NCBI BLAST+ 2.17.0,
// c++/src/util/lmdb/, which carries this notice:
//
// Copyright 2011-2019 Howard Chu, Symas Corp.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted only as authorized by the OpenLDAP
// Public License.
//
// A copy of this license is available in the file LICENSE in the
// top-level directory of the distribution or, alternatively, at
// <http://www.OpenLDAP.org/license.html>.
//
// This code is derived from btree.c written by Martin Hedenfalk.
//
// Copyright (c) 2009, 2010 Martin Hedenfalk <martin@bzero.se>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//
// (In LOSAT, the OpenLDAP Public License 2.8 is in THIRD_PARTY_NOTICES, section 2 (LMDB).)

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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp:60-67
/// Why a reader over [`ReadAt`] (this LMDB reader, `tax_files`'s lookup files, `tax_sqlite`) gives
/// up: the file's bytes are not what it reads (`Format`, LOSAT's text of a damaged or foreign
/// file, which the callers turn into their own error or NCBI's text), or the system does not give
/// them (`Unreadable`). The second is a typed boundary so that a read that fails after the file
/// was opened is not a format refusal (line WQA D18c, 2026-11-21): NCBI maps the LMDB and lookup
/// files, where such a failure is a SIGBUS, so it stays [`super::SeqDbError::Read`], LOSAT's own
/// "cannot be read"; it reads `taxonomy4blast.sqlite3` through SQLite, where it is SQLite's read
/// error, which `tax_sqlite`'s callers write as NCBI does (`tax_sqlite::QueryReadFailure`; review
/// 1 of line WQA D18b, L3).
#[derive(Debug)]
pub enum ReadError {
    /// The bytes are not those of the expected format, or a range passes the end of the source.
    Format(String),
    /// The system failed to give bytes of a file after its first byte ([`super::SeqDbError::Read`]:
    /// where NCBI reads its mapping, a SIGBUS, which `report::db_subject::LookupFailure::of_error`
    /// takes as NCBI's crash, line WQA D18 round 3; of `taxonomy4blast.sqlite3`, SQLite's read
    /// error, `tax_sqlite::QueryReadFailure`).
    Unreadable(super::SeqDbError),
}

impl ReadError {
    /// The database error of this failure: the `Unreadable` error as it is, the format message
    /// through `format` (each caller's own error or NCBI's text for it).
    pub fn into_seqdb(self, format: impl FnOnce(String) -> super::SeqDbError) -> super::SeqDbError {
        match self {
            ReadError::Format(message) => format(message),
            ReadError::Unreadable(error) => error,
        }
    }

    /// This failure with `prefix` before a format message (an `Unreadable` error keeps its text).
    pub fn prefixed(self, prefix: &str) -> Self {
        match self {
            ReadError::Format(message) => ReadError::Format(format!("{prefix}{message}")),
            unreadable => unreadable,
        }
    }
}

impl From<String> for ReadError {
    fn from(message: String) -> Self {
        ReadError::Format(message)
    }
}

impl From<&str> for ReadError {
    fn from(message: &str) -> Self {
        ReadError::Format(message.to_string())
    }
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Format(message) => f.write_str(message),
            ReadError::Unreadable(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReadError {}

/// Positioned reads over a database file ([`super::file::DbFile`]) or a buffer.
pub trait ReadAt {
    /// Total length in bytes.
    fn len(&self) -> u64;
    /// Fill `buf` from `offset`; `Err`: [`ReadError::Format`] when the range is not inside the
    /// source, [`ReadError::Unreadable`] when the system does not give the bytes.
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError>;
}

impl ReadAt for [u8] {
    fn len(&self) -> u64 {
        <[u8]>::len(self) as u64
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
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
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        self.as_slice().read_exact_at(buf, offset)
    }
}

impl ReadAt for super::file::DbFile {
    fn len(&self) -> u64 {
        super::file::DbFile::len(self)
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        super::file::DbFile::read_exact_at(self, buf, offset).map_err(|error| match error {
            super::SeqDbError::Read { .. } => ReadError::Unreadable(error),
            other => ReadError::Format(other.to_string()),
        })
    }
}

/// Plain files, for the tests (the library reads through [`super::file::DbFile`]).
#[cfg(test)]
impl ReadAt for std::fs::File {
    fn len(&self) -> u64 {
        self.metadata().map(|m| m.len()).unwrap_or(0)
    }
    #[cfg(unix)]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        std::os::unix::fs::FileExt::read_exact_at(self, buf, offset)
            .map_err(|e| format!("read of {} bytes at offset {}: {e}", buf.len(), offset).into())
    }
    #[cfg(not(any(unix, windows)))]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        use std::io::{Read, Seek, SeekFrom};
        let mut f = self;
        f.seek(SeekFrom::Start(offset))
            .and_then(|_| f.read_exact(buf))
            .map_err(|e| format!("read of {} bytes at offset {}: {e}", buf.len(), offset).into())
    }
    #[cfg(windows)]
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
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
                )
                .into());
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

fn slice(b: &[u8], off: usize, n: usize) -> Result<&[u8], ReadError> {
    off.checked_add(n)
        .and_then(|end| b.get(off..end))
        .ok_or_else(|| {
            ReadError::Format(format!(
                "LMDB: {n} bytes at offset {off} are outside a block of {} bytes",
                b.len()
            ))
        })
}

fn rd_u16(b: &[u8], off: usize) -> Result<u16, ReadError> {
    Ok(u16::from_le_bytes(
        slice(b, off, 2)?
            .try_into()
            .map_err(|_| "LMDB: short read")?,
    ))
}

fn rd_u32(b: &[u8], off: usize) -> Result<u32, ReadError> {
    Ok(u32::from_le_bytes(
        slice(b, off, 4)?
            .try_into()
            .map_err(|_| "LMDB: short read")?,
    ))
}

fn rd_u64(b: &[u8], off: usize) -> Result<u64, ReadError> {
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
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1005-1014
    fn parse(raw: &[u8]) -> Result<Db, ReadError> {
        if raw.len() < DB_RECORD_SZ {
            return Err(format!(
                "LMDB: database record of {} bytes, expected {DB_RECORD_SZ}",
                raw.len()
            )
            .into());
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

    /// True when the tree has no page (`md_root` is `P_INVALID`).
    #[cfg(test)]
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

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1036-1054
/// Read and validate the meta page that starts at `off` (header 16 bytes, then the meta).
fn read_meta<R: ReadAt + ?Sized>(src: &R, off: u64) -> Result<Meta, ReadError> {
    let end = off
        .checked_add(META_BYTES as u64)
        .filter(|&e| e <= src.len());
    if end.is_none() {
        return Err(format!(
            "LMDB: file of {} bytes is too small for a meta page at offset {off}",
            src.len()
        )
        .into());
    }
    let mut b = [0u8; META_BYTES];
    src.read_exact_at(&mut b, off)?;
    // The checks of mdb_env_read_header: P_META flag, magic, version.
    if rd_u16(&b, 10)? & P_META == 0 {
        return Err(format!("LMDB: page at offset {off} is not a meta page (MDB_INVALID)").into());
    }
    if rd_u32(&b, 16)? != MDB_MAGIC {
        return Err("LMDB: meta page has invalid magic (MDB_INVALID)"
            .to_string()
            .into());
    }
    let version = rd_u32(&b, 20)?;
    if version != MDB_DATA_VERSION {
        return Err(format!(
            "LMDB: database is version {version}, expected version {MDB_DATA_VERSION} (MDB_VERSION_MISMATCH)"
        ).into());
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
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:846-852
    fn new(buf: Rc<Vec<u8>>) -> Result<Page, ReadError> {
        let pad = rd_u16(&buf, 8)?;
        let flags = rd_u16(&buf, 10)?;
        let lower = rd_u16(&buf, 12)? as usize;
        if lower < PAGEHDRSZ || lower > buf.len() {
            return Err(format!(
                "LMDB: page header lower={lower} is outside a page of {} bytes",
                buf.len()
            )
            .into());
        }
        Ok(Page {
            buf,
            flags,
            pad,
            n: (lower - PAGEHDRSZ) >> 1,
        })
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:898-960
    fn node(&self, i: usize) -> Result<Node, ReadError> {
        if self.flags & P_LEAF2 != 0 || self.flags & (P_BRANCH | P_LEAF) == 0 {
            return Err(format!("LMDB: page with flags {:#x} has no nodes", self.flags).into());
        }
        if i >= self.n {
            return Err(format!("LMDB: node {i} is outside a page of {} nodes", self.n).into());
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

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:995-995
    fn leaf2_key(&self, i: usize) -> Result<&[u8], ReadError> {
        let ks = self.pad as usize;
        if ks == 0 {
            return Err("LMDB: LEAF2 page with key size 0".to_string().into());
        }
        if i >= self.n {
            return Err(format!("LMDB: LEAF2 key {i} is outside a page of {} keys", self.n).into());
        }
        slice(&self.buf, PAGEHDRSZ + i * ks, ks)
    }

    /// Key of entry `i` of a leaf page of a sub-database or sub-page: LEAF2 keys are the
    /// fixed-size values, in other leaves the value of a duplicate is the node key.
    fn dup_item(&self, i: usize) -> Result<&[u8], ReadError> {
        if self.flags & P_LEAF2 != 0 {
            self.leaf2_key(i)
        } else {
            let n = self.node(i)?;
            Ok(self.key(&n))
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5276-5365
    /// `mdb_node_search`: index of the first entry >= `key` and whether it is equal.
    fn node_search(&self, key: &[u8], db_flags: u16) -> Result<(usize, bool), ReadError> {
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

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5197-5217
fn cmp_cint(a: &[u8], b: &[u8]) -> Ordering {
    // Keys of one database have one size (4 bytes here); NCBI reads `a->mv_size` bytes of both.
    if a.len() != b.len() {
        return a.len().cmp(&b.len());
    }
    a.iter().rev().cmp(b.iter().rev())
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5225-5241
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
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3682-3738
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3943-3949
    /// Open an environment: validate both meta pages (magic, version), take the page size from
    /// `FREE_DBI.md_pad`, use the meta page with the larger transaction id (tie: meta 0).
    pub fn open(src: &'r R) -> Result<Self, ReadError> {
        let m0 = read_meta(src, 0)?;
        let psize = m0.free.pad as usize;
        if !(512..=32768).contains(&psize) || !psize.is_power_of_two() {
            return Err(format!("LMDB: invalid page size {psize} in the meta page").into());
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

    /// Number of pages that can be read: one more than `mm_last_pg` of the meta page in use, at
    /// most the whole pages of the file.
    pub fn num_pages(&self) -> u64 {
        self.npages
    }

    /// Transaction id of the meta page in use.
    #[cfg(test)]
    pub fn txnid(&self) -> u64 {
        self.txnid
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5465-5479
    fn page(&self, pgno: u64) -> Result<Page, ReadError> {
        if pgno >= self.npages {
            return Err(format!(
                "LMDB: page {pgno} not found (the database has {} pages) (MDB_PAGE_NOTFOUND)",
                self.npages
            )
            .into());
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

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5484-5562
    /// Leaf page whose key range covers `key` (`None` for an empty tree).
    fn leaf_for(&self, db: &Db, key: &[u8]) -> Result<Option<Page>, ReadError> {
        self.descend(db, key, |_, _| {})
    }

    /// `leaf_for`, handing each branch page on the way down to `on_branch` with the index of
    /// the child taken (the cursor's stack, `mc_pg` and `mc_ki`).
    fn descend(
        &self,
        db: &Db,
        key: &[u8],
        mut on_branch: impl FnMut(Page, usize),
    ) -> Result<Option<Page>, ReadError> {
        if db.root == P_INVALID {
            return Ok(None);
        }
        let mut page = self.page(db.root)?;
        for _ in 0..=MAX_DEPTH {
            if page.flags & P_BRANCH != 0 {
                if page.n == 0 {
                    return Err("LMDB: empty branch page (MDB_CORRUPTED)".to_string().into());
                }
                let (mut i, exact) = page.node_search(key, db.flags)?;
                if i >= page.n {
                    i = page.n - 1;
                } else if !exact {
                    if i == 0 {
                        return Err(
                            "LMDB: branch search ended before the first node (MDB_CORRUPTED)"
                                .to_string()
                                .into(),
                        );
                    }
                    i -= 1;
                }
                let child = page.node(i)?.child;
                on_branch(page, i);
                page = self.page(child)?;
            } else if page.flags & P_LEAF != 0 && page.flags & P_LEAF2 == 0 {
                return Ok(Some(page));
            } else {
                return Err(format!(
                    "LMDB: index points to a page with flags {:#x}, expected branch or leaf (MDB_CORRUPTED)",
                    page.flags
                ).into());
            }
        }
        Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)"
            .to_string()
            .into())
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6032-6047
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6148-6153
    /// `mdb_cursor_set` in mode `MDB_SET`: leaf page and node of an exact key match.
    fn find(&self, db: &Db, key: &[u8]) -> Result<Option<(Page, Node)>, ReadError> {
        if key.is_empty() {
            return Err(bad_valsize().into());
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

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6148-6160
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5821-5845
    /// `mdb_cursor_set` in mode `MDB_SET_RANGE`: the first key of `db` that is not before `key` in
    /// its order (`None`: every key is before it). Past the last node of the leaf that covers
    /// `key`, the first key of the next leaf (`mdb_cursor_sibling`: the deepest branch with a node
    /// to the right of the one taken, then the first node of each page down).
    pub fn first_key_at_least(&self, db: &Db, key: &[u8]) -> Result<Option<Vec<u8>>, ReadError> {
        if key.is_empty() {
            return Err(bad_valsize().into());
        }
        let mut path = Vec::new();
        let Some(leaf) = self.descend(db, key, |page, i| path.push((page, i)))? else {
            return Ok(None);
        };
        let (i, _) = leaf.node_search(key, db.flags)?;
        if i < leaf.n {
            return Ok(Some(leaf.key(&leaf.node(i)?).to_vec()));
        }
        while let Some((branch, i)) = path.pop() {
            if i + 1 >= branch.n {
                continue;
            }
            let mut page = self.page(branch.node(i + 1)?.child)?;
            for _ in path.len()..=MAX_DEPTH {
                if page.flags & P_BRANCH != 0 {
                    if page.n == 0 {
                        return Err("LMDB: empty branch page (MDB_CORRUPTED)".to_string().into());
                    }
                    page = self.page(page.node(0)?.child)?;
                } else if page.flags & P_LEAF != 0 && page.flags & P_LEAF2 == 0 {
                    if page.n == 0 {
                        return Err("LMDB: empty leaf page past the root (MDB_CORRUPTED)"
                            .to_string()
                            .into());
                    }
                    return Ok(Some(page.key(&page.node(0)?).to_vec()));
                } else {
                    return Err(format!(
                        "LMDB: index points to a page with flags {:#x}, expected branch or leaf (MDB_CORRUPTED)",
                        page.flags
                    ).into());
                }
            }
            return Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)"
                .to_string()
                .into());
        }
        Ok(None)
    }

    /// Whether `db` has a key that starts with `prefix` (`MDB_SET_RANGE` at `prefix`).
    pub fn has_key_starting_with(&self, db: &Db, prefix: &[u8]) -> Result<bool, ReadError> {
        Ok(self
            .first_key_at_least(db, prefix)?
            .is_some_and(|key| key.starts_with(prefix)))
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:9784-9799
    /// Look a named database up in the main database (`mdb_dbi_open` of an existing database).
    /// `Ok(None)` when there is no such name (NCBI: `MDB_NOTFOUND`).
    pub fn open_db(&self, name: &[u8]) -> Result<Option<Db>, ReadError> {
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
                "MDB_INCOMPATIBLE: Operation and DB incompatible, or DB type changed"
                    .to_string()
                    .into(),
            );
        }
        Ok(Some(Db::parse(page.data(&node))?))
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5755-5779
    /// Data of a node: inline, or read from the overflow pages for `F_BIGDATA`.
    fn node_value(&self, page: &Page, node: &Node) -> Result<Vec<u8>, ReadError> {
        if node.flags & F_BIGDATA == 0 {
            return Ok(page.data(node).to_vec());
        }
        let pgno = rd_u64(page.data(node), 0)?;
        if pgno >= self.npages {
            return Err(format!("LMDB: overflow page {pgno} not found (MDB_PAGE_NOTFOUND)").into());
        }
        let mut hdr = [0u8; PAGEHDRSZ];
        let base = pgno * self.psize as u64;
        self.src.read_exact_at(&mut hdr, base)?;
        if rd_u16(&hdr, 10)? & P_OVERFLOW == 0 {
            return Err(
                format!("LMDB: page {pgno} is not an overflow page (MDB_CORRUPTED)").into(),
            );
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
            )
            .into());
        }
        let mut v = vec![0u8; node.dsize];
        self.src.read_exact_at(&mut v, base + PAGEHDRSZ as u64)?;
        Ok(v)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7531-7571
    /// Shape of the values of one leaf node of `db`.
    fn node_dups(&self, db: &Db, page: &Page, node: &Node) -> Result<Dups, ReadError> {
        if node.flags & F_DUPDATA == 0 {
            return Ok(Dups::One(self.node_value(page, node)?));
        }
        if db.flags & MDB_DUPSORT == 0 || node.flags & F_BIGDATA != 0 {
            return Err("LMDB: duplicate node in a database without MDB_DUPSORT, or with overflow data (MDB_CORRUPTED)".to_string().into());
        }
        if node.flags & F_SUBDATA != 0 {
            return Ok(Dups::SubDb(Db::parse(page.data(node))?));
        }
        let sub = Page::new(Rc::new(page.data(node).to_vec()))?;
        if sub.flags & P_SUBP == 0 || sub.flags & P_LEAF == 0 {
            return Err(format!("LMDB: duplicate sub-page has flags {:#x}, expected a leaf sub-page (MDB_CORRUPTED)", sub.flags).into());
        }
        Ok(Dups::SubPage(sub))
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6260-6269
    /// First value in the order of the duplicates (`mdb_cursor_first` of the sub-cursor).
    fn first_dup(&self, dups: &Dups) -> Result<Vec<u8>, ReadError> {
        match dups {
            Dups::One(v) => Ok(v.clone()),
            Dups::SubPage(p) => {
                if p.n == 0 {
                    return Err("LMDB: empty duplicate sub-page (MDB_CORRUPTED)"
                        .to_string()
                        .into());
                }
                Ok(p.dup_item(0)?.to_vec())
            }
            Dups::SubDb(sub) => {
                if sub.root == P_INVALID {
                    return Err("LMDB: empty duplicate sub-database (MDB_CORRUPTED)"
                        .to_string()
                        .into());
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
                                    .to_string()
                                    .into(),
                            );
                        }
                        return Ok(page.dup_item(0)?.to_vec());
                    } else {
                        return Err(format!(
                            "LMDB: sub-database page with flags {:#x} (MDB_CORRUPTED)",
                            page.flags
                        )
                        .into());
                    }
                }
                Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)"
                    .to_string()
                    .into())
            }
        }
    }

    /// Call `f` for every value of `dups`, in LMDB duplicate order (the sub-cursor walk of
    /// `values` and `for_each`, cited there).
    fn each_dup(
        &self,
        dups: &Dups,
        f: &mut dyn FnMut(&[u8]) -> Result<(), ReadError>,
    ) -> Result<(), ReadError> {
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

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5810-5857
    /// Leaf pages of the tree under `root` in key order (the page sequence of
    /// `mdb_cursor_first` + `mdb_cursor_sibling`). A tree that visits more pages than the
    /// database has is a cycle.
    fn walk_leaves(
        &self,
        root: u64,
        f: &mut dyn FnMut(&Page) -> Result<(), ReadError>,
    ) -> Result<(), ReadError> {
        let mut budget = self.npages;
        self.walk_rec(root, 0, &mut budget, f)
    }

    fn walk_rec(
        &self,
        pgno: u64,
        depth: usize,
        budget: &mut u64,
        f: &mut dyn FnMut(&Page) -> Result<(), ReadError>,
    ) -> Result<(), ReadError> {
        if depth > MAX_DEPTH {
            return Err("LMDB: tree is deeper than 32 levels (MDB_CORRUPTED)"
                .to_string()
                .into());
        }
        if *budget == 0 {
            return Err(
                "LMDB: tree visits more pages than the database has (MDB_CORRUPTED)"
                    .to_string()
                    .into(),
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
            )
            .into())
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6333-6359
    /// `MDB_SET` + `MDB_GET_CURRENT`: the first duplicate of `key` (the only value of a plain
    /// key). Works for sub-page and sub-database keys. `Ok(None)`: no such key.
    pub fn first_value(&self, db: &Db, key: &[u8]) -> Result<Option<Vec<u8>>, ReadError> {
        let Some((page, node)) = self.find(db, key)? else {
            return Ok(None);
        };
        let dups = self.node_dups(db, &page, &node)?;
        Ok(Some(self.first_dup(&dups)?))
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5860-5896
    /// `MDB_SET` + `MDB_GET_CURRENT` + `MDB_NEXT_DUP` until the end: every value of `key` in
    /// LMDB order (for `MDB_DUPFIXED` without `MDB_INTEGERDUP`: `memcmp` of the value bytes).
    /// Empty when the key is absent.
    pub fn values(&self, db: &Db, key: &[u8]) -> Result<Vec<Vec<u8>>, ReadError> {
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

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6232-6273
    /// `MDB_FIRST` then `MDB_NEXT` until the end: every `(key, value)` pair in tree order (a key
    /// with several duplicates appears once per duplicate). Used for `GetDBTaxIds`.
    pub fn for_each(
        &self,
        db: &Db,
        mut f: impl FnMut(&[u8], &[u8]) -> Result<(), ReadError>,
    ) -> Result<(), ReadError> {
        if db.root == P_INVALID {
            return Ok(());
        }
        self.walk_leaves(db.root, &mut |leaf: &Page| {
            if leaf.flags & P_LEAF2 != 0 {
                return Err("LMDB: LEAF2 page in a main tree (MDB_CORRUPTED)"
                    .to_string()
                    .into());
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
#[path = "lmdb_read_tests.rs"]
mod tests;
