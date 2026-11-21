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

//! Port of the write path of the LMDB copy in NCBI BLAST+ (`c++/src/util/lmdb/mdb.c`,
//! `midl.c`; "LMDB 0.9.24: (July 24, 2019)", `lmdb.h:199-213`) that makeblastdb 2.17.0 uses
//! for the `.pdb`/`.ndb` and `.ptf`/`.ntf` files of a version 5 database.
//!
//! The port is literal at the level of page buffers: LMDB never clears freed space, so the
//! gaps of copied pages, the slack of sub-pages and the reserved freelist records carry bytes
//! whose values depend on which malloc'd page buffer (`me_dpages`), scratch page (the split's
//! `copy`) or DUPSORT scratch area (`me_pbuf`) LMDB used before. Every function here keeps the
//! buffer history of `mdb.c` so that the files are byte-identical to makeblastdb's.
//!
//! What makeblastdb's operation sequence never reaches is not ported and fails with
//! [`LmdbError::Unreachable`] (survey A, DECISIONS 9 of the makedb-v5 line): read
//! transactions and the reader table, nested transactions, page spilling, loose pages,
//! `mdb_node_move`, `mdb_page_merge`, `mdb_update_key`, the rebalance of non-root pages,
//! `MDB_WRITEMAP`, `MDB_NOMEMINIT`, `MDB_MULTIPLE`, `MDB_REVERSEKEY`, `MDB_INTEGERDUP`,
//! sub-pages of non-`MDB_DUPFIXED` databases and the split of a LEAF2 page without
//! `MDB_APPEND`. The memory map and the syncs change no byte of the data file and are left out:
//! clean pages are read back from the file instead of the map. The lock file is created and
//! removed as LMDB and the writers do, with the bytes of a lock file that was never used
//! (`write_lock_file`).

mod cursor;
mod midl;
mod page;

#[cfg(test)]
mod tests;

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use super::platform::LmdbPlatform;
use midl::*;
use page::*;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/lmdb.h:312-349
pub(crate) const MDB_REVERSEKEY: u32 = 0x02;
pub(crate) const MDB_DUPSORT: u32 = 0x04;
pub(crate) const MDB_INTEGERKEY: u32 = 0x08;
pub(crate) const MDB_DUPFIXED: u32 = 0x10;
pub(crate) const MDB_INTEGERDUP: u32 = 0x20;
pub(crate) const MDB_REVERSEDUP: u32 = 0x40;
pub(crate) const MDB_CREATE: u32 = 0x40000;
pub(crate) const MDB_NOOVERWRITE: u32 = 0x10;
pub(crate) const MDB_NODUPDATA: u32 = 0x20;
pub(crate) const MDB_CURRENT: u32 = 0x40;
pub(crate) const MDB_RESERVE: u32 = 0x10000;
pub(crate) const MDB_APPEND: u32 = 0x20000;
pub(crate) const MDB_APPENDDUP: u32 = 0x40000;
pub(crate) const MDB_MULTIPLE: u32 = 0x80000;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/lmdb.h:287
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1266
const MDB_NOSUBDIR: u32 = 0x4000;
const MDB_ENV_ACTIVE: u32 = 0x2000_0000;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/lmdb.h:403-450
const MDB_SUCCESS: i32 = 0;
const MDB_KEYEXIST: i32 = -30799;
const MDB_NOTFOUND: i32 = -30798;
const MDB_PAGE_NOTFOUND: i32 = -30797;
const MDB_CORRUPTED: i32 = -30796;
const MDB_MAP_FULL: i32 = -30792;
const MDB_DBS_FULL: i32 = -30791;
const MDB_TXN_FULL: i32 = -30788;
const MDB_CURSOR_FULL: i32 = -30787;
const MDB_PAGE_FULL: i32 = -30786;
const MDB_MAP_RESIZED: i32 = -30785;
const MDB_INCOMPATIBLE: i32 = -30784;
const MDB_BAD_TXN: i32 = -30782;
const MDB_BAD_VALSIZE: i32 = -30781;
const MDB_LAST_ERRCODE: i32 = -30780;
const EINVAL: i32 = 22;
const EIO: i32 = 5;
const ENOMEM: i32 = 12;
/// Not an LMDB code: a path of mdb.c that makeblastdb does not reach and that is not
/// ported; the name of the path is in `Env::unreachable`.
const MDB_UNREACHABLE: i32 = -1_000_000;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1455-1477 (mdb_errstr)
const MDB_ERRSTR: [&str; 20] = [
    "MDB_KEYEXIST: Key/data pair already exists",
    "MDB_NOTFOUND: No matching key/data pair found",
    "MDB_PAGE_NOTFOUND: Requested page not found",
    "MDB_CORRUPTED: Located page was wrong type",
    "MDB_PANIC: Update of meta page failed or environment had fatal error",
    "MDB_VERSION_MISMATCH: Database environment version mismatch",
    "MDB_INVALID: File is not an LMDB file",
    "MDB_MAP_FULL: Environment mapsize limit reached",
    "MDB_DBS_FULL: Environment maxdbs limit reached",
    "MDB_READERS_FULL: Environment maxreaders limit reached",
    "MDB_TLS_FULL: Thread-local storage keys full - too many environments open",
    "MDB_TXN_FULL: Transaction has too many dirty pages - transaction too big",
    "MDB_CURSOR_FULL: Internal error - cursor stack limit reached",
    "MDB_PAGE_FULL: Internal error - page has no more space",
    "MDB_MAP_RESIZED: Database contents grew beyond environment mapsize",
    "MDB_INCOMPATIBLE: Operation and DB incompatible, or DB flags changed",
    "MDB_BAD_RSLOT: Invalid reuse of reader locktable slot",
    "MDB_BAD_TXN: Transaction must abort, has a child, or is invalid",
    "MDB_BAD_VALSIZE: Unsupported size of key/DB name/data, or wrong DUPFIXED size",
    "MDB_BAD_DBI: The specified DBI handle was closed/changed unexpectedly",
];

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:520-564, 1015-1030
const MAX_PAGESIZE: usize = 0x8000;
const MDB_MINKEYS: usize = 2;
const MDB_MAGIC: u32 = 0xBEEF_C0DE;
const MDB_DATA_VERSION: u32 = 1;
const MDB_MAXKEYSIZE: usize = 511;
const MAXDATASIZE: usize = 0xffff_ffff;
const P_INVALID: u64 = !0;
const DEFAULT_MAPSIZE: u64 = 1_048_576;
const MDB_VALID: u16 = 0x8000;
/// `(0xffff & ~(MDB_VALID))` in 16 bits.
const PERSISTENT_FLAGS: u16 = !MDB_VALID;
const FREE_DBI: usize = 0;
const MAIN_DBI: usize = 1;
const CORE_DBS: usize = 2;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1121-1170
const DB_DIRTY: u8 = 0x01;
const DB_STALE: u8 = 0x02;
const DB_NEW: u8 = 0x04;
const DB_VALID: u8 = 0x08;
const DB_USRVALID: u8 = 0x10;
const DB_DUPDATA: u8 = 0x20;
const MDB_TXN_FINISHED: u32 = 0x01;
const MDB_TXN_ERROR: u32 = 0x02;
const MDB_TXN_DIRTY: u32 = 0x04;
const MDB_TXN_SPILLS: u32 = 0x08;
const MDB_TXN_HAS_CHILD: u32 = 0x10;
const MDB_TXN_BLOCKED: u32 = MDB_TXN_FINISHED | MDB_TXN_ERROR | MDB_TXN_HAS_CHILD;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1353-1362
const MDB_END_UPDATE: u32 = 0x10;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1173-1195
const CURSOR_STACK: usize = 32;
const C_INITIALIZED: u32 = 0x01;
const C_EOF: u32 = 0x02;
const C_SUB: u32 = 0x04;
const C_DEL: u32 = 0x08;

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1369-1377, 6537
const MDB_PS_MODIFY: u32 = 1;
const MDB_PS_ROOTONLY: u32 = 2;
const MDB_PS_FIRST: u32 = 4;
const MDB_PS_LAST: u32 = 8;
const MDB_SPLIT_REPLACE: u32 = MDB_APPENDDUP;
const MDB_NOSPILL: u32 = 0x8000;

/// An error of the LMDB writer.
#[derive(Debug)]
pub(crate) enum LmdbError {
    /// An LMDB return code, with the name of the LMDB function the NCBI layer called
    /// (lmdb++ raises `lmdb::error` with `"<function>: <mdb_strerror>"`).
    Mdb { func: &'static str, code: i32 },
    /// The operating system failed to create, read or write the file.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A path of mdb.c that makeblastdb's operation sequence does not reach (not ported).
    Unreachable(&'static str),
}

impl LmdbError {
    /// `mdb_strerror` of an LMDB code.
    pub(crate) fn strerror(code: i32) -> String {
        if code == MDB_SUCCESS {
            return "Successful return: 0".to_string();
        }
        if (MDB_KEYEXIST..=MDB_LAST_ERRCODE).contains(&code) {
            return MDB_ERRSTR[(code - MDB_KEYEXIST) as usize].to_string();
        }
        // glibc's `strerror(ENOMEM)`, the same on every platform.
        if code == ENOMEM {
            return "Cannot allocate memory".to_string();
        }
        // `std::io::Error` prints a system error as the system's text (`strerror` on Unix) and
        // ` (os error N)`; `mdb_strerror` gives the text alone.
        let text = std::io::Error::from_raw_os_error(code).to_string();
        match text.strip_suffix(&format!(" (os error {code})")) {
            Some(strerror) => strerror.to_string(),
            None => text,
        }
    }

    /// The LMDB code, when this is an LMDB error.
    pub(crate) fn code(&self) -> Option<i32> {
        match self {
            LmdbError::Mdb { code, .. } => Some(*code),
            _ => None,
        }
    }
}

impl std::fmt::Display for LmdbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LmdbError::Mdb { func, code } => write!(f, "{func}: {}", Self::strerror(*code)),
            LmdbError::Io { path, source } => write!(f, "{}: {source}", path.display()),
            LmdbError::Unreachable(what) => {
                write!(f, "LMDB writer: unported LMDB path reached: {what}")
            }
        }
    }
}

impl std::error::Error for LmdbError {}

/// A comparison function of a database (`MDB_dbx.md_cmp`, `md_dcmp`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cmp {
    None,
    /// `mdb_cmp_memn`.
    Memn,
    /// `mdb_cmp_long`, `mdb_cmp_int`, `mdb_cmp_cint`: unsigned integers, little-endian.
    Int,
}

impl Cmp {
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5174-5241
    // Only the sign of a comparison is used by LMDB.
    fn cmp(self, a: &[u8], b: &[u8]) -> i32 {
        match self {
            Cmp::Memn => {
                let len = a.len().min(b.len());
                match a[..len].cmp(&b[..len]) {
                    std::cmp::Ordering::Less => -1,
                    std::cmp::Ordering::Greater => 1,
                    std::cmp::Ordering::Equal => match a.len().cmp(&b.len()) {
                        std::cmp::Ordering::Less => -1,
                        std::cmp::Ordering::Greater => 1,
                        std::cmp::Ordering::Equal => 0,
                    },
                }
            }
            Cmp::Int => {
                // mdb_cmp_cint compares 16-bit words from the most significant end, over the
                // size of `a`; for keys of one size it orders as unsigned integers.
                let n = a.len();
                let mut k = n;
                while k >= 2 {
                    let u = get_u16(a, k - 2) as i32;
                    let c = get_u16(b, k - 2) as i32;
                    let x = u - c;
                    if x != 0 {
                        return x.signum();
                    }
                    k -= 2;
                }
                0
            }
            Cmp::None => 0,
        }
    }
}

/// A page a cursor points at (`MDB_page *`): a page by number (the dirty page of that number
/// in this transaction, or the page in the file), a sub-page inside the data of a node of a
/// page, or a scratch page of `mdb_page_split`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PageRef {
    Null,
    Pg(u64),
    Sub(u64, u16),
    Scratch(u16),
}

/// Where the data of a node is (`MDB_val.mv_data` handed back for `MDB_RESERVE`).
#[derive(Clone, Copy, Debug)]
struct Loc {
    page: PageRef,
    off: usize,
}

/// A key or data item found by a cursor operation (`MDB_val` pointing into a page).
#[derive(Clone, Copy, Debug)]
struct Slot {
    page: PageRef,
    off: usize,
    len: usize,
}

/// Where `MDB_cursor.mc_db` (and `mc_dbflag`) point: the transaction's record of a
/// database, or the `mx_db` of the sub-cursor at an index of the cursor arena.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DbLoc {
    Txn(usize),
    X(usize),
}

/// `MDB_cursor` (and, for a sub-cursor, the `mx_db`/`mx_dbflag` of its `MDB_xcursor`).
#[derive(Clone, Debug)]
struct Cursor {
    dbi: usize,
    db: DbLoc,
    dbflag: DbLoc,
    cmp: Cmp,
    dcmp: Cmp,
    flags: u32,
    snum: usize,
    top: usize,
    pg: [PageRef; CURSOR_STACK],
    ki: [u16; CURSOR_STACK],
    /// `mc_xcursor`: the arena index of the sub-cursor.
    xc: Option<usize>,
    /// `mx_db` when this cursor is the `mx_cursor` of an xcursor.
    xdb: MdbDb,
    /// `mx_dbflag` when this cursor is the `mx_cursor` of an xcursor.
    xdbflag: u8,
}

/// An entry of `txn->mt_cursors[dbi]`: a cursor, or the dummy cursor that
/// `WITH_CURSOR_TRACKING` links for a sub-cursor copy (`mc_xcursor` = that copy).
#[derive(Clone, Copy, Debug)]
enum Track {
    Cursor(usize),
    Dummy(usize),
}

/// An LMDB environment open for writing with its one write transaction (`MDB_env`,
/// `MDB_txn` `me_txn0`).
pub(crate) struct Env {
    file: File,
    path: PathBuf,
    /// `<path>-lock`, removed when the environment is dropped. `None` until `create` has
    /// finished: an environment whose open failed has no destructor to remove the file.
    lock: Option<PathBuf>,
    psize: usize,
    me_flags: u32,
    mapsize: u64,
    maxpg: u64,
    nodemax: usize,
    maxfree_1pg: i64,
    metas: [Meta; 2],
    /// `me_txns->mti_txnid`.
    ti_txnid: u64,
    pgoldest: u64,
    pghead: Option<Vec<u64>>,
    pglast: u64,
    /// `me_dpages`: the LIFO of malloc'd page buffers (top = last element).
    dpages: Vec<Box<[u8]>>,
    /// `me_pbuf`.
    pbuf: Box<[u8]>,
    maxdbs: usize,
    numdbs: usize,
    dbx_name: Vec<Vec<u8>>,
    dbx_cmp: Vec<Cmp>,
    dbx_dcmp: Vec<Cmp>,
    env_dbflags: Vec<u16>,
    t_txnid: u64,
    t_next_pgno: u64,
    t_flags: u32,
    t_free_pgs: Vec<u64>,
    /// `mt_u.dirty_list`, sorted by page number.
    t_dirty: Vec<(u64, Box<[u8]>)>,
    t_dirty_room: usize,
    t_dbs: Vec<MdbDb>,
    t_dbflags: Vec<u8>,
    t_numdbs: usize,
    /// Pages of the file read in this transaction (the map in LMDB).
    clean: HashMap<u64, Box<[u8]>>,
    cursors: Vec<Cursor>,
    tracked: Vec<Vec<Track>>,
    /// The `copy` pages of the `mdb_page_split` calls in progress.
    scratch: Vec<Box<[u8]>>,
    wbuf: Vec<u8>,
    io_error: Option<std::io::Error>,
    unreachable: &'static str,
    /// LMDB sets the length of the data file to the map size (Windows, `mdb_env_map`).
    file_length_is_map_size: bool,
    /// The largest map size the platform's `size_t` holds.
    max_map_size: u64,
}

/// A database handle (`MDB_dbi`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Dbi(usize);

impl Env {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:78-112
    /// `CBlastLMDBManager::CBlastEnv` for writing: `mdb_env_create`,
    /// `mdb_env_set_maxdbs(3)`, `mdb_env_set_mapsize(map_size)` unless `map_size` is 0, and
    /// `mdb_env_open(path, MDB_NOSUBDIR, 0664)` on a new file, on the platform `platform`
    /// (its page size, the file's length on Windows, its `size_t`; [`super::platform`]).
    pub(crate) fn create(
        path: &Path,
        map_size: u64,
        platform: &LmdbPlatform,
    ) -> Result<Env, LmdbError> {
        let io = |source| LmdbError::Io {
            path: path.to_path_buf(),
            source,
        };
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3969,4400-4402
        let psize = platform.os_page_size.min(MAX_PAGESIZE);
        // lmdb++ (NCBI 598d8ae6 copy): c++/include/util/lmdbxx/lmdb++.h:1355-1358
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4025-4029 (mdb_env_map)
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4995-5004 (mdb_env_open)
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4160-4164, 4749
        // LMDB opens (creates) its lock file `<path>-lock` before the data file; one that cannot
        // be created stops `mdb_env_open` there (a name the file system cannot hold,
        // `ENAMETOOLONG`: with makeblastdb's `NAME.pdb`, an output name of 247 bytes and more
        // where a file name has at most 255; a directory of that name, which `open` refuses
        // with `EISDIR`; makeblastdb removes a directory of that name with its contents before
        // (`DeleteBlastDb`, `crate::writedb::dir_entry`), so the directory is one that could
        // not be removed). The file stays when a later step of `mdb_env_open` fails (the data
        // file that cannot be created): only the destructor of an environment that was opened
        // removes it (`Env::drop`).
        let lock = lock_path(path);
        write_lock_file(&lock).map_err(|source| LmdbError::Io {
            path: lock.clone(),
            source,
        })?;
        // The map size (`Uint8` in NCBI's layer) reaches LMDB as a `size_t`. A size the
        // platform's `size_t` cannot hold (a 32-bit platform) would be cut by that conversion;
        // LOSAT refuses it with the error the map of such a size gets (`ENOMEM`) instead, before
        // it makes the data file. LMDB would fail in `mdb_env_map`, after it has created the
        // data file: the message is the same, and only which error is reported when the data
        // file cannot be made either differs (PD-LOSAT-DATABASE-SEARCH lists the cut map size
        // as NCBI undefined behaviour that is not reproduced).
        if map_size > platform.max_map_size {
            return Err(LmdbError::Mdb {
                func: "mdb_env_open",
                code: ENOMEM,
            });
        }
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3951-3972, 4094-4101
        let maxdbs = 3 + CORE_DBS;
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4235-4325 (mdb_fopen)
        // makeblastdb deletes the files of an older database before it writes
        // (`build_db.cpp:1111`, `DeleteLMDBFiles`), so LMDB always creates a new file; a file
        // left at the path is replaced here.
        let mut options = std::fs::OpenOptions::new();
        options.read(true).write(true).create(true).truncate(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o664);
        }
        let file = options.open(path).map_err(io)?;
        let mut env = Env {
            file,
            path: path.to_path_buf(),
            lock: None,
            psize,
            me_flags: MDB_NOSUBDIR | MDB_ENV_ACTIVE,
            mapsize: 0,
            maxpg: 0,
            nodemax: 0,
            maxfree_1pg: 0,
            metas: [Meta::default(); 2],
            ti_txnid: 0,
            pgoldest: 0,
            pghead: None,
            pglast: 0,
            dpages: Vec::new(),
            pbuf: vec![0u8; psize].into_boxed_slice(),
            maxdbs,
            numdbs: CORE_DBS,
            dbx_name: vec![Vec::new(); maxdbs],
            dbx_cmp: vec![Cmp::None; maxdbs],
            dbx_dcmp: vec![Cmp::None; maxdbs],
            env_dbflags: vec![0; maxdbs],
            t_txnid: 0,
            t_next_pgno: 0,
            t_flags: MDB_TXN_FINISHED,
            t_free_pgs: midl_new(),
            t_dirty: Vec::new(),
            t_dirty_room: 0,
            t_dbs: vec![MdbDb::default(); maxdbs],
            t_dbflags: vec![0; maxdbs],
            t_numdbs: 0,
            clean: HashMap::new(),
            cursors: Vec::new(),
            tracked: vec![Vec::new(); maxdbs],
            scratch: Vec::new(),
            wbuf: Vec::new(),
            io_error: None,
            unreachable: "",
            file_length_is_map_size: platform.file_length_is_map_size,
            max_map_size: platform.max_map_size,
        };
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4060-4091 (mdb_env_set_mapsize)
        if map_size != 0 {
            env.mapsize = map_size;
        }
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4993
        env.dbx_cmp[FREE_DBI] = Cmp::Int;
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4400-4456 (mdb_env_open2)
        let mut meta = env.init_meta0();
        meta.mapsize = DEFAULT_MAPSIZE;
        if env.mapsize == 0 {
            env.mapsize = meta.mapsize;
        }
        let minsize = (meta.last_pg + 1) * meta.dbs[FREE_DBI].pad as u64;
        if env.mapsize < minsize {
            env.mapsize = minsize;
        }
        meta.mapsize = env.mapsize;
        env.init_meta(&meta).map_err(io)?;
        env.map_file_size().map_err(io)?;
        env.maxfree_1pg = ((psize - PAGEHDRSZ) / 8) as i64 - 1;
        env.nodemax = (((psize - PAGEHDRSZ) / MDB_MINKEYS) & !1) - 2;
        env.maxpg = env.mapsize / psize as u64;
        // `mdb_env_open` has succeeded: from here on the writer's destructor removes the file.
        env.lock = Some(lock);
        Ok(env)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3740-3753 (mdb_env_init_meta0)
    fn init_meta0(&self) -> Meta {
        let mut meta = Meta {
            magic: MDB_MAGIC,
            version: MDB_DATA_VERSION,
            mapsize: self.mapsize,
            last_pg: 1,
            ..Meta::default()
        };
        meta.dbs[FREE_DBI].pad = self.psize as u32;
        meta.dbs[FREE_DBI].flags = (self.me_flags & 0xffff) as u16 | MDB_INTEGERKEY as u16;
        meta.dbs[FREE_DBI].root = P_INVALID;
        meta.dbs[MAIN_DBI].root = P_INVALID;
        meta
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3759-3807 (mdb_env_init_meta)
    fn init_meta(&mut self, meta: &Meta) -> std::io::Result<()> {
        let psize = self.psize;
        let mut p = vec![0u8; 2 * psize];
        let mb = meta.to_bytes();
        for i in 0..2 {
            let q = &mut p[i * psize..(i + 1) * psize];
            set_pgno(q, i as u64);
            set_flags(q, P_META);
            q[PAGEHDRSZ..PAGEHDRSZ + META_SIZE].copy_from_slice(&mb);
        }
        self.file.seek(SeekFrom::Start(0))?;
        self.file.write_all(&p)?;
        self.metas = [*meta; 2];
        Ok(())
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3944-3948 (mdb_env_pick_meta)
    fn pick_meta(&self) -> &Meta {
        &self.metas[(self.metas[0].txnid < self.metas[1].txnid) as usize]
    }

    /// `mdb_env_stat(env).ms_psize`.
    pub(crate) fn psize(&self) -> usize {
        self.psize
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:9688-9705 (mdb_env_info)
    /// `mdb_env_info(env).me_last_pgno`.
    pub(crate) fn last_pgno(&self) -> u64 {
        self.pick_meta().last_pg
    }

    /// `mdb_env_info(env).me_mapsize`.
    pub(crate) fn mapsize(&self) -> u64 {
        self.mapsize
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4060-4091 (mdb_env_set_mapsize)
    /// `mdb_env_set_mapsize` on the open environment (the map itself is not kept).
    pub(crate) fn set_mapsize(&mut self, size: u64) -> Result<(), LmdbError> {
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            return Err(LmdbError::Mdb {
                func: "mdb_env_set_mapsize",
                code: EINVAL,
            });
        }
        // A size the platform's `size_t` cannot hold, as in `create`.
        if size > self.max_map_size {
            return Err(LmdbError::Mdb {
                func: "mdb_env_set_mapsize",
                code: ENOMEM,
            });
        }
        let meta = *self.pick_meta();
        let mut size = if size == 0 { meta.mapsize } else { size };
        let minsize = (meta.last_pg + 1) * self.psize as u64;
        if size < minsize {
            size = minsize;
        }
        self.mapsize = size;
        self.map_file_size().map_err(|source| LmdbError::Io {
            path: self.path.clone(),
            source,
        })?;
        self.maxpg = self.mapsize / self.psize as u64;
        Ok(())
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3979-4003 (mdb_env_map)
    // On Windows LMDB sets the length of the data file to the map size whenever it maps the
    // file: when the environment opens (after the meta pages are written, mdb.c:4421-4436)
    // and when `mdb_env_set_mapsize` maps it again with a new size (mdb.c:4083-4087). The
    // pages are written inside that length, so a file of makeblastdb on Windows is as long
    // as the last map size (500000 bytes by default, `DEFAULT_LMDB_MAP_SIZE`). Elsewhere
    // the file is not extended (`MDB_WRITEMAP` is not set). Which applies is the platform's
    // (`LmdbPlatform::file_length_is_map_size`). Not checked against a makeblastdb.exe.
    /// The length `mdb_env_map` gives the data file: the map size on Windows.
    fn map_file_size(&mut self) -> std::io::Result<()> {
        if self.file_length_is_map_size {
            self.file.set_len(self.mapsize)?;
        }
        Ok(())
    }

    /// The error of an internal return code, for the NCBI layer's function `func`.
    fn error(&mut self, func: &'static str, rc: i32) -> LmdbError {
        if rc == MDB_UNREACHABLE {
            return LmdbError::Unreachable(self.unreachable);
        }
        if let Some(source) = self.io_error.take() {
            return LmdbError::Io {
                path: self.path.clone(),
                source,
            };
        }
        LmdbError::Mdb { func, code: rc }
    }

    /// Fail on a path that is not ported.
    fn unreachable(&mut self, what: &'static str) -> i32 {
        self.unreachable = what;
        self.t_flags |= MDB_TXN_ERROR;
        MDB_UNREACHABLE
    }

    // ---- transactions ----------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2747-2809
    // (mdb_txn_renew0, write txn)
    /// `mdb_txn_begin(env, NULL, 0)` (lmdb++ `txn::begin`).
    pub(crate) fn txn_begin(&mut self) -> Result<(), LmdbError> {
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            return Err(self.error("mdb_txn_begin", EINVAL));
        }
        self.t_txnid = self.ti_txnid;
        let meta = self.metas[(self.t_txnid & 1) as usize];
        self.t_txnid += 1;
        self.t_dirty_room = MDB_IDL_UM_MAX;
        self.t_dirty.clear();
        self.t_free_pgs.clear();
        self.t_free_pgs.push(0);
        self.t_dbs[FREE_DBI] = meta.dbs[FREE_DBI];
        self.t_dbs[MAIN_DBI] = meta.dbs[MAIN_DBI];
        self.t_next_pgno = meta.last_pg + 1;
        self.t_flags = 0;
        self.t_numdbs = self.numdbs;
        for i in CORE_DBS..self.t_numdbs {
            let x = self.env_dbflags[i];
            self.t_dbs[i].flags = x & PERSISTENT_FLAGS;
            self.t_dbflags[i] = if x & MDB_VALID != 0 {
                DB_VALID | DB_USRVALID | DB_STALE
            } else {
                0
            };
        }
        self.t_dbflags[MAIN_DBI] = DB_VALID | DB_USRVALID;
        self.t_dbflags[FREE_DBI] = DB_VALID;
        if self.maxpg < self.t_next_pgno {
            self.txn_end(0);
            return Err(self.error("mdb_txn_begin", MDB_MAP_RESIZED));
        }
        Ok(())
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2943-2970 (mdb_dbis_update)
    fn dbis_update(&mut self, keep: bool) {
        let n = self.t_numdbs;
        for i in (CORE_DBS..n).rev() {
            if self.t_dbflags[i] & DB_NEW != 0 {
                if keep {
                    self.env_dbflags[i] = self.t_dbs[i].flags | MDB_VALID;
                } else if !self.dbx_name[i].is_empty() {
                    self.dbx_name[i].clear();
                    self.env_dbflags[i] = 0;
                }
            }
        }
        if keep && self.numdbs < n {
            self.numdbs = n;
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2976-3045 (mdb_txn_end)
    fn txn_end(&mut self, mode: u32) {
        self.dbis_update(mode & MDB_END_UPDATE != 0);
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1831-1843 (mdb_dlist_free)
            let dl = std::mem::take(&mut self.t_dirty);
            for (_, dp) in dl {
                self.dpage_free(dp);
            }
            self.t_numdbs = 0;
            self.t_flags = MDB_TXN_FINISHED;
            self.pghead = None;
            self.pglast = 0;
        }
        self.clean.clear();
        self.cursors.clear();
        for t in &mut self.tracked {
            t.clear();
        }
        self.scratch.clear();
    }

    /// Abort the transaction after an error (`mdb_txn_abort`; lmdb++'s `txn` destructor).
    pub(crate) fn txn_abort(&mut self) {
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            self.txn_end(0);
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3448-3674 (mdb_txn_commit)
    // `mdb_env_sync` (fdatasync) changes no byte of the file and is left out.
    /// `mdb_txn_commit` (lmdb++ `txn::commit`).
    pub(crate) fn txn_commit(&mut self) -> Result<(), LmdbError> {
        let rc = self.txn_commit0();
        if rc != MDB_SUCCESS {
            self.txn_abort();
            return Err(self.error("mdb_txn_commit", rc));
        }
        Ok(())
    }

    fn txn_commit0(&mut self) -> i32 {
        if self.t_flags & (MDB_TXN_FINISHED | MDB_TXN_ERROR) != 0 {
            return MDB_BAD_TXN;
        }
        if self.t_dirty.is_empty() && self.t_flags & (MDB_TXN_DIRTY | MDB_TXN_SPILLS) == 0 {
            self.txn_end(MDB_END_UPDATE);
            return MDB_SUCCESS;
        }
        if self.t_numdbs > CORE_DBS {
            let mc = self.cursor_init(MAIN_DBI);
            for i in CORE_DBS..self.t_numdbs {
                if self.t_dbflags[i] & DB_DIRTY != 0 {
                    let data = self.t_dbs[i].to_bytes();
                    let name = self.dbx_name[i].clone();
                    let (rc, _) = self.cursor_put(mc, &name, Val::Bytes(&data), F_SUBDATA as u32);
                    if rc != MDB_SUCCESS {
                        return rc;
                    }
                }
            }
            self.cursors.truncate(mc);
        }
        let rc = self.freelist_save();
        if rc != MDB_SUCCESS {
            return rc;
        }
        self.pghead = None;
        let rc = self.page_flush();
        if rc != MDB_SUCCESS {
            return rc;
        }
        let rc = self.write_meta();
        if rc != MDB_SUCCESS {
            return rc;
        }
        self.txn_end(MDB_END_UPDATE);
        MDB_SUCCESS
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3812-3937 (mdb_env_write_meta, no
    // MDB_WRITEMAP)
    fn write_meta(&mut self) -> i32 {
        let toggle = (self.t_txnid & 1) as usize;
        let mut mapsize = self.metas[toggle ^ 1].mapsize;
        if mapsize < self.mapsize {
            mapsize = self.mapsize;
        }
        let mut meta = self.metas[toggle];
        meta.mapsize = mapsize;
        meta.dbs[FREE_DBI] = self.t_dbs[FREE_DBI];
        meta.dbs[MAIN_DBI] = self.t_dbs[MAIN_DBI];
        meta.last_pg = self.t_next_pgno - 1;
        meta.txnid = self.t_txnid;
        let bytes = meta.to_bytes();
        let off = (toggle * self.psize + PAGEHDRSZ + META_MAPSIZE_OFF) as u64;
        let r = self
            .file
            .seek(SeekFrom::Start(off))
            .and_then(|_| self.file.write_all(&bytes[META_MAPSIZE_OFF..]));
        if let Err(e) = r {
            self.io_error = Some(e);
            return EIO;
        }
        self.metas[toggle] = meta;
        self.ti_txnid = self.t_txnid;
        MDB_SUCCESS
    }

    // ---- page buffers ----------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1768-1805 (mdb_page_malloc)
    // A recycled buffer keeps all its bytes. A new buffer is zeroed here in full: LMDB zeroes
    // all of a single page but the header (whose fields its callers then set) and the last
    // page of a multi-page buffer; the earlier pages of an overflow buffer are covered by the
    // overflow data, which always reaches into the last page (`OVPAGES`), before they are
    // written.
    fn page_malloc(&mut self, num: usize) -> Box<[u8]> {
        if num == 1 {
            if let Some(b) = self.dpages.pop() {
                return b;
            }
        }
        vec![0u8; self.psize * num].into_boxed_slice()
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1807-1829 (mdb_page_free,
    // mdb_dpage_free)
    // `mp_next` (bytes 0-7) is always overwritten by `mp_pgno` before a buffer is written, so
    // the pointer is not stored.
    fn page_free(&mut self, mp: Box<[u8]>) {
        self.dpages.push(mp);
    }

    fn dpage_free(&mut self, dp: Box<[u8]>) {
        if !is_overflow(&dp) || pages(&dp) == 1 {
            self.page_free(dp);
        }
    }

    /// The index of page `pgno` in the dirty list (`mdb_mid2l_search`).
    fn dirty_index(&self, pgno: u64) -> Option<usize> {
        self.t_dirty.binary_search_by(|e| e.0.cmp(&pgno)).ok()
    }

    /// The bytes of a page the cursors point at. A page by number is the dirty page of that
    /// number, or the page read from the file by `page_get`.
    fn page(&self, r: PageRef) -> &[u8] {
        match r {
            PageRef::Pg(p) => self.page_by_pgno(p),
            PageRef::Sub(h, off) => &self.page_by_pgno(h)[off as usize..],
            PageRef::Scratch(k) => &self.scratch[k as usize],
            PageRef::Null => panic!("LMDB writer: NULL page"),
        }
    }

    fn page_by_pgno(&self, p: u64) -> &[u8] {
        match self.dirty_index(p) {
            Some(i) => &self.t_dirty[i].1,
            None => self
                .clean
                .get(&p)
                .map(|b| &b[..])
                .expect("LMDB writer: page not read"),
        }
    }

    /// The bytes of a writable page (a dirty page, a sub-page of one, or a scratch page).
    fn page_mut(&mut self, r: PageRef) -> &mut [u8] {
        match r {
            PageRef::Pg(p) => {
                let i = self.dirty_index(p).expect("LMDB writer: page not dirty");
                &mut self.t_dirty[i].1
            }
            PageRef::Sub(h, off) => {
                let i = self.dirty_index(h).expect("LMDB writer: page not dirty");
                &mut self.t_dirty[i].1[off as usize..]
            }
            PageRef::Scratch(k) => &mut self.scratch[k as usize],
            PageRef::Null => panic!("LMDB writer: NULL page"),
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5428-5479 (mdb_page_get)
    /// `mdb_page_get`; also says whether the page is dirty in this transaction (`level` 1)
    /// or in the file (`level` 0). A page of the file is read once per transaction.
    fn page_get(&mut self, pgno: u64) -> Result<(PageRef, bool), i32> {
        if self.dirty_index(pgno).is_some() {
            return Ok((PageRef::Pg(pgno), true));
        }
        if pgno >= self.t_next_pgno {
            self.t_flags |= MDB_TXN_ERROR;
            return Err(MDB_PAGE_NOTFOUND);
        }
        if !self.clean.contains_key(&pgno) {
            match self.read_page(pgno) {
                Ok(b) => {
                    self.clean.insert(pgno, b);
                }
                Err(e) => {
                    self.io_error = Some(e);
                    self.t_flags |= MDB_TXN_ERROR;
                    return Err(EIO);
                }
            }
        }
        Ok((PageRef::Pg(pgno), false))
    }

    /// Read page `pgno` from the file, with all its pages when it is an overflow page.
    fn read_page(&mut self, pgno: u64) -> std::io::Result<Box<[u8]>> {
        let psize = self.psize;
        let mut b = vec![0u8; psize];
        self.file.seek(SeekFrom::Start(pgno * psize as u64))?;
        self.file.read_exact(&mut b)?;
        if is_overflow(&b) && pages(&b) > 1 {
            let n = pages(&b) as usize;
            b.resize(n * psize, 0);
            self.file.read_exact(&mut b[psize..])?;
        }
        Ok(b.into_boxed_slice())
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2099-2116 (mdb_find_oldest)
    // makeblastdb opens no read transaction on the environment it writes: the reader table
    // is empty.
    fn find_oldest(&self) -> u64 {
        self.t_txnid - 1
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2118-2135 (mdb_page_dirty)
    fn page_dirty(&mut self, pgno: u64, np: Box<[u8]>) {
        let x = match self.t_dirty.binary_search_by(|e| e.0.cmp(&pgno)) {
            Ok(_) => panic!("LMDB writer: page {pgno} dirty twice"),
            Err(x) => x,
        };
        self.t_dirty.insert(x, (pgno, np));
        self.clean.remove(&pgno);
        self.t_dirty_room -= 1;
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2137-2322 (mdb_page_alloc)
    // No loose pages exist (they are made only by `mdb_page_merge`); `MDB_PARANOID` is off.
    /// `mdb_page_alloc`: the number of the new dirty page.
    fn page_alloc(&mut self, num: usize) -> Result<u64, i32> {
        let mut retry: i64 = num as i64 * 60;
        let n2 = num - 1;
        let mut mop_len = self.pghead.as_ref().map_or(0, |m| m[0] as usize);
        let mut oldest: u64 = 0;
        let mut last: u64 = 0;
        let mut found_old = false;
        let mut first = true;
        let mut set_range = false;
        let mut m2: Option<usize> = None;
        if self.t_dirty_room == 0 {
            self.t_flags |= MDB_TXN_ERROR;
            return Err(MDB_TXN_FULL);
        }
        let mut found: Option<(usize, u64)> = None;
        'fetch: loop {
            if mop_len > n2 {
                let mop = self.pghead.as_ref().unwrap();
                let mut i = mop_len;
                loop {
                    let pgno = mop[i];
                    if mop[i - n2] == pgno + n2 as u64 {
                        found = Some((i, pgno));
                        break 'fetch;
                    }
                    i -= 1;
                    if i <= n2 {
                        break;
                    }
                }
                retry -= 1;
                if retry < 0 {
                    break;
                }
            }
            if first {
                last = self.pglast;
                oldest = self.pgoldest;
                m2 = Some(self.cursor_init(FREE_DBI));
                if last != 0 {
                    set_range = true;
                }
            }
            last += 1;
            if oldest <= last {
                if !found_old {
                    oldest = self.find_oldest();
                    self.pgoldest = oldest;
                    found_old = true;
                }
                if oldest <= last {
                    break;
                }
            }
            let c = m2.unwrap();
            let (rc, key) = if first && set_range {
                // `key.mv_data = &last`: the search key is `last` after its increment.
                let k = last.to_le_bytes();
                let (rc, key, _) = self.cursor_set(c, &k, false, CursorOp::SetRange, None);
                (rc, key)
            } else if first {
                let (rc, key, _) = self.cursor_first(c, false);
                (rc, key)
            } else {
                let (rc, key, _) = self.cursor_next(c, false);
                (rc, key)
            };
            first = false;
            if rc != MDB_SUCCESS {
                if rc == MDB_NOTFOUND {
                    break;
                }
                self.t_flags |= MDB_TXN_ERROR;
                return Err(rc);
            }
            let key = key.unwrap();
            last = get_u64(self.page(key.page), key.off);
            if oldest <= last {
                if !found_old {
                    oldest = self.find_oldest();
                    self.pgoldest = oldest;
                    found_old = true;
                }
                if oldest <= last {
                    break;
                }
            }
            let cur = &self.cursors[c];
            let (np, ki) = (cur.pg[cur.top], cur.ki[cur.top] as usize);
            let leaf = nodeptr(self.page(np), ki);
            let data = match self.node_read(np, leaf) {
                Ok(d) => d,
                Err(rc) => {
                    self.t_flags |= MDB_TXN_ERROR;
                    return Err(rc);
                }
            };
            let idl: Vec<u64> = {
                let b = &self.page(data.page)[data.off..data.off + data.len];
                let n = get_u64(b, 0) as usize;
                (0..=n).map(|k| get_u64(b, 8 * k)).collect()
            };
            let mop = self.pghead.get_or_insert_with(midl_new);
            self.pglast = last;
            midl_xmerge(mop, &idl);
            mop_len = mop[0] as usize;
        }
        if let Some(c) = m2 {
            self.cursors.truncate(c);
        }
        let (i, pgno) = match found {
            Some(f) => f,
            None => {
                let pgno = self.t_next_pgno;
                if pgno + num as u64 >= self.maxpg {
                    self.t_flags |= MDB_TXN_ERROR;
                    return Err(MDB_MAP_FULL);
                }
                (0, pgno)
            }
        };
        let mut np = self.page_malloc(num);
        if i != 0 {
            let mop = self.pghead.as_mut().unwrap();
            mop_len -= num;
            mop[0] = mop_len as u64;
            let mut j = i - num;
            let mut i = i;
            while j < mop_len {
                j += 1;
                i += 1;
                mop[j] = mop[i];
            }
        } else {
            self.t_next_pgno = pgno + num as u64;
        }
        set_pgno(&mut np, pgno);
        self.page_dirty(pgno, np);
        Ok(pgno)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2324-2346 (mdb_page_copy)
    fn page_copy(dst: &mut [u8], src: &[u8], psize: usize) {
        let upper = upper(src);
        let lower = lower(src);
        let unused = upper.wrapping_sub(lower) & !7u16;
        if unused != 0 && !is_leaf2(src) {
            let upper = (upper & !7u16) as usize;
            let lo = (lower as usize + 7) & !7usize;
            dst[..lo].copy_from_slice(&src[..lo]);
            dst[upper..psize].copy_from_slice(&src[upper..psize]);
        } else {
            let n = psize - unused as usize;
            dst[..n].copy_from_slice(&src[..n]);
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3295-3445
    // (mdb_page_flush, keep 0, no MDB_WRITEMAP)
    // No page is P_LOOSE or P_KEEP (no loose pages, no spilling).
    fn page_flush(&mut self) -> i32 {
        let psize = self.psize;
        let dl = std::mem::take(&mut self.t_dirty);
        let mut wbuf = std::mem::take(&mut self.wbuf);
        let mut wpos: u64 = 0;
        let mut next_pos: u64 = u64::MAX;
        let mut n = 0usize;
        let mut result: std::io::Result<()> = Ok(());
        let mut dl = dl;
        for (pgno, dp) in dl.iter_mut() {
            let f = flags(dp) & !P_DIRTY;
            set_flags(dp, f);
            let pos = *pgno * psize as u64;
            let size = if is_overflow(dp) {
                psize * pages(dp) as usize
            } else {
                psize
            };
            if pos != next_pos || n == 64 || wbuf.len() + size > 0x4000_0000 {
                if n > 0 && result.is_ok() {
                    result = self.write_at(wpos, &wbuf);
                }
                wbuf.clear();
                n = 0;
                wpos = pos;
            }
            next_pos = pos + size as u64;
            wbuf.extend_from_slice(&dp[..size]);
            n += 1;
        }
        if n > 0 && result.is_ok() {
            result = self.write_at(wpos, &wbuf);
        }
        wbuf.clear();
        self.wbuf = wbuf;
        let count = dl.len();
        for (_, dp) in dl {
            self.dpage_free(dp);
        }
        self.t_dirty_room += count;
        if let Err(e) = result {
            self.io_error = Some(e);
            return EIO;
        }
        MDB_SUCCESS
    }

    fn write_at(&mut self, pos: u64, b: &[u8]) -> std::io::Result<()> {
        self.file.seek(SeekFrom::Start(pos))?;
        self.file.write_all(b)
    }

    /// Close the environment (`mdb_env_close`, which writes nothing); the lock file goes
    /// with it (`Env::drop`), as `CWriteDB_LMDB` and `CWriteDB_TaxID` remove it after it.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72
    pub(crate) fn close(mut self) -> Result<(), LmdbError> {
        self.txn_abort();
        self.file.flush().map_err(|source| LmdbError::Io {
            path: self.path.clone(),
            source,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72,437-444
// The destructors of `CWriteDB_LMDB` and `CWriteDB_TaxID` remove the lock file whenever the
// writer is destroyed, also when an exception of the build unwinds the stack; a writer whose
// constructor did not finish (the data file that `mdb_env_open` could not create) has no
// destructor, and its lock file stays (a data file that cannot be made, a header that cannot
// be written, a file that cannot be given its length: probes of 2026-11-11 and 2026-11-12). An
// `Env` that `create` has not finished (`lock` is `None`) is that writer: it removes nothing.
// The removal is `CFile(name).Remove()`: `CDirEntry::Remove(eRecursive)`; its result is not
// looked at. Where makeblastdb aborts instead of unwinding (a destructor that throws), it leaves
// the lock files; LOSAT, which ends with an error there, removes them (PD-LOSAT-DATABASE-SEARCH).
impl Drop for Env {
    fn drop(&mut self) {
        if let Some(lock) = &self.lock {
            let _ = crate::writedb::dir_entry::remove(lock);
        }
    }
}

/// `<path>-lock`: LMDB's name of the lock file of a `MDB_NOSUBDIR` environment.
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4160-4164
fn lock_path(path: &Path) -> PathBuf {
    let mut lock = path.as_os_str().to_owned();
    lock.push("-lock");
    PathBuf::from(lock)
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4739-4800,4880-4915
// (mdb_env_setup_locks)
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:540,780-784
// (MDB_MAGIC, MDB_LOCK_FORMAT)
/// The size of a new lock file: `(126 - 1) * sizeof(MDB_reader) + sizeof(MDB_txninfo)`,
/// 125 reader slots of 64 bytes and a header of three 64-byte cache lines.
const LOCK_FILE_SIZE: usize = 8192;

/// The bytes of a lock file that no transaction has used, as makeblastdb 2.17.0 of the Linux
/// x86-64 release (glibc, `MDB_USE_POSIX_MUTEX` with robust process-shared mutexes) leaves one
/// whose environment did not finish opening (the data file or its header could not be written)
/// and the one of an environment that no write transaction began (the same bytes in every run
/// of the probes of 2026-11-11 and 2026-11-12). The lock file of an environment that was opened
/// and then aborted holds more (`0x44` = 1 once a write transaction began, `0x30` = 1 after
/// a commit); LOSAT ends an aborting build with an error and removes the file instead, so it
/// never leaves those. The bytes are `mti_magic` and `mti_format` (`MDB_LOCK_VERSION` 1 in the
/// low half, `MDB_PIDLOCK` in bit 16), and the `__kind` of the read mutex (offset 0x18) and of the write mutex (0x50), which is
/// `PTHREAD_MUTEX_ROBUST_NORMAL_NP | PTHREAD_PROCESS_SHARED` (0x10 + 0x80). Every other byte is 0.
/// Other platforms keep other things in the file (a macOS semaphore name, a Windows mutex
/// name): their bytes are not compared (PD-LOSAT-DATABASE-SEARCH).
fn lock_file_bytes() -> [u8; LOCK_FILE_SIZE] {
    let mut bytes = [0u8; LOCK_FILE_SIZE];
    bytes[0..4].copy_from_slice(&0xBEEF_C0DEu32.to_le_bytes());
    bytes[4..8].copy_from_slice(&0x0001_0001u32.to_le_bytes());
    bytes[0x18] = 0x90;
    bytes[0x50] = 0x90;
    bytes
}

/// `mdb_fopen(env, fname, MDB_O_LOCKS, 0664, ...)` and the initialisation of the lock region:
/// the file is created (or an existing one opened) and holds the bytes of a new lock file.
fn write_lock_file(lock: &Path) -> std::io::Result<()> {
    let mut options = std::fs::OpenOptions::new();
    options.read(true).write(true).create(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o664);
    }
    let mut file = options.open(lock)?;
    file.set_len(LOCK_FILE_SIZE as u64)?;
    file.write_all(&lock_file_bytes())
}

/// The data of a put: bytes, or the size of an `MDB_RESERVE` item whose bytes the caller
/// writes afterwards.
#[derive(Clone, Copy, Debug)]
enum Val<'a> {
    Bytes(&'a [u8]),
    Reserve(usize),
}

impl Val<'_> {
    fn size(&self) -> usize {
        match self {
            Val::Bytes(b) => b.len(),
            Val::Reserve(n) => *n,
        }
    }
}

/// The cursor operations the writer uses (`MDB_cursor_op`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CursorOp {
    Set,
    SetRange,
}
