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
//! `MDB_APPEND`. The lock file, the memory map and the syncs change no byte of the data file
//! and are left out: clean pages are read back from the file instead of the map.

mod cursor;
mod midl;
mod page;

#[cfg(test)]
mod tests;

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use midl::*;
use page::*;

// NCBI reference (598d8ae6): c++/src/util/lmdb/lmdb.h:312-349
// ```c
// #define MDB_DUPSORT		0x04
// #define MDB_INTEGERKEY	0x08
// #define MDB_DUPFIXED	0x10
// #define MDB_INTEGERDUP	0x20
// #define MDB_CREATE		0x40000
// #define MDB_NOOVERWRITE	0x10
// #define MDB_NODUPDATA	0x20
// #define MDB_CURRENT	0x40
// #define MDB_RESERVE	0x10000
// #define MDB_APPEND	0x20000
// #define MDB_APPENDDUP	0x40000
// #define MDB_MULTIPLE	0x80000
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/lmdb.h:287 and mdb.c:1316-1319
// ```c
// #define MDB_NOSUBDIR	0x4000
// #define	MDB_ENV_ACTIVE	0x20000000U
// ```
const MDB_NOSUBDIR: u32 = 0x4000;
const MDB_ENV_ACTIVE: u32 = 0x2000_0000;

// NCBI reference (598d8ae6): c++/src/util/lmdb/lmdb.h:403-450
// ```c
// #define MDB_SUCCESS	 0
// #define MDB_KEYEXIST	(-30799)
// #define MDB_NOTFOUND	(-30798)
// ...
// #define MDB_BAD_DBI		(-30780)
// ```
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
/// Not an LMDB code: a path of mdb.c that makeblastdb does not reach and that is not
/// ported; the name of the path is in `Env::unreachable`.
const MDB_UNREACHABLE: i32 = -1_000_000;

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1455-1477 (mdb_errstr)
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:520-564, 1015-1030
// ```c
// #define MAX_PAGESIZE	 (PAGEBASE ? 0x10000 : 0x8000)
// #define MDB_MINKEYS	 2
// #define MDB_MAGIC	 0xBEEFC0DE
// #define MDB_DATA_VERSION	 ((MDB_DEVEL) ? 999 : 1)
// #define MDB_MAXKEYSIZE	 ((MDB_DEVEL) ? 0 : 511)
// #define MAXDATASIZE	0xffffffffUL
// #define P_INVALID	 (~(pgno_t)0)
// #define DEFAULT_MAPSIZE	1048576
// #define MDB_VALID	0x8000		/**< DB handle is valid, for me_dbflags */
// #define PERSISTENT_FLAGS	(0xffff & ~(MDB_VALID))
// #define	FREE_DBI	0
// #define	MAIN_DBI	1
// #define CORE_DBS	2
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1121-1170
// ```c
// #define DB_DIRTY	0x01		/**< DB was written in this txn */
// #define DB_STALE	0x02		/**< Named-DB record is older than txnID */
// #define DB_NEW		0x04		/**< Named-DB handle opened in this txn */
// #define DB_VALID	0x08		/**< DB handle is valid, see also #MDB_VALID */
// #define DB_USRVALID	0x10		/**< As #DB_VALID, but not set for #FREE_DBI */
// #define DB_DUPDATA	0x20		/**< DB is #MDB_DUPSORT data */
// #define MDB_TXN_FINISHED	0x01
// #define MDB_TXN_ERROR		0x02
// #define MDB_TXN_DIRTY		0x04
// #define MDB_TXN_SPILLS		0x08
// #define MDB_TXN_HAS_CHILD	0x10
// #define MDB_TXN_BLOCKED		(MDB_TXN_FINISHED|MDB_TXN_ERROR|MDB_TXN_HAS_CHILD)
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1353-1362
// ```c
// #define MDB_END_UPDATE	0x10	/**< update env state (DBIs) */
// ```
const MDB_END_UPDATE: u32 = 0x10;

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1173-1195
// ```c
// #define CURSOR_STACK		 32
// #define C_INITIALIZED	0x01	/**< cursor has been initialized and is valid */
// #define C_EOF	0x02			/**< No more data */
// #define C_SUB	0x04			/**< Cursor is a sub-cursor */
// #define C_DEL	0x08			/**< last op was a cursor_del */
// ```
const CURSOR_STACK: usize = 32;
const C_INITIALIZED: u32 = 0x01;
const C_EOF: u32 = 0x02;
const C_SUB: u32 = 0x04;
const C_DEL: u32 = 0x08;

// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1369-1377, 6537
// ```c
// #define MDB_PS_MODIFY	1
// #define MDB_PS_ROOTONLY	2
// #define MDB_PS_FIRST	4
// #define MDB_PS_LAST		8
// #define MDB_SPLIT_REPLACE	MDB_APPENDDUP	/**< newkey is not new */
// #define MDB_NOSPILL	0x8000
// ```
const MDB_PS_MODIFY: u32 = 1;
const MDB_PS_ROOTONLY: u32 = 2;
const MDB_PS_FIRST: u32 = 4;
const MDB_PS_LAST: u32 = 8;
const MDB_SPLIT_REPLACE: u32 = MDB_APPENDDUP;
const MDB_NOSPILL: u32 = 0x8000;

/// The page size LMDB takes for a new environment: the operating system's page size
/// (`GET_PAGESIZE`), at most `MAX_PAGESIZE` (DECISIONS 8 of the makedb-v5 line).
// NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:395-397, 3969, 4400-4402
// ```c
// #define	GET_PAGESIZE(x)	((x) = sysconf(_SC_PAGE_SIZE))
// 	GET_PAGESIZE(e->me_os_psize);
// 		env->me_psize = env->me_os_psize;
// 		if (env->me_psize > MAX_PAGESIZE)
// 			env->me_psize = MAX_PAGESIZE;
// ```
// On Windows (`mdb.c:383-384`) `GET_PAGESIZE` is `GetSystemInfo(&si); (x) = si.dwPageSize`.
// wasm32 has no operating-system page; LOSAT takes 4096, the size of the platforms the Web
// app's users mostly build on (DECISIONS 8).
pub(crate) fn os_page_size() -> usize {
    #[cfg(unix)]
    let os = rustix::param::page_size();
    #[cfg(windows)]
    let os = {
        use windows_sys::Win32::System::SystemInformation::{GetSystemInfo, SYSTEM_INFO};
        // SAFETY: GetSystemInfo fills the SYSTEM_INFO it is given and cannot fail.
        let mut si: SYSTEM_INFO = unsafe { std::mem::zeroed() };
        unsafe { GetSystemInfo(&mut si) };
        si.dwPageSize as usize
    };
    #[cfg(not(any(unix, windows)))]
    let os = 4096usize;
    os.min(MAX_PAGESIZE)
}

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
        std::io::Error::from_raw_os_error(code).to_string()
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
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5174-5241
    // ```c
    // static int mdb_cmp_long(const MDB_val *a, const MDB_val *b)
    // { return (*(size_t *)a->mv_data < *(size_t *)b->mv_data) ? -1 :
    //		*(size_t *)a->mv_data > *(size_t *)b->mv_data; }
    // static int mdb_cmp_cint(const MDB_val *a, const MDB_val *b)
    // {	unsigned short *u, *c; int x;
    //	u = (unsigned short *) ((char *) a->mv_data + a->mv_size);
    //	c = (unsigned short *) ((char *) b->mv_data + a->mv_size);
    //	do { x = *--u - *--c; } while(!x && u > (unsigned short *)a->mv_data);
    //	return x; }
    // static int mdb_cmp_memn(const MDB_val *a, const MDB_val *b)
    // {	int diff; ssize_t len_diff; unsigned int len;
    //	len = a->mv_size; len_diff = (ssize_t) a->mv_size - (ssize_t) b->mv_size;
    //	if (len_diff > 0) { len = b->mv_size; len_diff = 1; }
    //	diff = memcmp(a->mv_data, b->mv_data, len);
    //	return diff ? diff : len_diff<0 ? -1 : len_diff; }
    // ```
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
}

/// A database handle (`MDB_dbi`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Dbi(usize);

impl Env {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:78-112
    // ```c++
    // 		m_Filename(fname), m_FileType(file_type),m_Env(lmdb::env::create()), m_Count(1), m_ReadOnly(read_only)
    // {
    // 	const MDB_dbi num_db(3);
    // 	m_Env.set_max_dbs(num_db);
    // 	...
    // 	else {
    // 		LOG_POST(Info <<"Initial Map Size: " << map_size);
    // 		/// map_size 0 means use lmdb default
    // 		if(map_size != 0) {
    // 			m_Env.set_mapsize(map_size);
    // 		}
    //         try {
    //             m_Env.open(m_Filename.c_str(), MDB_NOSUBDIR , 0664);
    // ```
    /// `CBlastLMDBManager::CBlastEnv` for writing: `mdb_env_create`,
    /// `mdb_env_set_maxdbs(3)`, `mdb_env_set_mapsize(map_size)` unless `map_size` is 0, and
    /// `mdb_env_open(path, MDB_NOSUBDIR, 0664)` on a new file. `psize` is the page size LMDB
    /// takes from the operating system ([`os_page_size`]); tests pass 4096.
    pub(crate) fn create(path: &Path, map_size: u64, psize: usize) -> Result<Env, LmdbError> {
        let io = |source| LmdbError::Io {
            path: path.to_path_buf(),
            source,
        };
        // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3951-3972, 4094-4101
        // ```c
        // 	e->me_maxreaders = DEFAULT_READERS;
        // 	e->me_maxdbs = e->me_numdbs = CORE_DBS;
        // 	...
        // 	env->me_maxdbs = dbs + CORE_DBS;
        // ```
        let maxdbs = 3 + CORE_DBS;
        // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4235-4325 (mdb_fopen)
        // ```c
        // 	fd = open(fname->mn_val, flags, mode);
        // ```
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
        };
        // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4060-4091 (mdb_env_set_mapsize)
        // ```c
        // 	env->me_mapsize = size;
        // 	if (env->me_psize)
        // 		env->me_maxpg = env->me_mapsize / env->me_psize;
        // ```
        if map_size != 0 {
            env.mapsize = map_size;
        }
        // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4993
        // ```c
        // 	env->me_dbxs[FREE_DBI].md_cmp = mdb_cmp_long; /* aligned MDB_INTEGERKEY */
        // ```
        env.dbx_cmp[FREE_DBI] = Cmp::Int;
        // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4400-4456 (mdb_env_open2)
        // ```c
        // 		newenv = 1;
        // 		env->me_psize = env->me_os_psize;
        // 		if (env->me_psize > MAX_PAGESIZE)
        // 			env->me_psize = MAX_PAGESIZE;
        // 		memset(&meta, 0, sizeof(meta));
        // 		mdb_env_init_meta0(env, &meta);
        // 		meta.mm_mapsize = DEFAULT_MAPSIZE;
        // 	...
        // 	if (!env->me_mapsize) {
        // 		env->me_mapsize = meta.mm_mapsize;
        // 	}
        // 	{
        // 		/* Make sure mapsize >= committed data size.  Even when using
        // 		 * mm_mapsize, which could be broken in old files (ITS#7789).
        // 		 */
        // 		size_t minsize = (meta.mm_last_pg + 1) * meta.mm_psize;
        // 		if (env->me_mapsize < minsize)
        // 			env->me_mapsize = minsize;
        // 	}
        // 	meta.mm_mapsize = env->me_mapsize;
        // 	if (newenv && !(flags & MDB_FIXEDMAP)) {
        // 		rc = mdb_env_init_meta(env, &meta);
        // 	...
        // 	env->me_maxfree_1pg = (env->me_psize - PAGEHDRSZ) / sizeof(pgno_t) - 1;
        // 	env->me_nodemax = (((env->me_psize - PAGEHDRSZ) / MDB_MINKEYS) & -2)
        // 		- sizeof(indx_t);
        // 	env->me_maxpg = env->me_mapsize / env->me_psize;
        // ```
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
        Ok(env)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3740-3753 (mdb_env_init_meta0)
    // ```c
    // 	meta->mm_magic = MDB_MAGIC;
    // 	meta->mm_version = MDB_DATA_VERSION;
    // 	meta->mm_mapsize = env->me_mapsize;
    // 	meta->mm_psize = env->me_psize;
    // 	meta->mm_last_pg = NUM_METAS-1;
    // 	meta->mm_flags = env->me_flags & 0xffff;
    // 	meta->mm_flags |= MDB_INTEGERKEY; /* this is mm_dbs[FREE_DBI].md_flags */
    // 	meta->mm_dbs[FREE_DBI].md_root = P_INVALID;
    // 	meta->mm_dbs[MAIN_DBI].md_root = P_INVALID;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3759-3807 (mdb_env_init_meta)
    // ```c
    // 	p = calloc(NUM_METAS, psize);
    // 	p->mp_pgno = 0;
    // 	p->mp_flags = P_META;
    // 	*(MDB_meta *)METADATA(p) = *meta;
    // 	q = (MDB_page *)((char *)p + psize);
    // 	q->mp_pgno = 1;
    // 	q->mp_flags = P_META;
    // 	*(MDB_meta *)METADATA(q) = *meta;
    // 	DO_PWRITE(rc, env->me_fd, p, psize * NUM_METAS, len, 0);
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3944-3948 (mdb_env_pick_meta)
    // ```c
    // 	return metas[ metas[0]->mm_txnid < metas[1]->mm_txnid ];
    // ```
    fn pick_meta(&self) -> &Meta {
        &self.metas[(self.metas[0].txnid < self.metas[1].txnid) as usize]
    }

    /// `mdb_env_stat(env).ms_psize`.
    pub(crate) fn psize(&self) -> usize {
        self.psize
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:9688-9705 (mdb_env_info)
    // ```c
    // 	meta = mdb_env_pick_meta(env);
    // 	arg->me_mapaddr = meta->mm_address;
    // 	arg->me_last_pgno = meta->mm_last_pg;
    // 	arg->me_last_txnid = meta->mm_txnid;
    // 	arg->me_mapsize = env->me_mapsize;
    // ```
    /// `mdb_env_info(env).me_last_pgno`.
    pub(crate) fn last_pgno(&self) -> u64 {
        self.pick_meta().last_pg
    }

    /// `mdb_env_info(env).me_mapsize`.
    pub(crate) fn mapsize(&self) -> u64 {
        self.mapsize
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:4060-4091 (mdb_env_set_mapsize)
    // ```c
    // 	if (env->me_map) {
    // 		...
    // 		if (env->me_txn)
    // 			return EINVAL;
    // 		meta = mdb_env_pick_meta(env);
    // 		if (!size)
    // 			size = meta->mm_mapsize;
    // 		{
    // 			/* Silently round up to minimum if the size is too small */
    // 			size_t minsize = (meta->mm_last_pg + 1) * env->me_psize;
    // 			if (size < minsize)
    // 				size = minsize;
    // 		}
    // 		munmap(env->me_map, env->me_mapsize);
    // 		env->me_mapsize = size;
    // 		...
    // 	}
    // 	env->me_mapsize = size;
    // 	if (env->me_psize)
    // 		env->me_maxpg = env->me_mapsize / env->me_psize;
    // ```
    /// `mdb_env_set_mapsize` on the open environment (the map itself is not kept).
    pub(crate) fn set_mapsize(&mut self, size: u64) -> Result<(), LmdbError> {
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            return Err(LmdbError::Mdb {
                func: "mdb_env_set_mapsize",
                code: EINVAL,
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3979-4003 (mdb_env_map)
    // ```c
    // #ifdef _WIN32
    // 	...
    // 	} else {
    // 		msize = env->me_mapsize;
    // 		sizelo = msize & 0xffffffff;
    // 		sizehi = msize >> 16 >> 16; /* only needed on Win64 */
    //
    // 		/* Windows won't create mappings for zero length files.
    // 		 * and won't map more than the file size.
    // 		 * Just set the maxsize right now.
    // 		 */
    // 		if (!(flags & MDB_WRITEMAP) && (SetFilePointer(env->me_fd, sizelo, &sizehi, 0) != (DWORD)sizelo
    // 			|| !SetEndOfFile(env->me_fd)
    // 			|| SetFilePointer(env->me_fd, 0, NULL, 0) != 0))
    // 			return ErrCode();
    // 	}
    // ```
    // On Windows LMDB sets the length of the data file to the map size whenever it maps the
    // file: when the environment opens (after the meta pages are written, mdb.c:4421-4436)
    // and when `mdb_env_set_mapsize` maps it again with a new size (mdb.c:4083-4087). The
    // pages are written inside that length, so a file of makeblastdb on Windows is as long
    // as the last map size (500000 bytes by default, `DEFAULT_LMDB_MAP_SIZE`). Elsewhere
    // the file is not extended (`MDB_WRITEMAP` is not set). Not checked against a
    // makeblastdb.exe.
    /// The length `mdb_env_map` gives the data file: the map size on Windows.
    fn map_file_size(&mut self) -> std::io::Result<()> {
        #[cfg(windows)]
        {
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2747-2809 (mdb_txn_renew0, write txn)
    // ```c
    // 		if (ti) {
    // 			if (LOCK_MUTEX(rc, env, env->me_wmutex))
    // 				return rc;
    // 			txn->mt_txnid = ti->mti_txnid;
    // 			meta = env->me_metas[txn->mt_txnid & 1];
    // 		} else { ... }
    // 		txn->mt_txnid++;
    // 		txn->mt_child = NULL;
    // 		txn->mt_loose_pgs = NULL;
    // 		txn->mt_loose_count = 0;
    // 		txn->mt_dirty_room = MDB_IDL_UM_MAX;
    // 		txn->mt_u.dirty_list = env->me_dirty_list;
    // 		txn->mt_u.dirty_list[0].mid = 0;
    // 		txn->mt_free_pgs = env->me_free_pgs;
    // 		txn->mt_free_pgs[0] = 0;
    // 		txn->mt_spill_pgs = NULL;
    // 		env->me_txn = txn;
    // 	...
    // 	/* Copy the DB info and flags */
    // 	memcpy(txn->mt_dbs, meta->mm_dbs, CORE_DBS * sizeof(MDB_db));
    // 	/* Moved to here to avoid a data race in read TXNs */
    // 	txn->mt_next_pgno = meta->mm_last_pg+1;
    // 	txn->mt_flags = flags;
    // 	/* Setup db info */
    // 	txn->mt_numdbs = env->me_numdbs;
    // 	for (i=CORE_DBS; i<txn->mt_numdbs; i++) {
    // 		x = env->me_dbflags[i];
    // 		txn->mt_dbs[i].md_flags = x & PERSISTENT_FLAGS;
    // 		txn->mt_dbflags[i] = (x & MDB_VALID) ? DB_VALID|DB_USRVALID|DB_STALE : 0;
    // 	}
    // 	txn->mt_dbflags[MAIN_DBI] = DB_VALID|DB_USRVALID;
    // 	txn->mt_dbflags[FREE_DBI] = DB_VALID;
    // 	if (env->me_flags & MDB_FATAL_ERROR) {
    // 		DPUTS("environment had fatal error, must shutdown!");
    // 		rc = MDB_PANIC;
    // 	} else if (env->me_maxpg < txn->mt_next_pgno) {
    // 		rc = MDB_MAP_RESIZED;
    // 	} else {
    // 		return MDB_SUCCESS;
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2943-2970 (mdb_dbis_update)
    // ```c
    // 	for (i = n; --i >= CORE_DBS;) {
    // 		if (tdbflags[i] & DB_NEW) {
    // 			if (keep) {
    // 				env->me_dbflags[i] = txn->mt_dbs[i].md_flags | MDB_VALID;
    // 			} else {
    // 				char *ptr = env->me_dbxs[i].md_name.mv_data;
    // 				if (ptr) {
    // 					env->me_dbxs[i].md_name.mv_data = NULL;
    // 					env->me_dbxs[i].md_name.mv_size = 0;
    // 					env->me_dbflags[i] = 0;
    // 					env->me_dbiseqs[i]++;
    // 					free(ptr);
    // 				}
    // 			}
    // 		}
    // 	}
    // 	if (keep && env->me_numdbs < n)
    // 		env->me_numdbs = n;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2976-3045 (mdb_txn_end)
    // ```c
    // 	mdb_dbis_update(txn, mode & MDB_END_UPDATE);
    // 	...
    // 	} else if (!F_ISSET(txn->mt_flags, MDB_TXN_FINISHED)) {
    // 		pgno_t *pghead = env->me_pghead;
    // 		if (!(mode & MDB_END_UPDATE)) /* !(already closed cursors) */
    // 			mdb_cursors_close(txn, 0);
    // 		if (!(env->me_flags & MDB_WRITEMAP)) {
    // 			mdb_dlist_free(txn);
    // 		}
    // 		txn->mt_numdbs = 0;
    // 		txn->mt_flags = MDB_TXN_FINISHED;
    // 		if (!txn->mt_parent) {
    // 			mdb_midl_shrink(&txn->mt_free_pgs);
    // 			env->me_free_pgs = txn->mt_free_pgs;
    // 			/* me_pgstate: */
    // 			env->me_pghead = NULL;
    // 			env->me_pglast = 0;
    // 			env->me_txn = NULL;
    // 	...
    // 		mdb_midl_free(pghead);
    // ```
    fn txn_end(&mut self, mode: u32) {
        self.dbis_update(mode & MDB_END_UPDATE != 0);
        if self.t_flags & MDB_TXN_FINISHED == 0 {
            // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1831-1843 (mdb_dlist_free)
            // ```c
            // 	for (i = 1; i <= n; i++) {
            // 		mdb_dpage_free(env, dl[i].mptr);
            // 	}
            // 	dl[0].mid = 0;
            // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3448-3674 (mdb_txn_commit)
    // ```c
    // 	if (txn->mt_flags & (MDB_TXN_FINISHED|MDB_TXN_ERROR)) {
    // 		...
    // 		rc = MDB_BAD_TXN;
    // 		goto fail;
    // 	}
    // 	...
    // 	mdb_cursors_close(txn, 0);
    // 	if (!txn->mt_u.dirty_list[0].mid &&
    // 		!(txn->mt_flags & (MDB_TXN_DIRTY|MDB_TXN_SPILLS)))
    // 		goto done;
    // 	/* Update DB root pointers */
    // 	if (txn->mt_numdbs > CORE_DBS) {
    // 		MDB_cursor mc;
    // 		MDB_dbi i;
    // 		MDB_val data;
    // 		data.mv_size = sizeof(MDB_db);
    // 		mdb_cursor_init(&mc, txn, MAIN_DBI, NULL);
    // 		for (i = CORE_DBS; i < txn->mt_numdbs; i++) {
    // 			if (txn->mt_dbflags[i] & DB_DIRTY) {
    // 				...
    // 				data.mv_data = &txn->mt_dbs[i];
    // 				rc = mdb_cursor_put(&mc, &txn->mt_dbxs[i].md_name, &data,
    // 					F_SUBDATA);
    // 				if (rc)
    // 					goto fail;
    // 			}
    // 		}
    // 	}
    // 	rc = mdb_freelist_save(txn);
    // 	if (rc)
    // 		goto fail;
    // 	mdb_midl_free(env->me_pghead);
    // 	env->me_pghead = NULL;
    // 	mdb_midl_shrink(&txn->mt_free_pgs);
    // 	if ((rc = mdb_page_flush(txn, 0)) ||
    // 		(rc = mdb_env_sync(env, 0)) ||
    // 		(rc = mdb_env_write_meta(txn)))
    // 		goto fail;
    // 	end_mode = MDB_END_COMMITTED|MDB_END_UPDATE;
    // done:
    // 	mdb_txn_end(txn, end_mode);
    // 	return MDB_SUCCESS;
    // fail:
    // 	mdb_txn_abort(txn);
    // 	return rc;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3812-3937 (mdb_env_write_meta, no
    // MDB_WRITEMAP)
    // ```c
    // 	toggle = txn->mt_txnid & 1;
    // 	mp = env->me_metas[toggle];
    // 	mapsize = env->me_metas[toggle ^ 1]->mm_mapsize;
    // 	/* Persist any increases of mapsize config */
    // 	if (mapsize < env->me_mapsize)
    // 		mapsize = env->me_mapsize;
    // 	...
    // 	meta.mm_mapsize = mapsize;
    // 	meta.mm_dbs[FREE_DBI] = txn->mt_dbs[FREE_DBI];
    // 	meta.mm_dbs[MAIN_DBI] = txn->mt_dbs[MAIN_DBI];
    // 	meta.mm_last_pg = txn->mt_next_pgno - 1;
    // 	meta.mm_txnid = txn->mt_txnid;
    // 	off = offsetof(MDB_meta, mm_mapsize);
    // 	ptr = (char *)&meta + off;
    // 	len = sizeof(MDB_meta) - off;
    // 	off += (char *)mp - env->me_map;
    // 	...
    // 	rc = pwrite(mfd, ptr, len, off);
    // 	...
    // 	if (env->me_txns)
    // 		env->me_txns->mti_txnid = txn->mt_txnid;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1768-1805 (mdb_page_malloc)
    // ```c
    // 	MDB_page *ret = env->me_dpages;
    // 	size_t psize = env->me_psize, sz = psize, off;
    // 	if (num == 1) {
    // 		if (ret) {
    // 			VGMEMP_ALLOC(env, ret, sz);
    // 			VGMEMP_DEFINED(ret, sizeof(ret->mp_next));
    // 			env->me_dpages = ret->mp_next;
    // 			return ret;
    // 		}
    // 		psize -= off = PAGEHDRSZ;
    // 	} else {
    // 		sz *= num;
    // 		off = sz - psize;
    // 	}
    // 	if ((ret = malloc(sz)) != NULL) {
    // 		VGMEMP_ALLOC(env, ret, sz);
    // 		if (!(env->me_flags & MDB_NOMEMINIT)) {
    // 			memset((char *)ret + off, 0, psize);
    // 			ret->mp_pad = 0;
    // 		}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1807-1829 (mdb_page_free,
    // mdb_dpage_free)
    // ```c
    // 	mp->mp_next = env->me_dpages;
    // 	VGMEMP_FREE(env, mp);
    // 	env->me_dpages = mp;
    // ...
    // 	if (!IS_OVERFLOW(dp) || dp->mp_pages == 1) {
    // 		mdb_page_free(env, dp);
    // 	} else {
    // 		/* large pages just get freed directly */
    // 		VGMEMP_FREE(env, dp);
    // 		free(dp);
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5428-5479 (mdb_page_get)
    // ```c
    // 	if (! (txn->mt_flags & (MDB_TXN_RDONLY|MDB_TXN_WRITEMAP))) {
    // 		...
    // 			if (dl[0].mid) {
    // 				unsigned x = mdb_mid2l_search(dl, pgno);
    // 				if (x <= dl[0].mid && dl[x].mid == pgno) {
    // 					p = dl[x].mptr;
    // 					goto done;
    // 				}
    // 			}
    // 	...
    // 	if (pgno < txn->mt_next_pgno) {
    // 		level = 0;
    // 		p = (MDB_page *)(env->me_map + env->me_psize * pgno);
    // 	} else {
    // 		DPRINTF(("page %"Z"u not found", pgno));
    // 		txn->mt_flags |= MDB_TXN_ERROR;
    // 		return MDB_PAGE_NOTFOUND;
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2099-2116 (mdb_find_oldest)
    // ```c
    // 	txnid_t mr, oldest = txn->mt_txnid - 1;
    // 	if (txn->mt_env->me_txns) {
    // 		MDB_reader *r = txn->mt_env->me_txns->mti_readers;
    // 		for (i = txn->mt_env->me_txns->mti_numreaders; --i >= 0; ) {
    // ```
    // makeblastdb opens no read transaction on the environment it writes: the reader table
    // is empty.
    fn find_oldest(&self) -> u64 {
        self.t_txnid - 1
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2118-2135 (mdb_page_dirty)
    // ```c
    // 	mid.mid = mp->mp_pgno;
    // 	mid.mptr = mp;
    // 	rc = insert(txn->mt_u.dirty_list, &mid);
    // 	mdb_tassert(txn, rc == 0);
    // 	txn->mt_dirty_room--;
    // ```
    fn page_dirty(&mut self, pgno: u64, np: Box<[u8]>) {
        let x = match self.t_dirty.binary_search_by(|e| e.0.cmp(&pgno)) {
            Ok(_) => panic!("LMDB writer: page {pgno} dirty twice"),
            Err(x) => x,
        };
        self.t_dirty.insert(x, (pgno, np));
        self.clean.remove(&pgno);
        self.t_dirty_room -= 1;
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2137-2322 (mdb_page_alloc)
    // ```c
    // 	int rc, retry = num * 60;
    // 	MDB_txn *txn = mc->mc_txn;
    // 	MDB_env *env = txn->mt_env;
    // 	pgno_t pgno, *mop = env->me_pghead;
    // 	unsigned i, j, mop_len = mop ? mop[0] : 0, n2 = num-1;
    // 	...
    // 	for (op = MDB_FIRST;; op = MDB_NEXT) {
    // 		MDB_val key, data;
    // 		MDB_node *leaf;
    // 		pgno_t *idl;
    // 		/* Seek a big enough contiguous page range. Prefer
    // 		 * pages at the tail, just truncating the list.
    // 		 */
    // 		if (mop_len > n2) {
    // 			i = mop_len;
    // 			do {
    // 				pgno = mop[i];
    // 				if (mop[i-n2] == pgno+n2)
    // 					goto search_done;
    // 			} while (--i > n2);
    // 			if (--retry < 0)
    // 				break;
    // 		}
    // 		if (op == MDB_FIRST) {	/* 1st iteration */
    // 			/* Prepare to fetch more and coalesce */
    // 			last = env->me_pglast;
    // 			oldest = env->me_pgoldest;
    // 			mdb_cursor_init(&m2, txn, FREE_DBI, NULL);
    // 			if (last) {
    // 				op = MDB_SET_RANGE;
    // 				key.mv_data = &last; /* will look up last+1 */
    // 				key.mv_size = sizeof(last);
    // 			}
    // 			if (Paranoid && mc->mc_dbi == FREE_DBI)
    // 				retry = -1;
    // 		}
    // 		if (Paranoid && retry < 0 && mop_len)
    // 			break;
    // 		last++;
    // 		/* Do not fetch more if the record will be too recent */
    // 		if (oldest <= last) {
    // 			if (!found_old) {
    // 				oldest = mdb_find_oldest(txn);
    // 				env->me_pgoldest = oldest;
    // 				found_old = 1;
    // 			}
    // 			if (oldest <= last)
    // 				break;
    // 		}
    // 		rc = mdb_cursor_get(&m2, &key, NULL, op);
    // 		if (rc) {
    // 			if (rc == MDB_NOTFOUND)
    // 				break;
    // 			goto fail;
    // 		}
    // 		last = *(txnid_t*)key.mv_data;
    // 		if (oldest <= last) {
    // 			... (same as above)
    // 				break;
    // 		}
    // 		np = m2.mc_pg[m2.mc_top];
    // 		leaf = NODEPTR(np, m2.mc_ki[m2.mc_top]);
    // 		if ((rc = mdb_node_read(&m2, leaf, &data)) != MDB_SUCCESS)
    // 			goto fail;
    // 		idl = (MDB_ID *) data.mv_data;
    // 		i = idl[0];
    // 		if (!mop) {
    // 			if (!(env->me_pghead = mop = mdb_midl_alloc(i))) {
    // 		...
    // 		env->me_pglast = last;
    // 		/* Merge in descending sorted order */
    // 		mdb_midl_xmerge(mop, idl);
    // 		mop_len = mop[0];
    // 	}
    // 	/* Use new pages from the map when nothing suitable in the freeDB */
    // 	i = 0;
    // 	pgno = txn->mt_next_pgno;
    // 	if (pgno + num >= env->me_maxpg) {
    // 			DPUTS("DB size maxed out");
    // 			rc = MDB_MAP_FULL;
    // 			goto fail;
    // 	}
    // search_done:
    // 	...
    // 		if (!(np = mdb_page_malloc(txn, num))) {
    // 	...
    // 	if (i) {
    // 		mop[0] = mop_len -= num;
    // 		/* Move any stragglers down */
    // 		for (j = i-num; j < mop_len; )
    // 			mop[++j] = mop[++i];
    // 	} else {
    // 		txn->mt_next_pgno = pgno + num;
    // 	}
    // 	np->mp_pgno = pgno;
    // 	mdb_page_dirty(txn, np);
    // 	*mp = np;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2324-2346 (mdb_page_copy)
    // ```c
    // 	enum { Align = sizeof(pgno_t) };
    // 	indx_t upper = src->mp_upper, lower = src->mp_lower, unused = upper-lower;
    // 	/* If page isn't full, just copy the used portion. Adjust
    // 	 * alignment so memcpy may copy words instead of bytes.
    // 	 */
    // 	if ((unused &= -Align) && !IS_LEAF2(src)) {
    // 		upper = (upper + PAGEBASE) & -Align;
    // 		memcpy(dst, src, (lower + PAGEBASE + (Align-1)) & -Align);
    // 		memcpy((pgno_t *)((char *)dst+upper), (pgno_t *)((char *)src+upper),
    // 			psize - upper);
    // 	} else {
    // 		memcpy(dst, src, psize - unused);
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3295-3445 (mdb_page_flush, keep 0, no
    // MDB_WRITEMAP)
    // ```c
    // 	/* Write the pages */
    // 	for (;;) {
    // 		if (++i <= pagecount) {
    // 			dp = dl[i].mptr;
    // 			...
    // 			pgno = dl[i].mid;
    // 			/* clear dirty flag */
    // 			dp->mp_flags &= ~P_DIRTY;
    // 			pos = pgno * psize;
    // 			size = psize;
    // 			if (IS_OVERFLOW(dp)) size *= dp->mp_pages;
    // 		}
    // 		... (pwrite/writev of runs of up to MDB_COMMIT_PAGES contiguous pages)
    // 	}
    // 	for (i = keep; ++i <= pagecount; ) {
    // 		dp = dl[i].mptr;
    // 		...
    // 		mdb_dpage_free(env, dp);
    // 	}
    // done:
    // 	i--;
    // 	txn->mt_dirty_room += i - j;
    // 	dl[0].mid = j;
    // ```
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

    /// Close the environment (`mdb_env_close`, which writes nothing) and remove the lock
    /// file as `CWriteDB_LMDB`/`CWriteDB_TaxID` do after it.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72
    // ```c++
    // CWriteDB_LMDB::~CWriteDB_LMDB()
    // {
    // 	x_CreateOidToSeqidsLookupFile();
    // 	x_CommitTransaction();
    //     CBlastLMDBManager::GetInstance().CloseEnv(m_Db);
    //     CFile(m_Db+"-lock").Remove();
    // ```
    // LMDB's lock file holds only the reader table and the process id; LOSAT does not
    // create it, and removes one left at the path so that the files are as makeblastdb
    // leaves them.
    pub(crate) fn close(mut self) -> Result<(), LmdbError> {
        self.txn_abort();
        let flushed = self.file.flush();
        let mut lock = self.path.clone().into_os_string();
        lock.push("-lock");
        let lock = PathBuf::from(lock);
        if lock.exists() {
            let _ = std::fs::remove_file(&lock);
        }
        flushed.map_err(|source| LmdbError::Io {
            path: self.path.clone(),
            source,
        })
    }
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
