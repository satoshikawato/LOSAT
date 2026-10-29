//! The LMDB files of a version 5 database and their lookup files: `CWriteDB_LMDB` (the
//! `.pdb`/`.ndb` file with the volume tables and `acc2oid`, and the `.pos`/`.nos` file) and
//! `CWriteDB_TaxID` (the `.ptf`/`.ntf` file with `taxid2offset`, and the `.pot`/`.not` and
//! `.pto`/`.nto` files), `c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp`.
//!
//! The LMDB files are written by the port of NCBI's LMDB in [`super::mdb`].

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use super::mdb::{
    Dbi, Env, LmdbError, MDB_APPENDDUP, MDB_CREATE, MDB_DUPFIXED, MDB_DUPSORT, MDB_INTEGERKEY,
};
use super::WriteDbError;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_lmdb.hpp:51-57
// ```c++
// #ifdef NCBI_OS_MSWIN
// #define DEFAULT_LMDB_MAP_SIZE 500000
// #define DEFAULT_TAXID_MAP_SIZE 500000
// #else
// #define DEFAULT_LMDB_MAP_SIZE 700000000
// #define DEFAULT_TAXID_MAP_SIZE 300000000
// #endif
// ```
#[cfg(windows)]
pub(crate) const DEFAULT_LMDB_MAP_SIZE: u64 = 500_000;
#[cfg(windows)]
pub(crate) const DEFAULT_TAXID_MAP_SIZE: u64 = 500_000;
#[cfg(not(windows))]
pub(crate) const DEFAULT_LMDB_MAP_SIZE: u64 = 700_000_000;
#[cfg(not(windows))]
pub(crate) const DEFAULT_TAXID_MAP_SIZE: u64 = 300_000_000;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:46-48
// ```c++
// #define DEFAULT_MAX_ENTRY_PER_TXN 40000
// #define DEFAULT_MIN_SPLIT_SORT_SIZE 500000000
// #define DEFAULT_MIN_SPLIT_CHUNK_SIZE 25000000
// ```
const DEFAULT_MAX_ENTRY_PER_TXN: u32 = 40_000;
const DEFAULT_MIN_SPLIT_SORT_SIZE: u32 = 500_000_000;
const DEFAULT_MIN_SPLIT_CHUNK_SIZE: u32 = 25_000_000;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdb_lmdb.hpp:196-199
// ```c++
// const string volinfo_str = "volinfo";
// const string volname_str = "volname";
// const string acc2oid_str = "acc2oid";
// const string taxid2offset_str = "taxid2offset";
// ```
const VOLINFO_STR: &str = "volinfo";
const VOLNAME_STR: &str = "volname";
const ACC2OID_STR: &str = "acc2oid";
const TAXID2OFFSET_STR: &str = "taxid2offset";

/// The environment variables of NCBI's LMDB layer, as makeblastdb reads them (`getenv`),
/// unparsed: each is parsed where NCBI parses it, with NCBI's rules.
#[derive(Clone, Debug, Default)]
pub(crate) struct LmdbEnvSettings {
    /// `MAX_LMDB_TXN_ENTRY` (both writers' constructors).
    pub max_lmdb_txn_entry: Option<String>,
    /// `LMDB_MIN_SPLIT_SIZE` (`CWriteDB_LMDB::x_CommitTransaction`).
    pub lmdb_min_split_size: Option<String>,
    /// `LMDB_SPLIT_CHUNK_SIZE` (`CWriteDB_LMDB::x_CommitTransaction`).
    pub lmdb_split_chunk_size: Option<String>,
}

impl LmdbEnvSettings {
    /// The variables of this process's environment.
    pub(crate) fn from_env() -> Self {
        let get = |name: &str| std::env::var_os(name).map(|v| v.to_string_lossy().into_owned());
        LmdbEnvSettings {
            max_lmdb_txn_entry: get("MAX_LMDB_TXN_ENTRY"),
            lmdb_min_split_size: get("LMDB_MIN_SPLIT_SIZE"),
            lmdb_split_chunk_size: get("LMDB_SPLIT_CHUNK_SIZE"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1136-1144
// ```c++
//     if(m_DbVersion == eBDB_Version5 && m_Lmdbdb.Empty()) {
//         const string lmdb_fname_w_path = BuildLMDBFileName(m_Dbname, m_Protein);
//         Uint8 map_size = 0;
//         char* map_sz_str = getenv("BLASTDB_LMDB_MAP_SIZE");
//         if (map_sz_str) {
//         	map_size = NStr::StringToUInt8(map_sz_str);
//         }
//         if(map_size > 0){
// ```
/// `BLASTDB_LMDB_MAP_SIZE` as `CWriteDB_Impl::x_Publish` reads it: 0 when unset (the
/// writers then take their default map sizes).
pub(crate) fn lmdb_map_size_from_env() -> Result<u64, WriteDbError> {
    match std::env::var_os("BLASTDB_LMDB_MAP_SIZE") {
        None => Ok(0),
        Some(v) => string_to_uint8(&v.to_string_lossy()).map_err(WriteDbError::Other),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-874 (NStr::StringToInt8, flags 0,
// base 10)
// ```c++
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     ...
//     while (char ch = str[pos]) {
//         ...
//         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
//             break;
//         }
//         // Overflow check
//         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
//             S2N_CONVERT_ERROR_OVERFLOW(Int8);
//         }
//         n *= base;
//         n += delta;
//         pos++;
//     }
//     // Last checks
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     ...
//     n = sign ? -n : n;
//     CHECK_ENDPTR(Int8);
// ```
// StringToUInt8 (876-985) takes no '-'; StringToInt and StringToUInt (635-660) convert with
// the 8-byte functions and fail with ", overflow" outside int / unsigned int.
fn parse_digits(s: &str, allow_minus: bool, ty: &str) -> Result<(bool, u128), String> {
    let invalid = || format!("Cannot convert string '{s}' to {ty}");
    let b = s.as_bytes();
    let mut pos = 0;
    let mut neg = false;
    if pos < b.len() && (b[pos] == b'+' || (allow_minus && b[pos] == b'-')) {
        neg = b[pos] == b'-';
        pos += 1;
    }
    let start = pos;
    let mut n: u128 = 0;
    while pos < b.len() && b[pos].is_ascii_digit() {
        n = n * 10 + (b[pos] - b'0') as u128;
        if n > u64::MAX as u128 + 1 {
            return Err(format!("{}, overflow", invalid()));
        }
        pos += 1;
    }
    if pos == start || pos != b.len() {
        return Err(invalid());
    }
    Ok((neg, n))
}

fn string_to_int8(s: &str) -> Result<i64, String> {
    let (neg, n) = parse_digits(s, true, "Int8")?;
    let v: i128 = if neg { -(n as i128) } else { n as i128 };
    if v < i64::MIN as i128 || v > i64::MAX as i128 {
        return Err(format!("Cannot convert string '{s}' to Int8, overflow"));
    }
    Ok(v as i64)
}

fn string_to_int(s: &str) -> Result<i32, String> {
    let v = string_to_int8(s)?;
    i32::try_from(v).map_err(|_| format!("Cannot convert string '{s}' to int, overflow"))
}

fn string_to_uint8(s: &str) -> Result<u64, String> {
    let (_, n) = parse_digits(s, false, "Uint8")?;
    u64::try_from(n).map_err(|_| format!("Cannot convert string '{s}' to Uint8, overflow"))
}

fn string_to_uint(s: &str) -> Result<u32, String> {
    let v = string_to_uint8(s)?;
    u32::try_from(v).map_err(|_| format!("Cannot convert string '{s}' to unsigned int, overflow"))
}

/// `m_MaxEntryPerTxn` from `MAX_LMDB_TXN_ENTRY`.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:58-62
// ```c++
// 	char* max_entry_str = getenv("MAX_LMDB_TXN_ENTRY");
// 	if (max_entry_str) {
// 		m_MaxEntryPerTxn = NStr::StringToInt(max_entry_str);
// 	    _TRACE("DEBUG: LMDB TXN Max Entry " << m_MaxEntryPerTxn);
//     }
// ```
// `m_MaxEntryPerTxn` is an `unsigned int`: a negative value wraps.
fn max_entry_per_txn(settings: &LmdbEnvSettings) -> Result<u32, WriteDbError> {
    match &settings.max_lmdb_txn_entry {
        None => Ok(DEFAULT_MAX_ENTRY_PER_TXN),
        Some(v) => string_to_int(v)
            .map(|n| n as u32)
            .map_err(WriteDbError::Other),
    }
}

fn lmdb_error(e: LmdbError) -> WriteDbError {
    match e {
        LmdbError::Io { path, source } => WriteDbError::Io { path, source },
        other => WriteDbError::Other(other.to_string()),
    }
}

/// `GetFileNameFromExistingLMDBFile`: the LMDB file name with its last two letters replaced.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:763-788
// ```c++
// 	string filename (lmdb_filename, 0, lmdb_filename.size() - 2);
// 	switch (file_type) {
// 		case eLMDB : filename += "db"; break;
// 		case eOid2SeqIds : filename +="os"; break;
// 		case eOid2TaxIds : filename +="ot"; break;
// 		case eTaxId2Offsets : filename += "tf"; break;
// 		case eTaxId2Oids : filename += "to"; break;
// ```
pub(crate) fn lmdb_file_name(lmdb_filename: &Path, suffix: &str) -> PathBuf {
    let mut s = lmdb_filename
        .as_os_str()
        .to_owned()
        .into_string()
        .unwrap_or_default();
    s.truncate(s.len().saturating_sub(2));
    s.push_str(suffix);
    PathBuf::from(s)
}

/// The bytes up to the first NUL (`strlen(id.c_str())`).
fn c_str(b: &[u8]) -> &[u8] {
    match b.iter().position(|&c| c == 0) {
        Some(n) => &b[..n],
        None => b,
    }
}

/// Run the batches of `x_CommitTransaction`: `[i, j)` with
/// `j = i + m_MaxEntryPerTxn` in `unsigned int` arithmetic, clipped to the list size.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:312-321
// ```c++
//     unsigned int j=0;
//     while (j < m_list.size()){
//     	lmdb::txn txn = lmdb::txn::begin(m_Env);
//     	...
//     	unsigned int i = j;
//     	j= i+m_MaxEntryPerTxn;
//     	if(j > m_list.size()) {
//     		j = m_list.size();
//     	}
// ```
// `MAX_LMDB_TXN_ENTRY=0` makes makeblastdb loop for ever (no batch ever ends the list), and a
// batch end that wraps past 4294967295 sends its loop back over entries already put: LOSAT
// refuses both.
fn txn_batches(len: usize, max: u32) -> Result<Vec<(usize, usize)>, WriteDbError> {
    if len > u32::MAX as usize {
        return Err(WriteDbError::Other(
            "LMDB writer: more than 4294967295 entries in one transaction list".to_string(),
        ));
    }
    if max == 0 && len > 0 {
        return Err(WriteDbError::Other(
            "MAX_LMDB_TXN_ENTRY=0: makeblastdb never finishes its LMDB transactions".to_string(),
        ));
    }
    let mut out = Vec::new();
    let mut j: u32 = 0;
    while (j as usize) < len {
        let i = j;
        j = i.wrapping_add(max);
        if j < i && j as usize <= len {
            return Err(WriteDbError::Other(format!(
                "MAX_LMDB_TXN_ENTRY={max}: the end of an LMDB transaction wraps past 4294967295"
            )));
        }
        if j as usize > len {
            j = len as u32;
        }
        out.push((i as usize, j as usize));
    }
    Ok(out)
}

/// One element of `CWriteDB_LMDB::m_list` (`SKeyValuePair`); the id is in `WriteDbLmdb::ids`.
#[derive(Clone, Copy, Debug)]
struct KeyValuePair {
    off: u64,
    len: u32,
    oid: u32,
    save_to_oid_list: bool,
}

/// `CWriteDB_LMDB`: the `.pdb`/`.ndb` file and the `.pos`/`.nos` file of a database.
pub(crate) struct WriteDbLmdb {
    db: PathBuf,
    env: Env,
    max_entry_per_txn: u32,
    total_ids_length: usize,
    settings: LmdbEnvSettings,
    ids: Vec<u8>,
    list: Vec<KeyValuePair>,
}

impl WriteDbLmdb {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:51-63
    // ```c++
    // CWriteDB_LMDB::CWriteDB_LMDB(const string& dbname,  Uint8 map_size, Uint8 capacity): m_Db(dbname),
    //                              m_Env(CBlastLMDBManager::GetInstance().GetWriteEnv(dbname, map_size)),
    //                              m_ListCapacity(capacity),
    //                              m_MaxEntryPerTxn(DEFAULT_MAX_ENTRY_PER_TXN),
    //                              m_TotalIdsLength(0)
    // {
    // 	m_list.reserve(m_ListCapacity);
    // 	char* max_entry_str = getenv("MAX_LMDB_TXN_ENTRY");
    // ```
    /// Create the LMDB file `dbname` (`<db>.pdb`/`<db>.ndb`) with the map size `map_size`
    /// and the page size `page_size` ([`super::mdb::os_page_size`] for makeblastdb's bytes on
    /// this platform).
    pub(crate) fn new(
        dbname: &Path,
        map_size: u64,
        page_size: usize,
        settings: &LmdbEnvSettings,
    ) -> Result<Self, WriteDbError> {
        let env = Env::create(dbname, map_size, page_size).map_err(lmdb_error)?;
        let max_entry_per_txn = max_entry_per_txn(settings)?;
        Ok(WriteDbLmdb {
            db: dbname.to_path_buf(),
            env,
            max_entry_per_txn,
            total_ids_length: 0,
            settings: settings.clone(),
            ids: Vec::new(),
            list: Vec::new(),
        })
    }

    /// Add one `SKeyValuePair` to `m_list`, as `CWriteDB_LMDB::x_InsertEntry` makes them from
    /// a Seq-id of OID `oid` (`writedb_lmdb.cpp:119-196`): the key string, and whether the id
    /// goes into the `.pos` file (`saveToOidList`). Called in OID order.
    pub(crate) fn insert_entry(&mut self, id: &[u8], oid: u32, save_to_oid_list: bool) {
        let off = self.ids.len() as u64;
        self.ids.extend_from_slice(id);
        self.list.push(KeyValuePair {
            off,
            len: id.len() as u32,
            oid,
            save_to_oid_list,
        });
    }

    fn id(&self, kv: &KeyValuePair) -> &[u8] {
        &self.ids[kv.off as usize..kv.off as usize + kv.len as usize]
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:74-94
    // ```c++
    // void CWriteDB_LMDB::InsertVolumesInfo(const vector<string> & vol_names, const vector<blastdb::TOid> & vol_num_oids)
    // {
    // 	x_IncreaseEnvMapSize(vol_names, vol_num_oids);
    //     lmdb::txn txn = lmdb::txn::begin(m_Env);
    //     lmdb::dbi volinfo = lmdb::dbi::open(txn, blastdb::volinfo_str.c_str(), MDB_CREATE | MDB_INTEGERKEY);
    //     lmdb::dbi volname = lmdb::dbi::open(txn, blastdb::volname_str.c_str(), MDB_CREATE | MDB_INTEGERKEY);
    // 	for (unsigned int i =0; i < vol_names.size(); i++) {
    // 		const lmdb::val key{&i, sizeof(i)};
    // 		lmdb::val value{vol_names[i].c_str(), strlen(vol_names[i].c_str())};
    // 		bool rc =lmdb::dbi_put(txn, volname.handle(), key, value);
    // 		if (!rc) {
    // 			NCBI_THROW( CSeqDBException, eArgErr, "VolNames error ");
    // 		}
    // 		rc = volinfo.put(txn, i, vol_num_oids[i]);
    // 		if (!rc) {
    // 			NCBI_THROW( CSeqDBException, eArgErr, "VolInfo error ");
    // 		}
    // 	}
    // 	txn.commit();
    // }
    // ```
    /// `InsertVolumesInfo`: the names of the volumes (`CDirEntry(volume name).GetName()`) and
    /// the number of OIDs in each.
    pub(crate) fn insert_volumes_info(
        &mut self,
        vol_names: &[String],
        vol_num_oids: &[u32],
    ) -> Result<(), WriteDbError> {
        self.increase_env_map_size_volumes(vol_names)?;
        self.env.txn_begin().map_err(lmdb_error)?;
        let r = self.insert_volumes_info0(vol_names, vol_num_oids);
        if r.is_err() {
            self.env.txn_abort();
        }
        r
    }

    fn insert_volumes_info0(
        &mut self,
        vol_names: &[String],
        vol_num_oids: &[u32],
    ) -> Result<(), WriteDbError> {
        let volinfo = self
            .env
            .dbi_open(VOLINFO_STR, MDB_CREATE | MDB_INTEGERKEY)
            .map_err(lmdb_error)?;
        let volname = self
            .env
            .dbi_open(VOLNAME_STR, MDB_CREATE | MDB_INTEGERKEY)
            .map_err(lmdb_error)?;
        for (i, name) in vol_names.iter().enumerate() {
            let key = (i as u32).to_le_bytes();
            if !self
                .env
                .put(volname, &key, c_str(name.as_bytes()), 0)
                .map_err(lmdb_error)?
            {
                return Err(WriteDbError::Other("VolNames error ".to_string()));
            }
            if !self
                .env
                .put(volinfo, &key, &vol_num_oids[i].to_le_bytes(), 0)
                .map_err(lmdb_error)?
            {
                return Err(WriteDbError::Other("VolInfo error ".to_string()));
            }
        }
        self.env.txn_commit().map_err(lmdb_error)
    }

    /// `set_mapsize` when `total_pages_needed + last_page_num > max_num_pages`.
    fn grow_map(&mut self, total_pages_needed: usize) -> Result<(), WriteDbError> {
        let page_size = self.env.psize();
        let last_page_num = self.env.last_pgno() as usize;
        let max_num_pages = self.env.mapsize() as usize / page_size;
        if total_pages_needed + last_page_num > max_num_pages {
            let new_map_size = (total_pages_needed + last_page_num) * page_size;
            self.env
                .set_mapsize(new_map_size as u64)
                .map_err(lmdb_error)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:198-225
    // ```c++
    // 	// 2 meta pages
    // 	const size_t MIN_PAGES = 3;
    // 	const size_t BRANCH_PAGES = 2;
    // 	// Each entry has 8 byte overhead + size of (key + entry)
    // 	size_t vol_name_size = (vol_names.front().size() + 24)* vol_names.size();
    // 	size_t vol_info_size = 24* vol_names.size();
    // 	...
    // 	size_t page_size = stat.ms_psize;
    // 	// For each page 16 byte header
    // 	size_t page_max_size = page_size -16;
    // 	size_t last_page_num = info.me_last_pgno;
    // 	size_t max_num_pages = info.me_mapsize/page_size;
    // 	size_t leaf_pages_needed = vol_name_size/page_max_size + vol_info_size/page_max_size + 2;
    // 	size_t total_pages_needed = MIN_PAGES + BRANCH_PAGES + leaf_pages_needed;
    // 	if( (total_pages_needed + last_page_num) > max_num_pages ) {
    // 		size_t newMapSize = (total_pages_needed + last_page_num) * page_size;
    // 		m_Env.set_mapsize(newMapSize);
    // ```
    fn increase_env_map_size_volumes(&mut self, vol_names: &[String]) -> Result<(), WriteDbError> {
        let front = vol_names.first().map_or(0, |n| n.len());
        let vol_name_size = (front + 24) * vol_names.len();
        let vol_info_size = 24 * vol_names.len();
        let page_max_size = self.env.psize() - 16;
        let leaf_pages_needed = vol_name_size / page_max_size + vol_info_size / page_max_size + 2;
        self.grow_map(3 + 2 + leaf_pages_needed)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:227-250
    // ```c++
    // 	size_t size = m_TotalIdsLength  + m_list.size() * 16;
    // 	size_t avg_id_length = m_TotalIdsLength/m_list.size();
    // 	...
    // 	size_t leaf_pages_needed = size/page_max_size + 1;
    // 	size_t dup_pages = (leaf_pages_needed > 200) ? 14: 7;
    // 	size_t branch_pages_needed = (avg_id_length + 16)* leaf_pages_needed/page_max_size + 1;
    // 	size_t total_pages_needed = leaf_pages_needed + branch_pages_needed + dup_pages;
    // ```
    fn increase_env_map_size(&mut self) -> Result<(), WriteDbError> {
        let size = self.total_ids_length + self.list.len() * 16;
        let avg_id_length = self.total_ids_length / self.list.len();
        let page_max_size = self.env.psize() - 16;
        let leaf_pages_needed = size / page_max_size + 1;
        let dup_pages = if leaf_pages_needed > 200 { 14 } else { 7 };
        let branch_pages_needed = (avg_id_length + 16) * leaf_pages_needed / page_max_size + 1;
        self.grow_map(leaf_pages_needed + branch_pages_needed + dup_pages)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:273-342
    // (x_CommitTransaction)
    // ```c++
    // 	if(m_list.size() == 0) {
    // 		return;
    // 	}
    // #ifdef _OPENMP
    // 	unsigned int min_split_size = DEFAULT_MIN_SPLIT_SORT_SIZE;
    // 	unsigned int chunk_size = DEFAULT_MIN_SPLIT_CHUNK_SIZE;
    // 	char* min_split_str = getenv("LMDB_MIN_SPLIT_SIZE");
    // 	char* chunk_str = getenv("LMDB_SPLIT_CHUNK_SIZE");
    // 	if (chunk_str) {
    // 		chunk_size = NStr::StringToUInt(chunk_str);
    // 	}
    // 	if (min_split_str) {
    // 		min_split_size = NStr::StringToUInt(min_split_str);
    // 	}
    // 	if((m_list.size() < min_split_size) || (m_list.size() < 2*chunk_size)) {
    // 		std::sort (m_list.begin(), m_list.end(), SKeyValuePair::cmp_key);
    // 	}
    // 	else {
    // 		... (x_Split: nth_element and OpenMP tasks that sort chunks)
    // 	}
    // #else
    // 	std::sort (m_list.begin(), m_list.end(), SKeyValuePair::cmp_key);
    // #endif
    // 	x_IncreaseEnvMapSize();
    //     unsigned int j=0;
    //     while (j < m_list.size()){
    //     	lmdb::txn txn = lmdb::txn::begin(m_Env);
    //     	lmdb::dbi dbi = lmdb::dbi::open(txn, blastdb::acc2oid_str.c_str(),
    //     			                        MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED);
    //     	unsigned int i = j;
    //     	j= i+m_MaxEntryPerTxn;
    //     	if(j > m_list.size()) {
    //     		j = m_list.size();
    //     	}
    //     	for(; i < j; i++){
    //     		if( i > 0) {
    //     			if ((m_list[i-1].id == m_list[i].id) &&
    //     				(m_list[i-1].oid == m_list[i].oid)){
    //     				continue;
    //     			}
    //     		}
    //     		blastdb::TOid & oid = m_list[i].oid;
    //     		string & id = m_list[i].id;
    // 			lmdb::val value{&oid, sizeof(oid)};
    // 			lmdb::val key{id.c_str(), strlen(id.c_str())};
    // 			bool rc = lmdb::dbi_put(txn, dbi.handle(), key, value, MDB_APPENDDUP);
    // 			if (!rc) {
    // 		 		NCBI_THROW( CSeqDBException, eArgErr, "acc2oid error for id " + id);
    // 			}
    // 		}
    //     	txn.commit();
    //     }
    // ```
    // makeblastdb 2.17.0 (bioconda) is built with OpenMP: it reads both variables. The split
    // sort gives the order `std::sort` gives (`cmp_key` orders by id, then by the bytes of the
    // OID; ties are equal in both), but with a chunk size of 0 or 1 `x_Split` recurses without
    // end (a range of one element is split into itself): LOSAT refuses those settings where
    // NCBI takes that path.
    fn commit_transaction(&mut self) -> Result<(), WriteDbError> {
        if self.list.is_empty() {
            return Ok(());
        }
        let mut min_split_size = DEFAULT_MIN_SPLIT_SORT_SIZE;
        let mut chunk_size = DEFAULT_MIN_SPLIT_CHUNK_SIZE;
        if let Some(v) = &self.settings.lmdb_split_chunk_size {
            chunk_size = string_to_uint(v).map_err(WriteDbError::Other)?;
        }
        if let Some(v) = &self.settings.lmdb_min_split_size {
            min_split_size = string_to_uint(v).map_err(WriteDbError::Other)?;
        }
        let n = self.list.len();
        let split = !(n < min_split_size as usize || n < chunk_size.wrapping_mul(2) as usize);
        if split && chunk_size <= 1 {
            return Err(WriteDbError::Other(format!(
                "LMDB_SPLIT_CHUNK_SIZE={chunk_size}: makeblastdb's split sort does not end"
            )));
        }
        // SKeyValuePair::cmp_key (writedb_lmdb.hpp:116-127): `v.id < k.id` (std::string,
        // bytewise), ties by the OID's bytes from the least significant one.
        let ids = std::mem::take(&mut self.ids);
        self.list.sort_unstable_by(|a, b| {
            let ia = &ids[a.off as usize..a.off as usize + a.len as usize];
            let ib = &ids[b.off as usize..b.off as usize + b.len as usize];
            ia.cmp(ib)
                .then_with(|| a.oid.to_le_bytes().cmp(&b.oid.to_le_bytes()))
        });
        self.ids = ids;
        self.increase_env_map_size()?;
        for (i0, j) in txn_batches(self.list.len(), self.max_entry_per_txn)? {
            self.env.txn_begin().map_err(lmdb_error)?;
            let r = self.commit_batch(i0, j);
            if r.is_err() {
                self.env.txn_abort();
                return r;
            }
        }
        Ok(())
    }

    fn commit_batch(&mut self, i0: usize, j: usize) -> Result<(), WriteDbError> {
        let dbi: Dbi = self
            .env
            .dbi_open(ACC2OID_STR, MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED)
            .map_err(lmdb_error)?;
        for i in i0..j {
            let kv = self.list[i];
            if i > 0 {
                let prev = self.list[i - 1];
                if self.id(&prev) == self.id(&kv) && prev.oid == kv.oid {
                    continue;
                }
            }
            let oid = kv.oid.to_le_bytes();
            let (s, e) = (kv.off as usize, kv.off as usize + kv.len as usize);
            let ok = self
                .env
                .put(dbi, c_str(&self.ids[s..e]), &oid, MDB_APPENDDUP)
                .map_err(lmdb_error)?;
            if !ok {
                return Err(WriteDbError::Other(format!(
                    "acc2oid error for id {}",
                    String::from_utf8_lossy(self.id(&kv))
                )));
            }
        }
        self.env.txn_commit().map_err(lmdb_error)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:344-415
    // ```c++
    // Uint4 s_WirteIds(CNcbiOfstream & os, vector<string> & ids)
    // {
    // 	const unsigned char byte_max = 0xFF;
    // 	Uint4 diff =0;
    // 	sort(ids.begin(), ids.end());
    // 	for(unsigned int j =0; j < ids.size(); j++) {
    // 		Uint4 id_len = ids[j].size();
    // 		if(id_len >= byte_max) {
    // 			os.write((char *)&byte_max, 1);
    // 			os.write((char *)&id_len, 4);
    // 			diff += 5;
    // 		}
    // 		else {
    // 			char l = byte_max & id_len;
    // 			os.write(&l,1);
    // 			diff += 1;
    // 		}
    // 		os.write(ids[j].c_str(), id_len);
    // 		diff += id_len;
    // 	}
    // 	return diff;
    // }
    // void CWriteDB_LMDB::x_CreateOidToSeqidsLookupFile()
    // {
    // 	if(m_list.size() == 0) {
    // 		return;
    // 	}
    // 	Uint8 total_num_oids = m_list.back().oid + 1;
    // 	string filename = GetFileNameFromExistingLMDBFile(m_Db, ELMDBFileType::eOid2SeqIds);
    // 	Uint8 offset = 0;
    // 	CNcbiOfstream os(filename.c_str(), IOS_BASE::out | IOS_BASE::binary);
    // 	vector<Uint4> offsets(total_num_oids, 0);
    // 	os.write((char *)&total_num_oids, 8);
    // 	for(unsigned int i=0; i < total_num_oids; i++) {
    // 		os.write((char *) &offset, 8);
    // 	}
    // 	os.flush();
    // 	blastdb::TOid count = 0;
    // 	vector<string> tmp_ids;
    // 	for(unsigned int i = 0; i < m_list.size(); i++) {
    // 		if(i > 0 && m_list[i].oid != m_list[i-1].oid ) {
    // 			if((m_list[i].oid - m_list[i-1].oid) != 1) {
    // 		 		NCBI_THROW( CSeqDBException, eArgErr, "Input id list not in ascending oid order");
    // 			}
    // 			offsets[count] = s_WirteIds(os, tmp_ids);
    // 			count++;
    // 			tmp_ids.clear();
    // 		}
    // 		m_TotalIdsLength +=m_list[i].id.size();
    // 		if(!m_list[i].saveToOidList) {
    // 			continue;
    // 		}
    // 		tmp_ids.push_back(m_list[i].id);
    // 	}
    // 	offsets[count] = s_WirteIds(os, tmp_ids);
    // 	_ASSERT(count == m_list.back().oid);
    // 	os.flush();
    // 	os.seekp(8);
    // 	for(unsigned int i = 0; i < total_num_oids; i++) {
    // 		offset += offsets[i];
    // 		os.write((char *) &offset, 8);
    // 	}
    // ```
    fn create_oid_to_seqids_lookup_file(&mut self) -> Result<(), WriteDbError> {
        if self.list.is_empty() {
            return Ok(());
        }
        let total_num_oids = self.list.last().unwrap().oid as u64 + 1;
        let filename = lmdb_file_name(&self.db, "os");
        let io = |source| WriteDbError::Io {
            path: filename.clone(),
            source,
        };
        let mut os = BufWriter::new(File::create(&filename).map_err(io)?);
        let mut offsets: Vec<u32> = vec![0; total_num_oids as usize];
        os.write_all(&total_num_oids.to_le_bytes()).map_err(io)?;
        for _ in 0..total_num_oids {
            os.write_all(&0u64.to_le_bytes()).map_err(io)?;
        }
        let mut count = 0usize;
        let mut tmp_ids: Vec<(u64, u32)> = Vec::new();
        for i in 0..self.list.len() {
            let kv = self.list[i];
            if i > 0 && kv.oid != self.list[i - 1].oid {
                if kv.oid.wrapping_sub(self.list[i - 1].oid) != 1 {
                    os.flush().map_err(io)?;
                    return Err(WriteDbError::Other(
                        "Input id list not in ascending oid order".to_string(),
                    ));
                }
                offsets[count] = self.write_ids(&mut os, &mut tmp_ids).map_err(io)?;
                count += 1;
                tmp_ids.clear();
            }
            self.total_ids_length += kv.len as usize;
            if !kv.save_to_oid_list {
                continue;
            }
            tmp_ids.push((kv.off, kv.len));
        }
        offsets[count] = self.write_ids(&mut os, &mut tmp_ids).map_err(io)?;
        os.seek(SeekFrom::Start(8)).map_err(io)?;
        let mut offset: u64 = 0;
        for o in &offsets {
            offset += *o as u64;
            os.write_all(&offset.to_le_bytes()).map_err(io)?;
        }
        os.flush().map_err(io)
    }

    fn write_ids(&self, os: &mut impl Write, ids: &mut [(u64, u32)]) -> std::io::Result<u32> {
        let get = |&(off, len): &(u64, u32)| &self.ids[off as usize..off as usize + len as usize];
        ids.sort_unstable_by(|a, b| get(a).cmp(get(b)));
        let mut diff: u32 = 0;
        for e in ids.iter() {
            let id = get(e);
            let id_len = id.len() as u32;
            if id_len >= 0xFF {
                os.write_all(&[0xFF])?;
                os.write_all(&id_len.to_le_bytes())?;
                diff += 5;
            } else {
                os.write_all(&[id_len as u8])?;
                diff += 1;
            }
            os.write_all(id)?;
            diff += id_len;
        }
        Ok(diff)
    }

    /// `~CWriteDB_LMDB`: the `.pos` file, the `acc2oid` transactions, the environment closed
    /// and the lock file removed.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72
    // ```c++
    // CWriteDB_LMDB::~CWriteDB_LMDB()
    // {
    // 	x_CreateOidToSeqidsLookupFile();
    // 	x_CommitTransaction();
    //     CBlastLMDBManager::GetInstance().CloseEnv(m_Db);
    //     CFile(m_Db+"-lock").Remove();
    // }
    // ```
    pub(crate) fn close(mut self) -> Result<(), WriteDbError> {
        self.create_oid_to_seqids_lookup_file()?;
        self.commit_transaction()?;
        self.env.close().map_err(lmdb_error)
    }
}

/// `CWriteDB_TaxID`: the `.ptf`/`.ntf` file and the `.pot`/`.not` and `.pto`/`.nto` files.
pub(crate) struct WriteDbTaxId {
    db: PathBuf,
    env: Env,
    max_entry_per_txn: u32,
    /// `m_TaxId2OidList`: (tax id, OID).
    taxid2oid: Vec<(i32, u32)>,
    /// `m_TaxId2OffsetsList`: (tax id, offset in the `.pto` file).
    taxid2offsets: Vec<(i32, u64)>,
}

impl WriteDbTaxId {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:425-435
    // ```c++
    // CWriteDB_TaxID::CWriteDB_TaxID(const string& dbname,  Uint8 map_size, Uint8 capacity): m_Db(dbname),
    //                                m_Env(CBlastLMDBManager::GetInstance().GetWriteEnv(dbname, map_size)),
    //                                m_ListCapacity(capacity), m_MaxEntryPerTxn(DEFAULT_MAX_ENTRY_PER_TXN)
    // {
    // 	m_TaxId2OidList.reserve(m_ListCapacity);
    // 	char* max_entry_str = getenv("MAX_LMDB_TXN_ENTRY");
    // 	if (max_entry_str) {
    // 		m_MaxEntryPerTxn = NStr::StringToInt(max_entry_str);
    // ```
    /// Create the LMDB file `dbname` (`<db>.ptf`/`<db>.ntf`).
    pub(crate) fn new(
        dbname: &Path,
        map_size: u64,
        page_size: usize,
        settings: &LmdbEnvSettings,
    ) -> Result<Self, WriteDbError> {
        let env = Env::create(dbname, map_size, page_size).map_err(lmdb_error)?;
        let max_entry_per_txn = max_entry_per_txn(settings)?;
        Ok(WriteDbTaxId {
            db: dbname.to_path_buf(),
            env,
            max_entry_per_txn,
            taxid2oid: Vec::new(),
            taxid2offsets: Vec::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:446-464
    // ```c++
    // int CWriteDB_TaxID::InsertEntries(const set<TTaxId> & tax_ids, const blastdb::TOid oid)
    // {
    //     int count = 0;
    //     if(tax_ids.size() == 0) {
    //     	x_Resize();
    //     	SKeyValuePair<blastdb::TOid>  kv(ZERO_TAX_ID, oid);
    //     	m_TaxId2OidList.push_back(kv);
    //     	return 1;
    //     }
    //     ITERATE(set<TTaxId>, itr, tax_ids) {
    //     	x_Resize();
    //     	SKeyValuePair<blastdb::TOid> kv(*itr, oid);
    //     	m_TaxId2OidList.push_back(kv);
    //     	count++;
    //     }
    //     return count;
    // }
    // ```
    /// `InsertEntries`: the tax ids of OID `oid` in the order of a `set<TTaxId>` (ascending,
    /// each once). Called in OID order.
    pub(crate) fn insert_entries(&mut self, tax_ids: &[i32], oid: u32) {
        debug_assert!(tax_ids.windows(2).all(|w| w[0] < w[1]));
        if tax_ids.is_empty() {
            self.taxid2oid.push((0, oid));
            return;
        }
        for &t in tax_ids {
            self.taxid2oid.push((t, oid));
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:523-576
    // ```c++
    // Uint4 s_WirteTaxIds(CNcbiOfstream & os, vector<TTaxId> & tax_ids)
    // {
    // 	for(unsigned int j =0; j < tax_ids.size(); j++) {
    //         Int4 tid = TAX_ID_TO(Int4, tax_ids[j]);
    // 		os.write((char *)&tid, 4);
    // 	}
    // 	return tax_ids.size();
    // }
    // void CWriteDB_TaxID::x_CreateOidToTaxIdsLookupFile()
    // {
    // 	if(m_TaxId2OidList.size() == 0) {
    //  		NCBI_THROW( CSeqDBException, eArgErr, "No tax info for any oid");
    // 	}
    // 	Uint8 total_num_oids = m_TaxId2OidList.back().value + 1;
    // 	string filename = GetFileNameFromExistingLMDBFile(m_Db, ELMDBFileType::eOid2TaxIds);
    // 	Uint8 offset = 0;
    // 	CNcbiOfstream os(filename.c_str(), IOS_BASE::out | IOS_BASE::binary);
    // 	vector<Uint4> offsets(total_num_oids, 0);
    // 	os.write((char *)&total_num_oids, 8);
    // 	for(unsigned int i=0; i < total_num_oids; i++) {
    // 		os.write((char *) &offset, 8);
    // 	}
    // 	os.flush();
    // 	blastdb::TOid count = 0;
    // 	vector<TTaxId> tmp_tax_ids;
    // 	for(unsigned int i = 0; i < m_TaxId2OidList.size(); i++) {
    // 		if(i > 0 && m_TaxId2OidList[i].value != m_TaxId2OidList[i-1].value ) {
    // 			if((m_TaxId2OidList[i].value - m_TaxId2OidList[i-1].value) != 1) {
    // 		 		NCBI_THROW( CSeqDBException, eArgErr, "Input id list not in ascending oid order");
    // 			}
    // 			offsets[count] = s_WirteTaxIds(os, tmp_tax_ids);
    // 			count++;
    // 			tmp_tax_ids.clear();
    // 		}
    // 		tmp_tax_ids.push_back(m_TaxId2OidList[i].tax_id);
    // 	}
    // 	offsets[count] = s_WirteTaxIds(os, tmp_tax_ids);
    // 	_ASSERT(count == m_TaxId2OidList.back().value);
    // 	os.flush();
    // 	os.seekp(8);
    // 	for(unsigned int i = 0; i < total_num_oids; i++) {
    // 		offset += offsets[i];
    // 		os.write((char *) &offset, 8);
    // 	}
    // ```
    fn create_oid_to_taxids_lookup_file(&mut self) -> Result<(), WriteDbError> {
        if self.taxid2oid.is_empty() {
            return Err(WriteDbError::Other("No tax info for any oid".to_string()));
        }
        let total_num_oids = self.taxid2oid.last().unwrap().1 as u64 + 1;
        let filename = lmdb_file_name(&self.db, "ot");
        let io = |source| WriteDbError::Io {
            path: filename.clone(),
            source,
        };
        let mut os = BufWriter::new(File::create(&filename).map_err(io)?);
        let mut offsets: Vec<u32> = vec![0; total_num_oids as usize];
        os.write_all(&total_num_oids.to_le_bytes()).map_err(io)?;
        for _ in 0..total_num_oids {
            os.write_all(&0u64.to_le_bytes()).map_err(io)?;
        }
        let mut count = 0usize;
        let mut start = 0usize;
        for i in 0..self.taxid2oid.len() {
            let (_, oid) = self.taxid2oid[i];
            if i > 0 && oid != self.taxid2oid[i - 1].1 {
                if oid.wrapping_sub(self.taxid2oid[i - 1].1) != 1 {
                    os.flush().map_err(io)?;
                    return Err(WriteDbError::Other(
                        "Input id list not in ascending oid order".to_string(),
                    ));
                }
                for &(t, _) in &self.taxid2oid[start..i] {
                    os.write_all(&t.to_le_bytes()).map_err(io)?;
                }
                offsets[count] = (i - start) as u32;
                count += 1;
                start = i;
            }
        }
        for &(t, _) in &self.taxid2oid[start..] {
            os.write_all(&t.to_le_bytes()).map_err(io)?;
        }
        offsets[count] = (self.taxid2oid.len() - start) as u32;
        os.seek(SeekFrom::Start(8)).map_err(io)?;
        let mut offset: u64 = 0;
        for o in &offsets {
            offset += *o as u64;
            os.write_all(&offset.to_le_bytes()).map_err(io)?;
        }
        os.flush().map_err(io)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:578-613
    // ```c++
    // Uint4 s_WirteOids(CNcbiOfstream & os, vector<blastdb::TOid> & oids)
    // {
    // 	blastdb::SortAndUnique <blastdb::TOid> (oids);
    // 	Uint4 num_oids = oids.size();
    // 	os.write((char *)&num_oids, 4);
    // 	for(unsigned int j =0; j < num_oids; j++) {
    // 		os.write((char *)&oids[j], 4);
    // 	}
    // 	return ((num_oids +1) * 4);
    // }
    // void CWriteDB_TaxID::x_CreateTaxIdToOidsLookupFile()
    // {
    //     sort (m_TaxId2OidList.begin(), m_TaxId2OidList.end(), SKeyValuePair<blastdb::TOid>::cmp_key);
    // 	string filename = GetFileNameFromExistingLMDBFile(m_Db, ELMDBFileType::eTaxId2Oids);
    // 	CNcbiOfstream os(filename.c_str(), IOS_BASE::out | IOS_BASE::binary);
    // 	Uint8 offset =0;
    // 	vector<blastdb::TOid> tmp_oids;
    // 	for(unsigned int i = 0; i < m_TaxId2OidList.size(); i++) {
    // 		if(i > 0 && m_TaxId2OidList[i].tax_id != m_TaxId2OidList[i-1].tax_id ) {
    // 			SKeyValuePair<Uint8> p(m_TaxId2OidList[i-1].tax_id, offset);
    // 			offset += s_WirteOids(os, tmp_oids);
    // 			m_TaxId2OffsetsList.push_back(p);
    // 			tmp_oids.clear();
    // 		}
    // 		tmp_oids.push_back(m_TaxId2OidList[i].value);
    // 	}
    //    	SKeyValuePair<Uint8> kv(m_TaxId2OidList.back().tax_id, offset);
    //    	s_WirteOids(os, tmp_oids);
    // 	m_TaxId2OffsetsList.push_back(kv);
    // ```
    // `cmp_key` orders by tax id (signed), then by OID; the pairs that tie are equal.
    fn create_taxid_to_oids_lookup_file(&mut self) -> Result<(), WriteDbError> {
        self.taxid2oid.sort_unstable();
        let filename = lmdb_file_name(&self.db, "to");
        let io = |source| WriteDbError::Io {
            path: filename.clone(),
            source,
        };
        let mut os = BufWriter::new(File::create(&filename).map_err(io)?);
        let mut offset: u64 = 0;
        let mut tmp_oids: Vec<u32> = Vec::new();
        let write_oids = |os: &mut BufWriter<File>, oids: &mut Vec<u32>| -> std::io::Result<u64> {
            oids.sort_unstable();
            oids.dedup();
            os.write_all(&(oids.len() as u32).to_le_bytes())?;
            for o in oids.iter() {
                os.write_all(&o.to_le_bytes())?;
            }
            Ok((oids.len() as u64 + 1) * 4)
        };
        for i in 0..self.taxid2oid.len() {
            if i > 0 && self.taxid2oid[i].0 != self.taxid2oid[i - 1].0 {
                let p = (self.taxid2oid[i - 1].0, offset);
                offset += write_oids(&mut os, &mut tmp_oids).map_err(io)?;
                self.taxid2offsets.push(p);
                tmp_oids.clear();
            }
            tmp_oids.push(self.taxid2oid[i].1);
        }
        let kv = (self.taxid2oid.last().unwrap().0, offset);
        write_oids(&mut os, &mut tmp_oids).map_err(io)?;
        self.taxid2offsets.push(kv);
        os.flush().map_err(io)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:466-521
    // ```c++
    // void CWriteDB_TaxID::x_IncreaseEnvMapSize()
    // {
    // 	const size_t MIN_PAGES = 4;
    // 	...
    // 	size_t size = m_TaxId2OffsetsList.size()*32;
    // 	size_t page_size = stat.ms_psize;
    // 	size_t page_max_size = stat.ms_psize - 16;
    // 	size_t last_page_num = info.me_last_pgno;
    // 	size_t max_num_pages = info.me_mapsize/page_size;
    // 	size_t leaf_pages_needed = size/page_max_size + 1;
    // 	size_t branch_pages_needed = 24 * leaf_pages_needed/page_max_size + 1;
    // 	size_t total_pages_needed = leaf_pages_needed + branch_pages_needed + MIN_PAGES;
    // 	if( (total_pages_needed + last_page_num) > max_num_pages) {
    // 		size_t newMapSize = (total_pages_needed + last_page_num) * page_size;
    // 		m_Env.set_mapsize(newMapSize);
    // ...
    // void CWriteDB_TaxID::x_CommitTransaction()
    // {
    // 	_ASSERT(m_TaxId2OffsetsList.size());
    //     sort (m_TaxId2OffsetsList.begin(), m_TaxId2OffsetsList.end(), SKeyValuePair<Uint8>::cmp_key);
    //     x_IncreaseEnvMapSize();
    //     unsigned int j=0;
    //     while (j < m_TaxId2OffsetsList.size()){
    //     	lmdb::txn txn = lmdb::txn::begin(m_Env);
    //     	lmdb::dbi dbi = lmdb::dbi::open(txn, blastdb::taxid2offset_str.c_str(),
    //     			                        MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED);
    //     	unsigned int i = j;
    //     	j= i+m_MaxEntryPerTxn;
    //     	if(j > m_TaxId2OffsetsList.size()) {
    //     		j = m_TaxId2OffsetsList.size();
    //     	}
    //     	for(; i < j; i++){
    //     		Uint8 & offset = m_TaxId2OffsetsList[i].value;
    //             TTaxId & tax_id = m_TaxId2OffsetsList[i].tax_id;
    // 			lmdb::val value{&offset, sizeof(offset)};
    // 			lmdb::val key{&tax_id, sizeof(tax_id)};
    // 			bool rc = lmdb::dbi_put(txn, dbi.handle(), key, value, MDB_APPENDDUP);
    // 			if (!rc) {
    // 		 		NCBI_THROW( CSeqDBException, eArgErr, "taxid2offset error for tax id " + NStr::NumericToString(tax_id));
    // 			}
    // 		}
    //     	txn.commit();
    //     }
    // ```
    fn commit_transaction(&mut self) -> Result<(), WriteDbError> {
        self.taxid2offsets.sort_unstable();
        let size = self.taxid2offsets.len() * 32;
        let page_size = self.env.psize();
        let page_max_size = page_size - 16;
        let last_page_num = self.env.last_pgno() as usize;
        let max_num_pages = self.env.mapsize() as usize / page_size;
        let leaf_pages_needed = size / page_max_size + 1;
        let branch_pages_needed = 24 * leaf_pages_needed / page_max_size + 1;
        let total_pages_needed = leaf_pages_needed + branch_pages_needed + 4;
        if total_pages_needed + last_page_num > max_num_pages {
            let new_map_size = (total_pages_needed + last_page_num) * page_size;
            self.env
                .set_mapsize(new_map_size as u64)
                .map_err(lmdb_error)?;
        }
        for (i0, j) in txn_batches(self.taxid2offsets.len(), self.max_entry_per_txn)? {
            self.env.txn_begin().map_err(lmdb_error)?;
            let r = self.commit_batch(i0, j);
            if r.is_err() {
                self.env.txn_abort();
                return r;
            }
        }
        Ok(())
    }

    fn commit_batch(&mut self, i0: usize, j: usize) -> Result<(), WriteDbError> {
        let dbi = self
            .env
            .dbi_open(TAXID2OFFSET_STR, MDB_CREATE | MDB_DUPSORT | MDB_DUPFIXED)
            .map_err(lmdb_error)?;
        for i in i0..j {
            let (tax_id, offset) = self.taxid2offsets[i];
            let ok = self
                .env
                .put(
                    dbi,
                    &tax_id.to_le_bytes(),
                    &offset.to_le_bytes(),
                    MDB_APPENDDUP,
                )
                .map_err(lmdb_error)?;
            if !ok {
                return Err(WriteDbError::Other(format!(
                    "taxid2offset error for tax id {tax_id}"
                )));
            }
        }
        self.env.txn_commit().map_err(lmdb_error)
    }

    /// `~CWriteDB_TaxID`: the `.pot` and `.pto` files, the `taxid2offset` transactions, the
    /// environment closed and the lock file removed.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:437-444
    // ```c++
    // CWriteDB_TaxID::~CWriteDB_TaxID()
    // {
    // 	x_CreateOidToTaxIdsLookupFile();
    // 	x_CreateTaxIdToOidsLookupFile();
    // 	x_CommitTransaction();
    //     CBlastLMDBManager::GetInstance().CloseEnv(m_Db);
    //     CFile(m_Db+"-lock").Remove();
    // }
    // ```
    pub(crate) fn close(mut self) -> Result<(), WriteDbError> {
        self.create_oid_to_taxids_lookup_file()?;
        self.create_taxid_to_oids_lookup_file()?;
        self.commit_transaction()?;
        self.env.close().map_err(lmdb_error)
    }
}

#[cfg(test)]
#[path = "lmdb_tests.rs"]
mod tests;
