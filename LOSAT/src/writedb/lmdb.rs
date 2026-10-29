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
use super::platform::LmdbPlatform;
use super::WriteDbError;
use crate::seqid::{label, ParsedSeqId};

// The default map sizes (`DEFAULT_LMDB_MAP_SIZE`, `DEFAULT_TAXID_MAP_SIZE`) depend on the
// platform: `super::platform::Platform`.

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

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_lmdb.hpp:111-129
// ```c++
//     struct SKeyValuePair {
//     	string id;
//     	blastdb::TOid oid;
//     	bool saveToOidList;
//     	SKeyValuePair() : id(kEmptyStr), oid(kSeqDBEntryNotFound), saveToOidList(false) {}
//     	static bool cmp_key(const SKeyValuePair & v, const SKeyValuePair & k) {
//    			if(v.id == k.id) {
// 				blastdb::TOid mask = 0xFF;
//    				for(unsigned int i=0; i < sizeof(blastdb::TOid); i++) {
//    					if( (v.oid & mask) != (k.oid & mask)) {
//    						return (v.oid & mask) < (k.oid & mask);
//    					}
//    					mask =  mask << 8;
//    				}
//    			}
//     	    return v.id < k.id;
//     	}
//     };
//     vector<SKeyValuePair> m_list;
// ```
// NCBI's element is a `std::string` (32 bytes, its characters on the heap past 15 of them)
// beside the OID and the flag: 40 bytes plus the heap block, about 72 bytes for an id of
// 16-24 characters, kept for every key of the database until `~CWriteDB_LMDB`. LOSAT keeps
// the same list in 12 bytes per element: the OID and one word that locates the id in one
// byte arena (`WriteDbLmdb::ids`, the bytes of the ids one after another). The keys that
// `x_InsertEntry` derives from one Seq-id share its bytes: the id without the version is the
// start of the id with it, a PDB molecule the start of the molecule with its chain, and an
// upper-case copy (`NStr::ToUpper`) is the same bytes with a flag. The order of the list and
// every comparison are those of the strings (`id_cmp`).
/// One element of `CWriteDB_LMDB::m_list` (`SKeyValuePair`): the OID and a word with the
/// place of the id in `WriteDbLmdb::ids` (bits 0-39 the offset, bits 40-61 the length), bit
/// 62 when the id is the upper-case copy of those bytes and bit 63 `saveToOidList`.
#[derive(Clone, Copy, Debug)]
#[repr(C, packed(4))]
struct KeyValuePair {
    word: u64,
    oid: u32,
}

const _: () = assert!(std::mem::size_of::<KeyValuePair>() == 12);

const KV_OFF_BITS: u32 = 40;
const KV_MAX_OFF: u64 = (1 << KV_OFF_BITS) - 1;
const KV_MAX_LEN: u64 = (1 << 22) - 1;
const KV_UPPER: u64 = 1 << 62;
const KV_SAVE: u64 = 1 << 63;

impl KeyValuePair {
    /// The element for `len` bytes of the arena at `off`; an id past what the word holds (an
    /// arena of 1 TiB, an id of 4 MiB) is refused.
    fn new(
        off: usize,
        len: usize,
        oid: u32,
        upper: bool,
        save: bool,
    ) -> Result<Self, WriteDbError> {
        let (off, len) = (off as u64, len as u64);
        if off > KV_MAX_OFF || len > KV_MAX_LEN {
            return Err(WriteDbError::Other(format!(
                "LOSAT's makeblastdb holds Seq-id keys of at most {KV_MAX_LEN} bytes and {} bytes \
                 of keys in all (a key of {len} bytes at {off})",
                KV_MAX_OFF + 1
            )));
        }
        let mut word = off | len << KV_OFF_BITS;
        if upper {
            word |= KV_UPPER;
        }
        if save {
            word |= KV_SAVE;
        }
        Ok(KeyValuePair { word, oid })
    }

    fn off(self) -> usize {
        (self.word & KV_MAX_OFF) as usize
    }

    fn len(self) -> usize {
        ((self.word >> KV_OFF_BITS) & KV_MAX_LEN) as usize
    }

    fn upper(self) -> bool {
        self.word & KV_UPPER != 0
    }

    fn save_to_oid_list(self) -> bool {
        self.word & KV_SAVE != 0
    }

    fn oid(self) -> u32 {
        self.oid
    }

    /// The bytes of the id before the upper-case copy is made.
    fn raw(self, ids: &[u8]) -> &[u8] {
        &ids[self.off()..self.off() + self.len()]
    }

    /// The id (`SKeyValuePair::id`) into `out`.
    fn id_into(self, ids: &[u8], out: &mut Vec<u8>) {
        out.clear();
        out.extend_from_slice(self.raw(ids));
        if self.upper() {
            out.make_ascii_uppercase();
        }
    }
}

/// `a.id` against `b.id` as `std::string` compares them (`char_traits<char>::compare`: the
/// bytes as `unsigned char`, then the shorter first).
fn id_cmp(ids: &[u8], a: KeyValuePair, b: KeyValuePair) -> std::cmp::Ordering {
    let (x, y) = (a.raw(ids), b.raw(ids));
    let (ua, ub) = (a.upper(), b.upper());
    if !ua && !ub {
        return x.cmp(y);
    }
    for (&p, &q) in x.iter().zip(y) {
        let p = if ua { p.to_ascii_uppercase() } else { p };
        let q = if ub { q.to_ascii_uppercase() } else { q };
        if p != q {
            return p.cmp(&q);
        }
    }
    x.len().cmp(&y.len())
}

/// `SKeyValuePair::cmp_key`: by id, then by the bytes of the OID from the least significant
/// one (an OID is below 2^31, so its top byte is not negative in NCBI's `int` mask).
fn kv_cmp(ids: &[u8], a: KeyValuePair, b: KeyValuePair) -> std::cmp::Ordering {
    id_cmp(ids, a, b).then_with(|| a.oid().to_le_bytes().cmp(&b.oid().to_le_bytes()))
}

// NCBI reference (598d8ae6): c++/include/util/lmdbxx/lmdb++.h:850-861
// ```c++
// static inline bool
// lmdb::dbi_put(MDB_txn* const txn,
//               const MDB_dbi dbi,
//               const MDB_val* const key,
//               const MDB_val* const data,
//               const unsigned int flags = 0) {
//   const int rc = ::mdb_put(txn, dbi, const_cast<MDB_val*>(key), const_cast<MDB_val*>(data), flags);
//   if (rc != MDB_SUCCESS && rc != MDB_KEYEXIST) {
//     error::raise("mdb_put", rc);
//   }
//   return (rc == MDB_SUCCESS);
// }
// ```
// The error of a put in `x_CommitTransaction` (an id of more than 511 bytes, which LMDB refuses
// with MDB_BAD_VALSIZE) is raised in `~CWriteDB_LMDB`, where it ends makeblastdb 2.17.0 with
// SIGABRT after its `.pos` file (S3c x_key512). LOSAT ends with lmdb++'s text and the key it
// could not put (makedb-v5 DECISIONS 10: no NCBI abort reproduced).
/// lmdb++'s error of a put into `acc2oid`, with the key.
fn acc2oid_put_error(func: &'static str, code: i32, key: &[u8]) -> WriteDbError {
    WriteDbError::Other(format!(
        "{} (the acc2oid key '{}', {} bytes; makeblastdb 2.17.0 aborts here)",
        LmdbError::Mdb { func, code },
        String::from_utf8_lossy(key),
        key.len()
    ))
}

/// `CWriteDB_LMDB`: the `.pdb`/`.ndb` file and the `.pos`/`.nos` file of a database.
pub(crate) struct WriteDbLmdb {
    db: PathBuf,
    env: Env,
    max_entry_per_txn: u32,
    total_ids_length: usize,
    settings: LmdbEnvSettings,
    /// The bytes of the ids of `list`.
    ids: Vec<u8>,
    /// `m_list`.
    list: Vec<KeyValuePair>,
    /// A key of the current Seq-id (`x_InsertEntry`) or of a put.
    scratch: Vec<u8>,
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
    // `m_list.reserve` (and `x_Resize`, which doubles the reserve) changes no output; LOSAT's
    // list grows as it is filled.
    /// Create the LMDB file `dbname` (`<db>.pdb`/`<db>.ndb`) with the map size `map_size`
    /// on the platform `platform` (its page size and file length; makeblastdb's bytes there).
    pub(crate) fn new(
        dbname: &Path,
        map_size: u64,
        platform: &LmdbPlatform,
        settings: &LmdbEnvSettings,
    ) -> Result<Self, WriteDbError> {
        let env = Env::create(dbname, map_size, platform).map_err(lmdb_error)?;
        let max_entry_per_txn = max_entry_per_txn(settings)?;
        Ok(WriteDbLmdb {
            db: dbname.to_path_buf(),
            env,
            max_entry_per_txn,
            total_ids_length: 0,
            settings: settings.clone(),
            ids: Vec::new(),
            list: Vec::new(),
            scratch: Vec::new(),
        })
    }

    /// Add one `SKeyValuePair` to `m_list`: the key string, the OID and whether the id goes
    /// into the `.pos` file (`saveToOidList`). Called in OID order. (`insert_entries` makes
    /// the elements of the Seq-ids of a sequence; this one takes a key as it is, for the
    /// tests that rebuild makeblastdb's lists from its files.)
    #[cfg(test)]
    pub(crate) fn insert_entry(
        &mut self,
        id: &[u8],
        oid: u32,
        save_to_oid_list: bool,
    ) -> Result<(), WriteDbError> {
        let off = self.ids.len();
        self.ids.extend_from_slice(id);
        self.push(off, id.len(), oid, false, save_to_oid_list)
    }

    fn push(
        &mut self,
        off: usize,
        len: usize,
        oid: u32,
        upper: bool,
        save: bool,
    ) -> Result<(), WriteDbError> {
        self.list
            .push(KeyValuePair::new(off, len, oid, upper, save)?);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:108-117
    // ```c++
    // int CWriteDB_LMDB::InsertEntries(const vector<CRef<CSeq_id>> & seqids, const blastdb::TOid oid)
    // {
    //     int count = 0;
    //     ITERATE(vector<CRef<CSeq_id>>, itr, seqids) {
    //        x_InsertEntry((*itr), oid);
    //     	count++;
    //     }
    //     return count;
    // }
    // ```
    /// `InsertEntries(m_Ids, oid)`: the keys of the Seq-ids of the sequence of OID `oid`.
    /// Returns how many elements were added to `m_list` (0 when every id is a GI).
    pub(crate) fn insert_entries(
        &mut self,
        seqids: &[ParsedSeqId],
        oid: u32,
    ) -> Result<usize, WriteDbError> {
        let before = self.list.len();
        for seqid in seqids {
            self.insert_seq_id(seqid, oid)?;
        }
        Ok(self.list.len() - before)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:119-196
    // ```c++
    // void CWriteDB_LMDB::x_InsertEntry(const CRef<CSeq_id> &seqid, const blastdb::TOid oid)
    // {
    //     if(seqid->IsGi()) {
    //     	return;
    //     }
    //
    //     x_Resize();
    //
    //     if(seqid->IsPir() || seqid->IsPrf()) {
    //     	SKeyValuePair kv_pir;
    //     	kv_pir.id = seqid->AsFastaString();
    //     	kv_pir.oid = oid;
    //     	kv_pir.saveToOidList = true;
    //     	m_list.push_back(kv_pir);
    //     	return;
    //     }
    //
    //     if(seqid->IsPdb()) {
    //     	SKeyValuePair kv_pdb_mol;
    //     	kv_pdb_mol.id = seqid->GetPdb().GetMol().Get();
    //     	kv_pdb_mol.oid = oid;
    //     	m_list.push_back(kv_pdb_mol);
    //     	// mol code should be case insensitive but c++ tooklit
    //     	// is not converting it all to uppercase now
    //     	string id_upper = kv_pdb_mol.id;
    //     	NStr::ToUpper(id_upper);
    //     	if(kv_pdb_mol.id != id_upper) {
    //     		SKeyValuePair kv_u;
    //        		kv_u.id = id_upper;
    //        		kv_u.oid = oid;
    //        		m_list.push_back(kv_u);
    //        	}
    //     	SKeyValuePair kv;
    //     	kv.id = seqid->GetSeqIdString(true);
    //     	kv.oid = oid;
    //     	kv.saveToOidList = true;
    //     	m_list.push_back(kv);
    //     	return;
    //     }
    //
    //     if(seqid->GetTextseq_Id() != NULL) {
    //     	SKeyValuePair kv;
    //     	kv.id = seqid->GetSeqIdString(false);
    //     	kv.oid = oid;
    //     	string id_v = seqid->GetSeqIdString(true);
    //     	if( kv.id != id_v) {
    //     		m_list.push_back(kv);
    //     		SKeyValuePair kv_v;
    //     		kv_v.id = id_v;
    //     		kv_v.oid = oid;
    //     		kv_v.saveToOidList = true;
    //     		m_list.push_back(kv_v);
    //     	}
    //     	else {
    //     		kv.saveToOidList = true;
    //     		m_list.push_back(kv);
    //     	}
    //    	}
    //     else {
    //     	SKeyValuePair kv;
    //     	kv.id = seqid->GetSeqIdString(true);
    //     	kv.oid = oid;
    //     	kv.saveToOidList = true;
    //     	m_list.push_back(kv);
    //     	string id_upper = kv.id;
    //     	NStr::ToUpper(id_upper);
    //     	if(kv.id != id_upper) {
    //     		SKeyValuePair kv_u;
    //     		kv_u.id = id_upper;
    //     		kv_u.oid = oid;
    //     		m_list.push_back(kv_u);
    //
    //     	}
    //
    //     }
    //     return;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1949-1954
    // ```c++
    //         if ( !str.empty() ) {
    //             if ( (flags & CSeq_id::fLabel_Version)  &&  tsid->IsSetVersion()) {
    //                 str += "." + NStr::IntToString(tsid->GetVersion());
    //             }
    //         }
    //         *label += str;
    // ```
    // `NStr::ToUpper` changes the letters a-z only (the C locale; the ids of the FASTA reader
    // hold no byte past 0x7f, makedb-v5 DECISIONS 28), so a copy differs exactly when the id
    // has such a letter. Each Seq-id's bytes are stored once (`KeyValuePair`): the id with the
    // version is written and the id without it is its start (`GetSeqIdString(false)` is the
    // start of `GetSeqIdString(true)`, above); a PDB molecule is the start of
    // `GetSeqIdString(true)` (the molecule, then `_` and the chain: `label::content`). Where
    // a start is not the shorter id, it is stored on its own.
    /// `x_InsertEntry`: the elements of `m_list` of one Seq-id of OID `oid`.
    fn insert_seq_id(&mut self, seqid: &ParsedSeqId, oid: u32) -> Result<(), WriteDbError> {
        if matches!(seqid, ParsedSeqId::Gi(_)) {
            return Ok(());
        }
        let off = self.ids.len();
        match seqid {
            ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => {
                label::fasta(seqid, &mut self.ids);
                self.push(off, self.ids.len() - off, oid, false, true)
            }
            ParsedSeqId::Pdb(pdb) => {
                label::content(seqid, true, &mut self.ids);
                let full = self.ids.len() - off;
                let mol = if self.ids[off..].starts_with(&pdb.mol) {
                    off
                } else {
                    let at = self.ids.len();
                    self.ids.extend_from_slice(&pdb.mol);
                    at
                };
                let len = pdb.mol.len();
                self.push(mol, len, oid, false, false)?;
                if pdb.mol.iter().any(u8::is_ascii_lowercase) {
                    self.push(mol, len, oid, true, false)?;
                }
                self.push(off, full, oid, false, true)
            }
            _ if seqid.textseq().is_some() => {
                label::content(seqid, true, &mut self.ids);
                let full = self.ids.len() - off;
                self.scratch.clear();
                label::content(seqid, false, &mut self.scratch);
                if self.scratch[..] != self.ids[off..] {
                    let without = if self.ids[off..].starts_with(&self.scratch) {
                        off
                    } else {
                        let at = self.ids.len();
                        self.ids.extend_from_slice(&self.scratch);
                        at
                    };
                    self.push(without, self.scratch.len(), oid, false, false)?;
                    self.push(off, full, oid, false, true)
                } else {
                    self.push(off, full, oid, false, true)
                }
            }
            _ => {
                label::content(seqid, true, &mut self.ids);
                let len = self.ids.len() - off;
                self.push(off, len, oid, false, true)?;
                if self.ids[off..].iter().any(u8::is_ascii_lowercase) {
                    self.push(off, len, oid, true, false)?;
                }
                Ok(())
            }
        }
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

    /// `set_mapsize` when `total_pages_needed + last_page_num > max_num_pages`. NCBI's
    /// `size_t` arithmetic is done in `u64` on every platform (a 32-bit `usize` would cut
    /// it); a size the platform cannot map is refused by `Env::set_mapsize`.
    fn grow_map(&mut self, total_pages_needed: u64) -> Result<(), WriteDbError> {
        let page_size = self.env.psize() as u64;
        let last_page_num = self.env.last_pgno();
        let max_num_pages = self.env.mapsize() / page_size;
        if total_pages_needed + last_page_num > max_num_pages {
            let new_map_size = (total_pages_needed + last_page_num) * page_size;
            self.env.set_mapsize(new_map_size).map_err(lmdb_error)?;
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
        let front = vol_names.first().map_or(0, |n| n.len()) as u64;
        let vol_name_size = (front + 24) * vol_names.len() as u64;
        let vol_info_size = 24 * vol_names.len() as u64;
        let page_max_size = self.env.psize() as u64 - 16;
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
        let size = self.total_ids_length as u64 + self.list.len() as u64 * 16;
        let avg_id_length = self.total_ids_length as u64 / self.list.len() as u64;
        let page_max_size = self.env.psize() as u64 - 16;
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
        // SKeyValuePair::cmp_key (writedb_lmdb.hpp:116-127, `kv_cmp`): `v.id < k.id`
        // (std::string, bytewise), ties by the OID's bytes from the least significant one.
        // Elements equal under it are equal in every byte that reaches a file (their
        // `saveToOidList` is read before the sort only), so an unstable sort gives NCBI's list.
        let ids = &self.ids;
        self.list.sort_unstable_by(|a, b| kv_cmp(ids, *a, *b));
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
        let mut key = std::mem::take(&mut self.scratch);
        let r = self.put_batch(dbi, i0, j, &mut key);
        self.scratch = key;
        r?;
        self.env.txn_commit().map_err(lmdb_error)
    }

    fn put_batch(
        &mut self,
        dbi: Dbi,
        i0: usize,
        j: usize,
        key: &mut Vec<u8>,
    ) -> Result<(), WriteDbError> {
        for i in i0..j {
            let kv = self.list[i];
            if i > 0 {
                let prev = self.list[i - 1];
                if id_cmp(&self.ids, prev, kv).is_eq() && prev.oid() == kv.oid() {
                    continue;
                }
            }
            let oid = kv.oid().to_le_bytes();
            kv.id_into(&self.ids, key);
            let ok = match self.env.put(dbi, c_str(key), &oid, MDB_APPENDDUP) {
                Ok(ok) => ok,
                Err(LmdbError::Mdb { func, code }) => {
                    return Err(acc2oid_put_error(func, code, key));
                }
                Err(e) => return Err(lmdb_error(e)),
            };
            if !ok {
                return Err(WriteDbError::Other(format!(
                    "acc2oid error for id {}",
                    String::from_utf8_lossy(key)
                )));
            }
        }
        Ok(())
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
        let total_num_oids = self.list.last().unwrap().oid() as u64 + 1;
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
        let mut tmp_ids: Vec<KeyValuePair> = Vec::new();
        let mut id = Vec::new();
        for i in 0..self.list.len() {
            let kv = self.list[i];
            let prev_oid = if i > 0 {
                self.list[i - 1].oid()
            } else {
                kv.oid()
            };
            if i > 0 && kv.oid() != prev_oid {
                if kv.oid().wrapping_sub(prev_oid) != 1 {
                    os.flush().map_err(io)?;
                    return Err(WriteDbError::Other(
                        "Input id list not in ascending oid order".to_string(),
                    ));
                }
                offsets[count] = self.write_ids(&mut os, &mut tmp_ids, &mut id).map_err(io)?;
                count += 1;
                tmp_ids.clear();
            }
            self.total_ids_length += kv.len();
            if !kv.save_to_oid_list() {
                continue;
            }
            tmp_ids.push(kv);
        }
        offsets[count] = self.write_ids(&mut os, &mut tmp_ids, &mut id).map_err(io)?;
        os.seek(SeekFrom::Start(8)).map_err(io)?;
        let mut offset: u64 = 0;
        for o in &offsets {
            offset += *o as u64;
            os.write_all(&offset.to_le_bytes()).map_err(io)?;
        }
        os.flush().map_err(io)
    }

    /// `s_WirteIds`: the ids of one OID, sorted, each after its length; `id` is scratch.
    fn write_ids(
        &self,
        os: &mut impl Write,
        ids: &mut [KeyValuePair],
        id: &mut Vec<u8>,
    ) -> std::io::Result<u32> {
        ids.sort_unstable_by(|a, b| id_cmp(&self.ids, *a, *b));
        let mut diff: u32 = 0;
        for e in ids.iter() {
            e.id_into(&self.ids, id);
            let id = &id[..];
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
    /// Create the LMDB file `dbname` (`<db>.ptf`/`<db>.ntf`) on the platform `platform`.
    pub(crate) fn new(
        dbname: &Path,
        map_size: u64,
        platform: &LmdbPlatform,
        settings: &LmdbEnvSettings,
    ) -> Result<Self, WriteDbError> {
        let env = Env::create(dbname, map_size, platform).map_err(lmdb_error)?;
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
        // NCBI's `size_t` arithmetic in `u64`, as in `WriteDbLmdb::grow_map`.
        let size = self.taxid2offsets.len() as u64 * 32;
        let page_size = self.env.psize() as u64;
        let page_max_size = page_size - 16;
        let last_page_num = self.env.last_pgno();
        let max_num_pages = self.env.mapsize() / page_size;
        let leaf_pages_needed = size / page_max_size + 1;
        let branch_pages_needed = 24 * leaf_pages_needed / page_max_size + 1;
        let total_pages_needed = leaf_pages_needed + branch_pages_needed + 4;
        if total_pages_needed + last_page_num > max_num_pages {
            let new_map_size = (total_pages_needed + last_page_num) * page_size;
            self.env.set_mapsize(new_map_size).map_err(lmdb_error)?;
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

#[cfg(test)]
#[path = "lmdb_keys_tests.rs"]
mod keys_tests;
