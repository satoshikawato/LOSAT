//! NCBI layer over the v5 database lookup files of one LMDB entry.
//!
//! One entry is a `.ndb/.pdb` (accession and volume tables) with its `.ntf/.ptf`
//! (taxid -> offset), `.nto/.pto` (offset -> OID list), `.not/.pot` (OID -> taxids) and
//! `.nos/.pos` (OID -> seqids). All of them share the database's global OID numbering: OIDs
//! run over all volumes of the entry. The mapping of OIDs onto an alias or a set of entries
//! (`x_AdjustOidsOffset`) is not done here.
//!
//! Read-only. Integers are little-endian. Nothing is held in memory: every function reads the bytes it
//! needs by position. A value that NCBI would read out of bounds (undefined behaviour there) is
//! an `Err` here.

use std::collections::BTreeSet;

use super::lmdb_read::{Env, ReadAt};

/// NCBI text when the `.ndb` has no `volname`/`volinfo` table.
pub const NO_VOLUME_INFO: &str = "DB contains no vol info.";
/// NCBI text when the `.ndb` has no `acc2oid` table (database made without `-parse_seqids`).
pub const NO_ACCESSION_INFO: &str = "DB contains no accession info.";

/// Text of NCBI's catch block of `GetOid` / `GetOids` for any LMDB error; `dbname` is the base
/// name of the `.ndb` file. (`MDB_NOTFOUND` there gives "Seqid list specified but no accession
/// table is found in <dbname>", which cannot happen once the table exists.)
pub fn accessions_lookup_error(dbname: &str) -> String {
    format!("Accessions to Oids lookup error in {dbname}")
}

fn le_u32(v: &[u8], what: &str) -> Result<u32, String> {
    let a: [u8; 4] = v
        .try_into()
        .map_err(|_| format!("{what}: value of {} bytes, expected 4", v.len()))?;
    Ok(u32::from_le_bytes(a))
}

fn le_u64(v: &[u8], what: &str) -> Result<u64, String> {
    let a: [u8; 8] = v
        .try_into()
        .map_err(|_| format!("{what}: value of {} bytes, expected 8", v.len()))?;
    Ok(u64::from_le_bytes(a))
}

fn read_u64_at<R: ReadAt + ?Sized>(src: &R, off: u64) -> Result<u64, String> {
    let mut b = [0u8; 8];
    src.read_exact_at(&mut b, off)?;
    Ok(u64::from_le_bytes(b))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:124-146
// ```c++
// MDB_dbi CBlastLMDBManager::CBlastEnv::GetDbi(EDbiType dbi_type)
// {
//     if(m_dbis[dbi_type] == UINT_MAX) {
//         string err = "DB contains no ";
//         switch (dbi_type) {
//         case eDbiVolinof:
//         case eDbiVolname:
//             err += "vol info.";
//             break;
//         case eDbiAcc2oid:
//             err += "accession info.";
//             break;
//         case eDbiTaxid2offset:
//             err += "tax id info";
//             break;
//         default:
//             NCBI_THROW( CSeqDBException, eArgErr, "Invalid dbi type");
//         }
//         NCBI_THROW( CSeqDBException, eArgErr, err);
//
//     }
//     return m_dbis[dbi_type];
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:306-357
// ```c++
// void CSeqDBLMDB::GetVolumesInfo(vector<string> & vol_names, vector<blastdb::TOid> & vol_num_oids)
// {
//     MDB_dbi db_volname_handle;
//     MDB_dbi db_volinfo_handle;
//     lmdb::env & env = CBlastLMDBManager::GetInstance().GetReadEnvVol(m_LMDBFile, db_volname_handle, db_volinfo_handle);
//     vol_names.clear();
//     vol_num_oids.clear();
//     {
//     auto txn = lmdb::txn::begin(env, nullptr, MDB_RDONLY);
//     lmdb::dbi db_volname(db_volname_handle);
//     lmdb::dbi db_volinfo(db_volinfo_handle);
//     MDB_stat volinfo_stat, volname_stat;
//     lmdb::dbi_stat(txn, db_volinfo, &volinfo_stat);
//     lmdb::dbi_stat(txn, db_volname, &volname_stat);
//     if(volinfo_stat.ms_entries != volname_stat.ms_entries) {
//            NCBI_THROW( CSeqDBException, eArgErr, "Volinfo error ");
//     }
//
//     vol_names.resize(volinfo_stat.ms_entries);
//     vol_num_oids.resize(volinfo_stat.ms_entries);
//
//     auto cursor_volname = lmdb::cursor::open(txn, db_volname);
//     auto cursor_volinfo = lmdb::cursor::open(txn, db_volinfo);
//     for (unsigned int i=0; i < volinfo_stat.ms_entries; i++) {
//         lmdb::val data2find(&i, sizeof(Int4));
//         if (cursor_volname.get(data2find, MDB_SET)) {
//             {
//                 lmdb::val k, val;
//                 cursor_volname.get(k, val, MDB_GET_CURRENT);
//                 vol_names[i].assign(val.data(), val.size());
//             }
//             if (cursor_volinfo.get(data2find, MDB_SET)) {
//                 lmdb::val k, val;
//                 cursor_volinfo.get(k, val, MDB_GET_CURRENT);
//                 const char* d = val.data();
//                 vol_num_oids[i] = (((d[3] << 24)&0xFF000000) | ((d[2] << 16) & 0xFF0000) | ((d[1] << 8) & 0xFF00) | (d[0]&0xFF));
//             }
//             else {
//                 NCBI_THROW( CSeqDBException, eArgErr, "No volinfo for " + vol_names[i]);
//             }
//          }
//        }
//     cursor_volname.close();
//     cursor_volinfo.close();
//     txn.reset();
//     }
//     CBlastLMDBManager::GetInstance().CloseEnv(m_LMDBFile);
//
//     for(unsigned int j=0; j < vol_num_oids.size(); j++){
//         m_NumOids += vol_num_oids[j];
//     }
// }
// ```
/// Volumes of the entry in index order: `(volume base name, number of OIDs in that volume)`.
/// The OIDs of volume `i` start at the sum of the counts before it.
///
/// Errors: [`NO_VOLUME_INFO`], "Volinfo error " (the two tables have different entry counts),
/// "No volinfo for <name>".
pub fn volumes<R: ReadAt + ?Sized>(ndb: &Env<R>) -> Result<Vec<(Vec<u8>, u32)>, String> {
    let volname = ndb
        .open_db(b"volname")?
        .ok_or_else(|| NO_VOLUME_INFO.to_string())?;
    let volinfo = ndb
        .open_db(b"volinfo")?
        .ok_or_else(|| NO_VOLUME_INFO.to_string())?;
    if volinfo.entries != volname.entries {
        return Err("Volinfo error ".to_string());
    }
    let n = u32::try_from(volinfo.entries)
        .map_err(|_| format!("Volinfo error : {} volumes", volinfo.entries))?;
    let mut out = Vec::with_capacity((n as usize).min(1 << 16));
    for i in 0..n {
        let key = i.to_le_bytes();
        // A missing volname leaves the name empty and the count 0, as the resized vectors do.
        match ndb.first_value(&volname, &key)? {
            None => out.push((Vec::new(), 0)),
            Some(name) => {
                let Some(info) = ndb.first_value(&volinfo, &key)? else {
                    return Err(format!("No volinfo for {}", String::from_utf8_lossy(&name)));
                };
                out.push((name, le_u32(&info, "volinfo")?));
            }
        }
    }
    Ok(out)
}

/// `m_NumOids` of NCBI: the sum of the OID counts of all volumes of the entry.
pub fn total_oids(vols: &[(Vec<u8>, u32)]) -> u64 {
    vols.iter().map(|v| v.1 as u64).sum()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:359-400
// ```c++
// void
// CSeqDBLMDB::GetOids(const vector<string>& accessions, vector<blastdb::TOid>& oids) const
// {
//     try {
//     oids.clear();
//     oids.resize(accessions.size(), kSeqDBEntryNotFound);
//
//     MDB_dbi dbi_handle;
//     lmdb::env & env = CBlastLMDBManager::GetInstance().GetReadEnvAcc(m_LMDBFile, dbi_handle, &m_LMDBFileOpened);
//     {
//     lmdb::dbi dbi(dbi_handle);
//     auto txn = lmdb::txn::begin(env, nullptr, MDB_RDONLY);
//
//     auto cursor = lmdb::cursor::open(txn, dbi);
//
//     unsigned int i=0;
//     for (i=0; i < accessions.size(); i++) {
//         string acc = accessions[i];
//         lmdb::val data2find(acc);
//         if (cursor.get(data2find, MDB_SET)) {
//             lmdb::val k, val;
//             cursor.get(k, val, MDB_GET_CURRENT);
//             const char* d = val.data();
//             oids[i] = (((d[3] << 24)&0xFF000000) | ((d[2] << 16) & 0xFF0000) | ((d[1] << 8) & 0xFF00) | (d[0]&0xFF));
//         }
//     }
//
//     cursor.close();
//     txn.reset();
//     }
//     CBlastLMDBManager::GetInstance().CloseEnv(m_LMDBFile);
//     } catch (lmdb::error & e) {
//            string dbname;
//            CSeqDB_Path(m_LMDBFile).FindBaseName().GetString(dbname);
//            if(e.code() == MDB_NOTFOUND) {
//             NCBI_THROW( CSeqDBException, eArgErr, "Seqid list specified but no accession table is found in " + dbname);
//            }
//            else {
//             NCBI_THROW( CSeqDBException, eArgErr, "Accessions to Oids lookup error in " + dbname);
//            }
//     }
// }
// ```
/// First OID of an accession key (`GetOids`: `MDB_SET`, `MDB_GET_CURRENT`; the first duplicate
/// in `memcmp` order of the 4 value bytes). The key is the id string exactly as stored
/// (`AAA00001`, `AAA00001.1`, `myprot9`, `pdb|...`): no case folding.
///
/// `Err(NO_ACCESSION_INFO)` when the table is missing; other `Err` strings are LMDB or format
/// errors, which NCBI reports as [`accessions_lookup_error`]. `Ok(None)`: id not in the table.
pub fn acc_first_oid<R: ReadAt + ?Sized>(ndb: &Env<R>, acc: &[u8]) -> Result<Option<u32>, String> {
    let db = ndb
        .open_db(b"acc2oid")?
        .ok_or_else(|| NO_ACCESSION_INFO.to_string())?;
    ndb.first_value(&db, acc)?
        .map(|v| le_u32(&v, "acc2oid"))
        .transpose()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:261-303
// ```c++
// void
// CSeqDBLMDB::GetOid(const string & accession, vector<blastdb::TOid> & oids, const bool allow_dup) const
// {
//     try {
//     oids.clear();
//     {
//     MDB_dbi dbi_handle;
//     lmdb::env & env = CBlastLMDBManager::GetInstance().GetReadEnvAcc(m_LMDBFile, dbi_handle, &m_LMDBFileOpened);
//     lmdb::dbi dbi(dbi_handle);
//     auto txn = lmdb::txn::begin(env, nullptr, MDB_RDONLY);
//     auto cursor = lmdb::cursor::open(txn, dbi);
//
//     string acc = accession;
//     lmdb::val data2find(acc);
//
//     if (cursor.get(data2find, MDB_SET)) {
//         lmdb::val k, val;
//         cursor.get(k, val, MDB_GET_CURRENT);
//         const char* d = val.data();
//         oids.push_back(((d[3] << 24)&0xFF000000) | ((d[2] << 16) & 0xFF0000) | ((d[1] << 8) & 0xFF00) | (d[0]&0xFF));
//
//         if(allow_dup) {
//             while (cursor.get(k,val, MDB_NEXT_DUP)) {
//                 d = val.data();
//                 oids.push_back(((d[3] << 24)&0xFF000000) | ((d[2] << 16) & 0xFF0000) | ((d[1] << 8) & 0xFF00) | (d[0]&0xFF));
//             }
//         }
//     }
//     cursor.close();
//     txn.reset();
//     }
//     CBlastLMDBManager::GetInstance().CloseEnv(m_LMDBFile);
//     } catch (lmdb::error & e) {
//            string dbname;
//            CSeqDB_Path(m_LMDBFile).FindBaseName().GetString(dbname);
//            if(e.code() == MDB_NOTFOUND) {
//             NCBI_THROW( CSeqDBException, eArgErr, "Seqid list specified but no accession table is found in " + dbname);
//            }
//            else {
//             NCBI_THROW( CSeqDBException, eArgErr, "Accessions to Oids lookup error in " + dbname);
//            }
//     }
// }
// ```
/// All OIDs of an accession key in LMDB duplicate order (`GetOid(.., allow_dup = true)`):
/// `memcmp` order of the 4 LE bytes, so {1, 256, 512} come as 256, 512, 1. Empty when absent.
/// Works for keys with a sub-database (more than about 500 OIDs), where NCBI fails.
pub fn acc_oids<R: ReadAt + ?Sized>(ndb: &Env<R>, acc: &[u8]) -> Result<Vec<u32>, String> {
    let db = ndb
        .open_db(b"acc2oid")?
        .ok_or_else(|| NO_ACCESSION_INFO.to_string())?;
    ndb.values(&db, acc)?
        .iter()
        .map(|v| le_u32(v, "acc2oid"))
        .collect()
}

/// Failure of a taxid lookup, with the pieces for NCBI's message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaxError {
    /// The `.ntf` has no `taxid2offset` table (`MDB_NOTFOUND` of `mdb_dbi_open`). NCBI:
    /// "No taxonomy info found in <dbname>".
    NoTaxInfo,
    /// Any other LMDB failure while reading the `.ntf` (corrupt file, bad value size). NCBI:
    /// "Taxonomy Id to Oids lookup error in <dbname>"; the string is the detail.
    Lookup(String),
    /// The `.nto` / `.not` data is unreadable or points outside the file or the OID range
    /// (a crash or a `CFileException` in NCBI). The string is the full message.
    Data(String),
}

impl TaxError {
    /// The message NCBI throws (`CSeqDBException`, `eArgErr`) for this failure; `dbname` is the
    /// base name of the `.ntf` file, e.g. `nt.ntf`. `Data` has no NCBI text and gives its detail.
    pub fn message(&self, dbname: &str) -> String {
        match self {
            TaxError::NoTaxInfo => format!("No taxonomy info found in {dbname}"),
            TaxError::Lookup(_) => format!("Taxonomy Id to Oids lookup error in {dbname}"),
            TaxError::Data(m) => m.clone(),
        }
    }

    /// Detail of the failure for logs (no NCBI counterpart).
    pub fn detail(&self) -> &str {
        match self {
            TaxError::NoTaxInfo => "no taxid2offset table",
            TaxError::Lookup(d) | TaxError::Data(d) => d,
        }
    }
}

impl std::fmt::Display for TaxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxError::NoTaxInfo => write!(f, "No taxonomy info found"),
            TaxError::Lookup(d) => write!(f, "Taxonomy Id to Oids lookup error: {d}"),
            TaxError::Data(d) => f.write_str(d),
        }
    }
}

impl std::error::Error for TaxError {}

/// Result of [`get_oids_for_tax_ids`] and [`negative_tax_ids_to_oids`].
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TaxOids {
    /// Selected OIDs, ascending and unique (global OIDs of the entry).
    pub oids: Vec<u32>,
    /// The requested taxids that have a key in the `.ntf` (`tax_ids_found`), in the order of the
    /// request (ascending as `set<Int4>`).
    pub tax_ids_found: Vec<i32>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:558-589
// ```c++
// void CSeqDBLMDB::GetDBTaxIds(vector<TTaxId> & tax_ids) const
// {
//
//     tax_ids.clear();
//     try {
//        MDB_dbi dbi_handle;
//     lmdb::env & env = CBlastLMDBManager::GetInstance().GetReadEnvTax(m_TaxId2OffsetsFile, dbi_handle);
//     {
//         auto txn = lmdb::txn::begin(env, nullptr, MDB_RDONLY);
//         auto dbi(dbi_handle);
//         auto cursor = lmdb::cursor::open(txn, dbi);
//         lmdb::val key;
//         while (cursor.get(key, MDB_NEXT)) {
//             TTaxId taxid = TAX_ID_FROM(Int4, *((Int4 *)key.data()));
//             tax_ids.push_back(taxid);
//         }
//         cursor.close();
//         txn.reset();
//     }
//     }
//     catch (lmdb::error & e) {
//            string dbname;
//            CSeqDB_Path(m_TaxId2OffsetsFile).FindBaseName().GetString(dbname);
//            if(e.code() == MDB_NOTFOUND) {
//             NCBI_THROW( CSeqDBException, eArgErr, "No taxonomy info found in " + dbname);
//            }
//            else {
//             NCBI_THROW( CSeqDBException, eArgErr, "Taxonomy Id to Oids lookup error in " + dbname);
//            }
//     }
//     CBlastLMDBManager::GetInstance().CloseEnv(m_TaxId2OffsetsFile);
// }
// ```
/// All taxids with a key in the `.ntf` (`GetDBTaxIds`), in key order (`memcmp` of the LE bytes,
/// i.e. unsorted as numbers).
pub fn get_db_tax_ids<R: ReadAt + ?Sized>(ntf: &Env<R>) -> Result<Vec<i32>, TaxError> {
    let db = ntf
        .open_db(b"taxid2offset")
        .map_err(TaxError::Lookup)?
        .ok_or(TaxError::NoTaxInfo)?;
    let mut out = Vec::new();
    ntf.for_each(&db, |k, _| {
        let a: [u8; 4] = k
            .try_into()
            .map_err(|_| format!("taxid2offset key of {} bytes, expected 4", k.len()))?;
        // MDB_NEXT visits one key once per duplicate; the writer makes none.
        out.push(i32::from_le_bytes(a));
        Ok(())
    })
    .map_err(TaxError::Lookup)?;
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:591-672
// ```c++
// void CSeqDBLMDB::GetOidsForTaxIds(const set<TTaxId> & tax_ids, vector<blastdb::TOid>& oids, vector<TTaxId> & tax_ids_found) const
// {
//     try {
//     oids.clear();
//     tax_ids_found.clear();
//     vector<Uint8> offsets;
//        MDB_dbi dbi_handle;
//     lmdb::env & env = CBlastLMDBManager::GetInstance().GetReadEnvTax(m_TaxId2OffsetsFile, dbi_handle);
//     {
//     auto txn = lmdb::txn::begin(env, nullptr, MDB_RDONLY);
//        lmdb::dbi dbi(dbi_handle);
//     auto cursor = lmdb::cursor::open(txn, dbi);
//     ITERATE(set<TTaxId>, itr, tax_ids) {
//         Int4 tax_id = TAX_ID_TO(Int4, *itr);
//         lmdb::val data2find(tax_id);
//
//         if (cursor.get(data2find, MDB_SET)) {
//             lmdb::val k, val;
//             cursor.get(k, val, MDB_GET_CURRENT);
//             const char* d = val.data();
//             offsets.push_back((((Uint8) d[7] << 56) &0xFF00000000000000) | (((Uint8) d[6] << 48) & 0xFF000000000000) |
//                               (((Uint8) d[5] << 40) &0xFF0000000000)     | (((Uint8) d[4] << 32) & 0xFF00000000)     |
//                               (((Uint8) d[3] << 24) &0xFF000000)         | (((Uint8) d[2] << 16) & 0xFF0000)         |
//                               (((Uint8) d[1] << 8)  &0xFF00)             | ((Uint8)  d[0]&0xFF));
//             while (cursor.get(k,val, MDB_NEXT_DUP)) {
//                 d = val.data();
//                 offsets.push_back((((Uint8) d[7] << 56) &0xFF00000000000000) | (((Uint8) d[6] << 48) & 0xFF000000000000) |
//                                   (((Uint8) d[5] << 40) &0xFF0000000000)     | (((Uint8) d[4] << 32) & 0xFF00000000)     |
//                                   (((Uint8) d[3] << 24) &0xFF000000)         | (((Uint8) d[2] << 16) & 0xFF0000)         |
//                                   (((Uint8) d[1] << 8)  &0xFF00)             | ((Uint8)  d[0]&0xFF));
//             }
//             tax_ids_found.push_back(*itr);
//         }
//     }
//     cursor.close();
//     txn.reset();
//     }
//     CBlastLMDBManager::GetInstance().CloseEnv(m_TaxId2OffsetsFile);
//     vector<bool> oids_set(m_NumOids, false);
//
//        CMemoryFile oid_file(m_TaxId2OidsFile);
//        const char * start_ptr = (char *) oid_file.GetPtr();
//        for (unsigned int i=0; i < offsets.size(); i++) {
//            Uint4 * list_ptr = (Uint4 * ) (start_ptr + offsets[i]);
//            Uint4 num_of_oids = *list_ptr;
//            Uint4 count = 0 ;
//            list_ptr ++;
//            while(count < num_of_oids) {
//                if(!oids_set[*list_ptr]) {
//                    oids.push_back(*list_ptr);
//                    oids_set[*list_ptr] = true;
//                }
//                count++;
//                list_ptr++;
//            }
//        }
//
//        int oids_sz = oids.size();
//        if((oids_sz*log(oids_sz)) < m_NumOids) {
//            sort(oids.begin(), oids.end());
//        }
//        else {
//            oids.resize(0);
//            oids.reserve(oids_sz);
//            for (int i=0; i < oids_set.size(); i++) {
//                if(oids_set[i]) {
//                 oids.push_back(i);
//             }
//         }
//        }
//
//     } catch (lmdb::error & e) {
//            string dbname;
//            CSeqDB_Path(m_TaxId2OffsetsFile).FindBaseName().GetString(dbname);
//            if(e.code() == MDB_NOTFOUND) {
//             NCBI_THROW( CSeqDBException, eArgErr, "No taxonomy info found in " + dbname);
//            }
//            else {
//             NCBI_THROW( CSeqDBException, eArgErr, "Taxonomy Id to Oids lookup error in " + dbname);
//            }
//     }
// }
// ```
/// `CSeqDBLMDB::GetOidsForTaxIds`: the union of the OID lists of the requested taxids.
///
/// * `ntf`: the `.ntf`/`.ptf` environment (table `taxid2offset`, key = Int4 taxid in 4 LE
///   bytes, value(s) = u64 offset into the `.nto`).
/// * `nto`: the `.nto`/`.pto` bytes: `[u32 n][u32 oid x n]` at each offset.
/// * `num_oids`: `m_NumOids`, the OID count of the entry (`total_oids(volumes(..))`); an OID at
///   or above it in a list is `TaxError::Data` (undefined behaviour in NCBI).
/// * `tax_ids`: the requested taxids, a sorted set as NCBI's `set<Int4>`.
///
/// Every duplicate of a taxid key counts (the writer makes one value per key). The OIDs come out
/// ascending and unique, whichever of NCBI's two sort strategies applies.
pub fn get_oids_for_tax_ids<R1: ReadAt + ?Sized, R2: ReadAt + ?Sized>(
    ntf: &Env<R1>,
    nto: &R2,
    num_oids: u32,
    tax_ids: &BTreeSet<i32>,
) -> Result<TaxOids, TaxError> {
    let db = ntf
        .open_db(b"taxid2offset")
        .map_err(TaxError::Lookup)?
        .ok_or(TaxError::NoTaxInfo)?;
    let mut offsets: Vec<u64> = Vec::new();
    let mut tax_ids_found: Vec<i32> = Vec::new();
    for &tax_id in tax_ids {
        let values = ntf
            .values(&db, &tax_id.to_le_bytes())
            .map_err(TaxError::Lookup)?;
        if values.is_empty() {
            continue;
        }
        for v in &values {
            offsets.push(le_u64(v, "taxid2offset").map_err(TaxError::Lookup)?);
        }
        tax_ids_found.push(tax_id);
    }

    let mut oids_set = vec![0u64; (num_oids as usize).div_ceil(64)];
    let mut oids: Vec<u32> = Vec::new();
    let nto_len = nto.len();
    for &offset in &offsets {
        let count = {
            let mut b = [0u8; 4];
            nto.read_exact_at(&mut b, offset)
                .map_err(|e| TaxError::Data(format!("taxid-to-oids file: {e}")))?;
            u32::from_le_bytes(b) as u64
        };
        let list_end = count.checked_mul(4).and_then(|n| n.checked_add(offset + 4));
        if list_end.is_none_or(|e| e > nto_len) {
            return Err(TaxError::Data(format!(
                "taxid-to-oids file: list of {count} OIDs at offset {offset} runs past the end of the file ({nto_len} bytes)"
            )));
        }
        // Read the list in blocks of up to 16384 OIDs.
        let mut done = 0u64;
        let mut buf = vec![0u8; 4 * (count.min(16384) as usize)];
        while done < count {
            let take = (count - done).min(16384) as usize;
            nto.read_exact_at(&mut buf[..take * 4], offset + 4 + done * 4)
                .map_err(|e| TaxError::Data(format!("taxid-to-oids file: {e}")))?;
            for c in buf[..take * 4].chunks_exact(4) {
                let oid = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                if oid >= num_oids {
                    return Err(TaxError::Data(format!(
                        "taxid-to-oids file: OID {oid} is not below the {num_oids} sequences of the database"
                    )));
                }
                let (w, bit) = (oid as usize / 64, oid % 64);
                if (oids_set[w] >> bit) & 1 == 0 {
                    oids.push(oid);
                    oids_set[w] |= 1 << bit;
                }
            }
            done += take as u64;
        }
    }

    let oids_sz = oids.len();
    // NCBI: sort when `oids_sz * log(oids_sz) < m_NumOids`, else scan the bit vector. (For 0 the
    // product is NaN and the scan finds nothing.) Both give ascending unique OIDs.
    if (oids_sz as f64) * (oids_sz as f64).ln() < num_oids as f64 {
        oids.sort_unstable();
    } else {
        oids.clear();
        oids.reserve(oids_sz);
        for (w, &word) in oids_set.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let bit = bits.trailing_zeros();
                oids.push((w * 64) as u32 + bit);
                bits &= bits - 1;
            }
        }
    }
    Ok(TaxOids {
        oids,
        tax_ids_found,
    })
}

/// Start of the data and the cumulative-end index of an `.not` / `.nos` file:
/// `u64 n_oids`, `n_oids x u64` ends, then the data.
struct OidIndex {
    n_oids: u64,
    data_start: u64,
    file_len: u64,
}

impl OidIndex {
    fn open<R: ReadAt + ?Sized>(src: &R, what: &str) -> Result<OidIndex, String> {
        let file_len = src.len();
        if file_len == 0 {
            return Err(format!("Failed to open oid-to-{what} lookup file"));
        }
        let n_oids = read_u64_at(src, 0).map_err(|e| format!("oid-to-{what} lookup file: {e}"))?;
        let data_start = n_oids
            .checked_add(1)
            .and_then(|n| n.checked_mul(8))
            .filter(|&s| s <= file_len)
            .ok_or_else(|| format!("oid-to-{what} lookup file: {n_oids} OIDs do not fit a file of {file_len} bytes"))?;
        Ok(OidIndex {
            n_oids,
            data_start,
            file_len,
        })
    }

    /// `[begin, end)` of an OID in the units of the file (Int4 for `.not`, bytes for `.nos`).
    fn range<R: ReadAt + ?Sized>(
        &self,
        src: &R,
        oid: u32,
        what: &str,
    ) -> Result<(u64, u64), String> {
        if oid as u64 >= self.n_oids {
            return Err(format!(
                "oid-to-{what} lookup file: OID {oid} is not below the {} OIDs of the file",
                self.n_oids
            ));
        }
        let at = 8 + 8 * oid as u64;
        let (begin, end) = if oid == 0 {
            (0, read_u64_at(src, at)?)
        } else {
            let mut b = [0u8; 16];
            src.read_exact_at(&mut b, at - 8)?;
            (
                u64::from_le_bytes(b[..8].try_into().map_err(|_| "short read")?),
                u64::from_le_bytes(b[8..].try_into().map_err(|_| "short read")?),
            )
        };
        if begin > end {
            return Err(format!(
                "oid-to-{what} lookup file: ends of OID {oid} go backwards ({begin} > {end})"
            ));
        }
        Ok((begin, end))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:675-706
// ```c++
// class CLookupTaxIds
// {
// public:
//     CLookupTaxIds(CMemoryFile & file): m_IndexStart((Uint8*) file.GetPtr()), m_DataStart((Int4 *) file.GetPtr()) {
//         if(m_IndexStart == NULL){
//             NCBI_THROW( CSeqDBException, eArgErr, "Failed to open oid-to-taxids lookup file");
//         }
//
//         Uint8 num_of_oids = *m_IndexStart;
//         m_IndexStart ++;
//         m_DataStart += (2* (num_of_oids + 1));
//     }
//
//     inline void GetTaxIdListForOid(blastdb::TOid oid, vector<TTaxId> & taxid_list);
// private:
//
//     Uint8 * m_IndexStart;
//     Int4 * m_DataStart;
// };
//
// void CLookupTaxIds::GetTaxIdListForOid(blastdb::TOid oid, vector<TTaxId> & taxid_list)
// {
//     taxid_list.clear();
//     Uint8 * index_ptr = m_IndexStart + oid;
//     Int4 * end = m_DataStart + (*index_ptr);
//     index_ptr--;
//     Int4 * begin = (oid == 0) ? m_DataStart:m_DataStart + (*index_ptr);
//     while (begin < end) {
//         taxid_list.push_back(TAX_ID_FROM(Int4, *begin));
//         begin++;
//     }
// }
// ```
/// `CLookupTaxIds` over a `.not`/`.pot` file read by position. Keeps only the OID count.
pub struct OidTaxIds<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    idx: OidIndex,
}

impl<'r, R: ReadAt + ?Sized> OidTaxIds<'r, R> {
    /// Open and check the header (an empty file is NCBI's "Failed to open oid-to-taxids lookup
    /// file").
    pub fn open(src: &'r R) -> Result<Self, String> {
        Ok(OidTaxIds {
            src,
            idx: OidIndex::open(src, "taxids")?,
        })
    }

    /// Number of OIDs in the file.
    pub fn num_oids(&self) -> u64 {
        self.idx.n_oids
    }

    /// Taxids of one OID, as stored (ascending; 0 = unknown).
    pub fn taxids(&self, oid: u32) -> Result<Vec<i32>, String> {
        let (begin, end) = self.idx.range(self.src, oid, "taxids")?;
        let len = (end - begin).checked_mul(4);
        let start = begin
            .checked_mul(4)
            .and_then(|b| b.checked_add(self.idx.data_start));
        let (Some(len), Some(start)) = (len, start) else {
            return Err(format!(
                "oid-to-taxids lookup file: taxids of OID {oid} are out of range"
            ));
        };
        if start.checked_add(len).is_none_or(|e| e > self.idx.file_len) {
            return Err(format!(
                "oid-to-taxids lookup file: taxids of OID {oid} run past the end of the file"
            ));
        }
        let mut buf = vec![0u8; len as usize];
        self.src.read_exact_at(&mut buf, start)?;
        Ok(buf
            .chunks_exact(4)
            .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect())
    }
}

/// Taxids of one OID from a `.not`/`.pot` file (`CLookupTaxIds::GetTaxIdListForOid`). For many
/// OIDs use [`OidTaxIds`], which reads the header once.
pub fn taxids_of_oid<R: ReadAt + ?Sized>(not: &R, oid: u32) -> Result<Vec<i32>, String> {
    OidTaxIds::open(not)?.taxids(oid)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:708-735
// ```c++
// void
// CSeqDBLMDB::NegativeTaxIdsToOids(const set<TTaxId>& tax_ids, vector<blastdb::TOid>& rv, vector<TTaxId> & tax_ids_found) const
// {
//     rv.clear();
//     vector<blastdb::TOid> oids;
//     GetOidsForTaxIds(tax_ids, oids, tax_ids_found);
//     CMemoryFile oid_file(m_Oid2TaxIdsFile);
//     set<TTaxId> tax_id_list(tax_ids_found.begin(), tax_ids_found.end());
//     CLookupTaxIds lookup(oid_file);
//     for(unsigned int i=0; i < oids.size(); i++) {
//         vector<TTaxId>  file_list;
//         lookup.GetTaxIdListForOid(oids[i], file_list);
//         if(file_list.size() > tax_ids.size()) {
//             continue;
//         }
//         else {
//             unsigned int j = 0;
//             for(; j < file_list.size(); j++) {
//                 if(tax_id_list.find(file_list[j]) == tax_id_list.end()) {
//                     break;
//                 }
//             }
//             if(j == file_list.size()) {
//                 rv.push_back(oids[i]);
//             }
//         }
//     }
// }
// ```
/// `CSeqDBLMDB::NegativeTaxIdsToOids`: the OIDs that `-negative_taxids` removes.
///
/// Candidates are the OIDs of [`get_oids_for_tax_ids`]. A candidate is selected (= excluded
/// from the search) only if it has at most `tax_ids.len()` taxids in the `.not` and every one of
/// them is among the found taxids. So an OID with taxids {A, B} survives `-negative_taxids A`.
/// Note that `tax_ids.len()` is the size of the request, not of the found set.
/// `oids` is ascending; `tax_ids_found` as in [`get_oids_for_tax_ids`].
pub fn negative_tax_ids_to_oids<R1: ReadAt + ?Sized, R2: ReadAt + ?Sized, R3: ReadAt + ?Sized>(
    ntf: &Env<R1>,
    nto: &R2,
    not: &R3,
    num_oids: u32,
    tax_ids: &BTreeSet<i32>,
) -> Result<TaxOids, TaxError> {
    let positive = get_oids_for_tax_ids(ntf, nto, num_oids, tax_ids)?;
    let lookup = OidTaxIds::open(not).map_err(TaxError::Data)?;
    let tax_id_list: BTreeSet<i32> = positive.tax_ids_found.iter().copied().collect();
    let mut rv = Vec::new();
    for &oid in &positive.oids {
        let file_list = lookup.taxids(oid).map_err(TaxError::Data)?;
        if file_list.len() > tax_ids.len() {
            continue;
        }
        if file_list.iter().all(|t| tax_id_list.contains(t)) {
            rv.push(oid);
        }
    }
    Ok(TaxOids {
        oids: rv,
        tax_ids_found: positive.tax_ids_found,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:415-459
// ```c++
// class CLookupSeqIds
// {
// public:
//     CLookupSeqIds(CMemoryFile & file): m_IndexStart((Uint8*) file.GetPtr()), m_DataStart((char *) file.GetPtr()) {
//         if(m_IndexStart == NULL){
//             NCBI_THROW( CSeqDBException, eArgErr, "Failed to open oid-to-seqid lookup file");
//         }
//
//         Uint8 num_of_oids = *m_IndexStart;
//         m_IndexStart ++;
//         m_DataStart += (8 * (num_of_oids + 1));
//     }
//
//     inline void GetSeqIdListForOid(blastdb::TOid oid, vector<string> & idlist);
// private:
//
//     Uint8 * m_IndexStart;
//     char * m_DataStart;
// };
//
// void CLookupSeqIds::GetSeqIdListForOid(blastdb::TOid oid, vector<string> & idlist)
// {
//     Uint8 * index_ptr = m_IndexStart + oid;
//     Char * end = m_DataStart + (*index_ptr);
//     index_ptr--;
//     Char * begin = (oid == 0) ? m_DataStart:m_DataStart + (*index_ptr);
//     while (begin < end) {
//         unsigned char id_len = *begin;
//         begin ++;
//         if(id_len == 0xFF) {
//             Uint4 long_id_len = *((Uint4 *) begin);
//             begin +=4;
//             string id;
//             id.assign(begin, long_id_len);
//             begin += long_id_len;
//             idlist.push_back(id);
//         }
//         else {
//             string id;
//             id.assign(begin, id_len);
//             begin += id_len;
//             idlist.push_back(id);
//         }
//     }
// }
// ```
/// `CLookupSeqIds` over a `.nos`/`.pos` file read by position.
pub struct OidSeqIds<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    idx: OidIndex,
}

impl<'r, R: ReadAt + ?Sized> OidSeqIds<'r, R> {
    /// Open and check the header (an empty file is NCBI's "Failed to open oid-to-seqid lookup
    /// file").
    pub fn open(src: &'r R) -> Result<Self, String> {
        Ok(OidSeqIds {
            src,
            idx: OidIndex::open(src, "seqid")?,
        })
    }

    /// Number of OIDs in the file.
    pub fn num_oids(&self) -> u64 {
        self.idx.n_oids
    }

    /// Seqids of one OID in stored order (bytewise ascending): the versioned form where the id
    /// has a version (`AAA00001.1`), the plain id otherwise; GI numbers and names are absent.
    /// Each is `[u8 len][bytes]`, or `0xFF` + `u32 len` + bytes when 255 or longer.
    pub fn seqids(&self, oid: u32) -> Result<Vec<Vec<u8>>, String> {
        let (begin, end) = self.idx.range(self.src, oid, "seqid")?;
        let start = self.idx.data_start.checked_add(begin);
        let len = end - begin;
        if start
            .and_then(|s| s.checked_add(len))
            .is_none_or(|e| e > self.idx.file_len)
        {
            return Err(format!(
                "oid-to-seqid lookup file: seqids of OID {oid} run past the end of the file"
            ));
        }
        let mut buf = vec![0u8; len as usize];
        self.src.read_exact_at(&mut buf, start.unwrap_or(0))?;
        let mut ids = Vec::new();
        let mut p = 0usize;
        while p < buf.len() {
            let mut len = buf[p] as usize;
            p += 1;
            if len == 0xFF {
                let b = buf.get(p..p + 4).ok_or_else(|| {
                    format!("oid-to-seqid lookup file: truncated id length of OID {oid}")
                })?;
                len = u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize;
                p += 4;
            }
            let id = p
                .checked_add(len)
                .and_then(|e| buf.get(p..e))
                .ok_or_else(|| format!("oid-to-seqid lookup file: id of {len} bytes of OID {oid} runs past its end"))?;
            ids.push(id.to_vec());
            p += len;
        }
        Ok(ids)
    }
}

/// Seqids of one OID from a `.nos`/`.pos` file (`CLookupSeqIds::GetSeqIdListForOid`). For many
/// OIDs use [`OidSeqIds`].
pub fn seqids_of_oid<R: ReadAt + ?Sized>(nos: &R, oid: u32) -> Result<Vec<Vec<u8>>, String> {
    OidSeqIds::open(nos)?.seqids(oid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;

    /// Fixture databases (made by NCBI makeblastdb 2.17.0, see `tests/seqdb/README.md`).
    const FX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb");
    /// Values written by `dump.py` (a Python reference reader) and by NCBI programs.
    const EXP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb/expected");

    /// Path prefix of the NCBI-distributed `18S_fungal_sequences` database
    /// (`LOSAT_TEST_18S_DB`, without extension) for the `#[ignore]` tests.
    fn s18() -> String {
        std::env::var("LOSAT_TEST_18S_DB").expect("LOSAT_TEST_18S_DB")
    }

    /// Path of an `expected/<name>.txt` file: the repository fixtures, else the directory named by
    /// `LOSAT_TEST_18S_EXPECTED` (the dumps of the 18S database are too big to commit).
    fn exp_path(name: &str) -> String {
        let repo = format!("{EXP}/{name}.txt");
        if std::path::Path::new(&repo).exists() {
            return repo;
        }
        let dir = std::env::var("LOSAT_TEST_18S_EXPECTED").unwrap_or_default();
        format!("{dir}/{name}.txt")
    }

    /// Values dumped by `gen/expect_tax.py` with `dump.py`.
    #[derive(Default)]
    struct Expected {
        vols: Vec<(String, u32)>,
        /// taxid -> OIDs (ascending)
        t: Vec<(i32, Vec<u32>)>,
        /// oid -> taxids
        o: Vec<Vec<i32>>,
        /// oid -> seqids
        s: Vec<Vec<String>>,
    }

    fn csv<T: std::str::FromStr>(s: &str) -> Vec<T>
    where
        T::Err: std::fmt::Debug,
    {
        s.split(',')
            .filter(|x| !x.is_empty())
            .map(|x| x.parse().unwrap())
            .collect()
    }

    fn expected(name: &str) -> Expected {
        let text = std::fs::read_to_string(exp_path(&format!("{name}.tax"))).unwrap();
        let mut e = Expected::default();
        for line in text.lines() {
            let f: Vec<&str> = line.splitn(5, ' ').collect();
            match f[0] {
                "V" => e.vols.push((f[1].to_string(), f[2].parse().unwrap())),
                "T" => e.t.push((f[1].parse().unwrap(), csv(f[3]))),
                "O" => {
                    assert_eq!(f[1].parse::<usize>().unwrap(), e.o.len());
                    e.o.push(csv(f.get(2).copied().unwrap_or("")));
                }
                "S" => {
                    assert_eq!(f[1].parse::<usize>().unwrap(), e.s.len());
                    let ids = f.get(2).copied().unwrap_or("");
                    e.s.push(
                        ids.split('|')
                            .filter(|x| !x.is_empty())
                            .map(String::from)
                            .collect(),
                    );
                }
                "X" => {}
                _ => panic!("bad line {line}"),
            }
        }
        e
    }

    /// The five files of one LMDB entry.
    struct Files {
        ndb: File,
        ntf: File,
        nto: File,
        not: File,
        nos: Option<File>,
    }

    fn files(prefix: &str, n: char) -> Files {
        let o = |ext: &str| File::open(format!("{prefix}.{n}{ext}")).unwrap();
        Files {
            ndb: o("db"),
            ntf: o("tf"),
            nto: o("to"),
            not: o("ot"),
            nos: File::open(format!("{prefix}.{n}os")).ok(),
        }
    }

    fn set(v: &[i32]) -> BTreeSet<i32> {
        v.iter().copied().collect()
    }

    /// (path prefix, p/n, expected-file name)
    fn all_dbs() -> Vec<(String, char, &'static str)> {
        vec![
            (format!("{FX}/n5"), 'n', "n5"),
            (format!("{FX}/p5"), 'p', "p5"),
            (format!("{FX}/mv"), 'p', "mv"),
            (format!("{FX}/dup"), 'p', "dup"),
            (format!("{FX}/tx"), 'p', "tx"),
        ]
    }

    fn check_volumes(prefix: &str, n: char, name: &str) {
        {
            let f = files(prefix, n);
            let env = Env::open(&f.ndb).unwrap();
            let vols = volumes(&env).unwrap();
            let want = expected(name).vols;
            let got: Vec<(String, u32)> = vols
                .iter()
                .map(|(a, b)| (String::from_utf8(a.clone()).unwrap(), *b))
                .collect();
            assert_eq!(got, want, "{name}");
            // the sum is the OID count of the .not file
            let oidtax = OidTaxIds::open(&f.not).unwrap();
            assert_eq!(total_oids(&vols), oidtax.num_oids(), "{name}");
        }
    }

    #[test]
    fn volumes_match_dump_py() {
        for (prefix, n, name) in all_dbs() {
            check_volumes(&prefix, n, name);
        }
        // literals from the notes: mv has 4 volumes 8,8,8,6
        let f = files(&format!("{FX}/mv"), 'p');
        let vols = volumes(&Env::open(&f.ndb).unwrap()).unwrap();
        assert_eq!(
            vols,
            vec![
                (b"mv.00".to_vec(), 8),
                (b"mv.01".to_vec(), 8),
                (b"mv.02".to_vec(), 8),
                (b"mv.03".to_vec(), 6)
            ]
        );
    }

    #[test]
    fn volumes_of_400_volume_db() {
        let f = File::open(format!("{FX}/mvbig.pdb")).unwrap();
        let env = Env::open(&f).unwrap();
        let vols = volumes(&env).unwrap();
        assert_eq!(vols.len(), 400);
        for (i, (name, _)) in vols.iter().enumerate() {
            assert_eq!(name, format!("mvbig.{i:03}").as_bytes());
        }
        // counts from dump.py's volinfo dump
        let text = std::fs::read_to_string(exp_path("mvbig.pdb")).unwrap();
        let mut in_info = false;
        let mut want = Vec::new();
        for l in text.lines() {
            if l.starts_with("DB ") {
                in_info = l.starts_with("DB volinfo ");
            } else if in_info && l.starts_with("K ") {
                let v = l.split(' ').nth(2).unwrap();
                want.push(u32::from_le_bytes([
                    u8::from_str_radix(&v[0..2], 16).unwrap(),
                    u8::from_str_radix(&v[2..4], 16).unwrap(),
                    u8::from_str_radix(&v[4..6], 16).unwrap(),
                    u8::from_str_radix(&v[6..8], 16).unwrap(),
                ]));
            }
        }
        assert_eq!(vols.iter().map(|v| v.1).collect::<Vec<_>>(), want);
    }

    #[test]
    fn missing_tables_give_ncbi_texts() {
        // a .ptf has no volume or accession table; a .pdb has no taxid table
        let ptf = File::open(format!("{FX}/p5.ptf")).unwrap();
        let env = Env::open(&ptf).unwrap();
        assert_eq!(volumes(&env).unwrap_err(), "DB contains no vol info.");
        assert_eq!(
            acc_first_oid(&env, b"x").unwrap_err(),
            "DB contains no accession info."
        );
        assert_eq!(acc_oids(&env, b"x").unwrap_err(), NO_ACCESSION_INFO);
        let pdb = File::open(format!("{FX}/p5.pdb")).unwrap();
        let env = Env::open(&pdb).unwrap();
        let nto: Vec<u8> = Vec::new();
        let e = get_oids_for_tax_ids(&env, &nto, 7, &set(&[1])).unwrap_err();
        assert_eq!(e, TaxError::NoTaxInfo);
        assert_eq!(e.message("p5.ptf"), "No taxonomy info found in p5.ptf");
        assert_eq!(get_db_tax_ids(&env).unwrap_err(), TaxError::NoTaxInfo);
        assert_eq!(
            accessions_lookup_error("p5.pdb"),
            "Accessions to Oids lookup error in p5.pdb"
        );
        let e = TaxError::Lookup("detail".into());
        assert_eq!(
            e.message("x.ntf"),
            "Taxonomy Id to Oids lookup error in x.ntf"
        );
        assert_eq!(e.detail(), "detail");
        // volinfo/volname with different counts and a missing volinfo: synthetic LMDB is covered in lmdb_read
    }

    #[test]
    fn accession_lookups() {
        let f = File::open(format!("{FX}/p5.pdb")).unwrap();
        let env = Env::open(&f).unwrap();
        assert_eq!(acc_first_oid(&env, b"AAA00001").unwrap(), Some(2));
        assert_eq!(acc_oids(&env, b"AAA00001").unwrap(), vec![2, 3]);
        assert_eq!(acc_first_oid(&env, b"AAA00001.2").unwrap(), Some(3));
        assert_eq!(acc_first_oid(&env, b"aaa00001").unwrap(), None);
        assert_eq!(acc_oids(&env, b"nokey").unwrap(), Vec::<u32>::new());
        // 700 OIDs under one key (sub-database): memcmp order of the LE bytes, NCBI fails here
        let f = File::open(format!("{FX}/dup.pdb")).unwrap();
        let env = Env::open(&f).unwrap();
        let oids = acc_oids(&env, b"DUP00001").unwrap();
        assert_eq!(oids.len(), 700);
        assert_eq!(&oids[..7], &[0, 256, 512, 1, 257, 513, 2]);
        assert_eq!(acc_first_oid(&env, b"DUP00001").unwrap(), Some(0));
        // global OIDs over volumes (notes: MVA00017.1 -> 17, volume 2 holds 16..23)
        let f = File::open(format!("{FX}/mv.pdb")).unwrap();
        let env = Env::open(&f).unwrap();
        assert_eq!(acc_first_oid(&env, b"MVA00017.1").unwrap(), Some(17));
    }

    /// Every acc2oid key of the dump of `path` (`name` names the `expected` file).
    fn check_accession_keys(path: &str, name: &str) {
        {
            let f = File::open(path).unwrap();
            let env = Env::open(&f).unwrap();
            let text = std::fs::read_to_string(exp_path(name)).unwrap();
            let mut in_acc = false;
            let mut n = 0;
            for l in text.lines() {
                if l.starts_with("DB ") {
                    in_acc = l.starts_with("DB acc2oid ");
                } else if in_acc && l.starts_with("K ") {
                    let p: Vec<&str> = l.split(' ').collect();
                    let key: Vec<u8> = (0..p[1].len() / 2)
                        .map(|i| u8::from_str_radix(&p[1][2 * i..2 * i + 2], 16).unwrap())
                        .collect();
                    let want: Vec<u32> = p[2]
                        .split(',')
                        .map(|v| {
                            u32::from_str_radix(
                                &format!("{}{}{}{}", &v[6..8], &v[4..6], &v[2..4], &v[0..2]),
                                16,
                            )
                            .unwrap()
                        })
                        .collect();
                    assert_eq!(
                        acc_first_oid(&env, &key).unwrap(),
                        Some(want[0]),
                        "{name} {}",
                        String::from_utf8_lossy(&key)
                    );
                    assert_eq!(
                        acc_oids(&env, &key).unwrap(),
                        want,
                        "{name} {}",
                        String::from_utf8_lossy(&key)
                    );
                    n += 1;
                }
            }
            assert!(n > 0, "{name}");
        }
    }

    /// Every acc2oid key of the dump (all fixtures).
    #[test]
    fn every_accession_key_matches_dump_py() {
        for name in ["p5.pdb", "n5.ndb", "mv.pdb", "dup.pdb", "thr501.pdb"] {
            check_accession_keys(&format!("{FX}/{name}"), name);
        }
    }

    #[test]
    fn n5_literals_from_the_notes() {
        let f = files(&format!("{FX}/n5"), 'n');
        let ntf = Env::open(&f.ntf).unwrap();
        let ndb = Env::open(&f.ndb).unwrap();
        let num = total_oids(&volumes(&ndb).unwrap()) as u32;
        assert_eq!(num, 5);
        // n5.nto: 562 -> {2,4}; 7227 -> {0,3,4}; 10090 -> {1}
        let r = |t: &[i32]| get_oids_for_tax_ids(&ntf, &f.nto, num, &set(t)).unwrap();
        assert_eq!(r(&[562]).oids, vec![2, 4]);
        assert_eq!(r(&[7227]).oids, vec![0, 3, 4]);
        assert_eq!(r(&[10090]).oids, vec![1]);
        let both = r(&[7227, 562, 9606]);
        assert_eq!(both.oids, vec![0, 2, 3, 4]);
        assert_eq!(both.tax_ids_found, vec![562, 7227]);
        assert_eq!(r(&[9606]), TaxOids::default());
        // n5.not: oid 4 has two taxids (two deflines)
        assert_eq!(taxids_of_oid(&f.not, 0).unwrap(), vec![7227]);
        assert_eq!(taxids_of_oid(&f.not, 4).unwrap(), vec![562, 7227]);
        assert!(taxids_of_oid(&f.not, 5).is_err());
        // n5.nos: OID 4 = DD000001.1, Y00002.1
        let nos = f.nos.as_ref().unwrap();
        assert_eq!(
            seqids_of_oid(nos, 4).unwrap(),
            vec![b"DD000001.1".to_vec(), b"Y00002.1".to_vec()]
        );
        assert_eq!(seqids_of_oid(nos, 0).unwrap(), vec![b"AB000001.1".to_vec()]);
        assert!(seqids_of_oid(nos, 5).is_err());
        // negative: {562}: OID 2 only (OID 4 also has 7227); {562, 7227}: all four
        let neg = |t: &[i32]| negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &set(t)).unwrap();
        assert_eq!(neg(&[562]).oids, vec![2]);
        assert_eq!(neg(&[562, 7227]).oids, vec![0, 2, 3, 4]);
        // size rule: 2 requested ids, found {562}; OID 4 has {562, 7227} (2 <= 2) but 7227 not found
        let n = neg(&[562, 1]);
        assert_eq!((n.oids, n.tax_ids_found), (vec![2], vec![562]));
        assert_eq!(neg(&[1]), TaxOids::default());
    }

    /// Every OID list, taxid list and seqid list of every fixture against the dump.py values.
    fn check_lookup_files(prefix: &str, n: char, name: &str) {
        {
            let f = files(&prefix, n);
            let ex = expected(name);
            let ndb = Env::open(&f.ndb).unwrap();
            let ntf = Env::open(&f.ntf).unwrap();
            let num = total_oids(&volumes(&ndb).unwrap()) as u32;
            // GetDBTaxIds
            let mut keys = get_db_tax_ids(&ntf).unwrap();
            let want_keys: Vec<i32> = ex.t.iter().map(|t| t.0).collect();
            assert_eq!(keys, want_keys, "{name} key order");
            // one taxid at a time
            for (tax, oids) in &ex.t {
                let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &set(&[*tax])).unwrap();
                assert_eq!(&r.oids, oids, "{name} taxid {tax}");
                assert_eq!(r.tax_ids_found, vec![*tax]);
            }
            // all taxids at once: the union (a big list takes the bit-vector branch)
            keys.sort_unstable();
            let all =
                get_oids_for_tax_ids(&ntf, &f.nto, num, &keys.iter().copied().collect()).unwrap();
            let mut union: Vec<u32> = ex.t.iter().flat_map(|t| t.1.iter().copied()).collect();
            union.sort_unstable();
            union.dedup();
            assert_eq!(all.oids, union, "{name} union");
            assert_eq!(all.tax_ids_found, keys);
            // per OID
            let oidtax = OidTaxIds::open(&f.not).unwrap();
            assert_eq!(oidtax.num_oids(), ex.o.len() as u64);
            for (oid, want) in ex.o.iter().enumerate() {
                assert_eq!(
                    &oidtax.taxids(oid as u32).unwrap(),
                    want,
                    "{name} taxids of {oid}"
                );
            }
            if let Some(nos) = &f.nos {
                let seq = OidSeqIds::open(nos).unwrap();
                assert_eq!(seq.num_oids(), ex.s.len() as u64);
                for (oid, want) in ex.s.iter().enumerate() {
                    let got: Vec<String> = seq
                        .seqids(oid as u32)
                        .unwrap()
                        .into_iter()
                        .map(|b| String::from_utf8(b).unwrap())
                        .collect();
                    assert_eq!(&got, want, "{name} seqids of {oid}");
                }
            }
        }
    }

    /// Every OID list, taxid list and seqid list of every fixture against the dump.py values.
    #[test]
    fn lookup_files_match_dump_py() {
        for (prefix, n, name) in all_dbs() {
            check_lookup_files(&prefix, n, name);
        }
    }

    fn parse_oracle(name: &str) -> Vec<(char, Vec<i32>, Option<Vec<u32>>)> {
        let text = std::fs::read_to_string(exp_path(name)).unwrap();
        text.lines()
            .map(|l| {
                let p: Vec<&str> = l.splitn(3, ' ').collect();
                let rest = p.get(2).copied().unwrap_or("");
                let oids = if rest.starts_with("Error") {
                    None
                } else {
                    Some(
                        rest.split_whitespace()
                            .map(|x| x.parse().unwrap())
                            .collect(),
                    )
                };
                (p[0].chars().next().unwrap(), csv(p[1]), oids)
            })
            .collect()
    }

    /// NCBI itself (blastdbcmd -taxids for the positive rule, blastp -negative_taxids with one
    /// query per OID for the negative rule), on a DB with OIDs that have 1, 2, 3 taxids, taxid 0
    /// and a taxid listed twice. `gen/oracle.sh` wrote the file.
    #[test]
    fn positive_and_negative_rules_match_ncbi_tx_fixture() {
        let f = files(&format!("{FX}/tx"), 'p');
        let ntf = Env::open(&f.ntf).unwrap();
        let num = total_oids(&volumes(&Env::open(&f.ndb).unwrap()).unwrap()) as u32;
        assert_eq!(num, 14);
        let mut checked = (0, 0);
        for (kind, taxids, oids) in parse_oracle("tx_oracle") {
            let want_set = set(&taxids);
            match kind {
                'P' => {
                    let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &want_set).unwrap();
                    // NCBI's blastdbcmd errors "Taxonomy ID(s) not found" when nothing matched
                    assert_eq!(Some(r.oids), oids.or(Some(vec![])), "P {taxids:?}");
                    checked.0 += 1;
                }
                'N' => {
                    // the oracle lists the OIDs still searched; the port returns the excluded ones
                    let excluded = negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &want_set)
                        .unwrap()
                        .oids;
                    let included: Vec<u32> = (0..num).filter(|o| !excluded.contains(o)).collect();
                    assert_eq!(included, oids.unwrap(), "N {taxids:?}");
                    checked.1 += 1;
                }
                _ => panic!(),
            }
        }
        assert_eq!(checked, (16, 16));
    }

    /// Needs the NCBI-distributed `18S_fungal_sequences` database (`LOSAT_TEST_18S_DB`, path without
    /// extension). The expected OIDs are what `blastdbcmd -taxids` printed.
    #[test]
    #[ignore = "needs LOSAT_TEST_18S_DB"]
    fn positive_rule_matches_ncbi_on_distributed_18s() {
        let f = files(&s18(), 'n');
        let ntf = Env::open(&f.ntf).unwrap();
        let num = total_oids(&volumes(&Env::open(&f.ndb).unwrap()).unwrap()) as u32;
        let rows = parse_oracle("18S_oracle");
        assert!(rows.len() >= 10);
        for (_, taxids, oids) in rows {
            let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &set(&taxids)).unwrap();
            assert_eq!(r.oids, oids.unwrap_or_default(), "{taxids:?}");
        }
    }

    /// The negative rule re-derived from the dump.py tables on random sets, for every fixture
    /// (independent code: candidates from the T table, per-OID lists from the O table).
    fn check_negative_rule(prefix: &str, n: char, name: &str, next: &mut impl FnMut() -> usize) {
        {
            let f = files(prefix, n);
            let ex = expected(name);
            let ntf = Env::open(&f.ntf).unwrap();
            let num = ex.o.len() as u32;
            for round in 0..60 {
                let size = 1 + next() % 5;
                let mut req = BTreeSet::new();
                for _ in 0..size {
                    // mostly present taxids, sometimes an absent one
                    if next() % 5 == 0 {
                        req.insert(1_000_000 + (next() % 50) as i32);
                    } else {
                        req.insert(ex.t[next() % ex.t.len()].0);
                    }
                }
                let found: BTreeSet<i32> =
                    ex.t.iter()
                        .map(|t| t.0)
                        .filter(|t| req.contains(t))
                        .collect();
                let mut cand: Vec<u32> =
                    ex.t.iter()
                        .filter(|t| found.contains(&t.0))
                        .flat_map(|t| t.1.iter().copied())
                        .collect();
                cand.sort_unstable();
                cand.dedup();
                let want: Vec<u32> = cand
                    .into_iter()
                    .filter(|&o| {
                        let l = &ex.o[o as usize];
                        l.len() <= req.len() && l.iter().all(|t| found.contains(t))
                    })
                    .collect();
                let got = negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &req).unwrap();
                assert_eq!(got.oids, want, "{name} round {round} {req:?}");
                assert_eq!(got.tax_ids_found, found.iter().copied().collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn negative_rule_on_random_sets() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        for (prefix, n, name) in all_dbs() {
            check_negative_rule(&prefix, n, name, &mut next);
        }
    }

    #[test]
    fn corrupt_lookup_files_are_errors() {
        let f = files(&format!("{FX}/n5"), 'n');
        let ntf = Env::open(&f.ntf).unwrap();
        // offsets past the end of the .nto
        let short: Vec<u8> = std::fs::read(format!("{FX}/n5.nto")).unwrap()[..20].to_vec();
        let e = get_oids_for_tax_ids(&ntf, &short, 5, &set(&[10090])).unwrap_err();
        assert!(
            matches!(&e, TaxError::Data(m) if m.contains("beyond the end")),
            "{e:?}"
        );
        // a list that starts inside the file but runs past its end
        let e = get_oids_for_tax_ids(&ntf, &short, 5, &set(&[7227])).unwrap_err();
        assert!(
            matches!(&e, TaxError::Data(m) if m.contains("runs past the end")),
            "{e:?}"
        );
        assert!(get_oids_for_tax_ids(&ntf, &Vec::<u8>::new(), 5, &set(&[562])).is_err());
        // an OID at or above the OID count (undefined behaviour in NCBI)
        let e = get_oids_for_tax_ids(&ntf, &f.nto, 4, &set(&[562])).unwrap_err();
        assert!(
            matches!(&e, TaxError::Data(m) if m.contains("OID 4 is not below the 4")),
            "{e:?}"
        );
        // a huge count
        let mut bad = std::fs::read(format!("{FX}/n5.nto")).unwrap();
        bad[0..4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(get_oids_for_tax_ids(&ntf, &bad, 5, &set(&[562])).is_err());
        // .not: empty file, truncated, bad header, ends going backwards
        assert_eq!(
            OidTaxIds::open(&Vec::<u8>::new()).err().unwrap(),
            "Failed to open oid-to-taxids lookup file"
        );
        assert_eq!(
            OidSeqIds::open(&Vec::<u8>::new()).err().unwrap(),
            "Failed to open oid-to-seqid lookup file"
        );
        let not = std::fs::read(format!("{FX}/n5.not")).unwrap();
        assert!(OidTaxIds::open(&not[..12]).is_err());
        let t = OidTaxIds::open(&not[..not.len() - 4]).unwrap();
        assert!(t.taxids(4).is_err());
        assert!(t.taxids(3).is_ok());
        let mut b = not.clone();
        b[8..16].copy_from_slice(&5u64.to_le_bytes());
        b[16..24].copy_from_slice(&1u64.to_le_bytes());
        assert!(OidTaxIds::open(&b)
            .unwrap()
            .taxids(1)
            .unwrap_err()
            .contains("go backwards"));
        let mut b = not.clone();
        b[0..8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(OidTaxIds::open(&b).is_err());
        // the negative rule reports a bad .not as Data
        let e = negative_tax_ids_to_oids(&ntf, &f.nto, &not[..40].to_vec(), 5, &set(&[562]))
            .unwrap_err();
        assert!(matches!(e, TaxError::Data(_)), "{e:?}");
    }

    #[test]
    fn seqid_list_forms() {
        // 3 OIDs: [two ids], [empty], [one id of 300 bytes written as 0xFF + u32 length]
        let long = vec![b'x'; 300];
        let mut data: Vec<u8> = Vec::new();
        data.extend([3, b'a', b'b', b'c', 2, b'd', b'e']);
        let end0 = data.len() as u64;
        let end1 = end0;
        data.push(0xFF);
        data.extend(300u32.to_le_bytes());
        data.extend(&long);
        let end2 = data.len() as u64;
        let mut file: Vec<u8> = Vec::new();
        file.extend(3u64.to_le_bytes());
        for e in [end0, end1, end2] {
            file.extend(e.to_le_bytes());
        }
        file.extend(&data);
        let s = OidSeqIds::open(&file).unwrap();
        assert_eq!(s.seqids(0).unwrap(), vec![b"abc".to_vec(), b"de".to_vec()]);
        assert_eq!(s.seqids(1).unwrap(), Vec::<Vec<u8>>::new());
        assert_eq!(s.seqids(2).unwrap(), vec![long]);
        assert!(s.seqids(3).is_err());
        // an id that runs past its OID's end
        let mut bad = file.clone();
        bad[8..16].copy_from_slice(&5u64.to_le_bytes());
        assert!(OidSeqIds::open(&bad)
            .unwrap()
            .seqids(0)
            .unwrap_err()
            .contains("runs past"));
        // truncated file
        assert!(OidSeqIds::open(&file[..file.len() - 10])
            .unwrap()
            .seqids(2)
            .is_err());
        // a 0xFF length cut off
        let mut cut = file.clone();
        cut[24..32].copy_from_slice(&(end0 + 3).to_le_bytes());
        cut[16..24].copy_from_slice(&end0.to_le_bytes());
        assert!(OidSeqIds::open(&cut)
            .unwrap()
            .seqids(2)
            .unwrap_err()
            .contains("truncated"));
    }

    /// Needs the NCBI-distributed `18S_fungal_sequences` database (`LOSAT_TEST_18S_DB`, path without
    /// extension).
    #[test]
    #[ignore = "needs LOSAT_TEST_18S_DB"]
    fn nos_of_distributed_18s_has_uncovered_tail() {
        // The NCBI-distributed 18S .nos has 12 bytes after the last end: the last OID reads empty.
        let f = File::open(format!("{}.nos", s18())).unwrap();
        let s = OidSeqIds::open(&f).unwrap();
        let n = s.num_oids() as u32;
        assert_eq!(s.seqids(n - 1).unwrap(), Vec::<Vec<u8>>::new());
        assert_eq!(s.seqids(0).unwrap(), vec![b"NG_065155.1".to_vec()]);
    }

    /// Volumes, every lookup list, every accession key and the negative rule of the NCBI-distributed
    /// `18S_fungal_sequences` database against the `dump.py` outputs. Needs `LOSAT_TEST_18S_DB` (path
    /// without extension) and `LOSAT_TEST_18S_EXPECTED` (directory with `18S.tax.txt` and
    /// `18S_fungal_sequences.ndb.txt`).
    #[test]
    #[ignore = "needs LOSAT_TEST_18S_DB and LOSAT_TEST_18S_EXPECTED"]
    fn ncbi_distributed_18s_lookup_files_match_dump_py() {
        let prefix = s18();
        check_volumes(&prefix, 'n', "18S");
        check_lookup_files(&prefix, 'n', "18S");
        check_accession_keys(&format!("{prefix}.ndb"), "18S_fungal_sequences.ndb");
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        check_negative_rule(&prefix, 'n', "18S", &mut next);
    }

    #[test]
    fn many_oids_in_one_list_cross_the_read_blocks() {
        // 20000 OIDs in one list (more than one 16384 block), unsorted input OIDs with a duplicate
        let n = 20000u32;
        let mut nto: Vec<u8> = Vec::new();
        nto.extend((n + 1).to_le_bytes());
        for i in (0..n).rev() {
            nto.extend(i.to_le_bytes());
        }
        nto.extend(7u32.to_le_bytes());
        // a minimal .ntf is needed: borrow p5.ptf's layout is not possible, so build the lookup by
        // hand through the dumped fixture dup.ptf (taxid 0 -> offset 0)
        let f = File::open(format!("{FX}/dup.ptf")).unwrap();
        let ntf = Env::open(&f).unwrap();
        let r = get_oids_for_tax_ids(&ntf, &nto, n, &set(&[0])).unwrap();
        assert_eq!(r.oids, (0..n).collect::<Vec<_>>());
        assert_eq!(r.tax_ids_found, vec![0]);
    }

    /// Read back a version 5 database written by LOSAT's own writer (`writedb::lmdb`): five proteins
    /// in two volumes, two taxids on some, none on one (the writer files it under taxid 0, as
    /// makeblastdb does for a defline without a map entry).
    #[test]
    fn reads_back_a_database_written_by_the_losat_writer() {
        use crate::writedb::lmdb::{lmdb_file_name, LmdbEnvSettings, WriteDbLmdb, WriteDbTaxId};
        use crate::writedb::platform::{LmdbPlatform, Platform};

        let linux = Platform::LINUX_X86_64;
        let platform = LmdbPlatform::LINUX_X86_64;
        let dir = std::env::temp_dir().join(format!("losat-tax-files-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // (accession, tax ids as written, tax ids as read back)
        let oids: [(&str, Vec<i32>); 5] = [
            ("NP_000001", vec![9606, 10090]),
            ("NP_000002", vec![9606]),
            ("NP_000003", vec![]),
            ("NP_000004", vec![562, 10090, 65536]),
            ("NP_000005", vec![562]),
        ];
        let pdb = dir.join("db.pdb");
        let ptf = lmdb_file_name(&pdb, "tf");
        let settings = LmdbEnvSettings::default();
        let mut lmdb =
            WriteDbLmdb::new(&pdb, linux.default_lmdb_map_size, &platform, &settings).unwrap();
        let mut tax =
            WriteDbTaxId::new(&ptf, linux.default_taxid_map_size, &platform, &settings).unwrap();
        for (oid, (acc, taxids)) in oids.iter().enumerate() {
            lmdb.insert_entry(acc.as_bytes(), oid as u32, false);
            lmdb.insert_entry(format!("{acc}.1").as_bytes(), oid as u32, true);
            tax.insert_entries(taxids, oid as u32);
        }
        lmdb.insert_volumes_info(&["db.00".to_string(), "db.01".to_string()], &[3, 2])
            .unwrap();
        lmdb.close().unwrap();
        tax.close().unwrap();

        let open = |ext: &str| File::open(dir.join(format!("db.{ext}"))).unwrap();
        let (ndb, ntf, nto, not, nos) = (
            open("pdb"),
            open("ptf"),
            open("pto"),
            open("pot"),
            open("pos"),
        );
        let db = Env::open(&ndb).unwrap();
        let tf = Env::open(&ntf).unwrap();

        let vols = volumes(&db).unwrap();
        assert_eq!(vols, vec![(b"db.00".to_vec(), 3), (b"db.01".to_vec(), 2)]);
        let num = total_oids(&vols) as u32;
        assert_eq!(num, 5);

        // what was written, as a read-back expects it: no tax ids = taxid 0
        let written: Vec<Vec<i32>> = oids
            .iter()
            .map(|o| if o.1.is_empty() { vec![0] } else { o.1.clone() })
            .collect();
        let mut all: Vec<i32> = written.iter().flatten().copied().collect();
        all.sort_unstable();
        all.dedup();
        let mut keys = get_db_tax_ids(&tf).unwrap();
        keys.sort_unstable();
        assert_eq!(keys, all);
        for &t in &all {
            let want: Vec<u32> = (0..num)
                .filter(|&o| written[o as usize].contains(&t))
                .collect();
            let got = get_oids_for_tax_ids(&tf, &nto, num, &set(&[t])).unwrap();
            assert_eq!((got.oids, got.tax_ids_found), (want, vec![t]), "taxid {t}");
        }
        let both = get_oids_for_tax_ids(&tf, &nto, num, &set(&[10090, 562, 7])).unwrap();
        assert_eq!(
            (both.oids, both.tax_ids_found),
            (vec![0, 3, 4], vec![562, 10090])
        );
        for (oid, want) in written.iter().enumerate() {
            assert_eq!(&taxids_of_oid(&not, oid as u32).unwrap(), want, "OID {oid}");
        }
        assert!(taxids_of_oid(&not, num).is_err());
        // negative rule: excluding 9606 drops OID 1 only (OID 0 also has 10090)
        assert_eq!(
            negative_tax_ids_to_oids(&tf, &nto, &not, num, &set(&[9606]))
                .unwrap()
                .oids,
            vec![1]
        );
        // accessions and seqids
        assert_eq!(acc_first_oid(&db, b"NP_000004.1").unwrap(), Some(3));
        assert_eq!(acc_oids(&db, b"NP_000002").unwrap(), vec![1]);
        for (oid, (acc, _)) in oids.iter().enumerate() {
            assert_eq!(
                seqids_of_oid(&nos, oid as u32).unwrap(),
                vec![format!("{acc}.1").into_bytes()],
                "OID {oid}"
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
