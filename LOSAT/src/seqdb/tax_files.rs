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

use crate::utils::glibc_fma_libm::GlibcLibm;
use std::collections::BTreeSet;

use super::lmdb_read::{Env, ReadAt, ReadError};

/// NCBI text when the `.ndb` has no `volname`/`volinfo` table.
pub const NO_VOLUME_INFO: &str = "DB contains no vol info.";
/// NCBI text when the `.ndb` has no `acc2oid` table (database made without `-parse_seqids`).
pub const NO_ACCESSION_INFO: &str = "DB contains no accession info.";

/// Text of NCBI's catch block of `GetOid` / `GetOids` for any LMDB error; `dbname` is the base
/// name of the `.ndb` file. (`MDB_NOTFOUND` there gives "Seqid list specified but no accession
/// table is found in `<dbname>`", which cannot happen once the table exists.)
pub fn accessions_lookup_error(dbname: &[u8]) -> Vec<u8> {
    [b"Accessions to Oids lookup error in ".as_slice(), dbname].concat()
}

fn le_u32(v: &[u8], what: &str) -> Result<u32, ReadError> {
    let a: [u8; 4] = v
        .try_into()
        .map_err(|_| format!("{what}: value of {} bytes, expected 4", v.len()))?;
    Ok(u32::from_le_bytes(a))
}

fn le_u64(v: &[u8], what: &str) -> Result<u64, ReadError> {
    let a: [u8; 8] = v
        .try_into()
        .map_err(|_| format!("{what}: value of {} bytes, expected 8", v.len()))?;
    Ok(u64::from_le_bytes(a))
}

fn read_u64_at<R: ReadAt + ?Sized>(src: &R, off: u64) -> Result<u64, ReadError> {
    let mut b = [0u8; 8];
    src.read_exact_at(&mut b, off)?;
    Ok(u64::from_le_bytes(b))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:124-146
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:306-357
/// Volumes of the entry in index order: `(volume base name, number of OIDs in that volume)`.
/// The OIDs of volume `i` start at the sum of the counts before it.
///
/// Errors: [`NO_VOLUME_INFO`], "Volinfo error " (the two tables have different entry counts),
/// "Volinfo error : N volumes" (LOSAT: a count larger than the tables' leaf pages can hold),
/// "No volinfo for `<name>`".
pub fn volumes<R: ReadAt + ?Sized>(ndb: &Env<R>) -> Result<Vec<(Vec<u8>, u32)>, ReadError> {
    let volname = ndb
        .open_db(b"volname")?
        .ok_or_else(|| NO_VOLUME_INFO.to_string())?;
    let volinfo = ndb
        .open_db(b"volinfo")?
        .ok_or_else(|| NO_VOLUME_INFO.to_string())?;
    if volinfo.entries != volname.entries {
        return Err("Volinfo error ".to_string().into());
    }
    // LOSAT: each entry is a leaf node of at least 14 bytes (8-byte node header, 4-byte key,
    // 2-byte pointer) on one of the table's leaf pages (`md_leaf_pages`, exact for these tables,
    // which have no duplicates; at most the pages of the file), and a page holds at most its
    // size less the 16-byte page header. A larger count is damage; NCBI would size `vol_names`
    // and `vol_num_oids` to it and look up every number.
    let leaf_pages = volname
        .leaf_pages
        .min(volinfo.leaf_pages)
        .min(ndb.num_pages());
    let room = leaf_pages.saturating_mul(ndb.page_size().saturating_sub(16) as u64) / 14;
    let n = u32::try_from(volinfo.entries)
        .ok()
        .filter(|&n| u64::from(n) <= room)
        .ok_or_else(|| format!("Volinfo error : {} volumes", volinfo.entries))?;
    let mut out = Vec::with_capacity((n as usize).min(1 << 16));
    for i in 0..n {
        let key = i.to_le_bytes();
        // A missing volname leaves the name empty and the count 0, as the resized vectors do.
        match ndb.first_value(&volname, &key)? {
            None => out.push((Vec::new(), 0)),
            Some(name) => {
                let Some(info) = ndb.first_value(&volinfo, &key)? else {
                    return Err(format!("No volinfo for {}", String::from_utf8_lossy(&name)).into());
                };
                out.push((name, le_u32(&info, "volinfo")?));
            }
        }
    }
    Ok(out)
}

/// Whether an error of [`volumes`] is one of NCBI's exceptions (the rest are errors of the file).
pub fn is_ncbi_volume_error(message: &str) -> bool {
    message == NO_VOLUME_INFO
        || message == "Volinfo error "
        || message.starts_with("No volinfo for ")
}

/// `m_NumOids` of NCBI: the sum of the OID counts of all volumes of the entry.
#[cfg(test)]
pub fn total_oids(vols: &[(Vec<u8>, u32)]) -> u64 {
    vols.iter().map(|v| v.1 as u64).sum()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:359-400
/// First OID of an accession key (`GetOids`: `MDB_SET`, `MDB_GET_CURRENT`; the first duplicate
/// in `memcmp` order of the 4 value bytes). The key is the id string exactly as stored
/// (`AAA00001`, `AAA00001.1`, `myprot9`, `pdb|...`): no case folding.
///
/// `Err(NO_ACCESSION_INFO)` when the table is missing; other `Err` strings are LMDB or format
/// errors, which NCBI reports as [`accessions_lookup_error`]. `Ok(None)`: id not in the table.
pub fn acc_first_oid<R: ReadAt + ?Sized>(
    ndb: &Env<R>,
    acc: &[u8],
) -> Result<Option<u32>, ReadError> {
    let db = ndb
        .open_db(b"acc2oid")?
        .ok_or_else(|| NO_ACCESSION_INFO.to_string())?;
    ndb.first_value(&db, acc)?
        .map(|v| le_u32(&v, "acc2oid"))
        .transpose()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:261-303
/// All OIDs of an accession key in LMDB duplicate order (`GetOid` with `allow_dup` true):
/// `memcmp` order of the 4 LE bytes, so {1, 256, 512} come as 256, 512, 1. Empty when absent.
/// Works for keys with a sub-database (more than about 500 OIDs), where NCBI fails.
pub fn acc_oids<R: ReadAt + ?Sized>(ndb: &Env<R>, acc: &[u8]) -> Result<Vec<u32>, ReadError> {
    acc_oids_of(ndb, &acc2oid(ndb)?, acc)
}

/// The `acc2oid` database of `ndb` (`mdb_dbi_open`; `NO_ACCESSION_INFO` where it has none), for
/// `acc_oids_of` and `acc_case_variant_oids`.
pub fn acc2oid<R: ReadAt + ?Sized>(ndb: &Env<R>) -> Result<crate::seqdb::lmdb_read::Db, ReadError> {
    ndb.open_db(b"acc2oid")?
        .ok_or_else(|| NO_ACCESSION_INFO.into())
}

/// `acc_oids` in the `acc2oid` database `db` of `ndb` (`acc2oid`).
pub fn acc_oids_of<R: ReadAt + ?Sized>(
    ndb: &Env<R>,
    db: &crate::seqdb::lmdb_read::Db,
    acc: &[u8],
) -> Result<Vec<u32>, ReadError> {
    ndb.values(db, acc)?
        .iter()
        .map(|v| le_u32(v, "acc2oid"))
        .collect()
}

// No NCBI counterpart: NCBI's lookup reads an accession key as given (`acc_oids`), and its scope
// takes text Seq-ids that differ in the case of their letters for one
// (`seqid::defline_set::scope_keys`); LOSAT lists the other keys to find the sequences that hold
// such an id (`seqdb::lookup::SeqidLookup::case_folded_oids`).
/// How many bytes after the fixed `from` bytes a prefix that `acc_case_variant_oids` keeps in its
/// map may have: at most 2 letters' worth, a few thousand prefixes.
const KEPT_CASE_PREFIX: usize = 2;

/// The OIDs of the keys of `acc2oid` other than `acc` that differ from it only in the case of
/// ASCII letters at or after byte `from`, in key order, each key's in LMDB duplicate order
/// (`acc_oids_of`, of the `acc2oid` database `db`), with the answers of `MDB_SET_RANGE` to whether
/// a key starts with a prefix of at most `KEPT_CASE_PREFIX` bytes after `from` kept in `prefixes`
/// for the next calls (`case_variant_keys`).
pub fn acc_case_variant_oids<R: ReadAt + ?Sized>(
    ndb: &Env<R>,
    db: &crate::seqdb::lmdb_read::Db,
    acc: &[u8],
    from: usize,
    prefixes: &mut rustc_hash::FxHashMap<Vec<u8>, bool>,
) -> Result<Vec<u32>, ReadError> {
    // Only the short prefixes recur from one name to the next (`pir||a` for every `pir||A12345`);
    // keeping the longer ones, which each name asks alone, would only grow the map.
    let kept_up_to = from.saturating_add(KEPT_CASE_PREFIX);
    let keys = case_variant_keys(acc, from, |prefix| -> Result<bool, ReadError> {
        if prefix.len() > kept_up_to {
            return ndb.has_key_starting_with(db, prefix);
        }
        if let Some(&known) = prefixes.get(prefix) {
            return Ok(known);
        }
        let found = ndb.has_key_starting_with(db, prefix)?;
        prefixes.insert(prefix.to_vec(), found);
        Ok(found)
    })?;
    let mut oids = Vec::new();
    for key in keys {
        oids.extend(acc_oids_of(ndb, db, &key)?);
    }
    Ok(oids)
}

/// The strings other than `acc` that differ from it only in the case of ASCII letters at or after
/// byte `from` and whose every prefix that leaves `acc`'s letters `has_prefix` finds (a key
/// starts with it), in byte order: the candidates for keys of a case variant of `acc`. Along
/// `acc`'s own letters nothing is asked: for each letter, the prefix that turns it to the other
/// case is asked once, and only where it is found does the search go on from there, each later
/// byte in both cases. So where no key differs in case it is asked once for each letter.
pub fn case_variant_keys<E>(
    acc: &[u8],
    from: usize,
    mut has_prefix: impl FnMut(&[u8]) -> Result<bool, E>,
) -> Result<Vec<Vec<u8>>, E> {
    /// The other case of an ASCII letter.
    fn flip(byte: u8) -> Option<u8> {
        if byte.is_ascii_uppercase() {
            Some(byte.to_ascii_lowercase())
        } else if byte.is_ascii_lowercase() {
            Some(byte.to_ascii_uppercase())
        } else {
            None
        }
    }
    /// The variants under `prefix`, which already leaves `acc`'s letters and is found.
    fn below<E>(
        acc: &[u8],
        prefix: &mut Vec<u8>,
        has_prefix: &mut dyn FnMut(&[u8]) -> Result<bool, E>,
        keys: &mut Vec<Vec<u8>>,
    ) -> Result<(), E> {
        let Some(&byte) = acc.get(prefix.len()) else {
            keys.push(prefix.clone());
            return Ok(());
        };
        for byte in [Some(byte), flip(byte)].into_iter().flatten() {
            prefix.push(byte);
            if has_prefix(prefix)? {
                below(acc, prefix, has_prefix, keys)?;
            }
            prefix.pop();
        }
        Ok(())
    }
    let mut keys = Vec::new();
    let mut prefix = Vec::new();
    for (at, &byte) in acc.iter().enumerate().skip(from) {
        let Some(other) = flip(byte) else {
            continue;
        };
        prefix.clear();
        prefix.extend_from_slice(&acc[..at]);
        prefix.push(other);
        if has_prefix(&prefix)? {
            below(acc, &mut prefix, &mut has_prefix, &mut keys)?;
        }
    }
    keys.sort_unstable();
    Ok(keys)
}

/// Failure of a taxid lookup, with the pieces for NCBI's message.
#[derive(Debug)]
pub enum TaxError {
    /// The `.ntf` has no `taxid2offset` table (`MDB_NOTFOUND` of `mdb_dbi_open`). NCBI:
    /// "No taxonomy info found in `<dbname>`".
    NoTaxInfo,
    /// Any other LMDB failure while reading the `.ntf` (corrupt file, bad value size). NCBI:
    /// "Taxonomy Id to Oids lookup error in `<dbname>`"; the string is the detail.
    Lookup(String),
    /// The `.nto` / `.not` data is unreadable or points outside the file or the OID range
    /// (a crash or a `CFileException` in NCBI). The string is the full message.
    Data(String),
    /// The system failed to give bytes of the `.ntf`, `.nto` or `.not` after it opened
    /// ([`ReadError::Unreadable`]): LOSAT's own error, as it is.
    Unreadable(super::SeqDbError),
}

impl TaxError {
    /// The failure of a read of the `.ntf` (`Lookup`), or the `Unreadable` error as it is.
    fn lookup(error: ReadError) -> Self {
        match error {
            ReadError::Format(message) => TaxError::Lookup(message),
            ReadError::Unreadable(error) => TaxError::Unreadable(error),
        }
    }

    /// The failure of a read of the `.nto` or `.not` (`Data`), or the `Unreadable` error as it is.
    fn data(error: ReadError) -> Self {
        match error {
            ReadError::Format(message) => TaxError::Data(message),
            ReadError::Unreadable(error) => TaxError::Unreadable(error),
        }
    }

    /// The message NCBI throws (`CSeqDBException`, `eArgErr`) for this failure; `dbname` is the
    /// base name of the `.ntf` file, e.g. `nt.ntf`. `Data` has no NCBI text and gives its detail.
    pub fn message(&self, dbname: &[u8]) -> Vec<u8> {
        match self {
            TaxError::NoTaxInfo => [b"No taxonomy info found in ".as_slice(), dbname].concat(),
            TaxError::Lookup(_) => {
                [b"Taxonomy Id to Oids lookup error in ".as_slice(), dbname].concat()
            }
            TaxError::Data(m) => m.clone().into_bytes(),
            TaxError::Unreadable(error) => error.to_string().into_bytes(),
        }
    }

    /// Detail of the failure for logs (no NCBI counterpart).
    #[cfg(test)]
    pub fn detail(&self) -> &str {
        match self {
            TaxError::NoTaxInfo => "no taxid2offset table",
            TaxError::Lookup(d) | TaxError::Data(d) => d,
            TaxError::Unreadable(_) => "unreadable",
        }
    }
}

impl std::fmt::Display for TaxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxError::NoTaxInfo => write!(f, "No taxonomy info found"),
            TaxError::Lookup(d) => write!(f, "Taxonomy Id to Oids lookup error: {d}"),
            TaxError::Data(d) => f.write_str(d),
            TaxError::Unreadable(error) => write!(f, "{error}"),
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
/// All taxids with a key in the `.ntf` (`GetDBTaxIds`), in key order (`memcmp` of the LE bytes,
/// i.e. unsorted as numbers).
pub fn get_db_tax_ids<R: ReadAt + ?Sized>(ntf: &Env<R>) -> Result<Vec<i32>, TaxError> {
    let db = ntf
        .open_db(b"taxid2offset")
        .map_err(TaxError::lookup)?
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
    .map_err(TaxError::lookup)?;
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:591-672
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
        .map_err(TaxError::lookup)?
        .ok_or(TaxError::NoTaxInfo)?;
    let mut offsets: Vec<u64> = Vec::new();
    let mut tax_ids_found: Vec<i32> = Vec::new();
    for &tax_id in tax_ids {
        let values = ntf
            .values(&db, &tax_id.to_le_bytes())
            .map_err(TaxError::lookup)?;
        if values.is_empty() {
            continue;
        }
        for v in &values {
            offsets.push(le_u64(v, "taxid2offset").map_err(TaxError::lookup)?);
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
                .map_err(|e| TaxError::data(e.prefixed("taxid-to-oids file: ")))?;
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
                .map_err(|e| TaxError::data(e.prefixed("taxid-to-oids file: ")))?;
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
    // NCBI: sort when `oids_sz` times the log of `oids_sz` is below `m_NumOids`, else scan the bit
    // vector. (For 0 the product is NaN and the scan finds nothing.) Both give ascending unique
    // OIDs.
    if (oids_sz as f64) * (oids_sz as f64).glibc_ln() < num_oids as f64 {
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
    fn open<R: ReadAt + ?Sized>(src: &R, what: &str) -> Result<OidIndex, ReadError> {
        let file_len = src.len();
        if file_len == 0 {
            return Err(format!("Failed to open oid-to-{what} lookup file").into());
        }
        let n_oids =
            read_u64_at(src, 0).map_err(|e| e.prefixed(&format!("oid-to-{what} lookup file: ")))?;
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
    ) -> Result<(u64, u64), ReadError> {
        if oid as u64 >= self.n_oids {
            return Err(format!(
                "oid-to-{what} lookup file: OID {oid} is not below the {} OIDs of the file",
                self.n_oids
            )
            .into());
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
            )
            .into());
        }
        Ok((begin, end))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:675-706
/// `CLookupTaxIds` over a `.not`/`.pot` file read by position. Keeps only the OID count.
pub struct OidTaxIds<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    idx: OidIndex,
}

impl<'r, R: ReadAt + ?Sized> OidTaxIds<'r, R> {
    /// Open and check the header (an empty file is NCBI's "Failed to open oid-to-taxids lookup
    /// file").
    pub fn open(src: &'r R) -> Result<Self, ReadError> {
        Ok(OidTaxIds {
            src,
            idx: OidIndex::open(src, "taxids")?,
        })
    }

    /// Number of OIDs in the file.
    #[cfg(test)]
    pub fn num_oids(&self) -> u64 {
        self.idx.n_oids
    }

    /// Taxids of one OID, as stored (ascending; 0 = unknown).
    pub fn taxids(&self, oid: u32) -> Result<Vec<i32>, ReadError> {
        let (begin, end) = self.idx.range(self.src, oid, "taxids")?;
        let len = (end - begin).checked_mul(4);
        let start = begin
            .checked_mul(4)
            .and_then(|b| b.checked_add(self.idx.data_start));
        let (Some(len), Some(start)) = (len, start) else {
            return Err(
                format!("oid-to-taxids lookup file: taxids of OID {oid} are out of range").into(),
            );
        };
        if start.checked_add(len).is_none_or(|e| e > self.idx.file_len) {
            return Err(format!(
                "oid-to-taxids lookup file: taxids of OID {oid} run past the end of the file"
            )
            .into());
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
#[cfg(test)]
pub fn taxids_of_oid<R: ReadAt + ?Sized>(not: &R, oid: u32) -> Result<Vec<i32>, ReadError> {
    OidTaxIds::open(not)?.taxids(oid)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:708-735
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
    let lookup = OidTaxIds::open(not).map_err(TaxError::data)?;
    let tax_id_list: BTreeSet<i32> = positive.tax_ids_found.iter().copied().collect();
    let mut rv = Vec::new();
    for &oid in &positive.oids {
        let file_list = lookup.taxids(oid).map_err(TaxError::data)?;
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
/// `CLookupSeqIds` over a `.nos`/`.pos` file read by position.
pub struct OidSeqIds<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    idx: OidIndex,
}

impl<'r, R: ReadAt + ?Sized> OidSeqIds<'r, R> {
    /// Open and check the header (an empty file is NCBI's "Failed to open oid-to-seqid lookup
    /// file").
    pub fn open(src: &'r R) -> Result<Self, ReadError> {
        Ok(OidSeqIds {
            src,
            idx: OidIndex::open(src, "seqid")?,
        })
    }

    /// Number of OIDs in the file.
    #[cfg(test)]
    pub fn num_oids(&self) -> u64 {
        self.idx.n_oids
    }

    /// Seqids of one OID in stored order (bytewise ascending): the versioned form where the id
    /// has a version (`AAA00001.1`), the plain id otherwise; GI numbers and names are absent.
    /// Each is `[u8 len][bytes]`, or `0xFF` + `u32 len` + bytes when 255 or longer.
    pub fn seqids(&self, oid: u32) -> Result<Vec<Vec<u8>>, ReadError> {
        let (begin, end) = self.idx.range(self.src, oid, "seqid")?;
        let start = self.idx.data_start.checked_add(begin);
        let len = end - begin;
        if start
            .and_then(|s| s.checked_add(len))
            .is_none_or(|e| e > self.idx.file_len)
        {
            return Err(format!(
                "oid-to-seqid lookup file: seqids of OID {oid} run past the end of the file"
            )
            .into());
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
#[cfg(test)]
pub fn seqids_of_oid<R: ReadAt + ?Sized>(nos: &R, oid: u32) -> Result<Vec<Vec<u8>>, ReadError> {
    OidSeqIds::open(nos)?.seqids(oid)
}

#[cfg(test)]
#[path = "tax_files_tests.rs"]
mod tests;
