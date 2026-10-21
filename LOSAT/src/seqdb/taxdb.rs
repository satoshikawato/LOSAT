//! Reader of NCBI's taxonomy name database (`taxdb.bti` and `taxdb.btd`).
//!
//! NCBI BLAST+ writes the scientific, common and BLAST names of a taxid, and
//! the super kingdom, from this pair of files (`CSeqDBTaxInfo::GetTaxNames`,
//! c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp): the `sciname` of the
//! XML2 and JSON reports (`-outfmt 13` to `16`), the tabular fields `sscinames`,
//! `scomnames`, `sblastnames`, `sskingdoms` (and their singular forms) and the
//! organism report (`-outfmt 18`). This module is the reader only; the
//! report code decides when a name is looked up and what it writes when the
//! lookup finds nothing.
//!
//! ## Files
//!
//! Both files are found as NCBI finds them, with [`SearchPath::resolve_file`]
//! (`SeqDB_ResolveDbPath("taxdb.bti")`: the first of the current directory,
//! `BLASTDB` and the registry value `[BLAST] BLASTDB` that holds a regular file of
//! that name; the data file is the same path ending in `d`). All integers
//! are big-endian 32-bit.
//!
//! * `taxdb.bti`: magic `0x8739`, the number of taxids, four reserved words,
//!   then that many pairs (taxid, offset into `taxdb.btd`), sorted by taxid.
//! * `taxdb.btd`: the records back to back, `scientific\tcommon\tblast\tkingdom`;
//!   a record ends where the next entry's offset begins (the last: at the end
//!   of the file). The kingdom is the rest of the record, so a stray byte at the
//!   end of a record belongs to it.
//!
//! ## What is reproduced and what is refused
//!
//! The search and the split into four fields are NCBI's, including its
//! results for damaged files that are defined in practice: a file that is
//! missing, too short, of another magic number or of a count that does not
//! match its size makes every lookup fail (the last two with a line NCBI
//! writes through `ERR_POST`, see [`TaxDb::open_messages`]); a record with
//! fewer than three tabs, or with an empty kingdom, is "not found"; an
//! entry whose record would end before it begins is "not found".
//!
//! What NCBI does with a damaged file without a defined result is refused with
//! [`SeqDbError::Unsupported`] (it reads outside the file or crashes): a
//! record that reaches past the last 4096-byte page of `taxdb.btd` (within
//! that page NCBI reads the zero padding of its mapping, and so does this
//! reader), which includes every record of an empty `taxdb.btd` (NCBI
//! dereferences a null pointer). Lookups
//! of such a file that read no record (a taxid outside the index, a record
//! that ends before it begins, the offset error of the last entry) give
//! NCBI's result.
//!
//! The reader opens nothing until a name is asked for ([`LazyTaxDb`]); it
//! reads the index by position (about 22 reads of 8 bytes per lookup), so it
//! holds no table in memory.

use std::path::PathBuf;
use std::sync::OnceLock;

use super::file::DbFile;
use super::{SearchPath, SeqDbError};

/// Header of `taxdb.bti`: magic, count, four reserved words (seqdbtax.cpp,
/// `data_start`).
const HEADER_LEN: u64 = 4 + 4 + 16;
/// One (taxid, offset) pair of `taxdb.bti` (`CSeqDBTaxId`).
const ENTRY_LEN: u64 = 8;
/// The page size of the oracle platform (Linux x86-64): NCBI's mapping of `taxdb.btd` reads zero
/// bytes from the end of the file to the end of its last page.
const MAP_PAGE: u64 = 4096;
/// seqdbtax.cpp: `TAX_DB_MAGIC_NUMBER`.
const MAGIC: u32 = 0x8739;

/// `ERR_POST` text (seqdbtax.cpp:`CTaxDBFileInfo::CTaxDBFileInfo`): the magic number is wrong.
pub const MSG_WRONG_MAGIC: &str = "Error: Tax database file has wrong magic number.";
/// `ERR_POST` text (seqdbtax.cpp:`GetTaxNames`): the last entry begins after the end of `taxdb.btd`.
pub const MSG_OFFSET_ERROR: &str = "Error: Offset error at end of taxdb file.";

/// The line a BLAST+ program writes to standard error for an `ERR_POST` of
/// `text` (observed with BLAST+ 2.17.0: `Error: [blastp] Error: Tax database
/// file has wrong magic number.`). `program` is the name of the program, as
/// it appears after `-version`, in lower case.
pub fn err_post_line(program: &str, text: &str) -> String {
    format!("Error: [{program}] {text}")
}

/// The names of one taxid: `SSeqDBTaxInfo` (seqdb.hpp). The names are bytes as
/// they stand in `taxdb.btd`; NCBI does not check the encoding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaxNames {
    /// The taxid asked for.
    pub taxid: i32,
    /// `scientific_name`.
    pub scientific_name: Vec<u8>,
    /// `common_name` (may be empty).
    pub common_name: Vec<u8>,
    /// `blast_name` (may be empty; `-` and `unclassified` are values of the file).
    pub blast_name: Vec<u8>,
    /// `s_kingdom`: the super kingdom, never empty.
    pub super_kingdom: Vec<u8>,
}

/// The two files, open.
#[derive(Debug)]
struct Files {
    index: DbFile,
    data: DbFile,
    /// `m_AllTaxidCount`: equals the number of entries the index has room for, at least 1.
    count: i32,
}

/// `CTaxDBFileInfo`: the opened pair of files, or "missing" (`m_MissingDB`).
#[derive(Debug)]
pub struct TaxDb {
    files: Option<Files>,
    messages: Vec<String>,
}

impl TaxDb {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:128-201
    // ```c++
    // CTaxDBFileInfo::CTaxDBFileInfo()
    //     : m_AllTaxidCount(0), m_IndexPtr(NULL), m_DataPtr(NULL), m_DataFileSize(0), m_MissingDB(false)
    // {
    //     // It is reasonable for this database to not exist.
    //     m_IndexFN =  SeqDB_ResolveDbPath("taxdb.bti");
    //
    //     if (m_IndexFN.size()) {
    //         m_DataFN = m_IndexFN;
    //         m_DataFN[m_DataFN.size()-1] = 'd';
    //     }
    //
    //     if (! (m_IndexFN.size() &&
    //            m_DataFN.size()  &&
    //            CFile(m_IndexFN).Exists() &&
    //            CFile(m_DataFN).Exists())) {
    //         m_MissingDB = true;
    //         return;
    //     }
    //     // Size for header data plus one taxid object.
    //     Uint4 data_start = (4 + 4 + 16);
    //     Uint4 idx_file_len = (Uint4) CFile(m_IndexFN).GetLength();
    //     if (idx_file_len < (data_start + sizeof(CSeqDBTaxId))) {
    //         m_MissingDB = true;
    //         return;
    //     }
    //     ...
    //     if (TAX_DB_MAGIC_NUMBER != SeqDB_GetStdOrd(magic_num_ptr ++)) {
    //         m_MissingDB = true;
    //         m_IndexFileMap.reset();
    //         ERR_POST("Error: Tax database file has wrong magic number.");
    //         return;
    //     }
    //     m_AllTaxidCount = SeqDB_GetStdOrd(magic_num_ptr ++);
    //     // Skip the four reserved fields
    //     magic_num_ptr += 4;
    //     int taxid_array_size = int((idx_file_len - data_start)/sizeof(CSeqDBTaxId));
    //     if (taxid_array_size != m_AllTaxidCount) {
    //         m_MissingDB = true;
    //         m_IndexFileMap.reset();
    //         ERR_POST("SeqDB: Taxid metadata indicates (" << m_AllTaxidCount
    //                    << ") entries but file has room for (" << taxid_array_size
    //                    << ").");
    //         ...
    //         return;
    //     }
    //     m_DataFileMap.reset(new CMemoryFile(m_DataFN));
    //     m_DataPtr = (char *) (m_DataFileMap->GetPtr());
    //     m_DataFileSize = m_DataFileMap->GetSize();
    //     m_IndexPtr = (CSeqDBTaxId*) magic_num_ptr;
    // }
    // ```
    // NCBI builds this object once per process, at the first lookup (`static
    // CTaxDBFileInfo t;` in `GetTaxNames`); a "missing" database stays missing
    // for the rest of the run. `CFile::Exists` and `GetLength` accept a regular
    // file only (`alias::exists`).
    /// Opens the taxonomy database NCBI would use with the search path `search`.
    ///
    /// A database that NCBI treats as missing (no `taxdb.bti` on the path, no
    /// `taxdb.btd` beside it, an index shorter than a header and one entry, a
    /// wrong magic number, a count that does not match the size) opens as a
    /// database that finds no name; the lines NCBI writes for it are in
    /// [`Self::open_messages`]. Err: a read error. A `taxdb.btd` of no bytes
    /// opens (NCBI maps it as a null pointer): lookups that need no record
    /// work, one that reads a record is refused by [`Self::names`].
    pub fn open(search: &SearchPath) -> Result<Self, SeqDbError> {
        let missing = |messages: Vec<String>| {
            Ok(Self {
                files: None,
                messages,
            })
        };
        let Some(index_name) = search.resolve_file("taxdb.bti") else {
            return missing(Vec::new());
        };
        let mut data_name = index_name.clone();
        data_name.pop();
        data_name.push('d');
        let index_path = PathBuf::from(&index_name);
        let data_path = PathBuf::from(&data_name);
        if !super::alias::exists(&index_path) || !super::alias::exists(&data_path) {
            return missing(Vec::new());
        }
        let index = DbFile::open(&index_path)?;
        if index.len() < HEADER_LEN + ENTRY_LEN {
            return missing(Vec::new());
        }
        let mut header = [0u8; HEADER_LEN as usize];
        index.read_exact_at(&mut header, 0)?;
        if be32(&header[0..4]) != MAGIC {
            return missing(vec![MSG_WRONG_MAGIC.to_string()]);
        }
        let count = be32(&header[4..8]) as i32;
        let room = ((index.len() - HEADER_LEN) / ENTRY_LEN) as i32;
        if room != count {
            return missing(vec![format!(
                "SeqDB: Taxid metadata indicates ({count}) entries but file has room for ({room})."
            )]);
        }
        let data = DbFile::open(&data_path)?;
        Ok(Self {
            files: Some(Files { index, data, count }),
            messages: Vec::new(),
        })
    }

    /// The `ERR_POST` texts NCBI writes when it opens this database (wrong
    /// magic number, count that does not match the size); empty for a good
    /// database and for one that is silently missing. NCBI writes them once
    /// per run, at the first lookup, before anything else that lookup writes;
    /// [`err_post_line`] gives the line.
    pub fn open_messages(&self) -> &[String] {
        &self.messages
    }

    /// Whether NCBI treats the database as missing (`IsMissingTaxInfo`): every lookup finds nothing.
    pub fn is_missing(&self) -> bool {
        self.files.is_none()
    }

    /// The `(taxid, offset)` pair `i` of the index (`CSeqDBTaxId`).
    fn entry(files: &Files, i: i32) -> Result<(i32, u32), SeqDbError> {
        if i < 0 || i >= files.count {
            return Err(SeqDbError::Range(format!(
                "taxdb index entry {i} outside 0..{}",
                files.count
            )));
        }
        let mut buf = [0u8; ENTRY_LEN as usize];
        files
            .index
            .read_exact_at(&mut buf, HEADER_LEN + ENTRY_LEN * i as u64)?;
        Ok((be32(&buf[0..4]) as i32, be32(&buf[4..8])))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-300
    // ```c++
    //     static CTaxDBFileInfo t;
    //     if (t.IsMissingTaxInfo()) return false;
    //
    //     Int4 low_index  = 0;
    //     Int4 high_index = t.GetTaxidCount() - 1;
    //     ...
    //     TTaxId low_taxid  = Index[low_index ].GetTaxId();
    //     TTaxId high_taxid = Index[high_index].GetTaxId();
    //
    //     if((tax_id < low_taxid) || (tax_id > high_taxid))
    //         return false;
    //
    //     Int4 new_index =  (low_index+high_index)/2;
    //     Int4 old_index = new_index;
    //
    //     while(1) {
    //         TTaxId curr_taxid = Index[new_index].GetTaxId();
    //
    //         if (tax_id < curr_taxid) {
    //             high_index = new_index;
    //         } else if (tax_id > curr_taxid){
    //             low_index = new_index;
    //         } else { /* Got it ! */
    //             break;
    //         }
    //
    //         new_index = (low_index+high_index)/2;
    //         if (new_index == old_index) {
    //             if (tax_id > curr_taxid) {
    //                 new_index++;
    //             }
    //             break;
    //         }
    //         old_index = new_index;
    //     }
    //
    //     if (tax_id == Index[new_index].GetTaxId()) {
    //         info.taxid = tax_id;
    //
    //         Uint4 begin_data(Index[new_index].GetOffset());
    //         Uint4 end_data(0);
    //
    //         if (new_index == high_index) {
    //             // Last index is special...
    //             end_data = Uint4(t.GetDataFileSize());
    //
    //             if (end_data < begin_data) {
    //                 // Should not happen.
    //                 ERR_POST( "Error: Offset error at end of taxdb file.");
    //                 return false;
    //             }
    //         } else {
    //             end_data = (Index[new_index+1].GetOffset());
    //         }
    //
    //         const char * start_ptr = &Data[begin_data];
    //
    //         CSeqDB_Substring buffer(start_ptr, start_ptr + (end_data - begin_data));
    //         CSeqDB_Substring sci, com, blast, king;
    //         bool rc1, rc2, rc3;
    //
    //         rc1 = SeqDB_SplitString(buffer, sci, '\t');
    //         rc2 = SeqDB_SplitString(buffer, com, '\t');
    //         rc3 = SeqDB_SplitString(buffer, blast, '\t');
    //         king = buffer;
    //
    //         if (rc1 && rc2 && rc3 && buffer.Size()) {
    //             sci   .GetString(info.scientific_name);
    //             ...
    //             return true;
    //         }
    //     }
    //     return false;
    // ```
    // The index is not assumed to be sorted, as in NCBI: the search is NCBI's
    // loop with `Int4` taxids and `Int4` index arithmetic, and "the last
    // entry" is whatever `high_index` ends on (for a sorted index only the
    // last entry; for a damaged one `end_data` may then be the file size).
    // `CSeqDB_Substring::Size()` is `int(m_End-m_Begin)` (seqdbgeneral.hpp:375-378),
    // so a record whose end precedes its beginning (`end_data - begin_data`
    // wraps to 2^32-k as `Uint4`) has a negative size, as has a record of 2 GiB
    // or more: `SeqDB_SplitString` finds nothing in it and the lookup fails.
    // A record that reaches past the end of `taxdb.btd` is read by NCBI from
    // beyond the file: zero bytes up to the end of the last page of its
    // mapping (observed: `offset_next_past_end`, `btd_truncated`), unmapped
    // memory (a crash) after that. The page is 4096 bytes on the oracle
    // platform; beyond it the lookup is refused.
    /// `CSeqDBTaxInfo::GetTaxNames` for `tax_id`: `Ok(Some(names))` when NCBI
    /// returns true, `Ok(None)` when it returns false (the taxid is not in the
    /// database, the record is not of four fields, the database is missing).
    ///
    /// The `ERR_POST` texts NCBI writes during this lookup are appended to
    /// `messages` (only [`MSG_OFFSET_ERROR`]; those of opening are in
    /// [`Self::open_messages`]). A record that ends after the file is padded with
    /// zero bytes up to the end of a 4096-byte page, as NCBI's mapping is. Err: a
    /// record that reaches past that page (NCBI reads unmapped memory).
    pub fn names(
        &self,
        tax_id: i32,
        messages: &mut Vec<String>,
    ) -> Result<Option<TaxNames>, SeqDbError> {
        let Some(files) = &self.files else {
            return Ok(None);
        };
        let mut low_index: i32 = 0;
        let mut high_index: i32 = files.count - 1;
        let low_taxid = Self::entry(files, low_index)?.0;
        let high_taxid = Self::entry(files, high_index)?.0;
        if tax_id < low_taxid || tax_id > high_taxid {
            return Ok(None);
        }
        let mut new_index = (low_index + high_index) / 2;
        let mut old_index = new_index;
        loop {
            let curr_taxid = Self::entry(files, new_index)?.0;
            if tax_id < curr_taxid {
                high_index = new_index;
            } else if tax_id > curr_taxid {
                low_index = new_index;
            } else {
                break;
            }
            new_index = (low_index + high_index) / 2;
            if new_index == old_index {
                if tax_id > curr_taxid {
                    new_index += 1;
                }
                break;
            }
            old_index = new_index;
        }
        let (found, begin_data) = Self::entry(files, new_index)?;
        if found != tax_id {
            return Ok(None);
        }
        let end_data: u32 = if new_index == high_index {
            // Last index is special...
            let end = files.data.len() as u32;
            if end < begin_data {
                messages.push(MSG_OFFSET_ERROR.to_string());
                return Ok(None);
            }
            end
        } else {
            Self::entry(files, new_index + 1)?.1
        };
        // `int(m_End - m_Begin)` of a pointer difference that is the wrapped
        // `Uint4` difference.
        let size = end_data.wrapping_sub(begin_data) as i32;
        if size <= 0 {
            return Ok(None);
        }
        let begin = u64::from(begin_data);
        let end = begin + size as u64;
        let len = files.data.len();
        // NCBI maps the file; the last page is padded with zero bytes. Past that page
        // it reads unmapped memory (a crash on the oracle platform).
        let mapped = len.div_ceil(MAP_PAGE) * MAP_PAGE;
        if end > mapped {
            return Err(SeqDbError::Unsupported {
                path: files.data.path().to_path_buf(),
                message: format!(
                    "the taxonomy record of taxid {tax_id} (bytes {begin}..{end}) reaches past \
                     the mapped end of the file ({len} bytes); NCBI BLAST+ reads unmapped memory"
                ),
            });
        }
        let mut record = files.data.read_range(begin.min(len), end.min(len))?;
        record.resize(size as usize, 0);
        Ok(split_record(tax_id, &record))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:113-127
// ```c++
// bool SeqDB_SplitString(CSeqDB_Substring & buffer,
//                        CSeqDB_Substring & front,
//                        char               delim)
// {
//     for(int i = 0; i < buffer.Size(); i++) {
//         if (buffer[i] == delim) {
//             front = buffer;
//
//             buffer.EraseFront(i + 1);
//             front.Resize(i);
//
//             return true;
//         }
//     }
//     return false;
// }
// ```
// Three splits at the first tab of what is left; a split that finds no tab
// leaves the buffer as it is and returns false. The kingdom is what remains,
// and must not be empty (`EraseFront` of the whole buffer empties it).
/// The split of one record into four fields (seqdbtax.cpp:279-298).
fn split_record(taxid: i32, record: &[u8]) -> Option<TaxNames> {
    let mut rest = record;
    let mut fields: [&[u8]; 3] = [&[]; 3];
    let mut all = true;
    for field in &mut fields {
        match rest.iter().position(|&b| b == b'\t') {
            Some(i) => {
                *field = &rest[..i];
                rest = &rest[i + 1..];
            }
            None => all = false,
        }
    }
    if all && !rest.is_empty() {
        Some(TaxNames {
            taxid,
            scientific_name: fields[0].to_vec(),
            common_name: fields[1].to_vec(),
            blast_name: fields[2].to_vec(),
            super_kingdom: rest.to_vec(),
        })
    } else {
        None
    }
}

fn be32(bytes: &[u8]) -> u32 {
    u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

/// The taxonomy database of one run, opened at the first name asked for and
/// not before (NCBI: the `static CTaxDBFileInfo` of `GetTaxNames`). A run that
/// never asks for names touches no file.
#[derive(Debug)]
pub struct LazyTaxDb {
    search: SearchPath,
    db: OnceLock<Result<TaxDb, (PathBuf, String)>>,
}

impl LazyTaxDb {
    /// A handle that will open the database through `search`.
    pub fn new(search: SearchPath) -> Self {
        Self {
            search,
            db: OnceLock::new(),
        }
    }

    /// `GetTaxNames` through the lazily opened database. The `ERR_POST` texts
    /// of opening, which NCBI writes once at the first lookup, are appended to
    /// `messages` by the call that opens the database, before those of the
    /// lookup itself. An open that was refused is the same error on every call.
    pub fn names(
        &self,
        tax_id: i32,
        messages: &mut Vec<String>,
    ) -> Result<Option<TaxNames>, SeqDbError> {
        let mut opened = false;
        let db = self.db.get_or_init(|| {
            opened = true;
            TaxDb::open(&self.search).map_err(|error| match error {
                SeqDbError::Unsupported { path, message } => (path, message),
                other => (PathBuf::from("taxdb.bti"), other.to_string()),
            })
        });
        match db {
            Ok(db) => {
                if opened {
                    messages.extend(db.open_messages().iter().cloned());
                }
                db.names(tax_id, messages)
            }
            Err((path, message)) => Err(SeqDbError::Unsupported {
                path: path.clone(),
                message: message.clone(),
            }),
        }
    }

    /// Whether the database has been opened (a name has been asked for).
    pub fn is_open(&self) -> bool {
        self.db.get().is_some()
    }
}

#[cfg(test)]
mod tests {
    //! Tests against taxdb files written by `make_taxdb.py` (the fixtures under
    //! `tests/seqdb/taxdb/<variant>/`, with what NCBI BLAST+ 2.17.0 returned
    //! for them in `tests/seqdb/taxdb/expected.tsv`) and against files built
    //! here for the cases NCBI cannot be asked about.

    use super::*;
    use std::path::Path;

    fn fixture_dir() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb")
    }

    fn temp_dir(tag: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-taxdb-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// A search path whose only useful directory is `dir` (`BLASTDB`), run from an empty current directory.
    fn path_of(dir: &Path) -> SearchPath {
        let cwd = temp_dir("cwd");
        SearchPath::new(&cwd, Some(dir.as_os_str()), None).unwrap()
    }

    fn open_fixture(variant: &str) -> Result<TaxDb, SeqDbError> {
        TaxDb::open(&path_of(&fixture_dir().join(variant)))
    }

    fn be(v: u32) -> [u8; 4] {
        v.to_be_bytes()
    }

    /// Builds a taxdb: `rows` are (taxid, record), laid out back to back in the order given.
    fn build(rows: &[(i32, &[u8])]) -> (Vec<u8>, Vec<u8>) {
        let mut index = Vec::new();
        index.extend(be(MAGIC));
        index.extend(be(rows.len() as u32));
        index.extend([0u8; 16]);
        let mut data = Vec::new();
        for (taxid, record) in rows {
            index.extend(be(*taxid as u32));
            index.extend(be(data.len() as u32));
            data.extend(*record);
        }
        (index, data)
    }

    fn write_db(index: &[u8], data: &[u8]) -> PathBuf {
        let dir = temp_dir("db");
        std::fs::write(dir.join("taxdb.bti"), index).unwrap();
        std::fs::write(dir.join("taxdb.btd"), data).unwrap();
        dir
    }

    fn names_of(db: &TaxDb, taxid: i32) -> (Option<TaxNames>, Vec<String>) {
        let mut messages = Vec::new();
        let names = db.names(taxid, &mut messages).unwrap();
        (names, messages)
    }

    fn names_of_lazy(lazy: &LazyTaxDb, taxid: i32) -> (Option<TaxNames>, Vec<String>) {
        let mut messages = Vec::new();
        let names = lazy.names(taxid, &mut messages).unwrap();
        (names, messages)
    }

    fn four(names: &TaxNames) -> [&[u8]; 4] {
        [
            &names.scientific_name,
            &names.common_name,
            &names.blast_name,
            &names.super_kingdom,
        ]
    }

    #[test]
    fn good_database_finds_every_entry_and_nothing_else() {
        let db = open_fixture("good").unwrap();
        assert!(!db.is_missing());
        assert!(db.open_messages().is_empty());
        let expected: [(i32, [&[u8]; 4]); 7] = [
            (
                562,
                [b"Escherichia coli", b"", b"enterobacteria", b"Bacteria"],
            ),
            (7777, [b"Quote \"x\", y; z", b"a,b", b"c\"d", b"K;L"]),
            (9606, [b"Homo sapiens", b"human", b"primates", b"Eukaryota"]),
            (
                10090,
                [b"Mus musculus", b"house mouse", b"rodents", b"Eukaryota"],
            ),
            (
                12908,
                [
                    b"unclassified sequences",
                    b"",
                    b"unclassified",
                    b"unclassified",
                ],
            ),
            (32630, [b"synthetic construct", b"", b"-", b"-"]),
            // The last entry: the record ends at the end of taxdb.btd.
            (
                99990,
                [
                    b"Last entry sp.",
                    b"last common",
                    b"lastblast",
                    b"Eukaryota",
                ],
            ),
        ];
        for (taxid, fields) in expected {
            let (names, messages) = names_of(&db, taxid);
            let names = names.unwrap_or_else(|| panic!("taxid {taxid}"));
            assert_eq!(names.taxid, taxid);
            assert_eq!(four(&names), fields, "taxid {taxid}");
            assert!(messages.is_empty());
        }
        // Between two entries, below the first, above the last, zero, negative.
        for taxid in [
            777,
            8000,
            50000,
            3,
            561,
            99991,
            100000,
            0,
            -1,
            i32::MIN,
            i32::MAX,
        ] {
            assert_eq!(names_of(&db, taxid).0, None, "taxid {taxid}");
        }
    }

    #[test]
    fn every_position_of_every_size_is_found() {
        // The search loop for index sizes 1..=40, every entry found, every gap not found.
        for n in 1..=40i32 {
            let rows: Vec<(i32, Vec<u8>)> = (0..n)
                .map(|i| (10 * i + 5, format!("s{i}\tc\tb\tk{i}").into_bytes()))
                .collect();
            let borrowed: Vec<(i32, &[u8])> =
                rows.iter().map(|(t, r)| (*t, r.as_slice())).collect();
            let (index, data) = build(&borrowed);
            let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
            for i in 0..n {
                let names = names_of(&db, 10 * i + 5)
                    .0
                    .unwrap_or_else(|| panic!("n={n} i={i}"));
                assert_eq!(names.scientific_name, format!("s{i}").into_bytes());
                assert_eq!(names.super_kingdom, format!("k{i}").into_bytes());
                assert_eq!(names_of(&db, 10 * i + 6).0, None);
                assert_eq!(names_of(&db, 10 * i + 4).0, None);
            }
        }
    }

    #[test]
    fn missing_and_unusable_databases_find_nothing() {
        // (variant, messages written when it is opened)
        let cases: [(&str, &[&str]); 11] = [
            ("no_bti", &[]),
            ("no_btd", &[]),
            ("index_short", &[]),
            ("index_header_only", &[]),
            ("bad_magic", &[MSG_WRONG_MAGIC]),
            ("bad_magic_le", &[MSG_WRONG_MAGIC]),
            (
                "count_larger",
                &["SeqDB: Taxid metadata indicates (8) entries but file has room for (7)."],
            ),
            (
                "count_smaller",
                &["SeqDB: Taxid metadata indicates (6) entries but file has room for (7)."],
            ),
            (
                "count_zero",
                &["SeqDB: Taxid metadata indicates (0) entries but file has room for (7)."],
            ),
            (
                "index_trailing_8",
                &["SeqDB: Taxid metadata indicates (7) entries but file has room for (8)."],
            ),
            ("nonexistent", &[]),
        ];
        for (variant, messages) in cases {
            let db = open_fixture(variant).unwrap();
            assert!(db.is_missing(), "{variant}");
            assert_eq!(db.open_messages(), messages, "{variant}");
            assert_eq!(names_of(&db, 9606), (None, vec![]), "{variant}");
        }
    }

    #[test]
    fn index_sizes_at_the_edges() {
        // Exactly header + one entry is valid; seven stray bytes are ignored (integer division).
        let db = open_fixture("index_one_entry").unwrap();
        assert!(!db.is_missing());
        assert!(names_of(&db, 9606).0.is_some());
        assert_eq!(names_of(&db, 9607).0, None);
        let db = open_fixture("index_trailing").unwrap();
        assert!(!db.is_missing());
        assert!(names_of(&db, 99990).0.is_some());
    }

    #[test]
    fn empty_data_file_is_refused_only_where_a_record_is_read() {
        let db = open_fixture("btd_empty").unwrap();
        assert!(!db.is_missing());
        let mut messages = Vec::new();
        // NCBI: segmentation fault.
        assert!(matches!(
            db.names(9606, &mut messages),
            Err(SeqDbError::Unsupported { .. })
        ));
        // NCBI: not found; and the offset error of the last entry.
        assert_eq!(names_of(&db, 777), (None, vec![]));
        assert_eq!(names_of(&db, 100000), (None, vec![]));
        assert_eq!(
            names_of(&db, 99990),
            (None, vec![MSG_OFFSET_ERROR.to_string()])
        );
    }

    #[test]
    fn the_index_path_decides_and_the_data_file_is_its_sibling() {
        // The first directory with a taxdb.bti wins, even when its taxdb.btd is missing.
        let good = fixture_dir().join("good");
        let nobtd = fixture_dir().join("no_btd");
        let joined = std::env::join_paths([&nobtd, &good]).unwrap();
        let cwd = temp_dir("cwd");
        let search = SearchPath::new(&cwd, Some(joined.as_os_str()), None).unwrap();
        assert!(TaxDb::open(&search).unwrap().is_missing());
        let joined = std::env::join_paths([&fixture_dir().join("no_bti"), &good]).unwrap();
        let search = SearchPath::new(&cwd, Some(joined.as_os_str()), None).unwrap();
        assert!(!TaxDb::open(&search).unwrap().is_missing());
        // The current directory comes before BLASTDB.
        let search = SearchPath::new(&good, Some(nobtd.as_os_str()), None).unwrap();
        assert!(!TaxDb::open(&search).unwrap().is_missing());
        // The registry value comes after BLASTDB.
        let search = SearchPath::new(&cwd, Some(nobtd.as_os_str()), good.to_str()).unwrap();
        assert!(TaxDb::open(&search).unwrap().is_missing());
        let search = SearchPath::new(&cwd, None, good.to_str()).unwrap();
        assert!(!TaxDb::open(&search).unwrap().is_missing());
        // No directory at all.
        let search = SearchPath::new(&cwd, None, None).unwrap();
        assert!(TaxDb::open(&search).unwrap().is_missing());
    }

    #[test]
    fn records_with_other_than_four_fields() {
        let rows: [(i32, &[u8]); 9] = [
            (1, b"a\tb\tc\td"),
            (2, b"a\tb\tc"),       // two tabs: no kingdom
            (3, b"a\tb\tc\t"),     // empty kingdom
            (4, b"a"),             // no tab
            (5, b"a\tb\tc\td\te"), // kingdom takes the rest
            (6, b"\t\t\tx"),       // all other fields empty
            (7, b"a\tb"),
            (8, b"a\tb\tc\td\n"), // a terminator belongs to the kingdom
            (9, b"a\tb\tc\td\0"),
        ];
        let (index, data) = build(&rows);
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        assert_eq!(four(&names_of(&db, 1).0.unwrap()), [b"a", b"b", b"c", b"d"]);
        assert_eq!(names_of(&db, 2).0, None);
        assert_eq!(names_of(&db, 3).0, None);
        assert_eq!(names_of(&db, 4).0, None);
        assert_eq!(
            four(&names_of(&db, 5).0.unwrap()),
            [b"a".as_slice(), b"b", b"c", b"d\te"]
        );
        assert_eq!(
            four(&names_of(&db, 6).0.unwrap()),
            [b"".as_slice(), b"", b"", b"x"]
        );
        assert_eq!(names_of(&db, 7).0, None);
        assert_eq!(four(&names_of(&db, 8).0.unwrap())[3], b"d\n");
        // The last entry: the record ends at the end of the file.
        assert_eq!(four(&names_of(&db, 9).0.unwrap())[3], b"d\0");
    }

    #[test]
    fn last_entry_record_shapes() {
        for (record, found) in [
            (&b"a\tb\tc\td"[..], true),
            (b"a\tb\tc\t", false),
            (b"a\tb\tc", false),
            (b"a", false),
            (b"", false),
        ] {
            let (index, data) = build(&[(1, b"x\ty\tz\tw"), (2, record)]);
            let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
            assert_eq!(names_of(&db, 2).0.is_some(), found, "{record:?}");
            assert!(names_of(&db, 1).0.is_some());
        }
    }

    #[test]
    fn last_offset_past_the_end_is_an_error_message_and_not_found() {
        let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th")]);
        let at = 24 + 8 + 4; // offset of the second entry
        index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        let (names, messages) = names_of(&db, 2);
        assert_eq!(names, None);
        assert_eq!(messages, [MSG_OFFSET_ERROR]);
        // The message comes at every lookup.
        assert_eq!(names_of(&db, 2).1, [MSG_OFFSET_ERROR]);
        // An offset equal to the size is an empty record: not found, no message.
        index[at..at + 4].copy_from_slice(&be(data.len() as u32));
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        assert_eq!(names_of(&db, 2), (None, vec![]));
    }

    #[test]
    fn record_ending_before_it_begins_is_not_found() {
        // Entry 1 starts after entry 2: the size of 1 is negative.
        let (mut index, data) =
            build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
        let off1 = 24 + 4;
        let off2 = 24 + 8 + 4;
        index[off1..off1 + 4].copy_from_slice(&be(8));
        index[off2..off2 + 4].copy_from_slice(&be(0));
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        assert_eq!(names_of(&db, 1), (None, vec![]));
        // Entry 2 is now [0, 16): both records, the kingdom takes the rest.
        let names = names_of(&db, 2).0.unwrap();
        assert_eq!(names.scientific_name, b"a");
        assert_eq!(
            names.super_kingdom,
            b"d".to_vec()
                .into_iter()
                .chain(*b"e\tf\tg\th")
                .collect::<Vec<u8>>()
        );
    }

    #[test]
    fn record_reaching_past_the_end_reads_page_padding_and_then_is_refused() {
        let (mut index, data) =
            build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
        let off3 = 24 + 16 + 4;
        // Entry 2 ends 10 bytes after the file: zero padding belongs to its kingdom.
        index[off3..off3 + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        let names = names_of(&db, 2).0.unwrap();
        let mut kingdom = b"h".to_vec();
        kingdom.extend(*b"i\tj\tk\tl");
        kingdom.extend([0u8; 10]);
        assert_eq!(names.super_kingdom, kingdom);
        // Entry 3 (the last) begins after the end: NCBI's message.
        assert_eq!(names_of(&db, 3), (None, vec![MSG_OFFSET_ERROR.to_string()]));
        // Past the page: refused. Entry 1 is fine.
        index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 1));
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        let mut messages = Vec::new();
        assert!(matches!(
            db.names(2, &mut messages),
            Err(SeqDbError::Unsupported { .. })
        ));
        assert!(names_of(&db, 1).0.is_some());
    }

    #[test]
    fn unsorted_index_follows_the_search_loop() {
        // A port of the loop gives these results for an unsorted index; the same loop in NCBI
        // finds an entry only where its probes happen to land.
        let rows: [(i32, &[u8]); 5] = [
            (10, b"a\tb\tc\td10"),
            (50, b"a\tb\tc\td50"),
            (20, b"a\tb\tc\td20"),
            (40, b"a\tb\tc\td40"),
            (60, b"a\tb\tc\td60"),
        ];
        let (index, data) = build(&rows);
        let db = TaxDb::open(&path_of(&write_db(&index, &data))).unwrap();
        let found: Vec<i32> = [10, 20, 40, 50, 60]
            .into_iter()
            .filter(|t| names_of(&db, *t).0.is_some())
            .collect();
        // 50: new=2 (20) -> low=2, new=3 (40) -> low=3, new=3 again: 50 > 40 so new=4 (60), which is not 50.
        assert_eq!(found, [10, 20, 40, 60]);
        // High-bit taxids are negative Int4 values, sorted after 99990 in the file: out of range.
        let db = open_fixture("taxid_high_bit").unwrap();
        assert_eq!(names_of(&db, 9606).0, None);
        assert_eq!(names_of(&db, i32::MIN + 5).0, None);
    }

    #[test]
    fn names_are_bytes() {
        let db = open_fixture("high_bytes_utf8").unwrap();
        assert_eq!(
            names_of(&db, 9606).0.unwrap().scientific_name,
            "Caf\u{e9} sapiens".as_bytes()
        );
        let db = open_fixture("high_bytes_latin1").unwrap();
        let names = names_of(&db, 10090).0.unwrap();
        assert_eq!(names.common_name, b"mus\xe9");
        assert_eq!(names.super_kingdom, b"Euk\xffaryota");
    }

    #[test]
    fn lazy_handle_opens_once_and_hands_over_the_open_messages_once() {
        let dir = fixture_dir().join("bad_magic");
        let lazy = LazyTaxDb::new(path_of(&dir));
        assert!(!lazy.is_open());
        let mut messages = Vec::new();
        assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
        assert!(lazy.is_open());
        assert_eq!(messages, [MSG_WRONG_MAGIC]);
        assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
        assert_eq!(messages, [MSG_WRONG_MAGIC]);
        // A missing database: opened at the first name, silent, finds nothing.
        let lazy = LazyTaxDb::new(path_of(&fixture_dir().join("nonexistent")));
        assert!(!lazy.is_open());
        for _ in 0..2 {
            assert_eq!(names_of_lazy(&lazy, 9606), (None, vec![]));
        }
        assert!(lazy.is_open());
    }

    #[test]
    fn err_post_line_is_what_blast_prints() {
        assert_eq!(
            err_post_line("blastp", MSG_WRONG_MAGIC),
            "Error: [blastp] Error: Tax database file has wrong magic number."
        );
    }

    /// `make_expected.py`'s escaping: bytes outside 0x20..0x7e and the backslash as `\xNN`.
    fn escape(bytes: &[u8]) -> String {
        bytes
            .iter()
            .map(|&b| {
                if (0x20..0x7f).contains(&b) && b != b'\\' {
                    (b as char).to_string()
                } else {
                    format!("\\x{b:02x}")
                }
            })
            .collect()
    }

    /// Every (variant, taxid) of `expected.tsv`: what `blastdbcmd -tax_info` of NCBI BLAST+ 2.17.0
    /// printed (`%T|%S|%L|%B|%K`, escaped), `N/A` for a taxid it did not find, `CRASH` for a
    /// signal;
    /// then the lines it wrote to standard error.
    #[test]
    fn matches_blastdbcmd_tax_info_on_every_variant() {
        let text = std::fs::read_to_string(fixture_dir().join("expected.tsv")).unwrap();
        let mut checked = 0;
        for line in text.lines().filter(|l| !l.is_empty()) {
            let cols: Vec<&str> = line.splitn(4, '\t').collect();
            let (variant, taxid, expect) = (cols[0], cols[1].parse::<i32>().unwrap(), cols[2]);
            let ncbi_messages: Vec<&str> =
                cols[3].split('\u{1f}').filter(|m| !m.is_empty()).collect();
            let what = format!("{variant} taxid {taxid}");
            let opened = open_fixture(variant);
            let mut messages: Vec<String> = Vec::new();
            let result = match &opened {
                Ok(db) => {
                    messages.extend(db.open_messages().iter().cloned());
                    db.names(taxid, &mut messages)
                }
                Err(error) => Err(SeqDbError::Range(error.to_string())),
            };
            match (expect, &result) {
                ("CRASH", Err(_)) => continue,
                ("N/A", Ok(None)) => {}
                (expect, Ok(Some(names))) => {
                    let mut raw = format!("{}|", names.taxid).into_bytes();
                    for (i, field) in four(names).iter().enumerate() {
                        if i > 0 {
                            raw.push(b'|');
                        }
                        raw.extend(*field);
                    }
                    assert_eq!(escape(&raw), expect, "{what}");
                }
                (expect, other) => panic!("{what}: NCBI {expect:?}, got {other:?}"),
            }
            assert_eq!(messages, ncbi_messages, "{what}: messages");
            checked += 1;
        }
        assert!(checked > 300, "{checked}");
    }
}
