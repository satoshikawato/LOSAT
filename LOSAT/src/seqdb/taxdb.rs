//! Reader of NCBI's taxonomy name database (`taxdb.bti` and `taxdb.btd`).
//!
//! NCBI BLAST+ writes the scientific, common and BLAST names of a taxid, and
//! the super kingdom, from this pair of files (`CSeqDBTaxInfo::GetTaxNames`,
//! c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp): the `sciname` of the
//! XML2 and JSON reports (`-outfmt 13` to `16`), the tabular fields `sscinames`,
//! `scomnames`, `sblastnames`, `sskingdoms` (and their singular forms) and the
//! organism report (`-outfmt 18`). The lookups also write NCBI's messages on
//! standard error when a query batch of outfmt 0, 5, 6, 7 or 17 prefetches the
//! subjects it reports (`CSeqDBVol::x_GetTaxonomy`; the texts of opening and
//! [`MSG_OFFSET_ERROR`], `algorithm::common::database::prefetch_messages`).
//! This module is the reader only; the report code decides when a name is
//! looked up and what it writes when the lookup finds nothing.
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
//! A record that reaches past the end of `taxdb.btd` is read as NCBI reads its
//! mapping: the zero bytes of the mapping's last 4096-byte page, and after that
//! page memory that does not belong to the file, whose bytes, or whether NCBI
//! faults reading them, depend on NCBI's memory map. Only an empty `taxdb.btd`
//! faults for certain (NCBI maps it as a null pointer): a lookup that reads a
//! record of it is refused with [`SeqDbError::Unsupported`]. A record past the
//! mapping of a file that is not empty gives NCBI's result where its bytes lie
//! inside the mapping: the names before the record's third tab when that tab is
//! inside, the kingdom never ([`TaxNames::super_kingdom`]); where no third tab
//! is inside, whether the lookup finds names at all is not defined
//! ([`TaxLookup::PastMapping`]). The callers refuse only what they write of it
//! (PD-LOSAT-NCBI-DEFECTS, rule 3 and approved exception 14). Lookups
//! that read no record (a taxid outside the index, a record that ends before it
//! begins, the offset error of the last entry) give NCBI's result.
//!
//! The reader opens nothing until a name is asked for ([`LazyTaxDb`]); it
//! reads the index by position (about 22 reads of 8 bytes per lookup), so it
//! holds no table in memory. The readers of the stored titles keep what they
//! looked up as NCBI's volume does, in 256 slots ([`TaxCache`]).

use std::path::PathBuf;
use std::sync::OnceLock;

use super::file::DbFile;
use super::{SearchPath, SeqDbError};

/// Header of `taxdb.bti`: magic, count, four reserved words (seqdbtax.cpp,
/// `data_start`).
const HEADER_LEN: u64 = 4 + 4 + 16;
/// One (taxid, offset) pair of `taxdb.bti` (`CSeqDBTaxId`).
const ENTRY_LEN: u64 = 8;
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
    /// `s_kingdom`: the super kingdom, never empty. `Err` where the record reaches past NCBI's
    /// mapping of `taxdb.btd`: NCBI copies the kingdom from memory past it, whose bytes (or
    /// whether NCBI faults there) depend on its memory map; a caller that writes the kingdom
    /// refuses it ([`PastMapping::error`]).
    pub super_kingdom: Result<Vec<u8>, PastMapping>,
}

/// What `CSeqDBTaxInfo::GetTaxNames` gives for a taxid ([`TaxDb::lookup`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaxLookup {
    /// False before the search finds the taxid (a missing database, a taxid outside the
    /// index or not in it): the caller's `info.taxid` is not set.
    NotInIndex,
    /// False after the search found it (`info.taxid = tax_id`; seqdbtax.cpp:258-259): the offset
    /// error of the last entry, a record of no positive size or not of four fields.
    Unnamed,
    /// True: the names.
    Named(TaxNames),
    /// The search found the taxid (`info.taxid = tax_id`) and its record reaches past NCBI's
    /// mapping of `taxdb.btd` before its third tab: NCBI looks for the tabs in memory past the
    /// mapping, so whether it finds names, and which, depends on its memory map (or it faults).
    /// Nothing of the lookup but its texts can be written.
    PastMapping(PastMapping),
}

impl TaxLookup {
    /// The names, `None` where NCBI's lookup returns false. Err: a lookup whose result depends on
    /// memory past NCBI's mapping of `taxdb.btd` ([`TaxLookup::PastMapping`]), for a caller that
    /// writes something of it.
    pub fn into_names(self) -> Result<Option<TaxNames>, SeqDbError> {
        match self {
            TaxLookup::Named(names) => Ok(Some(names)),
            TaxLookup::NotInIndex | TaxLookup::Unnamed => Ok(None),
            TaxLookup::PastMapping(past) => Err(past.error("names")),
        }
    }
}

/// A record of `taxdb.btd` that reaches past NCBI's mapping of the file (the end of its last
/// 4096-byte page), which NCBI reads there: memory that does not belong to the file, whose bytes,
/// or whether NCBI faults reading them, depend on its memory map ([`TaxLookup::PastMapping`],
/// [`TaxNames::super_kingdom`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PastMapping {
    /// `taxdb.btd`.
    pub path: PathBuf,
    /// The taxid looked up.
    pub taxid: i32,
    /// The record, `begin..end` in bytes of the file.
    pub begin: u64,
    pub end: u64,
    /// The size of the file.
    pub len: u64,
}

impl PastMapping {
    /// LOSAT's refusal of a caller that writes `what` ("names", "super kingdom") of this lookup:
    /// those bytes come from past NCBI's mapping of `taxdb.btd`.
    pub fn error(&self, what: &str) -> SeqDbError {
        SeqDbError::Unsupported {
            path: self.path.clone(),
            message: format!(
                "the {what} of taxid {} in the taxonomy database, read from its record (bytes \
                 {}..{}) past the end of the file ({} bytes) and NCBI BLAST+'s mapping of it: \
                 NCBI reads memory past its mapping there, whose bytes, or whether NCBI faults, \
                 depend on its memory map",
                self.taxid, self.begin, self.end, self.len
            ),
        }
    }
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
    #[cfg(test)]
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
    // mapping (observed: `offset_next_past_end`, `btd_truncated`), then memory
    // past the mapping, whose bytes, or whether NCBI faults reading them, depend
    // on its memory map (review of line RG item T5b: a record of 266-byte
    // `taxdb.btd` that ends at 5,000, 9,000 or 20,000 bytes completes, one that
    // ends at 1 MiB faults). The page is 4096 bytes on the oracle platform.
    // `SeqDB_SplitString` reads only up to the tab it finds and `GetString`
    // copies a field whole, so a third tab inside the mapping gives the first
    // three names and a kingdom that reaches past it. An empty `taxdb.btd` is
    // mapped as a null pointer (`m_DataPtr`): every record read faults.
    /// `CSeqDBTaxInfo::GetTaxNames` for `tax_id`: `Ok(Some(names))` when NCBI
    /// returns true, `Ok(None)` when it returns false (the taxid is not in the
    /// database, the record is not of four fields, the database is missing).
    ///
    /// The `ERR_POST` texts NCBI writes during this lookup are appended to
    /// `messages` (only [`MSG_OFFSET_ERROR`]; those of opening are in
    /// [`Self::open_messages`]). A record that ends after the file is padded with
    /// zero bytes up to the end of a 4096-byte page, as NCBI's mapping is. Err: a
    /// record of an empty `taxdb.btd` (NCBI faults), and a record past the
    /// mapping before its third tab ([`TaxLookup::into_names`]).
    pub fn names(
        &self,
        tax_id: i32,
        messages: &mut Vec<String>,
    ) -> Result<Option<TaxNames>, SeqDbError> {
        self.lookup(tax_id, messages)?.into_names()
    }

    /// [`Self::names`], telling apart the two ways in which NCBI's lookup returns false: before
    /// its search finds the taxid in the index ([`TaxLookup::NotInIndex`]: `info.taxid` keeps
    /// its value) and after ([`TaxLookup::Unnamed`]: `info.taxid = tax_id`, the names stay empty),
    /// which blastdbcmd's `-tax_info` writes apart (`N/A` or empty names), and giving a record
    /// past NCBI's mapping of `taxdb.btd` as it is ([`TaxLookup::PastMapping`],
    /// [`TaxNames::super_kingdom`]): a caller that writes nothing of the names (the texts of
    /// `GetBioseq`'s lookups) goes on as NCBI does where it completes. Err: a read error, or a
    /// record of an empty `taxdb.btd` (NCBI faults for certain).
    pub fn lookup(&self, tax_id: i32, messages: &mut Vec<String>) -> Result<TaxLookup, SeqDbError> {
        let Some(files) = &self.files else {
            return Ok(TaxLookup::NotInIndex);
        };
        let mut low_index: i32 = 0;
        let mut high_index: i32 = files.count - 1;
        let low_taxid = Self::entry(files, low_index)?.0;
        let high_taxid = Self::entry(files, high_index)?.0;
        if tax_id < low_taxid || tax_id > high_taxid {
            return Ok(TaxLookup::NotInIndex);
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
            return Ok(TaxLookup::NotInIndex);
        }
        let end_data: u32 = if new_index == high_index {
            // Last index is special...
            let end = files.data.len() as u32;
            if end < begin_data {
                messages.push(MSG_OFFSET_ERROR.to_string());
                return Ok(TaxLookup::Unnamed);
            }
            end
        } else {
            Self::entry(files, new_index + 1)?.1
        };
        // `int(m_End - m_Begin)` of a pointer difference that is the wrapped
        // `Uint4` difference.
        let size = end_data.wrapping_sub(begin_data) as i32;
        if size <= 0 {
            return Ok(TaxLookup::Unnamed);
        }
        let begin = u64::from(begin_data);
        let end = begin + size as u64;
        let len = files.data.len();
        if len == 0 {
            // NCBI maps an empty file as a null pointer: it faults reading the record.
            return Err(SeqDbError::Unsupported {
                path: files.data.path().to_path_buf(),
                message: format!(
                    "the taxonomy record of taxid {tax_id} (bytes {begin}..{end}) reaches past \
                     the mapped end of the file (0 bytes); NCBI BLAST+ reads unmapped memory (it \
                     maps the empty file as a null pointer)"
                ),
            });
        }
        // NCBI maps the file; the last page is padded with zero bytes.
        let mapped = super::file::mapped_length(len);
        let mut record = files.data.read_range(begin.min(len), end.min(len))?;
        record.resize((end.min(mapped).saturating_sub(begin)) as usize, 0);
        if end <= mapped {
            return Ok(split_record(tax_id, &record).map_or(TaxLookup::Unnamed, TaxLookup::Named));
        }
        // The part of the record inside the mapping (none where it begins past it).
        let past = PastMapping {
            path: files.data.path().to_path_buf(),
            taxid: tax_id,
            begin,
            end,
            len,
        };
        Ok(match split_names(&record) {
            Some((fields, _)) => TaxLookup::Named(TaxNames {
                taxid: tax_id,
                scientific_name: fields[0].to_vec(),
                common_name: fields[1].to_vec(),
                blast_name: fields[2].to_vec(),
                super_kingdom: Err(past),
            }),
            None => TaxLookup::PastMapping(past),
        })
    }

    /// The number of entries of the index (`m_AllTaxidCount`); 0 for a database NCBI treats as
    /// missing.
    pub fn entry_count(&self) -> u64 {
        self.files.as_ref().map_or(0, |files| files.count as u64)
    }

    /// The taxid whose lookup writes [`MSG_OFFSET_ERROR`] ([`Self::names`]): the last entry's,
    /// when its offset is past the end of `taxdb.btd` and the search for its taxid ends on it (in
    /// an index that is not sorted the search may end elsewhere). `None` when no lookup writes the
    /// message, and for a database NCBI treats as missing. No other lookup writes a text (those of
    /// opening are [`Self::open_messages`]), and this one writes it at every lookup of the taxid.
    /// Err: a read error, or the lookup of that taxid faults (an empty `taxdb.btd`).
    pub fn offset_error_taxid(&self) -> Result<Option<i32>, SeqDbError> {
        let Some(files) = &self.files else {
            return Ok(None);
        };
        let (taxid, offset) = Self::entry(files, files.count - 1)?;
        // `names`: `end_data < begin_data` with `end_data` the file size as `Uint4`.
        if offset <= files.data.len() as u32 {
            return Ok(None);
        }
        let mut messages = Vec::new();
        self.lookup(taxid, &mut messages)?;
        Ok(messages
            .iter()
            .any(|message| message == MSG_OFFSET_ERROR)
            .then_some(taxid))
    }

    // No NCBI counterpart: LOSAT's check, from the index alone, of whether `lookup` (the port of
    // `GetTaxNames` above) reads past NCBI's mapping of `taxdb.btd`; it changes no value NCBI
    // computes.
    /// Whether [`Self::lookup`] can fail or give a part of a record from past NCBI's mapping of
    /// `taxdb.btd` for some taxid ([`TaxLookup::PastMapping`], a kingdom `Err`), decided by one
    /// sequential pass over the index alone (the checks before a search,
    /// `algorithm::common::database::taxdb_may_refuse`): only a lookup that reads a record past the
    /// mapping (the end of the file's last 4096-byte page; every record of an empty file, which
    /// NCBI maps as a null pointer) does. Besides a read error, a lookup reads the record
    /// `offset[i]..offset[i+1]` of the entry `i` that its search ends on (below the search's last
    /// `high_index`; one that ends on it reads to the end of `taxdb.btd`, inside the mapping). So
    /// `Ok(false)`, whatever the order of the index, when no two consecutive entries make a record
    /// of positive size that reaches past the mapping; a database that NCBI treats as missing finds
    /// nothing. The index is read in blocks of 64 Ki entries.
    pub fn any_lookup_may_fail(&self) -> Result<bool, SeqDbError> {
        const INDEX_BLOCK: u64 = 1 << 16;
        let Some(files) = &self.files else {
            return Ok(false);
        };
        let mapped = super::file::mapped_length(files.data.len());
        let count = files.count as u64;
        let mut previous: Option<u32> = None;
        let mut block: Vec<u8> = Vec::new();
        let mut i = 0u64;
        while i < count {
            let n = INDEX_BLOCK.min(count - i);
            block.resize((n * ENTRY_LEN) as usize, 0);
            files
                .index
                .read_exact_at(&mut block, HEADER_LEN + ENTRY_LEN * i)?;
            for entry in block.chunks_exact(ENTRY_LEN as usize) {
                let offset = be32(&entry[4..8]);
                if let Some(begin) = previous {
                    // `lookup`: `int(m_End - m_Begin)` of the wrapped `Uint4` difference.
                    let size = offset.wrapping_sub(begin) as i32;
                    if size > 0 && u64::from(begin) + size as u64 > mapped {
                        return Ok(true);
                    }
                }
                previous = Some(offset);
            }
            i += n;
        }
        Ok(false)
    }

    /// Whether [`Self::lookup`] faults (`Err`) for some taxid (the report streams' check before
    /// the search, `algorithm::common::database::reports_may_refuse`, whose prefetch writes no
    /// name): only where `taxdb.btd` is empty (NCBI maps it as a null pointer, [`Self::maps_null`])
    /// and a lookup reads a record of positive size ([`Self::any_lookup_may_fail`], whose mapping is
    /// then empty); `Ok(false)` for a `taxdb.btd` that is not empty.
    pub fn any_lookup_faults(&self) -> Result<bool, SeqDbError> {
        if !self.maps_null() {
            return Ok(false);
        }
        self.any_lookup_may_fail()
    }

    /// Whether NCBI maps `taxdb.btd` as a null pointer (an empty file): a lookup that reads a
    /// record faults ([`Self::lookup`]). Only then can a lookup whose names reach no output (the
    /// texts of a `GetBioseq`) fail; otherwise it writes at most the offset error of
    /// [`Self::offset_error_taxid`].
    pub fn maps_null(&self) -> bool {
        self.files
            .as_ref()
            .is_some_and(|files| files.data.is_empty())
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
/// The three splits of `bytes` at their first tabs (seqdbtax.cpp:279-283): scientific name,
/// common name, BLAST name, and what is left after the third tab; `None` when `bytes` holds fewer
/// than three tabs.
fn split_names(bytes: &[u8]) -> Option<([&[u8]; 3], &[u8])> {
    let mut rest = bytes;
    let mut fields: [&[u8]; 3] = [&[]; 3];
    for field in &mut fields {
        let i = rest.iter().position(|&b| b == b'\t')?;
        *field = &rest[..i];
        rest = &rest[i + 1..];
    }
    Some((fields, rest))
}

/// The split of one record into four fields (seqdbtax.cpp:279-298): scientific name, common
/// name, BLAST name, kingdom; `None` when NCBI's lookup returns false for the record.
fn split_fields(record: &[u8]) -> Option<[&[u8]; 4]> {
    let ([sci, common, blast], kingdom) = split_names(record)?;
    (!kingdom.is_empty()).then_some([sci, common, blast, kingdom])
}

/// The split of one record into four fields (seqdbtax.cpp:279-298).
fn split_record(taxid: i32, record: &[u8]) -> Option<TaxNames> {
    split_fields(record).map(|fields| TaxNames {
        taxid,
        scientific_name: fields[0].to_vec(),
        common_name: fields[1].to_vec(),
        blast_name: fields[2].to_vec(),
        super_kingdom: Ok(fields[3].to_vec()),
    })
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
    /// Whether the next call hands over the texts of opening again ([`LazyTaxDb::restart`]).
    again: std::sync::atomic::AtomicBool,
}

impl LazyTaxDb {
    /// A handle that will open the database through `search`.
    pub fn new(search: SearchPath) -> Self {
        Self {
            search,
            db: OnceLock::new(),
            again: std::sync::atomic::AtomicBool::new(false),
        }
    }

    /// Back to a handle that has handed over nothing: the next call hands over the texts of
    /// opening, as a new handle's first would (a run that searches the batches before a failed
    /// one again, `algorithm::common::database::TaxdbRun::restart`). The database read is kept:
    /// a new handle opens the same files to the same answers.
    pub fn restart(&self) {
        self.again.store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Opens the database if no call has (`CTaxDBFileInfo`, at the first `GetTaxNames`) and returns
    /// the `ERR_POST` texts of opening when this call opened it, none when it was open already.
    /// An open that was refused is the same error on every call.
    pub fn open(&self) -> Result<Vec<String>, SeqDbError> {
        let mut messages = Vec::new();
        self.open_into(&mut messages)?;
        Ok(messages)
    }

    /// Opens the database like [`LazyTaxDb::open`], appending the texts of opening to `messages`.
    fn open_into(&self, messages: &mut Vec<String>) -> Result<&TaxDb, SeqDbError> {
        let mut opened = false;
        let db = self.db.get_or_init(|| {
            opened = true;
            TaxDb::open(&self.search).map_err(|error| match error {
                SeqDbError::Unsupported { path, message } => (path, message),
                // The file that failed, and what the system said of it.
                SeqDbError::Open { path, source } | SeqDbError::Io { path, source } => {
                    (path, source.to_string())
                }
                SeqDbError::Format { path, message } => (path, message),
                other => (PathBuf::from("taxdb.bti"), other.to_string()),
            })
        });
        match db {
            Ok(db) => {
                if opened || self.handed_again() {
                    messages.extend(db.open_messages().iter().cloned());
                }
                Ok(db)
            }
            Err((path, message)) => Err(SeqDbError::Unsupported {
                path: path.clone(),
                message: message.clone(),
            }),
        }
    }

    /// Whether this call hands over the texts of opening again after [`LazyTaxDb::restart`] (the
    /// first call after it, once).
    fn handed_again(&self) -> bool {
        use std::sync::atomic::Ordering;
        self.again.load(Ordering::Relaxed) && self.again.swap(false, Ordering::Relaxed)
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
        self.open_into(messages)?.names(tax_id, messages)
    }

    /// [`TaxDb::lookup`] through the lazily opened database, the texts as [`LazyTaxDb::names`]
    /// appends them.
    pub fn lookup(&self, tax_id: i32, messages: &mut Vec<String>) -> Result<TaxLookup, SeqDbError> {
        self.open_into(messages)?.lookup(tax_id, messages)
    }

    /// The lazily opened database (this opens it, dropping the texts of opening: ask it after the
    /// run's first lookup). An open that was refused is the same error on every call.
    pub fn database(&self) -> Result<&TaxDb, SeqDbError> {
        self.open_into(&mut Vec::new())
    }

    /// [`TaxDb::entry_count`] of the lazily opened database (this opens it).
    pub fn entry_count(&self) -> Result<u64, SeqDbError> {
        Ok(self.open_into(&mut Vec::new())?.entry_count())
    }

    /// [`TaxDb::any_lookup_faults`] of the lazily opened database (this opens it). An open that
    /// was refused or a file that cannot be read is an `Err`.
    pub fn any_lookup_faults(&self) -> Result<bool, SeqDbError> {
        self.open_into(&mut Vec::new())?.any_lookup_faults()
    }

    /// [`TaxDb::any_lookup_may_fail`] of the lazily opened database (this opens it). An open that
    /// was refused or a file that cannot be read is an `Err`.
    pub fn any_lookup_may_fail(&self) -> Result<bool, SeqDbError> {
        self.open_into(&mut Vec::new())?.any_lookup_may_fail()
    }

    /// Whether the database has been opened (a name has been asked for).
    #[cfg(test)]
    pub fn is_open(&self) -> bool {
        self.db.get().is_some()
    }
}

/// What the lookup of one taxid gives the source of a Bioseq (`CSeqDBVol::x_GetTaxonomy`'s
/// Org-ref, seqdbvol.cpp:1161-1169): the names of the taxonomy database, and the scientific name
/// that `CDeflineGenerator` takes as `m_Taxname`.
#[derive(Debug)]
pub struct SourceNames {
    names: Option<TaxNames>,
    taxname: Option<std::sync::Arc<[u8]>>,
}

impl SourceNames {
    fn new(names: Option<TaxNames>) -> Self {
        // An empty scientific name sets an empty `m_Taxname`, which makes the titles of a Bioseq
        // without a name (`report::defline::SourceTaxname::Unnamed`).
        let taxname = names
            .as_ref()
            .filter(|names| !names.scientific_name.is_empty())
            .map(|names| std::sync::Arc::from(names.scientific_name.as_slice()));
        Self { names, taxname }
    }

    /// The names (`TaxDb::names`); `None` where `GetTaxNames` returned false.
    pub fn names(&self) -> Option<&TaxNames> {
        self.names.as_ref()
    }

    /// The scientific name of the source (`m_Taxname`); `None` where the lookup found no names or
    /// the name is empty.
    pub fn taxname(&self) -> Option<&std::sync::Arc<[u8]>> {
        self.taxname.as_ref()
    }
}

/// `CSeqDBVol::m_TaxCache(256)` (seqdbvol.cpp:84).
const TAX_CACHE_SLOTS: usize = 256;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1105-1107, 1140-1182
// ```c++
//     const bool provide_new_taxonomy_info = true;
//     const bool use_taxinfo_cache         = (CThread::GetSelf() == 0 ? true : false);
//     ...
//         if (use_taxinfo_cache && m_TaxCache.Lookup(TAX_ID_TO(int, taxid)).NotEmpty()) {
//             have_org_desc = true;
//         }
//
//         SSeqDBTaxInfo tnames(taxid);
//         bool found_taxid_in_taxonomy_blastdb = true;
//
//         if ((! have_org_desc) && provide_new_taxonomy_info) {
//             try {
//                 found_taxid_in_taxonomy_blastdb = CSeqDBTaxInfo::GetTaxNames(taxid, tnames);
//             } catch (CSeqDBException &) {
//                 found_taxid_in_taxonomy_blastdb = false;
//             }
//         }
//         ...
//                 if (use_taxinfo_cache) {
//                     m_TaxCache.Lookup(TAX_ID_TO(int, taxid)) = desc;
//                 }
// ```
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbgeneral.hpp:1184-1196
// ```c++
//     TValue & Lookup(int key)
//     {
//         _ASSERT(IS_POWER_OF_TWO(m_Slots.size()));
//
//         TSlot & slot = m_Slots[key & (m_Slots.size()-1)];
//
//         if (slot.first != key) {
//             slot.first = key;
//             slot.second = TValue();
//         }
//
//         return slot.second;
//     }
// ```
// NCBI's cache keeps the source of every lookup on the main thread, found or not (a
// `CSeqDBException` of `GetTaxNames` counts as not found), so a volume asks the taxonomy database
// again only where the slot of the taxid holds another one; on other threads it asks at every
// load. The texts that a lookup writes come from `GetTaxNames` and so only where it is asked:
// `algorithm::common::database::TaxdbRun` keeps that per volume for the texts of the loads (the
// search's and the readers'). A lookup that faults (a record of an empty `taxdb.btd`, which NCBI
// maps as a null pointer) or whose names NCBI reads past its mapping of the file (no defined
// result) is refused at the first by the reader that writes them: there is nothing to keep.
/// The sources' names that a reader of a database's Bioseqs looked up, by taxid: NCBI's taxonomy
/// cache of a volume (`CSeqDBVol::m_TaxCache`, 256 slots, a taxid in the slot `taxid & 255`, which
/// holds the last taxid looked up there). The one memo of the lookups of the readers of the stored
/// titles (`algorithm::common::database::HeaderReader`, blastdbcmd's `retrieve::Source`, the
/// maskers' `seqmasks_io::bdb_reader::MaskBdbReader`), each with its own.
///
/// It keeps only a lookup that found its answer without writing a text and without an error,
/// which every later lookup of that taxid gives again (`TaxDb::names` reads nothing else). A
/// lookup that writes a text (those of opening the database, which come with the first lookup;
/// the offset error of the last entry) or fails is asked again each time, as NCBI writes the
/// offset error at each lookup where it asks: the texts are those of the loads, which
/// `algorithm::common::database::TaxdbRun` writes for the search and the readers alike, and the
/// readers' lookups of names leave them out. So its answers do not depend on what it holds, and
/// it holds at most 256 of them (nothing until the first lookup).
#[derive(Debug, Default)]
pub struct TaxCache {
    slots: std::cell::RefCell<Vec<Option<(i32, std::sync::Arc<SourceNames>)>>>,
}

impl TaxCache {
    /// The names of `taxid`: the slot's where it holds `taxid`, else `lookup` (`TaxDb::names` or
    /// `LazyTaxDb::names` of the reader's database, or a refusal of the reader's own, `E`), which
    /// appends its texts to `messages`.
    pub fn names<E>(
        &self,
        taxid: i32,
        messages: &mut Vec<String>,
        lookup: impl FnOnce(&mut Vec<String>) -> Result<Option<TaxNames>, E>,
    ) -> Result<std::sync::Arc<SourceNames>, E> {
        let slot = taxid as usize & (TAX_CACHE_SLOTS - 1);
        if let Some(Some((held, names))) = self.slots.borrow().get(slot) {
            if *held == taxid {
                return Ok(std::sync::Arc::clone(names));
            }
        }
        let written = messages.len();
        let names = std::sync::Arc::new(SourceNames::new(lookup(messages)?));
        if messages.len() == written {
            let mut slots = self.slots.borrow_mut();
            if slots.is_empty() {
                slots.resize(TAX_CACHE_SLOTS, None);
            }
            slots[slot] = Some((taxid, std::sync::Arc::clone(&names)));
        }
        Ok(names)
    }
}

#[cfg(test)]
mod tests {
    //! Tests against taxdb files written by `make_taxdb.py` (the fixtures under
    //! `tests/seqdb/taxdb/<variant>/`, with what NCBI BLAST+ 2.17.0 returned
    //! for them in `tests/seqdb/taxdb/expected.tsv`) and against files built
    //! here for the cases NCBI cannot be asked about.

    use super::super::file::MAP_PAGE;
    use super::*;
    use std::path::Path;

    fn fixture_dir() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb")
    }

    /// A directory made for one test, removed (with what is in it) when it is dropped: a test
    /// drops it after the `TaxDb` and the `SearchPath` that read it.
    struct TempDir(PathBuf);

    impl std::ops::Deref for TempDir {
        type Target = Path;

        fn deref(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn temp_dir(tag: &str) -> TempDir {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-taxdb-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }

    /// A search path with the empty current directory it was made with, which is removed with it.
    struct Searching {
        search: SearchPath,
        _cwd: TempDir,
    }

    impl std::ops::Deref for Searching {
        type Target = SearchPath;

        fn deref(&self) -> &SearchPath {
            &self.search
        }
    }

    /// A search path whose only useful directory is `dir` (`BLASTDB`), run from an empty current directory.
    fn path_of(dir: &Path) -> Searching {
        let cwd = temp_dir("cwd");
        let search = SearchPath::new(&cwd, Some(dir.as_os_str()), None).unwrap();
        Searching { search, _cwd: cwd }
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

    fn write_db(index: &[u8], data: &[u8]) -> TempDir {
        let dir = temp_dir("db");
        std::fs::write(dir.join("taxdb.bti"), index).unwrap();
        std::fs::write(dir.join("taxdb.btd"), data).unwrap();
        dir
    }

    /// Writes a taxdb and opens it: the directory is returned first, so that it is dropped after
    /// the open database.
    fn open_built(index: &[u8], data: &[u8]) -> (TempDir, TaxDb) {
        let dir = write_db(index, data);
        let db = TaxDb::open(&path_of(&dir)).unwrap();
        (dir, db)
    }

    fn names_of(db: &TaxDb, taxid: i32) -> (Option<TaxNames>, Vec<String>) {
        let mut messages = Vec::new();
        let names = db.names(taxid, &mut messages).unwrap();
        (names, messages)
    }

    // seqdbtax.cpp:258-259 (`info.taxid = tax_id` once the search finds the taxid, before the
    // lookup can still fail): blastdbcmd's `-tax_info` writes empty names then, `N/A` only for a
    // taxid the search does not find (probe of line RG T5b, `T/t5b/out/extra`).
    #[test]
    fn a_lookup_tells_a_taxid_not_found_from_one_found_without_names() {
        let lookup = |variant: &str, taxid: i32| {
            let mut messages = Vec::new();
            let found = open_fixture(variant)
                .unwrap()
                .lookup(taxid, &mut messages)
                .unwrap();
            (found, messages)
        };
        assert!(matches!(lookup("good", 9606), (TaxLookup::Named(_), m) if m.is_empty()));
        for taxid in [150, 0, 77777, 100000] {
            assert_eq!(lookup("good", taxid), (TaxLookup::NotInIndex, Vec::new()));
        }
        // A record without three tabs, an empty record of the last entry, the offset error.
        assert_eq!(lookup("tabs0", 9606), (TaxLookup::Unnamed, Vec::new()));
        assert_eq!(lookup("tabs1", 9606), (TaxLookup::Unnamed, Vec::new()));
        assert_eq!(
            lookup("offset_last_at_end", 99990),
            (TaxLookup::Unnamed, Vec::new())
        );
        assert_eq!(
            lookup("offset_last_past_end", 99990),
            (TaxLookup::Unnamed, vec![MSG_OFFSET_ERROR.to_string()])
        );
        // A database NCBI treats as missing finds no taxid.
        assert_eq!(
            lookup("bad_magic", 9606),
            (TaxLookup::NotInIndex, Vec::new())
        );
        assert_eq!(lookup("no_bti", 9606), (TaxLookup::NotInIndex, Vec::new()));
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
            names
                .super_kingdom
                .as_deref()
                .expect("a kingdom inside the mapping"),
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
            let (_dir, db) = open_built(&index, &data);
            for i in 0..n {
                let names = names_of(&db, 10 * i + 5)
                    .0
                    .unwrap_or_else(|| panic!("n={n} i={i}"));
                assert_eq!(names.scientific_name, format!("s{i}").into_bytes());
                assert_eq!(names.super_kingdom, Ok(format!("k{i}").into_bytes()));
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
        let (_dir, db) = open_built(&index, &data);
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
            let (_dir, db) = open_built(&index, &data);
            assert_eq!(names_of(&db, 2).0.is_some(), found, "{record:?}");
            assert!(names_of(&db, 1).0.is_some());
        }
    }

    #[test]
    fn last_offset_past_the_end_is_an_error_message_and_not_found() {
        let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th")]);
        let at = 24 + 8 + 4; // offset of the second entry
        index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let (_dir, db) = open_built(&index, &data);
        let (names, messages) = names_of(&db, 2);
        assert_eq!(names, None);
        assert_eq!(messages, [MSG_OFFSET_ERROR]);
        // The message comes at every lookup.
        assert_eq!(names_of(&db, 2).1, [MSG_OFFSET_ERROR]);
        // An offset equal to the size is an empty record: not found, no message.
        index[at..at + 4].copy_from_slice(&be(data.len() as u32));
        let (_dir, db) = open_built(&index, &data);
        assert_eq!(names_of(&db, 2), (None, vec![]));
    }

    #[test]
    fn the_taxid_of_the_offset_error_is_the_last_entry_s_where_its_lookup_writes_it() {
        // Every variant: the taxid for which `blastdbcmd -tax_info` of NCBI BLAST+ 2.17.0 wrote
        // the message (`expected.tsv`), none where it wrote it for no taxid.
        let text = std::fs::read_to_string(fixture_dir().join("expected.tsv")).unwrap();
        let mut variants: Vec<String> = std::fs::read_dir(fixture_dir())
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.is_dir())
            .map(|path| path.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        variants.sort();
        let mut with_message = 0;
        for variant in &variants {
            let Ok(db) = open_fixture(variant) else {
                continue;
            };
            let ncbi: Vec<i32> = text
                .lines()
                .map(|line| line.splitn(4, '\t').collect::<Vec<&str>>())
                .filter(|cols| {
                    cols[0] == variant.as_str()
                        && cols.get(3).is_some_and(|m| m.contains(MSG_OFFSET_ERROR))
                })
                .map(|cols| cols[1].parse().unwrap())
                .collect();
            assert!(ncbi.len() <= 1, "{variant}");
            assert_eq!(
                db.offset_error_taxid().unwrap(),
                ncbi.first().copied(),
                "{variant}"
            );
            with_message += ncbi.len();
        }
        // offset_last_past_end, btd_truncated and btd_empty.
        assert_eq!(with_message, 3);
        // The last entry (50) begins past the end, but the search for 50 ends on entry 2 (90 >
        // 50, then 20 < 50): no lookup writes the message.
        let rows: [(i32, &[u8]); 5] = [
            (10, b"a\tb\tc\td"),
            (20, b"a\tb\tc\td"),
            (90, b"a\tb\tc\td"),
            (30, b"a\tb\tc\td"),
            (50, b"a\tb\tc\td"),
        ];
        let (mut index, data) = build(&rows);
        let at = 24 + 8 * 4 + 4;
        index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let (_dir, db) = open_built(&index, &data);
        assert_eq!(names_of(&db, 50), (None, vec![]));
        assert_eq!(db.offset_error_taxid().unwrap(), None);
        // Sorted, the same entry writes it.
        let rows: [(i32, &[u8]); 3] = [(10, b"a\tb\tc\td"), (20, b"a\tb\tc\td"), (50, b"x")];
        let (mut index, data) = build(&rows);
        let at = 24 + 8 * 2 + 4;
        index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 1));
        let (_dir, db) = open_built(&index, &data);
        assert_eq!(db.offset_error_taxid().unwrap(), Some(50));
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
        let (_dir, db) = open_built(&index, &data);
        assert_eq!(names_of(&db, 1), (None, vec![]));
        // Entry 2 is now [0, 16): both records, the kingdom takes the rest.
        let names = names_of(&db, 2).0.unwrap();
        assert_eq!(names.scientific_name, b"a");
        assert_eq!(
            names.super_kingdom,
            Ok(b"d"
                .to_vec()
                .into_iter()
                .chain(*b"e\tf\tg\th")
                .collect::<Vec<u8>>())
        );
    }

    #[test]
    fn record_reaching_past_the_end_reads_page_padding_and_then_memory_past_the_mapping() {
        let (mut index, data) =
            build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
        let off2 = 24 + 8 + 4;
        let off3 = 24 + 16 + 4;
        // Entry 2 ends 10 bytes after the file: zero padding belongs to its kingdom.
        index[off3..off3 + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let (_dir, db) = open_built(&index, &data);
        let names = names_of(&db, 2).0.unwrap();
        let mut kingdom = b"h".to_vec();
        kingdom.extend(*b"i\tj\tk\tl");
        kingdom.extend([0u8; 10]);
        assert_eq!(names.super_kingdom, Ok(kingdom));
        // Entry 3 (the last) begins after the end: NCBI's message.
        assert_eq!(names_of(&db, 3), (None, vec![MSG_OFFSET_ERROR.to_string()]));
        // Past the page: the three names before the third tab inside the mapping, a kingdom
        // NCBI reads past it. Entry 1 is fine.
        index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 1));
        let (_dir, db) = open_built(&index, &data);
        let past = PastMapping {
            path: db.files.as_ref().unwrap().data.path().to_path_buf(),
            taxid: 2,
            begin: 7,
            end: MAP_PAGE + 1,
            len: data.len() as u64,
        };
        let mut messages = Vec::new();
        assert_eq!(
            db.lookup(2, &mut messages).unwrap(),
            TaxLookup::Named(TaxNames {
                taxid: 2,
                scientific_name: b"e".to_vec(),
                common_name: b"f".to_vec(),
                blast_name: b"g".to_vec(),
                super_kingdom: Err(past.clone()),
            })
        );
        assert!(messages.is_empty());
        assert!(names_of(&db, 1).0.is_some());
        // A record past the mapping without its third tab inside it: no defined result, which
        // only a caller that writes the names refuses.
        index[off2..off2 + 4].copy_from_slice(&be(data.len() as u32 - 2));
        let (_dir, db) = open_built(&index, &data);
        let past = PastMapping {
            path: db.files.as_ref().unwrap().data.path().to_path_buf(),
            begin: data.len() as u64 - 2,
            ..past
        };
        assert_eq!(
            db.lookup(2, &mut messages).unwrap(),
            TaxLookup::PastMapping(past.clone())
        );
        assert!(matches!(
            db.names(2, &mut messages),
            Err(SeqDbError::Unsupported { .. })
        ));
        assert!(messages.is_empty());
        // A record that begins past the mapping.
        index[off2..off2 + 4].copy_from_slice(&be(MAP_PAGE as u32));
        index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 9));
        let (_dir, db) = open_built(&index, &data);
        assert!(matches!(
            db.lookup(2, &mut messages).unwrap(),
            TaxLookup::PastMapping(PastMapping {
                begin: 4096,
                end: 4105,
                ..
            })
        ));
        // Neither faults for certain, nor does the pass over the index say so; an empty
        // taxdb.btd does.
        assert!(!db.any_lookup_faults().unwrap());
        assert!(!db.maps_null());
        let empty = open_fixture("btd_empty").unwrap();
        assert!(empty.any_lookup_faults().unwrap());
        assert!(empty.maps_null());
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
        let (_dir, db) = open_built(&index, &data);
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
        assert_eq!(names.super_kingdom.as_deref(), Ok(&b"Euk\xffaryota"[..]));
    }

    #[test]
    fn lazy_handle_opens_once_and_hands_over_the_open_messages_once() {
        let dir = fixture_dir().join("bad_magic");
        let searching = path_of(&dir);
        let lazy = LazyTaxDb::new(searching.search.clone());
        assert!(!lazy.is_open());
        let mut messages = Vec::new();
        assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
        assert!(lazy.is_open());
        assert_eq!(messages, [MSG_WRONG_MAGIC]);
        assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
        assert_eq!(messages, [MSG_WRONG_MAGIC]);
        // A missing database: opened at the first name, silent, finds nothing.
        let searching = path_of(&fixture_dir().join("nonexistent"));
        let lazy = LazyTaxDb::new(searching.search.clone());
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

    /// What [`TaxDb::any_lookup_may_fail`] must agree with: the lookup of every taxid of the
    /// index, failing (`Err`) or giving a part of its record from past NCBI's mapping of
    /// `taxdb.btd` (`TaxLookup::PastMapping`, a kingdom `Err`); `names` then refuses the first.
    fn brute_force_may_fail(db: &TaxDb) -> bool {
        let Some(files) = &db.files else {
            return false;
        };
        (0..files.count).any(|i| {
            let taxid = TaxDb::entry(files, i).unwrap().0;
            !matches!(
                db.lookup(taxid, &mut Vec::new()),
                Ok(TaxLookup::NotInIndex
                    | TaxLookup::Unnamed
                    | TaxLookup::Named(TaxNames {
                        super_kingdom: Ok(_),
                        ..
                    }))
            )
        })
    }

    /// Whether the index of `db` is in the layout the pass relies on (strictly increasing taxids,
    /// non-decreasing offsets).
    fn regular_layout(db: &TaxDb) -> bool {
        let Some(files) = &db.files else {
            return true;
        };
        let entries: Vec<(i32, u32)> = (0..files.count)
            .map(|i| TaxDb::entry(files, i).unwrap())
            .collect();
        entries
            .windows(2)
            .all(|pair| pair[0].0 < pair[1].0 && pair[0].1 <= pair[1].1)
    }

    /// The pass over the index answers as the lookups do on every fixture variant: never `false`
    /// where the lookup of a taxid of the index, of a taxid next to one, or of the extremes fails;
    /// and exactly the lookups' answer on the taxids of the index when the index is sorted. Of the
    /// variants, `btd_empty` fails (every record past the mapped end); a record cut by the end of
    /// the file inside its last page does not (`btd_truncated`, `offset_next_past_end`).
    #[test]
    fn the_pass_over_the_index_agrees_with_the_failing_lookups_on_every_variant() {
        let (mut failing, mut clean) = (Vec::new(), 0);
        let mut dirs: Vec<_> = std::fs::read_dir(fixture_dir())
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.is_dir())
            .collect();
        dirs.sort();
        for dir in dirs {
            let variant = dir.file_name().unwrap().to_string_lossy().to_string();
            let Ok(db) = open_fixture(&variant) else {
                continue;
            };
            let pass = db.any_lookup_faults().unwrap();
            let fails = |taxid: i32| db.lookup(taxid, &mut Vec::new()).is_err();
            let indexed: Vec<i32> = db.files.as_ref().map_or(Vec::new(), |files| {
                (0..files.count)
                    .map(|i| TaxDb::entry(files, i).unwrap().0)
                    .collect()
            });
            let near = indexed
                .iter()
                .flat_map(|&taxid| [taxid.wrapping_sub(1), taxid.wrapping_add(1)])
                .chain([i32::MIN, -1, 0, 1, i32::MAX]);
            assert!(
                pass || !indexed.iter().copied().chain(near).any(fails),
                "{variant}: the pass says no, a lookup fails"
            );
            if regular_layout(&db) {
                assert_eq!(pass, indexed.iter().copied().any(fails), "{variant}");
            }
            if pass {
                failing.push(variant);
            } else {
                clean += 1;
            }
        }
        assert!(
            failing.iter().any(|variant| variant == "btd_empty") && clean >= 20,
            "{failing:?} {clean}"
        );
    }

    /// The one pass over the files answers as the lookup of every taxid does on every fixture
    /// variant: never `false` where a lookup refuses, and exactly the same when the index has the
    /// layout the pass relies on.
    #[test]
    fn one_pass_over_the_files_agrees_with_the_lookups_on_every_variant() {
        let mut refusing = 0;
        let mut clean = 0;
        let mut dirs: Vec<_> = std::fs::read_dir(fixture_dir())
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.is_dir())
            .collect();
        dirs.sort();
        for dir in dirs {
            let variant = dir.file_name().unwrap().to_string_lossy().to_string();
            let Ok(db) = open_fixture(&variant) else {
                continue;
            };
            let brute = brute_force_may_fail(&db);
            let walk = db.any_lookup_may_fail().unwrap();
            assert!(
                walk || !brute,
                "{variant}: the pass says no, a lookup refuses"
            );
            if regular_layout(&db) {
                assert_eq!(walk, brute, "{variant}");
            }
            if brute {
                refusing += 1;
            } else {
                clean += 1;
            }
        }
        // `btd_empty` refuses (an empty taxdb.btd maps no page); a scientific name of any bytes is
        // no refusal (the two `high_bytes_*` variants are clean); the test of the index's blocks
        // (`the_pass_counts_entries_and_reads_past_one_index_block`) refuses past the first.
        assert!(refusing >= 1 && clean >= 20, "{refusing} {clean}");
    }

    #[test]
    fn the_pass_counts_entries_and_reads_past_one_index_block() {
        assert_eq!(open_fixture("good").unwrap().entry_count(), 7);
        assert_eq!(open_fixture("no_btd").unwrap().entry_count(), 0);
        assert_eq!(open_fixture("bad_magic").unwrap().entry_count(), 0);
        // 70,000 entries (more than one block of the index); in the damaged copy the record of an
        // entry near the end reaches past NCBI's mapping of taxdb.btd (the two entries after it
        // begin past the file's last page), whose lookup gives its kingdom from past the mapping.
        let n = 70_000u32;
        let mut data = Vec::new();
        let mut offsets = Vec::new();
        for i in 0..n {
            offsets.push(data.len() as u32);
            data.extend(format!("name {i}\tc\tb\tK\n").as_bytes());
        }
        let past = (data.len().div_ceil(MAP_PAGE as usize) * MAP_PAGE as usize + 4096) as u32;
        let mut damaged = offsets.clone();
        damaged[n as usize - 2] = past;
        damaged[n as usize - 1] = past;
        for (tag, offsets, refuses) in [("walk", &offsets, false), ("walkpast", &damaged, true)] {
            let dir = temp_dir(tag);
            let mut index = Vec::new();
            index.extend(MAGIC.to_be_bytes());
            index.extend(n.to_be_bytes());
            index.extend([0u8; 16]);
            for (i, offset) in offsets.iter().enumerate() {
                index.extend((i as i32 + 1).to_be_bytes());
                index.extend(offset.to_be_bytes());
            }
            std::fs::write(dir.join("taxdb.bti"), &index).unwrap();
            std::fs::write(dir.join("taxdb.btd"), &data).unwrap();
            let db = TaxDb::open(&path_of(&dir)).unwrap();
            assert_eq!(db.entry_count(), u64::from(n));
            assert_eq!(db.any_lookup_may_fail().unwrap(), refuses, "{tag}");
            assert_eq!(brute_force_may_fail(&db), refuses, "{tag}");
        }
    }

    /// `TaxCache::names` of `taxid` through `db`, with the number of lookups it made.
    fn cached(
        cache: &TaxCache,
        db: &TaxDb,
        taxid: i32,
    ) -> (
        Result<std::sync::Arc<SourceNames>, SeqDbError>,
        Vec<String>,
        u32,
    ) {
        let mut messages = Vec::new();
        let mut asked = 0;
        let names = cache.names(taxid, &mut messages, |messages| {
            asked += 1;
            db.names(taxid, messages)
        });
        (names, messages, asked)
    }

    // NCBI's `CSeqDBIntCache` (seqdbgeneral.hpp:1184-1196): a slot per `taxid & 255` that holds the
    // last taxid looked up there. LOSAT keeps only the answers that a lookup gives again without
    // a text: those that wrote one, and the refusals, are asked again.
    #[test]
    fn the_cache_keeps_answers_without_texts_in_a_slot_per_taxid() {
        // 257 shares the slot of 1; 300 is in no entry (found nothing); 600 is the last entry,
        // whose offset is past the end (the offset error at each lookup).
        let (mut index, data) = build(&[
            (1, b"Bacillus\tb\tc\td"),
            (257, b"\tb\tc\td"),
            (600, b"e\tf\tg\th"),
        ]);
        let last = 24 + 16 + 4;
        index[last..last + 4].copy_from_slice(&be(data.len() as u32 + 10));
        let (_dir, db) = open_built(&index, &data);
        let cache = TaxCache::default();
        let (first, messages, asked) = cached(&cache, &db, 1);
        let first = first.unwrap();
        assert_eq!((messages.len(), asked), (0, 1));
        assert_eq!(
            first.taxname().map(|name| &name[..]),
            Some(&b"Bacillus"[..])
        );
        let (again, _, asked) = cached(&cache, &db, 1);
        assert!(std::sync::Arc::ptr_eq(&first, &again.unwrap()));
        assert_eq!(asked, 0);
        // An empty scientific name is no `m_Taxname`; the names are still there.
        let (other, _, asked) = cached(&cache, &db, 257);
        let other = other.unwrap();
        assert_eq!(
            (other.taxname(), other.names().is_some(), asked),
            (None, true, 1)
        );
        // 257 took the slot: 1 is asked again.
        assert_eq!(cached(&cache, &db, 1).2, 1);
        // Nothing found is kept.
        let (none, _, asked) = cached(&cache, &db, 300);
        assert!(none.unwrap().names().is_none());
        assert_eq!(asked, 1);
        assert_eq!(cached(&cache, &db, 300).2, 0);
        // The offset error comes at each lookup: never kept.
        for _ in 0..2 {
            let (names, messages, asked) = cached(&cache, &db, 600);
            assert!(names.unwrap().names().is_none());
            assert_eq!((messages, asked), (vec![MSG_OFFSET_ERROR.to_string()], 1));
        }
        // A refused lookup (a record that lies past NCBI's mapping of taxdb.btd, whose names are
        // not defined) is asked again.
        let (mut index, data) =
            build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
        let off2 = 24 + 8 + 4;
        let off3 = 24 + 16 + 4;
        index[off2..off2 + 4].copy_from_slice(&be(MAP_PAGE as u32));
        index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 9));
        let (_dir, db) = open_built(&index, &data);
        let cache = TaxCache::default();
        for _ in 0..2 {
            let (names, _, asked) = cached(&cache, &db, 2);
            assert!(matches!(names, Err(SeqDbError::Unsupported { .. })));
            assert_eq!(asked, 1);
        }
    }
}
