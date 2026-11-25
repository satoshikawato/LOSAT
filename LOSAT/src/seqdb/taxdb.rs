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
//! (`SeqDB_ResolveDbPath` of `"taxdb.bti"`: the first of the current directory,
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
/// A post of NCBI's main thread: its header alone after the stream failed
/// (`ncbi_exception::post_line`).
pub fn err_post_line(program: &str, text: &str) -> String {
    let mut line = String::from_utf8(crate::ncbi_exception::post_line(
        crate::ncbi_exception::Severity::Error,
        None,
        Some(program),
        text.as_bytes(),
    ))
    .expect("UTF-8 parts");
    line.pop();
    line
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
    /// False after the search found it (`info.taxid` is set to `tax_id`; seqdbtax.cpp:258-259): the
    /// offset error of the last entry, a record of no positive size or not of four fields.
    Unnamed,
    /// True: the names.
    Named(TaxNames),
    /// The search found the taxid (`info.taxid` is set to `tax_id`) and its record reaches past
    /// NCBI's mapping of `taxdb.btd` before its third tab: NCBI looks for the tabs in memory past
    /// the mapping, so whether it finds names, and which, depends on its memory map (or it faults).
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:125-203
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
        let Some(index_name) = search.resolve_file(b"taxdb.bti") else {
            return missing(Vec::new());
        };
        let mut data_name = index_name.clone();
        data_name.pop();
        data_name.push(b'd');
        let index_path = crate::utils::os_bytes::path_from_bytes(index_name);
        let data_path = crate::utils::os_bytes::path_from_bytes(data_name);
        if !super::alias::exists(&index_path) || !super::alias::exists(&data_path) {
            return missing(Vec::new());
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:154-161
        // NCBI takes the index's length from the file system before it maps the index: a short
        // index is missing even where it does not open.
        let index_len = std::fs::metadata(&index_path)
            .map_err(|source| SeqDbError::Io {
                path: index_path.clone(),
                source,
            })?
            .len();
        if index_len < HEADER_LEN + ENTRY_LEN {
            return missing(Vec::new());
        }
        let index = DbFile::open(&index_path).map_err(throws_at_every_lookup)?;
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
        let data = DbFile::open(&data_path).map_err(throws_at_every_lookup)?;
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
    // The index is not assumed to be sorted, as in NCBI: the search is NCBI's
    // loop with `Int4` taxids and `Int4` index arithmetic, and "the last
    // entry" is whatever `high_index` ends on (for a sorted index only the
    // last entry; for a damaged one `end_data` may then be the file size).
    // `CSeqDB_Substring::Size()` is `m_End` minus `m_Begin` cast to `int`
    // (seqdbgeneral.hpp:375-378), so a record whose end precedes its beginning
    // (`end_data` minus `begin_data` wraps to 2^32-k as `Uint4`) has a negative size, as has a
    // record of 2 GiB or more: `SeqDB_SplitString` finds nothing in it and the lookup fails.
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
    /// its value) and after ([`TaxLookup::Unnamed`]: `info.taxid` is set to `tax_id`, the names
    /// stay empty), which blastdbcmd's `-tax_info` writes apart (`N/A` or empty names), and giving
    /// a record past NCBI's mapping of `taxdb.btd` as it is ([`TaxLookup::PastMapping`],
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
        // `m_End` minus `m_Begin` cast to `int`: a pointer difference that is the wrapped
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
                    // `lookup`: `m_End` minus `m_Begin` cast to `int`, the wrapped `Uint4`
                    // difference.
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:161
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:199
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
// NCBI maps each file with `CMemoryFile`, whose exception leaves the constructor of the function's
// `static` taxonomy file object: the object stays unbuilt, so every later lookup builds it again
// and throws again (NCBI BLAST+ 2.17.0 with `taxdb.bti` or `taxdb.btd` of mode 000: every load of a
// subject with a taxid fails; line WQA D18b probe 12, 2026-11-21).
/// The error of a file of the taxonomy database that does not open or whose first byte cannot be
/// read: [`SeqDbError::TaxdbThrows`] where NCBI's `CMemoryFile` throws there (a system text), the
/// error as it is otherwise.
fn throws_at_every_lookup(error: SeqDbError) -> SeqDbError {
    let Some(exception) = error.memory_file_exception() else {
        return error;
    };
    match error {
        SeqDbError::Open { path, source } | SeqDbError::Map { path, source, .. } => {
            SeqDbError::TaxdbThrows {
                path,
                reason: crate::utils::os_error_text::ncbi_strerror(&source)
                    .map_or_else(|| source.to_string(), |text| text.into_owned()),
                exception: Box::new(exception),
            }
        }
        other => other,
    }
}

fn be32(bytes: &[u8]) -> u32 {
    u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

/// Why the taxonomy database of a run did not open, given again at every later call
/// ([`LazyTaxDb`]).
#[derive(Debug)]
enum OpenFailure {
    /// A file that NCBI's `CMemoryFile` throws on at every lookup ([`SeqDbError::TaxdbThrows`]).
    Throws {
        path: PathBuf,
        reason: String,
        exception: Box<crate::ncbi_exception::NcbiException>,
    },
    /// LOSAT's refusal of the database ([`SeqDbError::Unsupported`]): the file and the reason.
    Refused(PathBuf, String),
}

/// The taxonomy database of one run, opened at the first name asked for and
/// not before (NCBI: the `static CTaxDBFileInfo` of `GetTaxNames`). A run that
/// never asks for names touches no file.
#[derive(Debug)]
pub struct LazyTaxDb {
    search: SearchPath,
    db: OnceLock<Result<TaxDb, OpenFailure>>,
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
                SeqDbError::TaxdbThrows {
                    path,
                    reason,
                    exception,
                } => OpenFailure::Throws {
                    path,
                    reason,
                    exception,
                },
                SeqDbError::Unsupported { path, message } => OpenFailure::Refused(path, message),
                // The file that failed, and what the system said of it.
                SeqDbError::Open { path, source }
                | SeqDbError::Io { path, source }
                | SeqDbError::Map { path, source, .. }
                | SeqDbError::Read { path, source } => {
                    OpenFailure::Refused(path, source.to_string())
                }
                SeqDbError::Format { path, message } => OpenFailure::Refused(path, message),
                other => OpenFailure::Refused(PathBuf::from("taxdb.bti"), other.to_string()),
            })
        });
        match db {
            Ok(db) => {
                if opened || self.handed_again() {
                    messages.extend(db.open_messages().iter().cloned());
                }
                Ok(db)
            }
            Err(OpenFailure::Refused(path, message)) => Err(SeqDbError::Unsupported {
                path: path.clone(),
                message: message.clone(),
            }),
            Err(OpenFailure::Throws {
                path,
                reason,
                exception,
            }) => Err(SeqDbError::TaxdbThrows {
                path: path.clone(),
                reason: reason.clone(),
                exception: exception.clone(),
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

/// The size of `CSeqDBVol`'s `m_TaxCache`, 256 slots (seqdbvol.cpp:84).
const TAX_CACHE_SLOTS: usize = 256;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1105-1107, 1140-1182
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbgeneral.hpp:1184-1196
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
#[path = "taxdb_tests.rs"]
mod tests;
