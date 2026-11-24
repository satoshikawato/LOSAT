//! The string ISAM files of a version 4 database volume (`.psi`/`.psd` for protein,
//! `.nsi`/`.nsd` for nucleotide): the string index of the Seq-ids that `makeblastdb
//! -parse_seqids` writes (`crate::writedb::isam` is the writer of the same files), read-only.
//!
//! NCBI's reader is `CSeqDBIsam` (`c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp` and
//! `impl/seqdbisam.hpp`) with the string parts of `CSeqDBVol` (`seqdbvol.cpp`) and the
//! simplification of a user's id (`SeqDB_SimplifyAccession`, `seqdbcommon.cpp`). NCBI maps the
//! files; here the index file is read by position (`super::file::DbFile`) where a lookup needs
//! it and the data file one page at a time (the batch calls read the whole index file once, as
//! NCBI's `x_LoadIndex` loads every sample key, and the negative scan reads the data file in
//! chunks of whole pages). Pure `std`.
//!
//! # File layout (big-endian `Int4`s)
//!
//! Index file: nine fields (version `1`, type `2` = `eString`, length of the data file, number
//! of terms, number of samples, page size (`64`), longest line (`4096`), sparse flag, `0`), then
//! `samples + 1` offsets of the pages in the data file, then `samples + 1` offsets of the sample
//! keys in the index file, then the sample keys: the first record of every page less its line
//! feed (`key \x02 oid`) and a NUL. Data file: the records `key \x02 <decimal oid> \n`, sorted,
//! one per distinct (key, OID); the keys are lower-cased by the writer.
//!
//! # NCBI functions and what they became
//!
//! - `CSeqDBIsam` constructor, `x_InitSearch`, `x_FindIndexBounds`: [`StringIsam::open`] and
//!   [`StringIsam::from_files`]; `x_OutOfBounds` of a string: the private `out_of_bounds`;
//!   `GetIdBounds` of a string: `StringIsam::id_bounds`, which only the tests call.
//! - `x_StringSearch` with `x_DiffSample`, `x_DiffCharLease`, `x_DiffChar`, `x_LoadPage`,
//!   `x_ExtractAllData`, `x_ExtractPageData`, `x_ExtractData`: [`StringIsam::string_search_into`]
//!   (every line of the volume that matches a term, as `(key, data)` strings).
//! - `StringToOids`: [`StringIsam::string_to_oids`]; `CSeqDBVol::x_CheckVersions`:
//!   [`check_versions`]; the string case of `CSeqDBVol::x_StringToOids`:
//!   [`StringIsam::string_id_to_oids`].
//! - `IdsToOids` of a `CSeqDBGiList` (`x_TranslateGiList<string>`): [`StringIsam::ids_to_oids`];
//!   `CSeqDBGiList::InsureOrder` for the string list, with the list's order flag: [`insure_order`];
//!   `SiToOid`: `filter::find_si`.
//! - `IdsToOids` of a `CSeqDBNegativeList` (`x_SearchNegativeMultiSeq`):
//!   [`StringIsam::negative_scan`]; `CSeqDBNegativeList` and `CSeqDBBitVector`: [`NegativeOids`].
//! - `SeqDB_SimplifyAccession`, `SeqDB_SimplifySeqid`, `PreprocessIdsForISAMSiLookup`:
//!   [`seqdb_simplify_accession`], [`simplify_accession_for_isam`].
//!
//! The numeric ISAM of GIs (`.pni`, `.nni`) is [`super::isam_num`]. Not ported, because nothing
//! reaches it: the other numeric ISAMs, `x_GetIndexKeyOffset`
//! and `x_GetIndexString` (used only behind `short_match` and `follow_match`, two constants that
//! are `false`), `x_SparseStringToOids` (`_TROUBLE`), `HashToOids`, and the numeric
//! `x_SearchNegativeMulti`.
//!
//! # `CSeqDBNegativeList`: what the OID bit vectors mean
//!
//! A negative id list is applied to a volume by [`StringIsam::negative_scan`]. Every line of the
//! ISAM is read. A line whose key is in the sorted id list marks its OID *visible*
//! (`AddVisibleOid`); any other line marks it *included* (`AddIncludedOid`). `GetOidStatus` of
//! `oid` is true when the OID is included or is not visible: an OID is kept unless it was seen only
//! through listed keys, so a sequence with one listed key and one unlisted key (an accession and
//! its `.version` are two lines of one OID) stays. `GetNumOids()` is the larger of the two vectors'
//! sizes (the highest OID set plus one); `CSeqDBOIDList::x_ApplyNegativeList`
//! (`seqdboidlist.cpp:376-404`) clears every OID of the volume from `GetNumOids()` on and every OID
//! below it whose status is false. Both vectors are `CSeqDBBitVector`s: `SetBit` grows `Size` to
//! one more than `oid`, `GetBit` is false past `Size`.
//!
//! # Quirks of NCBI that are kept
//!
//! - `x_DiffChar` (`seqdbisam.cpp:589-632`) compares ignoring ASCII case and takes NUL, CR and LF
//!   of either side as equal (they end a line), and after the match it skips *spaces* of the
//!   file's key: a term `2xyy` finds the lines `2xyy` and `2xyy  ` (the writer adds `pdb|` keys
//!   whose chain is a blank). A term with trailing spaces matches only keys with at least those
//!   spaces. The scan of a page stops at the first non-matching line after a match.
//! - The sample search (`x_StringSearch`) reads `term[diff]` where `diff` may equal the term's
//!   length: NCBI reads the `std::string`'s terminating NUL; the Rust code reads `0`.
//! - `x_DiffCharLease` reads `term + 1` bytes of the index file past the sample key's start; the
//!   bytes past the end of the file read as NUL (a mapping's zero padding). When the
//!   sample key lies so near the end of the index file that the read would pass it, NCBI returns
//!   `file_length - offset` as the first difference even when the bytes match (`result`); that is
//!   kept.
//! - `StringToOids` looks for `gb|acc|` and `gb||acc` before `acc` (unless the id was simplified),
//!   strips a numeric `.version` of 1 to 3 digits only when nothing matched, and last tries
//!   `AsFastaString()` of a `CSeq_id` made from `acc` with the flags `fParse_RawText` and
//!   `fParse_AnyLocal`; the `OID`s come from `atoi` of the data strings. [`check_versions`] then
//!   keeps only OIDs whose Seq-ids hold the accession (the text after the last `|` of `acc` before
//!   its first `.`) with the version the text after that `.` starts with.
//! - The positive batch (`x_TranslateGiList<string>`) is not "first ISAM line wins". For a key
//!   that has several lines in one page, which line sets the OID of the equal list entries
//!   depends on the path: an entry reached by `s_AdvanceKeyList`, which stops at the *last* line
//!   with a key `<=` the id, takes that last line; an entry reached right after the previous
//!   line was used (the step of `index`, `s_AdvanceGiList`, `s_SetTranslation` with `keys[index]`)
//!   takes the line at `index`, the *first* line of the key when the previous list entry was a
//!   smaller key. Entries equal to a translated key later in the list are translated with it
//!   (`s_SetTranslation`, which overwrites an OID that was set), and entries that have an OID
//!   before the call are skipped by `s_AdvanceGiList` when it passes them.
//! - The negative scan (`x_SearchNegativeMultiSeq`) tests `vals[i] < vol_end` with the
//!   volume-local OID against the *global* end of the volume (`vol_end` is a parameter for that
//!   reason). `s_IsSameAccession`: when a key has no `.` and the *next* line of the same page
//!   is the same text plus a `.version` (`acc` found anywhere in it, split at its first `.`),
//!   the key's line is not marked included, and when it is listed the next line is skipped
//!   without being marked at all (it is assumed to be of the same OID).
//! - Keys are compared as `std::string`s, that is by unsigned bytes then length. The writer sorts
//!   the first six bytes of a record as signed `char`s; the two orders differ for keys with
//!   bytes of `0x80` or more, which Seq-ids do not hold.
//!
//! # Where LOSAT cannot or does not reproduce NCBI
//!
//! - A damaged file: NCBI reads past its mapping (undefined behaviour); every such read is an
//!   `Err(SeqDbError::Format)` here, and so are an index that holds no sample, a page size or
//!   count that does not fit the sample count, offsets outside the files, and an index
//!   shorter than its offset tables.
//! - `s_AdvanceKeyList` steps to `-1` where the walk of a list meets a key below the first
//!   sample key or below the first key of the page it loads: a list out of order (NCBI walks a
//!   list in the order its flag leaves it, [`insure_order`]) or a data file out of order. NCBI
//!   then reads the entry before the start of its vector of offsets or keys (undefined
//!   behaviour): LOSAT's refusal ([`walk_before_start`]).
//! - A volume whose files fail `x_InitSearch` (wrong version, type or data length) is opened as
//!   NCBI does, uninitialised: lookups find nothing, and the batch calls fail with NCBI's
//!   "Error: Unable to use ISAM index in batch mode.". NCBI's `x_FindIndexBounds` still reads the
//!   index of such a volume; its result is never used and is not computed.
//! - The warnings that `CSeq_id` posts while parsing a user's id (`ERR_POST` to the diagnostics
//!   stream) are dropped.
//! - Exceptions of `SeqDB_SimplifyAccession` and `SeqDB_SimplifySeqid` that are not a
//!   `CSeqIdException` (a `ti` tag that is not a number, a string tag of `BL_ORD_ID` or `PIG`)
//!   are [`SimplifyError::Exception`] with NCBI's exception (`ncbi_str::string_to_int8`'s, or
//!   `CObject_id::GetId()`'s `CInvalidChoiceSelection`, thrown in code that NCBI's `datatool`
//!   generates: `ThrowSite::generated`).
//! - The functions assume the C locale (`tolower`, `toupper`, `isdigit`, `atoi`).

use std::path::{Path, PathBuf};

use super::file::DbFile;
use super::SeqDbError;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::{string_to_int, string_to_int8, string_to_uint, ALLOW_TRAILING_SYMBOLS};
use crate::seqid::accguide::{parse_flags, SeqIdWarning, SEQID_NOT_SET};
use crate::seqid::label;
use crate::seqid::model::{ObjectId, ParsedSeqId, TextseqId};
use crate::seqid::parse::{parse_ids, seq_id_from_text, which_inverse_seq_id, SeqIdError};
use crate::seqid::rank::{best_rank_score, find_best};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:45-55
/// `ISAM_VERSION`.
const ISAM_VERSION: i32 = 1;
/// `CSeqDBIsam::eString`.
const ISAM_STRING_TYPE: i32 = 2;
/// `MEMORY_ONLY_PAGE_SIZE`.
const MEMORY_ONLY_PAGE_SIZE: i32 = 1;
/// `m_KeySampleOffset` (`x_InitSearch`): the nine header fields.
const KEY_SAMPLE_OFFSET: u64 = 9 * 4;
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:580-581
/// `ISAM_DATA_CHAR`: the byte that ends a key when data follows.
const ISAM_DATA_CHAR: u8 = 2;
/// The bytes of the data file the negative scan reads at once (whole pages; LOSAT: NCBI maps
/// the file).
const SCAN_CHUNK: u64 = 4 << 20;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:560-578
/// `s_SeqDBIsam_NullifyEOLs`.
fn nullify_eols(c: u8) -> u8 {
    if c == b'\n' || c == b'\r' {
        0
    } else {
        c
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:583-587
/// `ENDS_ISAM_KEY`.
fn ends_isam_key(c: u8) -> bool {
    c == ISAM_DATA_CHAR || nullify_eols(c) == 0
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:589-632
/// `CSeqDBIsam::x_DiffChar`: `begin` is `[begin, end)`; the index of the first difference, or -1
/// for a complete match.
fn diff_char(term_in: &[u8], begin: &[u8], ignore_case: bool) -> i64 {
    let bytes = begin.len();
    let mut i = 0;
    while i < bytes && i < term_in.len() {
        let mut ch1 = term_in[i];
        let mut ch2 = begin[i];
        if ch1 != ch2 {
            ch1 = nullify_eols(ch1);
            ch2 = nullify_eols(ch2);
            if ignore_case {
                ch1 = ch1.to_ascii_uppercase();
                ch2 = ch2.to_ascii_uppercase();
            }
            if ch1 != ch2 {
                break;
            }
        }
        i += 1;
    }
    let mut p = i;
    while p < bytes && begin[p] == b' ' {
        p += 1;
    }
    if (p == bytes || ends_isam_key(begin[p])) && i == term_in.len() {
        -1
    } else {
        i as i64
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:793-820
/// `CSeqDBIsam::x_ExtractData`: the key and the data of the line at the start of `line` (which
/// ends where the page does); nothing when the page ends inside the line.
fn extract_data(line: &[u8]) -> Option<IsamLine> {
    let mut data_ptr: Option<usize> = None;
    for (p, &c) in line.iter().enumerate() {
        match nullify_eols(c) {
            0 => {
                return Some(match data_ptr {
                    Some(d) => IsamLine {
                        key: line[..d].to_vec(),
                        data: line[d + 1..p].to_vec(),
                    },
                    None => IsamLine {
                        key: line[..p].to_vec(),
                        data: Vec::new(),
                    },
                });
            }
            ISAM_DATA_CHAR => data_ptr = Some(p),
            _ => {}
        }
    }
    None
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:634-686
/// `CSeqDBIsam::x_ExtractPageData`: every line of `page` that matches `term_in`, up to the
/// first line that does not once one matched.
fn extract_page_data(term_in: &[u8], page: &[u8], lines: &mut Vec<IsamLine>) {
    let ignore_case = true;
    let endp = page.len();
    let mut indexp = 0;
    let mut found_match = false;
    while indexp < endp {
        let diff = diff_char(term_in, &page[indexp..], ignore_case);
        if diff == -1 {
            found_match = true;
            if let Some(line) = extract_data(&page[indexp..]) {
                lines.push(line);
            }
        } else if found_match {
            break;
        }
        // Skip remainder of term, and any nulls after it.
        while indexp < endp && nullify_eols(page[indexp]) != 0 {
            indexp += 1;
        }
        while indexp < endp && nullify_eols(page[indexp]) == 0 {
            indexp += 1;
        }
    }
}

/// One line of the data file as `x_ExtractData` gives it: the key and the text after its
/// `\x02` (the decimal OID; empty for a line without one).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsamLine {
    pub key: Vec<u8>,
    pub data: Vec<u8>,
}

/// `x_StringSearch`'s `EErrorCode`: `eNoError`, `eNotFound` and the `eInitFailed` of an
/// uninitialised index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SearchResult {
    Found,
    NotFound,
    InitFailed,
}

/// The fields of the index header that `x_InitSearch` reads.
struct Header {
    num_terms: i32,
    num_samples: i32,
    page_size: i32,
    max_line_size: i32,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1201-1209
// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:4038-4042 and 3941-3944
/// `CSeqDBIsam::IndexExists` of the index file `index` and the data file `data` of an ISAM (the
/// string, PIG and trace-id ISAMs of a version 4 volume, `x_OpenStrFile`, `x_OpenPigFile`,
/// `x_OpenTiFile`): both are regular files, or symbolic links to one (`CFile::Exists` is
/// `IsFile`, [`super::alias::exists`]). A directory in their place is no ISAM: the volume then
/// has none, and its lookups find nothing there (review 4 of line DBID loadord, H1).
pub(super) fn index_exists(index: &Path, data: &Path) -> bool {
    super::alias::exists(index) && super::alias::exists(data)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1150-1160
// (cited at `StringIsam::open`)
/// The open of the two files of an ISAM by the `CSeqDBIsam` constructor: NCBI's `Error: Could
/// not open input file (<index>/<data>)` where one of them does not open
/// ([`DbFile::open`]: no regular file, or one that `open` refuses).
pub(super) fn open_files(index: &Path, data: &Path) -> Result<(DbFile, DbFile), SeqDbError> {
    let message = || {
        SeqDbError::NotFound(
            format!(
                "Error: Could not open input file ({}/{})",
                index.display(),
                data.display()
            )
            .into(),
        )
    };
    let index_file = DbFile::open(index).map_err(|_| message())?;
    let data_file = DbFile::open(data).map_err(|_| message())?;
    Ok((index_file, data_file))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:190-198
/// The files of the ISAM of a volume with sequences as `x_OpenStrFile` and `x_OpenPigFile` make
/// it (`index`, `data`): `None` where the volume has none ([`index_exists`]), the constructor's
/// exception where one of them does not open ([`open_files`]), else both open. The one owner of
/// "this volume has an ISAM and it opens", which the lookups (`StringIsam::open_existing`,
/// `NumericIsam::open_pig_volume`) and the search's question whether a lookup can stop on it
/// (`SeqDb::isams_may_fail_to_open`) both ask; neither reads a byte here.
pub(super) fn open_if_exists(
    index: &Path,
    data: &Path,
) -> Result<Option<(DbFile, DbFile)>, SeqDbError> {
    if !index_exists(index, data) {
        return Ok(None);
    }
    open_files(index, data).map(Some)
}

/// A string ISAM index of one volume (`CSeqDBIsam` of `eStringId`).
#[derive(Debug)]
pub struct StringIsam {
    index: DbFile,
    data: DbFile,
    initialized: bool,
    num_terms: i32,
    num_samples: i32,
    page_size: i32,
    max_line_size: i32,
    /// `m_FirstKey` and `m_LastKey` (lower-cased), set when the index is initialised.
    first_key: Option<Vec<u8>>,
    last_key: Option<Vec<u8>>,
    /// The bytes of the index file, where the ISAM keeps them (`keep_index`).
    index_bytes: Option<Vec<u8>>,
}

impl StringIsam {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1173-1210
    /// `CSeqDBIsam::x_MakeFilenames` (`file_ext_char` is `s`): `<dbname>.ps?` or `.ns?`, the index
    /// file ending in `i` and the data file in `d`. The volume's path `dbname` is the bytes of
    /// NCBI's `std::string`.
    pub fn file_names(dbname: &Path, protein: bool) -> Result<(PathBuf, PathBuf), SeqDbError> {
        if dbname.as_os_str().is_empty() {
            return Err(SeqDbError::NotFound(
                "Error: argument not valid".to_string().into(),
            ));
        }
        let mut index = dbname.as_os_str().to_os_string();
        index.push(if protein { ".ps" } else { ".ns" });
        let mut data = index.clone();
        index.push("i");
        data.push("d");
        Ok((PathBuf::from(index), PathBuf::from(data)))
    }

    /// `CSeqDBIsam::IndexExists` for the string index ([`index_exists`]: two regular files).
    pub fn index_exists(dbname: &Path, protein: bool) -> bool {
        match Self::file_names(dbname, protein) {
            Ok((index, data)) => index_exists(&index, &data),
            Err(_) => false,
        }
    }

    /// The string index of the volume `dbname` (its path without extension).
    pub fn open_volume(dbname: &Path, protein: bool) -> Result<Self, SeqDbError> {
        let (index, data) = Self::file_names(dbname, protein)?;
        Self::open(&index, &data)
    }

    /// `x_OpenStrFile`: the string index of the volume `dbname` (its path without extension, a
    /// volume with sequences), `None` where it has none ([`open_if_exists`]).
    pub fn open_existing(dbname: &Path, protein: bool) -> Result<Option<Self>, SeqDbError> {
        let (index, data) = Self::file_names(dbname, protein)?;
        open_if_exists(&index, &data)?
            .map(|(index, data)| Self::from_files(index, data))
            .transpose()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1102-1171
    /// The constructor: opens the two files. A file that cannot be opened is NCBI's "Error: Could
    /// not open input file (`<index>/<data>`)".
    pub fn open(index_path: &Path, data_path: &Path) -> Result<Self, SeqDbError> {
        let (index, data) = open_files(index_path, data_path)?;
        Self::from_files(index, data)
    }

    /// The constructor over two open files: `x_InitSearch`, then `x_FindIndexBounds`.
    pub fn from_files(index: DbFile, data: DbFile) -> Result<Self, SeqDbError> {
        let mut isam = Self {
            index,
            data,
            initialized: false,
            num_terms: 0,
            num_samples: 0,
            page_size: 0,
            max_line_size: 0,
            first_key: None,
            last_key: None,
            index_bytes: None,
        };
        if let Some(header) = isam.init_search()? {
            isam.num_terms = header.num_terms;
            isam.num_samples = header.num_samples;
            isam.page_size = header.page_size;
            isam.max_line_size = header.max_line_size;
            isam.initialized = true;
            isam.check_header()?;
            isam.find_index_bounds()?;
        }
        Ok(isam)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:58-121
    /// `x_InitSearch`: `None` for NCBI's `eWrongFile`, `eBadVersion` and `eBadType`, which leave
    /// the index uninitialised.
    fn init_search(&self) -> Result<Option<Header>, SeqDbError> {
        let info_needed = 10 * 4;
        if self.index.len() < info_needed {
            return Ok(None);
        }
        let mut raw = [0u8; 36];
        self.index.read_exact_at(&mut raw, 0)?;
        let field = |i: usize| {
            i32::from_be_bytes([raw[4 * i], raw[4 * i + 1], raw[4 * i + 2], raw[4 * i + 3]])
        };
        if field(0) != ISAM_VERSION {
            return Ok(None);
        }
        if field(1) != ISAM_STRING_TYPE {
            return Ok(None);
        }
        let header = Header {
            num_terms: field(3),
            num_samples: field(4),
            page_size: field(5),
            max_line_size: field(6),
        };
        if header.page_size != MEMORY_ONLY_PAGE_SIZE {
            // Special case of memory-only index
            let data_file_length = i64::from(field(2));
            if i64::try_from(self.data.len()).ok() != Some(data_file_length) {
                return Ok(None);
            }
        }
        Ok(Some(header))
    }

    /// LOSAT: the checks that keep a damaged header from being used (NCBI reads whatever the
    /// numbers point at).
    fn check_header(&self) -> Result<(), SeqDbError> {
        let bad = |message: String| SeqDbError::Format {
            path: self.index.path().to_path_buf(),
            message,
        };
        let (terms, samples, page) = (
            i64::from(self.num_terms),
            i64::from(self.num_samples),
            i64::from(self.page_size),
        );
        if samples < 1 || page < 1 || self.max_line_size < 0 || terms < 1 {
            return Err(bad(format!(
                "string ISAM header: {terms} terms, {samples} samples, page size {page}, longest \
                 line {}",
                self.max_line_size
            )));
        }
        if (samples - 1) * page >= terms || terms > samples * page {
            return Err(bad(format!(
                "string ISAM header: {terms} terms do not fill {samples} samples of {page}"
            )));
        }
        let tables = if self.page_size == MEMORY_ONLY_PAGE_SIZE {
            1
        } else {
            2
        };
        let needed = KEY_SAMPLE_OFFSET + 4 * (samples as u64 + 1) * tables;
        if self.index.len() < needed {
            return Err(bad(format!(
                "string ISAM index of {} bytes is shorter than its offset tables ({needed} bytes)",
                self.index.len()
            )));
        }
        Ok(())
    }

    /// Whether `x_InitSearch` accepted the files (`m_Initialized`).
    pub fn is_initialized(&self) -> bool {
        self.initialized
    }

    /// The number of lines (`m_NumTerms`).
    #[cfg(test)]
    pub fn num_terms(&self) -> i32 {
        self.num_terms
    }

    /// The number of pages (`m_NumSamples`).
    #[cfg(test)]
    pub fn num_samples(&self) -> i32 {
        self.num_samples
    }

    /// The lines per page (`m_PageSize`).
    #[cfg(test)]
    pub fn page_size(&self) -> i32 {
        self.page_size
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1646-1665
    /// `GetIdBounds` with two strings and an int, all by reference: the lower-cased first and last
    /// key and the number of lines; count 0 and empty keys for an uninitialised index.
    #[cfg(test)]
    pub fn id_bounds(&self) -> (Vec<u8>, Vec<u8>, i32) {
        if !self.initialized {
            return (Vec::new(), Vec::new(), 0);
        }
        match (&self.first_key, &self.last_key) {
            (Some(low), Some(high)) => (low.clone(), high.clone(), self.num_terms),
            _ => (Vec::new(), Vec::new(), 0),
        }
    }

    /// The four bytes at `at` of the index file.
    fn index_u32(&self, at: u64) -> Result<u32, SeqDbError> {
        let mut bytes = [0u8; 4];
        self.read_index(&mut bytes, at)?;
        Ok(u32::from_be_bytes(bytes))
    }

    /// The ISAM with the bytes of its index file kept in memory, which the binary search over
    /// the samples then reads (NCBI maps the file): for the search programs, which look up
    /// every reported sequence. No NCBI counterpart beyond the mapping.
    pub fn keep_index(mut self) -> Result<Self, SeqDbError> {
        self.index_bytes = Some(self.index.read_all()?);
        Ok(self)
    }

    /// `buf` from `at` of the index file (`DbFile::read_exact_at`, from the kept bytes where
    /// there are some, with the same error past the end).
    fn read_index(&self, buf: &mut [u8], at: u64) -> Result<(), SeqDbError> {
        let Some(bytes) = &self.index_bytes else {
            return self.index.read_exact_at(buf, at);
        };
        let range = usize::try_from(at)
            .ok()
            .and_then(|start| Some(start..start.checked_add(buf.len())?))
            .filter(|range| range.end <= bytes.len());
        match range {
            Some(range) => {
                buf.copy_from_slice(&bytes[range]);
                Ok(())
            }
            None => self.index.read_exact_at(buf, at),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:899-923
    /// `x_LoadPage` with `SampleNum1` and `SampleNum2`: the data file's bytes from the first record
    /// of page `s1` to the first of page `s2`.
    fn load_page(&self, s1: i64, s2: i64) -> Result<Vec<u8>, SeqDbError> {
        let begin_offset = KEY_SAMPLE_OFFSET + 4 * s1 as u64;
        let key_off1 = self.index_u32(begin_offset)?;
        let key_off2 = self.index_u32(begin_offset + 4 * (s2 - s1) as u64)?;
        if key_off2 < key_off1 {
            return Err(SeqDbError::Format {
                path: self.index.path().to_path_buf(),
                message: format!(
                    "string ISAM page offsets {key_off1} and {key_off2} run backwards"
                ),
            });
        }
        self.data
            .read_range(u64::from(key_off1), u64::from(key_off2))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1461-1583
    /// `x_FindIndexBounds` (the string branch): the lower-cased first key of the first page and
    /// last key of the last page.
    fn find_index_bounds(&mut self) -> Result<(), SeqDbError> {
        let corrupt = |isam: &Self, what: &str| SeqDbError::Format {
            path: isam.data.path().to_path_buf(),
            message: format!("string ISAM {what} page holds no complete line"),
        };
        let start = 0;
        let stop = i64::from(self.num_samples) - 1;
        // Load the first page, get the first term
        let page = self.load_page(start, start + 1)?;
        let first = extract_data(&page).ok_or_else(|| corrupt(self, "first"))?;
        // Load the last page, advance to the last item
        let page = self.load_page(stop, stop + 1)?;
        let endp = page.len();
        let mut lastp = None;
        let mut indexp = 0;
        while indexp < endp {
            // Remember our new "last term" value.
            lastp = Some(indexp);
            // Skip remainder of term, and any nulls after it.
            while indexp < endp && nullify_eols(page[indexp]) != 0 {
                indexp += 1;
            }
            while indexp < endp && nullify_eols(page[indexp]) == 0 {
                indexp += 1;
            }
        }
        let last = lastp
            .and_then(|at| extract_data(&page[at..]))
            .ok_or_else(|| corrupt(self, "last"))?;
        self.first_key = Some(first.key.to_ascii_lowercase());
        self.last_key = Some(last.key.to_ascii_lowercase());
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1604-1623
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1142-1148
    /// `x_OutOfBounds` of a string: the lower-cased `key` is before the first key or after the last
    /// one.
    fn out_of_bounds(&self, key: &[u8]) -> bool {
        let (Some(first), Some(last)) = (&self.first_key, &self.last_key) else {
            return false;
        };
        let key = key.to_ascii_lowercase();
        key < *first || key > *last
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:516-558
    /// `x_DiffCharLease` over the index file: the first difference of `term_in` and the bytes at
    /// `key_offset` (`x_DiffChar` over `term + 1` bytes, NUL past the end of the file), and
    /// those bytes (`KeyData`).
    fn diff_char_lease(
        &self,
        term_in: &[u8],
        at_least: u32,
        key_offset: u64,
    ) -> Result<(i64, Vec<u8>), SeqDbError> {
        let file_length = self.index.len() as i64;
        let offset_begin = key_offset as i64;
        // Add one to term_end to insure we don't consider "AA" and "AAB" as equal.
        let term_end = offset_begin + term_in.len() as i64 + 1;
        let map_end = term_end + i64::from(at_least);
        let mut result = -1;
        if map_end > file_length && term_end > file_length {
            result = file_length - offset_begin;
        }
        if offset_begin > file_length {
            return Err(SeqDbError::Format {
                path: self.index.path().to_path_buf(),
                message: format!(
                    "string ISAM key offset {offset_begin} is past the end of the file"
                ),
            });
        }
        let mut window = vec![0u8; term_in.len() + 1];
        let readable = window.len().min((file_length - offset_begin) as usize);
        self.read_index(&mut window[..readable], key_offset)?;
        let dc_result = diff_char(term_in, &window, true);
        if dc_result != -1 {
            return Ok((dc_result, window));
        }
        Ok((result, window))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:863-897
    /// `x_DiffSample`: the first difference of `term_in` and the key of sample `sample_num`,
    /// and the bytes of that key.
    fn diff_sample(&self, term_in: &[u8], sample_num: i64) -> Result<(i64, Vec<u8>), SeqDbError> {
        let mut sample_offset = KEY_SAMPLE_OFFSET;
        if self.page_size != MEMORY_ONLY_PAGE_SIZE {
            sample_offset += (self.num_samples as u64 + 1) * 4;
        }
        let offset_begin = sample_offset + 4 * sample_num as u64;
        let key_offset = u64::from(self.index_u32(offset_begin)?);
        let max_lines_2 = (self.max_line_size as u32).wrapping_mul(2);
        self.diff_char_lease(term_in, max_lines_2, key_offset)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:688-791
    /// `x_ExtractAllData`: the page of `sample_index`, known to hold a match, with the pages
    /// around it for as long as their first (before) and last (after) line match too.
    fn extract_all_data(
        &self,
        term_in: &[u8],
        sample_index: i64,
        lines: &mut Vec<IsamLine>,
    ) -> Result<(), SeqDbError> {
        let ignore_case = true;
        let num_samples = i64::from(self.num_samples);
        let mut pre_amt = 1;
        let mut post_amt = 1;
        let (mut done_b, mut done_e) = (false, false);
        let mut page = Vec::new();
        while !(done_b && done_e) {
            let beg_off = if sample_index < pre_amt {
                done_b = true;
                0
            } else {
                sample_index - pre_amt
            };
            let end_off = if (num_samples - sample_index) < post_amt {
                done_e = true;
                num_samples
            } else {
                sample_index + post_amt
            };
            page = self.load_page(beg_off, end_off)?;
            if !done_b {
                if diff_char(term_in, &page, ignore_case) != -1 {
                    done_b = true;
                } else {
                    pre_amt += 1;
                }
            }
            if !done_e {
                let mut last_term = None;
                let mut p = page.len() as i64 - 1;
                // Skip over any non-terminating junk at the end
                let mut end_nulls = true;
                while p > 0 {
                    let terminal = nullify_eols(page[p as usize]) == 0;
                    if end_nulls {
                        if !terminal {
                            end_nulls = false;
                        }
                    } else if terminal {
                        last_term = Some(p as usize + 1);
                        break;
                    }
                    p -= 1;
                }
                let last_term = last_term.unwrap_or(0);
                if diff_char(term_in, &page[last_term..], ignore_case) != -1 {
                    done_e = true;
                } else {
                    post_amt += 1;
                }
            }
        }
        extract_page_data(term_in, &page, lines);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:933-1100
    /// `x_StringSearch`: appends the lines that match `term_in` to `lines`.
    fn string_search_into(
        &self,
        term_in: &[u8],
        lines: &mut Vec<IsamLine>,
    ) -> Result<SearchResult, SeqDbError> {
        let preexisting_data_count = lines.len();
        if !self.initialized {
            return Ok(SearchResult::InitFailed);
        }
        if self.out_of_bounds(term_in) {
            return Ok(SearchResult::NotFound);
        }
        // search the sample file first
        let mut start: i64 = 0;
        let mut stop: i64 = i64::from(self.num_samples) - 1;
        let mut sample_num: i64 = -1;
        while stop >= start {
            sample_num = i64::from(((stop + start) as u32) >> 1);
            let (diff, key_data) = self.diff_sample(term_in, sample_num)?;
            // If this is an exact match, return the master term number.
            if diff == -1 {
                self.extract_all_data(term_in, sample_num, lines)?;
                return Ok(SearchResult::Found);
            }
            // Otherwise, search for the next sample. `term_in[diff]` of a `std::string` is its
            // NUL when `diff` is the length.
            let diff = diff as usize;
            let term_byte = term_in.get(diff).copied().unwrap_or(0);
            let key_byte = key_data.get(diff).copied().unwrap_or(0);
            if term_byte.to_ascii_lowercase() < key_byte.to_ascii_lowercase() {
                sample_num -= 1;
                stop = sample_num;
            } else {
                start = sample_num + 1;
            }
        }
        // If the term is out of range altogether, report not finding it.
        if sample_num < 0 || sample_num >= i64::from(self.num_samples) {
            return Ok(SearchResult::NotFound);
        }
        // Load the appropriate page of terms into memory.
        let page = self.load_page(sample_num, sample_num + 1)?;
        // Search the page for the term.
        extract_page_data(term_in, &page, lines);
        // For now the short and follow logic is not implemented.
        if preexisting_data_count == lines.len() {
            Ok(SearchResult::NotFound)
        } else {
            Ok(SearchResult::Found)
        }
    }

    /// `x_StringSearch` for one term: every line of the volume that matches it, ignoring ASCII
    /// case, in file order (see the quirks of `x_DiffChar` above). Empty for an uninitialised
    /// index.
    #[cfg(test)]
    pub fn string_search(&self, term_in: &[u8]) -> Result<Vec<IsamLine>, SeqDbError> {
        let mut lines = Vec::new();
        self.string_search_into(term_in, &mut lines)?;
        Ok(lines)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1236-1376
    /// `StringToOids`: the OIDs (local to the volume) of the lines found for `acc`. `adjusted`:
    /// `acc` was simplified (the `simpler` of `SeqDB_SimplifyAccession`), so `gb|acc|` and
    /// `gb||acc` are not tried; `strip_version`: the in/out `version_check` of NCBI, true when the
    /// caller will call [`check_versions`] if the result says so.
    pub fn string_to_oids(
        &self,
        acc: &[u8],
        adjusted: bool,
        strip_version: bool,
    ) -> Result<StringOids, SeqDbError> {
        let mut version_check = false;
        let mut result = StringOids {
            oids: Vec::new(),
            version_check,
        };
        if !self.initialized {
            return Ok(result);
        }
        let mut found = false;
        let mut accession = b"gb|".to_vec();
        accession.extend_from_slice(acc);
        accession.push(b'|');
        let mut locus_str = b"gb||".to_vec();
        locus_str.extend_from_slice(acc);
        let mut err = SearchResult::Found;
        let mut lines = Vec::new();
        if !adjusted {
            err = self.string_search_into(&accession, &mut lines)?;
            if err == SearchResult::Found {
                found = true;
            }
            if !found {
                err = self.string_search_into(&locus_str, &mut lines)?;
            }
            if err != SearchResult::NotFound {
                found = true;
            }
        }
        if !found {
            err = self.string_search_into(acc, &mut lines)?;
        }
        if err != SearchResult::NotFound {
            found = true;
        }
        if !found && strip_version {
            let pos = acc.iter().position(|&c| c == b'.');
            let mut is_version = false;
            if let Some(pos) = pos {
                let ver_len = (acc.len() - pos) as i64 - 1;
                is_version = (1..=3).contains(&ver_len);
                for &c in &acc[pos + 1..] {
                    if !c.is_ascii_digit() {
                        is_version = false;
                        break;
                    }
                }
            }
            if let (true, Some(pos)) = (is_version, pos) {
                err = self.string_search_into(&acc[..pos], &mut lines)?;
                if !lines.is_empty() {
                    version_check = true;
                }
            }
        }
        if err != SearchResult::NotFound {
            found = true;
        }
        if !found {
            // Use CSeq_id to parse the id string and build a replacement, FASTA type string.
            // This allows some IDs, such as PDBs with chains, such as '1qcfA' to be parsed.
            let mut id = Vec::new();
            if let Ok(seqid) = seq_id_from_text(
                acc,
                parse_flags::RAW_TEXT | parse_flags::ANY_LOCAL,
                &mut |_: SeqIdWarning| {},
            ) {
                label::fasta(&seqid, &mut id);
            }
            if !id.is_empty() {
                err = self.string_search_into(&id, &mut lines)?;
            }
        }
        if err != SearchResult::NotFound {
            found = true;
        }
        if found {
            for line in &lines {
                let oid = crate::ncbi_str::c_atoi(&line.data);
                result
                    .oids
                    .push(u32::try_from(oid).map_err(|_| SeqDbError::Format {
                        path: self.data.path().to_path_buf(),
                        message: format!(
                            "string ISAM line '{}' holds the OID {oid}",
                            String::from_utf8_lossy(&line.key)
                        ),
                    })?);
            }
        }
        result.version_check = version_check;
        Ok(result)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2638-2658
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2724-2733
    /// The string case of `CSeqDBVol::x_StringToOids`: `StringToOids` of the simplified id with
    /// `version_check` asked for, then `x_CheckVersions` of the user's `acc` when the lookup
    /// stripped a version. `seq_ids` is `GetSeqIDs` of `oid` of the volume.
    pub fn string_id_to_oids(
        &self,
        acc: &[u8],
        simplified: &SimplifiedAccession,
        seq_ids: &mut dyn FnMut(u32) -> Result<Vec<ParsedSeqId>, SeqDbError>,
    ) -> Result<Vec<u32>, SeqDbError> {
        let found = self.string_to_oids(&simplified.str_id, simplified.simpler, true)?;
        let mut oids = found.oids;
        if found.version_check {
            check_versions(acc, &mut oids, seq_ids)?;
        }
        Ok(oids)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1371-1404
    /// `x_LoadIndex<string>`: the offsets of the pages in the data file (`samples + 1` of them)
    /// and, when `with_keys`, the sample keys (the text before the `\x02`).
    fn load_index(&self, with_keys: bool) -> Result<(Vec<u64>, Vec<Vec<u8>>), SeqDbError> {
        let samples = self.num_samples as usize;
        let bytes = if with_keys {
            self.index.read_all()?
        } else {
            self.index
                .read_range(0, KEY_SAMPLE_OFFSET + 4 * (samples as u64 + 1))?
        };
        let word = |at: u64| {
            let at = at as usize;
            u32::from_be_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
        };
        let mut offs = Vec::with_capacity(samples + 1);
        for index in 0..=samples as u64 {
            offs.push(u64::from(word(KEY_SAMPLE_OFFSET + 4 * index)));
        }
        let mut keys = Vec::new();
        if with_keys {
            let sample_begin = KEY_SAMPLE_OFFSET + 4 * (samples as u64 + 1);
            keys.reserve(samples);
            for index in 0..samples as u64 {
                let at = word(sample_begin + 4 * index) as usize;
                let tail = bytes.get(at..).unwrap_or(&[]);
                let end = tail
                    .iter()
                    .position(|&c| c == ISAM_DATA_CHAR)
                    .ok_or_else(|| SeqDbError::Format {
                        path: self.index.path().to_path_buf(),
                        message: format!("string ISAM sample key {index} has no end"),
                    })?;
                keys.push(tail[..end].to_vec());
            }
        }
        Ok((offs, keys))
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1406-1425
    /// `x_LoadData<string>`: `num_keys` records of `page`: the keys and the OIDs.
    fn load_data<'a>(
        &self,
        page: &'a [u8],
        num_keys: usize,
    ) -> Result<(Vec<&'a [u8]>, Vec<i32>), SeqDbError> {
        let bad = |what: &str| SeqDbError::Format {
            path: self.data.path().to_path_buf(),
            message: format!("string ISAM data page: {what}"),
        };
        // LOSAT: `num_keys` comes from the index header; a record takes at least 3 bytes of the
        // page (the key end, a digit, the line end), so no more are reserved than fit.
        let room = num_keys.min(page.len() / 3);
        let mut keys = Vec::with_capacity(room);
        let mut vals = Vec::with_capacity(room);
        let mut at = 0;
        for _ in 0..num_keys {
            let key_end = page[at..]
                .iter()
                .position(|&c| c == ISAM_DATA_CHAR)
                .ok_or_else(|| bad("a record has no key end"))?
                + at;
            keys.push(&page[at..key_end]);
            at = key_end + 1;
            let val_end = page[at..]
                .iter()
                .position(|&c| c == b'\n')
                .ok_or_else(|| bad("a record has no line end"))?
                + at;
            let oid = string_to_uint(&page[at..val_end])
                .ok()
                .and_then(|v| i32::try_from(v).ok())
                .ok_or_else(|| {
                    bad(&format!(
                        "'{}' is not an OID",
                        String::from_utf8_lossy(&page[at..val_end])
                    ))
                })?;
            vals.push(oid);
            at = val_end + 1;
        }
        Ok((keys, vals))
    }

    /// The number of records of page `sample_index` (`num_keys` of the batch functions).
    fn page_num_keys(&self, sample_index: usize) -> usize {
        if sample_index + 1 == self.num_samples as usize {
            (i64::from(self.num_terms) - sample_index as i64 * i64::from(self.page_size)) as usize
        } else {
            self.page_size as usize
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2540-2550
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1388-1419
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:548-623
    // The list is walked in the order it has after `InsureOrder` (`insure_order`): sorted where
    // its order flag was not set yet, else as it is. A list whose ids changed after the flag was
    // set (the alias files' `SEQIDLIST`, simplified again for each volume, `alias::ListSet`) can
    // be out of order: the walk then skips the ids below a key it has passed, as NCBI's does,
    // and where it meets an id below the key it starts from it is refused (`walk_before_start`).
    /// `IdsToOids` with `vol_start`, `vol_end` and a `CSeqDBGiList` for the string list
    /// (`x_TranslateGiList<string>`): sets `oid` of every entry of `ids` whose `si` (lower-cased,
    /// simplified: [`simplify_accession_for_isam`]) is a key of the volume, to the OID of the
    /// line plus `vol_start`. An entry that has an `oid` already is left as it is. `sorted` is
    /// the list's order flag (`m_CurrentOrder` at `eGi`), which [`insure_order`] sets.
    pub fn ids_to_oids(
        &self,
        vol_start: u32,
        ids: &mut [SiOid],
        sorted: &mut bool,
    ) -> Result<(), SeqDbError> {
        let gilist_size = ids.len();
        if gilist_size == 0 {
            return Ok(());
        }
        insure_order(ids, sorted);
        if !self.initialized {
            return Err(batch_mode_error());
        }
        let (page_offs, sample_keys) = self.load_index(true)?;
        let num_samples = self.num_samples as usize;
        let order_error = || walk_before_start(self.data.path());
        let mut gilist_index = 0;
        let mut sample_index: i64 = 0;
        while gilist_index < gilist_size && (sample_index as usize) < num_samples {
            advance_gi_list(ids, &mut gilist_index, &sample_keys[sample_index as usize]);
            if gilist_index >= gilist_size {
                break;
            }
            advance_key_list(
                &sample_keys,
                &mut sample_index,
                num_samples,
                &ids[gilist_index].si,
            );
            if sample_index < 0 {
                return Err(order_error());
            }
            // Now we should be ready to search a data block.
            let sample = sample_index as usize;
            let num_keys = self.page_num_keys(sample);
            let page = self.data.read_range(
                page_offs[sample],
                page_offs[sample + 1].max(page_offs[sample]),
            )?;
            let (keys, vals) = self.load_data(&page, num_keys)?;
            let mut index: i64 = 0;
            while gilist_index < gilist_size && (index as usize) < num_keys {
                advance_key_list(&keys, &mut index, num_keys, &ids[gilist_index].si);
                if index < 0 {
                    return Err(order_error());
                }
                set_translation(
                    ids,
                    &mut gilist_index,
                    keys[index as usize],
                    vals[index as usize] as u32,
                    vol_start,
                );
                index += 1;
                if index as usize >= num_keys {
                    break;
                }
                advance_gi_list(ids, &mut gilist_index, keys[index as usize]);
                set_translation(
                    ids,
                    &mut gilist_index,
                    keys[index as usize],
                    vals[index as usize] as u32,
                    vol_start,
                );
            }
            // We could be finished here because we exhausted the GI list.
            // We must be done with that one by now..
            sample_index += 1;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2582-2591
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1421-1459
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:333-418
    /// `x_SearchNegativeMultiSeq` (through `IdsToOids` with `vol_start`, `vol_end` and a
    /// `CSeqDBNegativeList`): reads every line of the volume and marks its OID (plus `vol_start`)
    /// in `out`: visible if its key is in `ids`, included if not. `ids` are the lower-cased
    /// simplified ids; they are sorted here (`CSeqDBNegativeList::InsureOrder`). `vol_end` is the
    /// global end OID of the volume (see the quirks above).
    pub fn negative_scan(
        &self,
        vol_start: u32,
        vol_end: u32,
        ids: &mut [Vec<u8>],
        out: &mut NegativeOids,
    ) -> Result<(), SeqDbError> {
        ids.sort();
        let gilist_size = ids.len();
        if gilist_size == 0 {
            return Ok(());
        }
        if !self.initialized {
            // Most ordinary errors (missing GIs for example) are ignored for "multi" mode
            // searches.  But if a GI list is specified, and cannot be interpreted, it is an
            // error.
            return Err(batch_mode_error());
        }
        let (page_offs, _) = self.load_index(false)?;
        let num_samples = self.num_samples as usize;
        let overflow = || SeqDbError::Format {
            path: self.data.path().to_path_buf(),
            message: "string ISAM OID plus the start of the volume overflows".to_string(),
        };
        let mut gilist_index = 0;
        let mut chunk = ScanChunk::default();
        for sample_index in 0..num_samples {
            // Now we should be ready to search a data block.
            let num_keys = self.page_num_keys(sample_index);
            let page = chunk.page(
                &self.data,
                page_offs[sample_index],
                page_offs[sample_index + 1].max(page_offs[sample_index]),
                page_offs[num_samples],
            )?;
            let (keys, vals) = self.load_data(page, num_keys)?;
            let mut i = 0;
            while i < num_keys {
                // 2. Look for it in the negative id list.
                let mut found = false;
                if gilist_index < gilist_size {
                    found = find_in_negative_list(ids, &mut gilist_index, keys[i]);
                }
                // 3. If not found, add the OID to the negative ID list.
                if i64::from(vals[i]) < i64::from(vol_end) {
                    let oid = (vals[i] as u32)
                        .checked_add(vol_start)
                        .ok_or_else(overflow)?;
                    if i32::try_from(oid).is_err() {
                        return Err(overflow());
                    }
                    if found {
                        // OID is found, but may not be included yet.
                        out.add_visible_oid(oid);
                        // If next accession is the same as current, but with version
                        if is_same_accession_next(&keys, num_keys, i) {
                            i += 1; // skip next check - sequence already excluded
                        }
                    } else if !is_same_accession_next(&keys, num_keys, i) {
                        // OID is included for iteration. Only include next accession if it is
                        // not the same as current (but with version) because it may be in
                        // exclude list. The check will be done in the next step
                        out.add_included_oid(oid);
                    }
                }
                i += 1;
            }
            // Move to next data page.  Note that for a negative ID list processing, we don't
            // actually fetch any samples, because every ID->OID line needs to be examined anyway.
        }
        Ok(())
    }
}

/// The OIDs `StringToOids` found and its `version_check` out-parameter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StringOids {
    pub oids: Vec<u32>,
    /// A version was stripped from the id to find these: [`check_versions`] must filter them.
    pub version_check: bool,
}

/// The data pages the negative scan has read (LOSAT: whole pages, up to [`SCAN_CHUNK`] bytes
/// at a time).
#[derive(Default)]
struct ScanChunk {
    start: u64,
    bytes: Vec<u8>,
}

impl ScanChunk {
    /// The bytes `from..to` of the data file, read ahead when they are not held (`limit` is the
    /// end of the last page).
    fn page(&mut self, data: &DbFile, from: u64, to: u64, limit: u64) -> Result<&[u8], SeqDbError> {
        let held = from >= self.start && to <= self.start + self.bytes.len() as u64;
        if !held {
            let end = limit.max(to).min(from.saturating_add(SCAN_CHUNK).max(to));
            self.bytes = data.read_range(from, end)?;
            self.start = from;
        }
        let at = (from - self.start) as usize;
        Ok(&self.bytes[at..at + (to - from) as usize])
    }
}

/// LOSAT's refusal where the walk of a list over the string ISAM of `path` steps before the
/// first sample key or the first key of a page (`s_AdvanceKeyList` at -1): NCBI reads the entry
/// before the start of its vector there.
fn walk_before_start(path: &Path) -> SeqDbError {
    SeqDbError::Unsupported {
        path: path.to_path_buf(),
        message: "NCBI's walk of an id list over the string ISAM meets an id below the key it \
                  starts from (a list out of order, or keys out of order) and reads before the \
                  start of its keys"
            .to_string(),
    }
}

/// NCBI's `CSeqDBException` of kind `eArgErr` with the message
/// `"Error: Unable to use ISAM index in batch mode."`.
pub(super) fn batch_mode_error() -> SeqDbError {
    SeqDbError::NotFound(
        "Error: Unable to use ISAM index in batch mode."
            .to_string()
            .into(),
    )
}

/// One entry of a list of string ids (`CSeqDBGiList::SSiOid`): the id and its OID, set by
/// [`StringIsam::ids_to_oids`] (NCBI: -1 while unknown).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiOid {
    pub si: Vec<u8>,
    pub oid: Option<u32>,
}

impl SiOid {
    /// `SSiOid` of `si` with the OID -1.
    pub fn new(si: impl Into<Vec<u8>>) -> Self {
        Self {
            si: si.into(),
            oid: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:574-586
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:589-606
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:616-651
// The flag is NCBI's `m_CurrentOrder`: once it is `eGi`, `InsureOrder` leaves the list as it is,
// whatever changed it since. NCBI's `std::sort` is not stable, which no result here shows: the
// entries it sorts have no OID yet, and equal ids are then equal entries.
/// `CSeqDBGiList::InsureOrder` with `eGi` for the string list: where the flag `sorted` is not set,
/// sorts by `si` (`std::string` order: unsigned bytes, then length) unless it is sorted already,
/// and sets the flag.
pub fn insure_order(ids: &mut [SiOid], sorted: &mut bool) {
    if *sorted {
        return;
    }
    let already = ids.windows(2).all(|w| w[1].si >= w[0].si);
    if !already {
        ids.sort_by(|a, b| a.si.cmp(&b.si));
    }
    *sorted = true;
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:50-75
/// `s_AdvanceGiList` for strings: moves `gi_index` to the first entry that is `>= key` (by
/// jumps of doubling length), then past the entries already translated.
fn advance_gi_list(gis: &[SiOid], gi_index: &mut usize, key: &[u8]) {
    let gis_size = gis.len();
    while *gi_index < gis_size && gis[*gi_index].si.as_slice() < key {
        *gi_index += 1;
        let mut jump = 2;
        while *gi_index + jump < gis_size && gis[*gi_index + jump].si.as_slice() < key {
            *gi_index += jump;
            jump *= 2;
        }
    }
    // skipping translated elements
    while *gi_index < gis_size && gis[*gi_index].oid.is_some() {
        *gi_index += 1;
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:77-99
/// `s_AdvanceKeyList`: moves `index` to the last key that is `<= target`; -1 when none is.
fn advance_key_list<K: AsRef<[u8]>>(keys: &[K], index: &mut i64, num_keys: usize, target: &[u8]) {
    let num_keys = num_keys as i64;
    while *index < num_keys && keys[*index as usize].as_ref() <= target {
        *index += 1;
        let mut jump = 2;
        while *index + jump < num_keys && keys[(*index + jump) as usize].as_ref() <= target {
            *index += jump;
            jump *= 2;
        }
    }
    *index -= 1;
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:101-115
/// `s_SetTranslation` for strings: the OID of every entry from `gi_index` on that equals `key`.
fn set_translation(
    gis: &mut [SiOid],
    gi_index: &mut usize,
    key: &[u8],
    value: u32,
    vol_start: u32,
) {
    while *gi_index < gis.len() && gis[*gi_index].si.as_slice() == key {
        gis[*gi_index].oid = Some(value.wrapping_add(vol_start));
        *gi_index += 1;
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1460-1489
/// `x_FindInNegativeList` for strings: moves `index` past the ids `< key`; whether the id it
/// stops at is `key`.
fn find_in_negative_list(ids: &[Vec<u8>], index: &mut usize, key: &[u8]) -> bool {
    let ids_size = ids.len();
    // Skip any that are less than key.
    while *index < ids_size && ids[*index].as_slice() < key {
        *index += 1;
        let mut jump = 2;
        while *index + jump < ids_size && ids[*index + jump].as_slice() < key {
            *index += jump;
            jump += jump;
        }
    }
    // Check whether the GI or TI was found.
    *index < ids_size && ids[*index].as_slice() == key
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:304-319
/// `s_IsSameAccession` with `acc1` and `acc2`: `acc2` has no `.`, `acc1` has one and holds `acc2`,
/// and `acc2` is what `acc1` has before its first `.`.
fn is_same_accession(acc1: &[u8], acc2: &[u8]) -> bool {
    if !acc2.contains(&b'.') {
        // no version in acc2
        let holds = acc2.is_empty() || acc1.windows(acc2.len()).any(|w| w == acc2);
        if acc1.contains(&b'.') && holds {
            let dot = acc1.iter().position(|&c| c == b'.').unwrap_or(acc1.len());
            if acc2 == &acc1[..dot] {
                return true;
            }
        }
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:321-330
/// `s_IsSameAccession` with `keys`, `num_keys` and `currIndex`: whether line `curr_index + 1` of
/// the page is line `curr_index` with a version.
fn is_same_accession_next(keys: &[&[u8]], num_keys: usize, curr_index: usize) -> bool {
    curr_index + 1 < num_keys && is_same_accession(keys[curr_index + 1], keys[curr_index])
}

/// `CSeqDBBitVector`: `Size` is the highest bit set plus one; `GetBit` past it is false.
#[derive(Clone, Debug, Default)]
struct BitVector {
    bitmap: Vec<u32>,
    size: u32,
}

impl BitVector {
    fn set_bit(&mut self, oid: u32) {
        if oid >= self.size {
            self.resize(oid + 1);
        }
        self.bitmap[(oid / 32) as usize] |= 1 << (oid % 32);
    }

    fn get_bit(&self, oid: u32) -> bool {
        if oid >= self.size {
            return false;
        }
        self.bitmap[(oid / 32) as usize] & (1 << (oid % 32)) != 0
    }

    /// `x_Resize`.
    fn resize(&mut self, num: u32) {
        let need = u64::from(num).div_ceil(32);
        if (self.bitmap.len() as u64) < need {
            let mut new_size: u64 = 1024;
            while new_size < need {
                new_size *= 2;
            }
            self.bitmap.resize(new_size as usize, 0);
        }
        self.size = num;
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:976-1014
/// The OID sets of `CSeqDBNegativeList` that [`StringIsam::negative_scan`] fills (see the module
/// documentation for what they mean).
#[derive(Clone, Debug, Default)]
pub struct NegativeOids {
    included: BitVector,
    visible: BitVector,
}

impl NegativeOids {
    /// An empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// `AddIncludedOid`.
    pub fn add_included_oid(&mut self, oid: u32) {
        self.included.set_bit(oid);
    }

    /// `AddVisibleOid`.
    pub fn add_visible_oid(&mut self, oid: u32) {
        self.visible.set_bit(oid);
    }

    /// `GetOidStatus`: true for the OIDs that were included and for those never seen.
    pub fn get_oid_status(&self, oid: u32) -> bool {
        self.included.get_bit(oid) || !self.visible.get_bit(oid)
    }

    /// `GetNumOids`: the larger size of the two bit vectors.
    pub fn num_oids(&self) -> u32 {
        self.visible.size.max(self.included.size)
    }

    /// `GetBit` of `oid` in `m_Included`.
    #[cfg(test)]
    pub fn is_included(&self, oid: u32) -> bool {
        self.included.get_bit(oid)
    }

    /// `GetBit` of `oid` in `m_Visible`.
    #[cfg(test)]
    pub fn is_visible(&self, oid: u32) -> bool {
        self.visible.get_bit(oid)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2735-2789
/// `CSeqDBVol::x_CheckVersions`: removes the OIDs whose Seq-ids (`seq_ids` is `GetSeqIDs` of `oid`
/// of the volume) hold no Textseq-id with the accession and version of the user's `acc`.
pub fn check_versions(
    acc: &[u8],
    oids: &mut Vec<u32>,
    seq_ids: &mut dyn FnMut(u32) -> Result<Vec<ParsedSeqId>, SeqDbError>,
) -> Result<(), SeqDbError> {
    // If we resolved a string id by stripping the version off of the string, we need to check
    // (for each OID) if the real ID had a matching version.
    let pos = acc.iter().position(|&c| c == b'.');
    // The `find` of a "." in `acc` is asserted to succeed; without it `pos + 1` wraps to 0 and
    // `nover` is `acc`.
    let ver_str = &acc[pos.map_or(0, |p| p + 1)..];
    // `fConvErr_NoThrow`: 0 for a text NCBI cannot convert.
    let vernum = string_to_int(ver_str, ALLOW_TRAILING_SYMBOLS).unwrap_or(0);
    let mut nover = &acc[..pos.unwrap_or(acc.len())];
    while let Some(bar) = nover.iter().position(|&c| c == b'|') {
        nover = &nover[bar + 1..];
    }
    let mut kept = Vec::with_capacity(oids.len());
    for &oid in oids.iter() {
        let ids = seq_ids(oid)?;
        let found = ids.iter().any(|seqid| {
            seqid.textseq().is_some_and(|id| {
                id.accession.as_deref() == Some(nover) && id.version == Some(vernum)
            })
        });
        if found {
            kept.push(oid);
        }
    }
    *oids = kept;
    Ok(())
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:1951-1958
/// `ESeqDBIdType`, without `eHashId` (a Seq-id is never simplified to a hash).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqDbIdType {
    Gi,
    Ti,
    Pig,
    String,
    Oid,
}

/// The results of `SeqDB_SimplifyAccession` with `acc`, `num_id`, `str_id` and `simpler`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SimplifiedAccession {
    /// The returned `ESeqDBIdType`.
    pub id_type: SeqDbIdType,
    /// `num_id`: the GI, TI, PIG or OID (-1 cast to `Uint4` when none was set).
    pub num_id: i64,
    /// `str_id`: the key for the string index (not yet lower-cased: see
    /// [`simplify_accession_for_isam`]).
    pub str_id: Vec<u8>,
    /// `simpler`: the id was devolved to a number or to the text of a Seq-id.
    pub simpler: bool,
}

/// An exception of `SeqDB_SimplifyAccession`: NCBI lets it out of the call.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SimplifyError {
    /// A `CSeqIdException` that `CSeq_id` throws for `element` or for the `db:tag` id.
    SeqId(SeqIdError),
    /// Another exception of NCBI's: a string tag read as a number, an integer read from a string
    /// tag.
    Exception(NcbiException),
}

impl std::fmt::Display for SimplifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimplifyError::SeqId(error) => f.write_str(&String::from_utf8_lossy(&error.what())),
            SimplifyError::Exception(exception) => {
                f.write_str(&String::from_utf8_lossy(exception.get_msg()))
            }
        }
    }
}

impl std::error::Error for SimplifyError {}

// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:122-150
/// `CObject_id::GetId()`: a string tag is NCBI's `CInvalidChoiceSelection` (probe D11 of line
/// BDC), thrown by the code that NCBI's `datatool` generates for the CHOICE.
fn object_id_id(tag: &ObjectId) -> Result<i64, SimplifyError> {
    match tag {
        ObjectId::Id(n) => Ok(i64::from(*n)),
        ObjectId::Str(_) => Err(SimplifyError::Exception(NcbiException::new(
            ThrowSite::generated("CInvalidChoiceSelection::eFail"),
            "CObject_id::GetId(): Invalid choice selection: NCBI-General::Object-id.str",
        ))),
    }
}

fn fasta_string(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    label::fasta(id, &mut out);
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2264-2427
/// `SeqDB_SimplifySeqid` with `bestid`, the address of `acc`, `num_id`, `str_id` and `simpler` (the
/// call with no `acc` is [`seqdb_simplify_seqid`]).
fn simplify_seqid(bestid: &ParsedSeqId, acc: &[u8]) -> Result<SimplifiedAccession, SimplifyError> {
    let mut out = SimplifiedAccession {
        id_type: SeqDbIdType::String,
        num_id: i64::from(u32::MAX),
        str_id: Vec::new(),
        simpler: false,
    };
    let mut tsip: Option<&TextseqId> = None;
    let mut matched = true;
    match bestid {
        ParsedSeqId::Gi(gi) => {
            out.simpler = true;
            out.num_id = *gi;
            out.id_type = SeqDbIdType::Gi;
        }
        ParsedSeqId::Gibbsq(gibbsq) => {
            out.simpler = true;
            out.id_type = SeqDbIdType::String;
            out.str_id = (*gibbsq as u32).to_string().into_bytes();
        }
        ParsedSeqId::General(dbt) => 'general: {
            if dbt.db == b"BL_ORD_ID" {
                out.simpler = true;
                out.num_id = object_id_id(&dbt.tag)?;
                out.id_type = SeqDbIdType::Oid;
                break 'general;
            }
            if dbt.db == b"PIG" {
                out.simpler = true;
                out.num_id = object_id_id(&dbt.tag)?;
                out.id_type = SeqDbIdType::Pig;
                break 'general;
            }
            if dbt.db == b"ti" {
                out.simpler = true;
                out.num_id = match &dbt.tag {
                    ObjectId::Str(text) => {
                        string_to_int8(text, 0).map_err(SimplifyError::Exception)?
                    }
                    ObjectId::Id(n) => i64::from(*n),
                };
                out.id_type = SeqDbIdType::Ti;
                break 'general;
            }
            if dbt.db.eq_ignore_ascii_case(b"GNOMON") {
                out.str_id = fasta_string(bestid).to_ascii_lowercase();
                out.id_type = SeqDbIdType::String;
                break 'general;
            }
            if let ObjectId::Str(tag) = &dbt.tag {
                out.id_type = SeqDbIdType::String;
                out.str_id = tag.to_ascii_lowercase();
            } else {
                // Use the default logic.
                matched = false;
            }
        }
        ParsedSeqId::Local(objid) => {
            out.simpler = true;
            out.id_type = SeqDbIdType::String;
            match objid {
                // sparse version will leave "lcl|" off.
                ObjectId::Str(text) => out.str_id = text.to_ascii_lowercase(),
                // Local numeric ids are stored as strings.
                ObjectId::Id(n) => out.str_id = format!("lcl|{n}").into_bytes(),
            }
        }
        // tsip types
        ParsedSeqId::Embl(id)
        | ParsedSeqId::Ddbj(id)
        | ParsedSeqId::Genbank(id)
        | ParsedSeqId::Tpg(id)
        | ParsedSeqId::Tpe(id)
        | ParsedSeqId::Tpd(id)
        | ParsedSeqId::Other(id)
        | ParsedSeqId::Swissprot(id)
        | ParsedSeqId::Gpipe(id)
        | ParsedSeqId::Pir(id)
        | ParsedSeqId::Prf(id) => tsip = Some(id),
        _ => matched = false,
    }
    // Default: if we have a string, use it; if we only have seqid, create a string.  This
    // should not happen if the seqid matches one of the cases above, which currently
    // correspond to all the supported seqid types.
    if !matched {
        // (should not happen normally)
        out.simpler = false;
        out.id_type = SeqDbIdType::String;
        out.str_id = acc.to_ascii_lowercase();
    }
    if let Some(tsip) = tsip {
        let mut found = false;
        if let Some(accession) = &tsip.accession {
            out.str_id = accession.to_ascii_lowercase();
            found = true;
            if let Some(version) = tsip.version {
                out.str_id.push(b'.');
                out.str_id
                    .extend_from_slice((version as u32).to_string().as_bytes());
            }
        } else if let Some(name) = &tsip.name {
            out.str_id = name.to_ascii_lowercase();
            found = true;
        }
        if found {
            out.simpler = true;
            out.id_type = SeqDbIdType::String;
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2805-2817
/// `SeqDB_SimplifySeqid` with `seqid`, no accession, `num_id`, `str_id` and `simpler` of
/// `CSeqDBVol::SeqidToOids`: what a volume looks up for the Seq-id `seqid`. Without `acc`, NCBI's
/// default case labels the Seq-id itself (`GetLabel` with `eFasta`), which is its FASTA string, the
/// `acc` given here.
pub fn seqdb_simplify_seqid(seqid: &ParsedSeqId) -> Result<SimplifiedAccession, SimplifyError> {
    simplify_seqid(seqid, &fasta_string(seqid))
}

/// `WhichInverseSeqId` of `portion.c_str()`: the text up to its first NUL.
fn which_inverse(text: &[u8]) -> u32 {
    let end = text.iter().position(|&c| c == 0).unwrap_or(text.len());
    which_inverse_seq_id(&text[..end])
}

fn find_from(s: &[u8], byte: u8, from: usize) -> Option<usize> {
    s.get(from..)?
        .iter()
        .position(|&c| c == byte)
        .map(|p| p + from)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2441-2482
/// `s_SeqDB_EndOfFastaID`: the end of the first Seq-id of `s` from `pos`, `None` (NCBI's
/// `npos`) for a text that does not start with a Seq-id tag and `|`.
fn end_of_fasta_id(s: &[u8], pos: usize) -> Option<usize> {
    let vbar = find_from(s, b'|', pos)?;
    if which_inverse(&s[pos..vbar]) == SEQID_NOT_SET {
        return None; // bad
    }
    let mut vbar_prev = vbar;
    let end = loop {
        let Some(vbar) = find_from(s, b'|', vbar_prev + 1) else {
            break s.len();
        };
        let element = &s[vbar_prev + 1..vbar];
        if which_inverse(element) != SEQID_NOT_SET {
            break vbar_prev;
        }
        vbar_prev = vbar;
    };
    Some(end)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2496-2531
/// `s_SeqDB_ParseSeqIDs`: the Seq-ids at the start of `line` (`CSeq_id` of `element` of every
/// element up to the first text that is not one). A `CSeqIdException` is not caught by NCBI's
/// handler of `invalid_argument`.
fn parse_seq_ids(line: &[u8]) -> Result<Vec<ParsedSeqId>, SimplifyError> {
    let mut seqids = Vec::new();
    let mut pos = 0;
    while pos < line.len() {
        let Some(end) = end_of_fasta_id(line, pos) else {
            // We didn't get a clean parse -- ignore the data after this point, and return what
            // we have.
            break;
        };
        let element = &line[pos..end];
        let id = seq_id_from_text(element, parse_flags::ANY_RAW, &mut |_: SeqIdWarning| {})
            .map_err(SimplifyError::SeqId)?;
        seqids.push(id);
        pos = end + 1;
    }
    Ok(seqids)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2535-2608
/// `SeqDB_SimplifyAccession` with `acc`, `num_id`, `str_id` and `simpler`: what the volume looks up
/// for the user's id `acc`.
pub fn seqdb_simplify_accession(acc: &[u8]) -> Result<SimplifiedAccession, SimplifyError> {
    let seqid_set = parse_seq_ids(acc)?;
    if let Some(bestid) = find_best(seqid_set.iter(), best_rank_score) {
        // Something like SeqIdFindBest()
        return simplify_seqid(bestid, acc);
    }
    let mut out = SimplifiedAccession {
        id_type: SeqDbIdType::String,
        num_id: i64::from(u32::MAX),
        str_id: Vec::new(),
        simpler: false,
    };
    // Check for bare pdb accession with underscore (such as 12AS_A). These are not in the isam
    // index and need to be translated to the standard form (pdb|12AS|A).
    let seqids: Vec<ParsedSeqId> = match parse_ids(
        acc,
        parse_flags::RAW_TEXT | parse_flags::ANY_LOCAL,
        &mut |_: SeqIdWarning| {},
    ) {
        Ok(ids) => ids.into_iter().flatten().collect(),
        Err(_) => Vec::new(),
    };
    let front = seqids.first();
    if let (Some(pdb @ ParsedSeqId::Pdb(_)), true) = (front, acc.contains(&b'_')) {
        out.str_id = fasta_string(pdb).to_ascii_lowercase();
    } else if let Some(ParsedSeqId::Local(_)) = front {
        // Chec for gnl dbs
        if let Some(colon) = acc.iter().position(|&c| c == b':') {
            const GNL_DBS: [&[u8]; 4] = [b"CDD", b"SRA", b"TSA", b"GNOMON"];
            let (db_tag, gnl_id) = (&acc[..colon], &acc[colon + 1..]);
            if GNL_DBS.iter().any(|db| db.eq_ignore_ascii_case(db_tag)) {
                let mut str_id = b"gnl|".to_vec();
                str_id.extend_from_slice(db_tag);
                str_id.push(b'|');
                str_id.extend_from_slice(gnl_id);
                seq_id_from_text(&str_id, parse_flags::ANY_RAW, &mut |_: SeqIdWarning| {})
                    .map_err(SimplifyError::SeqId)?;
                out.str_id = str_id;
            } else {
                out.str_id = acc.to_vec();
            }
        } else if string_to_int8(acc, 0).is_ok() {
            out.str_id = [b"lcl|".as_slice(), acc].concat();
        } else {
            out.str_id = acc.to_vec();
        }
    } else {
        out.str_id = acc.to_vec();
    }
    out.id_type = SeqDbIdType::String;
    out.simpler = false;
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2610-2618
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:608-614
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1612-1618
/// `SeqDB_SimplifyAccession` of `acc` followed by `NStr::ToLower`
/// (`PreprocessIdsForISAMSiLookup`): the key of the string index for a user's id; empty for an id
/// that devolves to a number (a GI, TI, PIG or OID).
pub fn simplify_accession_for_isam(id: &[u8]) -> Result<Vec<u8>, SimplifyError> {
    let simplified = seqdb_simplify_accession(id)?;
    Ok(if simplified.id_type == SeqDbIdType::String {
        simplified.str_id.to_ascii_lowercase()
    } else {
        Vec::new()
    })
}

#[cfg(test)]
#[path = "isam_tests.rs"]
mod tests;
