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
//!   [`StringIsam::from_files`]; `x_OutOfBounds(string)`: the private `out_of_bounds`;
//!   `GetIdBounds(string)`: [`StringIsam::id_bounds`].
//! - `x_StringSearch` with `x_DiffSample`, `x_DiffCharLease`, `x_DiffChar`, `x_LoadPage`,
//!   `x_ExtractAllData`, `x_ExtractPageData`, `x_ExtractData`: [`StringIsam::string_search`]
//!   (every line of the volume that matches a term, as `(key, data)` strings).
//! - `StringToOids`: [`StringIsam::string_to_oids`]; `CSeqDBVol::x_CheckVersions`:
//!   [`check_versions`]; the string case of `CSeqDBVol::x_StringToOids`:
//!   [`StringIsam::string_id_to_oids`].
//! - `IdsToOids(CSeqDBGiList&)` (`x_TranslateGiList<string>`): [`StringIsam::ids_to_oids`];
//!   `CSeqDBGiList::InsureOrder` for the string list: [`insure_order`].
//! - `IdsToOids(CSeqDBNegativeList&)` (`x_SearchNegativeMultiSeq`): [`StringIsam::negative_scan`];
//!   `CSeqDBNegativeList` and `CSeqDBBitVector`: [`NegativeOids`].
//! - `SeqDB_SimplifyAccession`, `SeqDB_SimplifySeqid`, `PreprocessIdsForISAMSiLookup`:
//!   [`seqdb_simplify_accession`], [`simplify_accession_for_isam`].
//!
//! Not ported, because nothing reaches it: the numeric ISAM (`.pni`, `.nni`, ...), `x_GetIndexKeyOffset`
//! and `x_GetIndexString` (used only behind `short_match` and `follow_match`, two constants that
//! are `false`), `x_SparseStringToOids` (`_TROUBLE`), `HashToOids`, and the numeric
//! `x_SearchNegativeMulti`.
//!
//! # `CSeqDBNegativeList`: what the OID bit vectors mean
//!
//! A negative id list is applied to a volume by [`StringIsam::negative_scan`]. Every line of the
//! ISAM is read. A line whose key is in the sorted id list marks its OID *visible*
//! (`AddVisibleOid`); any other line marks it *included* (`AddIncludedOid`). `GetOidStatus(oid)` is
//! `included(oid) || !visible(oid)`: an OID is kept unless it was seen only through listed keys,
//! so a sequence with one listed key and one unlisted key (an accession and its `.version`
//! are two lines of one OID) stays. `GetNumOids()` is the larger of the two vectors' sizes (the
//! highest OID set plus one); `CSeqDBOIDList::x_ApplyNegativeList` (`seqdboidlist.cpp:376-404`)
//! clears every OID of the volume from `GetNumOids()` on and every OID below it whose status is
//! false. Both vectors are `CSeqDBBitVector`s: `SetBit` grows `Size` to `oid + 1`, `GetBit` is
//! false past `Size`.
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
//!   `CSeq_id(acc, fParse_RawText | fParse_AnyLocal).AsFastaString()`; the `OID`s come from
//!   `atoi` of the data strings. [`check_versions`] then keeps only OIDs whose Seq-ids hold the
//!   accession (the text after the last `|` of `acc` before its first `.`) with the version the
//!   text after that `.` starts with.
//! - The positive batch (`x_TranslateGiList<string>`) is not "first ISAM line wins". For a key
//!   that has several lines in one page, which line sets the OID of the equal list entries
//!   depends on the path: an entry reached by `s_AdvanceKeyList`, which stops at the *last* line
//!   with a key `<=` the id, takes that last line; an entry reached right after the previous
//!   line was used (`++index`, `s_AdvanceGiList`, `s_SetTranslation` with `keys[index]`) takes
//!   the line at `index`, the *first* line of the key when the previous list entry was a
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
//!   shorter than its offset tables. `s_AdvanceKeyList` can step to `-1` on an unsorted data
//!   file (NCBI indexes `keys[-1]`): `Err` here.
//! - A volume whose files fail `x_InitSearch` (wrong version, type or data length) is opened as
//!   NCBI does, uninitialised: lookups find nothing, and the batch calls fail with NCBI's
//!   "Error: Unable to use ISAM index in batch mode.". NCBI's `x_FindIndexBounds` still reads the
//!   index of such a volume; its result is never used and is not computed.
//! - The warnings that `CSeq_id` posts while parsing a user's id (`ERR_POST` to the diagnostics
//!   stream) are dropped.
//! - Exceptions of `SeqDB_SimplifyAccession` that are not NCBI's text (an out-of-range number in
//!   a `ti` tag, a string tag of `BL_ORD_ID`) are [`SimplifyError::Exception`] with LOSAT's text.
//! - The functions assume the C locale (`tolower`, `toupper`, `isdigit`, `atoi`).

use std::path::{Path, PathBuf};

use super::file::DbFile;
use super::SeqDbError;
use crate::seqid::accguide::{parse_flags, SeqIdWarning, SEQID_NOT_SET};
use crate::seqid::label;
use crate::seqid::model::{ObjectId, ParsedSeqId, TextseqId};
use crate::seqid::parse::{parse_ids, seq_id_from_text, which_inverse_seq_id, SeqIdError};
use crate::seqid::rank::{best_rank_score, find_best};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:45-55
// ```c++
// /// Format version of the ISAM files
// #define ISAM_VERSION 1
//
// /// Default page size for numeric indices
// #define DEFAULT_NISAM_SIZE 256
//
// /// Default page size for string indices
// #define DEFAULT_SISAM_SIZE 64
//
// /// Special page size value which indicates a memory-only string index
// #define MEMORY_ONLY_PAGE_SIZE 1
// ```
/// `ISAM_VERSION`.
const ISAM_VERSION: i32 = 1;
/// `CSeqDBIsam::eString`.
const ISAM_STRING_TYPE: i32 = 2;
/// `MEMORY_ONLY_PAGE_SIZE`.
const MEMORY_ONLY_PAGE_SIZE: i32 = 1;
/// `m_KeySampleOffset` (`x_InitSearch`): the nine header fields.
const KEY_SAMPLE_OFFSET: u64 = 9 * 4;
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:580-581
// ```c++
// /// The terminating character for string ISAM keys when data is present.
// const char ISAM_DATA_CHAR = (char) 2;
// ```
/// `ISAM_DATA_CHAR`: the byte that ends a key when data follows.
const ISAM_DATA_CHAR: u8 = 2;
/// The bytes of the data file the negative scan reads at once (whole pages; LOSAT: NCBI maps
/// the file).
const SCAN_CHUNK: u64 = 4 << 20;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:560-578
// ```c++
// /// Return NUL for nulls or EOL characters
// ///
// /// This function returns a NUL byte for any of NUL, CR, or NL.  This
// /// is done because these characters are used to terminate the
// /// variable length records in a string-based ISAM file.
// ///
// /// @param c
// ///   A character
// /// @return
// ///   NUL or the same character
// static inline char
// s_SeqDBIsam_NullifyEOLs(char c)
// {
//     if (SEQDB_ISEOL(c)) {
//         return 0;
//     } else {
//         return c;
//     }
// }
// ```
/// `s_SeqDBIsam_NullifyEOLs`.
fn nullify_eols(c: u8) -> u8 {
    if c == b'\n' || c == b'\r' {
        0
    } else {
        c
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:583-587
// ```c++
// /// Returns true if the character is a terminator for an ISAM key.
// static inline bool ENDS_ISAM_KEY(char P)
// {
//     return (P == ISAM_DATA_CHAR) || (s_SeqDBIsam_NullifyEOLs(P) == 0);
// }
// ```
/// `ENDS_ISAM_KEY`.
fn ends_isam_key(c: u8) -> bool {
    c == ISAM_DATA_CHAR || nullify_eols(c) == 0
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:589-632
// ```c++
// Int4 CSeqDBIsam::x_DiffChar(const string & term_in,
//                             const char   * begin,
//                             const char   * end,
//                             bool           ignore_case)
// {
//     int result(-1);
//     int i(0);
//
//     const char * file_data = begin;
//     int bytes = int(end - begin);
//
//     for(i = 0; (i < bytes) && i < (int) term_in.size(); i++) {
//         char ch1 = term_in[i];
//         char ch2 = file_data[i];
//
//         if (ch1 != ch2) {
//             ch1 = s_SeqDBIsam_NullifyEOLs(ch1);
//             ch2 = s_SeqDBIsam_NullifyEOLs(ch2);
//
//             if (ignore_case) {
//                 ch1 = toupper((unsigned char) ch1);
//                 ch2 = toupper((unsigned char) ch2);
//             }
//
//             if (ch1 != ch2) {
//                 break;
//             }
//         }
//     }
//
//     const char * p = file_data + i;
//
//     while((p < end) && ((*p) == ' ')) {
//         p++;
//     }
//
//     if (((p == end) || ENDS_ISAM_KEY(*p)) && (i == (int) term_in.size())) {
//         result = -1;
//     } else {
//         result = i;
//     }
//
//     return result;
// }
// ```
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
// ```c++
// void CSeqDBIsam::x_ExtractData(const char     * key_start,
//                                const char     * map_end,
//                                vector<string> & keys_out,
//                                vector<string> & data_out)
// {
//     const char * data_ptr(0);
//     const char * p(key_start);
//
//     while(p < map_end) {
//         switch(s_SeqDBIsam_NullifyEOLs(*p)) {
//         case 0:
//             if (data_ptr) {
//                 keys_out.push_back(string(key_start, data_ptr));
//                 data_out.push_back(string(data_ptr+1, p));
//             } else {
//                 keys_out.push_back(string(key_start, p));
//                 data_out.push_back("");
//             }
//             return;
//
//         case ISAM_DATA_CHAR:
//             data_ptr = p;
//
//         default:
//             p++;
//         }
//     }
// }
// ```
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
// ```c++
// void CSeqDBIsam::x_ExtractPageData(const string   & term_in,
//                                    TIndx            page_index,
//                                    const char     * beginp,
//                                    const char     * endp,
//                                    vector<TIndx>  & indices_out,
//                                    vector<string> & keys_out,
//                                    vector<string> & data_out)
// {
//     // Collect all 'good' data from the page.
//
//     bool ignore_case = true;
//
//     Uint4 TermNum(0);
//
//     const char * indexp(beginp);
//     bool found_match(false);
//
//     while (indexp < endp) {
//         Int4 Diff = x_DiffChar(term_in,
//                                indexp,
//                                endp,
//                                ignore_case);
//
//         if (Diff == -1) { // Complete match
//             found_match = true;
//
//             x_ExtractData(indexp,
//                           endp,
//                           keys_out,
//                           data_out);
//
//             indices_out.push_back(page_index + TermNum);
//         } else {
//             // If we found a match, but the current term doesn't
//             // match, then we are past the set of matching entries.
//
//             if (found_match) {
//                 break;
//             }
//         }
//
//         // Skip remainder of term, and any nulls after it.
//
//         while((indexp < endp) && s_SeqDBIsam_NullifyEOLs(*indexp)) {
//             indexp++;
//         }
//         while((indexp < endp) && (! s_SeqDBIsam_NullifyEOLs(*indexp))) {
//             indexp++;
//         }
//
//         TermNum++;
//     }
// }
// ```
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
}

impl StringIsam {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1173-1210
    // ```c++
    // void CSeqDBIsam::x_MakeFilenames(const string & dbname,
    //                                  char           prot_nucl,
    //                                  char           file_ext_char,
    //                                  string       & index_name,
    //                                  string       & data_name)
    // {
    //     if (dbname.empty() ||
    //         (! isalpha((unsigned char) prot_nucl)) ||
    //         (! isalpha((unsigned char) file_ext_char))) {
    //
    //         NCBI_THROW(CSeqDBException,
    //                    eArgErr,
    //                    "Error: argument not valid");
    //     }
    //
    //     index_name.reserve(dbname.size() + 4);
    //     data_name.reserve(dbname.size() + 4);
    //
    //     index_name = dbname;
    //     index_name += '.';
    //     index_name += prot_nucl;
    //     index_name += file_ext_char;
    //
    //     data_name = index_name;
    //     index_name += 'i';
    //     data_name  += 'd';
    // }
    //
    // bool CSeqDBIsam::IndexExists(const string & dbname,
    //                              char           prot_nucl,
    //                              char           file_ext_char)
    // {
    //     string iname, dname;
    //     x_MakeFilenames(dbname, prot_nucl, file_ext_char, iname, dname);
    //
    //     return CFile(iname).Exists() && CFile(dname).Exists();
    // }
    //
    // ```
    /// `CSeqDBIsam::x_MakeFilenames` (`file_ext_char` is `s`): `<dbname>.ps?` or `.ns?`, the index
    /// file ending in `i` and the data file in `d`.
    pub fn file_names(dbname: &str, protein: bool) -> Result<(PathBuf, PathBuf), SeqDbError> {
        if dbname.is_empty() {
            return Err(SeqDbError::NotFound(
                "Error: argument not valid".to_string(),
            ));
        }
        let mut index = format!("{dbname}.{}s", if protein { 'p' } else { 'n' });
        let mut data = index.clone();
        index.push('i');
        data.push('d');
        Ok((PathBuf::from(index), PathBuf::from(data)))
    }

    /// `CSeqDBIsam::IndexExists` for the string index.
    pub fn index_exists(dbname: &str, protein: bool) -> bool {
        match Self::file_names(dbname, protein) {
            Ok((index, data)) => index.exists() && data.exists(),
            Err(_) => false,
        }
    }

    /// The string index of the volume `dbname` (its path without extension).
    pub fn open_volume(dbname: &str, protein: bool) -> Result<Self, SeqDbError> {
        let (index, data) = Self::file_names(dbname, protein)?;
        Self::open(&index, &data)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1102-1171
    // ```c++
    // CSeqDBIsam::CSeqDBIsam(CSeqDBAtlas  & atlas,
    //                        const string & dbname,
    //                        char           prot_nucl,
    //                        char           file_ext_char,
    //                        ESeqDBIdType   ident_type)
    //     : m_Atlas          (atlas),
    //       m_IdentType      (ident_type),
    //       m_IndexLease     (atlas),
    //       m_DataLease      (atlas),
    //       m_Type           (eNumeric),
    //       m_NumTerms       (0),
    //       m_NumSamples     (0),
    //       m_PageSize       (0),
    //       m_MaxLineSize    (0),
    //       m_IdxOption      (0),
    //       m_Initialized    (false),
    //       m_KeySampleOffset(0),
    //       m_TestNonUnique  (true),
    //       m_FileStart      (0),
    //       m_FirstOffset    (0),
    //       m_LastOffset     (0),
    //       m_LongId         (false),
    //       m_TermSize       (8)
    // {
    //     // These are the types that readdb.c seems to use.
    //
    //     switch(ident_type) {
    //     case eGiId:
    //     case ePigId:
    //     case eTiId:
    //         m_Type = eNumeric;
    //         break;
    //
    //     case eStringId:
    //     case eHashId:
    //         m_Type = eString;
    //         break;
    //
    //     default:
    //         NCBI_THROW(CSeqDBException,
    //                    eArgErr,
    //                    "Error: ident type argument not valid");
    //     }
    //
    //     x_MakeFilenames(dbname,
    //                     prot_nucl,
    //                     file_ext_char,
    //                     m_IndexFname,
    //                     m_DataFname);
    //     try {
    //         m_IndexLease.Init(m_IndexFname);
    //         m_DataLease.Init(m_DataFname);
    //     } catch (CException & e){
    //         m_IndexLease.Clear();
    //         m_DataLease.Clear();
    //         string msg("Error: Could not open input file (");
    //         msg += m_IndexFname + "/" + m_DataFname + ")";
    //         NCBI_THROW(CSeqDBException, eFileErr, msg);
    //
    //     }
    //     if(m_Type == eNumeric) {
    //         m_PageSize = DEFAULT_NISAM_SIZE;
    //     } else {
    //         m_PageSize = DEFAULT_SISAM_SIZE;
    //     }
    //     if (eNoError !=x_InitSearch()) {
    //         m_Initialized = false;
    //     }
    //     x_FindIndexBounds();
    // }
    // ```
    /// The constructor: opens the two files. A file that cannot be opened is NCBI's "Error: Could
    /// not open input file (<index>/<data>)".
    pub fn open(index_path: &Path, data_path: &Path) -> Result<Self, SeqDbError> {
        let message = || {
            SeqDbError::NotFound(format!(
                "Error: Could not open input file ({}/{})",
                index_path.display(),
                data_path.display()
            ))
        };
        let index = DbFile::open(index_path).map_err(|_| message())?;
        let data = DbFile::open(data_path).map_err(|_| message())?;
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
    // ```c++
    // CSeqDBIsam::EErrorCode
    // CSeqDBIsam::x_InitSearch(void)
    // {
    //     if(m_Initialized == true)
    //         return eNoError;
    //
    //     TIndx info_needed = 10 * sizeof(Int4);
    //
    //     bool found_index_file =
    //         m_Atlas.GetFileSizeL(m_IndexFname, m_IndexFileLength);
    //
    //     if ((! found_index_file) || (m_IndexFileLength < info_needed)) {
    //         return eWrongFile;
    //     }
    //
    //
    //     Int4 * FileInfo = (Int4*)m_IndexLease.GetFileDataPtr(m_IndexFname,0);
    //
    //     // Check for consistence of files and parameters
    //
    //     Int4 Version = SeqDB_GetStdOrd(& FileInfo[0]);
    //
    //     if (Version != ISAM_VERSION)
    //         return eBadVersion;
    //
    //     Int4 IsamType = SeqDB_GetStdOrd(& FileInfo[1]);
    //
    //     if (IsamType == eNumericLongId && m_Type == eNumeric) {
    //         m_LongId = true;
    //         m_TermSize = 12;
    //         IsamType = eNumeric;
    //     }
    //
    //     if (IsamType != m_Type)
    //         return eBadType;
    //
    //     m_NumTerms    = SeqDB_GetStdOrd(& FileInfo[3]);
    //     m_NumSamples  = SeqDB_GetStdOrd(& FileInfo[4]);
    //     m_PageSize    = SeqDB_GetStdOrd(& FileInfo[5]);
    //     m_MaxLineSize = SeqDB_GetStdOrd(& FileInfo[6]);
    //
    //     if(m_PageSize != MEMORY_ONLY_PAGE_SIZE) {
    //         // Special case of memory-only index
    //         m_DataFileLength = SeqDB_GetStdOrd(& FileInfo[2]);
    //
    //         TIndx disk_file_length(0);
    //         bool found_data_file =
    //             m_Atlas.GetFileSizeL(m_DataFname, disk_file_length);
    //
    //         if ((! found_data_file) || (m_DataFileLength != disk_file_length)) {
    //             return eWrongFile;
    //         }
    //     }
    //
    //     // This space reserved for future use
    //
    //     m_IdxOption = SeqDB_GetStdOrd(& FileInfo[7]);
    //
    //     m_KeySampleOffset = (9 * sizeof(Int4));
    //
    //     m_Initialized = true;
    //
    //     return eNoError;
    // }
    // ```
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
    pub fn num_terms(&self) -> i32 {
        self.num_terms
    }

    /// The number of pages (`m_NumSamples`).
    pub fn num_samples(&self) -> i32 {
        self.num_samples
    }

    /// The lines per page (`m_PageSize`).
    pub fn page_size(&self) -> i32 {
        self.page_size
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1646-1665
    // ```c++
    // void CSeqDBIsam::GetIdBounds(string         & low_id,
    //                              string         & high_id,
    //                              int            & count)
    //
    //
    // {
    //     if(m_Initialized == false) {
    //         count = 0;
    //         return;
    //     }
    //
    //     if (! (m_FirstKey.IsSet() && m_LastKey.IsSet())) {
    //         count = 0;
    //         return;
    //     }
    //
    //     low_id = m_FirstKey.GetString();
    //     high_id = m_LastKey.GetString();
    //     count = m_NumTerms;
    // }
    // ```
    /// `GetIdBounds(string&, string&, int&)`: the lower-cased first and last key and the number
    /// of lines; count 0 and empty keys for an uninitialised index.
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
        self.index.read_exact_at(&mut bytes, at)?;
        Ok(u32::from_be_bytes(bytes))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:899-923
    // ```c++
    // void CSeqDBIsam::x_LoadPage(TIndx             SampleNum1,
    //                             TIndx             SampleNum2,
    //                             const char     ** beginp,
    //                             const char     ** endp)
    //
    // {
    //     // Load the appropriate page of terms into memory.
    //
    //     _ASSERT(SampleNum2 > SampleNum1);
    //
    //     TIndx begin_offset = m_KeySampleOffset + SampleNum1       * sizeof(Uint4);
    //     //TIndx end_offset   = m_KeySampleOffset + (SampleNum2 + 1) * sizeof(Uint4);
    //
    //     Uint4 * key_offsets((Uint4*)m_IndexLease.GetFileDataPtr(begin_offset));
    //
    //
    //     Uint4 key_off1 = SeqDB_GetStdOrd(& key_offsets[0]);
    //     Uint4 key_off2 = SeqDB_GetStdOrd(& key_offsets[SampleNum2 - SampleNum1]);
    //
    //     *beginp = (const char *) m_DataLease.GetFileDataPtr(m_DataFname,key_off1);
    //     *endp   = (const char *) m_DataLease.GetFileDataPtr(key_off2);
    // }
    //
    //
    // // ------------------------StringSearch--------------------------
    // ```
    /// `x_LoadPage(SampleNum1, SampleNum2)`: the data file's bytes from the first record of page
    /// `s1` to the first of page `s2`.
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
    // ```c++
    // void CSeqDBIsam::x_FindIndexBounds()
    // {
    //     Int4 Start (0);
    //     Int4 Stop  (m_NumSamples - 1);
    //
    //     //m_Atlas.Lock(locked);
    //
    //
    //     if (m_Type == eNumeric) {
    //         //
    //         // Get first key from data file
    //
    //         int num_elements(0);
    //         int start(0);
    //         const void * data_page(0);
    //
    //         x_MapDataPage(Start,
    //                       start,
    //                       num_elements,
    //                       & data_page);
    //
    //
    //         _ASSERT(num_elements);
    //
    //         int elem_index = 0;
    //
    //         Int8 data_gi(0);
    //         int data_oid(-1);
    //
    //         x_GetDataElement(data_page,
    //                          elem_index,
    //                          data_gi,
    //                          data_oid);
    //
    //         m_FirstKey.SetNumeric(data_gi);
    //
    //
    //         //
    //         // Get last key from data file
    //
    //         x_MapDataPage(Stop,
    //                       start,
    //                       num_elements,
    //                       & data_page);
    //
    //
    //         _ASSERT(num_elements);
    //
    //         elem_index = num_elements - 1;
    //
    //         x_GetDataElement(data_page,
    //                          elem_index,
    //                          data_gi,
    //                          data_oid);
    //
    //         m_LastKey.SetNumeric(data_gi);
    //     } else {
    //         //
    //         // Load the appropriate page of terms into memory.
    //
    //         const char * beginp(0);
    //         const char * endp(0);
    //
    //         //
    //         // Load the first page
    //
    //         x_LoadPage(Start, Start + 1, & beginp, & endp);
    //
    //         // Get first term
    //
    //         vector<string> keys_out;
    //         vector<string> data_out; // not used
    //
    //         x_ExtractData(beginp,
    //                       endp,
    //                       keys_out,
    //                       data_out);
    //
    //         x_Lower(keys_out.front());
    //         m_FirstKey.SetString(keys_out.front());
    //
    //
    //         //
    //         // Load the last page
    //
    //         x_LoadPage(Stop, Stop + 1, & beginp, & endp);
    //
    //         // Advance to last item
    //
    //         const char * lastp(0);
    //         const char * indexp(beginp);
    //
    //         while (indexp < endp) {
    //             // Remember our new "last term" value.
    //
    //             lastp = indexp;
    //
    //             // Skip remainder of term, and any nulls after it.
    //
    //             while((indexp < endp) && s_SeqDBIsam_NullifyEOLs(*indexp)) {
    //                 indexp++;
    //             }
    //             while((indexp < endp) && (! s_SeqDBIsam_NullifyEOLs(*indexp))) {
    //                 indexp++;
    //             }
    //         }
    //
    //         // Get the last key
    //
    //         _ASSERT(lastp);
    //
    //         keys_out.clear();
    //         data_out.clear();
    //
    //         x_ExtractData(lastp,
    //                       endp,
    //                       keys_out,
    //                       data_out);
    //
    //         x_Lower(keys_out.front());
    //         m_LastKey.SetString(keys_out.front());
    //     }
    // }
    // ```
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
    // ```c++
    // bool CSeqDBIsam::x_OutOfBounds(string key)
    // {
    //     if (! (m_FirstKey.IsSet() && m_LastKey.IsSet())) {
    //         return false;
    //     }
    //
    //     _ASSERT(m_Type == eString);
    //
    //     x_Lower(key);
    //
    //     if (m_FirstKey.OutsideFirstBound(key)) {
    //         return true;
    //     }
    //
    //     if (m_LastKey.OutsideLastBound(key)) {
    //         return true;
    //     }
    //
    //     return false;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1142-1148
    // ```c++
    //     /// Converts a string to lower case.
    //     static void x_Lower(string & s)
    //     {
    //         for(size_t i = 0; i < s.size(); i++) {
    //             s[i] = tolower(s[i]);
    //         }
    //     }
    // ```
    /// `x_OutOfBounds(string)`: the lower-cased `key` is before the first key or after the last
    /// one.
    fn out_of_bounds(&self, key: &[u8]) -> bool {
        let (Some(first), Some(last)) = (&self.first_key, &self.last_key) else {
            return false;
        };
        let key = key.to_ascii_lowercase();
        key < *first || key > *last
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:516-558
    // ```c++
    // int CSeqDBIsam::x_DiffCharLease(const string   & term_in,
    //                                 CSeqDBFileMemMap & lease,
    //                                 const string   & file_name,
    //                                 TIndx            file_length,
    //                                 Uint4            at_least,
    //                                 TIndx            KeyOffset,
    //                                 bool             ignore_case)
    //
    // {
    //     int result(-1);
    //
    //     //m_Atlas.Lock(locked);
    //
    //     // Add one to term_end to insure we don't consider "AA" and "AAB"
    //     // as equal.
    //
    //     TIndx offset_begin = KeyOffset;
    //     TIndx term_end     = KeyOffset + term_in.size() + 1;
    //     TIndx map_end      = term_end + at_least;
    //
    //     if (map_end > file_length) {
    //         map_end = file_length;
    //
    //         if (term_end > map_end) {
    //             term_end = map_end;
    //             result = int(file_length - offset_begin);
    //         }
    //     }
    //
    //     const char * file_data = (const char *)lease.GetFileDataPtr(file_name,offset_begin);
    //
    //     Int4 dc_result =
    //         x_DiffChar(term_in,
    //                    file_data,
    //                    file_data + term_in.size() + 1,
    //                    ignore_case);
    //
    //     if (dc_result != -1) {
    //         return dc_result;
    //     }
    //
    //     return result;
    // }
    // ```
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
        self.index
            .read_exact_at(&mut window[..readable], key_offset)?;
        let dc_result = diff_char(term_in, &window, true);
        if dc_result != -1 {
            return Ok((dc_result, window));
        }
        Ok((result, window))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:863-897
    // ```c++
    // int CSeqDBIsam::x_DiffSample(const string   & term_in,
    //                              Uint4            SampleNum,
    //                              TIndx          & KeyOffset)
    //
    // {
    //     // Meaning:
    //     // a. Compute SampleNum*4
    //     // b. Address this number into SamplePos (indexlease)
    //     // c. Swap this number to compute Key offset.
    //     // d. Add to beginning of file to get key data pointer.
    //
    //     bool ignore_case(true);
    //
    //     TIndx SampleOffset(m_KeySampleOffset);
    //
    //     if(m_PageSize != MEMORY_ONLY_PAGE_SIZE) {
    //         SampleOffset += (m_NumSamples + 1) * sizeof(Uint4);
    //     }
    //
    //     TIndx offset_begin = SampleOffset + (SampleNum * sizeof(Uint4));
    //     //TIndx offset_end   = offset_begin + sizeof(Uint4);
    //
    //     KeyOffset = SeqDB_GetStdOrd((Int4*)m_IndexLease.GetFileDataPtr(offset_begin));
    //
    //     Uint4 max_lines_2 = m_MaxLineSize * 2;
    //
    //     return x_DiffCharLease(term_in,
    //                            m_IndexLease,
    //                            m_IndexFname,
    //                            m_IndexFileLength,
    //                            max_lines_2,
    //                            KeyOffset,
    //                            ignore_case);
    //
    // }
    // ```
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
    // ```c++
    // void CSeqDBIsam::x_ExtractAllData(const string   & term_in,
    //                                   TIndx            sample_index,
    //                                   vector<TIndx>  & indices_out,
    //                                   vector<string> & keys_out,
    //                                   vector<string> & data_out)
    //
    // {
    //     // The object at sample_index is known to match; we will iterate
    //     // over the surrounding values to see if they match as well.  No
    //     // assumptions about how many keys can match are made here.
    //
    //     bool ignore_case = true;
    //
    //     int pre_amt  = 1;
    //     int post_amt = 1;
    //
    //     bool done_b(false), done_e(false);
    //
    //     const char * beginp(0);
    //     const char * endp(0);
    //
    //     TIndx beg_off(0);
    //     TIndx end_off(0);
    //
    //     while(! (done_b && done_e)) {
    //         if (sample_index < pre_amt) {
    //             beg_off = 0;
    //             done_b = true;
    //         } else {
    //             beg_off = sample_index - pre_amt;
    //         }
    //
    //         if ((m_NumSamples - sample_index) < post_amt) {
    //             end_off = m_NumSamples;
    //             done_e = true;
    //         } else {
    //             end_off = sample_index + post_amt;
    //         }
    //
    //         x_LoadPage(beg_off, end_off, & beginp, & endp);
    //
    //         if (! done_b) {
    //             Int4 diff_begin = x_DiffChar(term_in,
    //                                          beginp,
    //                                          endp,
    //                                          ignore_case);
    //
    //             if (diff_begin != -1) {
    //                 done_b = true;
    //             } else {
    //                 pre_amt ++;
    //             }
    //         }
    //
    //         if (! done_e) {
    //             const char * last_term(0);
    //             const char * p(endp-1);
    //
    //             // Skip over any non-terminating junk at the end
    //
    //             enum { eEndNulls, eLastTerm } search_stage = eEndNulls;
    //
    //             while(p > beginp) {
    //                 bool terminal = (0 == s_SeqDBIsam_NullifyEOLs(*p));
    //
    //                 if (search_stage == eEndNulls) {
    //                     if (! terminal) {
    //                         search_stage = eLastTerm;
    //                     }
    //                 } else {
    //                     if (terminal) {
    //                         last_term = p + 1;
    //                         break;
    //                     }
    //                 }
    //
    //                 p--;
    //             }
    //
    //             if (! last_term) {
    //                 last_term = beginp;
    //             }
    //
    //             Int4 diff_end = x_DiffChar(term_in,
    //                                        last_term,
    //                                        endp,
    //                                        ignore_case);
    //
    //             if (diff_end != -1) {
    //                 done_e = true;
    //             } else {
    //                 post_amt ++;
    //             }
    //         }
    //     }
    //
    //     x_ExtractPageData(term_in,
    //                       m_PageSize * beg_off,
    //                       beginp,
    //                       endp,
    //                       indices_out,
    //                       keys_out,
    //                       data_out);
    // }
    // ```
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
    // ```c++
    // CSeqDBIsam::EErrorCode
    // CSeqDBIsam::x_StringSearch(const string   & term_in,
    //                            vector<string> & terms_out,
    //                            vector<string> & values_out,
    //                            vector<TIndx>  & indices_out)
    //
    // {
    //     // These are always false; They may relate to the prior find_one /
    //     // expand_to_many method of getting multiple OIDs.
    //
    //     bool short_match(false);
    //     bool follow_match(false);
    //
    //     size_t preexisting_data_count = values_out.size();
    //
    //     if (m_Initialized == false) {
    //         return eInitFailed;
    //     }
    //
    //     if (x_OutOfBounds(term_in)) {
    //         return eNotFound;
    //     }
    //
    //     // We will set this option to avoid more complications
    //     bool ignore_case = true;
    //
    //     // search the sample file first
    //
    //     TIndx Start(0);
    //     TIndx Stop(m_NumSamples - 1);
    //
    //     int Length = (int) term_in.size();
    //
    //     TIndx SampleOffset(m_KeySampleOffset);
    //
    //     if(m_PageSize != MEMORY_ONLY_PAGE_SIZE) {
    //         SampleOffset += (m_NumSamples + 1) * sizeof(Uint4);
    //     }
    //
    //     int found_short(-1);
    //
    //     string short_term;
    //     int SampleNum(-1);
    //
    //     while(Stop >= Start) {
    //         SampleNum = ((Uint4)(Stop + Start)) >> 1;
    //
    //         TIndx KeyOffset(0);
    //
    //         int diff = x_DiffSample(term_in, SampleNum, KeyOffset);
    //
    //         // If this is an exact match, return the master term number.
    //
    //         const char * KeyData = (const char *)m_IndexLease.GetFileDataPtr(KeyOffset);
    //         TIndx BytesToEnd = m_IndexFileLength - KeyOffset;
    //
    //         Uint4 max_lines_2 = m_MaxLineSize * 2;
    //
    //         if (BytesToEnd > (TIndx) max_lines_2) {
    //             BytesToEnd = max_lines_2;
    //         }
    //
    //         if (diff == -1) {
    //             x_ExtractAllData(term_in,
    //                              SampleNum,
    //                              indices_out,
    //                              terms_out,
    //                              values_out);
    //
    //
    //             return eNoError;
    //         }
    //
    //         // If the key is a superset of the sample term, backup until
    //         // just before the term.
    //
    //         if (short_match && (diff >= Length)) {
    //             if (SampleNum > 0)
    //                 SampleNum--;
    //
    //             while(SampleNum > 0) {
    //                 TIndx key_offset =
    //                     x_GetIndexKeyOffset(SampleOffset,
    //                                         SampleNum);
    //
    //
    //                 string prefix;
    //                 x_GetIndexString(key_offset, Length, prefix, false);
    //
    //                 if (ignore_case) {
    //                     if (NStr::CompareNocase(prefix, term_in) != 0) {
    //                         break;
    //                     }
    //                 } else {
    //                     if (prefix != term_in) {
    //                         break;
    //                     }
    //                 }
    //
    //                 SampleNum--;
    //             }
    //
    //             found_short = SampleNum + 1;
    //
    //             TIndx key_offset =
    //                 x_GetIndexKeyOffset(SampleOffset,
    //                                     SampleNum + 1);
    //
    //
    //             string prefix;
    //             x_GetIndexString(key_offset, max_lines_2, short_term, true);
    //
    //             break;
    //         } else {
    //             // If preceding is desired, note the key.
    //
    //             if (follow_match) {
    //                 found_short = SampleNum;
    //
    //                 x_GetIndexString(KeyOffset, max_lines_2, short_term, true);
    //             }
    //         }
    //
    //         // Otherwise, search for the next sample.
    //
    //         if (ignore_case
    //             ? tolower((unsigned char) term_in[diff]) < tolower((unsigned char) KeyData[diff])
    //             : term_in[diff] < KeyData[diff]) {
    //             Stop = --SampleNum;
    //         } else {
    //             Start = SampleNum + 1;
    //         }
    //     }
    //
    //
    //     // If the term is out of range altogether, report not finding it.
    //
    //     if ( (SampleNum < 0) || (SampleNum >= m_NumSamples)) {
    //         return eNotFound;
    //     }
    //
    //     // Load the appropriate page of terms into memory.
    //
    //     const char * beginp(0);
    //     const char * endp(0);
    //
    //     x_LoadPage(SampleNum, SampleNum + 1, & beginp, & endp);
    //
    //     // Search the page for the term.
    //
    //     x_ExtractPageData(term_in,
    //                       m_PageSize * SampleNum,
    //                       beginp,
    //                       endp,
    //                       indices_out,
    //                       terms_out,
    //                       values_out);
    //
    //     // For now the short and follow logic is not implemented.
    //
    //     EErrorCode rv(eNoError);
    //
    //     if (preexisting_data_count == values_out.size()) {
    //         rv = eNotFound;
    //     }
    //
    //     return rv;
    // }
    // ```
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
    pub fn string_search(&self, term_in: &[u8]) -> Result<Vec<IsamLine>, SeqDbError> {
        let mut lines = Vec::new();
        self.string_search_into(term_in, &mut lines)?;
        Ok(lines)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1236-1376
    // ```c++
    // void CSeqDBIsam::StringToOids(const string   & acc,
    //                               vector<TOid>   & oids,
    //                               bool             adjusted,
    //                               bool           & version_check)
    //
    // {
    //     bool strip_version = version_check;
    //     version_check = false;
    //
    //     _ASSERT(m_IdentType == eStringId);
    //
    //     if(m_Initialized == false) {
    //             return;
    //     }
    //
    //     bool found = false;
    //
    //     string accession(string("gb|") + acc + "|");
    //     string locus_str(string("gb||") + acc);
    //
    //     EErrorCode err = eNoError;
    //
    //     vector<string> keys_out;
    //     vector<string> data_out;
    //     vector<TIndx>  indices_out;
    //
    //     if (! adjusted) {
    //         if ((err = x_StringSearch(accession,
    //                                   keys_out,
    //                                   data_out,
    //                                   indices_out)) < 0) {
    //             return;
    //         }
    //
    //         if (err == eNoError) {
    //             found = true;
    //         }
    //
    //         if ((! found) &&
    //             (err = x_StringSearch(locus_str,
    //                                   keys_out,
    //                                   data_out,
    //                                   indices_out)) < 0) {
    //
    //             return;
    //         }
    //
    //         if (err != eNotFound) {
    //             found = true;
    //         }
    //     }
    //
    //     if ((! found) &&
    //         (err = x_StringSearch(acc,
    //                               keys_out,
    //                               data_out,
    //                               indices_out)) < 0) {
    //
    //
    //         return;
    //     }
    //
    //     if (err != eNotFound) {
    //         found = true;
    //     }
    //
    //     if ((! found) && strip_version) {
    //         size_t pos = acc.find(".");
    //
    //         bool is_version = false;
    //
    //         if (pos != string::npos) {
    //             int ver_len = static_cast<int>(acc.size() - pos) - 1;
    //
    //             is_version = (ver_len <= 3 && ver_len >= 1);
    //
    //             for(size_t vp = pos+1; vp < acc.size(); vp++) {
    //                 if (! isdigit(acc[vp])) {
    //                     is_version = false;
    //                     break;
    //                 }
    //             }
    //         }
    //
    //         if (is_version) {
    //             string nover(acc, 0, pos);
    //
    //             err = x_StringSearch(nover,
    //                                  keys_out,
    //                                  data_out,
    //                                  indices_out);
    //
    //
    //             if (data_out.size()) {
    //                 version_check = true;
    //             }
    //
    //             if (err < 0) {
    //                 return;
    //             }
    //         }
    //     }
    //
    //     if (err != eNotFound) {
    //         found = true;
    //     }
    //
    //     if (! found) {
    //         // Use CSeq_id to parse the id string and build a replacement,
    //         // FASTA type string.  This allows some IDs, such as PDBs with
    //         // chains, such as '1qcfA' to be parsed.
    //
    //         string id;
    //
    //         try {
    //             CSeq_id seqid(acc, CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal);
    //             id = seqid.AsFastaString();
    //         }
    //         catch(CSeqIdException &) {
    //         }
    //
    //         if (id.size() &&
    //             ((err = x_StringSearch(id,
    //                                    keys_out,
    //                                    data_out,
    //                                    indices_out)) < 0)) {
    //
    //             return;
    //         }
    //     }
    //
    //     if (err != eNotFound) {
    //         found = true;
    //     }
    //
    //     if (found) {
    //         ITERATE(vector<string>, iter, data_out) {
    //             oids.push_back(atoi((*iter).c_str()));
    //         }
    //     }
    // }
    // ```
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
                let oid = atoi(&line.data);
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
    // ```c++
    // void CSeqDBVol::x_StringToOids(const string          & acc,
    //                                ESeqDBIdType            id_type,
    //                                Int8                    ident,
    //                                const string          & str_id,
    //                                bool                    simpler,
    //                                vector<int>           & oids) const
    // {
    //     bool vcheck (false);
    //     bool fits_in_four = (ident == -1) || ! (ident >> 32);
    //     bool needs_four = true;
    //
    //     switch(id_type) {
    //     case eStringId:
    //         x_OpenStrFile();
    //         if (! m_IsamStr.Empty()) {
    //             // Not simplified
    //             vcheck = true;
    //             m_IsamStr->StringToOids(str_id, oids, simpler, vcheck);
    //             x_UnleaseStrFile();
    //         }
    //         break;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2724-2733
    // ```c++
    //     if ((! fits_in_four) && needs_four) {
    //         NCBI_THROW(CSeqDBException,
    //                    eArgErr,
    //                    "ID overflows range of specified type.");
    //     }
    //
    //     if (vcheck) {
    //         x_CheckVersions(acc, oids);
    //     }
    // }
    // ```
    /// The string case of `CSeqDBVol::x_StringToOids`: `StringToOids` of the simplified id with
    /// `version_check` asked for, then `x_CheckVersions` of the user's `acc` when the lookup
    /// stripped a version. `seq_ids` is `GetSeqIDs(oid)` of the volume.
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
    // ```c++
    // template <> inline void
    // CSeqDBIsam::x_LoadIndex<string>(CSeqDBFileMemMap & lease,
    //                                 vector<string> & keys,
    //                                 vector<TIndx>  & offs)
    // {
    //     TIndx offset_begin = m_KeySampleOffset;
    //     TIndx sample_begin = offset_begin + sizeof(Uint4) * (m_NumSamples + 1);
    //
    //     // load offset array
    //     const Uint4 * offset = (const Uint4 *) lease.GetFileDataPtr(offset_begin);
    //     for (int index=0; index <= m_NumSamples; ++index, ++offset) {
    //         // Get the data_offsets
    //         offs.push_back(SeqDB_GetStdOrd((Uint4*) offset));
    //     }
    //
    //     // load sample array
    //     offset = (const Uint4 *) lease.GetFileDataPtr(sample_begin);
    //     for (int index=0; index < m_NumSamples; ++index, ++offset) {
    //         // Get the index_offsets
    //         offset_begin = SeqDB_GetStdOrd((Uint4*) offset);
    //
    //         // Lookup the samples
    //         const char * keydatap =  (const char *) lease.GetFileDataPtr(offset_begin) - 1;
    //
    //
    //         const char * key_begin = ++ keydatap;
    //         while (*keydatap != 0x02) ++keydatap;
    //         keys.push_back(string(key_begin, keydatap));
    //
    //         /* key_begin = ++keydatap;
    //         while (*keydatap != 0x00) ++keydatap;
    //         vals.push_back(NStr::StringToUInt(string(key_begin, keydatap))); */
    //     }
    // }
    // ```
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
    // ```c++
    // template <> inline void
    // CSeqDBIsam::x_LoadData<string>(CSeqDBFileMemMap & lease,
    //                                vector<string> & keys,
    //                                vector<int>    & vals,
    //                                int              num_keys,
    //                                TIndx            begin)
    // {
    //     const char * keydatap = (const char *) lease.GetFileDataPtr(begin) - 1;
    //
    //     for (int index=0; index < num_keys; ++index) {
    //
    //         const char * key_begin = ++keydatap;
    //         while (*keydatap != 0x02) ++keydatap;
    //         keys.push_back(string(key_begin, keydatap));
    //
    //         key_begin = ++keydatap;
    //         while (*keydatap != 0x0a) ++keydatap;
    //         vals.push_back(NStr::StringToUInt(string(key_begin, keydatap)));
    //     }
    // }
    // ```
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
        let mut keys = Vec::with_capacity(num_keys);
        let mut vals = Vec::with_capacity(num_keys);
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
    // ```c++
    //     if (ids.GetNumSis() && (GetLMDBFileName() == kEmptyStr)) {
    //         x_OpenStrFile();
    //         if (m_IsamStr.NotEmpty()) {
    //             m_IsamStr->IdsToOids(m_VolStart, m_VolEnd, ids);
    //             x_UnleaseStrFile();
    //         } else {
    //             NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "SI list specified but no ISAM file found for SI in " + m_VolName);
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1388-1419
    // ```c++
    // void CSeqDBIsam::IdsToOids(int              vol_start,
    //                            int              vol_end,
    //                            CSeqDBGiList   & ids)
    //
    //
    // {
    //     // The vol_start parameter is needed because translations in the
    //     // GI list should refer to global OIDs, not per-volume OIDs.
    //
    //     switch (m_IdentType) {
    //     case eGiId:
    //         x_TranslateGiList<TGi>(vol_start, ids);
    //         break;
    //
    //     case eTiId:
    //         x_TranslateGiList<TTi>(vol_start, ids);
    //         break;
    //
    //     case eStringId:
    //         x_TranslateGiList<string>(vol_start, ids);
    //         break;
    //
    //     case ePigId:
    //         x_TranslateGiList<TPig>(vol_start, ids);
    //         break;
    //
    //     default:
    //         NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "Error: Wrong type of idlist specified.");
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:548-623
    // ```c++
    //     template <class T>
    //     void x_TranslateGiList(int              vol_start,
    //                            CSeqDBGiList   & gis)
    //
    //     {
    //         int gilist_size = gis.GetSize<T>();
    //         if (! gilist_size) return;
    //
    //         gis.InsureOrder(CSeqDBGiList::eGi);
    //
    //         if(m_Initialized == false) {
    //                 NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "Error: Unable to use ISAM index in batch mode.");
    //         }
    //
    //         vector<T> sample_keys;
    //         vector<TIndx> page_offs;
    //         vector<T> keys;
    //         vector<int> vals;
    //
    //         sample_keys.reserve(m_NumSamples);
    //         page_offs.reserve(m_NumSamples + 1);
    //         keys.reserve(m_PageSize);
    //         vals.reserve(m_PageSize);
    //
    //         x_LoadIndex(m_IndexLease, sample_keys, page_offs);
    //
    //         int gilist_index = 0;
    //         int sample_index = 0;
    //
    //         while((gilist_index < gilist_size) && (sample_index < m_NumSamples)) {
    //
    //             s_AdvanceGiList<T>(gis, gilist_index, gilist_size,
    //                                sample_keys[sample_index]);
    //
    //             if (gilist_index >= gilist_size) break;
    //
    //             s_AdvanceKeyList<T>(sample_keys, sample_index, m_NumSamples,
    //                                 gis.GetKey<T>(gilist_index));
    //
    //             // Now we should be ready to search a data block.
    //             keys.clear();
    //             vals.clear();
    //
    //             int num_keys = m_PageSize;
    //             if (sample_index + 1 == m_NumSamples) {
    //                 num_keys = m_NumTerms - sample_index * m_PageSize;
    //             }
    //             x_LoadData(m_DataLease, keys, vals, num_keys, page_offs[sample_index]);
    //
    //             int index = 0;
    //
    //             while ((gilist_index < gilist_size) && (index < num_keys)) {
    //
    //                 s_AdvanceKeyList<T>(keys, index, num_keys,
    //                                     gis.GetKey<T>(gilist_index));
    //
    //                 s_SetTranslation<T>(gis, gilist_index, gilist_size,
    //                                     keys[index], vals[index] + vol_start);
    //
    //                 ++index;
    //                 if (index >= num_keys) break;
    //
    //                 s_AdvanceGiList<T>(gis, gilist_index, gilist_size, keys[index]);
    //
    //                 s_SetTranslation<T>(gis, gilist_index, gilist_size,
    //                                     keys[index], vals[index] + vol_start);
    //
    //             }
    //
    //             // We could be finished here because we exhausted the GI list
    //             // We must be done with that one by now..
    //             ++sample_index;
    //         }
    //     }
    // ```
    /// `IdsToOids(vol_start, vol_end, CSeqDBGiList&)` for the string list
    /// (`x_TranslateGiList<string>`): sets `oid` of every entry of `ids` whose `si` (lower-cased,
    /// simplified: [`simplify_accession_for_isam`]) is a key of the volume, to the OID of the
    /// line plus `vol_start`. An entry that has an `oid` already is left as it is.
    pub fn ids_to_oids(&self, vol_start: u32, ids: &mut [SiOid]) -> Result<(), SeqDbError> {
        let gilist_size = ids.len();
        if gilist_size == 0 {
            return Ok(());
        }
        insure_order(ids);
        if !self.initialized {
            return Err(batch_mode_error());
        }
        let (page_offs, sample_keys) = self.load_index(true)?;
        let num_samples = self.num_samples as usize;
        let order_error = || SeqDbError::Format {
            path: self.data.path().to_path_buf(),
            message: "string ISAM keys are not in order".to_string(),
        };
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
    // ```c++
    //     if (ids.GetNumSis()) {
    //         x_OpenStrFile();
    //         if (m_IsamStr.NotEmpty()) {
    //             m_IsamStr->IdsToOids(m_VolStart, m_VolEnd, ids);
    //             x_UnleaseStrFile();
    //         } else {
    //             NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "SI list specified but no ISAM file found for SI in " + m_VolName);
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1421-1459
    // ```c++
    // void CSeqDBIsam::IdsToOids(int                  vol_start,
    //                            int                  vol_end,
    //                            CSeqDBNegativeList & ids)
    //
    //
    // {
    //     // The vol_start parameter is needed because translations in the
    //     // GI list should refer to global OIDs, not per-volume OIDs.
    //
    //     _ASSERT(m_IdentType == eGiId || m_IdentType == eTiId || m_IdentType == eStringId);
    //
    //     //m_Atlas.Lock(locked);
    //
    //     ids.InsureOrder();
    //
    //     if ((m_IdentType == eGiId) && ids.GetNumGis()) {
    //         x_SearchNegativeMulti(vol_start,
    //                               vol_end,
    //                               ids,
    //                               false);
    //
    //     }
    //
    //     if ((m_IdentType == eTiId) && ids.GetNumTis()) {
    //         x_SearchNegativeMulti(vol_start,
    //                               vol_end,
    //                               ids,
    //                               true);
    //
    //     }
    //
    //     if(m_IdentType == eStringId && ids.GetNumSis()) {
    //         x_SearchNegativeMultiSeq(vol_start,
    //                               vol_end,
    //                               ids);
    //                               //true,
    //
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:333-418
    // ```c++
    // CSeqDBIsam::x_SearchNegativeMultiSeq(int              vol_start,
    //                                      int              vol_end,
    //                                      CSeqDBNegativeList & ids)
    //
    // {
    //         int gilist_size = ids.ListSize();
    //         if (! gilist_size) return;
    //
    //         if(m_Initialized == false) {
    //             // Most ordinary errors (missing GIs for example) are
    //             // ignored for "multi" mode searches.  But if a GI list is
    //             // specified, and cannot be interpreted, it is an error.
    //
    //             NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "Error: Unable to use ISAM index in batch mode.");
    //         }
    //
    //
    //         vector<string> sample_keys;
    //         vector<TIndx> page_offs;
    //         vector<string> keys;
    //         vector<int> vals;
    //
    //         sample_keys.reserve(m_NumSamples);
    //         page_offs.reserve(m_NumSamples + 1);
    //         keys.reserve(m_PageSize);
    //         vals.reserve(m_PageSize);
    //
    //
    //         x_LoadIndex(m_IndexLease, sample_keys, page_offs);
    //
    //         int gilist_index = 0;
    //         int sample_index = 0;
    //
    //         while(sample_index < m_NumSamples) {
    //
    //             // Now we should be ready to search a data block.
    //             keys.clear();
    //             vals.clear();
    //
    //             int num_keys = m_PageSize;
    //             if (sample_index + 1 == m_NumSamples) {
    //                 num_keys = m_NumTerms - sample_index * m_PageSize;
    //             }
    //
    //             x_LoadData(m_DataLease, keys, vals, num_keys, page_offs[sample_index]);
    //
    //             for(int i = 0; i < num_keys; i++) {
    //                 // 2. Look for it in the negative id list.
    //
    //                 bool found = false;
    //                 if (gilist_index < gilist_size) {
    //                     found = x_FindInNegativeList(ids,
    //                                              gilist_index,
    //                                              keys[i]);
    //
    //                 }
    //                 if (vals[i] < vol_end) {
    //                     if (found) {
    //                         // OID is found, but may not be included yet.
    //                         ids.AddVisibleOid(vals[i] + vol_start);
    //                         //If next accession is the same as current, but with version
    //                         if(s_IsSameAccession(keys, num_keys, i)) {
    //                             i++; //skip next check - sequence already excluded
    //                         }
    //                     } else {
    //                         // OID is included for iteration.
    //                         //Only include next accession if it is not the same as current (but with version)
    //                         //because it may be in exclude list. The check will be done in the next step
    //                         if(!s_IsSameAccession(keys, num_keys, i)) {
    //                             ids.AddIncludedOid(vals[i] + vol_start);
    //                         }
    //                     }
    //                 }
    //
    //             }
    //             // Move to next data page.  Note that for a negative ID list
    //             // processing, we don't actually fetch any samples, because
    //             // every ID->OID line needs to be examined anyway.
    //
    //             sample_index ++;
    //
    //         }
    // }
    //
    // ```
    /// `x_SearchNegativeMultiSeq` (through `IdsToOids(vol_start, vol_end, CSeqDBNegativeList&)`):
    /// reads every line of the volume and marks its OID (plus `vol_start`) in `out`: visible if
    /// its key is in `ids`, included if not. `ids` are the lower-cased simplified ids; they are
    /// sorted here (`CSeqDBNegativeList::InsureOrder`). `vol_end` is the global end OID of the
    /// volume (see the quirks above).
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

/// NCBI's `CSeqDBException(eArgErr, "Error: Unable to use ISAM index in batch mode.")`.
fn batch_mode_error() -> SeqDbError {
    SeqDbError::NotFound("Error: Unable to use ISAM index in batch mode.".to_string())
}

/// `atoi` of the C library (C locale): leading white space, a sign, the digits; the `long` it
/// reads, truncated to `int`.
fn atoi(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && (s[i] == b' ' || (b'\t'..=b'\r').contains(&s[i])) {
        i += 1;
    }
    let mut negative = false;
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        negative = s[i] == b'-';
        i += 1;
    }
    let mut value: i64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        value = value
            .saturating_mul(10)
            .saturating_add(i64::from(s[i] - b'0'));
        i += 1;
    }
    (if negative { -value } else { value }) as i32
}

/// `NStr::StringToUInt(str)` with no flags: an optional `+`, digits only, at most `UINT_MAX`.
fn string_to_uint(s: &[u8]) -> Option<u32> {
    let digits = s.strip_prefix(b"+").unwrap_or(s);
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    digits.iter().try_fold(0u32, |n, &c| {
        n.checked_mul(10)?.checked_add(u32::from(c - b'0'))
    })
}

/// One entry of a list of string ids (`CSeqDBGiList::SSiOid`): the id and its OID, set by
/// [`StringIsam::ids_to_oids`] (NCBI: -1 while unknown).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiOid {
    pub si: Vec<u8>,
    pub oid: Option<u32>,
}

impl SiOid {
    /// `SSiOid(si, -1)`.
    pub fn new(si: impl Into<Vec<u8>>) -> Self {
        Self {
            si: si.into(),
            oid: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:574-586
// ```c++
// /// Compare SSeqIdOid structs by SeqId.
// class CSeqDB_SortSiLessThan {
// public:
//     /// Test whether lhs is less than (occurs before) rhs.
//     /// @param lhs Left hand side of less-than operator. [in]
//     /// @param rhs Right hand side of less-than operator. [in]
//     /// @return True if lhs sorts before rhs by Seq-id.
//     int operator()(const CSeqDBGiList::SSiOid & lhs,
//                    const CSeqDBGiList::SSiOid & rhs)
//     {
//         return lhs.si < rhs.si;
//     }
// };
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:589-606
// ```c++
// template<class TCompare, class TVector>
// void s_InsureOrder(TVector & v)
// {
//     bool already = true;
//
//     TCompare compare_less;
//
//     for(int i = 1; i < (int) v.size(); i++) {
//         if (compare_less(v[i], v[i-1])) {
//             already = false;
//             break;
//         }
//     }
//
//     if (! already) {
//         sort(v.begin(), v.end(), compare_less);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:616-651
// ```c++
// void CSeqDBGiList::InsureOrder(ESortOrder order)
// {
//     // Code depends on OID order after translation, because various
//     // methods of SeqDB use this class for filtering purposes.
//     static CFastMutex mtx;
//     CFastMutexGuard mtx_gurad(mtx);
//     if ((order < m_CurrentOrder) || (order == eNone)) {
//         NCBI_THROW(CSeqDBException,
//                    eFileErr,
//                    "Out of sequence sort order requested.");
//     }
//
//     // Input is usually sorted by GI, so we first test for sortedness.
//     // If it will fail it will probably do so almost immediately.
//
//     if (order != m_CurrentOrder) {
//         switch(order) {
//         case eNone:
//             break;
//
//         case eGi:
//             s_InsureOrder<CSeqDB_SortGiLessThan>(m_GisOids);
//             s_InsureOrder<CSeqDB_SortTiLessThan>(m_TisOids);
//             s_InsureOrder<CSeqDB_SortSiLessThan>(m_SisOids);
//             s_InsureOrder<CSeqDB_SortPigLessThan>(m_PigsOids);
//             break;
//
//         default:
//             NCBI_THROW(CSeqDBException,
//                        eFileErr,
//                        "Unrecognized sort order requested.");
//         }
//
//         m_CurrentOrder = order;
//     }
// }
// ```
/// `CSeqDBGiList::InsureOrder(eGi)` for the string list: sorts by `si` (`std::string` order:
/// unsigned bytes, then length) unless it is sorted already.
pub fn insure_order(ids: &mut [SiOid]) {
    let already = ids.windows(2).all(|w| w[1].si >= w[0].si);
    if !already {
        ids.sort_by(|a, b| a.si.cmp(&b.si));
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:50-75
// ```c++
// // Use Parabolic Binary Search to find the first gilist that
// // is greater-than-or-equal-to the key
// template <class T> static inline void
// s_AdvanceGiList(CSeqDBGiList & gis,
//                 int          & gi_index,
//                 int            gis_size,
//                 const T      & key)
// {
//     while( (gi_index < gis_size)
//        &&  (gis.GetKey<T>(gi_index) < key)) {
//
//         ++gi_index;
//         int jump = 2;
//
//         while( (gi_index + jump < gis_size)
//            &&  (gis.GetKey<T>(gi_index + jump) < key)) {
//
//             gi_index += jump;
//             jump *= 2;
//         }
//     }
//
//     // skipping translated elements
//     while( (gi_index < gis_size)
//        &&  (gis.IsValueSet<T>(gi_index) )) ++gi_index;
// }
// ```
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
// ```c++
// // Use Parabolic Binary Search to find the largest sample that
// // is less-than-or-equal-to the first untranslated target GI.
// template <class T> static inline void
// s_AdvanceKeyList(const vector<T> & keys,
//                  int             & index,
//                  int               num_keys,
//                  const T         & target)
// {
//     while( (index < num_keys)
//        &&  (keys[index] <= target)) {
//
//         ++index;
//         int jump = 2;
//
//         while( (index + jump < num_keys)
//            &&  (keys[index + jump] <= target)) {
//
//             index += jump;
//             jump *= 2;
//         }
//     }
//     --index;
// }
// ```
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
// ```c++
// // apply the translation (if we have it) for those GIs.
// template <class T> static inline void
// s_SetTranslation(CSeqDBGiList & gis,
//                  int          & gi_index,
//                  int            gis_size,
//                  const T      & key,
//                  int            value)
// {
//     while( (gi_index < gis_size)
//        &&  (gis.GetKey<T>(gi_index) == key)) {
//
//         gis.SetValue<T>(gi_index, value);
//         ++gi_index;
//     }
// }
// ```
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
// ```c++
// inline bool
// CSeqDBIsam::x_FindInNegativeList(CSeqDBNegativeList & ids,
//                                  int                & index,
//                                  string                 key)
// {
//     bool found = false;
//
//     // Skip any that are less than key.
//     int ids_size = ids.ListSize();
//
//     while((index < ids_size) && (x_GetId(ids, index) < key)) {
//         index++;
//
//         int jump = 2;
//
//         while((index + jump) < ids_size &&
//               x_GetId(ids, index + jump) < key) {
//             index += jump;
//             jump += jump;
//         }
//     }
//
//     // Check whether the GI or TI was found.
//
//     if ((index < ids_size) && (x_GetId(ids,index) == key)) {
//         found = true;
//     }
//
//     return found;
// }
// ```
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
// ```c++
// //In case if acc2 does not have version and acc1 has version
// //check if it is the same accession
// static bool s_IsSameAccession(string acc1, string acc2)
// {
//     bool sameAccession = false;
//     if(NStr::Find(acc2,".") == NPOS) { // no version in acc2
//         if(NStr::Find(acc1,".") != NPOS && NStr::Find(acc1, acc2) != NPOS) {
//             string accession, version;
//             NStr::SplitInTwo(acc1,".", accession, version);
//             if(acc2 == accession) {
//                 sameAccession = true;
//             }
//         }
//     }
//     return sameAccession;
// }
// ```
/// `s_IsSameAccession(acc1, acc2)`: `acc2` has no `.`, `acc1` has one and holds `acc2`, and
/// `acc2` is what `acc1` has before its first `.`.
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
// ```c++
// //In case if keys[currIndex] does not have version and keys[currIndex + 1] has version
// //check if it is the same accession
// static bool s_IsSameAccession(vector <string> keys, int num_keys, int currIndex)
// {
//     bool sameAccession = false;
//     if(currIndex < num_keys - 1) {
//         sameAccession = s_IsSameAccession(keys[currIndex + 1], keys[currIndex]);
//     }
//     return sameAccession;
// }
// ```
/// `s_IsSameAccession(keys, num_keys, currIndex)`: whether line `curr_index + 1` of the page is
/// line `curr_index` with a version.
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
        let need = (u64::from(num) + 31) / 32;
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
// ```c++
//     /// Include an OID in the iteration.
//     ///
//     /// The OID will be included by SeqDB in the set returned to users
//     /// by OID iteration.
//     ///
//     /// @param oid The OID in question. [in]
//     void AddIncludedOid(int oid)
//     {
//         m_Included.SetBit(oid);
//     }
//
//     /// Indicate a visible OID.
//     ///
//     /// The OID will be marked as having been found in a GI or TI
//     /// ISAM index (but possibly not included for iteration).
//     ///
//     /// @param oid The OID in question. [in]
//     void AddVisibleOid(int oid)
//     {
//         m_Visible.SetBit(oid);
//     }
//
//     /// Get the inclusion status of an OID.
//     ///
//     /// This returns true for OIDs that were in the included set and
//     /// for OIDs that were not found in the ISAM file at all.
//     ///
//     /// @param oid The OID in question. [in]
//     /// @return True if the OID is included by SeqDB.
//     bool GetOidStatus(int oid)
//     {
//         return m_Included.GetBit(oid) || (! m_Visible.GetBit(oid));
//     }
//
//     /// Get the size of the OID array.
//     int GetNumOids()
//     {
//         return max(m_Visible.Size(), m_Included.Size());
//     }
// ```
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

    /// `m_Included.GetBit(oid)`.
    pub fn is_included(&self, oid: u32) -> bool {
        self.included.get_bit(oid)
    }

    /// `m_Visible.GetBit(oid)`.
    pub fn is_visible(&self, oid: u32) -> bool {
        self.visible.get_bit(oid)
    }
}

/// `NStr::StringToInt(str, fConvErr_NoThrow | fAllowTrailingSymbols)`: an optional sign and
/// digits, anything after them ignored; 0 when there are no digits or the number does not fit.
fn string_to_int_trailing(s: &[u8]) -> i32 {
    let mut pos = 0;
    let mut sign = false;
    match s.first() {
        Some(b'-') => {
            sign = true;
            pos = 1;
        }
        Some(b'+') => pos = 1,
        _ => {}
    }
    let pos0 = pos;
    let limdiv = i64::MAX / 10;
    let limoff = i64::MAX % 10 + i64::from(sign);
    let mut n: i64 = 0;
    while pos < s.len() && s[pos].is_ascii_digit() {
        let delta = i64::from(s[pos] - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return 0;
        }
        n = n * 10 + delta;
        pos += 1;
    }
    if pos == pos0 {
        return 0;
    }
    let n = if sign { -n } else { n };
    i32::try_from(n).unwrap_or(0)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2735-2789
// ```c++
// void CSeqDBVol::x_CheckVersions(const string         & acc,
//                                 vector<int>          & oids) const
// {
//     // If we resolved a string id by stripping the version off of the
//     // string, we need to check (for each OID) if the real ID had a
//     // matching version.
//
//     // This condition happens in two cases: where the database has a
//     // different version of the same ID, and in the case of sparse
//     // databases (which do not store the version).  In the latter
//     // case, we may get a list of OIDs where some of the OIDs pass
//     // this test, and others fail.
//
//     size_t pos = acc.find(".");
//     _ASSERT(pos != acc.npos);
//
//     string ver_str(acc, pos+1, acc.size()-(pos+1));
//     int vernum = NStr::StringToInt(ver_str,
//                                    NStr::fConvErr_NoThrow |
//                                    NStr::fAllowTrailingSymbols);
//
//     string nover(acc, 0, pos);
//
//     size_t pos2(0);
//     while((pos2 = nover.find("|")) != nover.npos) {
//         nover.erase(0, pos2+1);
//     }
//
//     NON_CONST_ITERATE(vector<int>, iter, oids) {
//         list< CRef<CSeq_id> > ids =
//             GetSeqIDs(*iter);
//
//         bool found = false;
//
//         ITERATE(list< CRef<CSeq_id> >, seqid, ids) {
//             const CTextseq_id * id = (*seqid)->GetTextseq_Id();
//
//             if (id &&
//                 id->CanGetAccession() &&
//                 id->GetAccession() == nover &&
//                 id->CanGetVersion() &&
//                 id->GetVersion() == vernum) {
//
//                 found = true;
//                 break;
//             }
//         }
//
//         if (! found) {
//             *iter = -1;
//         }
//     }
//
//     oids.erase(remove(oids.begin(), oids.end(), -1), oids.end());
// }
// ```
/// `CSeqDBVol::x_CheckVersions`: removes the OIDs whose Seq-ids (`seq_ids` is `GetSeqIDs(oid)`
/// of the volume) hold no Textseq-id with the accession and version of the user's `acc`.
pub fn check_versions(
    acc: &[u8],
    oids: &mut Vec<u32>,
    seq_ids: &mut dyn FnMut(u32) -> Result<Vec<ParsedSeqId>, SeqDbError>,
) -> Result<(), SeqDbError> {
    // If we resolved a string id by stripping the version off of the string, we need to check
    // (for each OID) if the real ID had a matching version.
    let pos = acc.iter().position(|&c| c == b'.');
    // `acc.find(".")` is asserted to exist; without it `pos + 1` wraps to 0 and `nover` is `acc`.
    let ver_str = &acc[pos.map_or(0, |p| p + 1)..];
    let vernum = string_to_int_trailing(ver_str);
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
// ```c++
// enum ESeqDBIdType {
//     eGiId,     /// Genomic ID is a relatively stable numeric identifier for sequences.
//     eTiId,     /// Trace ID is a numeric identifier for Trace sequences.
//     ePigId,    /// Each PIG identifier refers to exactly one protein sequence.
//     eStringId, /// Some sequence sources uses string identifiers.
//     eHashId,   /// Lookup from sequence hash values to OIDs.
//     eOID       /// The ordinal id indicates the order of the data in the volume's index file.
// };
// ```
/// `ESeqDBIdType`, without `eHashId` (a Seq-id is never simplified to a hash).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqDbIdType {
    Gi,
    Ti,
    Pig,
    String,
    Oid,
}

/// The results of `SeqDB_SimplifyAccession(acc, num_id, str_id, simpler)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SimplifiedAccession {
    /// The returned `ESeqDBIdType`.
    pub id_type: SeqDbIdType,
    /// `num_id`: the GI, TI, PIG or OID (`(Uint4)-1` when none was set).
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
    /// A `CSeqIdException` of `CSeq_id(element)` or of the `db:tag` id.
    SeqId(SeqIdError),
    /// Another exception (LOSAT's text): a string tag read as a number, an integer read from a
    /// string tag.
    Exception(String),
}

impl std::fmt::Display for SimplifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimplifyError::SeqId(error) => f.write_str(&error.what()),
            SimplifyError::Exception(message) => f.write_str(message),
        }
    }
}

impl std::error::Error for SimplifyError {}

/// `NStr::StringToInt8(str)` with no flags: an optional sign and digits, nothing else.
fn string_to_int8(s: &[u8]) -> Option<i64> {
    let (negative, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let mut n: i128 = 0;
    for &c in digits {
        n = n * 10 + i128::from(c - b'0');
        if n > i128::from(i64::MAX) + 1 {
            return None;
        }
    }
    let n = if negative { -n } else { n };
    i64::try_from(n).ok()
}

/// `CObject_id::GetId()`.
fn object_id_id(tag: &ObjectId) -> Result<i64, SimplifyError> {
    match tag {
        ObjectId::Id(n) => Ok(i64::from(*n)),
        ObjectId::Str(_) => Err(SimplifyError::Exception(
            "Object-id::GetId: wrong choice selection (str)".to_string(),
        )),
    }
}

fn fasta_string(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    label::fasta(id, &mut out);
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2264-2427
// ```c++
// ESeqDBIdType SeqDB_SimplifySeqid(CSeq_id       & bestid,
//                                  const string  * acc,
//                                  Int8          & num_id,
//                                  string        & str_id,
//                                  bool          & simpler)
// {
//     ESeqDBIdType result = eStringId;
//
//     const CTextseq_id * tsip = 0;
//
//     bool matched = true;
//
//     switch(bestid.Which()) {
//     case CSeq_id::e_Gi:
//         simpler = true;
//         num_id = GI_TO(Int8, bestid.GetGi());
//         result = eGiId;
//         break;
//
//     case CSeq_id::e_Gibbsq:    /* gibbseq */
//         simpler = true;
//         result = eStringId;
//         str_id = NStr::UIntToString(bestid.GetGibbsq());
//         break;
//
//     case CSeq_id::e_General:
//         {
//             const CDbtag & dbt = bestid.GetGeneral();
//
//             if (dbt.CanGetDb()) {
//                 if (dbt.GetDb() == "BL_ORD_ID") {
//                     simpler = true;
//                     num_id = dbt.GetTag().GetId();
//                     result = eOID;
//                     break;
//                 }
//
//                 if (dbt.GetDb() == "PIG") {
//                     simpler = true;
//                     num_id = dbt.GetTag().GetId();
//                     result = ePigId;
//                     break;
//                 }
//
//                 if (dbt.GetDb() == "ti") {
//                     simpler = true;
//                     num_id = (dbt.GetTag().IsStr()
//                               ? NStr::StringToInt8(dbt.GetTag().GetStr())
//                               : dbt.GetTag().GetId());
//
//                     result = eTiId;
//                     break;
//                 }
//
//
//                 if (NStr::CompareNocase(dbt.GetDb(), "GNOMON") == 0) {
//                     str_id = bestid.AsFastaString();
//                     str_id = NStr::ToLower(str_id);
//                     result = eStringId;
//                     break;
//                 }
//             }
//
//             if (dbt.CanGetTag() && dbt.GetTag().IsStr()) {
//                 result = eStringId;
//                 str_id = dbt.GetTag().GetStr();
//                 str_id = NStr::ToLower(str_id);
//             } else {
//                 // Use the default logic.
//                 matched = false;
//             }
//         }
//         break;
//
//     case CSeq_id::e_Local:     /* local */
//         simpler = true;
//         result = eStringId;
//         {
//             const CObject_id & objid = bestid.GetLocal();
//
//             if (objid.IsStr()) {
//                 // sparse version will leave "lcl|" off.
//                 str_id = objid.GetStr();
//                 str_id = NStr::ToLower(str_id);
//             } else {
//                 // Local numeric ids are stored as strings.
//                 str_id = "lcl|" + NStr::IntToString(objid.GetId());
//             }
//         }
//         break;
//
//         // tsip types
//
//     case CSeq_id::e_Embl:      /* embl */
//     case CSeq_id::e_Ddbj:      /* ddbj */
//     case CSeq_id::e_Genbank:   /* genbank */
//     case CSeq_id::e_Tpg:       /* Third Party Annot/Seq Genbank */
//     case CSeq_id::e_Tpe:       /* Third Party Annot/Seq EMBL */
//     case CSeq_id::e_Tpd:       /* Third Party Annot/Seq DDBJ */
//     case CSeq_id::e_Other:     /* other */
//     case CSeq_id::e_Swissprot: /* swissprot (now with versions) */
//     case CSeq_id::e_Gpipe:     /* internal NCBI genome pipeline */
//         tsip = bestid.GetTextseq_Id();
//         break;
//
//     case CSeq_id::e_Pir:       /* pir   */
//     case CSeq_id::e_Prf:       /* prf   */
//         tsip = bestid.GetTextseq_Id();
//         break;
//
//     default:
//         matched = false;
//     }
//
//     // Default: if we have a string, use it; if we only have seqid,
//     // create a string.  This should not happen if the seqid matches
//     // one of the cases above, which currently correspond to all the
//     // supported seqid types.
//
//     CSeq_id::ELabelFlags label_flags = (CSeq_id::ELabelFlags)
//         (CSeq_id::fLabel_GeneralDbIsContent | CSeq_id::fLabel_Version);
//
//     if (! matched) {
//         // (should not happen normally)
//
//         simpler = false;
//         result  = eStringId;
//
//         if (acc) {
//             str_id = *acc;
//             str_id = NStr::ToLower(str_id);
//         } else {
//             bestid.GetLabel(& str_id, CSeq_id::eFasta, label_flags);
//             str_id = NStr::ToLower(str_id);
//         }
//     }
//
//     if (tsip) {
//         bool found = false;
//
//         if (tsip->CanGetAccession()) {
//             str_id = tsip->GetAccession();
//             str_id = NStr::ToLower(str_id);
//             found = true;
//
//             if (tsip->CanGetVersion()) {
//                 str_id += ".";
//                 str_id += NStr::UIntToString(tsip->GetVersion());
//             }
//         } else if (tsip->CanGetName()) {
//             str_id = tsip->GetName();
//             str_id = NStr::ToLower(str_id);
//             found = true;
//         }
//
//         if (found) {
//             simpler = true;
//             result = eStringId;
//         }
//     }
//
//     return result;
// }
//
// ```
/// `SeqDB_SimplifySeqid(bestid, &acc, num_id, str_id, simpler)` (the call with `acc == NULL`,
/// made by `CSeqDBGiList::FindId`, is not ported).
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
                    ObjectId::Str(text) => string_to_int8(text).ok_or_else(|| {
                        SimplifyError::Exception(format!(
                            "Cannot convert string '{}' to Int8",
                            String::from_utf8_lossy(text)
                        ))
                    })?,
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

/// `WhichInverseSeqId(portion.c_str())`: the text up to its first NUL.
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
// ```c++
// static size_t
// s_SeqDB_EndOfFastaID(const string & str, size_t pos)
// {
//     // (Derived from s_EndOfFastaID()).
//
//     size_t vbar = str.find('|', pos);
//
//     if (vbar == string::npos) {
//         return string::npos; // bad
//     }
//
//     string portion(str, pos, vbar - pos);
//
//     CSeq_id::E_Choice choice =
//         CSeq_id::WhichInverseSeqId(portion.c_str());
//
//     if (choice != CSeq_id::e_not_set) {
//         size_t vbar_prev = vbar;
//         int count;
//         for (count=0; ; ++count, vbar_prev = vbar) {
//             vbar = str.find('|', vbar_prev + 1);
//
//             if (vbar == string::npos) {
//                 break;
//             }
//
//             int start_pt = int(vbar_prev + 1);
//             string element(str, start_pt, vbar - start_pt);
//
//             choice = CSeq_id::WhichInverseSeqId(element.c_str());
//
//             if (choice != CSeq_id::e_not_set) {
//                 vbar = vbar_prev;
//                 break;
//             }
//         }
//     } else {
//         return string::npos; // bad
//     }
//
//     return (vbar == string::npos) ? str.size() : vbar;
// }
// ```
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
// ```c++
// static bool
// s_SeqDB_ParseSeqIDs(const string              & line,
//                     vector< CRef< CSeq_id > > & seqids)
// {
//     // (Derived from s_ParseFastaDefline()).
//
//     seqids.clear();
//     size_t pos = 0;
//
//     while (pos < line.size()) {
//         size_t end = s_SeqDB_EndOfFastaID(line, pos);
//
//         if (end == string::npos) {
//             // We didn't get a clean parse -- ignore the data after
//             // this point, and return what we have.
//             break;
//         }
//
//         string element(line, pos, end - pos);
//
//         CRef<CSeq_id> id;
//
//         try {
//             id = new CSeq_id(element);
//         }
//         catch(invalid_argument &) {
//             // Maybe this should be done: "seqids.clear();"
//             break;
//         }
//
//         seqids.push_back(id);
//         pos = end + 1;
//     }
//
//     return ! seqids.empty();
// }
// ```
/// `s_SeqDB_ParseSeqIDs`: the Seq-ids at the start of `line` (`CSeq_id(element)` of every
/// element up to the first text that is not one). A `CSeqIdException` is not caught by NCBI's
/// `catch(invalid_argument &)`.
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
// ```c++
// ESeqDBIdType SeqDB_SimplifyAccession(const string & acc,
//                                      Int8         & num_id,
//                                      string       & str_id,
//                                      bool         & simpler)
// {
//     ESeqDBIdType result = eStringId;
//     num_id = (Uint4)-1;
//
//     vector< CRef< CSeq_id > > seqid_set;
//
//     if (s_SeqDB_ParseSeqIDs(acc, seqid_set)) {
//         // Something like SeqIdFindBest()
//         CRef<CSeq_id> bestid =
//             FindBestChoice(seqid_set, CSeq_id::BestRank);
//
//         result = SeqDB_SimplifySeqid(*bestid, & acc, num_id, str_id, simpler);
//     } else {
//
//         // Check for bare pdb accession with underscore (such as 12AS_A). These
//         // are not in the isam index and need to be translated to the
//         // standard form (pdb|12AS|A).
//         list< CRef<CSeq_id> > seqids;
//         try {
//             CSeq_id::ParseFastaIds(seqids, acc, false);
//         }
//         catch (...) {
//             seqids.clear();
//         }
//
//         if (!seqids.empty() && seqids.front()->IsPdb() &&
//             acc.find("_") != string::npos) {
//
//             str_id = seqids.front()->AsFastaString();
//             str_id = NStr::ToLower(str_id);
//         }
//         else if (!seqids.empty() && seqids.front()->IsLocal()) {
//             // Chec for gnl dbs
//             if( acc.find(":") != string::npos) {
//                 static const char* GNL_DBs[] = {"CDD", "SRA", "TSA", "GNOMON", NULL};
//                 string db_tag, gnl_id;
//                 NStr::SplitInTwo(acc, ":", db_tag, gnl_id);
//                 const char** p = GNL_DBs;
//                 for (;  p  &&  *p;  ++p) {
//                     if(NStr::EqualNocase(*p, db_tag.c_str())) {
//                         str_id  = "gnl|" + db_tag + "|"  + gnl_id;
//                         seqids.front().Reset();
//                         CRef<CSeq_id> new_id(new CSeq_id(str_id));
//                         seqids.front() = new_id;
//                         break;
//                     }
//                 }
//                 if(*p == NULL) {
//                     str_id =  acc;
//                 }
//             }
//             else {
//                 Int8 n_id = 0;
//                 if (NStr::StringToNumeric<Int8>(acc,&n_id, NStr::fConvErr_NoThrow)) {
//                     str_id = "lcl|" + acc;
//                 }
//                 else {
//                     str_id = acc;
//                 }
//             }
//         }
//         else {
//             str_id = acc;
//         }
//         result = eStringId;
//         simpler = false;
//     }
//
//     return result;
// }
// ```
/// `SeqDB_SimplifyAccession(acc, num_id, str_id, simpler)`: what the volume looks up for the
/// user's id `acc`.
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
        } else if string_to_int8(acc).is_some() {
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
// ```c++
// const string SeqDB_SimplifyAccession(const string &acc)
// {
//     Int8 num_id;
//     string str_id;
//     bool simpler(false);
//     ESeqDBIdType result = SeqDB_SimplifyAccession(acc, num_id, str_id, simpler);
//     if (result == eStringId) return str_id;
//     else return "";
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:608-614
// ```c++
// void CSeqDBGiList::PreprocessIdsForISAMSiLookup()
// {
//     NON_CONST_ITERATE(vector<CSeqDBGiList::SSiOid>, itr, m_SisOids) {
//         string str_id = SeqDB_SimplifyAccession(itr->si);
//         itr->si = NStr::ToLower(str_id);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1612-1618
// ```c++
// void CSeqDBNegativeList::PreprocessIdsForISAMSiLookup()
// {
//     NON_CONST_ITERATE(vector<string>, itr, m_Sis) {
//         string str_id = SeqDB_SimplifyAccession(*itr);
//         *itr = NStr::ToLower(str_id);
//     }
// }
// ```
/// `SeqDB_SimplifyAccession(acc)` followed by `NStr::ToLower`
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
