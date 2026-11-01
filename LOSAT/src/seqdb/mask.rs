//! The mask data of a BLAST database (`makeblastdb -mask_data`): the filtering algorithms
//! that were applied to the sequences, and the masked ranges of each sequence.
//!
//! NCBI reads them in `CSeqDBImpl` (`seqdbimpl.cpp`: `x_BuildMaskAlgorithmList`,
//! `GetMaskAlgorithmDetails`, `GetAvailableMaskAlgorithmDescriptions`, `GetMaskData`) with
//! the id table `CSeqDB_IdRemapper`. The data is in the column titled `BlastDb/MaskData` of
//! each volume ([`super::column`]): its metadata maps the volume's algorithm ids to
//! descriptions; the blob of a sequence holds, per algorithm, a list of ranges.
//!
//! The reader of `losat blastdbcmd`; LOSAT's searches do not call it. A database is opened with
//! [`SeqDb::open_in_with_masks`] to read it.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use super::column::{Blob, MASK_DATA_TITLE};
use super::volume::Volume;
use super::{SeqDb, SeqDbError};

/// A filtering algorithm of the database (`GetMaskAlgorithmDetails`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskAlgorithm {
    /// The database-wide id (the `algo_id` of [`SeqDb::mask_ranges`]).
    pub id: i32,
    /// The program: a `Blast-filter-program` number, or a name.
    pub program: String,
    /// The name of the program.
    pub program_name: String,
    /// The options the program ran with; empty for the defaults.
    pub options: String,
}

/// What [`SeqDb::mask_ranges`] found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MaskRanges {
    /// The masked ranges `(start, end)` of the algorithm in the sequence (`end` is
    /// exclusive), in the order the database stores them; empty when there are none.
    Ranges(Vec<(u32, u32)>),
    /// The volume of the sequence has a mask column and the sequence has mask data, but the
    /// volume does not know the algorithm. NCBI answers from a stale cache here (it may
    /// read the data of another algorithm), so a caller should refuse the search.
    AlgorithmNotInVolume {
        /// The volume (its path without extension).
        volume: PathBuf,
    },
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
// ```c++
// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Int8 value = StringToInt8(str, flags, base);
//     if ( value < kMin_Int  ||  value > kMax_Int ) {
//         S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
//     }
//     return (int) value;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-873
// ```c++
// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD(flags);
//
//     // Current position in the string
//     SIZE_TYPE pos = 0;
//
//     // Skip allowed leading symbols
//     if (flags & fAllowLeadingSymbols) {
//         bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
//         s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
//     }
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     default:
//         if (flags & fMandatorySign) {
//             S2N_CONVERT_ERROR_INVAL(Int8);
//         }
//         break;
//     }
//     SIZE_TYPE pos0 = pos;
//     // Check radix base
//     if ( !s_CheckRadix(str, pos, base) ) {
//         S2N_CONVERT_ERROR_RADIX(Int8, "bad numeric base '" +
//                                 NStr::IntToString(base)+ "'");
//     }
//
//     // Begin conversion
//     Int8 n = 0;
//     Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
//     Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);
//
//     // Number of symbols between two commas. '-1' means -- no comma yet.
//     int       comma  = -1;
//     SIZE_TYPE numpos = pos;
//
//     while (char ch = str[pos]) {
//         int  delta;   // corresponding numeric value of 'ch'
//
//         // Check on possible commas
//         CHECK_COMMAS;
//         // Sanity check
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
//
//     // Last checks
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     // Skip allowed trailing symbols
//     if (flags & fAllowTrailingSymbols) {
//         bool spaces = ((flags & fAllowTrailingSymbols) == fAllowTrailingSpaces);
//         s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
//     }
//     // Assign sign before the end pointer check
//     n = sign ? -n : n;
//     CHECK_ENDPTR(Int8);
//
//     return n;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2138-2141
// ```c++
// static bool s_IsNumericId(const string &id) {
//     Int4 nid(-1);
//     return NStr::StringToNumeric(id, &nid, NStr::fConvErr_NoThrow, 10);
// }
// ```
/// `NStr::StringToInt(s)` (flags 0, base 10): an optional sign and at least one digit up
/// to the first NUL, nothing else, within the range of an `int`.
pub(super) fn string_to_int(s: &str) -> Option<i32> {
    let bytes = s.as_bytes();
    let bytes = &bytes[..bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len())];
    let (negative, digits) = match bytes.first() {
        Some(b'-') => (true, &bytes[1..]),
        Some(b'+') => (false, &bytes[1..]),
        _ => (false, bytes),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let mut value: i128 = 0;
    for &digit in digits {
        value = (value * 10 + i128::from(digit - b'0')).min(1 << 70);
    }
    let value = if negative { -value } else { value };
    i32::try_from(value).ok()
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-303
// ```c++
//     void Do(TContainer&     target,
//             TPosContainer&  token_pos,
//             const TString&  empty_str = TString())
//     {
//         auto target_initial_size = target.size();
//
//         // Special cases
//         if (m_Str.empty()) {
//             return;
//         } else if (m_Delim.empty()) {
//             target.push_back(m_Str);
//             token_pos.push_back(0);
//             return;
//         }
//         // Do target space reservation (if applicable)
//         TReserveTrait::Reserve(*this, target, token_pos);
//
//         // Tokenization
//
//         CTempStringList part_collector(m_Storage);
//         SIZE_TYPE prev_pos;
//         SIZE_TYPE delim_pos = NPOS;
//         m_Pos = 0;
//         do {
//             Advance(&part_collector, &prev_pos, &delim_pos);
//             target.push_back(empty_str); // reserve space for value added next line
//             part_collector.Join(&target.back());
//             part_collector.Clear();
//             token_pos.push_back(prev_pos);
//         } while ( !AtEnd() );
//
//         if ( (m_Flags & NStr::fSplit_Truncate_End) == 0 ) {
//             // account training delimiter
//             if (delim_pos != NPOS) {
//                 // add empty token after last delimiter
//                 target.push_back(empty_str);
//                 token_pos.push_back(delim_pos + 1);
//             }
//         }
//         else {
// ```
/// `NStr::Split(s, ":", items)` with flags 0: every `:` separates, empty tokens are kept
/// (so `"30:"` is `["30", ""]`), and the empty string has no tokens.
pub(super) fn ncbi_split(s: &str) -> Vec<&str> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.split(':').collect()
    }
}

/// The name of a `Blast-filter-program` value (`blastdb.asn:51-59`).
fn filter_program_name(value: i32) -> Option<&'static str> {
    Some(match value {
        0 => "not-set",
        10 => "dust",
        20 => "seg",
        30 => "windowmasker",
        40 => "repeat",
        100 => "other",
        255 => "max",
        _ => return None,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2143-2165
// ```c++
// static const string* s_CheckUniqueValues(const map<string, string> & m)
// {
//     typedef map<string, string> TStringMap;
//
//     set<string> seen;
//
//     ITERATE(TStringMap, iter, m) {
//         string v = iter->second;
//         vector<string> items;
//         NStr::Split(v, ":", items);
//
//         if (items.size() == 4) {
//             v = items[2];
//         }
//
//         if (s_Contains(seen, v)) {
//             return & iter->second;
//         }
//
//         seen.insert(v);
//     }
//
//     return NULL;
// ```
/// `s_CheckUniqueValues`: the first value (in key order) whose description repeats an
/// earlier one.
fn first_duplicate<'a>(metadata: &[(&str, &'a str)]) -> Option<&'a str> {
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    for &(_, value) in metadata {
        let items = ncbi_split(value);
        let description = if items.len() == 4 { items[2] } else { value };
        if !seen.insert(description) {
            return Some(value);
        }
    }
    None
}

/// A metadata key as text. A key that is not UTF-8 cannot be an integer, so NCBI's
/// `StringToInt` fails on it as well: a damaged file.
fn utf8_key<'a>(bytes: &'a [u8], path: &Path) -> Result<&'a str, SeqDbError> {
    std::str::from_utf8(bytes).map_err(|_| SeqDbError::Format {
        path: path.to_path_buf(),
        message: "a metadata key of the mask column is not UTF-8 text".to_string(),
    })
}

/// A metadata value (the stored algorithm description) as text. NCBI keeps any bytes
/// (`-mask_desc` and `-mask_id` typed in Latin-1, say), so this is a limit of LOSAT, not a
/// damaged file.
fn utf8_value<'a>(bytes: &'a [u8], path: &Path) -> Result<&'a str, SeqDbError> {
    std::str::from_utf8(bytes).map_err(|_| SeqDbError::Unsupported {
        path: path.to_path_buf(),
        message: "a stored mask algorithm description is not UTF-8 text and LOSAT does not read it"
            .to_string(),
    })
}

/// Pads `text` with spaces on the right to `width` bytes (`setw(width) << left`; text that
/// is longer is not cut).
fn push_padded(out: &mut String, text: &str, width: usize) {
    out.push_str(text);
    for _ in text.len()..width {
        out.push(' ');
    }
}

/// The filtering algorithms of a database: `CSeqDB_IdRemapper` with the descriptions
/// (`GetAvailableMaskAlgorithms`, `GetMaskAlgorithmDetails`, `GetMaskAlgorithmId`).
///
/// The ids that volumes use for an algorithm may differ from volume to volume; the table
/// gives each description one database-wide id and remembers each volume's own.
#[derive(Debug)]
pub struct MaskAlgorithms {
    /// `m_NextId`: where the search for an unused id starts.
    next_id: i32,
    /// `m_IdToDesc`, with the file that first stated the description (for messages).
    id_to_desc: BTreeMap<i32, (String, PathBuf)>,
    /// `m_DescToId`, keyed by the description (item 2 of a four-item one).
    desc_to_id: BTreeMap<String, i32>,
    /// `m_RealIdToVolumeId`: by volume index, database-wide id to the volume's id.
    real_to_volume: BTreeMap<usize, BTreeMap<i32, i32>>,
}

impl MaskAlgorithms {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2167-2171
    // ```c++
    //
    // CSeqDB_IdRemapper::CSeqDB_IdRemapper()
    //     : m_NextId(100), m_Empty(true), m_CacheRealAlgo(-1)
    // {
    // }
    // ```
    fn new() -> Self {
        Self {
            next_id: 100,
            id_to_desc: BTreeMap::new(),
            desc_to_id: BTreeMap::new(),
            real_to_volume: BTreeMap::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2182-2232
    // ```c++
    // void CSeqDB_IdRemapper::AddMapping(int vol_id, int id, const string & desc)
    // {
    //     string real_desc = desc;
    //     vector<string> items;
    //     NStr::Split(desc, ":", items);
    //     if (items.size() == 4) {
    //         real_desc = items[2];
    //     }
    //     bool found_desc = s_Contains(m_DescToId, real_desc);
    //     bool found_id   = s_Contains(m_IdToDesc, id);
    //
    //     int real_id = id;
    //
    //     if (found_desc) {
    //         if ((! found_id) || (m_DescToId[real_desc] != id)) {
    //             // This description is mapped to a different ID.
    //             real_id = m_DescToId[real_desc];
    //         }
    //     } else {
    //         // New description.
    //
    //         if (found_id) {
    //             // Pick a 'synthetic' ID for this description,
    //             // i.e. one that is not actually used by any of
    //             // the existing volumes so far.
    //
    //             while(s_Contains(m_IdToDesc, m_NextId)) {
    //                 m_NextId++;
    //             }
    //
    //             real_id = m_NextId;
    //         }
    //
    //         // Add the new description.
    //
    //         m_IdToDesc[real_id] = desc;
    //         m_DescToId[real_desc] = real_id;
    //     }
    //
    //     m_RealIdToVolumeId[vol_id][real_id] = id;
    // }
    //
    // bool CSeqDB_IdRemapper::GetDesc(int algorithm_id, string & desc)
    // {
    //     if (! s_Contains(m_IdToDesc, algorithm_id)) {
    //        return false;
    //     }
    //
    //     desc = m_IdToDesc[algorithm_id];
    //     return true;
    // }
    // ```
    fn add_mapping(&mut self, volume: usize, id: i32, description: &str, source: &Path) {
        let items = ncbi_split(description);
        let real_desc = if items.len() == 4 {
            items[2]
        } else {
            description
        };
        let known_id = self.desc_to_id.get(real_desc).copied();
        let found_id = self.id_to_desc.contains_key(&id);

        let mut real_id = id;

        if let Some(known) = known_id {
            if !found_id || known != id {
                // This description is mapped to a different ID.
                real_id = known;
            }
        } else {
            // New description.
            if found_id {
                // Pick a 'synthetic' ID for this description, i.e. one that is not
                // actually used by any of the existing volumes so far.
                while self.id_to_desc.contains_key(&self.next_id) {
                    self.next_id += 1;
                }
                real_id = self.next_id;
            }
            self.id_to_desc
                .insert(real_id, (description.to_string(), source.to_path_buf()));
            self.desc_to_id.insert(real_desc.to_string(), real_id);
        }

        self.real_to_volume
            .entry(volume)
            .or_default()
            .insert(real_id, id);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2414-2467
    // ```c++
    // void CSeqDBImpl::x_BuildMaskAlgorithmList(CSeqDBLockHold & locked)
    // {
    //     m_Atlas.Lock(locked);
    //
    //     if (! m_AlgorithmIds.Empty()) {
    //         return;
    //     }
    //
    //     int col_id = x_GetMaskDataColumn(locked);
    //
    //     if (col_id < 0) {
    //         // No mask data column exists, therefore, the algorithms list
    //         // is empty, and we are done.
    //         return;
    //     }
    //
    //     CSeqDB_ColumnEntry & entry = *m_ColumnInfo[col_id];
    //
    //     typedef map<string,string> TStringMap;
    //
    //     // Results needed:
    //     // 1. Map global ids to desc.
    //     // 2. Map local id+vol -> global id
    //
    //     for(int vol_idx = 0; vol_idx < m_VolSet.GetNumVols(); vol_idx++) {
    //         // Get volume column #.
    //         int vol_col_id = entry.GetVolumeIndex(vol_idx);
    //
    //         if (vol_col_id < 0) {
    //             continue;
    //         }
    //
    //         CSeqDBVol * volp = m_VolSet.GetVolNonConst(vol_idx);
    //         const TStringMap & volmap =
    //             volp->GetColumnMetaData(vol_col_id, locked);
    //
    //         // Check for identical algorithm descriptions (should not happen.)
    //
    //         const string * dup = s_CheckUniqueValues(volmap);
    //
    //         if (dup != NULL) {
    //             ostringstream oss;
    //             oss << "Error: volume (" << volp->GetVolName()
    //                 << ") mask data has duplicates value (" << *dup <<  ")";
    //
    //             NCBI_THROW(CSeqDBException, eArgErr, oss.str());
    //         }
    //
    //         ITERATE(TStringMap, iter, volmap) {
    //             int id1 = NStr::StringToInt(iter->first);
    //             const string & desc1 = iter->second;
    //
    //             m_AlgorithmIds.AddMapping(vol_idx, id1, desc1);
    //         }
    // ```
    /// Builds the table from the mask columns of `volumes` (in database order). Every
    /// volume opens its columns first, as `x_GetColumnId` does, so an error in a later
    /// volume wins over one in an earlier volume's data.
    pub(super) fn build(volumes: &[Volume]) -> Result<Self, SeqDbError> {
        for volume in volumes {
            volume.columns()?;
        }
        refuse_other_first_columns(volumes)?;
        let mut table = Self::new();
        for (index, volume) in volumes.iter().enumerate() {
            let Some(column) = volume.mask_column()? else {
                continue;
            };
            let path = column.index_path();
            let mut metadata: Vec<(&str, &str)> = Vec::new();
            for (key, value) in column.metadata() {
                metadata.push((utf8_key(key, path)?, utf8_value(value, path)?));
            }

            // Check for identical algorithm descriptions (should not happen.)
            if let Some(duplicate) = first_duplicate(&metadata) {
                return Err(SeqDbError::Format {
                    path: path.to_path_buf(),
                    message: format!(
                        "Error: volume ({}) mask data has duplicates value ({})",
                        volume.base().display(),
                        duplicate
                    ),
                });
            }

            for (key, description) in metadata {
                let id = string_to_int(key).ok_or_else(|| SeqDbError::Format {
                    path: path.to_path_buf(),
                    message: format!("the mask algorithm id {key:?} is not an integer"),
                })?;
                table.add_mapping(index, id, description, path);
            }
        }
        Ok(table)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2173-2180
    // ```c++
    // void CSeqDB_IdRemapper::GetIdList(vector<int> & algorithms)
    // {
    //     typedef map<int,string> TIdMap;
    //
    //     ITERATE(TIdMap, iter, m_IdToDesc) {
    //         algorithms.push_back(iter->first);
    //     }
    // }
    // ```
    /// The database-wide algorithm ids, ascending.
    pub fn ids(&self) -> Vec<i32> {
        self.id_to_desc.keys().copied().collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2261-2268
    // ```c++
    // int CSeqDB_IdRemapper::GetAlgoId(const string & id)
    // {
    //     if (! s_Contains(m_DescToId, id)) {
    //         NCBI_THROW(CSeqDBException, eArgErr,
    //                    "Cannot find string algorithm id in algorithm map.");
    //     }
    //
    //     return m_DescToId[id];
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2292-2307
    // ```c++
    // int CSeqDBImpl::GetMaskAlgorithmId(const string & algo)
    // {
    //     if (m_UseGiMask) {
    //         return m_GiMask->GetAlgorithmId(algo);
    //     }
    //
    //     CHECK_MARKER();
    //     CSeqDBLockHold locked(m_Atlas);
    //     m_Atlas.Lock(locked);
    //
    //     if (m_AlgorithmIds.Empty()) {
    //         x_BuildMaskAlgorithmList(locked);
    //     }
    //
    //     return m_AlgorithmIds.GetAlgoId(algo);
    // }
    // ```
    /// `GetMaskAlgorithmId`: the id of the algorithm with this description (item 2 of a
    /// four-item description; otherwise the whole description, for example `"30:"`), or
    /// `None` where NCBI throws "Cannot find string algorithm id in algorithm map.".
    pub fn id_of(&self, name: &str) -> Option<i32> {
        self.desc_to_id.get(name).copied()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2344-2378
    // ```c++
    // static const string s_RestoreColon(const string &in) {
    //     const char l = 0x1;
    //     return NStr::Replace(in, string(l,1), ":");
    // }
    //
    // static
    // void s_GetDetails(const string          & desc,
    //                   string                & program,
    //                   string                & program_name,
    //                   string                & algo_opts)
    // {
    //     static const CEnumeratedTypeValues* enum_type_vals = NULL;
    //     if (enum_type_vals == NULL) {
    //         enum_type_vals = GetTypeInfo_enum_EBlast_filter_program();
    //     }
    //     _ASSERT(enum_type_vals);
    //
    //     vector<string> items;
    //     NStr::Split(desc, ":", items);
    //
    //     if (items.size() == 2) {
    //         EBlast_filter_program
    //            pid = (EBlast_filter_program) NStr::StringToInt(items[0]);
    //         program.assign(items[0]);
    //         program_name.assign(enum_type_vals->FindName(pid, false));
    //         algo_opts.assign(s_RestoreColon(items[1]));
    //     } else if (items.size() == 4) {
    //         program.assign(s_RestoreColon(items[0]));
    //         program_name.assign(s_RestoreColon(items[2]));
    //         algo_opts.assign(s_RestoreColon(items[1]));
    //     } else {
    //         NCBI_THROW(CSeqDBException, eArgErr,
    //                    "Error in stored mask algorithm description data.");
    //     }
    // }
    // ```
    fn parse_details(
        id: i32,
        description: &str,
        source: &Path,
    ) -> Result<MaskAlgorithm, SeqDbError> {
        let bad = |message: String| SeqDbError::Format {
            path: source.to_path_buf(),
            message,
        };
        // `s_RestoreColon`: 0x01 stands for a ':' inside an item.
        let restore = |item: &str| item.replace('\u{1}', ":");
        let items = ncbi_split(description);
        match items.len() {
            2 => {
                // `NStr::StringToInt` throws a CStringException on a text that is no
                // integer; `FindName(pid, false)` throws on a value the enum lacks.
                let number = string_to_int(items[0]).ok_or_else(|| {
                    bad(format!(
                        "the program {:?} of the mask algorithm description {description:?} \
                         is not an integer",
                        items[0]
                    ))
                })?;
                let name = filter_program_name(number).ok_or_else(|| {
                    bad(format!(
                        "the program {number} of the mask algorithm description \
                         {description:?} is not a Blast-filter-program value"
                    ))
                })?;
                Ok(MaskAlgorithm {
                    id,
                    program: items[0].to_string(),
                    program_name: name.to_string(),
                    options: restore(items[1]),
                })
            }
            4 => Ok(MaskAlgorithm {
                id,
                program: restore(items[0]),
                program_name: restore(items[2]),
                options: restore(items[1]),
            }),
            _ => Err(bad(
                "Error in stored mask algorithm description data.".to_string()
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2380-2412
    // ```c++
    // void CSeqDBImpl::GetMaskAlgorithmDetails(int                 algorithm_id,
    //                          string            & program,
    //                          string            & program_name,
    //                          string            & algo_opts)
    // {
    //     CHECK_MARKER();
    //     CSeqDBLockHold locked(m_Atlas);
    //     m_Atlas.Lock(locked);
    //
    //     string s;
    //     bool found;
    //
    //     if (m_UseGiMask) {
    //        // TODO:  to get description s
    //        s = m_GiMask->GetDesc(algorithm_id, locked);
    //        found = true;
    //     } else {
    //        if (m_AlgorithmIds.Empty()) {
    //            x_BuildMaskAlgorithmList(locked);
    //        }
    //        found = m_AlgorithmIds.GetDesc(algorithm_id, s);
    //     }
    //
    //     if (found == false) {
    //         CNcbiOstrstream oss;
    //         oss << "Filtering algorithm ID " << algorithm_id
    //             << " is not supported." << endl;
    //         oss << GetAvailableMaskAlgorithmDescriptions();
    //         NCBI_THROW(CSeqDBException, eArgErr, CNcbiOstrstreamToString(oss));
    //     }
    //
    //     s_GetDetails(s, program, program_name, algo_opts);
    // }
    // ```
    /// `GetMaskAlgorithmDetails`. An id the database does not have gives NCBI's text
    /// ("Filtering algorithm ID <id> is not supported." and the list of algorithms) in
    /// [`SeqDbError::Range`].
    pub fn details(&self, id: i32) -> Result<MaskAlgorithm, SeqDbError> {
        let Some((description, source)) = self.id_to_desc.get(&id) else {
            return Err(SeqDbError::Range(format!(
                "Filtering algorithm ID {id} is not supported.\n{}",
                self.descriptions()?
            )));
        };
        Self::parse_details(id, description, source)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2309-2342
    // ```c++
    // string CSeqDBImpl::GetAvailableMaskAlgorithmDescriptions()
    // {
    //     vector<int> algorithms;
    //     GetAvailableMaskAlgorithms(algorithms);
    //     if (algorithms.empty()) {
    //         return kEmptyStr;
    //     }
    //
    //     CNcbiOstrstream retval;
    //     retval << endl
    //            << "Available filtering algorithms applied to database sequences:"
    //            << endl << endl;
    //
    //     retval << setw(13) << left << "Algorithm ID"
    //            << setw(40) << left << "Algorithm name"
    //            << setw(40) << left << "Algorithm options" << endl;
    //     ITERATE(vector<int>, algo_id, algorithms) {
    //         string algo, algo_opts, algo_name;
    //         GetMaskAlgorithmDetails(*algo_id, algo, algo_name, algo_opts);
    //         if (algo_opts.empty()) {
    //             algo_opts.assign("default options used");
    //         }
    //         if (s_IsNumericId(algo)) {
    //             retval << setw(13) << left << (*algo_id)
    //                    << setw(40) << left << algo_name
    //                    << setw(40) << left << algo_opts << endl;
    //         } else {
    //             retval << setw(13) << left << (*algo_id)
    //                    << setw(40) << left << algo
    //                    << setw(40) << left << algo_opts << endl;
    //         }
    //     }
    //     return CNcbiOstrstreamToString(retval);
    // }
    // ```
    /// `GetAvailableMaskAlgorithmDescriptions`: the table `blastdbcmd -list`-style tools
    /// print; empty when the database has no mask data.
    pub fn descriptions(&self) -> Result<String, SeqDbError> {
        let ids = self.ids();
        if ids.is_empty() {
            return Ok(String::new());
        }
        let mut out = String::new();
        out.push_str("\nAvailable filtering algorithms applied to database sequences:\n\n");
        push_padded(&mut out, "Algorithm ID", 13);
        push_padded(&mut out, "Algorithm name", 40);
        push_padded(&mut out, "Algorithm options", 40);
        out.push('\n');
        for id in ids {
            let details = self.details(id)?;
            let options = if details.options.is_empty() {
                "default options used"
            } else {
                details.options.as_str()
            };
            let name = if string_to_int(&details.program).is_some() {
                &details.program_name
            } else {
                &details.program
            };
            push_padded(&mut out, &id.to_string(), 13);
            push_padded(&mut out, name, 40);
            push_padded(&mut out, options, 40);
            out.push('\n');
        }
        Ok(out)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2244-2259
    // ```c++
    // int CSeqDB_IdRemapper::RealToVol(int vol_idx, int algo_id)
    // {
    //     if (! s_Contains(m_RealIdToVolumeId, vol_idx)) {
    //         NCBI_THROW(CSeqDBException, eArgErr,
    //                    "Cannot find volume in algorithm map.");
    //     }
    //
    //     map<int,int> & trans = m_RealIdToVolumeId[vol_idx];
    //
    //     if (! s_Contains(trans, algo_id)) {
    //         NCBI_THROW(CSeqDBException, eArgErr,
    //                    "Cannot find volume algorithm in algorithm map.");
    //     }
    //
    //     return trans[algo_id];
    // }
    // ```
    /// `RealToVol`: the id the volume `volume` (index in the database) uses for the
    /// database-wide `algo_id`; `None` where NCBI throws ("Cannot find volume in algorithm
    /// map." or "Cannot find volume algorithm in algorithm map.").
    fn volume_algorithm(&self, volume: usize, algo_id: i32) -> Option<i32> {
        self.real_to_volume.get(&volume)?.get(&algo_id).copied()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2470-2506
// ```c++
//     m_AlgorithmIds.SetNotEmpty();
// }
//
// struct SReadInt4 {
//     enum { numeric_size = 4 };
//
//     static int Read(CBlastDbBlob & blob)
//     {
//         return blob.ReadInt4();
//     }
//
//     static void Read(CBlastDbBlob & blob, int n,
//                      CSeqDB::TSequenceRanges & ranges)
//     {
//         const void * src = (const void *) blob.ReadRaw(n*8);
//         ranges.append(src, n);
//     }
// };
//
// template<class TRead>
// void s_ReadRanges(int                       vol_algo,
//                   CSeqDB::TSequenceRanges & ranges,
//                   CBlastDbBlob            & blob)
// {
//     int num_ranges = TRead::Read(blob);
//
//     for(int rng = 0; rng < num_ranges; rng++) {
//         int algo = TRead::Read(blob);
//         int num_pairs = TRead::Read(blob);
//         if (algo == vol_algo) {
//             TRead::Read(blob, num_pairs, ranges);
//             break;
//         }
//         int skip_amt = num_pairs * 2 * TRead::numeric_size;
//         blob.SeekRead(blob.GetReadOffset() + skip_amt);
//     }
// }
// ```
/// `s_ReadRanges<SReadInt4>`: the ranges of `vol_algo` in a mask blob. The pairs are
/// little-endian `u32` (the `.?xc` file). Every read is checked against the blob.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3276-3292
// ```c++
// void CSeqDBVol::GetColumnBlob(int              col_id,
//                               int              oid,
//                               CBlastDbBlob   & blob,
//                               bool             keep,
//                               CSeqDBLockHold & locked)
// {
//     if (! m_HaveColumns) {
//         x_OpenAllColumns(locked);
//     }
//
//     if ((int)m_Columns.size() == 0 || m_Columns[col_id].Empty())
// 	return;
//
//     _ASSERT(col_id >= 0);
//
//     m_Columns[col_id]->GetBlob(oid, blob, keep, & locked);
// }
// ```
// `col_id` is the database's number of the mask column (0), not the volume's: NCBI reads
// the first column of a volume as its masks. `makeblastdb` writes the mask column alone.
/// Refuses a database with masks that has a volume whose first column is another one.
fn refuse_other_first_columns(volumes: &[Volume]) -> Result<(), SeqDbError> {
    let mut masked = false;
    for volume in volumes {
        masked |= volume.mask_column()?.is_some();
    }
    if !masked {
        return Ok(());
    }
    for volume in volumes {
        if let Some(first) = volume.columns()?.first() {
            if first.title() != MASK_DATA_TITLE {
                return Err(SeqDbError::Unsupported {
                    path: first.index_path().to_path_buf(),
                    message: "the first column of the volume is not its mask data, which NCBI \
                              reads as the masks"
                        .to_string(),
                });
            }
        }
    }
    Ok(())
}

fn read_ranges(blob: &[u8], vol_algo: i32, path: &Path) -> Result<Vec<(u32, u32)>, SeqDbError> {
    let format = |message: String| SeqDbError::Format {
        path: path.to_path_buf(),
        message,
    };
    let mut blob = Blob::new(blob, 0);
    let mut ranges = Vec::new();
    let num_ranges = blob.read_i32().map_err(format)?;
    // A skipped entry with a negative `num_pairs` moves the read offset backwards, and NCBI
    // then reads the same entries again for up to 2^31 passes. A pass that starts at an
    // offset where an earlier pass started cannot match or fail differently (that pass
    // would have returned), so the loop ends there with no ranges, as NCBI's does.
    let mut pass_starts = std::collections::HashSet::new();
    for _ in 0..num_ranges.max(0) {
        if !pass_starts.insert(blob.read_offset()) {
            break;
        }
        let algo = blob.read_i32().map_err(format)?;
        let num_pairs = blob.read_i32().map_err(format)?;
        if algo == vol_algo {
            // `ReadRaw(n*8)`: a negative count is NCBI's "hit end of data".
            let bytes = blob.read_raw(i64::from(num_pairs) * 8).map_err(format)?;
            for pair in bytes.chunks_exact(8) {
                let start = u32::from_le_bytes([pair[0], pair[1], pair[2], pair[3]]);
                let end = u32::from_le_bytes([pair[4], pair[5], pair[6], pair[7]]);
                ranges.push((start, end));
            }
            break;
        }
        let skip = i64::from(num_pairs) * 8;
        blob.seek_read(blob.read_offset() + skip);
    }
    Ok(ranges)
}

impl SeqDb {
    /// The filtering algorithms of the database, built from the mask columns of its
    /// volumes the first time they are asked for (`x_BuildMaskAlgorithmList`). Empty for a
    /// database without mask data. A database opened with [`SeqDb::open`] has no columns;
    /// open it with [`SeqDb::open_in_with_masks`].
    pub fn mask_algorithms(&self) -> Result<&MaskAlgorithms, SeqDbError> {
        if let Some(table) = self.mask_algorithms.get() {
            return Ok(table);
        }
        let table = MaskAlgorithms::build(&self.volumes)?;
        Ok(self.mask_algorithms.get_or_init(|| table))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2508-2547
    // ```c++
    // void CSeqDBImpl::GetMaskData(int                       oid,
    //                              int                       algo_id,
    //                              CSeqDB::TSequenceRanges & ranges)
    // {
    //     CHECK_MARKER();
    //
    //     // This reads the data written by CWriteDB_Impl::SetMaskData
    //     ranges.clear();
    //
    //     CSeqDBLockHold locked(m_Atlas);
    //     m_Atlas.Lock(locked);
    //
    //     if (m_UseGiMask) {
    //         m_GiMask->GetMaskData(algo_id, x_GetSeqGI(oid, locked), ranges, locked);
    //         return;
    //     }
    //
    //     if (m_AlgorithmIds.Empty()) {
    //         x_BuildMaskAlgorithmList(locked);
    //     }
    //
    //     int vol_oid = 0, vol_idx = -1;
    //
    //     CSeqDBVol * vol = const_cast<CSeqDBVol*>
    //         (m_VolSet.FindVol(oid, vol_oid, vol_idx));
    //
    //     if (! vol) {
    //         NCBI_THROW(CSeqDBException, eArgErr, CSeqDB::kOidNotFound);
    //     }
    //
    //     // Get the data.
    //
    //     CBlastDbBlob blob;
    //     vol->GetColumnBlob(x_GetMaskDataColumn(locked), vol_oid, blob, false, locked);
    //
    //     if (blob.Size() != 0) {
    //         // If there actually is mask data, then we need to do the
    //         // algorithm translation.
    //
    //         int vol_algo_id = -1;
    // ```
    /// `GetMaskData`: the ranges that the algorithm `algo_id` (a database-wide id, see
    /// [`MaskAlgorithms::ids`]) masked in the sequence `ordinal`.
    ///
    /// No ranges are returned (as in NCBI) when the sequence's volume has no mask column or
    /// the sequence has no mask data. A sequence with mask data in a volume that does not
    /// know the algorithm is [`MaskRanges::AlgorithmNotInVolume`]: NCBI's volume-algorithm
    /// cache (`GetVolAlgo`) keeps a stale value there, which LOSAT does not reproduce.
    pub fn mask_ranges(&self, ordinal: u32, algo_id: i32) -> Result<MaskRanges, SeqDbError> {
        let algorithms = self.mask_algorithms()?;
        if ordinal >= self.num_seqs() {
            return Err(SeqDbError::Range("OID not found".to_string()));
        }
        let (volume_index, volume_oid) = self.find_volume(ordinal)?;
        let volume = &self.volumes[volume_index];
        let Some(column) = volume.mask_column()? else {
            return Ok(MaskRanges::Ranges(Vec::new()));
        };
        let blob = column.blob(volume_oid)?;
        if blob.is_empty() {
            return Ok(MaskRanges::Ranges(Vec::new()));
        }
        // If there actually is mask data, then we need to do the algorithm translation.
        let Some(volume_algo) = algorithms.volume_algorithm(volume_index, algo_id) else {
            return Ok(MaskRanges::AlgorithmNotInVolume {
                volume: volume.base().to_path_buf(),
            });
        };
        read_ranges(&blob, volume_algo, column.data_path()).map(MaskRanges::Ranges)
    }
}
