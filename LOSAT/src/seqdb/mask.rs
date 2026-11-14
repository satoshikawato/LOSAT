//! The mask data of a BLAST database (`makeblastdb -mask_data`): the filtering algorithms
//! that were applied to the sequences, and the masked ranges of each sequence.
//!
//! NCBI reads them in `CSeqDBImpl` (`seqdbimpl.cpp`: `x_BuildMaskAlgorithmList`,
//! `GetMaskAlgorithmDetails`, `GetAvailableMaskAlgorithmDescriptions`, `GetMaskData`) with
//! the id table `CSeqDB_IdRemapper`. The data is in the column titled `BlastDb/MaskData` of
//! each volume ([`super::column`]): its metadata maps the volume's algorithm ids to
//! descriptions; the blob of a sequence holds, per algorithm, a list of ranges.
//!
//! The reader of `losat blastdbcmd`, and of BLASTN's `-db_soft_mask` (`blastinput::db_mask`):
//! NCBI's searches read the columns only with a masking algorithm named by `-db_soft_mask`,
//! `-db_hard_mask` or an imported search strategy (LOSAT refuses the last two), and so do
//! LOSAT's. The columns are opened on first use.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::ncbi_str::string_to_int;

use super::column::{Blob, MASK_DATA_TITLE};
use super::volume::Volume;
use super::{BlockRange, SeqDb, SeqDbError, SeqLocation};

/// A filtering algorithm of the database (`GetMaskAlgorithmDetails`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskAlgorithm {
    /// The database-wide id (the `algo_id` of [`SeqDb::mask_ranges_at`]).
    pub id: i32,
    /// The program: a `Blast-filter-program` number, or a name.
    pub program: String,
    /// The name of the program.
    pub program_name: String,
    /// The options the program ran with; empty for the defaults.
    pub options: String,
}

/// What [`SeqDb::mask_ranges_at`] found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MaskRanges {
    /// The masked ranges `(start, end)` of the algorithm in the sequence (`end` is
    /// exclusive), in the order the database stores them; empty when there are none. A mask
    /// column read from its big-endian data file `.?xb` alone gives each number byte-reversed,
    /// as NCBI does, so the ranges can be unsorted, inverted (`start > end`) or past the end
    /// of the sequence. A consumer that indexes letters with them must check them: blastdbcmd
    /// `-outfmt %s` refuses them, its FASTA path clips them, and a port of the searches'
    /// `-db_soft_mask` / `-db_hard_mask`, where NCBI's use of such ranges is undefined, must
    /// refuse them (line NP, M/DECISIONS 9).
    Ranges(Vec<(u32, u32)>),
    /// The volume of the sequence has a mask column and the sequence has mask data, but the
    /// volume does not know the algorithm. NCBI answers from a stale cache here (it may
    /// read the data of another algorithm), so a caller should refuse the search.
    AlgorithmNotInVolume {
        /// The volume (its path without extension).
        volume: PathBuf,
    },
}

/// Why [`SeqDb::name_mask_details`] has no algorithm of a name of the `-db` list.
#[derive(Debug)]
pub enum NameMaskError<'a> {
    /// The error of reading the table of the name's database (`NameRecord::mask_algorithms`).
    Table(&'a SeqDbError),
    /// `GetMaskAlgorithmDetails` of the name's table ([`MaskAlgorithms::details`]): NCBI's error
    /// of an id that it lacks, or the error of the algorithm's stored description.
    Details(SeqDbError),
    /// The records were read without their tables (`read_name_records` without `masks`): a
    /// LOSAT caller's error.
    NotRead,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2138-2141
// ```c++
// static bool s_IsNumericId(const string &id) {
//     Int4 nid(-1);
//     return NStr::StringToNumeric(id, &nid, NStr::fConvErr_NoThrow, 10);
// }
// ```
// `StringToNumeric(str, &int, flags, base)` is `NStr::StringToInt` that returns whether it
// converted the text (`*value || !errno`; a converted 0 leaves `errno` at 0).
/// `s_IsNumericId(id)` (`blastinput::db_mask::FilteringAlgorithm::parse` is the same test of
/// the searches, `uniform_search.cpp:149-153`).
fn is_numeric_id(id: &str) -> bool {
    string_to_int(id.as_bytes(), 0).is_ok()
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

// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-59
// ```text
// Blast-filter-program ::= INTEGER {
//     not-set             (0),
//     dust                (10),
//     seg                 (20),
//     windowmasker        (30),
//     repeat              (40),
//     other               (100),
//     max                 (255)
// }
// ```
/// The named values of `Blast-filter-program`, with their names as NCBI prints them: the one
/// table of LOSAT (the mask writers of the masking tools, makeblastdb's mask reader and
/// registry, and this reader of the column), made from [`FilterProgram`], so a name or value
/// belongs to its variant whatever the order of the variants.
pub(crate) const FILTER_PROGRAMS: [(&str, i32); FilterProgram::ALL.len()] = {
    let mut table = [("", 0); FilterProgram::ALL.len()];
    let mut i = 0;
    while i < table.len() {
        let program = FilterProgram::ALL[i];
        table[i] = (program.name(), program.value());
        i += 1;
    }
    table
};

/// A named value of `Blast-filter-program`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FilterProgram {
    NotSet,
    Dust,
    Seg,
    Windowmasker,
    Repeat,
    Other,
    Max,
}

// Every variant is in `FilterProgram::ALL`, once, in the order of the enum: the discriminant of
// the last one gives the count, and each position holds the variant of that discriminant.
const _: () = {
    assert!(FilterProgram::ALL.len() == FilterProgram::Max as usize + 1);
    let mut i = 0;
    while i < FilterProgram::ALL.len() {
        assert!(FilterProgram::ALL[i] as usize == i);
        i += 1;
    }
};

impl FilterProgram {
    /// Every variant, in blastdb.asn's order (the rows of [`FILTER_PROGRAMS`]).
    pub(crate) const ALL: [FilterProgram; 7] = [
        FilterProgram::NotSet,
        FilterProgram::Dust,
        FilterProgram::Seg,
        FilterProgram::Windowmasker,
        FilterProgram::Repeat,
        FilterProgram::Other,
        FilterProgram::Max,
    ];

    pub(crate) const fn name(self) -> &'static str {
        match self {
            FilterProgram::NotSet => "not-set",
            FilterProgram::Dust => "dust",
            FilterProgram::Seg => "seg",
            FilterProgram::Windowmasker => "windowmasker",
            FilterProgram::Repeat => "repeat",
            FilterProgram::Other => "other",
            FilterProgram::Max => "max",
        }
    }

    pub(crate) const fn value(self) -> i32 {
        match self {
            FilterProgram::NotSet => 0,
            FilterProgram::Dust => 10,
            FilterProgram::Seg => 20,
            FilterProgram::Windowmasker => 30,
            FilterProgram::Repeat => 40,
            FilterProgram::Other => 100,
            FilterProgram::Max => 255,
        }
    }
}

/// The name of a `Blast-filter-program` value.
fn filter_program_name(value: i32) -> Option<&'static str> {
    FILTER_PROGRAMS
        .iter()
        .find(|&&(_, number)| number == value)
        .map(|&(name, _)| name)
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
#[derive(Debug, Clone)]
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
                let id = string_to_int(key.as_bytes(), 0).map_err(|_| SeqDbError::Format {
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
                let number = string_to_int(items[0].as_bytes(), 0).map_err(|_| {
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
    /// ("Filtering algorithm ID `<id>` is not supported." and the list of algorithms) in
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
    /// `GetAvailableMaskAlgorithmDescriptions`: the table of `blastdbcmd -info`; empty when
    /// the database has no mask data.
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
            let name = if is_numeric_id(&details.program) {
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

/// `s_ReadRanges<SReadInt4>`: the ranges of `vol_algo` in a mask blob. The pairs are copied as
/// they are into `TSeqPos` (`ranges.append`, a `memcpy`), so they read in the host's order:
/// little-endian on every host of NCBI BLAST+ 2.17.0 (Linux x64 and arm64, macOS x64 and
/// arm64, Windows x64) and of LOSAT (x86-64, aarch64, wasm32). There a blob of the `.?xc` file
/// (`WriteInt4_LE`) gives the stored ranges, and one of a `.?xb` file read alone
/// (`column::open_all`; `WriteInt4`) gives each number with its bytes reversed, as NCBI reads
/// it. Every read is checked against the blob. The ranges go into `ranges` (appended to), with
/// `pass_starts` the set of the offsets where the passes started (cleared first; a caller that
/// reads many blobs keeps one set).
fn read_ranges_into(
    blob: &[u8],
    vol_algo: i32,
    path: &Path,
    ranges: &mut Vec<(u32, u32)>,
    pass_starts: &mut std::collections::HashSet<i64>,
) -> Result<(), SeqDbError> {
    let format = |message: String| SeqDbError::Format {
        path: path.to_path_buf(),
        message,
    };
    let mut blob = Blob::new(blob, 0);
    let num_ranges = blob.read_i32().map_err(format)?;
    // A skipped entry with a negative `num_pairs` moves the read offset backwards, and NCBI
    // then reads the same entries again for up to 2^31 passes. A pass that starts at an
    // offset where an earlier pass started cannot match or fail differently (that pass
    // would have returned), so the loop ends there with no ranges, as NCBI's does.
    pass_starts.clear();
    for _ in 0..num_ranges.max(0) {
        if !pass_starts.insert(blob.read_offset()) {
            break;
        }
        let algo = blob.read_i32().map_err(format)?;
        let num_pairs = blob.read_i32().map_err(format)?;
        if algo == vol_algo {
            // `ReadRaw(n*8)`: a negative count is NCBI's "hit end of data".
            let bytes = blob.read_raw(i64::from(num_pairs) * 8).map_err(format)?;
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1418-1423
            // ```c++
            //             ITERATE(TPairVector, r2, r1->offsets) {
            //                 blob.WriteInt4(r2->first);
            //                 blob.WriteInt4(r2->second);
            //                 blob2.WriteInt4_LE(r2->first);
            //                 blob2.WriteInt4_LE(r2->second);
            //             }
            // ```
            // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:302-306
            // ```c++
            //         void append(const void *src, size_type num_elements) {
            //             reserve(_size + num_elements);
            //             memcpy(&_data[1+ 2*_size], src, num_elements * sizeof(value_type));
            //             _size += num_elements;
            //         }
            // ```
            // `blob` goes to `.?xb`, `blob2` to `.?xc`. The copy keeps the bytes in host order;
            // on a little-endian host they are the numbers of `blob2`, and those of `blob`
            // byte-reversed (`blastdbcmd -outfmt %m` of a column without its `.?xc`, line BDC's
            // probe M784: 120-161 prints as 2013265920-2701131776).
            for pair in bytes.chunks_exact(8) {
                let start = u32::from_ne_bytes([pair[0], pair[1], pair[2], pair[3]]);
                let end = u32::from_ne_bytes([pair[4], pair[5], pair[6], pair[7]]);
                ranges.push((start, end));
            }
            break;
        }
        let skip = i64::from(num_pairs) * 8;
        blob.seek_read(blob.read_offset() + skip);
    }
    Ok(())
}

impl SeqDb {
    /// The filtering algorithms of the database, built from the mask columns of its
    /// volumes the first time they are asked for (`x_BuildMaskAlgorithmList`). Empty for a
    /// database without mask data.
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
    /// [`MaskAlgorithms::ids`]) masked in the sequence at `at` (`SeqDb::locate_oid`: NCBI's
    /// `FindVol`, whose `OID not found` NCBI throws here after it reads the algorithms, which
    /// a caller has read before it finds an OID).
    ///
    /// No ranges are returned (as in NCBI) when the sequence's volume has no mask column or
    /// the sequence has no mask data. A sequence with mask data in a volume that does not
    /// know the algorithm is [`MaskRanges::AlgorithmNotInVolume`]: NCBI's volume-algorithm
    /// cache (`GetVolAlgo`) keeps a stale value there, which LOSAT does not reproduce.
    ///
    /// The ranges of a mask column read from its `.?xb` alone are byte-reversed, as NCBI
    /// reads them (line SW, SW8 DECISIONS 40): unsorted, inverted or past the sequence, so
    /// a caller that indexes letters must check them ([`MaskRanges::Ranges`]).
    pub fn mask_ranges_at(&self, at: SeqLocation, algo_id: i32) -> Result<MaskRanges, SeqDbError> {
        let SeqLocation {
            volume, vol_oid, ..
        } = at;
        let mut block = MaskBlock::default();
        self.read_mask_block(
            BlockRange {
                volume,
                first: vol_oid,
                end: vol_oid + 1,
                first_ordinal: at.id,
            },
            algo_id,
            &mut block,
        );
        match block.ranges(vol_oid) {
            Ok(BlockMasks::Ranges(ranges)) => return Ok(MaskRanges::Ranges(ranges.to_vec())),
            Ok(BlockMasks::AlgorithmNotInVolume) => {
                return Ok(MaskRanges::AlgorithmNotInVolume {
                    volume: self.volumes[volume].base().to_path_buf(),
                })
            }
            Err(BlockMaskError::Ranges(error)) => return Err(error),
            Err(BlockMaskError::Read(_)) => {}
        }
        Err(block.into_read_error(vol_oid))
    }

    // No NCBI counterpart (reads only): NCBI maps the column files and its `GetMaskData` points into
    // them for each sequence (seqdbimpl.cpp:2508-2547, quoted at `mask_ranges_at`).
    /// The masks of the sequences of `range` for the filtering algorithm `algo_id` into `block`
    /// (unless it holds them already), as a search reads them with the sequences of its blocks
    /// (`SeqDb::read_block`): one read of the mask column's offset array for the block and one of
    /// its blobs (`Column::read_blobs`); [`MaskBlock::ranges`] gives each sequence's, with
    /// [`SeqDb::mask_ranges_at`]'s checks and errors in its order (the algorithms' table and the
    /// volume's columns: the error of every sequence of the block).
    pub fn read_mask_block(&self, range: BlockRange, algo_id: i32, block: &mut MaskBlock) {
        if block.read == Some((range, algo_id)) {
            return;
        }
        block.read = Some((range, algo_id));
        block.error = None;
        block.column = false;
        block.volume_algo = None;
        let algorithms = match self.mask_algorithms() {
            Ok(algorithms) => algorithms,
            Err(error) => {
                block.error = Some(error);
                return;
            }
        };
        let column = match self.volumes[range.volume].mask_column() {
            Ok(Some(column)) => column,
            Ok(None) => return,
            Err(error) => {
                block.error = Some(error);
                return;
            }
        };
        block.column = true;
        // If there actually is mask data, then we need to do the algorithm translation (the
        // volume's id, read for the blobs that are not empty).
        block.volume_algo = algorithms.volume_algorithm(range.volume, algo_id);
        if block.path.as_path() != column.data_path() {
            block.path = column.data_path().to_path_buf();
        }
        column.read_blobs(range.first, range.end, &mut block.blobs);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:479-492
    // ```c++
    //     vector<CTempString> dbs;
    //     SeqDB_SplitQuoted(blastdb_names, dbs, true);
    //     retval.reserve(dbs.size());
    //
    //     ITERATE(vector<CTempString>, i, dbs) {
    //         CAlignFormatUtil::SDbInfo info;
    //         info.is_protein = is_protein;
    //         bool success = false;
    // 	// Unsafe OK as kDbName only used in this loop.
    //         const string kDbName = NStr::TruncateSpaces_Unsafe(*i);
    //         if (kDbName.empty())
    //             continue;
    //
    //         success = s_FillDbInfoLocally(kDbName, info, dbfilt_algorithm);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:376-406
    // ```c++
    //     CRef<CSeqDB> seqdb(new CSeqDB(dbname, info.is_protein
    //                           ? CSeqDB::eProtein : CSeqDB::eNucleotide));
    //     ...
    //     if (dbfilt_algorithm == -1) {
    //         return true;
    //     }
    //     ...
    //     string filtering_algorithm;
    //     seqdb->GetMaskAlgorithmDetails(dbfilt_algorithm,
    //                                    filtering_algorithm,
    //                                    info.filt_algorithm_name,
    //                                    info.filt_algorithm_options);
    // ```
    /// What the formatter of a search with the filtering algorithm `algo_id` (not -1) reads of
    /// each name of the `-db` list, in the order of the list: that name's database's own
    /// `GetMaskAlgorithmDetails`, from the table its record keeps (`NameRecord::mask_algorithms`,
    /// read by [`SeqDb::read_name_records`] with `masks`: the names are not opened again). The
    /// first name whose table could not be read gives that error; the first whose database does
    /// not have the id gives NCBI's (`MaskAlgorithms::details`: "Filtering algorithm ID ... is
    /// not supported." and that database's table). One algorithm per record, so the reports'
    /// lines of the names (`name_records`) and these pair one to one.
    pub fn name_mask_details(&self, algo_id: i32) -> Result<Vec<MaskAlgorithm>, NameMaskError<'_>> {
        self.name_records
            .iter()
            .map(|record| match &record.mask_algorithms {
                Some(Ok(table)) => table.details(algo_id).map_err(NameMaskError::Details),
                Some(Err(error)) => Err(NameMaskError::Table(error)),
                None => Err(NameMaskError::NotRead),
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2234-2242
    // ```c++
    // int CSeqDB_IdRemapper::GetVolAlgo(int vol_idx, int algo_id)
    // {
    //     if (vol_idx != m_CacheVolIndex || algo_id != m_CacheRealAlgo) {
    //         m_CacheVolIndex = vol_idx;
    //         m_CacheRealAlgo = algo_id;
    //         m_CacheVolAlgo = RealToVol(vol_idx, algo_id);
    //     }
    //     return m_CacheVolAlgo;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3345-3358
    // (quoted at `column::open_all`): a little-endian host reads a column from its `.?xc` file,
    // or from its `.?xb` file when it has no `.?xc`.
    /// LOSAT's refusals of a search masked with the filtering algorithm `algo_id` (line NP,
    /// M/DECISIONS 9; decision D5 of M4a-2), before the search writes anything: a volume with a
    /// mask column that does not know the algorithm (NCBI translates the id of each of its
    /// sequences with mask data through the volume-algorithm cache above, which keeps the value of
    /// another volume, or an unset one on its first call), and a volume whose mask column is read
    /// from its big-endian `.?xb` file alone (NCBI's search then reads each number byte-reversed:
    /// ranges past the end of the sequence or inverted, which its scan uses as they are). A volume
    /// without a mask column has no masks, as in NCBI.
    pub fn check_search_masks(&self, algo_id: i32) -> Result<(), SeqDbError> {
        let algorithms = self.mask_algorithms()?;
        for (index, volume) in self.volumes.iter().enumerate() {
            let Some(column) = volume.mask_column()? else {
                continue;
            };
            let refusal = if algorithms.volume_algorithm(index, algo_id).is_none() {
                SeqDbError::Unsupported {
                    path: volume.base().to_path_buf(),
                    message: format!(
                        "the volume has mask data without the filtering algorithm {algo_id}, \
                         for which NCBI BLAST+ 2.17.0 reads the masks of another volume's \
                         algorithm (or an unset one)"
                    ),
                }
            } else if column.reads_byte_reversed() {
                SeqDbError::Unsupported {
                    path: column.data_path().to_path_buf(),
                    message: "the mask data of the volume has its big-endian file alone, whose \
                              ranges NCBI BLAST+ 2.17.0 reads byte-reversed on this host"
                        .to_string(),
                }
            } else {
                continue;
            };
            if self.searches_a_blob_of(index, column) {
                return Err(refusal);
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321 (quoted at
    // `MaskBlock::ranges`): the search reads the masks of the sequences it searches only.
    /// Whether the search reads a blob of the mask column `column` of volume `volume`: one of its
    /// searched sequences (those the list and the alias files' ranges include) has one, by the
    /// offset array, before a sequence whose entries cannot be read (the search fails there
    /// first; `Column::first_blob_in`).
    fn searches_a_blob_of(&self, volume: usize, column: &super::column::Column) -> bool {
        let start = self.volume_starts[volume];
        let mut ordinal = self.ordinal_of(start);
        let end_ordinal = self.ordinal_of(self.volume_starts[volume + 1]);
        while let Some(run) = self.x_block_range(ordinal, end_ordinal) {
            if run.volume != volume {
                break;
            }
            match column.first_blob_in(run.first, run.end) {
                Some(found) => return found,
                None => ordinal = run.end_ordinal(),
            }
        }
        false
    }
}

/// The masks of the sequences of one block of a volume for one filtering algorithm
/// ([`SeqDb::read_mask_block`]), kept with their buffers from one block to the next (a search's
/// worker keeps one). No NCBI counterpart: NCBI keeps one `TSequenceRanges` per search thread
/// (`SSeqDB_SeqSrcData::seq_ranges`) and points into its mapping of the column.
#[derive(Debug, Default)]
pub struct MaskBlock {
    /// The block and the algorithm read.
    read: Option<(BlockRange, i32)>,
    /// The error of every sequence of the block: the algorithms' table or the volume's columns.
    error: Option<SeqDbError>,
    /// Whether the volume has a mask column (none: no sequence has masks).
    column: bool,
    /// The volume's own id of the algorithm (`None`: the volume does not have it).
    volume_algo: Option<i32>,
    /// The column's data file, which the errors of the ranges name.
    path: PathBuf,
    /// The block's blobs.
    blobs: super::column::ColumnBlobs,
    /// The ranges of the last sequence asked for, and the passes of `read_ranges_into`.
    ranges: Vec<(u32, u32)>,
    pass_starts: std::collections::HashSet<i64>,
}

/// The masks of a sequence of a [`MaskBlock`].
#[derive(Debug)]
pub enum BlockMasks<'a> {
    /// The ranges, as stored (none for a sequence without mask data).
    Ranges(&'a [(u32, u32)]),
    /// The sequence has mask data but its volume does not know the algorithm
    /// ([`MaskRanges::AlgorithmNotInVolume`]).
    AlgorithmNotInVolume,
}

/// Why a [`MaskBlock`] has no masks of a sequence.
#[derive(Debug)]
pub enum BlockMaskError<'a> {
    /// The read of the block or of the sequence's blob (the block keeps it).
    Read(&'a SeqDbError),
    /// The ranges of the sequence's blob.
    Ranges(SeqDbError),
}

impl BlockMaskError<'_> {
    /// The error.
    pub fn error(&self) -> &SeqDbError {
        match self {
            BlockMaskError::Read(error) => error,
            BlockMaskError::Ranges(error) => error,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2531-2546
// ```c++
//     CBlastDbBlob blob;
//     vol->GetColumnBlob(x_GetMaskDataColumn(locked), vol_oid, blob, false, locked);
//
//     if (blob.Size() != 0) {
//         // If there actually is mask data, then we need to do the
//         // algorithm translation.
//
//         int vol_algo_id = -1;
//
//         try {
//             vol_algo_id = m_AlgorithmIds.GetVolAlgo(vol_idx, algo_id);
//         }
//         catch(CSeqDBException &) {
//             return;
//         }
//
//         s_ReadRanges<SReadInt4>(vol_algo_id, ranges, blob);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321
// ```c++
//     if (datap->mask_type != eNoSubjMasking) {
//         ASSERT(datap->mask_algo_id != -1);
//         seqdb.GetMaskData(oid, datap->mask_algo_id, datap->seq_ranges);
//     }
// ```
impl MaskBlock {
    /// `GetMaskData` of the sequence `vol_oid` (volume-relative) of the block read: its ranges
    /// of the algorithm in the order they are stored, valid until the next call.
    pub fn ranges(&mut self, vol_oid: u32) -> Result<BlockMasks<'_>, BlockMaskError<'_>> {
        if let Some(error) = &self.error {
            return Err(BlockMaskError::Read(error));
        }
        if !self.column {
            return Ok(BlockMasks::Ranges(&[]));
        }
        let blob = self.blobs.blob(vol_oid).map_err(BlockMaskError::Read)?;
        if blob.is_empty() {
            return Ok(BlockMasks::Ranges(&[]));
        }
        let Some(volume_algo) = self.volume_algo else {
            return Ok(BlockMasks::AlgorithmNotInVolume);
        };
        self.ranges.clear();
        read_ranges_into(
            blob,
            volume_algo,
            &self.path,
            &mut self.ranges,
            &mut self.pass_starts,
        )
        .map_err(BlockMaskError::Ranges)?;
        Ok(BlockMasks::Ranges(&self.ranges))
    }

    /// The read error that [`MaskBlock::ranges`] gives the sequence `vol_oid` (`Read`), owned.
    fn into_read_error(self, vol_oid: u32) -> SeqDbError {
        match self.error {
            Some(error) => error,
            None => self.blobs.into_error(vol_oid),
        }
    }
}
