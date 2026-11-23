//! The mask data of a BLAST database (`makeblastdb -mask_data`): the filtering algorithms
//! that were applied to the sequences, and the masked ranges of each sequence.
//!
//! NCBI reads them in `CSeqDBImpl` (`c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp`: `x_BuildMaskAlgorithmList`,
//! `GetMaskAlgorithmDetails`, `GetAvailableMaskAlgorithmDescriptions`, `GetMaskData`) with
//! the id table `CSeqDB_IdRemapper`. The data is in the column titled `BlastDb/MaskData` of
//! each volume ([`super::column`]): its metadata maps the volume's algorithm ids to
//! descriptions; the blob of a sequence holds, per algorithm, a list of ranges.
//!
//! A database whose alias file names GI-indexed masks (`MASKLIST`, `makeblastdb -gi_mask`) has
//! those instead (`gi_mask`, NCBI's `m_UseGiMask`): [`MaskAlgorithms`] dispatches between the two
//! once, and the reads below take the masks of either.
//!
//! The reader of `losat blastdbcmd`, and of BLASTN's `-db_soft_mask` and `-db_hard_mask`
//! (`blastinput::db_mask`): NCBI's searches read the columns only with a masking algorithm named
//! by `-db_soft_mask`, `-db_hard_mask` or an imported search strategy (LOSAT refuses the last),
//! and so do LOSAT's. The columns are opened on first use.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ncbi_str::string_to_int;

use super::column::{Blob, MASK_DATA_TITLE};
use super::gi_mask::GiMasks;
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
// `StringToNumeric` (with `str`, the address of an `int`, `flags` and `base`) is
// `NStr::StringToInt` that returns whether it converted the text (true when the value is not 0 or
// `errno` is 0; a converted 0 leaves `errno` at 0).
/// `s_IsNumericId` of `id` (`blastinput::db_mask::FilteringAlgorithm::parse` is the same test of
/// the searches, `uniform_search.cpp:149-153`).
fn is_numeric_id(id: &str) -> bool {
    string_to_int(id.as_bytes(), 0).is_ok()
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-303
/// `NStr::Split` of `s` at `":"` into `items` with flags 0: every `:` separates, empty tokens are
/// kept (so `"30:"` is `["30", ""]`), and the empty string has no tokens.
pub(super) fn ncbi_split(s: &str) -> Vec<&str> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.split(':').collect()
    }
}

// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-59
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

/// Pads `text` with spaces on the right to `width` bytes (`setw` of `width` with `left` alignment;
/// text that is longer is not cut): the one helper of the names tables of the mask columns and of
/// the GI masks (`gi_mask`).
pub(super) fn push_padded(out: &mut String, text: &str, width: usize) {
    out.push_str(text);
    for _ in text.len()..width {
        out.push(' ');
    }
}

/// The filtering algorithms of a database (`GetAvailableMaskAlgorithms`,
/// `GetMaskAlgorithmDetails`, `GetMaskAlgorithmId`): those of the volumes' mask columns, or the
/// GI masks that its alias file names (NCBI's `m_UseGiMask`), never both.
#[derive(Debug, Clone)]
pub struct MaskAlgorithms {
    kind: AlgorithmsKind,
}

#[derive(Debug, Clone)]
enum AlgorithmsKind {
    /// The mask columns of the volumes.
    Columns(ColumnAlgorithms),
    /// The GI masks of `MASKLIST`.
    Gi(Arc<GiMasks>),
}

/// The filtering algorithms of the mask columns: `CSeqDB_IdRemapper` with the descriptions.
///
/// The ids that volumes use for an algorithm may differ from volume to volume; the table
/// gives each description one database-wide id and remembers each volume's own.
#[derive(Debug, Clone)]
struct ColumnAlgorithms {
    /// `m_NextId`: where the search for an unused id starts.
    next_id: i32,
    /// `m_IdToDesc`, with the file that first stated the description (for messages).
    id_to_desc: BTreeMap<i32, (String, PathBuf)>,
    /// `m_DescToId`, keyed by the description (item 2 of a four-item one).
    desc_to_id: BTreeMap<String, i32>,
    /// `m_RealIdToVolumeId`: by volume index, database-wide id to the volume's id.
    real_to_volume: BTreeMap<usize, BTreeMap<i32, i32>>,
}

impl ColumnAlgorithms {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2167-2171
    fn new() -> Self {
        Self {
            next_id: 100,
            id_to_desc: BTreeMap::new(),
            desc_to_id: BTreeMap::new(),
            real_to_volume: BTreeMap::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2182-2232
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
    /// Builds the table from the mask columns of `volumes` (in database order). Every
    /// volume opens its columns first, as `x_GetColumnId` does, so an error in a later
    /// volume wins over one in an earlier volume's data.
    fn build(volumes: &[Volume]) -> Result<Self, SeqDbError> {
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2244-2259
    /// `RealToVol`: the id the volume `volume` (index in the database) uses for the
    /// database-wide `algo_id`; `None` where NCBI throws ("Cannot find volume in algorithm
    /// map." or "Cannot find volume algorithm in algorithm map.").
    fn volume_algorithm(&self, volume: usize, algo_id: i32) -> Option<i32> {
        self.real_to_volume.get(&volume)?.get(&algo_id).copied()
    }
}

impl MaskAlgorithms {
    /// The table of the mask columns of `volumes` (`x_BuildMaskAlgorithmList`,
    /// `ColumnAlgorithms::build`).
    pub(super) fn build(volumes: &[Volume]) -> Result<Self, SeqDbError> {
        Ok(Self {
            kind: AlgorithmsKind::Columns(ColumnAlgorithms::build(volumes)?),
        })
    }

    /// The algorithms of the GI masks `masks` (`m_UseGiMask`).
    pub(super) fn gi(masks: Arc<GiMasks>) -> Self {
        Self {
            kind: AlgorithmsKind::Gi(masks),
        }
    }

    /// The GI masks, where the database has them instead of mask columns.
    pub(super) fn gi_masks(&self) -> Option<&GiMasks> {
        match &self.kind {
            AlgorithmsKind::Columns(_) => None,
            AlgorithmsKind::Gi(masks) => Some(masks),
        }
    }

    /// `RealToVol` of the mask columns (`ColumnAlgorithms::volume_algorithm`); `None` for GI masks.
    fn volume_algorithm(&self, volume: usize, algo_id: i32) -> Option<i32> {
        match &self.kind {
            AlgorithmsKind::Columns(table) => table.volume_algorithm(volume, algo_id),
            AlgorithmsKind::Gi(_) => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2173-2180
    /// The database-wide algorithm ids, ascending (for GI masks 0 to N-1, `GiMasks::ids`).
    pub fn ids(&self) -> Vec<i32> {
        match &self.kind {
            AlgorithmsKind::Columns(table) => table.id_to_desc.keys().copied().collect(),
            AlgorithmsKind::Gi(masks) => masks.ids(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2261-2268
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2292-2307
    /// `GetMaskAlgorithmId`: the id of the algorithm with this description (item 2 of a
    /// four-item description; otherwise the whole description, for example `"30:"`), or the
    /// text of NCBI's exception: "Cannot find string algorithm id in algorithm map."; for GI
    /// masks the id of the `MASKLIST` name (`GiMasks::algorithm_id`, whose text lists them).
    pub fn id_of(&self, name: &str) -> Result<i32, String> {
        match &self.kind {
            AlgorithmsKind::Columns(table) => table
                .desc_to_id
                .get(name)
                .copied()
                .ok_or_else(|| "Cannot find string algorithm id in algorithm map.".to_string()),
            AlgorithmsKind::Gi(masks) => masks.algorithm_id(name),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2380-2412
    /// `GetMaskAlgorithmDetails`. An id the database does not have gives NCBI's text
    /// ("Filtering algorithm ID `<id>` is not supported." and the list of algorithms) in
    /// [`SeqDbError::Range`]. For GI masks the description is the one of the algorithm's index
    /// file, which is opened here (`GiMasks::description`: its errors, and for an id outside 0 to
    /// N-1 "Filtering algorithm ID `<id>` does not exist." and the list of names).
    pub fn details(&self, id: i32) -> Result<MaskAlgorithm, SeqDbError> {
        match &self.kind {
            AlgorithmsKind::Columns(table) => {
                let Some((description, source)) = table.id_to_desc.get(&id) else {
                    return Err(SeqDbError::Range(format!(
                        "Filtering algorithm ID {id} is not supported.\n{}",
                        self.descriptions()?
                    )));
                };
                parse_details(id, description, source)
            }
            AlgorithmsKind::Gi(masks) => {
                let (description, source) = masks.description(id)?;
                parse_details(id, &description, &source)
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2309-2342
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
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2344-2378
/// `s_GetDetails` of the stored description `description` of the algorithm `id` (`source`: the
/// file that holds it, for LOSAT's errors).
fn parse_details(id: i32, description: &str, source: &Path) -> Result<MaskAlgorithm, SeqDbError> {
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
            // integer; `FindName` with `pid` and false throws on a value the enum lacks.
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2470-2506
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3276-3292
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
            // `ReadRaw` of 8 times `n` bytes: a negative count is NCBI's "hit end of data".
            let bytes = blob.read_raw(i64::from(num_pairs) * 8).map_err(format)?;
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1418-1423
            // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:302-306
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
    /// database without mask data. With GI masks, their table, or LOSAT's refusal of a
    /// `MASKLIST` name (`alias::gi_mask_list`).
    pub fn mask_algorithms(&self) -> Result<&MaskAlgorithms, SeqDbError> {
        if let Some(table) = self.mask_algorithms.get() {
            return Ok(table);
        }
        if let Some(refusal) = &self.gi_mask_refusal {
            return Err(refusal.duplicate());
        }
        let table = MaskAlgorithms::build(&self.volumes)?;
        Ok(self.mask_algorithms.get_or_init(|| table))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2508-2547
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
        self.mask_ranges_in(at, algo_id, &mut MaskBlock::default())
    }

    /// [`SeqDb::mask_ranges_at`] through the caller's `block`, whose buffers (and the name of the
    /// column's file) are kept from one sequence to the next: the reader of a program that asks
    /// for the masks of one sequence after another (`blastdbcmd`'s `%m` and
    /// `-mask_sequence_with`, `makembindex -db_mask`).
    pub fn mask_ranges_in(
        &self,
        at: SeqLocation,
        algo_id: i32,
        block: &mut MaskBlock,
    ) -> Result<MaskRanges, SeqDbError> {
        let SeqLocation {
            volume, vol_oid, ..
        } = at;
        self.read_mask_block(
            BlockRange {
                volume,
                first: vol_oid,
                end: vol_oid + 1,
                first_ordinal: at.id,
            },
            algo_id,
            block,
        );
        match block.ranges(vol_oid) {
            Ok(BlockMasks::NoMaskData) => return Ok(MaskRanges::Ranges(Vec::new())),
            Ok(BlockMasks::Ranges(ranges)) => return Ok(MaskRanges::Ranges(ranges.to_vec())),
            Ok(BlockMasks::AlgorithmNotInVolume) => {
                return Ok(MaskRanges::AlgorithmNotInVolume {
                    volume: self.volumes[volume].base().to_path_buf(),
                })
            }
            Err(BlockMaskError::Ranges(error)) => return Err(error),
            Err(BlockMaskError::Read(_)) => {}
        }
        Err(std::mem::take(block).into_read_error(vol_oid))
    }

    // No NCBI counterpart (reads only): NCBI maps the column files and its `GetMaskData` points into
    // them for each sequence (seqdbimpl.cpp:2508-2547, cited at `mask_ranges_at`).
    /// The masks of the sequences of `range` for the filtering algorithm `algo_id` into `block`
    /// (unless it holds them already), as a search reads them with the sequences of its blocks
    /// (`SeqDb::read_block`): one read of the mask column's offset array for the block and one of
    /// its blobs (`Column::read_blobs`); [`MaskBlock::ranges`] gives each sequence's, with
    /// [`SeqDb::mask_ranges_at`]'s checks and errors in its order (the algorithms' table and the
    /// volume's columns: the error of every sequence of the block).
    #[inline]
    pub fn read_mask_block(&self, range: BlockRange, algo_id: i32, block: &mut MaskBlock) {
        if block.read != Some((range, algo_id)) {
            self.fill_mask_block(range, algo_id, block);
        }
    }

    /// [`SeqDb::read_mask_block`] of a block that `block` does not hold: the algorithms' table,
    /// then the volume's mask column or the GI masks (apart, so that the check of the block held
    /// stays a short call in the caller).
    #[inline(never)]
    fn fill_mask_block(&self, range: BlockRange, algo_id: i32, block: &mut MaskBlock) {
        block.read = Some((range, algo_id));
        block.error = None;
        block.source = BlockSource::None;
        block.volume_algo = None;
        let algorithms = match self.mask_algorithms() {
            Ok(algorithms) => algorithms,
            Err(error) => {
                block.error = Some(error);
                return;
            }
        };
        let table = match &algorithms.kind {
            AlgorithmsKind::Columns(table) => table,
            AlgorithmsKind::Gi(masks) => {
                self.read_gi_mask_block(masks, range, algo_id, block);
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
        block.source = BlockSource::Column;
        // If there actually is mask data, then we need to do the algorithm translation (the
        // volume's id, read for the blobs that are not empty).
        block.volume_algo = table.volume_algorithm(range.volume, algo_id);
        if block.path.as_path() != column.data_path() {
            block.path = column.data_path().to_path_buf();
        }
        column.read_blobs(range.first, range.end, &mut block.blobs);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2520-2523
    // For each sequence the GI comes first, then `x_Open` of the algorithm, then the lookup.
    // NCBI keeps the files once opened and opens them again at the next sequence after a failed
    // opening, which fails again on the same files: the block opens them at its first sequence
    // whose GI is read, and a failure is the error of every later sequence whose GI is read.
    /// [`SeqDb::read_mask_block`] of a database with GI masks: the ranges of every sequence of
    /// `range`, or its error, into `block`.
    fn read_gi_mask_block(
        &self,
        masks: &GiMasks,
        range: BlockRange,
        algo_id: i32,
        block: &mut MaskBlock,
    ) {
        block.source = BlockSource::Gi;
        block.ranges.clear();
        block.gi_spans.clear();
        let volume = &self.volumes[range.volume];
        let mut files = None;
        for vol_oid in range.first..range.end {
            let start = block.ranges.len();
            let read = masks
                .seq_gi_key(range.volume, volume, vol_oid)
                .and_then(
                    |key| match files.get_or_insert_with(|| masks.open(algo_id)) {
                        Ok(files) => {
                            files.ranges_into(key, &mut block.ranges, &mut block.gi_window)
                        }
                        Err(error) => Err(error.duplicate()),
                    },
                );
            let span = match read {
                Ok(()) => Ok(start..block.ranges.len()),
                Err(error) => {
                    block.ranges.truncate(start);
                    Err(error)
                }
            };
            block.gi_spans.push(span);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:476-492
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:376-406
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3345-3358
    // (cited at `column::open_all`): a little-endian host reads a column from its `.?xc` file,
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
        if algorithms.gi_masks().is_some() {
            // GI masks: no volume algorithm to translate, no byte order to choose.
            return Ok(());
        }
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321 (cited at
    // `MaskBlock::ranges`): the search reads the masks of the sequences it searches only.
    /// The first sequence that the search reads (by ordinal id: those the list and the alias
    /// files' ranges include) whose masks of the filtering algorithm `algo_id`, as stored, `test`
    /// finds a refusal in, with its NCBI OID (as `blastdbcmd`'s `%o`) and its volume's path.
    /// `test` gets the masks (never empty), the bounds of the sequence's length that the index
    /// gives (`Volume::nucleotide_length_bounds`, or a protein's exact length,
    /// `Volume::x_protein_lengths`, read once per block at its first sequence with
    /// masks; `None` where they give none or cannot be read) and a function that reads the
    /// sequence's length (the sequence's index and its last byte; `None` when they cannot be
    /// read, and `test` then gives `Err(())`), which `test` calls only where the bounds do not
    /// decide (review 2 S4 of M4a-2b: three reads per sequence were most of the walk). No NCBI
    /// counterpart: a check of the masks before the search for LOSAT's refusals
    /// (`blastinput::db_mask`, series NP, M4a-2b), with `-db_hard_mask` only. The masks are read
    /// 4,096 sequences at a time (`SeqDb::read_mask_block`); a sequence whose masks or length
    /// cannot be read ends the walk without a finding (the search meets that read error first,
    /// as `searches_a_blob_of`).
    pub fn first_searched_masks<T>(
        &self,
        algo_id: i32,
        mut test: impl FnMut(
            &[(u32, u32)],
            Option<(usize, usize)>,
            &dyn Fn() -> Option<usize>,
        ) -> Result<Option<T>, ()>,
    ) -> Option<(u32, PathBuf, T)> {
        const WALK_BLOCK: u32 = 4096;
        let protein = self.seqtype() == super::SeqType::Protein;
        let mut block = MaskBlock::default();
        let mut ordinal = 0u32;
        let end_ordinal = self.num_seqs();
        while let Some(run) =
            self.x_block_range(ordinal, end_ordinal.min(ordinal.saturating_add(WALK_BLOCK)))
        {
            self.read_mask_block(run, algo_id, &mut block);
            // The block's length bounds, read at its first sequence with masks (none where the
            // index cannot be read: each length is then read alone, with its own error).
            let mut bounds: Option<Vec<Option<(usize, usize)>>> = None;
            for vol_oid in run.first..run.end {
                let seq_ordinal = run.first_ordinal + (vol_oid - run.first);
                let masks = match block.ranges(vol_oid) {
                    Ok(BlockMasks::Ranges(masks)) => masks,
                    // No masks; an algorithm not in the volume is refused before
                    // (`check_search_masks`).
                    Ok(BlockMasks::NoMaskData | BlockMasks::AlgorithmNotInVolume) => continue,
                    Err(_) => return None,
                };
                if masks.is_empty() {
                    continue;
                }
                let volume = &self.volumes[run.volume];
                let length_bounds = bounds
                    .get_or_insert_with(|| {
                        if protein {
                            // A protein's length is exact in its index (series NP, M4b).
                            volume
                                .x_protein_lengths(run.first, run.end)
                                .map(|lengths| {
                                    lengths
                                        .into_iter()
                                        .map(|length| Some((length as usize, length as usize)))
                                        .collect()
                                })
                                .unwrap_or_default()
                        } else {
                            volume
                                .nucleotide_length_bounds(run.first, run.end)
                                .unwrap_or_default()
                        }
                    })
                    .get((vol_oid - run.first) as usize)
                    .copied()
                    .flatten();
                let length = || {
                    if protein {
                        self.protein_length(seq_ordinal)
                            .ok()
                            .map(|length| length as usize)
                    } else {
                        self.nucleotide_length(seq_ordinal).ok()
                    }
                };
                match test(masks, length_bounds, &length) {
                    Ok(None) => {}
                    Ok(Some(found)) => {
                        return Some((
                            self.volume_starts[run.volume] + vol_oid,
                            self.volumes[run.volume].base().to_path_buf(),
                            found,
                        ))
                    }
                    Err(()) => return None,
                }
            }
            ordinal = run.end_ordinal();
        }
        None
    }

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
    /// Where the masks were read from: none (no sequence has masks), the volume's mask column, or
    /// the GI masks.
    source: BlockSource,
    /// The volume's own id of the algorithm (`None`: the volume does not have it).
    volume_algo: Option<i32>,
    /// The column's data file, which the errors of the ranges name.
    path: PathBuf,
    /// The block's blobs.
    blobs: super::column::ColumnBlobs,
    /// The ranges of the last sequence asked for, and the passes of `read_ranges_into`.
    ranges: Vec<(u32, u32)>,
    pass_starts: std::collections::HashSet<i64>,
    /// With GI masks (`SeqDb::read_gi_mask_block`), each sequence's part of `ranges` (which holds
    /// the ranges of all the block's sequences), or its error.
    gi_spans: Vec<Result<std::ops::Range<usize>, SeqDbError>>,
    /// With GI masks, the bytes of the data file read last (`gi_mask::DataWindow`).
    gi_window: super::gi_mask::DataWindow,
}

/// Where a [`MaskBlock`] read its masks from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum BlockSource {
    /// Nowhere: the volume has no mask column (no sequence has masks), or the read failed.
    #[default]
    None,
    /// The volume's mask column.
    Column,
    /// The GI masks of `MASKLIST`.
    Gi,
}

/// The masks of a sequence of a [`MaskBlock`].
#[derive(Debug)]
pub enum BlockMasks<'a> {
    /// The sequence has no mask data (its volume has no mask column, or the sequence's blob is
    /// empty): no ranges, and NCBI does not look the algorithm up in the volume.
    NoMaskData,
    /// The ranges of the sequence's mask data, as stored; with GI masks the ranges of its GI
    /// (none for a GI without masks: NCBI never looks the algorithm up in a volume then).
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
    /// The ranges of the sequence's blob; or a LOSAT internal error ([`SeqDbError::Internal`]):
    /// a sequence outside the block read.
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2531-2555
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321
impl MaskBlock {
    /// `GetMaskData` of the sequence `vol_oid` (volume-relative) of the block read: its ranges
    /// of the algorithm in the order they are stored, valid until the next call.
    pub fn ranges(&mut self, vol_oid: u32) -> Result<BlockMasks<'_>, BlockMaskError<'_>> {
        if let Some(error) = &self.error {
            return Err(BlockMaskError::Read(error));
        }
        match self.source {
            BlockSource::Column => {}
            BlockSource::None => return Ok(BlockMasks::NoMaskData),
            BlockSource::Gi => {
                let first = self.read.map_or(0, |(range, _)| range.first);
                return match self.gi_spans.get(vol_oid.wrapping_sub(first) as usize) {
                    Some(Ok(span)) => Ok(BlockMasks::Ranges(&self.ranges[span.clone()])),
                    Some(Err(error)) => Err(BlockMaskError::Read(error)),
                    None => Err(BlockMaskError::Ranges(SeqDbError::Internal(format!(
                        "LOSAT internal error: the masks of sequence {vol_oid} were asked for \
                         outside the block read"
                    )))),
                };
            }
        }
        let blob = self.blobs.blob(vol_oid).map_err(BlockMaskError::Read)?;
        if blob.is_empty() {
            return Ok(BlockMasks::NoMaskData);
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
    pub(crate) fn into_read_error(self, vol_oid: u32) -> SeqDbError {
        if let Some(error) = self.error {
            return error;
        }
        if self.source == BlockSource::Gi {
            let first = self.read.map_or(0, |(range, _)| range.first);
            let at = vol_oid.wrapping_sub(first) as usize;
            if let Some(Err(error)) = self.gi_spans.into_iter().nth(at) {
                return error;
            }
            return SeqDbError::Internal(format!(
                "LOSAT internal error: no GI-mask read error for sequence {vol_oid}"
            ));
        }
        self.blobs.into_error(vol_oid)
    }
}
