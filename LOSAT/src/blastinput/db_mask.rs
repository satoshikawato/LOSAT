//! `-db_soft_mask` and `-db_hard_mask`: the masks of a BLAST database's sequences that a search
//! applies (NCBI's `CBlastDatabaseArgs` with `SetDatabaseMaskingSupport` set to true and the filtering
//! algorithm of `CSearchDatabase`), from the option to the algorithm the search and its reports
//! use. One resolver for NCBI's whole sequence (decision D2 of M4a-2), in NCBI's order:
//!
//! 1. the command line keeps the option's text (`DbMaskArgs::option`): an id when it is an `Int4`
//!    in NCBI's form, else a name;
//! 2. the database is opened: an id above 0 that the database does not have is NCBI's options
//!    error (`DbMaskOption::validate_at_open`);
//! 3. after `Query is Empty!`, the formatter: a name that the database does not know and the id -1
//!    are NCBI's warning and a search without masks; any other id must be in the table of each
//!    name of the `-db` list, else NCBI's database error (`DbMaskOption::resolve`). The result
//!    (`ResolvedMask`) is the search's algorithm, the program options of the BLAST archive and of
//!    the search strategy, and the `Masked using:` lines of the pairwise report.
//!
//! LOSAT's BLASTN (series NP, M4a-2a and M4a-2b), BLASTP and BLASTX (series NP, M4b) search with
//! `-db_soft_mask` and `-db_hard_mask`, as do TBLASTN and TBLASTX (M4b PR 2); they refuse the
//! databases of `SeqDb::check_search_masks` and, with `-db_hard_mask`, a searched sequence of a
//! `HardMaskShape` (`ResolvedMask::check_losat_limits`). psiblast, rpsblast and rpstblastn have
//! neither option, as in NCBI.

use crate::wording::{DATABASE_ERROR, OPTIONS_ERROR};
use anyhow::Result;

use std::ffi::OsString;

use crate::algorithm::common::hard_chunk_walk::{HardStarts, NcbiStart, NcbiWalk};
use crate::cli::NativeError;
use crate::core::subject_split::HardSplitState;
use crate::report::archive::fields;
use crate::report::archive::options::{Param, ParamValue};
use crate::seqdb::{MaskAlgorithm, NameMaskError, SeqDb, SeqDbError};
use crate::utils::os_bytes::os_str_bytes;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2342-2356
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2226-2229,2367-2382
// NCBI's argument reader reports the exclusions (`cli::reader`).
/// The two options of a program whose database arguments support masking (BLASTN).
#[derive(Debug, Clone, Default)]
pub struct DbMaskArgs {
    /// `-db_soft_mask`.
    pub db_soft_mask: Option<OsString>,
    /// `-db_hard_mask` (refused for the masks of `HardMaskShape`: `ResolvedMask::check_losat_limits`).
    pub db_hard_mask: Option<OsString>,
}

impl DbMaskArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, the arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(DbMaskArgs {
            db_soft_mask: args.given_option("db_soft_mask")?,
            db_hard_mask: args.given_option("db_hard_mask")?,
        })
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:234-239
/// `ESubjectMaskingType` of a masked search (`eNoSubjMasking` is no `DbMaskOption`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaskKind {
    /// `eSoftSubjMasking` (1): the masks keep words of the database sequences out of the scan.
    Soft = 1,
    /// `eHardSubjMasking` (2): the masked letters are cut out of the search.
    Hard = 2,
}

/// What the option's text names (`CSearchDatabase::SetFilteringAlgorithm` with the string and the type).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FilteringAlgorithm {
    /// `m_FilteringAlgorithmId`: a text that is an `Int4` in NCBI's form (`s_IsNumericId`).
    Id(i32),
    /// `m_FilteringAlgorithmString`: any other text, looked up by the database
    /// (`GetMaskAlgorithmId`, case-sensitive).
    Name(Vec<u8>),
}

/// The masking of a `-db` search that the command line asks for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DbMaskOption {
    /// The masking type.
    pub kind: MaskKind,
    /// The option's algorithm.
    pub algorithm: FilteringAlgorithm,
}

impl DbMaskArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2517-2523
    /// The masking that the options ask for (the soft one first; the parser refuses both), for a
    /// search of `-db` (the parser refuses them with `-subject` and `-subject_loc`).
    pub fn option(&self) -> Option<DbMaskOption> {
        let (kind, text) = match (&self.db_soft_mask, &self.db_hard_mask) {
            (Some(text), _) => (MaskKind::Soft, text),
            (None, Some(text)) => (MaskKind::Hard, text),
            (None, None) => return None,
        };
        let bytes = os_str_bytes(text);
        Some(DbMaskOption {
            kind,
            algorithm: FilteringAlgorithm::parse(&bytes),
        })
    }
}

impl FilteringAlgorithm {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:149-171
    // `x_ValidateMaskingAlgorithm` returns at once here: the database is not open yet.
    /// The algorithm that the option's text names: `+11` and `011` are the id 11; ` 11`, `11 `,
    /// `11x`, `2147483648` and the empty text are names.
    pub fn parse(text: &[u8]) -> Self {
        // `s_IsNumericId` is `StringToInt` that converts (`seqdb::mask::is_numeric_id`).
        match crate::ncbi_str::string_to_int(text, 0) {
            Ok(id) => FilteringAlgorithm::Id(id),
            Err(_) => FilteringAlgorithm::Name(text.to_vec()),
        }
    }
}

/// The masking of a search as the formatter resolved it (`DbMaskOption::resolve`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedMask {
    /// The masking type.
    pub kind: MaskKind,
    /// The database-wide id of the algorithm (`GetFilteringAlgorithm`).
    pub id: i32,
    /// `GetFilteringAlgorithmKey`: the option's text where it was a name.
    pub key: Option<Vec<u8>>,
    /// The algorithm in the database of each name of the `-db` list, in the order of the list
    /// (`SDbInfo::filt_algorithm_name`, `filt_algorithm_options`).
    pub shown: Vec<MaskAlgorithm>,
    /// With `-db_hard_mask` in BLASTP, BLASTX, TBLASTN and TBLASTX: the searched sequences that
    /// the search starts on a new table at NCBI's one-thread offset (or apart for walked
    /// translated sequences), from the walk of `ResolvedMask::check_losat_limits` (shared by the
    /// clones of the masking).
    pub(crate) starts: std::sync::Arc<std::sync::OnceLock<HardStarts>>,
}

impl DbMaskOption {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:295-318
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-220
    /// NCBI's check of an id when it opens the database `db` of the `-db` value `db_name`
    /// (`x_InitializeDb`): an id above 0 must be one of the database's algorithms (all its
    /// volumes'), else NCBI's options error, exit 1. A name is not looked up here, nor an id of
    /// 0 or less. `program` names the program in LOSAT's refusals.
    pub fn validate_at_open(
        &self,
        db: &SeqDb,
        db_name: &[u8],
        nucleotide: bool,
        program: &str,
    ) -> Result<()> {
        let FilteringAlgorithm::Id(id) = self.algorithm else {
            return Ok(());
        };
        if id <= 0 {
            return Ok(());
        }
        let algorithms = db
            .mask_algorithms()
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        if !algorithms.ids().contains(&id) {
            let kind = if nucleotide { "nucleotide" } else { "protein" };
            let mut message =
                format!("{OPTIONS_ERROR} Masking algorithm ID {id} is not supported in {kind} '")
                    .into_bytes();
            message.extend_from_slice(db_name);
            message.extend_from_slice(b"' BLAST database\n");
            return Err(NativeError { exit: 1, message }.into());
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:195-231
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:228-242
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:196-200
    /// The formatter's resolution, after `Query is Empty!`: the id of a name (by the table of all
    /// the database's volumes, `MaskAlgorithms::id_of`). A name that the database does not have
    /// (or that NCBI's reader of the table fails on) and the id -1 write NCBI's warning to
    /// `warning` and give `None`: the search is not masked. Any other id is looked up in the
    /// database of each name of the `-db` list (`SeqDb::name_mask_details`, the tables that the
    /// database's name records keep): an id that one of them does not have is NCBI's database
    /// error, exit 2 (the id 0, an id below -1, or a name's id that a database of the list
    /// lacks). `app` is the running application, the prefix of the warning.
    pub fn resolve(
        &self,
        db: &SeqDb,
        db_name: &[u8],
        program: &str,
        app: &str,
        warning: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<Option<ResolvedMask>> {
        let (id, key) = match &self.algorithm {
            FilteringAlgorithm::Id(id) => (Some(*id), None),
            FilteringAlgorithm::Name(name) => {
                let id = match db.mask_algorithms() {
                    Ok(algorithms) => algorithms.id_of(name).ok(),
                    // NCBI catches its `CSeqDBException` and answers -1.
                    Err(error) if error.ncbi_message().is_some() => None,
                    Err(error) => {
                        return Err(crate::algorithm::common::database::database_error(
                            error, program,
                        ))
                    }
                };
                (id, Some(name.clone()))
            }
        };
        let Some(id) = id.filter(|&id| id != -1) else {
            let mut line = format!("Warning: [{app}] Subject mask not found in ").into_bytes();
            line.extend_from_slice(db_name);
            line.extend_from_slice(b", proceeding without subject masking.\n");
            warning(&line)?;
            return Ok(None);
        };
        let shown = db
            .name_mask_details(id)
            .map_err(|error| name_mask_error(&error, program))?;
        Ok(Some(ResolvedMask {
            kind: self.kind,
            id,
            key,
            shown,
            starts: Default::default(),
        }))
    }
}

impl DbMaskOption {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:644-651
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:54-58
    // The search strategy reads the command line's `CSearchDatabase`; the formatter's warning
    // resets the database adapter's copy (blast_format.cpp:146), which the archive reads. So
    // after the warning the strategy still writes the name as its key (search_strategy.cpp:500-505,
    // cited at `ResolvedMask::archive_params`); an empty name and the id -1 write nothing.
    /// The search strategy's program options of the masking (`-export_search_strategy`): the
    /// archive's (`ResolvedMask::archive_params`) when the formatter resolved it (`resolved`), else,
    /// after its warning, `DbFilteringAlgorithmKey` and `SubjectMaskingType` of a name that is not
    /// empty.
    pub fn strategy_params(&self, resolved: Option<&ResolvedMask>) -> Vec<Param> {
        match (resolved, &self.algorithm) {
            (Some(mask), _) => mask.archive_params(),
            (None, FilteringAlgorithm::Name(name)) if !name.is_empty() => vec![
                Param::new(
                    fields::DB_FILTERING_ALGORITHM_KEY,
                    ParamValue::Str(name.clone()),
                ),
                Param::new(
                    fields::SUBJECT_MASKING_TYPE,
                    ParamValue::Integer(self.kind as i32),
                ),
            ],
            (None, _) => Vec::new(),
        }
    }
}

impl ResolvedMask {
    /// LOSAT's refusals of the masked search of `db`, where NCBI starts the search (after its own
    /// checks; `program` names the program), before anything is written: the databases whose masks
    /// NCBI reads from undefined or reversed data (`SeqDb::check_search_masks`; decision D5 of
    /// M4a-2), then with `-db_hard_mask` the first searched sequence whose masks NCBI's chunks
    /// do not handle (`HardMaskShape`; series NP, M4a-2b decisions H3-H5, and M4a-2c: a mask that
    /// runs to the end only below `TAIL_MASK_SEARCHED_WORD`, in BLASTN; series NP, M4b-3: no such
    /// mask in the protein and translated programs), for a search whose scan is `scan`
    /// (`HardMaskShape::of_scan`). The protein programs' walk visits every searched sequence and
    /// keeps the sequences that the search starts on a new table (`ResolvedMask::hard_starts`).
    pub fn check_losat_limits(&self, db: &SeqDb, program: &str, scan: HardScan) -> Result<()> {
        db.check_search_masks(self.id)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        if self.kind == MaskKind::Hard {
            let found = match scan {
                HardScan::Protein {
                    word,
                    window,
                    zero_words,
                } => {
                    let (starts, found) = protein_walk(db, self.id, word, window, zero_words);
                    if found.is_none() {
                        // A second check of the same masking finds the same table.
                        let _ = self.starts.set(starts);
                    }
                    found
                }
                HardScan::Translated {
                    word,
                    window,
                    zero_words,
                } => {
                    enum FirstPass {
                        Refused(HardMaskShape),
                        Walked,
                    }
                    let first = db.first_searched_masks(self.id, |masks, bounds, length| {
                        if let Some(shape) =
                            HardMaskShape::of_scan(masks, bounds, length, scan, None)?
                        {
                            return Ok(Some(FirstPass::Refused(shape)));
                        }
                        if translated_walked(masks) {
                            return Ok(Some(FirstPass::Walked));
                        }
                        Ok(None)
                    });
                    match first {
                        None => None,
                        Some((oid, volume, FirstPass::Refused(shape))) => {
                            Some((oid, volume, shape))
                        }
                        Some((_, _, FirstPass::Walked)) => {
                            let (starts, found) =
                                translated_walk(db, self.id, word, window, zero_words);
                            if found.is_none() {
                                let _ = self.starts.set(starts);
                            }
                            found
                        }
                    }
                }
                HardScan::Nucleotide { .. } => {
                    db.first_searched_masks(self.id, |masks, bounds, length| {
                        HardMaskShape::of_scan(masks, bounds, length, scan, None)
                    })
                }
            };
            if let Some((oid, volume, shape)) = found {
                anyhow::bail!(
                    "-db_hard_mask: the sequence of OID {oid} (volume {}) {}; LOSAT's {program} \
                     does not search such a database with -db_hard_mask",
                    volume.display(),
                    shape.text(),
                );
            }
        }
        Ok(())
    }

    /// The searched sequences that a BLASTP or BLASTX search of `db` with `-db_hard_mask`, the word
    /// `word` and the window `window` starts on a new table at NCBI's one-thread offset
    /// (`algorithm::common::hard_chunk_walk::NcbiWalk`): the table of `check_losat_limits`, which
    /// the search's start runs, or the same walk where it did not run (a search that a caller
    /// set up without it).
    pub(crate) fn hard_starts(&self, db: &SeqDb, word: i32, window: i32) -> &HardStarts {
        self.starts
            .get_or_init(|| protein_walk(db, self.id, word, window, false).0)
    }

    /// The searched sequences of a TBLASTN or TBLASTX `-db_hard_mask` search that run in order on
    /// their own diagonal table (`translated_walk`): `None` when no sequence in the database had
    /// overlapping masks (no walk ran).
    pub(crate) fn translated_starts(&self) -> Option<&HardStarts> {
        self.starts.get()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:499-516
    /// The program options of the BLAST archive and of the search strategy (after the list's).
    /// An empty name is no key (only an `algo_key` other than `kEmptyStr` is one): the id it was
    /// translated to is written.
    pub fn archive_params(&self) -> Vec<Param> {
        let algorithm = match self.key.as_deref() {
            Some(key) if !key.is_empty() => Param::new(
                fields::DB_FILTERING_ALGORITHM_KEY,
                ParamValue::Str(key.to_vec()),
            ),
            _ => Param::new(
                fields::DB_FILTERING_ALGORITHM_ID,
                ParamValue::Integer(self.id),
            ),
        };
        vec![
            algorithm,
            Param::new(
                fields::SUBJECT_MASKING_TYPE,
                ParamValue::Integer(self.kind as i32),
            ),
        ]
    }
}

/// The error of `SeqDb::name_mask_details`, as the formatter's `GetBlastDbInfo` ends the run with
/// it (`CATCH_ALL`, cited at `DbMaskOption::resolve`): NCBI's `Filtering algorithm ID ... is not
/// supported.` of a name's database that lacks the id, or the error of reading its table, as
/// `LOSAT database error: <text>` with exit status 2 (`common::database::database_error_of`).
/// The one mapping of these errors: the search programs' database masking (`resolve`) and
/// blast_formatter's database of an archive (`blast_formatter::archive`).
pub(crate) fn name_mask_error(error: &NameMaskError<'_>, program: &str) -> anyhow::Error {
    let error = match error {
        NameMaskError::Table(error) => *error,
        NameMaskError::Details(error) => error,
        NameMaskError::NotRead => {
            return anyhow::anyhow!(
                "LOSAT internal error: the database's names were read without their mask algorithms"
            )
        }
    };
    match error {
        SeqDbError::Range(text) => NativeError {
            exit: 2,
            message: format!("{DATABASE_ERROR} {text}\n").into(),
        }
        .into(),
        error => crate::algorithm::common::database::database_error_of(error, program),
    }
}

/// The walk of a BLASTP or BLASTX search of `db` with `-db_hard_mask` (the filtering algorithm
/// `algo_id`, the word `word`, the window `window`; `zero_words` as in `HardScan::Protein`) over
/// every sequence that it searches, in the order of the search (`SeqDb::searched_sequences`):
/// NCBI's one-thread offsets (`NcbiWalk`) with the sequences that the search starts on a new
/// table and whether a sequence is kept apart (`HardStarts::keeps_apart`), and the first sequence
/// whose masks LOSAT refuses (`HardMaskShape::of_scan`, at the place that the walk gives it),
/// after which the table is not complete (the search does not run). A sequence whose masks or
/// length cannot be read ends the walk without a finding (the search meets that read error first).
/// With a negative window (a search strategy's; the command line takes none) the search itself runs
/// on NCBI's one-thread table, in NCBI's order (BLASTP's `db_masked::neg_window_chunks`, BLASTX's
/// `search_db_masked.rs`): no sequence starts on a new table or is kept apart, and the walk only
/// finds the refusals.
fn protein_walk(
    db: &SeqDb,
    algo_id: i32,
    word: i32,
    window: i32,
    zero_words: bool,
) -> (HardStarts, Option<(u32, std::path::PathBuf, HardMaskShape)>) {
    let scan = HardScan::Protein {
        word,
        window,
        zero_words,
    };
    let carried = window < 0;
    let mut walk = NcbiWalk::new(word, window);
    let mut starts = HardStarts::default();
    let mut plan = crate::algorithm::common::db_mask_search::SubjectPlan::default();
    let mut chunks = Vec::new();
    let found = db
        .searched_sequences(algo_id, |ordinal, length, masks| {
            let masks = masks.read()?;
            if masks.is_empty() {
                if let Some(start) = walk.plain(length).filter(|_| !carried) {
                    starts.push(ordinal, start);
                }
                return Ok(None);
            }
            // The search's chunker (`db_mask_search::SubjectPlan`), so that the walk and the search
            // read the same chunks.
            plan.start_hard(masks, length);
            plan.hard_chunks(&mut chunks);
            let (start, reads_back) = walk.masked(length, &chunks);
            let (start, reads_back) = if carried {
                (None, false)
            } else {
                (start, reads_back)
            };
            if let Some(start) = start {
                starts.push(ordinal, start);
            } else if reads_back {
                starts.push_reads_back();
            }
            // The length is exact: never read.
            Ok(
                HardMaskShape::of_scan(
                    masks,
                    Some((length, length)),
                    &|| Some(length),
                    scan,
                    start,
                )
                .ok()
                .flatten(),
            )
        })
        .ok()
        .flatten();
    (starts, found)
}

/// The shortest translated subject searched for a word size in residues (three times the word
/// size plus two bases, NCBI's `min_subj_seq_length` in `blast_engine.c:1318-1331`).
fn translated_min_length(word: i32) -> usize {
    usize::try_from(word).unwrap_or(0) * 3 + 2
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 805-830
// Each frame maps nucleotide mask boundaries by dividing by 3, so two consecutive stored masks
// yield a negative frame chunk in every reading frame iff the first mask's floor-divided end
// exceeds the next mask's floor-divided start.
/// Whether a translated sequence with stored masks `masks` has a negative frame chunk (a walked
/// sequence: two consecutive masks whose frame-mapped end and start overlap).
pub(crate) fn translated_walked(masks: &[(u32, u32)]) -> bool {
    masks.windows(2).any(|pair| pair[0].1 / 3 > pair[1].0 / 3)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 478-500, 804-855,
// 1318-1331, 1368-1369, 1429
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107, 162-173
/// The one-thread walk of a TBLASTN or TBLASTX search of `db` with `-db_hard_mask` (series NP,
/// M4b-3c part 2): runs only when the first pass of `ResolvedMask::check_losat_limits` finds a
/// walked sequence (`translated_walked`). Walks every searched sequence in database order with
/// exact nucleotide lengths (`SeqDb::searched_sequences`), checks the first-pass shapes and the
/// walk's refusals in order (`of_scan`, `WordPastEnd`, `NoChunkingWithOverlaps`), and records the
/// `HardStarts` for every walked or New sequence (approved exception 29: a sequence whose hits
/// would reach another sequence's entries runs on a table of its own at NCBI's offset where it
/// meets unwritten entries, else apart at `window + below_start`).
fn translated_walk(
    db: &SeqDb,
    algo_id: i32,
    word: i32,
    window: i32,
    zero_words: bool,
) -> (HardStarts, Option<(u32, std::path::PathBuf, HardMaskShape)>) {
    use crate::algorithm::common::translated_chunks::{frame_length, translated_db_mask_chunks};
    use crate::core::subject_split::complement_masks;
    const NCBI_FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

    let scan = HardScan::Translated {
        word,
        window,
        zero_words,
    };
    let min_length = translated_min_length(word);
    let carried = window < 0;
    let mut walk = NcbiWalk::new(word, window);
    let mut starts = HardStarts::default();
    let mut ranges = Vec::new();
    let mut split = HardSplitState::default();
    let mut items = Vec::new();

    let found = db
        .searched_sequences(algo_id, |ordinal, length, masks| {
            let masks = masks.read()?;
            if let Some(shape) =
                HardMaskShape::of_scan(masks, Some((length, length)), &|| Some(length), scan, None)
                    .ok()
                    .flatten()
            {
                return Ok(Some(shape));
            }
            if length < min_length {
                return Ok(None);
            }
            let walked = translated_walked(masks);
            complement_masks(masks, length, &mut ranges);
            let neg_first = frame_length(length, -1).unwrap_or(0);
            items.clear();
            let mut no_chunking: Option<(i8, u32, u32)> = None;
            for frame in NCBI_FRAMES {
                let Some(fl) = frame_length(length, frame) else {
                    continue;
                };
                let Ok(frame_chunks) =
                    translated_db_mask_chunks(&ranges, true, frame, fl, neg_first, &mut split)
                else {
                    return Ok(None);
                };
                for (chunk_idx, (chunk, _)) in frame_chunks.iter().enumerate() {
                    if walked && no_chunking.is_none() && chunk.full_range {
                        let mask_idx = chunk_idx.min(masks.len() - 1);
                        let (start, end) = if frame > 0 {
                            masks[mask_idx]
                        } else {
                            masks[masks.len() - 1 - mask_idx]
                        };
                        no_chunking = Some((frame, start, end));
                    }
                    items.push((fl, chunk.offset, chunk.ncbi_length));
                }
            }
            let step = walk.translated_masked(&items, zero_words);
            if let Some(offset) = step.word_past_end {
                return Ok(Some(HardMaskShape::WordPastEnd { offset }));
            }
            if let Some((frame, start, end)) = no_chunking {
                return Ok(Some(HardMaskShape::NoChunkingWithOverlaps {
                    frame,
                    start,
                    end,
                }));
            }
            if !carried {
                if let Some(start) = step.ncbi {
                    starts.push(ordinal, start);
                } else if walked {
                    let offset =
                        i32::try_from(i64::from(window) + step.below_start).unwrap_or(i32::MAX);
                    starts.push(
                        ordinal,
                        NcbiStart {
                            offset,
                            cleared: false,
                        },
                    );
                }
            }
            Ok(None)
        })
        .ok()
        .flatten();
    (starts, found)
}

#[cfg(test)]
#[path = "db_mask_tests.rs"]
mod db_mask_tests;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:447-448, 478-479
// The chunker is told whether the program's subjects are nucleotide; for a protein subject it
// takes no residual, so each chunk starts at its hard range's left.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-500 and 762-766 (cited
// at `blastp::blast_engine::x_soft_subject::x_for_each_hit`): the two-hit scan of a protein
// subject raises the end of its first range to its start, so it reads one word of a chunk shorter
// than a word, running past the chunk; the one-hit scan does not, and reads no word there.
/// The scan of a hard-masked search, which decides the chunks of its sequences and which of
/// their shapes LOSAT refuses (`HardMaskShape::of_scan`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardScan {
    /// BLASTN: chunks from a multiple of `COMPRESSION_RATIO` (NCBI's `is_nucleotide`), whose scan
    /// reads words of `scanned_word` letters (the word size, or the discontiguous template's
    /// length).
    Nucleotide {
        /// The letters of a scanned word.
        scanned_word: usize,
    },
    /// blastp and blastx (protein subjects, series NP, M4b): chunks from their hard range's left
    /// (`residual` 0), in a search whose words have `word` letters and whose two-hit window is
    /// `window` (0: the one-hit search, `-window_size 0`). A mask that runs to the end after the
    /// chunk from the start (after empty masks at 0 too) is searched over that chunk alone, at
    /// every word (`HardMaskShape::TailMask` is never refused here: the same NCBI defect as
    /// BLASTN's [`TAIL_MASK_SEARCHED_WORD`] (PD-LOSAT-NCBI-DEFECTS exception 13), approved for
    /// these programs at every word as exception 21; the Owner's answer Q1 = B of M4b,
    /// 2026-11-21). Masks may overlap or nest, and be out of the order of their start (series NP,
    /// M4b-3b: the Owner's answer of 2026-11-21 on overlapping masks): every chunk is searched as
    /// NCBI searches it, those between overlapping masks of a negative length too, and the
    /// diagonal offset moves by each chunk's length plus the window
    /// (`algorithm::common::hard_chunk_walk::DiagonalWalk`); the one-hit scan reads no word of a
    /// chunk shorter than the word, the two-hit scan the one word at its start. Refused: a mask
    /// that ends before it starts (`HardMaskShape::Inverted`), a part between masks that can hold
    /// alignments and starts before the end of an earlier such part, whose alignments NCBI's
    /// merge joins (`HardMaskShape::PartsOutOfOrder`, not yet ported), a hard range longer than
    /// `MAX_DBSEQ_LEN`, which NCBI cuts into chunks that its merge joins across their overlap
    /// (`HardMaskShape::LongRange`, not yet ported), a mask that ends past the sequence
    /// (`HardMaskShape::PastEnd`), a mask that runs to the end after the chunk from the start with
    /// masks stored after it (`HardMaskShape::TailMaskWithLaterMasks`), and the one word of a
    /// chunk that runs past the sequence's sentinel where it can change a result
    /// (`HardMaskShape::WordPastEnd`). Chunks whose word hits would reach the diagonal entries of
    /// other sequences are searched with those entries out of reach (PD-LOSAT-NCBI-DEFECTS
    /// exception 29, the Owner's answer Q5 = B of M4b, 2026-11-22: `algorithm::common::
    /// hard_chunk_walk::DiagonalWalk::start`), and those whose hits meet entries that no sequence
    /// wrote in NCBI's one-thread search start on such a table at NCBI's offset (`algorithm::
    /// common::hard_chunk_walk::NcbiWalk`).
    Protein {
        /// The letters of a word.
        word: i32,
        /// The two-hit window (0: the one-hit search).
        window: i32,
        /// Whether the query's lookup table can hold a word with the residue 0, the sequence's
        /// sentinel byte (BLASTP: a query that holds the gap `-`; BLASTX never, its translation
        /// gives `X` for a codon with a gap): `common::hard_chunk_walk::DiagonalWalk::new`.
        zero_words: bool,
    },
    /// tblastn and tblastx (translated nucleotide subjects; series NP, M4b PR 2 and M4b-3c part 2):
    /// the complement of a sequence's stored masks is mapped to each of its six frames (dividing
    /// by 3, the minus frames reflected on frame -1's length; `common::translated_chunks::
    /// translated_frame_ranges`) and each frame is chunked as a protein subject (residual 0). Masks
    /// that overlap or nest in the order of their start are searched as NCBI searches them (every
    /// frame's chunk moves the diagonal table by its `Int4` length plus the window, a chunk shorter
    /// than the word gets its one word in the two-hit scan, and all six frames of a walked or New
    /// sequence run in order on one table); masks that end before they start or start before the
    /// start of an earlier mask stay refused (`HardMaskShape::Unsorted`), as do a mask ending past
    /// the sequence (`HardMaskShape::PastEnd`), a one word past a frame's sentinel when
    /// `zero_words` is true (`HardMaskShape::WordPastEnd`), and a walked sequence with a frame
    /// taking NCBI's no-chunking branch (`HardMaskShape::NoChunkingWithOverlaps`). A sequence whose
    /// hits would reach entries written by earlier sequences is searched on a table of its own
    /// (PD-LOSAT-NCBI-DEFECTS approved exception 29).
    /// For masks apart, `TailMask` read per frame (exception 21) and an empty first part before the
    /// no-chunking branch (exception 27) are searched as before. Both engines join the chunks of a
    /// frame longer than `MAX_DBSEQ_LEN` across their overlap, as for an unmasked subject (no
    /// `LongRange`). A sequence shorter than three times `word` plus two bases is not searched
    /// (NCBI's subject loop skips it), so it is never refused.
    Translated {
        /// The letters of a word.
        word: i32,
        /// The two-hit window (0: the one-hit search).
        window: i32,
        /// Whether the query's lookup table can hold a word with the residue 0, the frame's
        /// sentinel byte (TBLASTN: a query with `-` or a PSSM whose column 0 scores above
        /// `BLAST_SCORE_MIN`; TBLASTX never, its translation has no residue 0).
        zero_words: bool,
    },
}

impl HardScan {
    /// NCBI's `is_nucleotide` of the chunks (`HardSplitState::start`): a translated subject's
    /// frames are chunked as protein subjects.
    fn nucleotide(self) -> bool {
        matches!(self, HardScan::Nucleotide { .. })
    }
}

/// A searched sequence's masks that `-db_hard_mask` refuses (`ResolvedMask::check_losat_limits`;
/// series NP, M4a-2b, ncbi-parity-20261019/M/m4a2/SURVEY.md §7.1-7.3): NCBI's hard chunks
/// (`s_GetNextSubjectChunk`, blast_engine.c:221-310) give no result that LOSAT can check for them.
/// The shapes are read off the chunks themselves (`core::subject_split::HardSplitState`, the
/// chunker that the search runs; review S1 of M4a-2b), except `Unsorted` and `Inverted`, checks of
/// the stored masks. The first shape is refused only in BLASTN at a scanned word below
/// [`TAIL_MASK_SEARCHED_WORD`] (blastn-short): at that word and above, and at every word in the
/// protein and translated programs, LOSAT searches it over the chunk before the mask, an approved
/// NCBI defect (series NP, M4a-2c: the Owner's answer of 2026-11-15 to question Q1; M4b-3: the
/// Owner's answers of 2026-11-21 to questions Q1 and Q2 of M4b).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardMaskShape {
    /// The mask `[start, end)` runs to the subject's end after a chunk from the subject's start
    /// (only empty masks at 0 before it), which takes NCBI's "no chunking" branch and so scans
    /// the full range past the chunk, in a subject that the scan reads a word of, at a scanned
    /// word below [`TAIL_MASK_SEARCHED_WORD`] (BLASTN only).
    TailMask { start: u32, end: u32 },
    /// A mask that ends before it starts, or that starts before the end of the mask before it in
    /// BLASTN (`HardScan::Nucleotide`), or before the start of the mask before it in TBLASTN and
    /// TBLASTX (`HardScan::Translated`: masks that overlap or nest in the order of their start are
    /// searched, while masks out of start order fail NCBI's partial-fetch slice scan and merge;
    /// the protein programs, `HardScan::Protein`, read their masks' order off the chunks:
    /// `HardMaskShape::Inverted`, `HardMaskShape::PartsOutOfOrder`).
    Unsorted {
        before: (u32, u32),
        mask: (u32, u32),
    },
    /// The last mask `[start, end)`, after which the subject's final chunk (one chunk for the
    /// last hard range) starts before the end of an earlier chunk of another hard range, both
    /// chunks scanning a word (so both can hold HSPs): the mask is shorter than the `end % 4`
    /// letters that the final chunk starts with, and the overlap that NCBI reads past its list of
    /// ranges decides which of the final chunk's HSPs merge with the earlier chunk's across the
    /// mask.
    ShortLastMask { start: u32, end: u32 },
    /// A hard range `[left, right)` between masks (or a sequence end) longer than
    /// `MAX_DBSEQ_LEN`, in a protein program (`HardScan::Protein`; series NP, M4b): NCBI cuts it
    /// into chunks that overlap and joins their HSPs across the overlap, which LOSAT's protein
    /// programs do not port.
    LongRange { left: i32, right: i32 },
    /// A mask `[start, end)` that ends past the sequence's `length` letters, in a protein program
    /// (`HardScan::Protein`; series NP, M4b; review L1 of PR 1) or a translated one
    /// (`HardScan::Translated`, PR 2: `length` bases). makeblastdb does not write such
    /// a mask (it refuses one, writedb_impl.cpp:1376-1384), so only a database written by another
    /// tool, or a damaged one, has it; where such a mask starts past the end, NCBI's chunk before
    /// it runs past the sequence's letters (blast_engine.c:254-258).
    PastEnd { start: u32, end: u32, length: usize },
    /// A chunk at `offset` shorter than the word, of a two-hit search in a protein or translated
    /// program (`HardScan::Protein`, `HardScan::Translated`), whose one word runs past the
    /// sequence's or frame's sentinel where its hits can change a result
    /// (`algorithm::common::hard_chunk_walk::DiagonalWalk::word_past_end`).
    WordPastEnd { offset: usize },
    /// A mask `[start, end)` that ends before it starts, in a protein program
    /// (`HardScan::Protein`; series NP, M4b-3b). makeblastdb does not write such a mask (it
    /// refuses one, writedb_impl.cpp:1376-1384), so only a database written by another tool, or a
    /// damaged one, has it.
    Inverted { start: u32, end: u32 },
    /// A chunk that can hold alignments (of `word` letters or more) and starts before the end of
    /// an earlier such chunk, in a protein program (`HardScan::Protein`; series NP, M4b-3b):
    /// `mask` is the mask before the later chunk and `before` the mask after the earlier one, so
    /// that `mask` ends before `before` starts though it is stored after it. NCBI's merge after the
    /// later chunk takes the earlier chunks' HSPs that end past its offset to join them with its
    /// own (`Blast_HSPListsMerge`, blast_hits.c:2888-2911), which LOSAT's protein programs do not
    /// do yet (with `LongRange`'s merge, item M4b-3d of series NP). Masks out of the order of
    /// their start whose chunks do not meet so are searched (the chunks between them shorter
    /// than the word, or apart: segmasker writes such masks); where the HSPs of a chunk shorter
    /// than the word meet another chunk's so, the search refuses the sequence when it finds them
    /// (`common::hard_chunk_scan::short_chunk_merge_refusal`).
    PartsOutOfOrder {
        before: (u32, u32),
        mask: (u32, u32),
    },
    /// A mask `[start, end)` that runs to the end after the chunk from the start, with masks
    /// stored after it, in a protein program (`HardScan::Protein`; series NP, M4b-3b, review L1):
    /// NCBI's chunk loop ends at that chunk (its "no chunking" branch scans the full range past
    /// it, PD-LOSAT-NCBI-DEFECTS exception 21) and never applies the later masks in its search,
    /// while its alignments read their letters as masked (probe of M/m4b3b/probe/fix1: rows that
    /// differ from those of the masks without the later one, and searches that never end in the
    /// one-hit search). Exception 21 defines LOSAT's result of a tail mask by the input with that
    /// mask one letter shorter, whose chunks differ here, so no near input defines it.
    TailMaskWithLaterMasks { start: u32, end: u32 },
    /// A walked translated sequence (`HardScan::Translated` with overlapping masks giving a
    /// negative frame chunk) one of whose reading frames has a chunk taking NCBI's no-chunking
    /// branch (`HardChunk::full_range`), which scans the whole frame past the later masks
    /// (exceptions 21 and 27 define a result only for masks apart).
    NoChunkingWithOverlaps { frame: i8, start: u32, end: u32 },
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-268 (cited at
// `core::subject_split::HardSplitState::next_chunk`)
/// The shortest scanned word (the word size, or the discontiguous template's length) at which
/// `-db_hard_mask` searches a sequence of the shape [`HardMaskShape::TailMask`] instead of refusing
/// it: over the chunk before the mask alone, scanned from its start as NCBI scans its other hard
/// chunks (`core::subject_split::HardSplitState::next_chunk`). An approved NCBI defect
/// (PD-LOSAT-NCBI-DEFECTS exception 13, version 1.35; the Owner's decision of 2026-11-15, 「M4a
/// Q1=B」; series NP, M4a-2c): NCBI BLAST+ 2.17.0 takes the "no chunking" branch for that chunk,
/// so it scans the subject's stored letters over the full range (its soft ranges,
/// blast_engine.c:162-164), the masked ones included, while its extensions and its interval tree
/// take the chunk's length (blast_gapalign.c:3830-3831). In the 900 runs of the Owner's question
/// at words of 11 and more (megablast, blastn, dc-megablast) its output was the same as for the
/// input with the mask one letter shorter, which its normal chunks search, while at blastn-short's
/// 7 it differed in 29 of 300 runs and crashed in 71, so the refusal stays below this word. At 11
/// and more NCBI can still add a row from a word across the mask's start, lose the hits of the
/// subject it searches next on the diagonals of the masked words (its diagonal table is left
/// beyond the offset that `Blast_ExtendWordExit` moves on by the chunk's length; with threads
/// which subject comes next varies by run), or run out of memory in its interval tree and crash.
/// LOSAT's traceback reads the whole mask as `B`, as NCBI's does; everywhere else LOSAT gives
/// what NCBI's normal chunks give for the near input (the exception's text). The protein and
/// translated programs search the shape at every word (`HardScan::Protein`, `HardScan::Translated`;
/// PD-LOSAT-NCBI-DEFECTS exception 21, the Owner's answers Q1 = B and Q2 = B of M4b, 2026-11-21:
/// NCBI BLAST+ 2.17.0 crashes on some of these sequences and gives rows that differ from the near
/// input's on others, in tblastx the same row many times, M/m4b/SURVEY.md §4).
pub const TAIL_MASK_SEARCHED_WORD: usize = 11;

impl HardMaskShape {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276 and 579-584
    // (cited at `core::subject_split::HardSplitState::next_chunk`): the chunks, the "no
    // chunking" branch, and the overlap after each chunk, read one range past the list after the
    // last one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1668
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354, 1368-1370
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500, 565
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2873-2911, 2956
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1527-1537
    // (cited at blastn's `database::mask_hard_letters`): the partial fetch's slices are masked by
    // a scan that assumes sorted masks.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:233-242
    // A protein database is never fetched in part.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1326-1333
    /// The shape of the stored masks `masks` (`[start, end)`, as stored, not empty) of a sequence
    /// whose length `length` gives, that `-db_hard_mask` refuses in a search whose scan is `scan`
    /// (in the protein programs where the search starts it: on a new table at `ncbi`, else
    /// after sequences whose entries do not reach it, `algorithm::common::hard_chunk_walk::
    /// NcbiWalk`)
    /// (BLASTN's `HardScan::Nucleotide`, whose words have `scanned_word` letters: the word size,
    /// or the discontiguous template's length), or
    /// `None`; `Err(())` when the length cannot be read (`SeqDb::first_searched_masks`). With
    /// `bounds`, the lowest and highest length that the index gives (`Volume::
    /// nucleotide_length_bounds`; the sequence's length is one of the four), `length` is read
    /// only where the masks end at or after the lowest (a mask that may run to the end) or the
    /// shapes at the two ends differ (review 2 S4 of M4a-2b). With every mask ending before
    /// the lowest, the chunks of the ranges before the last are the same at every length, no
    /// chunk shorter than the sequence takes the "no chunking" branch, and the last range gives
    /// `ShortLastMask` only while it is one chunk (fewer letters; after a full chunk of its own
    /// range its final chunk starts past every other range's chunks) and that chunk scans a
    /// word (more letters): the two conditions change within the bounds' three letters only for
    /// a word of more than `MAX_DBSEQ_LEN - 6` letters, and BLASTN's words have at most
    /// `DBSEQ_CHUNK_OVERLAP` (blast_options.c:1326-1333; templates 16 to 21), so the shape at
    /// both ends is the shape at every length between them:
    /// - masks out of order (a check of the masks, before the chunks): hard ranges that go back
    ///   (negative chunk lengths, which move the diagonal offset back) and slices that NCBI's
    ///   partial fetch masks by a scan of sorted masks (§7.2); in the protein programs, which
    ///   NCBI never fetches in part (seqsrc_seqdb.cpp:233-242), only a mask that ends before it
    ///   starts (`Inverted`), the order of the others read off the chunks (`PartsOutOfOrder`);
    /// - else the sequence's chunks (`HardSplitState` on the complement of the masks, the
    ///   search's own chunker, with `DBSEQ_CHUNK_OVERLAP`): a chunk that takes the "no chunking"
    ///   branch while it is shorter than the sequence (so a mask with letters follows it) scans
    ///   the full range `[0, length)` past the chunk, which reads a word when `length >=
    ///   scanned_word` (`min_subj_seq_length` is 1 for BLASTN, so every subject is scanned):
    ///   NCBI reads seeds past the chunk (§7.1: rows of the input with the mask one letter
    ///   shorter, rows of its own, or a crash; PD-LOSAT-NCBI-DEFECTS exception 13), refused at a
    ///   `scanned_word` below [`TAIL_MASK_SEARCHED_WORD`] (at that word and above LOSAT scans the
    ///   chunk alone, the approved NCBI defect; the chunk is the sequence's last, so the next rule
    ///   never reads it);
    /// - the final chunk, after which NCBI reads the overlap one range past the list (another
    ///   sequence's mask end, or unset memory, §7.3), where that overlap can decide a merge: the
    ///   final chunk scans a word (else its list is empty and NCBI goes on to the next chunk
    ///   before it reads the overlap, blast_engine.c:500, 565), and an earlier chunk of another
    ///   hard range that scans a word ends after the final chunk's offset (its HSPs can then end
    ///   there, the merge's `hspcnt1`; a chunk that scans no word holds none, review 2 F1). A
    ///   chunk's scan reads a word when its range `[residual, length)` (the full range in the
    ///   "no chunking" branch) holds `scanned_word` letters (`s_DetermineScanningOffsets`). Where
    ///   only full chunks of the final chunk's own long range (a last range over
    ///   `MAX_DBSEQ_LEN` letters) end after its offset, the merge is within the range and
    ///   LOSAT's overlap (`dbseq_chunk_overlap`) is kept, as the unmasked search keeps it
    ///   (decision B-d).
    ///
    /// The protein programs' scan (`HardScan::Protein`, series NP, M4b) has chunks that start at
    /// their range's left, so that `ShortLastMask` (a residual before the final chunk's range)
    /// never occurs, a mask that runs to the end after the chunk from the start is searched over
    /// that chunk (series NP, M4b-3), a mask that ends past the sequence is `PastEnd` (before its
    /// chunks: NCBI's chunk before a mask that starts past the end reads past the sequence), a
    /// hard range longer than `MAX_DBSEQ_LEN` is `LongRange`, a chunk that can hold alignments
    /// starting before the end of an earlier one is `PartsOutOfOrder`, and the walk of the chunks'
    /// diagonal offsets gives `WordPastEnd` (series NP, M4b-3b, `HardMaskShape::of_protein_chunks`). Its shape is read at the sequence's length only, which
    /// a protein database's index gives exactly (`SeqDb::first_searched_masks`): whether the one
    /// word of a short chunk near the end covers the sentinel changes with the length.
    pub fn of_scan(
        masks: &[(u32, u32)],
        bounds: Option<(usize, usize)>,
        length: &dyn Fn() -> Option<usize>,
        scan: HardScan,
        ncbi: Option<NcbiStart>,
    ) -> Result<Option<HardMaskShape>, ()> {
        if masks.is_empty() {
            return Ok(None);
        }
        // The protein programs search masks in any order that end where they start or after, the
        // order read off their chunks (`of_protein_chunks`); TBLASTN and TBLASTX accept masks in
        // the order of their start (overlapping or nested included) and refuse masks out of start
        // order or inverted; BLASTN needs the masks apart and sorted.
        if matches!(scan, HardScan::Protein { .. }) {
            if let Some(&(start, end)) = masks.iter().find(|&&(start, end)| end < start) {
                return Ok(Some(HardMaskShape::Inverted { start, end }));
            }
        } else {
            let start_order_only = matches!(scan, HardScan::Translated { .. });
            let mut before: Option<(u32, u32)> = None;
            for &mask in masks {
                let out_of_order = before.is_some_and(|(start, end)| {
                    mask.0 < if start_order_only { start } else { end }
                });
                if mask.1 < mask.0 || out_of_order {
                    return Ok(Some(HardMaskShape::Unsorted {
                        before: before.unwrap_or(mask),
                        mask,
                    }));
                }
                before = Some(mask);
            }
        }
        if let HardScan::Translated { word, .. } = scan {
            // Only a mask that ends past the sequence is refused here (`HardMaskShape::of_frames`;
            // overlapping masks may nest, so the last mask need not have the largest end): where
            // every mask ends at or before the shortest length within the bounds none does, and
            // the length is not read.
            let max_end = masks
                .iter()
                .map(|&(_, end)| end as usize)
                .max()
                .unwrap_or(0);
            if bounds.is_some_and(|(low, _)| max_end <= low) {
                return Ok(None);
            }
            let length = length().ok_or(())?;
            return Ok(Self::of_frames(masks, length, translated_min_length(word)));
        }
        let mut split = HardSplitState::default();
        split.ranges.reserve(masks.len() + 1);
        if matches!(scan, HardScan::Protein { .. }) {
            // The exact length, from the index or read.
            let length = match bounds {
                Some((low, high)) if low == high => low,
                _ => length().ok_or(())?,
            };
            return Ok(Self::of_chunks(masks, length, scan, ncbi, &mut split));
        }
        let last_end = masks[masks.len() - 1].1 as usize;
        if let Some((low, high)) = bounds.filter(|&(low, _)| last_end < low) {
            let shape = Self::of_chunks(masks, low, scan, ncbi, &mut split);
            if Self::of_chunks(masks, high, scan, ncbi, &mut split) == shape {
                return Ok(shape);
            }
        }
        let length = length().ok_or(())?;
        Ok(Self::of_chunks(masks, length, scan, ncbi, &mut split))
    }

    /// The shape of the masks `masks` of a sequence of `length` letters (sorted and apart except in
    /// a protein program, where the search starts it on a new table at `ncbi`, or not), read off
    /// its chunks (`HardMaskShape::of_scan`); `split` is the chunker, whose buffer of ranges is
    /// reused.
    fn of_chunks(
        masks: &[(u32, u32)],
        length: usize,
        scan: HardScan,
        ncbi: Option<NcbiStart>,
        split: &mut HardSplitState,
    ) -> Option<HardMaskShape> {
        use crate::algorithm::common::hard_chunk_walk::DiagonalWalk;
        use crate::core::subject_split::{
            complement_masks, HardChunk, DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN,
        };
        let protein = matches!(scan, HardScan::Protein { .. });
        if protein {
            if let Some(&(start, end)) = masks.iter().find(|&&(_, end)| end as usize > length) {
                return Some(HardMaskShape::PastEnd { start, end, length });
            }
        }
        let mut ranges = std::mem::take(&mut split.ranges);
        complement_masks(masks, length, &mut ranges);
        if protein {
            // The masks end within the sequence; a range between masks that overlap is inverted
            // (no letters, never long).
            if let Some(&(left, right)) = ranges
                .iter()
                .find(|&&(left, right)| right - left > MAX_DBSEQ_LEN as i32)
            {
                split.ranges = ranges;
                return Some(HardMaskShape::LongRange { left, right });
            }
        }
        split.start(length, ranges, scan.nucleotide());
        let (scanned_word, tail_refused) = match scan {
            HardScan::Nucleotide { scanned_word } => {
                (scanned_word, scanned_word < TAIL_MASK_SEARCHED_WORD)
            }
            HardScan::Protein {
                word,
                window,
                zero_words,
            } => {
                let walk = DiagonalWalk::new(word, window, length, zero_words);
                return Self::of_protein_chunks(masks, length, walk, ncbi, split);
            }
            HardScan::Translated { .. } => (0, false),
        };
        let scans_word = |chunk: &HardChunk| {
            if chunk.full_range {
                length >= scanned_word
            } else {
                chunk.letters().saturating_sub(chunk.residual) >= scanned_word
            }
        };
        // The farthest end of the chunks that scan a word: of the ranges before the chunk's
        // range, and of the chunk's own range so far.
        let (mut before_range, mut in_range, mut range) = (0usize, 0usize, 0usize);
        while let Some(chunk) = split.next_chunk(DBSEQ_CHUNK_OVERLAP) {
            if chunk.range != range {
                (before_range, in_range, range) = (before_range.max(in_range), 0, chunk.range);
            }
            let scans = scans_word(&chunk);
            // The chunk's range ends where the mask of the same index starts.
            if chunk.full_range && chunk.letters() < length && scans && tail_refused {
                let (start, end) = masks[chunk.range];
                return Some(HardMaskShape::TailMask { start, end });
            }
            if chunk.overlap_past_list && scans && before_range > chunk.offset {
                let (start, end) = masks[chunk.range - 1];
                return Some(HardMaskShape::ShortLastMask { start, end });
            }
            if scans {
                in_range = in_range.max(chunk.offset + chunk.letters());
            }
        }
        None
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276 (cited at
    // `core::subject_split::HardSplitState::next_chunk`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107, 110-155,
    // 162-173 (cited at `algorithm::common::hard_chunk_walk::DiagonalWalk`): a thread's table is
    // made with its entries 0 and its offset the window, a clear sets every entry to the window
    // less than 0, and the exit of each chunk's word finder moves the offset.
    /// The shape of a protein program's chunks of a sequence of `length` letters (`split`, started
    /// on the hard ranges of its masks `masks`), walked with `walk` (a new walk of the search's
    /// word and window) from where the search starts it (`ncbi`: on a new table at NCBI's
    /// one-thread offset, `NcbiWalk`): refused where the chunk from the start before a mask that
    /// runs to the end has masks stored after that mask (`TailMaskWithLaterMasks`), where a chunk
    /// that can hold alignments starts before the end of an earlier one (`PartsOutOfOrder`, the
    /// first such chunk), else by the walk of their diagonal offsets that the search runs
    /// (`DiagonalWalk`) where the one word of a short chunk runs past its sentinel with an effect
    /// (`WordPastEnd`; on a new table the walk starts after its entries, which a chunk below the
    /// start can meet, `DiagonalWalk::start_on_ncbi_table`; other sequences' entries never reach
    /// the hits of the others, `DiagonalWalk::start_after_others`). The chunks start at their
    /// range's left, so that `ShortLastMask` never occurs.
    fn of_protein_chunks(
        masks: &[(u32, u32)],
        length: usize,
        mut walk: crate::algorithm::common::hard_chunk_walk::DiagonalWalk,
        ncbi: Option<NcbiStart>,
        split: &mut HardSplitState,
    ) -> Option<HardMaskShape> {
        use crate::core::subject_split::DBSEQ_CHUNK_OVERLAP;
        let word = walk.word();
        if let Some(start) = ncbi {
            walk.start_on_ncbi_table(start);
        }
        // The farthest end of the earlier chunks that can hold alignments (the first of them on a
        // tie) and its range. A chunk's HSPs end within its letters; those of a chunk shorter
        // than the word come from its one word and the search checks them as it finds them
        // (`common::hard_chunk_scan::short_chunk_merge_refusal`). NCBI's merge after a chunk with
        // HSPs
        // takes those of the earlier chunks that end past its offset.
        let mut held: Option<(usize, usize)> = None;
        while let Some(chunk) = split.next_chunk(DBSEQ_CHUNK_OVERLAP) {
            // The chunk's range ends where the mask of the same index starts: a mask that runs to
            // the end, after which the loop ends.
            if chunk.full_range && chunk.letters() < length && chunk.range + 1 < masks.len() {
                let (start, end) = masks[chunk.range];
                return Some(HardMaskShape::TailMaskWithLaterMasks { start, end });
            }
            if chunk.length >= word {
                let end = chunk.offset + chunk.letters();
                match held {
                    Some((held_end, range)) if held_end > chunk.offset => {
                        // The earlier chunk's range ends where the mask of its index starts, and
                        // this chunk's range starts where the mask before its index ends.
                        return Some(HardMaskShape::PartsOutOfOrder {
                            before: masks[range],
                            mask: masks[chunk.range - 1],
                        });
                    }
                    Some((held_end, _)) if held_end >= end => {}
                    _ => held = Some((end, chunk.range)),
                }
            }
            walk.chunk(chunk.offset, chunk.length);
        }
        walk.word_past_end()
            .map(|offset| HardMaskShape::WordPastEnd { offset })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:254-258 (cited at
    // `HardMaskShape::PastEnd`)
    /// The refused shape of the sorted masks `masks` of a translated subject of `length` bases
    /// (`HardScan::Translated`): a mask that ends past the sequence (`PastEnd`), in a sequence that
    /// the search reads (`min_length` bases or more). Every other shape of sorted masks is
    /// searched as its frames' chunks give it (`common::translated_chunks::
    /// translated_db_mask_chunks`): a frame's first part that NCBI's "no chunking" branch scans
    /// past (exception 21, series NP, M4b-3) and one after an empty first part, on which NCBI
    /// crashes (exception 27, M4b-3c).
    fn of_frames(masks: &[(u32, u32)], length: usize, min_length: usize) -> Option<HardMaskShape> {
        if length < min_length {
            return None;
        }
        // Masks are in start order but may overlap or nest, so scan all ends for one past the
        // sequence. Past the sequence, its nucleotide range would be inverted.
        masks
            .iter()
            .find(|&&(_, end)| end as usize > length)
            .map(|&(start, end)| HardMaskShape::PastEnd { start, end, length })
    }

    /// The refusal's words for the shape.
    pub fn text(&self) -> String {
        match *self {
            HardMaskShape::TailMask { start, end } => format!(
                "has a mask, [{start}, {end}), running to its end after the chunk from its \
                 start, for which NCBI BLAST+ 2.17.0 scans the masked letters past that chunk \
                 while it extends within the chunk (at a word of fewer than \
                 {TAIL_MASK_SEARCHED_WORD} letters: rows that differ from a near input's, or a \
                 crash)"
            ),
            HardMaskShape::Unsorted { before, mask } => format!(
                "has masks out of order ([{}, {}) after [{}, {})), whose hard masking NCBI BLAST+ \
                 2.17.0 ties to the slices of its partial fetch",
                mask.0, mask.1, before.0, before.1
            ),
            HardMaskShape::ShortLastMask { start, end } => format!(
                "has a last mask, [{start}, {end}), shorter than the {} letters before its end \
                 that the last chunk starts with, after which NCBI BLAST+ 2.17.0 merges that \
                 chunk's alignments with an earlier chunk's by an overlap read past its list of \
                 ranges",
                end % 4
            ),
            HardMaskShape::LongRange { left, right } => format!(
                "has a range between its masks, [{left}, {right}), longer than the 5,000,000 \
                 letters of NCBI BLAST+ 2.17.0's chunks, whose overlapping chunks NCBI joins by a \
                 merge across their overlap"
            ),
            HardMaskShape::PastEnd { start, end, length } => format!(
                "has a mask, [{start}, {end}), ending past its {length} letters, which \
                 makeblastdb does not write (NCBI BLAST+ 2.17.0's chunk before a mask that \
                 starts past the end reads past the sequence)"
            ),
            HardMaskShape::WordPastEnd { offset } => format!(
                "has a part between masks at {offset} shorter than the word, whose one word NCBI \
                 BLAST+ 2.17.0 reads past the sequence's end, from the next sequence's letters, \
                 where its hits can change the rows"
            ),
            HardMaskShape::Inverted { start, end } => format!(
                "has a mask, [{start}, {end}), that ends before it starts, which makeblastdb does \
                 not write"
            ),
            HardMaskShape::PartsOutOfOrder { before, mask } => format!(
                "has masks out of order ([{}, {}) after [{}, {})): its part after [{}, {}) starts \
                 before the end of its part before [{}, {}), and NCBI BLAST+ 2.17.0 joins the two \
                 parts' alignments by a merge across that overlap, which LOSAT does not do yet",
                mask.0, mask.1, before.0, before.1, mask.0, mask.1, before.0, before.1
            ),
            HardMaskShape::TailMaskWithLaterMasks { start, end } => format!(
                "has a mask, [{start}, {end}), running to its end after the part from its start, \
                 and masks stored after it, which NCBI BLAST+ 2.17.0 does not apply when it \
                 searches that part past its end but applies to the letters of its alignments"
            ),
            HardMaskShape::NoChunkingWithOverlaps { frame, start, end } => format!(
                "has masks that overlap and reading frame {frame}'s part from its start running \
                 to the frame's end at the mask [{start}, {end}), which NCBI BLAST+ 2.17.0 \
                 searches as the whole frame past the later masks"
            ),
        }
    }
}
