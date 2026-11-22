//! `-db_soft_mask` and `-db_hard_mask`: the masks of a BLAST database's sequences that a search
//! applies (NCBI's `CBlastDatabaseArgs` with `SetDatabaseMaskingSupport(true)` and the filtering
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

use crate::cli::NativeError;
use crate::core::subject_split::HardSplitState;
use crate::report::archive::fields;
use crate::report::archive::options::{Param, ParamValue};
use crate::seqdb::{MaskAlgorithm, NameMaskError, SeqDb, SeqDbError};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2342-2356
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2226-2229,2367-2382
// The exclusions are errors of the argument parser: LOSAT's parser reports them with its own
// text (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1).
/// The two options of a program whose database arguments support masking (BLASTN).
#[derive(clap::Args, Debug, Clone, Default)]
pub struct DbMaskArgs {
    /// `-db_soft_mask`.
    #[arg(
        long = "db_soft_mask",
        value_name = "filtering_algorithm",
        conflicts_with_all = ["subject", "subject_loc", "db_hard_mask"],
        help = "Filtering algorithm ID to apply to the BLAST database as soft masking"
    )]
    pub db_soft_mask: Option<String>,
    /// `-db_hard_mask` (refused for the masks of `HardMaskShape`: `ResolvedMask::check_losat_limits`).
    #[arg(
        long = "db_hard_mask",
        value_name = "filtering_algorithm",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "Filtering algorithm ID to apply to the BLAST database as hard masking"
    )]
    pub db_hard_mask: Option<String>,
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

/// What the option's text names (`CSearchDatabase::SetFilteringAlgorithm(string, type)`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FilteringAlgorithm {
    /// `m_FilteringAlgorithmId`: a text that is an `Int4` in NCBI's form (`s_IsNumericId`).
    Id(i32),
    /// `m_FilteringAlgorithmString`: any other text, looked up by the database
    /// (`GetMaskAlgorithmId`, case-sensitive).
    Name(String),
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
        Some(DbMaskOption {
            kind,
            algorithm: FilteringAlgorithm::parse(text),
        })
    }
}

impl FilteringAlgorithm {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:149-171
    // `x_ValidateMaskingAlgorithm` returns at once here: the database is not open yet.
    /// The algorithm that the option's text names: `+11` and `011` are the id 11; ` 11`, `11 `,
    /// `11x`, `2147483648` and the empty text are names.
    pub fn parse(text: &str) -> Self {
        // `s_IsNumericId` is `StringToInt` that converts (`seqdb::mask::is_numeric_id`).
        match crate::ncbi_str::string_to_int(text.as_bytes(), 0) {
            Ok(id) => FilteringAlgorithm::Id(id),
            Err(_) => FilteringAlgorithm::Name(text.to_string()),
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
    pub key: Option<String>,
    /// The algorithm in the database of each name of the `-db` list, in the order of the list
    /// (`SDbInfo::filt_algorithm_name`, `filt_algorithm_options`).
    pub shown: Vec<MaskAlgorithm>,
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
        db_name: &str,
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
            return Err(NativeError {
                exit: 1,
                message: format!(
                    "{OPTIONS_ERROR} Masking algorithm ID {id} is not supported in {} '{db_name}' BLAST database\n",
                    if nucleotide { "nucleotide" } else { "protein" }
                ).into(),
            }
            .into());
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
    /// lacks).
    pub fn resolve(
        &self,
        db: &SeqDb,
        db_name: &str,
        program: &str,
        warning: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<Option<ResolvedMask>> {
        let (id, key) = match &self.algorithm {
            FilteringAlgorithm::Id(id) => (Some(*id), None),
            FilteringAlgorithm::Name(name) => {
                let id = match db.mask_algorithms() {
                    Ok(algorithms) => algorithms.id_of(name),
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
            warning(
                format!(
                    "Warning: [{}] Subject mask not found in {db_name}, proceeding without subject masking.\n",
                    program.to_ascii_lowercase()
                )
                .as_bytes(),
            )?;
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
    /// (`HardMaskShape::of_scan`).
    pub fn check_losat_limits(&self, db: &SeqDb, program: &str, scan: HardScan) -> Result<()> {
        db.check_search_masks(self.id)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        if self.kind == MaskKind::Hard {
            let found = db.first_searched_masks(self.id, |masks, bounds, length| {
                HardMaskShape::of_scan(masks, bounds, length, scan)
            });
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:499-516
    /// The program options of the BLAST archive and of the search strategy (after the list's).
    /// An empty name is no key (`algo_key != kEmptyStr`): the id it was translated to is written.
    pub fn archive_params(&self) -> Vec<Param> {
        let algorithm = match self.key.as_deref() {
            Some(key) if !key.is_empty() => Param::new(
                fields::DB_FILTERING_ALGORITHM_KEY,
                ParamValue::Str(key.to_string()),
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
    /// (`residual` 0). A mask that runs to the end after the chunk from the start (after empty
    /// masks at 0 too) is searched over that chunk alone, at every word (`HardMaskShape::TailMask`
    /// is never refused here: the same NCBI defect as BLASTN's [`TAIL_MASK_SEARCHED_WORD`]
    /// (PD-LOSAT-NCBI-DEFECTS exception 13), approved for these programs at every word as exception
    /// 18; the Owner's answer Q1 = B of M4b, 2026-11-21). A hard range longer than `MAX_DBSEQ_LEN`,
    /// which NCBI cuts into chunks that its merge joins across their overlap, is refused
    /// (`HardMaskShape::LongRange`), as is a mask that ends past the sequence
    /// (`HardMaskShape::PastEnd`). A chunk shorter than the word is searched as NCBI searches it,
    /// in both scans: the one-hit scan (`-window_size 0`) reads no word of it, and the one word
    /// that the two-hit scan reads at its start has no other hit on its diagonal, so it extends
    /// nothing (`algorithm::common::db_mask_search::SubjectPlan::next_searched_chunk`).
    Protein,
    /// tblastn and tblastx (translated nucleotide subjects; series NP, M4b PR 2): the complement
    /// of a sequence's stored masks is mapped to each of its six frames (`x / 3`, the minus
    /// frames from frame -1's length; `common::translated_chunks::translated_frame_ranges`) and
    /// each frame is chunked as a protein subject (residual 0). The shapes refused are a mask
    /// that ends past the sequence (`HardMaskShape::PastEnd`, as for the protein programs: its
    /// last range would be inverted) and a frame's "no chunking" chunk after an earlier chunk
    /// (`HardMaskShape::FrameStartMask`, on which NCBI crashes). `TailMask` read per frame (the
    /// chunk from the frame's start, shorter than the frame, in NCBI's "no chunking" branch) is
    /// searched over that chunk alone (`common::translated_chunks::translated_db_mask_chunks`;
    /// PD-LOSAT-NCBI-DEFECTS exception 21, the Owner's answers Q1 = B and Q2 = B of M4b,
    /// 2026-11-21). Both engines join the chunks of a frame longer than `MAX_DBSEQ_LEN` across
    /// their overlap, as for an unmasked subject (no `LongRange`), and a chunk shorter than a word
    /// gives no hit in either scan, as NCBI's (the one-hit scan reads no word of it,
    /// aa_ungapped.c:762-766; the two-hit scan reads one, :495-500, which no second hit on its
    /// diagonal follows within the window). A sequence shorter than `min_length` bases is not
    /// searched (NCBI's subject loop skips it), so it is never refused.
    Translated {
        /// The shortest subject searched (`word * 3 + 2` bases).
        min_length: usize,
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
/// chunker that the search runs; review S1 of M4a-2b), except `Unsorted`, a check of the stored
/// masks. The first shape is refused only in BLASTN at a scanned word below
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
    /// A mask that ends before it starts, or that starts before the end of the mask before it.
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
    /// A chunk of frame `frame` of a translated subject (`HardScan::Translated`) that takes NCBI's
    /// "no chunking" branch after an earlier chunk of the frame: after the mapping to the frame
    /// the frame's first hard range has no letters and the range after it starts at 0 too (and
    /// is not cut at `MAX_DBSEQ_LEN`), with no range after it or one at the frame's end.
    /// `[start, end)` is the stored mask between those two ranges in the frame's direction: a
    /// mask within the sequence's first two bases (the plus frames) or within the one or two
    /// bases after its last whole codon from the start (the minus frames, whose mapping reads
    /// `x / 3` from frame -1's length). The earlier chunk allocates NCBI's list of ranges, the
    /// "no chunking" branch points the list at the backup's full range, and NCBI BLAST+ 2.17.0
    /// frees that pointer when it restores the subject (a crash in every run, M/m4b/NOTES_M4B2.md).
    FrameStartMask { frame: i8, start: u32, end: u32 },
    /// A mask `[start, end)` that ends past the sequence's `length` letters, in a protein program
    /// (`HardScan::Protein`; series NP, M4b; review L1 of PR 1) or a translated one
    /// (`HardScan::Translated`, PR 2: `length` bases). makeblastdb does not write such
    /// a mask (it refuses one, writedb_impl.cpp:1376-1384), so only a database written by another
    /// tool, or a damaged one, has it; where such a mask starts past the end, NCBI's chunk before
    /// it runs past the sequence's letters (blast_engine.c:254-258).
    PastEnd { start: u32, end: u32, length: usize },
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1326-1333
    /// The shape of the stored masks `masks` (`[start, end)`, as stored, not empty) of a sequence
    /// whose length `length` gives, that `-db_hard_mask` refuses in a search whose scan is `scan`
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
    ///   partial fetch masks by a scan of sorted masks (§7.2);
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
    /// chunks: its last hard range would be inverted), and a hard range
    /// longer than `MAX_DBSEQ_LEN` is `LongRange`.
    pub fn of_scan(
        masks: &[(u32, u32)],
        bounds: Option<(usize, usize)>,
        length: &dyn Fn() -> Option<usize>,
        scan: HardScan,
    ) -> Result<Option<HardMaskShape>, ()> {
        if masks.is_empty() {
            return Ok(None);
        }
        let mut before: Option<(u32, u32)> = None;
        for &mask in masks {
            if mask.1 < mask.0 || before.is_some_and(|(_, end)| mask.0 < end) {
                return Ok(Some(HardMaskShape::Unsorted {
                    before: before.unwrap_or(mask),
                    mask,
                }));
            }
            before = Some(mask);
        }
        let mut split = HardSplitState::default();
        split.ranges.reserve(masks.len() + 1);
        if let HardScan::Translated { min_length } = scan {
            // A frame's shape changes with the length's remainder by 3 and with masks near
            // either end, so every length within the bounds is read (the sequence's is one of
            // them; a mask that ends past the shortest makes no frame ranges), or the length
            // itself.
            let last_end = masks[masks.len() - 1].1 as usize;
            if let Some((low, high)) = bounds.filter(|&(low, _)| last_end <= low) {
                let shape = Self::of_frames(masks, low, min_length, &mut split);
                if (low + 1..=high)
                    .all(|length| Self::of_frames(masks, length, min_length, &mut split) == shape)
                {
                    return Ok(shape);
                }
            }
            let length = length().ok_or(())?;
            return Ok(Self::of_frames(masks, length, min_length, &mut split));
        }
        // The masks are sorted: the last one ends last.
        let last_end = masks[masks.len() - 1].1 as usize;
        if let Some((low, high)) = bounds.filter(|&(low, _)| last_end < low) {
            let shape = Self::of_chunks(masks, low, scan, &mut split);
            if Self::of_chunks(masks, high, scan, &mut split) == shape {
                return Ok(shape);
            }
        }
        let length = length().ok_or(())?;
        Ok(Self::of_chunks(masks, length, scan, &mut split))
    }

    /// The shape of the sorted masks `masks` of a sequence of `length` letters, read off its
    /// chunks (`HardMaskShape::of_scan`); `split` is the chunker, whose buffer of ranges is reused.
    fn of_chunks(
        masks: &[(u32, u32)],
        length: usize,
        scan: HardScan,
        split: &mut HardSplitState,
    ) -> Option<HardMaskShape> {
        use crate::core::subject_split::{
            complement_masks, HardChunk, DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN,
        };
        if matches!(scan, HardScan::Protein) {
            if let Some(&(start, end)) = masks.iter().find(|&&(_, end)| end as usize > length) {
                return Some(HardMaskShape::PastEnd { start, end, length });
            }
        }
        let mut ranges = std::mem::take(&mut split.ranges);
        complement_masks(masks, length, &mut ranges);
        if matches!(scan, HardScan::Protein) {
            // The masks are sorted and end within the sequence: no range is inverted.
            if let Some(&(left, right)) = ranges
                .iter()
                .find(|&&(left, right)| (right - left) as usize > MAX_DBSEQ_LEN)
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
            HardScan::Protein | HardScan::Translated { .. } => (0, false),
        };
        let scans_word = |chunk: &HardChunk| {
            if chunk.full_range {
                length >= scanned_word
            } else {
                chunk.length.saturating_sub(chunk.residual) >= scanned_word
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
            if chunk.full_range && chunk.length < length && scans && tail_refused {
                let (start, end) = masks[chunk.range];
                return Some(HardMaskShape::TailMask { start, end });
            }
            if chunk.overlap_past_list && scans && before_range > chunk.offset {
                let (start, end) = masks[chunk.range - 1];
                return Some(HardMaskShape::ShortLastMask { start, end });
            }
            if scans {
                in_range = in_range.max(chunk.offset + chunk.length);
            }
        }
        None
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-830
    // The frames' ranges are mapped once per strand from the nucleotide ranges (contexts 0 and
    // 3) and each frame is searched by `s_BlastSearchEngineOneContext` (:422-600), whose chunks
    // are those of a protein subject (`s_GetNextSubjectChunk`, :221-311).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1045-1110
    // (`BLAST_GetAllTranslations`: frame `f` of `L` bases has `(L - |f| + 1) / 3` residues).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:145-215, 264-276
    // (`s_BackupSubject`, `s_AllocateSeqRange`, `s_RestoreSubject` and the "no chunking" branch:
    // `HardMaskShape::FrameStartMask`)
    /// The first refused shape of the sorted masks `masks` of a translated subject of `length`
    /// bases (`HardScan::Translated`), read off the chunks of its six frames in NCBI's order of
    /// contexts: the frame's "no chunking" chunk after an earlier chunk (`FrameStartMask`, NCBI
    /// crashes). As the frame's first chunk, shorter than the frame, the same chunk is searched
    /// over its own letters (`common::translated_chunks::translated_db_mask_chunks`; series NP,
    /// M4b-3).
    fn of_frames(
        masks: &[(u32, u32)],
        length: usize,
        min_length: usize,
        split: &mut HardSplitState,
    ) -> Option<HardMaskShape> {
        use crate::algorithm::common::translated_chunks::translated_frame_ranges;
        use crate::core::subject_split::{complement_masks, DBSEQ_CHUNK_OVERLAP};
        if length < min_length {
            return None;
        }
        // The masks are sorted: the last one ends last. Past the sequence, its nucleotide range
        // would be inverted, and so would a frame's chunk of it.
        if let Some(&(start, end)) = masks.last().filter(|&&(_, end)| end as usize > length) {
            return Some(HardMaskShape::PastEnd { start, end, length });
        }
        let frame_length =
            |frame: i8| (length + 1).saturating_sub(usize::from(frame.unsigned_abs())) / 3;
        let mut nucleotide = std::mem::take(&mut split.ranges);
        complement_masks(masks, length, &mut nucleotide);
        let mut found = None;
        'frames: for frame in [1i8, 2, 3, -1, -2, -3] {
            let frame_len = frame_length(frame);
            if frame_len == 0 {
                continue;
            }
            split.start(
                frame_len,
                translated_frame_ranges(&nucleotide, frame, frame_length(-1)),
                false,
            );
            // The plus frames' ranges are in the order of the masks, the minus frames' in
            // reverse; range `r` ends where its mask starts (the last range has none).
            let mask_after = |range: usize| {
                let index = if frame > 0 {
                    range
                } else {
                    (masks.len() - 1).saturating_sub(range)
                };
                masks[index.min(masks.len() - 1)]
            };
            let mut first = true;
            while let Some(chunk) = split.next_chunk(DBSEQ_CHUNK_OVERLAP) {
                if chunk.full_range && !first {
                    // Only ranges at 0 come before it, so it is not the frame's first range.
                    let (start, end) = mask_after(chunk.range.saturating_sub(1));
                    found = Some(HardMaskShape::FrameStartMask { frame, start, end });
                    break 'frames;
                }
                first = false;
            }
        }
        split.end();
        nucleotide.clear();
        split.ranges = nucleotide;
        found
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
            HardMaskShape::FrameStartMask { frame, start, end } => format!(
                "has a mask, [{start}, {end}), that leaves reading frame {frame:+} an empty part \
                 at its start, followed by a part from its start that NCBI BLAST+ 2.17.0 searches \
                 as the whole frame: NCBI then frees memory that it did not allocate (a crash in \
                 every run)"
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
        }
    }
}
