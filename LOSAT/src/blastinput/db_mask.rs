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
//! LOSAT's BLASTN searches with `-db_soft_mask` and `-db_hard_mask` (series NP, M4a-2a and M4a-2b);
//! it refuses the databases of `SeqDb::check_search_masks` and, with `-db_hard_mask`, a searched
//! sequence of a `HardMaskShape` (`ResolvedMask::check_losat_limits`).

use anyhow::Result;

use crate::cli::NativeError;
use crate::core::subject_split::HardSplitState;
use crate::report::archive::options::{Param, ParamValue};
use crate::seqdb::{MaskAlgorithm, NameMaskError, SeqDb, SeqDbError};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2342-2356
// ```c++
//     // Masking of database
//     if (m_SupportsDatabaseMasking) {
//         arg_desc.AddOptionalKey(kArgDbSoftMask,
//                 "filtering_algorithm",
//                 "Filtering algorithm ID to apply to the BLAST database as soft "
//                 "masking",
//                 CArgDescriptions::eString);
//         arg_desc.SetDependency(kArgDbSoftMask, CArgDescriptions::eExcludes,
//                            kArgDbHardMask);
//
//         arg_desc.AddOptionalKey(kArgDbHardMask,
//                 "filtering_algorithm",
//                 "Filtering algorithm ID to apply to the BLAST database as hard "
//                 "masking",
//                 CArgDescriptions::eString);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2226-2229,2367-2382
// ```c++
//     if (m_SupportsDatabaseMasking) {
//         database_args.push_back(kArgDbSoftMask);
//         database_args.push_back(kArgDbHardMask);
//     }
//     ...
//         ITERATE(vector<string>, dbarg, database_args) {
//             arg_desc.SetDependency(kArgSubject, CArgDescriptions::eExcludes,
//                                    *dbarg);
//         }
//     ...
//         ITERATE(vector<string>, dbarg, database_args) {
//             arg_desc.SetDependency(kArgSubjectLocation,
//                                    CArgDescriptions::eExcludes,
//                                    *dbarg);
//         }
// ```
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
// ```c
// /** Define the possible subject masking types */
// typedef enum ESubjectMaskingType {
//     eNoSubjMasking,
//     eSoftSubjMasking,
//     eHardSubjMasking
// } ESubjectMaskingType;
// ```
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
    // ```c++
    //         if (args.Exist(kArgDbSoftMask) && args[kArgDbSoftMask]) {
    //             m_SearchDb->SetFilteringAlgorithm(args[kArgDbSoftMask].AsString(), eSoftSubjMasking);
    //         } else if (args.Exist(kArgDbHardMask) && args[kArgDbHardMask]) {
    //             m_SearchDb->SetFilteringAlgorithm(args[kArgDbHardMask].AsString(), eHardSubjMasking);
    //         }
    // ```
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
    // ```c++
    // static bool s_IsNumericId(const string &id)
    // {
    //     Int4 nid(-1);
    //     return NStr::StringToNumeric(id, &nid, NStr::fConvErr_NoThrow, 10);
    // }
    //
    // void
    // CSearchDatabase::SetFilteringAlgorithm(const string &filt_algorithm,
    //                                        ESubjectMaskingType mask_type)
    // {
    //     m_MaskType = mask_type;
    //     m_FilteringAlgorithmString = "";
    //     if (mask_type == eNoSubjMasking) {
    //         m_FilteringAlgorithmId = -1;
    //         return;
    //     }
    //     if (s_IsNumericId(filt_algorithm)) {
    //         m_FilteringAlgorithmId = NStr::StringToInt(filt_algorithm);
    //         x_ValidateMaskingAlgorithm();
    //     } else {
    //         m_FilteringAlgorithmString = filt_algorithm;
    //         m_NeedsFilteringTranslation = true;
    //     }
    // }
    // ```
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:296-318
    // ```c++
    //     x_ValidateMaskingAlgorithm();
    //     _ASSERT(m_SeqDb.NotEmpty());
    //     m_DbInitialized = true;
    // }
    //
    // void
    // CSearchDatabase::x_ValidateMaskingAlgorithm() const
    // {
    //     if (m_FilteringAlgorithmId <= 0 || m_SeqDb.Empty()) {
    //         return;
    //     }
    //
    //     vector<int> supported_algorithms;
    //     m_SeqDb->GetAvailableMaskAlgorithms(supported_algorithms);
    //     if (find(supported_algorithms.begin(),
    //              supported_algorithms.end(),
    //              m_FilteringAlgorithmId) == supported_algorithms.end()) {
    //         CNcbiOstrstream oss;
    //         oss << "Masking algorithm ID " << m_FilteringAlgorithmId << " is "
    //             << "not supported in " <<
    //             (IsProtein() ? "protein" : "nucleotide") << " '"
    //             << GetDatabaseName() << "' BLAST database";
    //         string msg = CNcbiOstrstreamToString(oss);
    //         NCBI_THROW(CBlastException, eInvalidOptions, msg);
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-220
    // ```c++
    //     catch (const blast::CBlastException& e) {                               \
    //         const string& msg = e.GetMsg();                                     \
    //         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
    //             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
    //             exit_code = BLAST_INPUT_ERROR;                                  \
    // ```
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
                    "BLAST options error: Masking algorithm ID {id} is not supported in {} '{db_name}' BLAST database\n",
                    if nucleotide { "nucleotide" } else { "protein" }
                ),
            }
            .into());
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:195-231
    // ```c++
    // int
    // CSearchDatabase::GetFilteringAlgorithm() const
    // {
    //     if (m_MaskType == eNoSubjMasking) return -1;
    //     try {
    //         if (m_NeedsFilteringTranslation) {
    //             x_TranslateFilteringAlgorithm();
    //         }
    //         return m_FilteringAlgorithmId;
    //     } catch (const CSeqDBException&) {
    //         return -1;
    //     }
    // }
    //
    // string
    // CSearchDatabase::GetFilteringAlgorithmKey() const
    // {
    //     if (m_MaskType == eNoSubjMasking) return "";
    //     return m_FilteringAlgorithmString;
    // }
    // ...
    // void
    // CSearchDatabase::x_TranslateFilteringAlgorithm() const
    // {
    //     if (m_MaskType == eNoSubjMasking) return;
    //     if (!m_DbInitialized) {
    //         x_InitializeDb();
    //     }
    //     m_FilteringAlgorithmId =
    //         m_SeqDb->GetMaskAlgorithmId(m_FilteringAlgorithmString);
    //     m_NeedsFilteringTranslation = false;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153
    // ```c++
    //         int filteringAlgorithmId = db_adapter.GetFilteringAlgorithm();
    //         if(filteringAlgorithmId == -1) {
    //             CRef <CSearchDatabase> db_Info = db_adapter.GetSearchDatabase();
    //             if (db_Info && db_Info.NotEmpty()) {
    //                 ESubjectMaskingType maskType = db_Info->GetMaskType();
    //                 if(maskType != eNoSubjMasking) {
    //                     db_Info->SetFilteringAlgorithm(-1, eNoSubjMasking);
    //                     ERR_POST(Warning << "Subject mask not found in " + m_DbName +", proceeding without subject masking.");
    //                 }
    //             }
    //         }
    //         CBlastFormatUtil::GetBlastDbInfo(m_DbInfo, m_DbName, m_DbIsAA,
    //                                    dbfilt_algorithm, is_remote_search);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:228-242
    // ```c++
    //         CBlastFormat formatter(opt, *db_adapter,
    //     ...
    //                                m_CmdLineArgs->ExecuteRemotely(),
    //                                db_adapter->GetFilteringAlgorithm(),
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:196-200
    // ```c++
    //         else {                                                              \
    //         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
    //         	exit_code = BLAST_DATABASE_ERROR;                               \
    //         }                                                                   \
    // ```
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
        let shown = db.name_mask_details(id).map_err(|error| {
            let error = match &error {
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
                    message: format!("BLAST Database error: {text}\n"),
                }
                .into(),
                error => crate::algorithm::common::database::database_error_of(error, program),
            }
        })?;
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
    // ```c++
    //     CNcbiOstream* out = cmdline_args->GetExportSearchStrategyStream(args);
    //     if ( !out ) {
    //         return;
    //     }
    //
    //     s_ExportSearchStrategy(out, queries, opts_hndl,
    //                            cmdline_args->GetBlastDatabaseArgs(),
    //                            pssm, num_iters);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:54-58
    // ```c++
    // CLocalDbAdapter::CLocalDbAdapter(const CSearchDatabase& dbinfo)
    //     : m_SeqSrc(0), m_SeqInfoSrc(0), m_DbName(dbinfo.GetDatabaseName()), m_DbScanMode(false)
    // {
    //     m_DbInfo.Reset(new CSearchDatabase(dbinfo));
    // }
    // ```
    // The search strategy reads the command line's `CSearchDatabase`; the formatter's warning
    // resets the database adapter's copy (blast_format.cpp:146), which the archive reads. So
    // after the warning the strategy still writes the name as its key (search_strategy.cpp:500-505,
    // quoted at `ResolvedMask::archive_params`); an empty name and the id -1 write nothing.
    /// The search strategy's program options of the masking (`-export_search_strategy`): the
    /// archive's (`ResolvedMask::archive_params`) when the formatter resolved it (`resolved`), else,
    /// after its warning, `DbFilteringAlgorithmKey` and `SubjectMaskingType` of a name that is not
    /// empty.
    pub fn strategy_params(&self, resolved: Option<&ResolvedMask>) -> Vec<Param> {
        match (resolved, &self.algorithm) {
            (Some(mask), _) => mask.archive_params(),
            (None, FilteringAlgorithm::Name(name)) if !name.is_empty() => vec![
                Param {
                    name: "DbFilteringAlgorithmKey",
                    value: ParamValue::Str(name.clone()),
                },
                Param {
                    name: "SubjectMaskingType",
                    value: ParamValue::Integer(self.kind as i32),
                },
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
    /// runs to the end only below `TAIL_MASK_SEARCHED_WORD`), for a search whose scan reads words
    /// of `scanned_word` letters (`HardMaskShape::of`).
    pub fn check_losat_limits(&self, db: &SeqDb, program: &str, scanned_word: usize) -> Result<()> {
        db.check_search_masks(self.id)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        if self.kind == MaskKind::Hard {
            let found = db.first_searched_masks(self.id, |masks, bounds, length| {
                HardMaskShape::of(masks, bounds, length, scanned_word)
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
    // ```c++
    //     // Set the filtering algorithms
    //     string algo_key = db->GetFilteringAlgorithmKey();
    //     if (algo_key != kEmptyStr)
    //     {
    //     	int mask_type = (int) db->GetMaskType();
    //        	x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_DbFilteringAlgorithmKey), algo_key);
    //        	x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_SubjectMaskingType), mask_type);
    //     }
    //     else {
    //     	// Set the filtering algorithms
    //     	int algo_id = db->GetFilteringAlgorithm();
    //     	if (algo_id != -1)
    //     	{
    //     		int mask_type = (int) db->GetMaskType();
    //     		x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_DbFilteringAlgorithmId), algo_id);
    //     		x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_SubjectMaskingType), mask_type);
    //     	}
    //     }
    // ```
    /// The program options of the BLAST archive and of the search strategy (after the list's).
    /// An empty name is no key (`algo_key != kEmptyStr`): the id it was translated to is written.
    pub fn archive_params(&self) -> Vec<Param> {
        let algorithm = match self.key.as_deref() {
            Some(key) if !key.is_empty() => Param {
                name: "DbFilteringAlgorithmKey",
                value: ParamValue::Str(key.to_string()),
            },
            _ => Param {
                name: "DbFilteringAlgorithmId",
                value: ParamValue::Integer(self.id),
            },
        };
        vec![
            algorithm,
            Param {
                name: "SubjectMaskingType",
                value: ParamValue::Integer(self.kind as i32),
            },
        ]
    }
}

#[cfg(test)]
#[path = "db_mask_tests.rs"]
mod db_mask_tests;

/// A searched sequence's masks that `-db_hard_mask` refuses (`ResolvedMask::check_losat_limits`;
/// series NP, M4a-2b, ncbi-parity-20261019/M/m4a2/SURVEY.md §7.1-7.3): NCBI's hard chunks
/// (`s_GetNextSubjectChunk`, blast_engine.c:221-310) give no result that LOSAT can check for them.
/// The shapes are read off the chunks themselves (`core::subject_split::HardSplitState`, the
/// chunker that the search runs; review S1 of M4a-2b), except `Unsorted`, a check of the stored
/// masks. The first shape is refused only at a scanned word below [`TAIL_MASK_SEARCHED_WORD`]
/// (blastn-short): at that word and above LOSAT searches it over the chunk before the mask, an
/// approved NCBI defect (series NP, M4a-2c; the Owner's answer of 2026-11-15 to question Q1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardMaskShape {
    /// The mask `[start, end)` runs to the subject's end after a chunk from the subject's start
    /// (only empty masks at 0 before it), which takes NCBI's "no chunking" branch and so scans
    /// the full range past the chunk, in a subject that the scan reads a word of, at a scanned
    /// word below [`TAIL_MASK_SEARCHED_WORD`].
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
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-268 (quoted at
// `core::subject_split::HardSplitState::next_chunk`)
// ```c
//     /* if no chunking is performed */
//     if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
//         subject->seq_ranges = backup->soft_ranges;
//         subject->num_seq_ranges = backup->num_soft_ranges;
//         return SUBJECT_SPLIT_OK;
//     }
// ```
/// The shortest scanned word (the word size, or the discontiguous template's length) at which
/// `-db_hard_mask` searches a sequence of the shape [`HardMaskShape::TailMask`] instead of refusing
/// it: over the chunk before the mask alone, scanned from its start as NCBI scans its other hard
/// chunks (`core::subject_split::HardSplitState::next_chunk`). An approved NCBI defect
/// (PD-LOSAT-NCBI-DEFECTS exception 13, version 1.34; the Owner's decision of 2026-11-15, 「M4a
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
/// what NCBI's normal chunks give for the near input (the exception's text).
pub const TAIL_MASK_SEARCHED_WORD: usize = 11;

impl HardMaskShape {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276 and 579-584
    // (quoted at `core::subject_split::HardSplitState::next_chunk`): the chunks, the "no
    // chunking" branch, and the overlap after each chunk, read one range past the list after the
    // last one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1668
    // ```c
    //     /* if sequence is masked, fall back to generic scanner and extender */
    //     if (subject->mask_type != eNoSubjMasking) {
    //     ...
    //         /* generic scanner permits any (non-aligned) starting offset */
    //         scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
    //         scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354, 1368-1370
    // ```c
    //     Int4 min_subj_seq_length = 1;
    //     ...
    //     if (Blast_SubjectIsTranslated(program_number)) {
    //         min_subj_seq_length = s_GetMinimumSubjSeqLen(lookup_wrap);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
    // ```c
    // s_DetermineScanningOffsets(const BLAST_SequenceBlk* subject,
    //                            Int4  word_length,
    //                            Int4  lut_word_length,
    //                            Int4* range)
    // {
    //     ...
    //     while (range[1] > range[2]) {
    //         range[0]++;
    //         if (range[0] >= (Int4)subject->num_seq_ranges) {
    //             return FALSE;
    //         }
    //         range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
    //         range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
    //     }
    //     return TRUE;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500, 565
    // ```c
    //         if (aux_struct->WordFinder) {
    //             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
    //                                    ...);
    //
    //             if (init_hitlist->total == 0) continue;
    //         }
    //     ...
    //         if (hsp_list->hspcnt == 0) continue;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2873-2911, 2956
    // ```c
    //    if (!hsp_list || hsp_list->hspcnt == 0)
    //       return 0;
    //
    //    /* If no previous HSP list, just return a copy of the new one. */
    //    if (!combined_hsp_list) {
    //    ...
    //    if (contexts_per_query < 0) {      /* subject seq is split */
    //       for (index1 = 0; index1 < combined_hsp_list->hspcnt; index1++) {
    //          hsp1 = combined_hsp_list->hsp_array[index1];
    //          if (hsp1->subject.end > split_offsets[0]) {
    //       ...
    //       for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
    //          hsp2 = hsp_list->hsp_array[index2];
    //          if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) {
    //    ...
    //    if (hspcnt1 > 0 && hspcnt2 > 0) {
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1527-1537
    // (quoted at blastn's `database::mask_hard_letters`): the partial fetch's slices are masked by
    // a scan that assumes sorted masks.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1326-1333
    // ```c
    //     } else if (program_number == eBlastTypeBlastn &&
    //                options->word_size > DBSEQ_CHUNK_OVERLAP) {
    //         ...
    //                   "Word-size must be less than or equal to %d", DBSEQ_CHUNK_OVERLAP);
    //         ...
    //         return BLASTERR_OPTION_VALUE_INVALID;
    // ```
    /// The shape of the stored masks `masks` (`[start, end)`, as stored, not empty) of a sequence
    /// whose length `length` gives, that `-db_hard_mask` refuses in a search whose scan reads words
    /// of `scanned_word` letters (the word size, or the discontiguous template's length), or
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
    pub fn of(
        masks: &[(u32, u32)],
        bounds: Option<(usize, usize)>,
        length: &dyn Fn() -> Option<usize>,
        scanned_word: usize,
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
        // The masks are sorted: the last one ends last.
        let last_end = masks[masks.len() - 1].1 as usize;
        if let Some((low, high)) = bounds.filter(|&(low, _)| last_end < low) {
            let shape = Self::of_chunks(masks, low, scanned_word, &mut split);
            if Self::of_chunks(masks, high, scanned_word, &mut split) == shape {
                return Ok(shape);
            }
        }
        let length = length().ok_or(())?;
        Ok(Self::of_chunks(masks, length, scanned_word, &mut split))
    }

    /// The shape of the sorted masks `masks` of a sequence of `length` letters, read off its
    /// chunks (`HardMaskShape::of`); `split` is the chunker, whose buffer of ranges is reused.
    fn of_chunks(
        masks: &[(u32, u32)],
        length: usize,
        scanned_word: usize,
        split: &mut HardSplitState,
    ) -> Option<HardMaskShape> {
        use crate::core::subject_split::{complement_masks, HardChunk, DBSEQ_CHUNK_OVERLAP};
        let mut ranges = std::mem::take(&mut split.ranges);
        complement_masks(masks, length, &mut ranges);
        split.start(length, ranges);
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
            if chunk.full_range
                && chunk.length < length
                && scans
                && scanned_word < TAIL_MASK_SEARCHED_WORD
            {
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

    /// The refusal's words for the shape.
    pub fn text(&self) -> String {
        match *self {
            HardMaskShape::TailMask { start, end } => format!(
                "has a mask, [{start}, {end}), running to its end after the chunk from its start, \
                 for which NCBI BLAST+ 2.17.0 scans the masked letters past that chunk while it \
                 extends within the chunk (at a word of fewer than {TAIL_MASK_SEARCHED_WORD} \
                 letters: rows that differ from a near input's, or a crash)"
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
        }
    }
}
