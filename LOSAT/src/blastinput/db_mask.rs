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
    /// do not handle (`HardMaskShape`; series NP, M4a-2b decisions H3-H5), for a search whose scan
    /// reads words of `scanned_word` letters (`HardMaskShape::of`).
    pub fn check_losat_limits(&self, db: &SeqDb, program: &str, scanned_word: usize) -> Result<()> {
        db.check_search_masks(self.id)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        if self.kind == MaskKind::Hard {
            let found = db.first_searched_masks(self.id, |masks, length| {
                HardMaskShape::of(masks, length, scanned_word)
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
/// The one place of these rules: the answer of Owner question Q1 for the first shape (an approved
/// NCBI defect instead of the refusal) removes `TailMask` here and changes the chunker's "no
/// chunking" branch for hard ranges (`SubjectSplitState::next_chunk`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardMaskShape {
    /// The first mask `[start, end)` runs to the subject's end and starts before it, at or below
    /// `MAX_DBSEQ_LEN`, in a subject that the scan reads a word of.
    TailMask { start: u32, end: u32 },
    /// A mask that ends before it starts, or that starts before the end of the mask before it.
    Unsorted {
        before: (u32, u32),
        mask: (u32, u32),
    },
    /// The last mask `[start, end)` ends before the subject's end and is shorter than the
    /// `end % 4` letters that the last chunk starts with.
    ShortLastMask { start: u32, end: u32 },
}

impl HardMaskShape {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:241-268
    // ```c
    //     if (backup->offset + MAX_DBSEQ_LEN <
    //         backup->hard_ranges[backup->hm_index].right) {
    //
    //         subject->length = MAX_DBSEQ_LEN;
    //         backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    //
    //     } else {
    //
    //         subject->length = backup->hard_ranges[backup->hm_index].right
    //                         - backup->offset;
    //         backup->hm_index++;
    //         backup->next = (backup->hm_index < backup->num_hard_ranges) ?
    //                         backup->hard_ranges[backup->hm_index].left :
    //                         backup->full_range.right;
    //     }
    //
    //     (subject->chunk)++;
    //
    //     /* if no chunking is performed */
    //     if (backup->offset == 0 && residual == 0 && backup->next == backup->full_range.right) {
    //         subject->seq_ranges = backup->soft_ranges;
    //         subject->num_seq_ranges = backup->num_soft_ranges;
    //         return SUBJECT_SPLIT_OK;
    //     }
    // ```
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1318-1331, 1429-1432
    // ```c
    // static Int4 s_GetMinimumSubjSeqLen(LookupTableWrap* lookup_wrap)
    // {
    //     Int4 word_length = 1;
    //     ...
    //     return word_length * 3 + 2;
    // }
    //     ...
    //        if (seq_arg.seq->length < min_subj_seq_length) {
    //            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
    //            continue;
    //        }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:579-584
    // ```c
    //         Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
    //         overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
    //                   0 : dbseq_chunk_overlap;
    //         status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
    //                      kHspNumMax, &(backup.offset), INT4_MIN,
    //                      overlap, score_options->gapped_calculation,
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1527-1537
    // (quoted at blastn's `database::mask_hard_letters`): the partial fetch's slices are masked by
    // a scan that assumes sorted masks.
    /// The shape of the stored masks `masks` (`[start, end)`, as stored, not empty) of a sequence
    /// whose length `length` gives, that `-db_hard_mask` refuses in a search whose scan reads words
    /// of `scanned_word` letters (the word size, or the discontiguous template's length), or
    /// `None`; `Err(())` when the length cannot be read (`SeqDb::first_searched_masks`). The hard
    /// ranges are the complement in stored order (`BlastSeqBlkSetSeqRanges`), whose first starts
    /// at 0:
    /// - the first chunk (offset 0, residual 0) ends where the first mask starts unless that is
    ///   above `MAX_DBSEQ_LEN`; when the first mask runs to the end, `next` is the end and the
    ///   "no chunking" branch scans the full range `[0, length)` of a chunk of `start` letters:
    ///   NCBI reads seeds past the chunk (SURVEY §7.1: rows like the input with the mask one letter
    ///   shorter, or a corrupt heap; an empty last mask after it changes nothing). A subject
    ///   shorter than the word (the scan of `[0, length)` reads none) or than 5 letters (the
    ///   engine skips it, `min_subj_seq_length` of a nucleotide lookup) is not scanned;
    /// - masks out of order give hard ranges that go back (negative chunk lengths, which move the
    ///   diagonal offset back) and slices that NCBI's partial fetch masks by a scan of sorted
    ///   masks (§7.2);
    /// - after the last hard range NCBI reads the overlap from one range past the list (another
    ///   sequence's mask end, or unset memory, §7.3); it decides the merge when the previous
    ///   chunk's HSPs can end after the last chunk's offset `end - end % 4`, i.e. when the last
    ///   mask is shorter than `end % 4` (the last hard range is not empty: `end` before the
    ///   subject's end).
    pub fn of(
        masks: &[(u32, u32)],
        length: &dyn Fn() -> Option<usize>,
        scanned_word: usize,
    ) -> Result<Option<HardMaskShape>, ()> {
        const MAX_DBSEQ_LEN: u32 = 5_000_000;
        let Some(&(first_start, first_end)) = masks.first() else {
            return Ok(None);
        };
        // The length, read once and only for the shapes that need it.
        let mut known: Option<usize> = None;
        let mut length = || -> Result<usize, ()> {
            if known.is_none() {
                known = Some(length().ok_or(())?);
            }
            Ok(known.unwrap_or_default())
        };
        // A first mask that runs to the end leaves only empty masks after it (or masks out of
        // order): read the length only then.
        let tail_candidate = first_start <= MAX_DBSEQ_LEN
            && first_start < first_end
            && masks[1..]
                .iter()
                .all(|&(start, end)| start == first_end && end == first_end);
        if tail_candidate && {
            let length = length()?;
            length == first_end as usize && length >= scanned_word.max(5)
        } {
            return Ok(Some(HardMaskShape::TailMask {
                start: first_start,
                end: first_end,
            }));
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
        let (last_start, last_end) = masks[masks.len() - 1];
        if last_end - last_start < last_end % 4 && (last_end as usize) < length()? {
            return Ok(Some(HardMaskShape::ShortLastMask {
                start: last_start,
                end: last_end,
            }));
        }
        Ok(None)
    }

    /// The refusal's words for the shape.
    pub fn text(&self) -> String {
        match *self {
            HardMaskShape::TailMask { start, end } => format!(
                "has its only mask, [{start}, {end}), running to its end, for which NCBI BLAST+ \
                 2.17.0 scans the masked letters past the chunk before the mask (it reads past \
                 its buffer: rows of a near input, or a corrupt heap)"
            ),
            HardMaskShape::Unsorted { before, mask } => format!(
                "has masks out of order ([{}, {}) after [{}, {})), whose hard masking NCBI BLAST+ \
                 2.17.0 ties to the slices of its partial fetch",
                mask.0, mask.1, before.0, before.1
            ),
            HardMaskShape::ShortLastMask { start, end } => format!(
                "has a last mask, [{start}, {end}), shorter than the {} letters before its end \
                 that the last chunk starts with, after which NCBI BLAST+ 2.17.0 reads the \
                 chunks' overlap past its list of ranges",
                end % 4
            ),
        }
    }
}
