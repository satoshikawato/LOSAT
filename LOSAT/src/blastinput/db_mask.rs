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
//! LOSAT's BLASTN searches with `-db_soft_mask`; it refuses a search that `-db_hard_mask` would
//! mask (series NP, M4a-2b) and the databases of `SeqDb::check_search_masks`.

use anyhow::Result;

use crate::cli::NativeError;
use crate::report::archive::options::{Param, ParamValue};
use crate::seqdb::{MaskAlgorithm, SeqDb, SeqDbError};

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
    /// `-db_hard_mask` (refused where it would mask a search: `DbMaskOption::check_losat_limits`).
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

impl MaskKind {
    /// The option's name, as LOSAT's refusals name it.
    fn option(self) -> &'static str {
        match self {
            MaskKind::Soft => "db_soft_mask",
            MaskKind::Hard => "db_hard_mask",
        }
    }
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
    /// database of each name of the `-db` list (`SeqDb::name_mask_details`): an id that one of
    /// them does not have is NCBI's database error, exit 2 (the id 0, an id below -1, or a
    /// name's id that a database of the list lacks).
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
            .name_mask_details(db_name, id)
            .map_err(|error| match error {
                SeqDbError::Range(text) => NativeError {
                    exit: 2,
                    message: format!("BLAST Database error: {text}\n"),
                }
                .into(),
                error => crate::algorithm::common::database::database_error(error, program),
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
    /// checks; `program` names the program): `-db_hard_mask` (`cli::run_refusal`; series NP,
    /// M4a-2b), and the databases whose masks NCBI reads from undefined or reversed data
    /// (`SeqDb::check_search_masks`; decision D5 of M4a-2).
    pub fn check_losat_limits(&self, db: &SeqDb, program: &str) -> Result<()> {
        if self.kind == MaskKind::Hard {
            anyhow::bail!("{}", crate::cli::run_refusal(program, self.kind.option()));
        }
        db.check_search_masks(self.id)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))
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
