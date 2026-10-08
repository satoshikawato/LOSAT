//! Coordination module for BLASTN pipeline setup
//!
//! This module handles:
//! - Task-specific parameter configuration (word size, scoring, etc.)
//! - Sequence reading and preprocessing
//! - DUST masking
//! - Lookup table building

use super::args::BlastnArgs;
use super::constants::{
    LUT_WIDTH_11_THRESHOLD_10, LUT_WIDTH_11_THRESHOLD_8, MAX_DIRECT_LOOKUP_WORD_SIZE,
    MIN_DIAG_SEPARATION_BLASTN, MIN_DIAG_SEPARATION_MEGABLAST, MIN_UNGAPPED_SCORE_BLASTN,
    MIN_UNGAPPED_SCORE_MEGABLAST, SCAN_RANGE_BLASTN, SCAN_RANGE_MEGABLAST, X_DROP_GAPPED_FINAL,
    X_DROP_GAPPED_GREEDY, X_DROP_GAPPED_NUCL,
};
use super::disc_lookup::DiscWordType;
use super::lookup::{
    build_db_word_counts, build_disc_mb_lookup, build_na_lookup, build_pv_direct_lookup,
    build_two_stage_lookup, NaLookupTable, PvDirectLookup, TwoStageLookup,
};
use crate::blastinput::fasta_reader::{FastaRecord, InputRecord};
use crate::utils::dust::{DustMasker, MaskedInterval};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LookupTableKind {
    Small,
    Mb,
    Na,
}

/// Task-specific configuration for BLASTN
pub struct TaskConfig {
    pub effective_word_size: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_ungapped_score: i32,
    pub use_dp: bool,
    pub scan_step: usize,
    pub use_direct_lookup: bool,
    pub use_two_stage: bool,
    pub lut_word_length: usize,
    /// Whether NCBI chose its small-query lookup table (`eSmallNaLookupTable`), whose
    /// word extension reads the compressed query.
    pub small_na_lookup: bool,
    /// Whether NCBI chose its megablast lookup table (`eMBLookupTable`), whose chains
    /// list query offsets newest first; the small and standard tables list them in
    /// ascending order.
    pub mb_lookup: bool,
    pub x_drop_gapped: i32, // Task-specific gapped X-dropoff (blastn: 30, megablast: 25)
    pub x_drop_final: i32,  // Final traceback X-dropoff (100 for all nucleotide tasks)
    pub scan_range: usize,  // Off-diagonal scan range (BLAST_SCAN_RANGE_NUCL: 0 for every task)
    pub min_diag_separation: i32, // NCBI reference: blast_nucl_options.cpp:239,259 (blastn: 50, megablast: 6)
    /// The two-hit window (`window_size`): 40 for dc-megablast, 0 (one hit) otherwise.
    pub window_size: usize,
    /// The discontiguous template (`mb_template_type`, `mb_template_length`); a length
    /// of 0 is a contiguous word.
    pub mb_template_type: DiscWordType,
    pub mb_template_length: u8,
}

/// Lookup tables for seed finding
pub struct LookupTables {
    pub two_stage_lookup: Option<TwoStageLookup>,
    pub pv_direct_lookup: Option<PvDirectLookup>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:131-156
    // ```c
    // typedef struct BlastNaLookupTable {
    //     ...
    // } BlastNaLookupTable;
    // ```
    pub na_lookup: Option<NaLookupTable>,
}

/// Subject metadata needed for search setup.
#[derive(Clone)]
pub struct SubjectMetadata {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    pub db_len_total: usize,
    pub db_num_seqs: usize,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    pub subject_ids: Vec<String>,
}

/// Sequence data and metadata (`R`: the query records, `InputRecord`)
pub struct SequenceData<R = FastaRecord> {
    pub queries: Vec<R>,
    pub query_ids: Vec<String>,
    pub query_masks: Vec<Vec<MaskedInterval>>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    pub db_len_total: usize,
    pub db_num_seqs: usize,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    pub subject_ids: Vec<String>,
}

/// The option values of a task, which the command line options replace (NCBI
/// `CBlastOptionsFactory::CreateTask`). The CLI accepts the tasks megablast, blastn,
/// dc-megablast and blastn-short.
///
/// NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:344-380
/// ```c
///     if (!NStr::CompareNocase(task, "blastn") ||
///         !NStr::CompareNocase(task, "blastn-short") ||
///         // -RMH-
///         !NStr::CompareNocase(task, "rmblastn") ||
///         !NStr::CompareNocase(task, "vecscreen"))
///     {
///         CBlastNucleotideOptionsHandle* opts =
///              dynamic_cast<CBlastNucleotideOptionsHandle*>
///                 (CBlastOptionsFactory::Create(eBlastn, locality));
///         _ASSERT(opts);
///         if (!NStr::CompareNocase(task, "blastn-short"))
///         {
///              opts->SetMatchReward(1);
///              opts->SetMismatchPenalty(-3);
///              opts->SetEvalueThreshold(1000);
///              opts->SetWordSize(7);
///              opts->ClearFilterOptions();
///         }
///         ...
///     }
///     else if (!NStr::CompareNocase(task, "megablast"))
///     {
///          retval = CBlastOptionsFactory::Create(eMegablast, locality);
///     }
///     else if (!NStr::CompareNocase(task, "dc-megablast"))
///     {
///          retval = CBlastOptionsFactory::Create(eDiscMegablast, locality);
///     }
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_nucl_options.cpp:137-150,153-170,173-194,198-221,223-262
/// ```c
/// CBlastNucleotideOptionsHandle::SetLookupTableDefaults()
/// {
///     SetLookupTableType(eNaLookupTable);
///     SetWordSize(BLAST_WORDSIZE_NUCL);
/// ...
/// CBlastNucleotideOptionsHandle::SetMBLookupTableDefaults()
/// {
///     SetLookupTableType(eMBLookupTable);
///     SetWordSize(BLAST_WORDSIZE_MEGABLAST);
/// ...
/// CBlastNucleotideOptionsHandle::SetQueryOptionDefaults()
/// {
///     SetDustFiltering(true);
/// ...
/// CBlastNucleotideOptionsHandle::SetInitialWordOptionsDefaults()
/// {
///     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL);
///     SetWindowSize(BLAST_WINDOW_SIZE_NUCL);
///     SetOffDiagonalRange(BLAST_SCAN_RANGE_NUCL);
/// }
/// ...
/// CBlastNucleotideOptionsHandle::SetGappedExtensionDefaults()
/// {
///     SetGapXDropoff(BLAST_GAP_X_DROPOFF_NUCL);
///     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_NUCL);
///     SetGapTrigger(BLAST_GAP_TRIGGER_NUCL);
///     SetGapExtnAlgorithm(eDynProgScoreOnly);
///     SetGapTracebackAlgorithm(eDynProgTbck);
/// }
/// ...
/// CBlastNucleotideOptionsHandle::SetMBGappedExtensionDefaults()
/// {
///     SetGapXDropoff(BLAST_GAP_X_DROPOFF_GREEDY);
///     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_NUCL);
///     SetGapTrigger(BLAST_GAP_TRIGGER_NUCL);
///     SetGapExtnAlgorithm(eGreedyScoreOnly);
///     SetGapTracebackAlgorithm(eGreedyTbck);
/// }
/// CBlastNucleotideOptionsHandle::SetScoringOptionsDefaults()
/// {
///     SetMatrixName(NULL);
///     SetGapOpeningCost(BLAST_GAP_OPEN_NUCL);
///     SetGapExtensionCost(BLAST_GAP_EXTN_NUCL);
///     SetMatchReward(2);
///     SetMismatchPenalty(-3);
/// ...
/// CBlastNucleotideOptionsHandle::SetMBScoringOptionsDefaults()
/// {
///     SetMatrixName(NULL);
///     SetGapOpeningCost(BLAST_GAP_OPEN_MEGABLAST);
///     SetGapExtensionCost(BLAST_GAP_EXTN_MEGABLAST);
///     SetMatchReward(1);
///     SetMismatchPenalty(-2);
/// ...
/// CBlastNucleotideOptionsHandle::SetHitSavingOptionsDefaults()
/// {
///     SetHitlistSize(500);
///     SetEvalueThreshold(BLAST_EXPECT_VALUE);
///     ...
///     SetMinDiagSeparation(50);
/// ...
/// CBlastNucleotideOptionsHandle::SetMBHitSavingOptionsDefaults()
/// {
///     SetHitlistSize(500);
///     SetEvalueThreshold(BLAST_EXPECT_VALUE);
///     ...
///     SetMinDiagSeparation(6);
/// ```
/// dc-megablast is the megablast handle with the discontiguous overrides; the base
/// constructor's calls bind to the base class, and the derived constructor's
/// `SetDefaults` then calls the overrides (the hit-saving defaults are not overridden).
///
/// NCBI reference: c++/src/algo/blast/api/disc_nucl_options.cpp:47-90
/// ```c
/// CDiscNucleotideOptionsHandle::CDiscNucleotideOptionsHandle(EAPILocality locality)
///     : CBlastNucleotideOptionsHandle(locality)
/// {
///     SetDefaults();
///     m_Opts->SetProgram(eDiscMegablast);
/// }
///
/// void
/// CDiscNucleotideOptionsHandle::SetMBLookupTableDefaults()
/// {
///     CBlastNucleotideOptionsHandle::SetMBLookupTableDefaults();
///     bool defaults_mode = m_Opts->GetDefaultsMode();
///     m_Opts->SetDefaultsMode(false);
///     SetTemplateType(0);
///     SetTemplateLength(18);
///     SetWordSize(BLAST_WORDSIZE_NUCL);
///     m_Opts->SetDefaultsMode(defaults_mode);
/// }
///
/// void
/// CDiscNucleotideOptionsHandle::SetMBInitialWordOptionsDefaults()
/// {
///     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL);
///     bool defaults_mode = m_Opts->GetDefaultsMode();
///     m_Opts->SetDefaultsMode(false);
///     SetWindowSize(BLAST_WINDOW_SIZE_DISC);
///     m_Opts->SetDefaultsMode(defaults_mode);
/// }
///
/// void
/// CDiscNucleotideOptionsHandle::SetMBGappedExtensionDefaults()
/// {
///     SetGapXDropoff(BLAST_GAP_X_DROPOFF_NUCL);
///     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_NUCL);
///     SetGapTrigger(BLAST_GAP_TRIGGER_NUCL);
///     SetGapExtnAlgorithm(eDynProgScoreOnly);
///     SetGapTracebackAlgorithm(eDynProgTbck);
/// }
///
/// void
/// CDiscNucleotideOptionsHandle::SetMBScoringOptionsDefaults()
/// {
///     CBlastNucleotideOptionsHandle::SetScoringOptionsDefaults();
/// }
/// ```
/// NCBI reference: c++/include/algo/blast/core/blast_options.h:58-62,67-68,86-87,94-95,130-133,158
/// ```c
/// #define BLAST_WINDOW_SIZE_NUCL 0   /**< default window size (blastn) */
/// #define BLAST_WINDOW_SIZE_MEGABLAST 0   /**< default window size
///                                           (contiguous megablast) */
/// #define BLAST_WINDOW_SIZE_DISC 40  /**< default window size
///                                           (discontiguous megablast) */
/// ...
/// #define BLAST_WORDSIZE_NUCL 11   /**< default word size (blastn) */
/// #define BLAST_WORDSIZE_MEGABLAST 28   /**< default word size (contiguous
/// ...
/// #define BLAST_GAP_OPEN_NUCL 5 /**< default gap open penalty (blastn) */
/// #define BLAST_GAP_OPEN_MEGABLAST 0 /**< default gap open penalty (megablast
/// ...
/// #define BLAST_GAP_EXTN_NUCL 2 /**< default gap open penalty (blastn) */
/// #define BLAST_GAP_EXTN_MEGABLAST 0 /**< default gap open penalty (megablast)
/// ...
/// #define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
///                                          nucleotide gapped extensions */
/// #define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
///                                          nucleotide gapped extensions */
/// ...
/// #define BLAST_EXPECT_VALUE 10.0 /**< by default, alignments whose expect
/// ```
struct TaskDefaults {
    word_size: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    evalue: f64,
    dust: bool,
    greedy: bool,
    x_drop_gapped: i32,
    min_diag_separation: i32,
    window_size: usize,
    template_type: DiscWordType,
    template_length: u8,
}

fn task_defaults(task: &str) -> TaskDefaults {
    match task {
        "megablast" => TaskDefaults {
            word_size: 28,
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
            evalue: 10.0,
            dust: true,
            greedy: true,
            x_drop_gapped: X_DROP_GAPPED_GREEDY,
            min_diag_separation: MIN_DIAG_SEPARATION_MEGABLAST,
            window_size: 0,
            template_type: DiscWordType::Coding,
            template_length: 0,
        },
        "dc-megablast" => TaskDefaults {
            word_size: 11,
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            evalue: 10.0,
            dust: true,
            greedy: false,
            x_drop_gapped: X_DROP_GAPPED_NUCL,
            min_diag_separation: MIN_DIAG_SEPARATION_MEGABLAST,
            window_size: 40,
            template_type: DiscWordType::Coding,
            template_length: 18,
        },
        "blastn-short" => TaskDefaults {
            word_size: 7,
            reward: 1,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            evalue: 1000.0,
            dust: false,
            greedy: false,
            x_drop_gapped: X_DROP_GAPPED_NUCL,
            min_diag_separation: MIN_DIAG_SEPARATION_BLASTN,
            window_size: 0,
            template_type: DiscWordType::Coding,
            template_length: 0,
        },
        _ => TaskDefaults {
            word_size: 11,
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            evalue: 10.0,
            dust: true,
            greedy: false,
            x_drop_gapped: X_DROP_GAPPED_NUCL,
            min_diag_separation: MIN_DIAG_SEPARATION_BLASTN,
            window_size: 0,
            template_type: DiscWordType::Coding,
            template_length: 0,
        },
    }
}

/// The e-value threshold: the given one, or the task's (1000 for blastn-short).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:142-146,253-255
/// ```c
///     	string des = "Expectation value (E) threshold for saving hits. Default = 10";
///     	if(m_IsBlastn) {
///     		des += " (1000 for blastn-short)";
///     	}
///         arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);
/// ...
///     if (args.Exist(kArgEvalue) && args[kArgEvalue]) {
///         opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());
///     }
/// ```
pub fn determine_evalue(args: &BlastnArgs) -> f64 {
    args.evalue
        .unwrap_or_else(|| task_defaults(&args.task).evalue)
}

/// Whether the task filters queries with DUST when `-dust` is not given (not for
/// blastn-short, whose task clears the filtering options).
///
/// NCBI reference: c++/src/algo/blast/api/blast_options_cxx.cpp:1022-1031
/// ```c
/// CBlastOptions::ClearFilterOptions()
/// {
///     SetDustFiltering(false);
///     SetSegFiltering(false);
///     SetRepeatFiltering(false);
///     SetMaskAtHash(false);
///     SetWindowMaskerTaxId(0);
///     SetWindowMaskerDatabase(NULL);
///     return;
/// }
/// ```
pub fn task_dust_by_default(task: &str) -> bool {
    task_defaults(task).dust
}

/// Whether the task's program splits queries in chunks of 5,000,000 (megablast and
/// dc-megablast) rather than 1,000,000 (blastn, blastn-short).
///
/// NCBI reference: c++/src/algo/blast/api/local_blast.cpp:64-72
/// ```c
///         switch (program) {
///         case eBlastn:
///             retval = 1000000;
///             break;
///         case eMegablast:
///         case eDiscMegablast:
///         case eMapper:
///             retval = 5000000;
///             break;
/// ```
pub fn task_uses_megablast_chunks(task: &str) -> bool {
    matches!(task, "megablast" | "dc-megablast")
}

/// The discontiguous template: the given `-template_type` and `-template_length`, which
/// NCBI applies to any task, or the task's.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:735-763
/// ```c
///     if (args[kArgDMBTemplateType]) {
///         const string& type = args[kArgDMBTemplateType].AsString();
///         EDiscWordType temp_type = eMBWordCoding;
///
///         if (type == kTemplType_Coding) {
///             temp_type = eMBWordCoding;
///         } else if (type == kTemplType_Optimal) {
///             temp_type = eMBWordOptimal;
///         } else if (type == kTemplType_CodingAndOptimal) {
///             temp_type = eMBWordTwoTemplates;
///         } else {
///             abort();
///         }
///         options.SetMBTemplateType(static_cast<unsigned char>(temp_type));
///     }
///
///     if (args[kArgDMBTemplateLength]) {
///         unsigned char tlen =
///             static_cast<unsigned char>(args[kArgDMBTemplateLength].AsInteger());
///         options.SetMBTemplateLength(tlen);
///     }
/// ```
pub fn determine_template(args: &BlastnArgs) -> (DiscWordType, u8) {
    let defaults = task_defaults(&args.task);
    (
        args.template_type.unwrap_or(defaults.template_type),
        args.template_length.unwrap_or(defaults.template_length),
    )
}

/// The word size: the given one, or the task's.
pub fn determine_effective_word_size(args: &BlastnArgs) -> usize {
    args.word_size
        .unwrap_or_else(|| task_defaults(&args.task).word_size)
}

/// The reward, penalty and gap costs: each given one, or the task's. NCBI keeps the reward
/// and the penalty in 16 bits, so a given value outside that range wraps before any check
/// or use.
///
/// NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:1628-1643
/// ```c
/// CBlastOptionsLocal::SetMatchReward(int r)
/// {
///     m_ScoringOpts->reward = r;
/// ...
/// CBlastOptionsLocal::SetMismatchPenalty(int p)
/// {
///     m_ScoringOpts->penalty = p;
/// }
/// ```
/// NCBI reference: c++/include/algo/blast/core/blast_options.h:465-466
/// ```c
///    Int2 reward;      /**< Reward for a match */
///    Int2 penalty;     /**< Penalty for a mismatch */
/// ```
pub fn determine_scoring_params(args: &BlastnArgs) -> (i32, i32, i32, i32) {
    let defaults = task_defaults(&args.task);
    (
        i32::from(args.reward.unwrap_or(defaults.reward) as i16),
        i32::from(args.penalty.unwrap_or(defaults.penalty) as i16),
        args.gap_open.unwrap_or(defaults.gap_open),
        args.gap_extend.unwrap_or(defaults.gap_extend),
    )
}

/// Calculate initial scan step based on word size
pub fn calculate_initial_scan_step(effective_word_size: usize, user_scan_step: usize) -> usize {
    if user_scan_step > 0 {
        user_scan_step
    } else {
        if effective_word_size >= 16 {
            4
        } else if effective_word_size >= 11 {
            2
        } else {
            1
        }
    }
}

/// Configure task-specific parameters
pub fn configure_task(args: &BlastnArgs) -> TaskConfig {
    let effective_word_size = determine_effective_word_size(args);
    let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(args);

    let min_ungapped_score = match args.task.as_str() {
        "megablast" => MIN_UNGAPPED_SCORE_MEGABLAST,
        _ => MIN_UNGAPPED_SCORE_BLASTN,
    };

    let defaults = task_defaults(&args.task);
    // NCBI BLAST algorithm selection (`task_defaults`):
    // - megablast: eGreedyScoreOnly (greedy alignment)
    // - blastn, blastn-short, dc-megablast: eDynProgScoreOnly (dynamic programming)
    // Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:182, 192;
    // disc_nucl_options.cpp:76-84
    let use_dp = !defaults.greedy;

    // NCBI BLAST: Task-specific gapped X-dropoff (`task_defaults`)
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122-148
    // DP (blastn, blastn-short, dc-megablast): 30, greedy (megablast): 25
    let x_drop_gapped = defaults.x_drop_gapped;

    // NCBI BLAST: Final traceback X-dropoff (100 for all nucleotide tasks)
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:146
    // BLAST_GAP_X_DROPOFF_FINAL_NUCL = 100
    // Used in traceback phase to extend alignments further than preliminary extension
    let x_drop_final = X_DROP_GAPPED_FINAL; // 100

    let scan_step = calculate_initial_scan_step(effective_word_size, args.scan_step);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:163-168
    // ```c
    // SetWindowSize(BLAST_WINDOW_SIZE_NUCL);
    // SetOffDiagonalRange(BLAST_SCAN_RANGE_NUCL);
    // ```
    let scan_range = match args.task.as_str() {
        "megablast" => SCAN_RANGE_MEGABLAST, // 0
        _ => SCAN_RANGE_BLASTN,              // 0
    };

    // NCBI reference: blast_nucl_options.cpp:239, 259
    // Minimum diagonal separation for HSP containment checking
    // Used in MB_HSP_CLOSE macro (blast_gapalign_priv.h:123-124)
    // megablast and dc-megablast: 6; blastn and blastn-short: 50 (`task_defaults`).
    let min_diag_separation = defaults.min_diag_separation;
    let (mb_template_type, mb_template_length) = determine_template(args);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
    // ```c
    // lut_type = BlastChooseNaLookupTable(lookup_options,
    //                                     approx_table_entries,
    //                                     max_q_off,
    //                                     &lut_width);
    // ```
    let use_two_stage = false;
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;
    let lut_word_length = effective_word_size;

    TaskConfig {
        effective_word_size,
        reward,
        penalty,
        gap_open,
        gap_extend,
        min_ungapped_score,
        use_dp,
        scan_step,
        use_direct_lookup,
        use_two_stage,
        lut_word_length,
        small_na_lookup: false,
        mb_lookup: false,
        x_drop_gapped,
        x_drop_final,
        scan_range,
        min_diag_separation,
        window_size: defaults.window_size,
        mb_template_type,
        mb_template_length,
    }
}

/// Choose lookup table kind and LUT width using NCBI's BlastChooseNaLookupTable logic.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
/// ```c
/// if (lookup_options->mb_template_length > 0) {
///    *lut_width = lookup_options->word_size;
///    return eMBLookupTable;
/// }
/// switch(lookup_options->word_size) {
/// case 4: case 5: case 6:
///    lut_type = eSmallNaLookupTable;
///    *lut_width = lookup_options->word_size;
///    break;
/// case 7:
///    lut_type = eSmallNaLookupTable;
///    if (approx_table_entries < 250) *lut_width = 6;
///    else *lut_width = 7;
///    break;
/// case 8:
///    lut_type = eSmallNaLookupTable;
///    if (approx_table_entries < 8500) *lut_width = 7;
///    else *lut_width = 8;
///    break;
/// case 9:
///    if (approx_table_entries < 1250) { *lut_width = 7; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 21000) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else { *lut_width = 9; lut_type = eMBLookupTable; }
///    break;
/// case 10:
///    if (approx_table_entries < 1250) { *lut_width = 7; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 18000) { *lut_width = 9; lut_type = eMBLookupTable; }
///    else { *lut_width = 10; lut_type = eMBLookupTable; }
///    break;
/// case 11:
///    if (approx_table_entries < 12000) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 180000) { *lut_width = 10; lut_type = eMBLookupTable; }
///    else { *lut_width = 11; lut_type = eMBLookupTable; }
///    break;
/// case 12:
///    if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 18000) { *lut_width = 9; lut_type = eMBLookupTable; }
///    else if (approx_table_entries < 60000) { *lut_width = 10; lut_type = eMBLookupTable; }
///    else if (approx_table_entries < 900000) { *lut_width = 11; lut_type = eMBLookupTable; }
///    else { *lut_width = 12; lut_type = eMBLookupTable; }
///    break;
/// default:
///    if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 300000) { *lut_width = 11; lut_type = eMBLookupTable; }
///    else { *lut_width = 12; lut_type = eMBLookupTable; }
///    break;
/// }
/// if (lut_type == eSmallNaLookupTable &&
///     (approx_table_entries >= 32767 || max_q_off >= 32768)) {
///    lut_type = eNaLookupTable;
/// }
/// ```
fn choose_na_lookup_table(
    word_size: usize,
    approx_table_entries: usize,
    max_q_off: usize,
    discontig_template: bool,
) -> (LookupTableKind, usize) {
    debug_assert!(word_size >= 4);

    if discontig_template {
        return (LookupTableKind::Mb, word_size);
    }

    let (mut lut_kind, mut lut_width) = match word_size {
        4 | 5 | 6 => (LookupTableKind::Small, word_size),
        7 => {
            if approx_table_entries < 250 {
                (LookupTableKind::Small, 6)
            } else {
                (LookupTableKind::Small, 7)
            }
        }
        8 => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 7)
            } else {
                (LookupTableKind::Small, 8)
            }
        }
        9 => {
            if approx_table_entries < 1_250 {
                (LookupTableKind::Small, 7)
            } else if approx_table_entries < 21_000 {
                (LookupTableKind::Small, 8)
            } else {
                (LookupTableKind::Mb, 9)
            }
        }
        10 => {
            if approx_table_entries < 1_250 {
                (LookupTableKind::Small, 7)
            } else if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 18_000 {
                (LookupTableKind::Mb, 9)
            } else {
                (LookupTableKind::Mb, 10)
            }
        }
        11 => {
            if approx_table_entries < LUT_WIDTH_11_THRESHOLD_8 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < LUT_WIDTH_11_THRESHOLD_10 {
                (LookupTableKind::Mb, 10)
            } else {
                (LookupTableKind::Mb, 11)
            }
        }
        12 => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 18_000 {
                (LookupTableKind::Mb, 9)
            } else if approx_table_entries < 60_000 {
                (LookupTableKind::Mb, 10)
            } else if approx_table_entries < 900_000 {
                (LookupTableKind::Mb, 11)
            } else {
                (LookupTableKind::Mb, 12)
            }
        }
        _ => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 300_000 {
                (LookupTableKind::Mb, 11)
            } else {
                (LookupTableKind::Mb, 12)
            }
        }
    };

    if lut_kind == LookupTableKind::Small && (approx_table_entries >= 32_767 || max_q_off >= 32_768)
    {
        lut_kind = LookupTableKind::Na;
    }

    (lut_kind, lut_width)
}

/// Finalize task configuration with query-dependent parameters.
/// Must be called after queries are loaded to enable adaptive lookup table selection.
///
/// NCBI reference: blast_nalookup.c (BlastChooseNaLookupTable)
/// - Uses total query length (approx_table_entries) and max query offset
/// - Handles discontiguous template forcing MB lookup
pub fn finalize_task_config(
    config: &mut TaskConfig,
    total_query_length: usize,
    max_query_length: usize,
    discontig_template: bool,
) {
    let max_q_off = max_query_length.saturating_sub(1);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
    // ```c
    // lut_type = BlastChooseNaLookupTable(lookup_options,
    //                                     approx_table_entries,
    //                                     max_q_off,
    //                                     &lut_width);
    // ```
    let (lut_kind, lut_width) = choose_na_lookup_table(
        config.effective_word_size,
        total_query_length,
        max_q_off,
        discontig_template,
    );

    config.lut_word_length = lut_width;
    config.small_na_lookup = lut_kind == LookupTableKind::Small;
    config.mb_lookup = lut_kind == LookupTableKind::Mb;
    config.use_two_stage =
        lut_kind == LookupTableKind::Mb || config.lut_word_length < config.effective_word_size;
    config.use_direct_lookup =
        !config.use_two_stage && config.lut_word_length <= MAX_DIRECT_LOOKUP_WORD_SIZE;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:397,1334
    // ```c
    // lookup->scan_step = lookup->word_length - lookup->lut_word_length + 1;
    // mb_lt->scan_step = mb_lt->word_length - mb_lt->lut_word_length + 1;
    // ```
    config.scan_step = (config.effective_word_size - config.lut_word_length + 1).max(1);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1413
// ```c
// db_length = BlastSeqSrcGetTotLen(seq_src);
// itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
// ```
pub fn subject_metadata_from_records<R: InputRecord>(records: &[R]) -> SubjectMetadata {
    let mut subject_ids = Vec::with_capacity(records.len());
    let mut db_len_total: usize = 0;
    for record in records {
        // The title up to its first space: the ID of a `bio` record (whose ID has no white
        // space), "unknown" when it is empty.
        let title = record.title_bytes();
        let end = title
            .iter()
            .position(|&byte| byte == b' ')
            .unwrap_or(title.len());
        subject_ids.push(if end == 0 {
            "unknown".to_string()
        } else {
            String::from_utf8_lossy(&title[..end]).into_owned()
        });
        db_len_total += record.seq().len();
    }
    SubjectMetadata {
        db_len_total,
        db_num_seqs: records.len(),
        subject_ids,
    }
}

// NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:856-874
// ```c
// case 'a': case 'b': case 'c': case 'd':
// case 'g': case 'h':
// ...
//     char_type = eCharType_MaskedNonGap;
//     break;
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:1079-1089
// ```c
// void CFastaReader::x_OpenMask(void)
// {
//     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
// }
// void CFastaReader::x_CloseMask(void)
// {
//     m_CurrentMask->SetPacked_int().AddInterval(...);
// }
// ```
pub fn collect_lowercase_masks(seq: &[u8]) -> Vec<MaskedInterval> {
    let mut masks = Vec::new();
    let mut current_start: Option<usize> = None;

    for (idx, &base) in seq.iter().enumerate() {
        let is_lower = base.is_ascii_lowercase();
        match (current_start, is_lower) {
            (None, true) => current_start = Some(idx),
            (Some(start), false) => {
                if start < idx {
                    masks.push(MaskedInterval::new(start, idx));
                }
                current_start = None;
            }
            _ => {}
        }
    }

    if let Some(start) = current_start {
        if start < seq.len() {
            masks.push(MaskedInterval::new(start, seq.len()));
        }
    }

    masks
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
// ```c
// const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking) ? ... : kDfltArgUseLCaseMasking;
// ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
// ```
fn collect_lowercase_masks_for_records<R: InputRecord>(records: &[R]) -> Vec<Vec<MaskedInterval>> {
    records
        .iter()
        .map(|record| collect_lowercase_masks(record.seq()))
        .collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-128
// ```c
// const int kTopFlags = CSeq_loc::fStrand_Ignore|CSeq_loc::fMerge_All|CSeq_loc::fSort;
// if (orig_query_mask.NotEmpty()) {
//     orig_query_mask->Add(*query_masks,  kTopFlags, 0);
// } else {
//     query_masks->Merge(kTopFlags, 0);
// }
// ```
fn merge_mask_intervals(mut intervals: Vec<MaskedInterval>) -> Vec<MaskedInterval> {
    if intervals.is_empty() {
        return intervals;
    }

    intervals.sort_by_key(|interval| interval.start);
    let mut merged: Vec<MaskedInterval> = Vec::with_capacity(intervals.len());
    let mut current = intervals[0].clone();

    for interval in intervals.into_iter().skip(1) {
        if interval.start <= current.end {
            current.end = current.end.max(interval.end);
        } else {
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);
    merged
}

/// Apply DUST filter to query sequences
pub fn apply_dust_masking<R: InputRecord>(
    args: &BlastnArgs,
    queries: &[R],
) -> Vec<Vec<MaskedInterval>> {
    // NCBI blast_args.cpp:418-420: opt.SetDustFilteringLevel(...);
    // opt.SetDustFilteringWindow(...); opt.SetDustFilteringLinker(...);
    if let Some((dust_level, dust_window, dust_linker)) = args.dust.params() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-128
        // ```c
        // orig_query_mask->Add(*query_masks,  kTopFlags, 0);
        // query_masks->Merge(kTopFlags, 0);
        // ```
        if args.verbose {
            eprintln!(
                "Applying DUST filter (level={}, window={}, linker={})...",
                dust_level, dust_window, dust_linker
            );
        }
        let masks: Vec<Vec<MaskedInterval>> = queries
            .iter()
            .map(|record| {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-102
                // ```c
                // void s_CombineDustMasksWithUserProvidedMasks(...)
                // {
                //     CSymDustMasker duster(level, window, linker);
                //     CRef<CPacked_seqint> masked_locations =
                //         duster.GetMaskedInts(*query_id, data);
                // ```
                // NCBI constructs a fresh CSymDustMasker for each query, which
                // resets the converter's CRandom state before scanning that query.
                DustMasker::new(dust_level, dust_window, dust_linker).mask_sequence(record.seq())
            })
            .collect();

        let total_masked: usize = masks
            .iter()
            .map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>())
            .sum();
        let total_bases: usize = queries.iter().map(|r| r.seq().len()).sum();
        if args.verbose && total_bases > 0 {
            eprintln!(
                "DUST masked {} bases ({:.2}%) across {} sequences",
                total_masked,
                100.0 * total_masked as f64 / total_bases as f64,
                queries.len()
            );
        }
        masks
    } else {
        vec![Vec::new(); queries.len()]
    }
}

/// Build lookup tables based on configuration
pub fn build_lookup_tables<R: InputRecord>(
    config: &TaskConfig,
    args: &BlastnArgs,
    queries_blastna: &[Vec<u8>],
    query_masks: &[Vec<MaskedInterval>],
    query_offsets: &[i32],
    subjects: &[R],
    subjects_packed: Option<&[Vec<u8>]>,
    approx_table_entries: usize,
) -> (LookupTables, usize) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1231-1344
    // ```c
    // Int4 BlastMBLookupTableNew(...);
    // Int4 BlastNaLookupTableNew(...);
    // ```
    if args.verbose {
        eprintln!(
            "Building lookup (Task: {}, Word: {}, TwoStage: {}, LUTWord: {}, Direct: {}, DUST: {})...",
            args.task,
            config.effective_word_size,
            config.use_two_stage,
            config.lut_word_length,
            config.use_direct_lookup,
            args.dust.is_enabled()
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // ```c
    // if (lookup_options->db_filter) {
    //    counts = (Uint1*)calloc(mb_lt->hashsize / 2, sizeof(Uint1));
    // }
    // if (lookup_options->db_filter) {
    //    s_FillPV(query, location, mb_lt, lookup_options);
    //    s_ScanSubjectForWordCounts(seqsrc, mb_lt, counts,
    //                               lookup_options->max_db_word_count);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // ```c
    // pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
    // mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
    // ```
    let db_word_counts = if args.limit_lookup {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // if (subj_is_na) {
        //     BlastSeqBlkSetSequence(subj, sequence.data.release(),
        //        ((sentinels == eSentinels) ? sequence.length - 2 :
        //         sequence.length));
        //     ...
        //     SBlastSequence compressed_seq =
        //         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
        //                                   eNa_strand_plus, eNoSentinels);
        //     BlastSeqBlkSetCompressedSequence(subj,
        //                               compressed_seq.data.release());
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:865-868
        // ```c
        // if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
        //     continue;
        // ```
        Some(build_db_word_counts(
            queries_blastna,
            query_masks,
            subjects,
            config.lut_word_length,
            config.effective_word_size,
            args.max_db_word_count,
            approx_table_entries,
            subjects_packed,
        ))
    } else {
        None
    };
    let db_word_counts_ref = db_word_counts.as_deref();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
    // ```c
    // /* Also add 1 to all indices, because lookup table indices count
    //    from 1. */
    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    // mb_lt->hashtable[ecode] = index;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1327-1341
    // ```c
    //    if (lookup_options->mb_template_length > 0) {
    //         /* discontiguous megablast */
    //         mb_lt->scan_step = 1;
    //         status = s_FillDiscMBTable(query, location, mb_lt, lookup_options);
    //    }
    //    else {
    //         /* contiguous megablast */
    //         mb_lt->scan_step = mb_lt->word_length - mb_lt->lut_word_length + 1;
    //         status = s_FillContigMBTable(query, location, mb_lt, lookup_options,
    //                                      counts);
    // ```
    let two_stage_lookup: Option<TwoStageLookup> = if config.mb_template_length > 0 {
        Some(build_disc_mb_lookup(
            queries_blastna,
            query_offsets,
            config.effective_word_size,
            config.mb_template_length as usize,
            config.mb_template_type,
            query_masks,
            approx_table_entries,
        ))
    } else if config.use_two_stage {
        Some(build_two_stage_lookup(
            queries_blastna,
            query_offsets,
            config.effective_word_size,
            config.lut_word_length,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
            approx_table_entries,
            !config.mb_lookup,
        ))
    } else {
        None
    };

    let pv_direct_lookup: Option<PvDirectLookup> =
        if !config.use_two_stage && config.use_direct_lookup {
            Some(build_pv_direct_lookup(
                queries_blastna,
                query_offsets,
                config.effective_word_size,
                config.effective_word_size,
                query_masks,
                db_word_counts_ref,
                args.max_db_word_count,
                approx_table_entries,
                false,
            ))
        } else {
            None
        };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:548-583
    // ```c
    // Int4 BlastNaLookupTableNew(BLAST_SequenceBlk* query,
    //                            BlastSeqLoc* locations,
    //                            BlastNaLookupTable * *lut,
    //                            const LookupTableOptions * opt,
    //                            const QuerySetUpOptions* query_options,
    //                            Int4 lut_width)
    // ```
    let na_lookup: Option<NaLookupTable> = if !config.use_two_stage && !config.use_direct_lookup {
        Some(build_na_lookup(
            queries_blastna,
            query_offsets,
            config.effective_word_size,
            config.lut_word_length,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
        ))
    } else {
        None
    };

    // Use scan_step from config - if user specified a value > 0, use it; otherwise use the configured default
    // For two-stage lookup, the default was already calculated in configure_task
    let scan_step = config.scan_step;

    if args.verbose && config.use_two_stage {
        eprintln!(
            "[INFO] Using two-stage lookup: lut_word_length={}, word_length={}, scan_step={}",
            config.lut_word_length, config.effective_word_size, scan_step
        );
    }

    (
        LookupTables {
            two_stage_lookup,
            pv_direct_lookup,
            na_lookup,
        },
        scan_step,
    )
}

/// The masks of the queries: their DUST masks and, with `-lcase_masking`, their lower-case
/// letters.
pub fn query_masks<R: InputRecord>(args: &BlastnArgs, queries: &[R]) -> Vec<Vec<MaskedInterval>> {
    let mut query_masks = apply_dust_masking(args, queries);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
    // ```c
    // const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking) ? ... : kDfltArgUseLCaseMasking;
    // ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-128
    // ```c
    // orig_query_mask->Add(*query_masks,  kTopFlags, 0);
    // query_masks->Merge(kTopFlags, 0);
    // ```
    if args.lcase_masking {
        let lcase_masks = collect_lowercase_masks_for_records(queries);
        for (dust_masks, lcase) in query_masks.iter_mut().zip(lcase_masks) {
            if !lcase.is_empty() {
                dust_masks.extend(lcase);
                let merged = merge_mask_intervals(std::mem::take(dust_masks));
                *dust_masks = merged;
            }
        }
    }
    query_masks
}

/// The masks of the query parts of a query chunk: the DUST masks of each part added to the
/// masks that NCBI keeps from its query (`query_split::restrict_masks`), which already have
/// the lower-case letters.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:171-180
/// ```c
///         CSeqVector data(*queries.GetQuerySeqLoc(i), *queries.GetScope(i),
///                         CBioseq_Handle::eCoding_Iupac);
///         ...
///         CRef<CSeq_loc> masks = queries.GetMasks(i);
///         s_CombineDustMasksWithUserProvidedMasks(data,
///                                                 queries.GetQuerySeqLoc(i),
///                                                 queries.GetScope(i), query_id,
///                                                 masks, level, window, linker);
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:121-128
/// ```c
///     const int kTopFlags = CSeq_loc::fStrand_Ignore|CSeq_loc::fMerge_All|CSeq_loc::fSort;
///     if (orig_query_mask.NotEmpty() && !orig_query_mask->IsNull()) {
///         CRef<CSeq_loc> tmp = orig_query_mask->Add(*query_masks,  kTopFlags, 0);
///         orig_query_mask.Reset(tmp);
///     } else {
///         query_masks->Merge(kTopFlags, 0);
///         orig_query_mask.Reset(query_masks);
///     }
/// ```
/// The part's own masks are only merged when DUST finds a region; the lookup table and the
/// search use the masked residues, which merging does not change.
pub fn chunk_query_masks<R: InputRecord>(
    args: &BlastnArgs,
    parts: &[R],
    restricted: Vec<Vec<MaskedInterval>>,
) -> Vec<Vec<MaskedInterval>> {
    apply_dust_masking(args, parts)
        .into_iter()
        .zip(restricted)
        .map(|(mut masks, restricted)| {
            masks.extend(restricted);
            merge_mask_intervals(masks)
        })
        .collect()
}

/// Prepare all sequence data and configuration
pub fn prepare_sequence_data<R>(
    queries: Vec<R>,
    query_ids: Vec<String>,
    query_masks: Vec<Vec<MaskedInterval>>,
    subject_metadata: SubjectMetadata,
) -> SequenceData<R> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    let db_len_total = subject_metadata.db_len_total;
    let db_num_seqs = subject_metadata.db_num_seqs;
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    let subject_ids = subject_metadata.subject_ids;

    SequenceData {
        queries,
        query_ids,
        query_masks,
        db_len_total,
        db_num_seqs,
        subject_ids,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:867-949
    // ```c
    // case 'a': case 'b': case 'c': case 'd': ...
    //     char_type = eCharType_MaskedNonGap;
    // ...
    // case eCharType_MaskedNonGap:
    //     OpenMask();
    // ```
    fn blastn_args(words: &[&str]) -> BlastnArgs {
        let argv = [
            "LOSAT", "blastn", "-query", "q", "-subject", "s", "-outfmt", "6",
        ];
        let cli: crate::cli::Cli =
            crate::cli::try_parse_from(argv.iter().chain(words)).expect("valid arguments");
        let crate::cli::Commands::Blastn(args) = cli.command else {
            panic!("blastn")
        };
        args
    }

    // NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:344-380
    // ```c
    //         if (!NStr::CompareNocase(task, "blastn-short"))
    //         {
    //              opts->SetMatchReward(1);
    //              opts->SetMismatchPenalty(-3);
    //              opts->SetEvalueThreshold(1000);
    //              opts->SetWordSize(7);
    //              opts->ClearFilterOptions();
    //         }
    // ```
    // NCBI reference: c++/src/algo/blast/api/disc_nucl_options.cpp:55-84
    // ```c
    //     SetTemplateType(0);
    //     SetTemplateLength(18);
    //     SetWordSize(BLAST_WORDSIZE_NUCL);
    // ...
    //     SetWindowSize(BLAST_WINDOW_SIZE_DISC);
    // ...
    //     SetGapExtnAlgorithm(eDynProgScoreOnly);
    // ```
    #[test]
    fn task_defaults_follow_ncbi_handles() {
        let config = configure_task(&blastn_args(&["-task", "dc-megablast"]));
        assert_eq!(
            (config.effective_word_size, config.reward, config.penalty),
            (11, 2, -3)
        );
        assert_eq!((config.gap_open, config.gap_extend), (5, 2));
        assert!(config.use_dp);
        assert_eq!((config.x_drop_gapped, config.x_drop_final), (30, 100));
        assert_eq!(config.min_diag_separation, 6);
        assert_eq!(config.window_size, 40);
        assert_eq!(
            (config.mb_template_type, config.mb_template_length),
            (DiscWordType::Coding, 18)
        );

        let args = blastn_args(&["-task", "blastn-short"]);
        let config = configure_task(&args);
        assert_eq!(
            (config.effective_word_size, config.reward, config.penalty),
            (7, 1, -3)
        );
        assert_eq!((config.gap_open, config.gap_extend), (5, 2));
        assert!(config.use_dp);
        assert_eq!(config.min_diag_separation, 50);
        assert_eq!((config.window_size, config.mb_template_length), (0, 0));
        assert_eq!(determine_evalue(&args), 1000.0);
        assert!(!task_dust_by_default("blastn-short"));
        assert!(!task_uses_megablast_chunks("blastn-short"));
        assert!(task_uses_megablast_chunks("dc-megablast"));

        // Given options replace the task's values; the template applies to any task.
        let args = blastn_args(&[
            "-task",
            "megablast",
            "-word_size",
            "12",
            "-evalue",
            "5",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "21",
        ]);
        let config = configure_task(&args);
        assert_eq!(config.effective_word_size, 12);
        assert!(!config.use_dp);
        assert_eq!(config.window_size, 0);
        assert_eq!(
            (config.mb_template_type, config.mb_template_length),
            (DiscWordType::TwoTemplates, 21)
        );
        assert_eq!(determine_evalue(&args), 5.0);
        for task in ["megablast", "blastn", "dc-megablast"] {
            assert!(task_dust_by_default(task));
            assert_eq!(determine_evalue(&blastn_args(&["-task", task])), 10.0);
        }
    }

    #[test]
    fn test_collect_lowercase_masks() {
        let masks = collect_lowercase_masks(b"AAaaBBbC");
        assert_eq!(
            masks,
            vec![MaskedInterval::new(2, 4), MaskedInterval::new(6, 7)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:121-127
    // ```c
    // const int kTopFlags = ... | CSeq_loc::fMerge_All | ...;
    // query_masks->Merge(kTopFlags, 0);
    // ```
    #[test]
    fn test_merge_mask_intervals() {
        let merged = merge_mask_intervals(vec![
            MaskedInterval::new(0, 2),
            MaskedInterval::new(2, 5),
            MaskedInterval::new(6, 7),
        ]);
        assert_eq!(
            merged,
            vec![MaskedInterval::new(0, 5), MaskedInterval::new(6, 7)]
        );
    }
}
