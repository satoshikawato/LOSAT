//! Coordination module for BLASTN pipeline setup
//!
//! This module handles:
//! - Task-specific parameter configuration (word size, scoring, etc.)
//! - Sequence reading and preprocessing
//! - DUST masking
//! - Lookup table building

use super::args::{BlastnArgs, QueryStrand};
use super::constants::{
    BLAST_PENALTY, BLAST_REWARD, LUT_WIDTH_11_THRESHOLD_10, LUT_WIDTH_11_THRESHOLD_8,
    MAX_DIRECT_LOOKUP_WORD_SIZE, MIN_DIAG_SEPARATION_BLASTN, MIN_DIAG_SEPARATION_MEGABLAST,
    MIN_UNGAPPED_SCORE_BLASTN, MIN_UNGAPPED_SCORE_MEGABLAST, SCAN_RANGE_BLASTN,
    SCAN_RANGE_MEGABLAST, X_DROP_GAPPED_FINAL, X_DROP_GAPPED_GREEDY, X_DROP_GAPPED_NUCL,
    X_DROP_UNGAPPED,
};
use super::disc_lookup::DiscWordType;
use super::lookup::{
    build_db_word_counts, build_disc_mb_lookup_kept, build_na_lookup, build_pv_direct_lookup,
    build_two_stage_lookup_kept, x_mb_lut_keep, NaLookupTable, PvDirectLookup, TwoStageLookup,
};
use crate::blastinput::fasta_reader::{FastaRecord, InputRecord};
use crate::utils::dust::MaskedInterval;
use crate::utils::threading::SearchPool;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LookupTableKind {
    Small,
    Mb,
    Na,
}

/// Task-specific configuration for BLASTN
pub struct TaskConfig {
    pub effective_word_size: usize,
    /// The option pair (`-reward`/`-penalty`, else the task's): the reports and the
    /// endpoint re-evaluation of the traceback read it.
    pub reward: i32,
    pub penalty: i32,
    /// NCBI's matrix-only scoring: the option pair is 0/0 (`-task rmblastn`).
    pub matrix_only: bool,
    /// The engine pair: the matrix, the score table and the Karlin-Altschul statistics
    /// score with it (1/-3 under matrix-only scoring, else the option pair).
    pub score_reward: i32,
    pub score_penalty: i32,
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
    /// The gapped X-dropoff in bits (`-xdrop_gap`, else the task's: blastn 30, megablast 25),
    /// converted to a raw score once per batch (`scoring.rs` `gap_x_dropoffs`).
    pub x_drop_gapped: f64,
    /// The final (traceback) X-dropoff in bits (`-xdrop_gap_final`, else 100 for every task).
    pub x_drop_final: f64,
    /// The ungapped X-dropoff in bits (`-xdrop_ungap`, else 20; 0 for megablast, whose raw
    /// X-drop is then the word cutoff of the context, blast_parameters.c:378-381).
    pub x_drop_ungapped: f64,
    /// The raw hit-saving cutoff that replaces the e-value's (`-min_raw_gapped_score` above 0).
    pub min_raw_gapped_score: Option<i32>,
    /// The off-diagonal range of the second hit (`-off_diagonal_range`, else 0).
    pub scan_range: usize,
    pub min_diag_separation: i32, // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:239,259 (blastn: 50, megablast: 6)
    /// The two-hit window (`-window_size`, else the task's: 40 for dc-megablast, 0 (one hit)
    /// otherwise).
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:131-156
    pub na_lookup: Option<NaLookupTable>,
}

/// Subject metadata needed for search setup.
#[derive(Clone)]
pub struct SubjectMetadata {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1409
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-906
    /// The statistics' total length, NCBI's `Int8` (`DbStatistics::setup_total` of a database,
    /// the sum of a subject set's records): an alias file's `LENGTH` declares it, so it is not
    /// bounded by the letters LOSAT holds (2^32 or more on a 32-bit target).
    pub db_len_total: i64,
    /// The statistics' number of sequences, NCBI's `Int4` (`DbStatistics::setup_count`: negative
    /// for an alias file's negative `NSEQ`, or one of 2^31 or more).
    pub db_num_seqs: i32,
}

/// Sequence data and metadata (`R`: the query records, `InputRecord`)
pub struct SequenceData<'a, R = FastaRecord> {
    pub queries: &'a [R],
    pub query_masks: Vec<Vec<MaskedInterval>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1409
    /// NCBI's `Int8` (`SubjectMetadata::db_len_total`).
    pub db_len_total: i64,
    /// The statistics' number of sequences, NCBI's `Int4` (`DbStatistics::setup_count`: negative
    /// for an alias file's negative `NSEQ`, or one of 2^31 or more).
    pub db_num_seqs: i32,
}

/// The option values of a task, which the command line options replace (NCBI
/// `CBlastOptionsFactory::CreateTask`). The CLI accepts the tasks megablast, blastn,
/// dc-megablast and blastn-short.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:344-380
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:137-150,153-170,173-194,198-221,223-262
/// dc-megablast is the megablast handle with the discontiguous overrides; the base
/// constructor's calls bind to the base class, and the derived constructor's
/// `SetDefaults` then calls the overrides (the hit-saving defaults are not overridden).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:47-90
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:58-62,67-68,86-87,94-95,130-133,158
struct TaskDefaults {
    word_size: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    evalue: f64,
    dust: bool,
    greedy: bool,
    x_drop_gapped: f64,
    x_drop_ungapped: f64,
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
            x_drop_gapped: f64::from(X_DROP_GAPPED_GREEDY),
            // `SetMBInitialWordOptionsDefaults` sets the window only; the calloc'd X-drop
            // (api/blast_options_local_priv.cpp:54) stays 0.
            x_drop_ungapped: 0.0,
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
            x_drop_gapped: f64::from(X_DROP_GAPPED_NUCL),
            x_drop_ungapped: f64::from(X_DROP_UNGAPPED),
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
            x_drop_gapped: f64::from(X_DROP_GAPPED_NUCL),
            x_drop_ungapped: f64::from(X_DROP_UNGAPPED),
            min_diag_separation: MIN_DIAG_SEPARATION_BLASTN,
            window_size: 0,
            template_type: DiscWordType::Coding,
            template_length: 0,
        },
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:365-370
        // The blastn handle with reward 0 and penalty 0 (matrix-only scoring).
        "rmblastn" => TaskDefaults {
            reward: 0,
            penalty: 0,
            ..task_defaults("blastn")
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
            x_drop_gapped: f64::from(X_DROP_GAPPED_NUCL),
            x_drop_ungapped: f64::from(X_DROP_UNGAPPED),
            min_diag_separation: MIN_DIAG_SEPARATION_BLASTN,
            window_size: 0,
            template_type: DiscWordType::Coding,
            template_length: 0,
        },
    }
}

/// The e-value threshold: the given one, or the task's (1000 for blastn-short).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:142-146,253-255
pub fn determine_evalue(args: &BlastnArgs) -> f64 {
    args.evalue
        .unwrap_or_else(|| task_defaults(&args.task).evalue)
}

/// Whether the task filters queries with DUST when `-dust` is not given (not for
/// blastn-short, whose task clears the filtering options).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1022-1031
pub fn task_dust_by_default(task: &str) -> bool {
    task_defaults(task).dust
}

/// Whether the task's program splits queries in chunks of 5,000,000 (megablast and
/// dc-megablast) rather than 1,000,000 (blastn, blastn-short).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:64-72
pub fn task_uses_megablast_chunks(task: &str) -> bool {
    matches!(task, "megablast" | "dc-megablast")
}

/// The discontiguous template: the given `-template_type` and `-template_length`, which
/// NCBI applies to any task, or the task's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:735-763
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1628-1643
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:465-466
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
    } else if effective_word_size >= 16 {
        4
    } else if effective_word_size >= 11 {
        2
    } else {
        1
    }
}

/// The two-hit window of a run: `-window_size` when given, else the task's (40 for
/// dc-megablast, 0 otherwise; `task_defaults`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:487-497
/// BLASTN has no matrix, so no suggested window replaces the task's.
pub fn determine_window_size(args: &BlastnArgs) -> usize {
    args.window_size
        .map_or(task_defaults(&args.task).window_size, |window| {
            window as usize
        })
}

/// The matrix-only flag and the engine pair of an option pair (reward, penalty).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:359-368
/// (`BLAST_REWARD`/`BLAST_PENALTY`: constants.rs.) The Karlin-Altschul blocks (blast_setup.c:90-96) and alpha/beta (blast_setup.c:796-799)
/// take the same 1/-3 under matrix-only scoring.
pub fn engine_scoring_pair(reward: i32, penalty: i32) -> (bool, i32, i32) {
    if matrix_only_scoring(reward, penalty) {
        (true, BLAST_REWARD, BLAST_PENALTY)
    } else {
        (false, reward, penalty)
    }
}

/// Whether an option pair (reward, penalty), as NCBI holds them (16-bit), is matrix-only scoring:
/// the one test of the concept, for the search (`engine_scoring_pair`), the option checks and the
/// report's introduction, which NCBI's formatter decides from the same two option values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:168-170
pub fn matrix_only_scoring(reward: i32, penalty: i32) -> bool {
    penalty == 0 && reward == 0
}

/// Configure task-specific parameters
pub fn configure_task(args: &BlastnArgs) -> TaskConfig {
    let effective_word_size = determine_effective_word_size(args);
    let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(args);
    let (matrix_only, score_reward, score_penalty) = engine_scoring_pair(reward, penalty);

    let min_ungapped_score = match args.task.as_str() {
        "megablast" => MIN_UNGAPPED_SCORE_MEGABLAST,
        _ => MIN_UNGAPPED_SCORE_BLASTN,
    };

    let defaults = task_defaults(&args.task);
    // NCBI BLAST algorithm selection (`task_defaults`):
    // - megablast: eGreedyScoreOnly (greedy alignment)
    // - blastn, blastn-short, dc-megablast: eDynProgScoreOnly (dynamic programming)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:182, 192;
    // disc_nucl_options.cpp:76-84
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:681-685
    let use_dp = !defaults.greedy || args.no_greedy;

    // NCBI BLAST: Task-specific gapped X-dropoff (`task_defaults`), in bits
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:122-148
    // DP (blastn, blastn-short, dc-megablast): 30, greedy (megablast): 25; `-xdrop_gap`
    // replaces it (blast_args.cpp:280-282).
    let x_drop_gapped = args.xdrop_gap.unwrap_or(defaults.x_drop_gapped);

    // NCBI BLAST: Final traceback X-dropoff (100 for all nucleotide tasks), in bits
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:146
    // Default `BLAST_GAP_X_DROPOFF_FINAL_NUCL` is 100; `-xdrop_gap_final` replaces it
    // (blast_args.cpp:284-286). The raw value is at least the gapped one (`gap_x_dropoffs`).
    let x_drop_final = args
        .xdrop_gap_final
        .unwrap_or(f64::from(X_DROP_GAPPED_FINAL));

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:276-278
    let x_drop_ungapped = args.xdrop_ungap.unwrap_or(defaults.x_drop_ungapped);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:876-877
    // 0 and below are the options' unset value (the handles' `SetCutoffScore` of 0,
    // blast_nucl_options.cpp:243,263).
    let min_raw_gapped_score = args.min_raw_gapped_score.filter(|&score| score > 0);

    let scan_step = calculate_initial_scan_step(effective_word_size, args.scan_step);

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:163-168
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:518-522
    let scan_range = args.off_diagonal_range.map_or(
        match args.task.as_str() {
            "megablast" => SCAN_RANGE_MEGABLAST, // 0
            _ => SCAN_RANGE_BLASTN,              // 0
        },
        |range| range as usize,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:239, 259
    // Minimum diagonal separation for HSP containment checking
    // Used in MB_HSP_CLOSE macro (blast_gapalign_priv.h:123-124)
    // megablast and dc-megablast: 6; blastn and blastn-short: 50 (`task_defaults`).
    let min_diag_separation = defaults.min_diag_separation;
    let (mb_template_type, mb_template_length) = determine_template(args);

    // No NCBI counterpart: placeholders; `finalize_task_config` sets the lookup kind and width.
    let use_two_stage = false;
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;
    let lut_word_length = effective_word_size;

    TaskConfig {
        effective_word_size,
        reward,
        penalty,
        matrix_only,
        score_reward,
        score_penalty,
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
        x_drop_ungapped,
        min_raw_gapped_score,
        scan_range,
        min_diag_separation,
        window_size: determine_window_size(args),
        mb_template_type,
        mb_template_length,
    }
}

/// Choose lookup table kind and LUT width using NCBI's BlastChooseNaLookupTable logic.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:45-185
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

    let (mut lut_kind, lut_width) = match word_size {
        4..=6 => (LookupTableKind::Small, word_size),
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c (BlastChooseNaLookupTable)
/// - Uses total query length (approx_table_entries) and max query offset
/// - Handles discontiguous template forcing MB lookup
pub fn finalize_task_config(
    config: &mut TaskConfig,
    total_query_length: usize,
    max_query_length: usize,
    discontig_template: bool,
) {
    let max_q_off = max_query_length.saturating_sub(1);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:128-130
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:397,1334
    config.scan_step = (config.effective_word_size - config.lut_word_length + 1).max(1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1413
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
/// `length(index)` is the number of letters of record `index` (its letters, or with
/// `LOSAT_X_BNSTORE` its packed subject's); the total is NCBI's `Int8` sum of them (a length
/// held in memory is at most `isize::MAX`, so each converts without loss).
pub fn subject_metadata_from_records<R: InputRecord>(
    records: &[R],
    length: impl Fn(usize) -> usize,
) -> SubjectMetadata {
    SubjectMetadata {
        db_len_total: (0..records.len()).map(|index| length(index) as i64).sum(),
        db_num_seqs: i32::try_from(records.len()).unwrap_or(i32::MAX),
    }
}

// No NCBI counterpart: the name of a subject in LOSAT's debug and trace lines
// (`LOSAT_DEBUG_BLASTN`, `blastn_trace`), made only when they are on; it changes no value
// NCBI computes.
/// The title of a subject record up to its first space, "unknown" when that is empty.
pub fn subject_trace_id(title: &[u8]) -> String {
    let end = title
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(title.len());
    if end == 0 {
        "unknown".to_string()
    } else {
        String::from_utf8_lossy(&title[..end]).into_owned()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:856-874
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1079-1089
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
fn collect_lowercase_masks_for_records<L: QueryLetters + ?Sized>(
    records: &L,
) -> Vec<Vec<MaskedInterval>> {
    (0..records.count())
        .map(|index| records.read(index, collect_lowercase_masks))
        .collect()
}

/// The letters of a query batch's records for the readers that set the batch up (WindowMasker,
/// DUST, the lower-case masks, the encoding), one record at a time. No NCBI counterpart (how
/// LOSAT holds the letters).
pub(crate) trait QueryLetters: Sync {
    /// The number of records.
    fn count(&self) -> usize;
    /// The number of letters of record `index`.
    fn length(&self, index: usize) -> usize;
    /// `read` of the letters of record `index`.
    fn read<T>(&self, index: usize, read: impl FnOnce(&[u8]) -> T) -> T;
}

impl<R: InputRecord + Sync> QueryLetters for [R] {
    fn count(&self) -> usize {
        self.len()
    }

    fn length(&self, index: usize) -> usize {
        self[index].seq().len()
    }

    fn read<T>(&self, index: usize, read: impl FnOnce(&[u8]) -> T) -> T {
        read(self[index].seq())
    }
}

impl QueryLetters for [&[u8]] {
    fn count(&self) -> usize {
        self.len()
    }

    fn length(&self, index: usize) -> usize {
        self[index].len()
    }

    fn read<T>(&self, index: usize, read: impl FnOnce(&[u8]) -> T) -> T {
        read(self[index])
    }
}

/// The reverse complements of `L`'s letters, each made for the reader of one record (`-strand
/// minus` DUST below).
struct ReverseComplements<'a, L: ?Sized>(&'a L);

impl<L: QueryLetters + ?Sized> QueryLetters for ReverseComplements<'_, L> {
    fn count(&self) -> usize {
        self.0.count()
    }

    fn length(&self, index: usize) -> usize {
        self.0.length(index)
    }

    fn read<T>(&self, index: usize, read: impl FnOnce(&[u8]) -> T) -> T {
        self.0.read(index, |letters| {
            read(&super::lookup::reverse_complement(letters))
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-128
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1173-1178
pub(crate) fn reverse_mask_intervals(
    masks: &[MaskedInterval],
    query_length: usize,
) -> Vec<MaskedInterval> {
    let mut out: Vec<MaskedInterval> = masks
        .iter()
        .map(|m| {
            let start = query_length.saturating_sub(m.end);
            let end = query_length.saturating_sub(m.start);
            MaskedInterval { start, end }
        })
        .collect();
    out.sort_by_key(|m| m.start);
    out
}

/// Apply DUST filter to query sequences
pub(crate) fn apply_dust_masking<L: QueryLetters + ?Sized>(
    args: &BlastnArgs,
    queries: &L,
    pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:418-420
    if let Some((dust_level, dust_window, dust_linker)) = args.dust.params() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-128
        if args.verbose {
            eprintln!(
                "Applying DUST filter (level={}, window={}, linker={})...",
                dust_level, dust_window, dust_linker
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-102
        // NCBI constructs a fresh CSymDustMasker for each query, which
        // resets the converter's CRandom state before scanning that query.
        // `DustMasker::new(..).mask_sequence(..)` is that, once per query. EXPERIMENT
        // (LOSAT_X_QPREPPAR): the queries are masked on the search pool, in query order
        // (`x_query_dust_par.rs`: the mask of a query is a function of its letters and the three
        // parameters only, so the results are the serial ones); without the switch, or without a
        // pool of two or more threads, the same map runs serially.
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:171-180
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:91-119
        // The vector is that of the query's Seq-loc with its strand: with `-strand minus` DUST
        // reads the reverse complement and its intervals are mapped back to the plus strand.
        // (The conversion to two bits reads letters other than C, G, T as A and N at random,
        // `symdust.hpp:76-84`, so on a query with ambiguity codes the intervals differ by strand.)
        // `plus` and `both` dust the letters as they are.
        let masks: Vec<Vec<MaskedInterval>> =
            if QueryStrand::from_option(&args.strand) == QueryStrand::Minus {
                let reverse = ReverseComplements(queries);
                super::x_query_dust_par::masks(&reverse, dust_level, dust_window, dust_linker, pool)
                    .into_iter()
                    .enumerate()
                    .map(|(index, masks)| reverse_mask_intervals(&masks, queries.length(index)))
                    .collect()
            } else {
                super::x_query_dust_par::masks(queries, dust_level, dust_window, dust_linker, pool)
            };

        let total_masked: usize = masks
            .iter()
            .map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>())
            .sum();
        let total_bases: usize = (0..queries.count())
            .map(|index| queries.length(index))
            .sum();
        if args.verbose && total_bases > 0 {
            eprintln!(
                "DUST masked {} bases ({:.2}%) across {} sequences",
                total_masked,
                100.0 * total_masked as f64 / total_bases as f64,
                queries.count()
            );
        }
        masks
    } else {
        vec![Vec::new(); queries.count()]
    }
}

/// Build lookup tables based on configuration. `x_keep`: the batch loop's keep of
/// LOSAT_X_MBLUTKEEP for the megablast table (`lookup::x_mb_lut_keep`), `None` elsewhere.
pub(crate) fn build_lookup_tables<R: InputRecord>(
    config: &TaskConfig,
    args: &BlastnArgs,
    queries_blastna: &[&[u8]],
    query_masks: &[Vec<MaskedInterval>],
    query_offsets: &[i32],
    subjects: &[R],
    subjects_packed: Option<&[Vec<u8>]>,
    approx_table_entries: usize,
    x_keep: Option<&std::sync::Arc<x_mb_lut_keep::Keep>>,
    x_scan_letters: Option<u64>,
) -> (LookupTables, usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1231-1344
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    let db_word_counts = if args.limit_lookup {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:865-868
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1027-1034
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1327-1341
    let two_stage_lookup: Option<TwoStageLookup> = if config.mb_template_length > 0 {
        Some(build_disc_mb_lookup_kept(
            queries_blastna,
            query_offsets,
            config.effective_word_size,
            config.mb_template_length as usize,
            config.mb_template_type,
            query_masks,
            approx_table_entries,
            x_keep,
        ))
    } else if config.use_two_stage {
        Some(build_two_stage_lookup_kept(
            queries_blastna,
            query_offsets,
            config.effective_word_size,
            config.lut_word_length,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
            approx_table_entries,
            !config.mb_lookup,
            x_keep,
            x_scan_letters,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:548-583
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

/// Why the masks of queries could not be made (`query_masks`, `chunk_query_masks`).
#[derive(Debug)]
pub(crate) enum MaskError {
    /// The counts file's or WindowMasker's error (`window_masker`).
    WindowMasker(crate::utils::winmask::WinMaskError),
    /// The repeat search's error (`-filtering_db`, `repeat_filter`).
    Repeats(super::repeat_filter::RepeatError),
}

/// The masks of the queries: their DUST masks and, with `-lcase_masking`, their lower-case
/// letters, then with `-filtering_db` the hits of the repeat search (`repeat_filter`); with
/// `-window_masker_db` the masks of WindowMasker, which replace DUST's and the lower case
/// (`window_masker`; the options exclude `-filtering_db`). With the masks, the repeat search's
/// hits of each query (`RepeatFilter::add_masks`; empty without `-filtering_db`).
#[allow(clippy::type_complexity)]
pub(crate) fn query_masks<L: QueryLetters + ?Sized>(
    args: &BlastnArgs,
    queries: &L,
    outside: &mut [super::repeat_filter::MasksOutside],
    pool: &SearchPool<'_>,
) -> Result<(Vec<Vec<MaskedInterval>>, Vec<Vec<(usize, usize)>>), MaskError> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:216-224
    // After DUST, and `SetMaskedRegions` of the query with `mqr` (winmask_filter.cpp:311) replaces the
    // query's masks: DUST's and the lower case are not made.
    if let Some(window_masker) = &args.window_masker {
        let minus = QueryStrand::from_option(&args.strand) == QueryStrand::Minus;
        return (0..queries.count())
            .map(|index| queries.read(index, |letters| window_masker.masks(letters, minus)))
            .collect::<Result<_, _>>()
            .map(|masks| (masks, Vec::new()))
            .map_err(MaskError::WindowMasker);
    }
    let mut query_masks = apply_dust_masking(args, queries, pool);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-128
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:206-214
    let repeat_hits = match &args.repeat_filter {
        Some(filter) => filter
            .add_masks(&args.strand, queries, &mut query_masks, outside, pool)
            .map_err(MaskError::Repeats)?,
        None => Vec::new(),
    };
    Ok((query_masks, repeat_hits))
}

/// The masks of the query parts of a query chunk: the DUST masks of each part added to the
/// masks that NCBI keeps from its query (`query_split::restrict_masks`), which already have
/// the lower-case letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:171-180
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:121-128
/// The part's own masks are only merged when DUST finds a region; the lookup table and the
/// search use the masked residues, which merging does not change. With `-filtering_db` the
/// chunk's query source searches the repeats of the parts' letters again and adds their hits
/// (`repeat_filter`). With `-window_masker_db` the chunk's `x_CalculateMasks` ends with
/// WindowMasker on the part's letters, which replaces every mask the part had (`query_masks`).
pub(crate) fn chunk_query_masks<R: InputRecord + Sync>(
    args: &BlastnArgs,
    parts: &[R],
    restricted: Vec<Vec<MaskedInterval>>,
    outside: &mut [super::repeat_filter::MasksOutside],
    pool: &SearchPool<'_>,
) -> Result<Vec<Vec<MaskedInterval>>, MaskError> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:876-880
    // A chunk's query data is made with the options (`MakeLocalQueryData` with `m_Options`), so its
    // query source calculates the masks of the chunk's parts again, WindowMasker last.
    if args.window_masker.is_some() {
        return query_masks(args, parts, &mut [], pool).map(|(masks, _)| masks);
    }
    let mut masks: Vec<Vec<MaskedInterval>> = apply_dust_masking(args, parts, pool)
        .into_iter()
        .zip(restricted)
        .map(|(mut masks, restricted)| {
            masks.extend(restricted);
            merge_mask_intervals(masks)
        })
        .collect();
    // The chunk's hits are its parts' own: the archive reads those of the batch (`query_masks`).
    if let Some(filter) = &args.repeat_filter {
        filter
            .add_masks(&args.strand, parts, &mut masks, outside, pool)
            .map_err(MaskError::Repeats)?;
    }
    Ok(masks)
}

/// Prepare all sequence data and configuration
pub fn prepare_sequence_data<'a, R>(
    queries: &'a [R],
    query_masks: Vec<Vec<MaskedInterval>>,
    subject_metadata: SubjectMetadata,
) -> SequenceData<'a, R> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1409
    let db_len_total = subject_metadata.db_len_total;
    let db_num_seqs = subject_metadata.db_num_seqs;

    SequenceData {
        queries,
        query_masks,
        db_len_total,
        db_num_seqs,
    }
}

#[cfg(test)]
#[path = "coordination_tests.rs"]
mod tests;
