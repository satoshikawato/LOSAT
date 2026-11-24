//! The scoring options of a BLASTN search: NCBI's checks of the given values, and the
//! Karlin-Altschul blocks of the queries.

use super::args::BlastnArgs;
use super::blast_engine::NO_DATA_MESSAGE;
use super::coordination::determine_scoring_params;
use super::pairwise::query_ungapped_karlin_blastna;
use crate::blastinput::report_counts::{run_report_counts_with, ReportCounts};
use crate::cli::NativeError;
use crate::config::NuclScoringSpec;
use crate::stats::{ungapped_beta, KarlinParams, NuclValues};
use crate::wording::{ENGINE_ERROR, QUERY_OPTIONS_ERROR};

fn scoring_spec(args: &BlastnArgs) -> NuclScoringSpec {
    let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(args);
    NuclScoringSpec {
        reward,
        penalty,
        gap_open,
        gap_extend,
    }
}

/// NCBI's checks of the options (`Validate`), before the query is read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:881-906
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1326-1333
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1247-1261,1399-1413
/// The template length and type are 16, 18 or 21 and 0, 1 or 2 by the argument
/// constraints (`-template_length`, `-template_type`); a caller of `run_local` can set another
/// length, which the length's check refuses (the type is an enum of the three).
/// The lookup table type of the options is the task's: megablast and dc-megablast
/// `eMBLookupTable`, blastn and blastn-short `eNaLookupTable`
/// (blast_nucl_options.cpp:137-150, disc_nucl_options.cpp:55-64).
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1699-1711
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1759-1776 (BLAST_ValidateOptions
/// runs the scoring, lookup table, hit saving and extension-scoring checks in this order)
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3636-3639
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
/// Only megablast extends greedily (dc-megablast and blastn-short use dynamic programming,
/// `coordination.rs` `task_defaults`). The reward and the penalty are NCBI's 16-bit values
/// (`determine_scoring_params`). LOSAT's own limits are `check_losat_limits`.
pub fn check_scoring_options(args: &BlastnArgs) -> anyhow::Result<()> {
    check_scoring_options_with(args, &args.report_counts, args.ungapped)
}

/// `check_scoring_options` with the hit list size of `counts` (the counts that
/// [`report_counts_of`] derives for a run; `args.report_counts` is set by the run), for an
/// ungapped (`-ungapped`) search when `ungapped`.
fn check_scoring_options_with(
    args: &BlastnArgs,
    counts: &ReportCounts,
    ungapped: bool,
) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    let word_size = super::coordination::determine_effective_word_size(args);
    // `-no_greedy` selects dynamic programming (blast_args.cpp:681-685).
    let greedy = args.task == "megablast" && !args.no_greedy;
    let mb_lookup_options = matches!(args.task.as_str(), "megablast" | "dc-megablast");
    let (_, template_length) = super::coordination::determine_template(args);
    let matrix_only = super::coordination::matrix_only_scoring(spec.reward, spec.penalty);
    let message = if !matrix_only && spec.penalty >= 0 {
        "BLASTN penalty must be negative".to_string()
    } else if !ungapped && spec.gap_open > 0 && spec.gap_extend == 0 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:900-904
        "BLASTN gap extension penalty cannot be 0".to_string()
    } else if word_size == 0 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1326
        // The command line's `-word_size` is at least 4 (`blastn_word_size`); a caller of
        // `run_local` sets the field itself.
        "Word-size must be greater than zero".to_string()
    } else if word_size < 4 {
        "Word-size must be 4 or greater for nucleotide comparison".to_string()
    } else if word_size > DBSEQ_CHUNK_OVERLAP {
        format!("Word-size must be less than or equal to {DBSEQ_CHUNK_OVERLAP}")
    } else if template_length > 0 && word_size != 11 && word_size != 12 {
        "Invalid discontiguous template parameters: word size must be either 11 or 12".to_string()
    } else if template_length > 0 && !matches!(template_length, 16 | 18 | 21) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1263-1269
        // The command line's `-template_length` is one of the three (`blastn_template_length`);
        // a caller of `run_local` sets the field itself.
        "Invalid discontiguous template parameters: template length must be 16, 18, or 21"
            .to_string()
    } else if template_length > 0 && args.strategy_template_type.is_some() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1271-1277
        // A search strategy's type above 2 (`MBTemplateType`; the command line names the three).
        "Invalid discontiguous template parameters: template type must be 0, 1, or 2".to_string()
    } else if template_length > 0 && !mb_lookup_options {
        "Invalid lookup table type for discontiguous Mega BLAST".to_string()
    } else if args.limit_lookup && word_size < 16 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1422-1427
        // `limit_lookup` has no command-line spelling in LOSAT; a caller of `run_local` sets it.
        "The limit_lookup option can only be used with word size >= 16".to_string()
    } else if args.off_diagonal_range.is_some_and(|range| range != 0)
        && super::coordination::determine_window_size(args) == 0
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:620-624,1583-1596
        // (`BlastInitialWordOptionsValidate`, after the lookup table's checks and before the
        // hit saving's)
        "off_diagonal_range is only useful in 2-hit algorithm".to_string()
    } else if counts.hitlist_size < 1 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1511-1517, `BlastHitSavingOptionsValidate`, whose hit list
        // check comes before the e-value's (the message of `ReportCounts::validate`, which
        // is not called because this check has its place among the others).
        "No hits are being saved".to_string()
    } else if super::coordination::determine_evalue(args) <= 0.0
        && args.min_raw_gapped_score.unwrap_or(0) <= 0
    {
        "expect value or cutoff score must be greater than zero".to_string()
    } else if let Some(message) = crate::blastinput::app::hsp_filtering_error(
        args.culling_limit,
        args.best_hit_overhang,
        args.best_hit_score_edge,
    ) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1543-1548 (the HSP
        // filtering options, the last check of the hit saving options)
        message.to_string()
    } else if spec.gap_open == 0
        && spec.gap_extend == 0
        && match args.strategy_algorithms {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1700-1712
            // A search strategy's pair: greedy where either algorithm is.
            Some((extension, traceback)) => extension != 1 && traceback != 1,
            None => !greedy,
        }
    {
        "Greedy extension must be used if gap existence and extension options are zero".to_string()
    } else if args.strategy_algorithms.is_some() {
        // NCBI's checks pass: NCBI runs the pair (its jumper, Smith-Waterman or a greedy
        // extension with dynamic programming's traceback), which LOSAT does not.
        return Err(crate::report::archive::read::import::Start::unported(
            "a search strategy with this pair of gapped extension algorithms",
            crate::blastinput::app::Program::Blastn.name(),
        ));
    } else {
        return Ok(());
    };
    // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
    let footer = crate::wording::USER_MANUAL_FOOTER;
    Err(NativeError {
        exit: 1,
        message: format!("{QUERY_OPTIONS_ERROR} {message}\n{footer}\n").into(),
    }
    .into())
}

/// LOSAT's limits on the options that NCBI accepts, checked where NCBI starts the search
/// (after its `Query is Empty!` success, before its Karlin-Altschul table error, whose
/// computation the range limit bounds): a reward of 0 or less with a penalty that is not 0
/// (no valid query; reward 0 with penalty 0 is matrix-only scoring, which runs) and a reward of
/// `BLAST_SCORE_MAX` or a penalty of `BLAST_SCORE_MIN` (the 16-bit values after NCBI's
/// conversion). The limit on greedy gap costs follows the table check
/// (`check_greedy_gap_costs`).
pub fn check_losat_limits(args: &BlastnArgs) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    // Reward 0 with penalty 0 is matrix-only scoring (`coordination::engine_scoring_pair`).
    // Another reward of 0 or less leaves every query without Karlin-Altschul values (NCBI:
    // warnings and no hits), an explicit rejection of PD-LOSAT-NCBI-DEFECTS (version 1.1).
    if spec.reward <= 0 && !super::coordination::matrix_only_scoring(spec.reward, spec.penalty) {
        anyhow::bail!(
            "a reward of {} with a penalty of {} (NCBI BLAST+'s 16-bit values), for which no query has Karlin-Altschul values, is not supported by LOSAT's BLASTN; reward 0 needs penalty 0 (matrix-only scoring, as -task rmblastn)",
            spec.reward,
            spec.penalty
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1506-1520
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2175-2179
    // NCBI leaves a reward of BLAST_SCORE_MAX (32767) or a penalty of BLAST_SCORE_MIN
    // (-32768) out of the score range of the ungapped Karlin-Altschul computation. It then
    // counts a reward of 32767 beyond the end of the frequency array (oracle: invalid
    // queries, or a crash), and without a penalty of -32768 every query is invalid (oracle:
    // the warnings and no hits for a 3000-base query, a crash for primer-length queries of
    // every task, SD audit (c)), which LOSAT's computation does not reproduce. LOSAT's
    // scores below these values were compared with NCBI's up to 24000/-30000 (E2g V1).
    if spec.reward >= i32::from(i16::MAX) || spec.penalty <= i32::from(i16::MIN) {
        anyhow::bail!(
            "a reward of 32767 or a penalty of -32768 (NCBI BLAST+'s BLAST_SCORE_MAX and BLAST_SCORE_MIN, which it leaves out of the score range of its Karlin-Altschul computation) is not supported by LOSAT's BLASTN; reward {} and penalty {} were given",
            spec.reward,
            spec.penalty
        );
    }
    Ok(())
}

/// LOSAT's limit on gap costs with greedy extension (`-task megablast`), for scoring that
/// NCBI's Karlin-Altschul tables accept: NCBI sizes the greedy distances in 32-bit
/// integers (`blast_gapalign.c:233-236`), which larger costs can overflow.
pub fn check_greedy_gap_costs(args: &BlastnArgs) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    if args.task == "megablast" && spec.gap_open.max(spec.gap_extend) > MAX_GREEDY_GAP_COST {
        anyhow::bail!(
            "gap costs above {MAX_GREEDY_GAP_COST} with greedy extension (-task megablast) are not supported by LOSAT's BLASTN"
        );
    }
    Ok(())
}

/// LOSAT's refusal of matrix-only scoring (reward 0 and penalty 0) with the greedy extension of
/// `-task megablast`, checked where `check_greedy_gap_costs` is (a batch with a valid query).
/// NCBI sets up the greedy alignment with the option reward and penalty: with gap costs 0 and 0
/// (megablast's own, also in an ungapped search, which sets it up too) it divides by their sum,
/// 0 (SIGFPE, exit 136); with other gap costs the set-up passes, and the alignment divides by the
/// half reward, 0, converts the infinite quotient to `Int4` (undefined behaviour) and does not
/// end. An ungapped search with other gap costs runs no greedy alignment, and `-no_greedy` runs
/// dynamic programming: both run (`TaskConfig::matrix_only`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:183,191-193,208-209
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:979-981 (`match_score_half`
/// is `match_score / 2`, :837)
pub fn check_greedy_matrix_only(args: &BlastnArgs, ungapped: bool) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    let greedy = args.task == "megablast" && !args.no_greedy;
    let matrix_only = super::coordination::matrix_only_scoring(spec.reward, spec.penalty);
    let zero_gap_costs = spec.gap_open == 0 && spec.gap_extend == 0;
    if greedy && matrix_only && (zero_gap_costs || !ungapped) {
        anyhow::bail!(
            "reward 0 and penalty 0 (matrix-only scoring) with the greedy extension of -task megablast, where NCBI BLAST+ divides by zero (gap costs 0 and 0, also with -ungapped) or does not end (other gap costs; -no_greedy or -ungapped runs them), is not supported by LOSAT's BLASTN"
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:192
use crate::core::subject_split::DBSEQ_CHUNK_OVERLAP;

/// The largest gap cost that LOSAT's BLASTN accepts with greedy extension.
const MAX_GREEDY_GAP_COST: i32 = 32767;

/// NCBI's error for scores or gap costs without Karlin-Altschul values, raised when the
/// first batch with a valid query is set up: every query of that batch gets the message,
/// and a query without letters (`true` in `empty_queries`, one item per query of the
/// batch) keeps its set-up warning after it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:94-102
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1013-1024
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:180-183
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640 (the query's
/// warning, before the score block: its id is the query's index)
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1042-1052 (after the error is
/// added, blast_aux_priv.cpp:107-108)
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:294-303
/// `operator<` is no order: the error is below the warning by id (`kBlastMessageNoContext`,
/// -1, against the query's index) and the warning below the error by severity. Of two
/// messages, libstdc++'s insertion sort moves the second one (the error, added after the
/// set-up warning) to the front, so the error comes first (TBLASTX's unknown matrix has the
/// same order, `run_impl.rs`).
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:170-184
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
pub(crate) fn karlin_error(
    message: &str,
    empty_queries: impl IntoIterator<Item = bool>,
) -> anyhow::Error {
    let mut text = format!("{ENGINE_ERROR} ");
    for empty in empty_queries {
        text.push_str("Error: ");
        text.push_str(message);
        text.push(' ');
        if empty {
            text.push_str("Warning: ");
            text.push_str(NO_DATA_MESSAGE);
            text.push(' ');
        }
    }
    text.push('\n');
    NativeError {
        exit: 3,
        message: text.into(),
    }
    .into()
}

/// The hit list size and the counts of the pairwise report of a run that writes the output
/// formats `outfmts` (`-max_target_seqs`, `-num_descriptions`, `-num_alignments` of `args`),
/// with the formatting handler's warnings: the one derivation of `process_options`.
pub(crate) fn report_counts_of(
    args: &BlastnArgs,
    outfmts: &[i32],
) -> anyhow::Result<(ReportCounts, Vec<u8>)> {
    run_report_counts_with(
        args.app().ncbi_name(),
        outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
        &args.formatting,
    )
}

/// The Karlin-Altschul blocks of one query context (a strand of a query).
#[derive(Debug, Clone, Copy)]
pub(crate) struct ContextKarlin {
    /// From the context's composition (NCBI `kbp_std`): the gap trigger and the ungapped
    /// X-drop.
    pub ungapped: KarlinParams,
    /// The block that NCBI's `gapped_calculation` selects for e-values, bit scores, cutoffs
    /// and the length adjustment: the gapped block (NCBI `kbp_gap`) with its alpha and beta,
    /// or, in an ungapped (`-ungapped`) search, the ungapped block (`kbp_std`) with the alpha
    /// and beta of `Blast_GetNuclAlphaBeta` (see [`context_blocks_ungapped`]).
    pub gapped: KarlinParams,
}

/// The ungapped block of every query context, from the residues of the context (a query
/// or its reverse complement); `None` for a context whose block cannot be computed, which
/// NCBI does not search (both strands of an invalid query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2762-2784
/// The two strands have the same composition in another order, so their blocks can
/// differ in the last bits.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2772
/// A context without letters (a query without letters) is invalid from its set-up
/// (`s_QueryInfo_SetContext`): it gets no block and no message here.
/// The contexts come in BLASTNA (the search holds them in no other encoding).
pub(crate) fn context_ungapped_blocks<'a>(
    contexts: impl Iterator<Item = &'a [u8]>,
    spec: &NuclScoringSpec,
) -> Vec<Option<KarlinParams>> {
    // Twin: the closure body below is copied in `parallel_blocks` of
    // blast_engine/x_query_prep_par.rs (LOSAT_X_QPREPPAR); change it there too.
    contexts
        .map(|context| {
            if context.is_empty() {
                None
            } else {
                query_ungapped_karlin_blastna(context, spec.reward, spec.penalty)
            }
        })
        .collect()
}

/// The blocks of every valid context, and whether scores are rounded down to even
/// numbers. With no valid context, NCBI computes no gapped block and raises no error.
/// The error is NCBI's message (`NuclValues`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:73-107
pub(crate) fn context_blocks(
    ungapped: &[Option<KarlinParams>],
    spec: &NuclScoringSpec,
) -> Result<(Vec<Option<ContextKarlin>>, bool), String> {
    if ungapped.iter().all(Option::is_none) {
        return Ok((vec![None; ungapped.len()], false));
    }
    let values = NuclValues::new(spec.reward, spec.penalty)?;
    let blocks = ungapped
        .iter()
        .map(|block| {
            block
                .map(|ungapped| {
                    values
                        .gapped(spec, &ungapped)
                        .map(|gapped| ContextKarlin { ungapped, gapped })
                })
                .transpose()
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok((blocks, values.round_down))
}

/// The blocks of every valid context of an ungapped (`-ungapped`) search, and whether scores
/// are rounded down to even numbers (never). NCBI builds no gapped block then, so no table is
/// asked for the gap costs and a pair of scores without a table row is no error: its alpha
/// and beta stay 0 and the length adjustment is 0.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:513-529
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:705,768,806-812,821-824
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3989-3992, 4024-4029
pub(crate) fn context_blocks_ungapped(
    ungapped: &[Option<KarlinParams>],
    spec: &NuclScoringSpec,
) -> (Vec<Option<ContextKarlin>>, bool) {
    let tabled = NuclValues::new(spec.reward, spec.penalty).is_ok();
    let blocks = ungapped
        .iter()
        .map(|block| {
            block.map(|ungapped| {
                let (alpha, beta) = if tabled {
                    (
                        ungapped.lambda / ungapped.h,
                        ungapped_beta(spec.reward, spec.penalty),
                    )
                } else {
                    (0.0, 0.0)
                };
                ContextKarlin {
                    ungapped,
                    gapped: KarlinParams {
                        alpha,
                        beta,
                        ..ungapped
                    },
                }
            })
        })
        .collect();
    (blocks, false)
}

/// The gapped X-drops of the search (raw scores of the preliminary and the final
/// extension), from the smallest gapped Lambda of the valid contexts of a batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:92-116
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
/// Gapped blocks from the tables are the same for every context; blocks copied from the
/// ungapped ones differ with the composition, so the X-drops depend on the contexts of the
/// batch. Without a valid context the batch is not searched. The conversion is the protein
/// programs' (`core::x_drop::gapped_x_dropoffs`): a value past `Int4` is `INT_MIN`, as the
/// x86-64 build of NCBI BLAST+ 2.17.0 converts it.
pub(crate) fn gap_x_dropoffs(
    context_karlin: &[Option<ContextKarlin>],
    x_drop_gapped_bits: f64,
    x_drop_final_bits: f64,
) -> (i32, i32) {
    let min_lambda = context_karlin
        .iter()
        .flatten()
        .map(|blocks| blocks.gapped.lambda)
        .reduce(f64::min)
        .unwrap_or(KarlinParams::default().lambda);
    crate::core::x_drop::gapped_x_dropoffs(x_drop_gapped_bits, x_drop_final_bits, min_lambda)
}

/// The refusal, before a batch's subjects are searched, of the gapped X-drops (raw,
/// [`gap_x_dropoffs`]) that NCBI's kernels cannot run with the search's scores:
/// - dynamic programming (`use_dp`): the final X-drop of `ALIGN_EX` past the sizes that NCBI keeps
///   in `Int4` (`core::x_drop::kernel_x_drop_fits`, as the protein programs refuse it). The
///   preliminary stage runs with the smaller of the X-drop and the HSP's ungapped score
///   (`s_BlastDynProgNtGappedAlignment`), so no preliminary X-drop is refused;
/// - greedy extension (`-task megablast`): an X-drop of either stage whose `xdrop_offset`
///   (`alignment::greedy::greedy_xdrop_offset`) is 0 or less, or above the offset of the larger
///   X-drop, by which NCBI sizes `max_score` (`d_diff`): NCBI then reads and writes outside its
///   `max_score` array or fails to allocate it (NCBI 2.17.0: `BLAST engine error: Unknown error
///   code -1`, a crash on a corrupted heap, or SIGSEGV, line RG X1b probe and review 1). A
///   batch that NCBI would search without any greedy extension (no seed) is refused too, as the
///   protein programs refuse their kernels' X-drops before the search.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:183-188, 207, 257
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2959-2963
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:496-498
///
/// Under matrix-only scoring (`matrix_only`) the X-drops are the options' values as raw scores
/// (`subject_env`, from blast_parameters.c:208-217 and 472-482; the final one is not raised to
/// the preliminary one), and the search extends by dynamic programming
/// (`check_greedy_matrix_only` refuses the greedy extension before). Only the final X-drop is
/// checked, as in the other DP searches: the preliminary DP caps its X-drop at the HSP's ungapped
/// score (blast_gapalign.c:2959-2963).
#[allow(clippy::too_many_arguments)]
pub(crate) fn check_gapped_x_dropoffs(
    use_dp: bool,
    matrix_only: bool,
    (gap_bits, final_bits): (f64, f64),
    (gapped, finale): (i32, i32),
    (reward, penalty): (i32, i32),
    (gap_open, gap_extend): (i32, i32),
) -> anyhow::Result<()> {
    // The final X-drop is the larger of its own and the preliminary one: when they are equal,
    // `-xdrop_gap` set it.
    let (final_option, final_bits) = if finale == gapped && !matrix_only {
        ("-xdrop_gap", gap_bits)
    } else {
        ("-xdrop_gap_final", final_bits)
    };
    if use_dp || matrix_only {
        if !crate::core::x_drop::kernel_x_drop_fits(finale, gap_open, gap_extend) {
            if matrix_only {
                return Err(crate::core::x_drop::kernel_x_drop_value_error(
                    "BLASTN",
                    final_option,
                    &format!("{final_bits:?} (a raw score with reward 0 and penalty 0)"),
                ));
            }
            return Err(crate::core::x_drop::kernel_x_drop_error(
                "BLASTN",
                final_option,
                final_bits,
            ));
        }
        return Ok(());
    }
    use super::alignment::greedy::greedy_xdrop_offset;
    // NCBI's `max_score` holds `max_d` plus 1 plus `d_diff` cells, `d_diff` the offset of the
    // larger raw X-drop (`s_BlastGreedyAlignMemAlloc`). With an odd reward the doubling wraps, so
    // the smaller X-drop can have the larger offset (`-xdrop_gap -3e9`: 174306649 against 19),
    // which NCBI zero-fills past the array (SIGSEGV, review 1 of line RG X1b).
    // A `d_diff` of 0 or less is the offset of one of the two stages: that stage is named.
    let d_diff = greedy_xdrop_offset(gapped.max(finale), reward, penalty, gap_open, gap_extend);
    let stages = [
        ("-xdrop_gap", gap_bits, gapped),
        (final_option, final_bits, finale),
    ];
    let offset = |x_dropoff| greedy_xdrop_offset(x_dropoff, reward, penalty, gap_open, gap_extend);
    let culprit = stages
        .iter()
        .find(|&&(_, _, x_dropoff)| offset(x_dropoff) <= 0)
        .or_else(|| {
            stages
                .iter()
                .find(|&&(_, _, x_dropoff)| offset(x_dropoff) > d_diff)
        });
    if let Some(&(option, bits, _)) = culprit {
        anyhow::bail!(
            "an {option} of {bits:?} bits, whose raw X-drop makes NCBI BLAST+'s greedy alignment (-task megablast) read and write outside its arrays or fail to allocate them, is not supported by LOSAT's BLASTN"
        );
    }
    Ok(())
}

#[cfg(test)]
#[path = "scoring_tests.rs"]
mod tests;
