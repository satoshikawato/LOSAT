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
    let matrix_only = spec.reward == 0 && spec.penalty == 0;
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
    } else if spec.gap_open == 0 && spec.gap_extend == 0 && !greedy {
        "Greedy extension must be used if gap existence and extension options are zero".to_string()
    } else {
        return Ok(());
    };
    // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
    let footer = crate::wording::USER_MANUAL_FOOTER;
    Err(NativeError {
        exit: 1,
        message: format!("BLAST query/options error: {message}\n{footer}\n"),
    }
    .into())
}

/// LOSAT's limits on the options that NCBI accepts, checked where NCBI starts the search
/// (after its `Query is Empty!` success, before its Karlin-Altschul table error, whose
/// computation the range limit bounds): a reward of 0 or less (NCBI's rmblastn matrix
/// scoring when the penalty is 0 too, otherwise no valid query), a hit list size whose
/// preliminary size overflows (`get_prelim_hitlist_size`),
/// and a reward of `BLAST_SCORE_MAX` or a penalty of `BLAST_SCORE_MIN` (the 16-bit values
/// after NCBI's conversion). The limit on greedy gap costs follows the table check
/// (`check_greedy_gap_costs`).
pub fn check_losat_limits(args: &BlastnArgs) -> anyhow::Result<()> {
    check_losat_limits_with(args, &args.report_counts, args.ungapped)
}

/// `check_losat_limits` with the hit list size of `counts`, see `check_scoring_options_with`.
/// An ungapped search keeps `hitlist_size` preliminary hits (`get_prelim_hitlist_size` with
/// `gapped_calculation` false), which cannot overflow.
fn check_losat_limits_with(
    args: &BlastnArgs,
    counts: &ReportCounts,
    ungapped: bool,
) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    if spec.reward <= 0 {
        anyhow::bail!(
            "a reward of {} (NCBI BLAST+'s 16-bit value; 0 is NCBI's matrix scoring of rmblastn, otherwise no query is valid) is not supported by LOSAT's BLASTN",
            spec.reward
        );
    }
    super::hsp::check_prelim_hitlist_size(
        counts.hitlist_size,
        args.max_target_seqs.is_some(),
        false,
        !ungapped,
        "BLASTN",
    )?;
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
    let mut text = String::from("BLAST engine error: ");
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
        message: text,
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
        "blastn",
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:701-705, 768-848
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3008-3010
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:496-498
pub(crate) fn check_gapped_x_dropoffs(
    use_dp: bool,
    (gap_bits, final_bits): (f64, f64),
    (gapped, finale): (i32, i32),
    (reward, penalty): (i32, i32),
    (gap_open, gap_extend): (i32, i32),
) -> anyhow::Result<()> {
    // The final X-drop is the larger of its own and the preliminary one: when they are equal,
    // `-xdrop_gap` set it.
    let (final_option, final_bits) = if finale == gapped {
        ("-xdrop_gap", gap_bits)
    } else {
        ("-xdrop_gap_final", final_bits)
    };
    if use_dp {
        if !crate::core::x_drop::kernel_x_drop_fits(finale, gap_open, gap_extend) {
            return Err(crate::core::x_drop::kernel_x_drop_error(
                "BLASTN",
                final_option,
                final_bits,
            ));
        }
        return Ok(());
    }
    use super::alignment::greedy::greedy_xdrop_offset;
    // NCBI's `max_score` holds `max_d + 1 + d_diff` cells, `d_diff` the offset of the larger raw
    // X-drop (`s_BlastGreedyAlignMemAlloc`). With an odd reward the doubling wraps, so the smaller
    // X-drop can have the larger offset (`-xdrop_gap -3e9`: 174306649 against 19), which NCBI
    // zero-fills past the array (SIGSEGV, review 1 of line RG X1b).
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
mod tests {
    use super::super::pairwise::query_ungapped_karlin;
    use super::*;

    fn args(words: &[&str]) -> BlastnArgs {
        let argv = [
            "losat", "blastn", "-query", "q", "-subject", "s", "-outfmt", "6",
        ];
        let cli: crate::cli::Cli =
            crate::cli::try_parse_from(argv.iter().chain(words)).expect("valid arguments");
        let crate::cli::Commands::Blastn(args) = cli.command else {
            panic!("blastn")
        };
        args
    }

    fn message(result: anyhow::Result<()>) -> String {
        result
            .expect_err("rejected")
            .downcast::<NativeError>()
            .expect("an NCBI error")
            .message
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1247-1261,1399-1413
    #[test]
    fn discontiguous_templates_are_checked_after_the_word_size_limit() {
        let word = "Invalid discontiguous template parameters: word size must be either 11 or 12";
        let table = "Invalid lookup table type for discontiguous Mega BLAST";
        assert!(check_scoring_options(&args(&["-task", "dc-megablast"])).is_ok());
        assert!(
            check_scoring_options(&args(&["-task", "dc-megablast", "-word_size", "12"])).is_ok()
        );
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-word_size",
            "13"
        ])))
        .contains(word));
        // The template applies to any task: megablast keeps its MB table, blastn and
        // blastn-short have the standard table.
        let template = ["-template_type", "coding", "-template_length", "18"];
        let with = |words: &[&str]| args(&[words, &template[..]].concat());
        assert!(check_scoring_options(&with(&["-word_size", "11"])).is_ok());
        assert!(message(check_scoring_options(&with(&[]))).contains(word));
        assert!(message(check_scoring_options(&with(&["-task", "blastn"]))).contains(table));
        assert!(message(check_scoring_options(&with(&["-task", "blastn-short"]))).contains(word));
        assert!(message(check_scoring_options(&with(&[
            "-task",
            "blastn-short",
            "-word_size",
            "12"
        ])))
        .contains(table));
        // The word size limit (100) comes first, the zero gap costs last.
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-word_size",
            "101"
        ])))
        .contains("Word-size must be less than or equal to 100"));
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-word_size",
            "13",
            "-gapopen",
            "0",
            "-gapextend",
            "0"
        ])))
        .contains(word));
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-gapopen",
            "0",
            "-gapextend",
            "0"
        ])))
        .contains("Greedy extension must be used"));
    }

    // DI-6: `-no_greedy` makes megablast's gap costs 0/0 invalid (blast_options.c:1699-1711),
    // and an off-diagonal range needs a window (blast_options.c:1583-1596), checked after the
    // lookup table's options.
    #[test]
    fn no_greedy_and_off_diagonal_range_are_checked() {
        assert!(message(check_scoring_options(&args(&["-no_greedy"])))
            .contains("Greedy extension must be used"));
        assert!(
            check_scoring_options(&args(&["-no_greedy", "-gapopen", "5", "-gapextend", "2"]))
                .is_ok()
        );
        assert!(check_scoring_options(&args(&["-task", "blastn", "-no_greedy"])).is_ok());
        assert_eq!(
            message(check_scoring_options(&args(&["-off_diagonal_range", "2"]))),
            "BLAST query/options error: off_diagonal_range is only useful in 2-hit algorithm\nPlease refer to the LOSAT user manual.\n"
        );
        assert!(check_scoring_options(&args(&["-off_diagonal_range", "0"])).is_ok());
        assert!(
            check_scoring_options(&args(&["-off_diagonal_range", "2", "-window_size", "40"]))
                .is_ok()
        );
        assert!(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-off_diagonal_range",
            "2"
        ]))
        .is_ok());
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "dc-megablast",
            "-window_size",
            "0",
            "-off_diagonal_range",
            "2"
        ])))
        .contains("off_diagonal_range is only useful"));
        // The lookup table's check comes first.
        assert!(message(check_scoring_options(&args(&[
            "-off_diagonal_range",
            "2",
            "-word_size",
            "101"
        ])))
        .contains("Word-size must be less than"));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518: the e-value check needs a cutoff score of 0 or below too.
    #[test]
    fn min_raw_gapped_score_above_zero_allows_an_evalue_of_zero() {
        const REFUSED: &str = "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the LOSAT user manual.\n";
        for words in [
            &["-evalue", "0"][..],
            &["-evalue", "0", "-min_raw_gapped_score", "0"],
            &["-evalue", "-1", "-min_raw_gapped_score", "-5"],
        ] {
            assert_eq!(message(check_scoring_options(&args(words))), REFUSED);
        }
        for words in [
            &["-evalue", "0", "-min_raw_gapped_score", "1"][..],
            &["-evalue", "-1", "-min_raw_gapped_score", "40"],
            &["-min_raw_gapped_score", "0"],
            &["-min_raw_gapped_score", "-5"],
        ] {
            assert!(check_scoring_options(&args(words)).is_ok(), "{words:?}");
        }
    }

    #[test]
    fn options_are_checked_in_ncbi_order() {
        assert!(check_scoring_options(&args(&[])).is_ok());
        assert!(check_scoring_options(&args(&["-task", "blastn"])).is_ok());
        assert_eq!(
            message(check_scoring_options(&args(&[
                "-task", "blastn", "-penalty", "0", "-gapopen", "0", "-gapextend", "0"
            ]))),
            "BLAST query/options error: BLASTN penalty must be negative\nPlease refer to the LOSAT user manual.\n"
        );
        assert!(message(check_scoring_options(&args(&["-gapopen", "3"])))
            .contains("BLASTN gap extension penalty cannot be 0"));
        assert!(message(check_scoring_options(&args(&[
            "-task",
            "blastn",
            "-gapopen",
            "0",
            "-gapextend",
            "0"
        ])))
        .contains("Greedy extension must be used"));
        // LOSAT's limits are separate: NCBI's checks accept these options.
        for words in [
            &["-reward", "0"][..],
            &["-reward", "32767", "-penalty", "-1"],
            &["-reward", "16384", "-penalty", "-32768"],
            // 98303 is 32767 as NCBI's 16-bit reward.
            &["-reward", "98303", "-penalty", "-2"],
        ] {
            let parsed = args(words);
            assert!(check_scoring_options(&parsed).is_ok(), "{words:?}");
            let error = check_losat_limits(&parsed).unwrap_err().to_string();
            assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
        }
        // Scores below the 16-bit limits are accepted (E2g V1), and so are the infinite
        // and NaN e-values that NCBI searches with (E2g).
        for words in [
            &["-reward", "5000", "-penalty", "-1"][..],
            &["-reward", "32766", "-penalty", "-32767"],
            &["-evalue", "+inf"],
            &["-evalue", "+nan"],
            &["-evalue", "-nan(7)"],
        ] {
            assert!(check_losat_limits(&args(words)).is_ok(), "{words:?}");
        }
    }

    /// The checks of the options that a run makes (`process_options`): the counts of the formats
    /// `[0]` from `report_counts_of`, then `check_scoring_options` and `check_losat_limits`.
    fn run_checks(words: &[&str]) -> anyhow::Result<()> {
        let mut args = args(words);
        args.report_counts = report_counts_of(&args, &[0]).unwrap().0;
        check_scoring_options(&args)?;
        check_losat_limits(&args)
    }

    /// `run_checks` for a run of the formats `outfmts`, gapped or not (`ungapped`).
    fn run_checks_with(words: &[&str], outfmts: &[i32], ungapped: bool) -> anyhow::Result<()> {
        let args = args(words);
        let (counts, _) = report_counts_of(&args, outfmts).unwrap();
        check_scoring_options_with(&args, &counts, ungapped)?;
        check_losat_limits_with(&args, &counts, ungapped)
    }

    #[test]
    fn the_counts_a_run_derives_give_the_refusals_of_the_options() {
        let error = message(run_checks(&[
            "-num_descriptions",
            "0",
            "-num_alignments",
            "0",
        ]));
        assert_eq!(
            error,
            "BLAST query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n"
        );
        // The hit list of format 0 is the larger of the two numbers.
        assert!(run_checks(&["-num_alignments", "0"]).is_ok());
        assert!(run_checks(&["-num_descriptions", "0"]).is_ok());
        assert!(run_checks(&["-num_descriptions", "1", "-num_alignments", "0"]).is_ok());
        // The preliminary hit list size overflows.
        let error = run_checks(&["-num_alignments", "2147483598"])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("(-num_descriptions / -num_alignments)")
                && error.contains("not supported by LOSAT's BLASTN"),
            "{error}"
        );
        let error = run_checks(&["-num_descriptions", "2147483598"])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("whose preliminary hit list size overflows"),
            "{error}"
        );
        // `-max_target_seqs` keeps its own wording (the default path, before O1b).
        let error = run_checks(&["-max_target_seqs", "2147483598"])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("a -max_target_seqs of 2147483598")
                && error.contains("whose preliminary hit list size overflows"),
            "{error}"
        );
        // Without the counts of a run (`args.report_counts` at its default), no hit list is empty.
        assert!(
            check_scoring_options(&args(&["-num_descriptions", "0", "-num_alignments", "0"]))
                .is_ok()
        );
    }

    #[test]
    fn table_errors_repeat_for_each_query_of_the_batch() {
        let error = karlin_error("Substitution scores 1 and -6 are not supported", [false; 2]);
        assert_eq!(
            error.downcast::<NativeError>().unwrap().message,
            "BLAST engine error: Error: Substitution scores 1 and -6 are not supported Error: Substitution scores 1 and -6 are not supported \n"
        );
        // A query without letters keeps its warning after the error (NCBI 2.17.0).
        let error = karlin_error(
            "Substitution scores 1 and -6 are not supported",
            [true, false],
        );
        assert_eq!(
            error.downcast::<NativeError>().unwrap().message,
            "BLAST engine error: Error: Substitution scores 1 and -6 are not supported Warning: Sequence contains no data Error: Substitution scores 1 and -6 are not supported \n"
        );
    }

    #[test]
    fn no_valid_query_needs_no_gapped_block() {
        let spec = NuclScoringSpec {
            reward: 1,
            penalty: -6,
            gap_open: 0,
            gap_extend: 0,
        };
        assert!(context_blocks(&[None, None], &spec).is_ok());
        let ungapped = query_ungapped_karlin(b"ACGTACGT", 1, -6);
        assert!(context_blocks(&[ungapped, ungapped], &spec).is_err());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3955-3963, 3989-3992, 4024-4029
    // (oracle `-ungapped`: 2/-3 beta -2, 1/-1 beta -2, 2/-2 and 4/-6 beta 0, 1/-6 and 3/-1 no table
    // row: alpha = beta = 0 and a length adjustment of 0, effective search space 377 x 373)
    #[test]
    fn ungapped_blocks_take_alpha_beta_from_the_ungapped_block() {
        let spec = |reward, penalty| NuclScoringSpec {
            reward,
            penalty,
            gap_open: 0,
            gap_extend: 0,
        };
        let block = KarlinParams {
            lambda: 0.634,
            k: 0.408,
            h: 0.912,
            alpha: 7.0,
            beta: 7.0,
        };
        for (reward, penalty, beta) in [
            (2, -3, -2.0),
            (1, -1, -2.0),
            (2, -2, 0.0),
            (4, -6, 0.0),
            (1, -2, 0.0),
        ] {
            let (blocks, round_down) =
                context_blocks_ungapped(&[Some(block), None], &spec(reward, penalty));
            assert!(!round_down);
            assert!(blocks[1].is_none());
            let blocks = blocks[0].unwrap();
            assert_eq!(blocks.ungapped.lambda, block.lambda);
            assert_eq!(blocks.gapped.lambda, block.lambda);
            assert_eq!(blocks.gapped.k, block.k);
            assert_eq!(blocks.gapped.h, block.h);
            assert_eq!(
                blocks.gapped.alpha,
                block.lambda / block.h,
                "{reward}/{penalty}"
            );
            assert_eq!(blocks.gapped.beta, beta, "{reward}/{penalty}");
        }
        // A pair with no table row (1/-6, 3/-1) is no error: alpha = beta = 0.
        for (reward, penalty) in [(1, -6), (3, -1)] {
            let (blocks, _) = context_blocks_ungapped(&[Some(block)], &spec(reward, penalty));
            let blocks = blocks[0].unwrap();
            assert_eq!((blocks.gapped.alpha, blocks.gapped.beta), (0.0, 0.0));
            assert!(context_blocks(&[Some(block)], &spec(reward, penalty)).is_err());
        }
        // No valid context: no blocks, as `context_blocks`.
        assert!(context_blocks_ungapped(&[None, None], &spec(1, -6))
            .0
            .iter()
            .all(Option::is_none));
    }

    /// Oracle (`t2`, 377 x 373, 1/-6 `-ungapped`): the effective search space is 377 x 373 =
    /// 140,621, so the length adjustment of alpha = beta = 0 is 0.
    #[test]
    fn ungapped_blocks_without_a_table_row_have_no_length_adjustment() {
        let spec = NuclScoringSpec {
            reward: 1,
            penalty: -6,
            gap_open: 0,
            gap_extend: 0,
        };
        let block = KarlinParams {
            lambda: 1.39,
            k: 0.749,
            h: 1.38,
            alpha: 0.0,
            beta: 0.0,
        };
        let (blocks, _) = context_blocks_ungapped(&[Some(block)], &spec);
        let adjustment = crate::stats::length_adjustment::compute_length_adjustment_ncbi(
            377,
            373,
            1,
            &blocks[0].unwrap().gapped,
        );
        assert_eq!(adjustment.length_adjustment, 0);
    }

    #[test]
    fn ungapped_option_checks_follow_gapped_calculation() {
        let zero_extension = args(&["-gapopen", "3"]);
        let zero_costs = args(&["-task", "blastn", "-gapopen", "0", "-gapextend", "0"]);
        let counts = ReportCounts::default();
        // "gap extension penalty cannot be 0" only in a gapped search.
        assert!(check_scoring_options_with(&zero_extension, &counts, false).is_err());
        assert!(check_scoring_options_with(&zero_extension, &counts, true).is_ok());
        // "Greedy extension must be used ..." is raised in both.
        assert!(
            message(check_scoring_options_with(&zero_costs, &counts, true))
                .contains("Greedy extension must be used")
        );
        // Gap costs without a table row (gapped: exit 3) and an unsupported pair are accepted.
        for words in [
            &["-reward", "1", "-penalty", "-6"][..],
            &["-reward", "3", "-penalty", "-1"],
            &["-task", "blastn", "-gapopen", "3", "-gapextend", "7"],
            &["-gapopen", "5", "-gapextend", "0"],
            &["-gapopen", "40000", "-gapextend", "40000"],
        ] {
            assert!(run_checks_with(words, &[6], true).is_ok(), "{words:?}");
        }
        // (A gapped search of 1/-6 stops at its table: `context_blocks`,
        // `no_valid_query_needs_no_gapped_block`.)
        // reward <= 0 stays refused.
        assert!(run_checks_with(&["-reward", "0"], &[6], true).is_err());
        // The preliminary hit list of an ungapped search is the hit list: 2147483598 is accepted.
        assert!(run_checks_with(&["-num_alignments", "2147483598"], &[0], true).is_ok());
        assert!(run_checks_with(&["-num_alignments", "2147483598"], &[0], false).is_err());
    }

    // NCBI 2.17.0's greedy offsets (`greedy_xdrop_offset`, x86-64 `Int4`) and the refusals of
    // the X-drops that NCBI cannot run (line RG X1b probe `b2.*`).
    #[test]
    fn gapped_x_drops_that_ncbis_kernels_cannot_run_are_refused() {
        use super::super::alignment::greedy::greedy_xdrop_offset;
        // megablast 1/-2, linear (0/0): the X-drop doubles, (2x + 1) / 6 + 1.
        assert_eq!(greedy_xdrop_offset(13, 1, -2, 0, 0), 5);
        assert_eq!(greedy_xdrop_offset(-1, 1, -2, 0, 0), 1);
        assert_eq!(greedy_xdrop_offset(-3, 1, -2, 0, 0), 1);
        // -8 bits = raw -4: the offset 0 reads `max_score` before NCBI writes it.
        assert_eq!(greedy_xdrop_offset(-4, 1, -2, 0, 0), 0);
        // INT_MIN (past Int4) doubles to 0; 2^30 doubles to INT_MIN; 2^31 - 1 doubles to -2.
        assert_eq!(greedy_xdrop_offset(i32::MIN, 1, -2, 0, 0), 1);
        assert!(greedy_xdrop_offset(1 << 30, 1, -2, 0, 0) <= 0);
        assert_eq!(greedy_xdrop_offset(i32::MAX, 1, -2, 0, 0), 1);
        // An even reward is not doubled: INT_MIN and INT_MAX (whose sum wraps) leave no offset.
        assert_eq!(greedy_xdrop_offset(-1, 2, -3, 0, 0), 1);
        assert!(greedy_xdrop_offset(i32::MIN, 2, -3, 0, 0) <= 0);
        assert!(greedy_xdrop_offset(i32::MAX, 2, -3, 0, 0) <= 0);
        // Affine 1/-2 2/2: doubled 2/4 4/4, common factor gdb3(6, 4, 5) = 1.
        assert_eq!(greedy_xdrop_offset(13, 1, -2, 2, 2), 28);
        assert_eq!(greedy_xdrop_offset(i32::MIN, 1, -2, 2, 2), 2);

        let check = |use_dp, (gapped, finale), (reward, penalty), (gap_open, gap_extend)| {
            check_gapped_x_dropoffs(
                use_dp,
                (25.0, 100.0),
                (gapped, finale),
                (reward, penalty),
                (gap_open, gap_extend),
            )
        };
        let refusal = |result: anyhow::Result<()>| result.expect_err("refused").to_string();
        // Greedy: the defaults, INT_MIN at either stage, and the top of Int4 run.
        assert!(check(false, (13, 54), (1, -2), (0, 0)).is_ok());
        assert!(check(false, (i32::MIN, 54), (1, -2), (0, 0)).is_ok());
        assert!(check(false, (13, i32::MIN), (1, -2), (0, 0)).is_ok());
        assert!(refusal(check(false, (-4, 54), (1, -2), (0, 0)))
            .starts_with("an -xdrop_gap of 25.0 bits"));
        assert!(refusal(check(false, (13, 1 << 30), (1, -2), (0, 0)))
            .starts_with("an -xdrop_gap_final of 100.0 bits"));
        // The final X-drop that the preliminary one set is named -xdrop_gap.
        assert!(refusal(check(false, (1 << 30, 1 << 30), (1, -2), (0, 0)))
            .starts_with("an -xdrop_gap of 25.0 bits"));
        assert!(check(false, (i32::MIN, 54), (2, -3), (0, 0)).is_err());
        // -xdrop_gap -3e9 (raw -1624563704) doubles to 1045839888: its offset 174306649 passes
        // the 19 cells that the final X-drop's 54 gives NCBI's array (NCBI: SIGSEGV). With a
        // final X-drop as large in NCBI's arithmetic (-3e9 or 1e9 bits), NCBI runs it.
        assert!(refusal(check(false, (-1_624_563_704, 54), (1, -2), (0, 0)))
            .starts_with("an -xdrop_gap of 25.0 bits"));
        assert!(check(false, (-1_624_563_704, -1_624_563_704), (1, -2), (0, 0)).is_ok());
        assert!(check(false, (-1_624_563_704, 541_521_234), (1, -2), (0, 0)).is_ok());
        // Dynamic programming: only the final X-drop, against the Int4 sizes of `ALIGN_EX`; the
        // preliminary one is never past the HSP's ungapped score, INT_MIN rises to the gap costs.
        assert!(check(true, (i32::MIN, i32::MIN), (2, -3), (5, 2)).is_ok());
        assert!(check(true, (-100, 2_147_483_441), (2, -3), (5, 2)).is_ok());
        assert!(check(true, (-100, 2_147_483_442), (2, -3), (5, 2)).is_err());
        assert!(check(true, (13, (1 << 30) - 104), (1, -2), (2, 1)).is_ok());
        assert!(refusal(check(true, (13, (1 << 30) - 103), (1, -2), (2, 1)))
            .contains("-xdrop_gap_final of 100.0 bits"));
    }
}
