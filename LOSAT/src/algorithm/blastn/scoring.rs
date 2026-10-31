//! The scoring options of a BLASTN search: NCBI's checks of the given values, and the
//! Karlin-Altschul blocks of the queries.

use super::args::BlastnArgs;
use super::coordination::determine_scoring_params;
use super::ncbi_cutoffs::NCBIMATH_LN2;
use super::pairwise::query_ungapped_karlin_blastna;
use crate::blastinput::report_counts::{run_report_counts, ReportCounts};
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
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:881-906
/// ```c
///            if ( ! ( options->penalty == 0 && options->reward == 0 ) )
///            {
/// 		if (options->penalty >= 0)
/// 		{
/// 			Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
///                             "BLASTN penalty must be negative");
/// 			return BLASTERR_OPTION_VALUE_INVALID;
/// 		}
/// ...
///              if (options->gapped_calculation && options->gap_open > 0 && options->gap_extend == 0)
///              {
///                      Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
///                         "BLASTN gap extension penalty cannot be 0");
///                      return BLASTERR_OPTION_VALUE_INVALID;
///              }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1326-1333
/// ```c
///     } else if (program_number == eBlastTypeBlastn &&
///                options->word_size > DBSEQ_CHUNK_OVERLAP) {
///         char buffer[256];
///         int bytes_written = snprintf(buffer, DIM(buffer),
///                   "Word-size must be less than or equal to %d", DBSEQ_CHUNK_OVERLAP);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1247-1261,1399-1413
/// ```c
/// static Boolean
/// s_DiscWordOptionsValidate(Int4 word_size, Uint1 template_length,
///                           Uint1 template_type,
///                           Blast_Message** blast_msg)
/// {
///    if (template_length == 0)
///       return TRUE;
///
///
///    if (word_size != 11 && word_size != 12) {
///       Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///                          "Invalid discontiguous template parameters: word "
///                          "size must be either 11 or 12");
///       return FALSE;
///    }
/// ...
///     if (/*program_number == eBlastTypeBlastn &&*/
///         Blast_ProgramIsNucleotide(program_number) &&
///         !Blast_QueryIsPattern(program_number) &&
///         options->mb_template_length > 0) {
///       if (!s_DiscWordOptionsValidate(options->word_size,
///               options->mb_template_length,
///               options->mb_template_type,
///               blast_msg)) {
///          return BLASTERR_OPTION_VALUE_INVALID;
///       } else if (options->lut_type != eMBLookupTable) {
///          Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///             "Invalid lookup table type for discontiguous Mega BLAST");
///          return BLASTERR_OPTION_VALUE_INVALID;
///       }
///    }
/// ```
/// The template length and type are 16, 18 or 21 and 0, 1 or 2 by the argument
/// constraints (`-template_length`, `-template_type`), so their two checks cannot fail.
/// The lookup table type of the options is the task's: megablast and dc-megablast
/// `eMBLookupTable`, blastn and blastn-short `eNaLookupTable`
/// (blast_nucl_options.cpp:137-150, disc_nucl_options.cpp:55-64).
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1518-1523
/// ```c
/// 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
/// 	{
/// 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///          "expect value or cutoff score must be greater than zero");
/// 		return BLASTERR_OPTION_VALUE_INVALID;
/// 	}
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1699-1711
/// ```c
///     if (program_number == eBlastTypeBlastn)
///     {
///         if (score_options->gap_open == 0 && score_options->gap_extend == 0)
///         {
///             if (ext_options->ePrelimGapExt != eGreedyScoreOnly &&
///                 ext_options->eTbackExt != eGreedyTbck)
///                 {
///                     Blast_MessageWrite(blast_msg, eBlastSevWarning,
///                                        kBlastMessageNoContext,
///                                        "Greedy extension must be used if gap existence and extension options are zero");
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1759-1776 (BLAST_ValidateOptions
/// runs the scoring, lookup table, hit saving and extension-scoring checks in this order)
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3636-3639
/// ```c
///     try { retval->Validate(); }
///     catch (const CBlastException& e) {
///         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
///     }
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:172-175
/// ```c
///     catch (const blast::CInputException& e) {                               \
///         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
///         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
///         exit_code = BLAST_INPUT_ERROR;                                      \
/// ```
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
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:900-904
        // ```c
        //  if (options->gapped_calculation && options->gap_open > 0 && options->gap_extend == 0)
        //  {
        //      Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
        //         "BLASTN gap extension penalty cannot be 0");
        // ```
        "BLASTN gap extension penalty cannot be 0".to_string()
    } else if word_size > DBSEQ_CHUNK_OVERLAP {
        format!("Word-size must be less than or equal to {DBSEQ_CHUNK_OVERLAP}")
    } else if template_length > 0 && word_size != 11 && word_size != 12 {
        "Invalid discontiguous template parameters: word size must be either 11 or 12".to_string()
    } else if template_length > 0 && !mb_lookup_options {
        "Invalid lookup table type for discontiguous Mega BLAST".to_string()
    } else if args.off_diagonal_range.is_some_and(|range| range != 0)
        && super::coordination::determine_window_size(args) == 0
    {
        // NCBI reference: c++/src/algo/blast/core/blast_options.c:1583-1596
        // (`BlastInitialWordOptionsValidate`, after the lookup table's checks and before the
        // hit saving's)
        // ```c
        //    if (program_number == eBlastTypeBlastn &&
        //        options->scan_range && !options->window_size)
        //    {
        //       Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //                             "off_diagonal_range is only useful in 2-hit algorithm");
        // ```
        "off_diagonal_range is only useful in 2-hit algorithm".to_string()
    } else if counts.hitlist_size < 1 {
        // NCBI blast_options.c:1511-1517, `BlastHitSavingOptionsValidate`, whose hit list
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
    Err(NativeError {
        exit: 1,
        message: format!(
            "BLAST query/options error: {message}\nPlease refer to the BLAST+ user manual.\n"
        ),
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
    let hitlist_size = counts.hitlist_size;
    if super::hsp::get_prelim_hitlist_size(hitlist_size, false, !ungapped) < 1 {
        let what = if args.max_target_seqs.is_some() {
            format!("a -max_target_seqs of {hitlist_size}")
        } else {
            format!("a hit list size of {hitlist_size} (-num_descriptions / -num_alignments)")
        };
        anyhow::bail!(
            "{what}, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTN"
        );
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1506-1520
    // ```c
    //     sbp->loscore = BLAST_SCORE_MAX;
    //     sbp->hiscore = BLAST_SCORE_MIN;
    //     matrix = sbp->matrix->data;
    //     for (index1=0; index1<sbp->alphabet_size; index1++)
    //     {
    //       for (index2=0; index2<sbp->alphabet_size; index2++)
    //       {
    //          score = matrix[index1][index2];
    //          if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
    //             continue;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2175-2179
    // ```c
    //          score = matrix[index1][index2];
    //          if (score >= sbp->loscore)
    //          {
    //             sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
    //          }
    // ```
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

/// NCBI reference: c++/include/algo/blast/core/blast_hits.h:192
/// ```c
/// #define DBSEQ_CHUNK_OVERLAP 100
/// ```
const DBSEQ_CHUNK_OVERLAP: usize = 100;

/// The largest gap cost that LOSAT's BLASTN accepts with greedy extension.
const MAX_GREEDY_GAP_COST: i32 = 32767;

/// NCBI's error for scores or gap costs without Karlin-Altschul values, raised when the
/// first batch with a valid query is set up: every query of that batch gets the message.
///
/// NCBI reference: c++/src/algo/blast/api/blast_aux_priv.cpp:94-102
/// ```c
///             // applies to all queries
///             CRef<CSearchMessage> sm(new CSearchMessage(blmsg->severity,
///                                                        kBlastMessageNoContext,
///                                                        msg));
///             NON_CONST_ITERATE(TSearchMessages, query_messages, messages) {
///                 query_messages->push_back(sm);
///             }
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:1013-1024
/// ```c
/// TSearchMessages::ToString() const
/// {
///     string retval;
///     ITERATE(vector<TQueryMessages>, qm, *this) {
///         if (qm->empty()) {
///             continue;
///         }
///         ITERATE(TQueryMessages, msg, *qm) {
///             retval += (*msg)->GetMessage() + " ";
///         }
/// ```
/// NCBI reference: c++/include/algo/blast/api/blast_types.hpp:180-183
/// ```c
///     string GetMessage(bool withSeverity = true) const
///     {
///     	if (withSeverity) {
///     		return GetSeverityString() + ": " + m_Message;
/// ```
/// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:170-184
/// ```c
///     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, search_messages);
///     if (status != 0 &&
///     ...
///             		msg = search_messages.ToString();
///     ...
///                 NCBI_THROW(CBlastException, eCoreBlastError, msg);
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
pub(crate) fn karlin_error(message: &str, batch_queries: usize) -> anyhow::Error {
    NativeError {
        exit: 3,
        message: format!(
            "BLAST engine error: {}\n",
            format!("Error: {message} ").repeat(batch_queries)
        ),
    }
    .into()
}

/// The hit list size and the counts of the pairwise report of a run that writes the output
/// formats `outfmts` (`-max_target_seqs`, `-num_descriptions`, `-num_alignments` of `args`),
/// with the formatting handler's warnings: the one derivation of `process_options` and of
/// `check_scoring`.
pub(crate) fn report_counts_of(
    args: &BlastnArgs,
    outfmts: &[i32],
) -> anyhow::Result<(ReportCounts, Vec<u8>)> {
    run_report_counts(
        "blastn",
        outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
    )
}

/// The checks of `check_scoring_options`, `check_losat_limits` and the Karlin-Altschul
/// tables, without the queries: the checks that a host can run before a search. A table
/// error is reported as for a batch of one query.
pub fn check_scoring(args: &BlastnArgs) -> anyhow::Result<()> {
    check_scoring_for_formats(args, &[0])
}

/// `check_scoring` for a search that writes the output formats `outfmts` (formats of one
/// hit list size): the hit list size is the one that the run of these formats derives. A
/// host that runs one search per hit list size checks each of them.
pub fn check_scoring_for_formats(args: &BlastnArgs, outfmts: &[i32]) -> anyhow::Result<()> {
    check_scoring_for_formats_with(args, outfmts, args.ungapped)
}

/// `check_scoring_for_formats` for an ungapped (`-ungapped`) search when `ungapped`: NCBI
/// builds no gapped Karlin-Altschul block then, so the gap costs are not looked up in the
/// tables and greedy gap costs are not limited.
///
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:513-529
/// ```c
///     status = Blast_ScoreBlkKbpUngappedCalc(program_number, sbp, query_blk->sequence,
///                                            query_info, blast_message);
///     ...
///     if (scoring_options->gapped_calculation) {
///         status = s_BlastScoreBlkGappedFill(sbp, scoring_options, program_number,
///                                            query_info, blast_message);
/// ```
pub(crate) fn check_scoring_for_formats_with(
    args: &BlastnArgs,
    outfmts: &[i32],
    ungapped: bool,
) -> anyhow::Result<()> {
    // The host has no database files (`run_local`).
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // `args.report_counts` is set by `process_options` of a run; a host that has not run
    // it derives the counts the same way, for the output formats of the search.
    let (counts, _) = report_counts_of(args, outfmts)?;
    check_scoring_options_with(args, &counts, ungapped)?;
    check_losat_limits_with(args, &counts, ungapped)?;
    if ungapped {
        return Ok(());
    }
    let spec = scoring_spec(args);
    NuclValues::new(spec.reward, spec.penalty)
        .and_then(|values| values.check_gaps(&spec))
        .map_err(|message| karlin_error(&message, 1))?;
    check_greedy_gap_costs(args)
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
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2762-2784
/// ```c
///    for (context = query_info->first_context;
///         context <= query_info->last_context; ++context) {
///     ...
///       query_length = contexts[context].query_length;
///       context_offset = contexts[context].query_offset;
///       buffer = &query[context_offset];
///
///       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
///       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
///       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
///       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
///       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
/// ```
/// The two strands have the same composition in another order, so their blocks can
/// differ in the last bits.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2771-2772
/// ```c
///       if ( !contexts[context].is_valid )
///           continue;
/// ```
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
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:73-107
/// ```c
///     for (index = query_info->first_context;
///          index <= query_info->last_context; index++) {
///
///         if ( !query_info->contexts[index].is_valid ) {
///             continue;
///         }
///
///         sbp->kbp_gap_std[index] = Blast_KarlinBlkNew();
///     ...
///               retval =
///                 Blast_KarlinBlkNuclGappedCalc(sbp->kbp_gap_std[index],
///                     scoring_options->gap_open, scoring_options->gap_extend,
///                     scoring_options->reward, scoring_options->penalty,
///                     sbp->kbp_std[index], &(sbp->round_down), error_return);
/// ```
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
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:513-529
/// ```c
///     if (scoring_options->gapped_calculation) { ... s_BlastScoreBlkGappedFill ... }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:701-705, 768-848
/// ```c
///     double alpha = 0, beta = 0; /* alpha and beta for new scoring system */
///     ...
///     kbp_ptr = sbp->kbp;   (ungapped)
///     ...
///             Blast_GetNuclAlphaBeta(scoring_options->reward, scoring_options->penalty,
///                 scoring_options->gap_open, scoring_options->gap_extend,
///                 sbp->kbp_std[index], scoring_options->gapped_calculation, &alpha, &beta);
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3989-3992, 4024-4029
/// ```c
///     if (status)
///        return status;       (s_GetNuclValuesArray: no table row for the pair; alpha, beta unset)
///     ...
///     if (!found) {
///         *alpha = kbp->Lambda/kbp->H;
///         *beta = s_GetUngappedBeta(reward, penalty);
///     }
/// ```
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
/// NCBI reference: c++/src/algo/blast/core/blast_parameters.c:92-116
/// ```c
///     double min_lambda = (double) INT4_MAX;
///     ...
///     for (i=query_info->first_context; i<=query_info->last_context; i++) {
///         ...
///         if (s_BlastKarlinBlkIsValid(kbp_in[i])) {
///             if (min_lambda > kbp_in[i]->Lambda)
///             {
///                 min_lambda = kbp_in[i]->Lambda;
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_parameters.c:455-463
/// ```c
///       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
///       params->gap_x_dropoff = (Int4)
///           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
///     ...
///       params->gap_x_dropoff_final = (Int4)
///           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
/// ```
/// Gapped blocks from the tables are the same for every context; blocks copied from the
/// ungapped ones differ with the composition, so the X-drops depend on the contexts of the
/// batch. Without a valid context the batch is not searched.
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
    let gapped = (x_drop_gapped_bits * NCBIMATH_LN2 / min_lambda) as i32;
    let finale = (x_drop_final_bits * NCBIMATH_LN2 / min_lambda).max(gapped as f64) as i32;
    (gapped, finale)
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

    // NCBI reference: c++/src/algo/blast/core/blast_options.c:1247-1261,1399-1413
    // ```c
    //    if (word_size != 11 && word_size != 12) {
    //       Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "Invalid discontiguous template parameters: word "
    //                          "size must be either 11 or 12");
    // ...
    //       } else if (options->lut_type != eMBLookupTable) {
    //          Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //             "Invalid lookup table type for discontiguous Mega BLAST");
    // ```
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
            "BLAST query/options error: off_diagonal_range is only useful in 2-hit algorithm\nPlease refer to the BLAST+ user manual.\n"
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

    // NCBI blast_options.c:1518: the e-value check needs a cutoff score of 0 or below too.
    #[test]
    fn min_raw_gapped_score_above_zero_allows_an_evalue_of_zero() {
        const REFUSED: &str = "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n";
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
            "BLAST query/options error: BLASTN penalty must be negative\nPlease refer to the BLAST+ user manual.\n"
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

    #[test]
    fn check_scoring_derives_the_counts_a_run_derives() {
        // `args.report_counts` stays at its default until a run sets it, so a host's
        // `check_scoring` must derive the counts from the options itself.
        let error = message(check_scoring(&args(&[
            "-num_descriptions",
            "0",
            "-num_alignments",
            "0",
        ])));
        assert_eq!(
            error,
            "BLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
        );
        // The hit list of format 0 is the larger of the two numbers.
        assert!(check_scoring(&args(&["-num_alignments", "0"])).is_ok());
        assert!(check_scoring(&args(&["-num_descriptions", "0"])).is_ok());
        assert!(check_scoring(&args(&["-num_descriptions", "1", "-num_alignments", "0"])).is_ok());
        // The preliminary hit list size overflows.
        let error = check_scoring(&args(&["-num_alignments", "2147483598"]))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("(-num_descriptions / -num_alignments)")
                && error.contains("not supported by LOSAT's BLASTN"),
            "{error}"
        );
        let error = check_scoring(&args(&["-num_descriptions", "2147483598"]))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("whose preliminary hit list size overflows"),
            "{error}"
        );
        // `-max_target_seqs` keeps its own wording (the default path, before O1b).
        let error = check_scoring(&args(&["-max_target_seqs", "2147483598"]))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("a -max_target_seqs of 2147483598")
                && error.contains("whose preliminary hit list size overflows"),
            "{error}"
        );
        // The same options set on the args by a run give the same refusals.
        let mut run_args = args(&["-num_descriptions", "0", "-num_alignments", "0"]);
        run_args.report_counts = report_counts_of(&run_args, &[0]).unwrap().0;
        assert!(check_scoring_options(&run_args)
            .unwrap_err()
            .to_string()
            .contains("No hits are being saved"));
    }

    #[test]
    fn table_errors_repeat_for_each_query_of_the_batch() {
        let error = karlin_error("Substitution scores 1 and -6 are not supported", 2);
        assert_eq!(
            error.downcast::<NativeError>().unwrap().message,
            "BLAST engine error: Error: Substitution scores 1 and -6 are not supported Error: Substitution scores 1 and -6 are not supported \n"
        );
        assert!(
            message(check_scoring(&args(&["-reward", "1", "-penalty", "-6"])))
                .starts_with("BLAST engine error: Error: Substitution scores 1 and -6")
        );
        assert!(check_scoring(&args(&["-gapopen", "5", "-gapextend", "2"])).is_ok());
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

    // NCBI reference: c++/src/algo/blast/core/blast_stat.c:3955-3963, 3989-3992, 4024-4029
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
        let adjustment =
            crate::stats::compute_length_adjustment_ncbi(377, 373, 1, &blocks[0].unwrap().gapped);
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
            let parsed = args(words);
            assert!(
                check_scoring_for_formats_with(&parsed, &[6], true).is_ok(),
                "{words:?}"
            );
        }
        assert!(check_scoring_for_formats_with(
            &args(&["-reward", "1", "-penalty", "-6"]),
            &[6],
            false
        )
        .is_err());
        // reward <= 0 stays refused.
        assert!(check_scoring_for_formats_with(&args(&["-reward", "0"]), &[6], true).is_err());
        // The preliminary hit list of an ungapped search is the hit list: 2147483598 is accepted.
        assert!(check_scoring_for_formats_with(
            &args(&["-num_alignments", "2147483598"]),
            &[0],
            true
        )
        .is_ok());
        assert!(check_scoring_for_formats_with(
            &args(&["-num_alignments", "2147483598"]),
            &[0],
            false
        )
        .is_err());
    }
}
