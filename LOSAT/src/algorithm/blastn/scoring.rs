//! The scoring options of a BLASTN search: NCBI's checks of the given values, and the
//! Karlin-Altschul blocks of the queries.

use super::args::BlastnArgs;
use super::coordination::determine_scoring_params;
use super::ncbi_cutoffs::NCBIMATH_LN2;
use super::pairwise::query_ungapped_karlin;
use crate::cli::NativeError;
use crate::config::NuclScoringSpec;
use crate::stats::{KarlinParams, NuclValues};

fn scoring_spec(args: &BlastnArgs) -> NuclScoringSpec {
    let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(args);
    NuclScoringSpec {
        reward,
        penalty,
        gap_open,
        gap_extend,
    }
}

/// NCBI's checks of the options, before the queries are read, then the limits of LOSAT.
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
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3609-3612
/// ```c
///         try { m_OptsHandle->Validate(); }
///         catch (const CBlastException& e) {
///             NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
///         }
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:172-175
/// ```c
///     catch (const blast::CInputException& e) {                               \
///         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
///         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
///         exit_code = BLAST_INPUT_ERROR;                                      \
/// ```
/// Only megablast extends greedily. The reward and the penalty are NCBI's 16-bit values
/// (`determine_scoring_params`).
///
/// LOSAT's limits come after NCBI's checks, for options that NCBI runs: a reward of 0 or
/// less (NCBI's rmblastn matrix scoring when the penalty is 0 too, otherwise no valid
/// query),
/// an infinite or NaN e-value, and scores whose range (reward - penalty) is above
/// `MAX_SCORE_RANGE`, where LOSAT's Karlin-Altschul computation has not been compared with
/// NCBI's; that range limit comes before NCBI's Karlin-Altschul table error, which needs
/// that computation. The limit on greedy gap costs follows the table check
/// (`check_greedy_gap_costs`).
pub fn check_scoring_options(args: &BlastnArgs) -> anyhow::Result<()> {
    let spec = scoring_spec(args);
    let word_size = super::coordination::determine_effective_word_size(args);
    let greedy = args.task == "megablast";
    let matrix_only = spec.reward == 0 && spec.penalty == 0;
    let message = if !matrix_only && spec.penalty >= 0 {
        "BLASTN penalty must be negative".to_string()
    } else if spec.gap_open > 0 && spec.gap_extend == 0 {
        "BLASTN gap extension penalty cannot be 0".to_string()
    } else if word_size > DBSEQ_CHUNK_OVERLAP {
        format!("Word-size must be less than or equal to {DBSEQ_CHUNK_OVERLAP}")
    } else if args.evalue <= 0.0 {
        "expect value or cutoff score must be greater than zero".to_string()
    } else if spec.gap_open == 0 && spec.gap_extend == 0 && !greedy {
        "Greedy extension must be used if gap existence and extension options are zero".to_string()
    } else {
        if spec.reward <= 0 {
            anyhow::bail!(
                "a reward of {} (NCBI BLAST+'s 16-bit value; 0 is NCBI's matrix scoring of rmblastn, otherwise no query is valid) is not supported by LOSAT's BLASTN",
                spec.reward
            );
        }
        if !args.evalue.is_finite() {
            anyhow::bail!("an infinite or NaN e-value is not supported by LOSAT's BLASTN");
        }
        let range = i64::from(spec.reward) - i64::from(spec.penalty);
        if range > MAX_SCORE_RANGE {
            anyhow::bail!(
                "reward {} and penalty {} span more than {MAX_SCORE_RANGE} score units, which is not supported by LOSAT's BLASTN",
                spec.reward,
                spec.penalty
            );
        }
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

/// The largest score range (reward - penalty) that LOSAT's BLASTN accepts; pairs up to
/// 1000/-2000 were compared with NCBI (S07+ independent audit).
const MAX_SCORE_RANGE: i64 = 3000;

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

/// The check of `check_scoring_options` and of the Karlin-Altschul tables, without the
/// queries: the checks that a host can run before a search. A table error is reported
/// as for a batch of one query.
pub fn check_scoring(args: &BlastnArgs) -> anyhow::Result<()> {
    check_scoring_options(args)?;
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
    /// The gapped block (NCBI `kbp_gap`) with the alpha and beta of the length adjustment:
    /// e-values, bit scores and cutoffs.
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
pub(crate) fn context_ungapped_blocks<'a>(
    contexts: impl Iterator<Item = &'a [u8]>,
    spec: &NuclScoringSpec,
) -> Vec<Option<KarlinParams>> {
    contexts
        .map(|context| query_ungapped_karlin(context, spec.reward, spec.penalty))
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
/// ungapped ones differ with the composition. LOSAT knows only the first batch
/// (`first_query_batch`), so contexts whose blocks give different X-drops are exact only
/// when the first batch holds every query (`one_batch`); otherwise this is the reason
/// for failing.
pub(crate) fn gap_x_dropoffs(
    context_karlin: &[Option<ContextKarlin>],
    x_drop_gapped_bits: i32,
    x_drop_final_bits: i32,
    one_batch: bool,
) -> Result<(i32, i32), String> {
    let x_drops = |lambda: f64| {
        let gapped = (x_drop_gapped_bits as f64 * NCBIMATH_LN2 / lambda) as i32;
        let finale = (x_drop_final_bits as f64 * NCBIMATH_LN2 / lambda).max(gapped as f64) as i32;
        (gapped, finale)
    };
    let lambdas: Vec<f64> = context_karlin
        .iter()
        .flatten()
        .map(|blocks| blocks.gapped.lambda)
        .collect();
    let Some(&first) = lambdas.first() else {
        return Ok(x_drops(KarlinParams::default().lambda));
    };
    if lambdas
        .iter()
        .all(|&lambda| x_drops(lambda) == x_drops(first))
    {
        return Ok(x_drops(first));
    }
    if !one_batch {
        return Err("the gapped X-drop of these gap costs depends on the queries of each of NCBI BLAST+'s adaptive query batches, which LOSAT does not reproduce (the queries do not fit in the first batch)".to_string());
    }
    Ok(x_drops(
        lambdas.into_iter().fold(f64::from(i32::MAX), f64::min),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(words: &[&str]) -> BlastnArgs {
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

    fn message(result: anyhow::Result<()>) -> String {
        result
            .expect_err("rejected")
            .downcast::<NativeError>()
            .expect("an NCBI error")
            .message
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
}
