//! The cutoffs of a TBLASTX query batch's contexts, set up before its subjects are searched: the
//! check of the average subject length, the cutoff score, length adjustment and effective search
//! space of every context, the smallest cutoff, the per-query statistics of the reports, and
//! the Karlin block of the smallest Lambda that the linking cutoffs read.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:887-999
//! (`BLAST_GapAlignSetUp`), which calls blast_setup.c:699-850 (`BLAST_CalcEffLengths`);
//! c++/src/algo/blast/core/blast_parameters.c:279-419 (`BlastInitialWordParametersUpdate`) and
//! :830-994 (`BlastHitSavingParametersUpdate`); c++/src/algo/blast/api/blast_results.cpp:72-115
//! (`CBlastAncillaryData`).
//!
//! This module is a child of `run_impl` for the error type it keeps private
//! (`ErrorAfterBatchReports`). Its function holds statements of `search_query_batch` moved here
//! in their order; what it reads is a struct of references or `Copy` values, and what it makes is
//! what `search_query_batch` then reads under the same names.

use super::{
    cutoff_score_for_update_tblastx, cutoff_score_max_for_tblastx, find_smallest_lambda_params,
    gap_trigger_raw_score, ErrorAfterBatchReports, KarlinParams, QueryContext, QueryFrame,
    SubjectSetNumbers, Subjects, TblastxQueryStats, BLAST_GAP_DECAY_RATE, GAP_TRIGGER_BIT_SCORE,
};
use anyhow::{Context, Result};

/// What `context_cutoffs` reads.
pub(super) struct CutoffInputs<'a> {
    pub(super) contexts_ref: &'a [QueryContext],
    pub(super) numbers: &'a SubjectSetNumbers,
    /// The subject set's letters (`SubjectSetNumbers::total`).
    pub(super) db_length_nucl: i64,
    /// The subject set's sequences (`SubjectSetNumbers::count`).
    pub(super) db_num_seqs: i64,
    pub(super) subjects: Subjects<'a>,
    pub(super) evalue_threshold: f64,
    pub(super) diag_enabled: bool,
    pub(super) debug_cutoffs_all: bool,
    pub(super) query_frames: &'a [Vec<QueryFrame>],
}

/// The cutoffs and statistics of the batch's contexts, which `search_query_batch` reads under the
/// same names.
pub(super) struct ContextCutoffs {
    /// One cutoff score per context of `contexts_ref` (`vec![0; contexts_ref.len()]`, filled).
    pub(super) cutoff_scores: Vec<i32>,
    pub(super) cutoff_score_min: i32,
    pub(super) length_adj_per_context: Vec<i64>,
    pub(super) eff_searchsp_per_context: Vec<i64>,
    pub(super) query_stats: Vec<TblastxQueryStats>,
    pub(super) linking_params_for_cutoff: KarlinParams,
}

/// The average subject length check (an error after the reports of the batches before), then the
/// cutoff score, length adjustment and effective search space of every context, the queries'
/// report statistics and the smallest-Lambda block of the linking cutoffs.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:887-999
/// (`BLAST_GapAlignSetUp`, the average length at :969-980) and :699-850 (`BLAST_CalcEffLengths`);
/// c++/src/algo/blast/core/blast_parameters.c:279-419 (`BlastInitialWordParametersUpdate`),
/// :830-994 (`BlastHitSavingParametersUpdate`) and :91-116 (`s_BlastFindSmallestLambda`).
pub(super) fn context_cutoffs(inputs: &CutoffInputs<'_>) -> Result<ContextCutoffs> {
    let &CutoffInputs {
        contexts_ref,
        numbers,
        db_length_nucl,
        db_num_seqs,
        subjects,
        evalue_threshold,
        diag_enabled,
        debug_cutoffs_all,
        query_frames,
    } = inputs;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:521-528
    // ```c
    //        } else {
    //           ASSERT(sbp->kbp_gap == NULL);
    //           /* for ungapped cases we do not have gbp filled */
    //           if (sbp->gbp) {
    //               sfree(sbp->gbp);
    //               sbp->gbp=NULL;
    //           }
    //        }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:908-939
    // ```c
    //    if (seq_src) {
    //       total_length = BlastSeqSrcGetTotLenStats(seq_src);
    //       if (total_length <= 0)
    //           total_length = BlastSeqSrcGetTotLen(seq_src);
    //       ...
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    //               total_length = -1;
    //               num_seqs = -1;
    //           }
    //           num_seqs = 1;
    //       }
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
    // ```c
    //    if (sbp->gbp) {
    //        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
    //        if (Blast_SubjectIsTranslated(program_number)) {
    //            min_subject_length/=3;
    //        }
    //    } else {
    //        min_subject_length = (Int4) (total_length/num_seqs);
    //    }
    //
    //    if(min_subject_length <=0) {
    // 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_message.c:215-219
    // ```c
    //     case BLASTERR_SUBJECT_LENGTH_INVALID:
    //         new_msg->message = strdup("The average subject length is too short");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:310-314
    // ```c
    //             retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
    //             if (retval) {
    //                 NCBI_THROW(CBlastException, eCoreBlastError,
    //                            BlastErrorCode2String(retval));
    //             }
    // ```
    // The set-up of the search comes after `CheckInternalData` (in `search_query_batch`, before
    // this call): a batch with a valid context stops NCBI when the subjects have fewer letters
    // than subjects (subjects without letters count with none; `AUTHORITY.md` §E3, BI-44).
    // TBLASTX is ungapped, so its Gumbel block is freed and the average length decides. The
    // error comes after the reports of the batches before (`ErrorAfterBatchReports`,
    // `run_in_pool`).
    let (total_length, num_seqs) = if db_length_nucl > 0 {
        (db_length_nucl, db_num_seqs)
    } else {
        (
            if subjects.count() > 0 {
                subjects.first_length() as i64
            } else {
                0
            },
            1,
        )
    };
    // `stats::setup::setup_subject_length`: `(Uint4)(Int4)(total / num_seqs) == 0` (a negative
    // `NSEQ` passes; one of 2^31 fails).
    if crate::stats::setup::setup_subject_length(total_length, num_seqs as i32, None).is_none() {
        return Err(ErrorAfterBatchReports(
            crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: The average subject length is too short\n"
                    .to_string(),
            }
            .into(),
        )
        .into());
    }
    // Precompute per-context cutoff scores using NCBI BLAST algorithm.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:280-419
    //
    // NCBI cutoff calculation for tblastx ungapped path:
    // 1. gap_trigger from ungapped params (kbp_std)
    // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
    // 3. Initial word cutoff: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
    // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
    //
    // The word cutoff uses the average subject nucleotide length. With E = 1e-300 (step 3) its
    // score exceeds gap_trigger whatever the lengths (above 2,000 against about 40), so step 4
    // never takes it and no output depends on this average; it is computed as NCBI computes it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
    // ```c
    //     if ((status = BlastInitialWordParametersNew(program_number,
    //                       word_options, hit_params, lookup_wrap, sbp, query_info,
    //                       BlastSeqSrcGetAvgSeqLen(seq_src), &word_params)) != 0)
    // ```
    // The word cutoffs read the sequence source's average (`Uint4 subj_length`,
    // blast_parameters.c:283): a database's `GetTotalLength / MAX(1, GetNumSeqs)` as NCBI's
    // `Int4` (`SubjectSetNumbers::average`, without the `STATS_*` keys and with a negative
    // `NSEQ` read as 1), a FASTA set's total over its records.
    let subject_len_nucl = match numbers.average {
        Some(average) => i64::from(average as u32),
        None => db_length_nucl / db_num_seqs,
    };
    // NCBI: cutoff scores are stored per query context (no subject-frame dimension).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:320-324
    // ```c
    // BlastUngappedCutoffs *curr_cutoffs = parameters->cutoffs + context;
    // ...
    // curr_cutoffs->cutoff_score = new_cutoff;
    // ```
    let mut cutoff_scores: Vec<i32> = vec![0; contexts_ref.len()];
    // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:401-403
    let mut cutoff_score_min = i32::MAX;

    // =======================================================================
    // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
    // =======================================================================
    // NCBI stores these in query_info->contexts[ctx].length_adjustment and
    // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths
    // (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:700-850).
    // Compute once for the complete subject set and share with all subject jobs.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846-847
    //   query_info->contexts[index].eff_searchsp = effective_search_space;
    //   query_info->contexts[index].length_adjustment = length_adjustment;
    let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
    let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());

    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:324-331;
        // blast_setup.c:774-847: invalid contexts keep zero effective lengths.
        // if (!query_info->contexts[context].is_valid) {
        //     curr_cutoffs->cutoff_score = INT4_MAX; continue;
        // }
        if !ctx.is_valid || ctx.aa_len == 0 {
            cutoff_scores[ctx_idx] = i32::MAX;
            length_adj_per_context.push(0);
            eff_searchsp_per_context.push(numbers.searchsp.unwrap_or(0));
            continue;
        }
        // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
        let ctx_params = &ctx.karlin_params;

        // NCBI: query_length = query_info->contexts[context].query_length
        let query_len_aa = ctx.aa_len as i64;

        // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);

        // =======================================================================
        // NCBI Parity: Use try_compute_eff_lengths_tblastx to get BOTH
        // length_adjustment and eff_searchsp from a single source of truth.
        // This mirrors BLAST_CalcEffLengths which computes and stores both values.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
        // =======================================================================
        // `-dbsize` replaces the database length here only (`SubjectSetNumbers`).
        let eff_lengths = numbers.context_lengths(
            query_len_aa,
            ctx_params, // tblastx uses per-context ungapped params (kbp_gap is NULL)
        )?;
        let eff_searchsp = eff_lengths.eff_searchsp;
        length_adj_per_context.push(eff_lengths.length_adjustment);
        eff_searchsp_per_context.push(eff_searchsp);

        // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
        // This uses the effective search space WITH length adjustment
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:942-946
        let cutoff_score_max = cutoff_score_max_for_tblastx(
            eff_searchsp,
            evalue_threshold, // User's E-value (typically 10.0)
            ctx_params,
        );

        // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
        // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374
        let cutoff = cutoff_score_for_update_tblastx(
            query_len_aa,
            subject_len_nucl, // NUCLEOTIDE length, NOT divided by 3!
            gap_trigger,
            cutoff_score_max,
            BLAST_GAP_DECAY_RATE, // 0.5
            ctx_params,
            1.0, // scale_factor (standard BLOSUM62)
        );

        // DEBUG: Print cutoff values for first context
        static PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if diag_enabled && ctx_idx == 0 && !PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
        {
            eprintln!(
                "[DEBUG CUTOFF] query_len_aa={}, subject_len_nucl={}",
                query_len_aa, subject_len_nucl
            );
            eprintln!("[DEBUG CUTOFF] eff_searchsp={}", eff_searchsp);
            eprintln!(
                "[DEBUG CUTOFF] length_adjustment={}",
                eff_lengths.length_adjustment
            );
            eprintln!("[DEBUG CUTOFF] cutoff_score_max={}", cutoff_score_max);
            eprintln!("[DEBUG CUTOFF] gap_trigger={}", gap_trigger);
            eprintln!("[DEBUG CUTOFF] final cutoff={}", cutoff);
        }
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:943-946
        // ```c
        // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
        // params->cutoffs[context].cutoff_score = new_cutoff;
        // params->cutoffs[context].cutoff_score_max = new_cutoff;
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:360-374
        // ```c
        // BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
        //               MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length),
        //               TRUE, gap_decay_rate);
        // new_cutoff = MIN(new_cutoff, gap_trigger);
        // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
        // ```
        // Diagnostic-only dump of the same per-context values that feed
        // word_params->cutoff_score_min and CalculateLinkHSPCutoffs.
        if debug_cutoffs_all {
            eprintln!(
                "[DEBUG CUTOFF_ALL] ctx_idx={} q_frame={} query_len_aa={} subject_len_nucl={} eff_searchsp={} length_adjustment={} lambda={:.12e} k={:.12e} h={:.12e} cutoff_score_max={} gap_trigger={} word_cutoff={}",
                ctx_idx,
                ctx.frame,
                query_len_aa,
                subject_len_nucl,
                eff_searchsp,
                eff_lengths.length_adjustment,
                ctx_params.lambda,
                ctx_params.k,
                ctx_params.h,
                cutoff_score_max,
                gap_trigger,
                cutoff
            );
        }

        // Track minimum cutoff for linking
        cutoff_score_min = cutoff_score_min.min(cutoff);

        // All subject frames use the same cutoff (NCBI: per-context cutoffs only).
        cutoff_scores[ctx_idx] = cutoff;
    }
    // NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-100
    // ```c++
    //     // find the first valid context corresponding to this query
    //     for (i = 0; i < context_per_query; i++) {
    //         BlastContextInfo *ctx = query_info->contexts +
    //                                 query_number * context_per_query + i;
    //         if (ctx->is_valid) {
    //             m_SearchSpace = ctx->eff_searchsp;
    // 	    m_LengthAdjustment = ctx->length_adjustment;
    //             break;
    //     ...
    //     const int ctx_index = query_number * context_per_query + i;
    //     if (sbp->kbp_std) {
    //         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
    // ```
    // The per-query footer of the reports shows the first valid context's block and
    // search space. The SEG intervals of each frame are the query masks of the report.
    let mut query_stats: Vec<TblastxQueryStats> = query_frames
        .iter()
        .map(|frames| TblastxQueryStats {
            karlin: None,
            eff_searchsp: 0,
            aa_len: 0,
            context0_search_space: 0,
            context0_aa_len: 0,
            seg_masks: frames
                .iter()
                .map(|frame| (frame.frame, frame.seg_masks.clone()))
                .collect(),
        })
        .collect();
    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        let stats = &mut query_stats[ctx.q_idx as usize];
        // Context 0 is the query's frame +1 (invalid contexts hold search space 0).
        if ctx.frame == 1 {
            stats.context0_search_space = eff_searchsp_per_context[ctx_idx];
            if ctx.is_valid {
                stats.context0_aa_len = ctx.aa_len as i64;
            }
        }
        if ctx.is_valid && stats.karlin.is_none() {
            stats.karlin = Some(ctx.karlin_params);
            stats.eff_searchsp = eff_searchsp_per_context[ctx_idx];
            stats.aa_len = ctx.aa_len;
        }
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:401-416
    // ```c
    // if (new_cutoff < cutoff_min) {
    //    cutoff_min = new_cutoff;
    // }
    // parameters->cutoff_score_min = cutoff_min;
    // ```
    if debug_cutoffs_all {
        eprintln!(
            "[DEBUG CUTOFF_ALL] word_params_cutoff_score_min={}",
            cutoff_score_min
        );
    }

    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:101-110
    // if (kbp[index] && query_info->contexts[index].is_valid &&
    //     kbp[index]->Lambda > 0.0 && kbp[index]->Lambda < min_lambda) { ... }
    let context_params: Vec<KarlinParams> = contexts_ref
        .iter()
        .filter(|ctx| ctx.is_valid)
        .map(|ctx| ctx.karlin_params)
        .collect();
    let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
        .context("no valid query statistical parameters")?;
    Ok(ContextCutoffs {
        cutoff_scores,
        cutoff_score_min,
        length_adj_per_context,
        eff_searchsp_per_context,
        query_stats,
        linking_params_for_cutoff,
    })
}
