//! The cutoffs of a TBLASTX query batch's contexts, set up before its subjects are searched: the
//! check of the average subject length, the cutoff score, length adjustment and effective search
//! space of every context, the smallest cutoff, the per-query statistics of the reports, and
//! the Karlin block of the smallest Lambda that the linking cutoffs read; and the statistics each
//! subject reads (`SubjectStatistics`): the batch's, or the subject's own where NCBI computes the
//! statistics of each subject (`per_subject_cutoffs`, `SubjectStatistics::of_subject`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:887-999
//! (`BLAST_GapAlignSetUp`), which calls blast_setup.c:699-850 (`BLAST_CalcEffLengths`);
//! c++/src/algo/blast/core/blast_parameters.c:279-419 (`BlastInitialWordParametersUpdate`) and
//! :830-994 (`BlastHitSavingParametersUpdate`); c++/src/algo/blast/api/blast_results.cpp:72-115
//! (`CBlastAncillaryData`); c++/src/algo/blast/core/blast_engine.c:1434-1444 (the statistics of
//! each subject, `BLAST_OneSubjectUpdateParameters`).
//!
//! This module is a child of `run_impl` for the error type it keeps private
//! (`ErrorAfterBatchReports`). Its function holds statements of `search_query_batch` moved here
//! in their order; what it reads is a struct of references or `Copy` values, and what it makes is
//! what `search_query_batch` then reads under the same names. `per_subject_cutoffs` and
//! `SubjectStatistics::of_subject` are statements of `process_subject` moved here.

use super::word_finder;
use super::{
    cutoff_score_for_update_tblastx, cutoff_score_max_for_tblastx, find_smallest_lambda_params,
    gap_trigger_raw_score, ErrorAfterBatchReports, KarlinParams, QueryContext, QueryFrame,
    SubjectSetNumbers, Subjects, TblastxQueryStats, BLAST_GAP_DECAY_RATE, GAP_TRIGGER_BIT_SCORE,
};
use crate::wording::ENGINE_ERROR;
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
    /// The preliminary stage's option spaces of the contexts, where they differ from the setup's
    /// or the statistics of each subject read them (`SubjectSetNumbers::preliminary_spaces`).
    pub(super) option_spaces: Option<Vec<i64>>,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:521-528
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-939
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-219
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
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
    // A count of 0 (of a positive total) is LOSAT's refusal before (`search_query_batch`). Else
    // `stats::setup::setup_subject_length`: `(Uint4)(Int4)(total / num_seqs) == 0` (a negative
    // `NSEQ` passes; one of 2^31 fails).
    if crate::stats::setup::setup_subject_length(total_length, num_seqs as i32, None).is_none() {
        return Err(ErrorAfterBatchReports(
            crate::cli::NativeError {
                exit: 3,
                message: format!(
                    "{ENGINE_ERROR} {}\n",
                    crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT
                )
                .into(),
            }
            .into(),
        )
        .into());
    }
    // Precompute per-context cutoff scores using NCBI BLAST algorithm.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:280-419
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
    // The word cutoffs read the sequence source's average (`Uint4 subj_length`,
    // blast_parameters.c:283): a database's `GetTotalLength / MAX(1, GetNumSeqs)` as NCBI's
    // `Int4` (`SubjectSetNumbers::average`, without the `STATS_*` keys and with a negative
    // `NSEQ` read as 1), a FASTA set's total over its records.
    let subject_len_nucl = match numbers.average {
        Some(average) => i64::from(average as u32),
        None => db_length_nucl / db_num_seqs,
    };
    // NCBI: cutoff scores are stored per query context (no subject-frame dimension).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:320-324
    let mut cutoff_scores: Vec<i32> = vec![0; contexts_ref.len()];
    // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:401-403
    let mut cutoff_score_min = i32::MAX;

    // =======================================================================
    // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
    // =======================================================================
    // NCBI stores these in query_info->contexts[ctx].length_adjustment and
    // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths
    // (c++/src/algo/blast/core/blast_setup.c:700-850).
    // Compute once for the complete subject set and share with all subject jobs.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:846-847
    let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
    let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220 (cited at
    // `stats::setup::preliminary_option_space`): the preliminary stage's option spaces, which
    // its setup reads in place of its own where they differ (`SubjectSetNumbers::
    // preliminary_spaces_differ`), and the statistics of each subject read (`subject_cutoffs`).
    let option_spaces = numbers.preliminary_spaces(contexts_ref.iter().map(|ctx| {
        let query_length = if ctx.is_valid { ctx.aa_len as i64 } else { 0 };
        (query_length, &ctx.karlin_params)
    }))?;
    let preliminary_spaces = option_spaces
        .as_deref()
        .filter(|_| numbers.preliminary_spaces_differ());
    // The setup's search spaces, which the reports print (`query_stats` below), where the
    // preliminary stage's differ from them; else they are `eff_searchsp_per_context`.
    let mut report_searchsp_per_context: Option<Vec<i64>> =
        preliminary_spaces.map(|_| Vec::with_capacity(contexts_ref.len()));

    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:324-331;
        // blast_setup.c:774-847: invalid contexts keep zero effective lengths.
        if !ctx.is_valid || ctx.aa_len == 0 {
            cutoff_scores[ctx_idx] = i32::MAX;
            length_adj_per_context.push(0);
            eff_searchsp_per_context.push(
                preliminary_spaces.map_or(numbers.searchsp.unwrap_or(0), |spaces| spaces[ctx_idx]),
            );
            if let Some(report) = report_searchsp_per_context.as_mut() {
                report.push(numbers.searchsp.unwrap_or(0));
            }
            continue;
        }
        // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
        let ctx_params = &ctx.karlin_params;

        // NCBI: query_length = query_info->contexts[context].query_length
        let query_len_aa = ctx.aa_len as i64;

        // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-345
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);

        // =======================================================================
        // NCBI Parity: Use try_compute_eff_lengths_tblastx to get BOTH
        // length_adjustment and eff_searchsp from a single source of truth.
        // This mirrors BLAST_CalcEffLengths which computes and stores both values.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-847
        // =======================================================================
        // `-dbsize` replaces the database length here only (`SubjectSetNumbers`).
        let eff_lengths = numbers.context_lengths(
            query_len_aa,
            ctx_params, // tblastx uses per-context ungapped params (kbp_gap is NULL)
        )?;
        if let Some(report) = report_searchsp_per_context.as_mut() {
            report.push(eff_lengths.eff_searchsp);
        }
        // The preliminary stage's lengths, with its option space where it differs.
        let eff_lengths = match preliminary_spaces {
            Some(spaces) => {
                numbers.context_lengths_with(query_len_aa, ctx_params, spaces[ctx_idx])?
            }
            None => eff_lengths,
        };
        let eff_searchsp = eff_lengths.eff_searchsp;
        length_adj_per_context.push(eff_lengths.length_adjustment);
        eff_searchsp_per_context.push(eff_searchsp);

        // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
        // This uses the effective search space WITH length adjustment
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:942-946
        let cutoff_score_max = cutoff_score_max_for_tblastx(
            eff_searchsp,
            evalue_threshold, // User's E-value (typically 10.0)
            ctx_params,
        );

        // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
        // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-374
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:943-946
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:360-374
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-100
    // The per-query footer of the reports shows the first valid context's block and
    // search space. The SEG intervals of each frame are the query masks of the report.
    let mut query_stats: Vec<TblastxQueryStats> = query_frames
        .iter()
        .map(|frames| TblastxQueryStats {
            karlin: None,
            eff_searchsp: 0,
            aa_len: 0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:157-161
            // (the request's search space is `-searchsp` when it was given)
            // A query without frame +1 among its contexts (`-strand minus`, or fewer than three
            // letters) keeps this value: `-searchsp`, else 0 (NCBI's context 0 then has no
            // length, or is not even visited, blast_setup.c:770-772 with `first_context` 3).
            context0_search_space: numbers.searchsp.unwrap_or(0),
            context0_aa_len: 0,
            seg_masks: frames
                .iter()
                .map(|frame| (frame.frame, frame.seg_masks.clone()))
                .collect(),
        })
        .collect();
    let report_searchsp = report_searchsp_per_context
        .as_deref()
        .unwrap_or(&eff_searchsp_per_context);
    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        let stats = &mut query_stats[ctx.q_idx as usize];
        // Context 0 is the query's frame +1 (invalid contexts hold search space 0).
        if ctx.frame == 1 {
            stats.context0_search_space = report_searchsp[ctx_idx];
            if ctx.is_valid {
                stats.context0_aa_len = ctx.aa_len as i64;
            }
        }
        if ctx.is_valid && stats.karlin.is_none() {
            stats.karlin = Some(ctx.karlin_params);
            stats.eff_searchsp = report_searchsp[ctx_idx];
            stats.aa_len = ctx.aa_len;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:401-416
    if debug_cutoffs_all {
        eprintln!(
            "[DEBUG CUTOFF_ALL] word_params_cutoff_score_min={}",
            cutoff_score_min
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:101-110
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
        option_spaces,
    })
}

/// The cutoffs, length adjustments and search spaces of the batch's contexts for the subject of
/// `subject_len_nucl` bases, where NCBI's preliminary search computes the statistics of each
/// subject (`SubjectSetNumbers::per_subject`): `BLAST_OneSubjectUpdateParameters` with the
/// preliminary stage's option spaces `option_spaces` (`ContextCutoffs::option_spaces`). The word
/// cutoff of a context is `context_cutoffs`' with the subject's length in place of the average.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1001-1024 (cited at
/// `stats::setup::subject_eff_lengths`); c++/src/algo/blast/core/blast_parameters.c:279-419
/// (`BlastInitialWordParametersUpdate`, `subj_length`) and :830-994
/// (`BlastHitSavingParametersUpdate`: an ungapped search without the Gumbel block).
pub(super) fn subject_cutoffs(
    contexts_ref: &[QueryContext],
    numbers: &SubjectSetNumbers,
    option_spaces: &[i64],
    subject_len_nucl: i64,
    evalue_threshold: f64,
) -> Result<SubjectCutoffs, crate::stats::setup::SearchSpaceOverflow> {
    let mut cutoffs = SubjectCutoffs {
        cutoff_scores: vec![i32::MAX; contexts_ref.len()],
        cutoff_score_min: i32::MAX,
        length_adj_per_context: vec![0; contexts_ref.len()],
        eff_searchsp_per_context: option_spaces.to_vec(),
    };
    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        if !ctx.is_valid || ctx.aa_len == 0 {
            continue;
        }
        let ctx_params = &ctx.karlin_params;
        let query_len_aa = ctx.aa_len as i64;
        let eff_lengths = numbers.subject_context_lengths(
            query_len_aa,
            subject_len_nucl,
            option_spaces[ctx_idx],
            ctx_params,
        )?;
        let cutoff_score_max =
            cutoff_score_max_for_tblastx(eff_lengths.eff_searchsp, evalue_threshold, ctx_params);
        let cutoff = cutoff_score_for_update_tblastx(
            query_len_aa,
            subject_len_nucl,
            gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params),
            cutoff_score_max,
            BLAST_GAP_DECAY_RATE,
            ctx_params,
            1.0,
        );
        cutoffs.cutoff_scores[ctx_idx] = cutoff;
        cutoffs.cutoff_score_min = cutoffs.cutoff_score_min.min(cutoff);
        cutoffs.length_adj_per_context[ctx_idx] = eff_lengths.length_adjustment;
        cutoffs.eff_searchsp_per_context[ctx_idx] = eff_lengths.eff_searchsp;
    }
    Ok(cutoffs)
}

/// The values of `subject_cutoffs`, which the subject's search reads under the names of
/// `ContextCutoffs`.
pub(super) struct SubjectCutoffs {
    pub(super) cutoff_scores: Vec<i32>,
    pub(super) cutoff_score_min: i32,
    pub(super) length_adj_per_context: Vec<i64>,
    pub(super) eff_searchsp_per_context: Vec<i64>,
}

/// What the rest of `process_subject` reads of the statistics after a subject's length check: the
/// hit-saving cutoffs (`cutoff_scores`, their minimum, the word finder's `WordFinderCx`) and the
/// effective lengths of each context, the batch's (`context_cutoffs::context_cutoffs`, made once
/// per batch) or, where NCBI computes the statistics of each subject, the subject's own
/// (`context_cutoffs::subject_cutoffs`). One pointer per subject, so that the default path reads
/// the batch's as it did before the statistics of each subject.
pub(super) struct SubjectStatistics<'a> {
    pub(super) cutoff_scores: &'a [i32],
    pub(super) cutoff_score_min: i32,
    pub(super) length_adj_per_context: &'a [i64],
    pub(super) eff_searchsp_per_context: &'a [i64],
    pub(super) word_finder_cx: word_finder::WordFinderCx<'a>,
}

impl<'a> SubjectStatistics<'a> {
    /// The statistics of a subject with its own cutoffs (`own`, from `per_subject_cutoffs`):
    /// `own`'s cutoffs, length adjustments and search spaces, and the batch's word finder values
    /// (`self`'s) with `own`'s `cutoff_scores`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1019-1022
    /// NCBI updates the word finder's cutoffs in place; the batch's stay for the other subjects
    /// here, and the subject's word finder reads its own.
    pub(super) fn of_subject<'s>(&self, own: &'s SubjectCutoffs) -> SubjectStatistics<'s>
    where
        'a: 's,
    {
        SubjectStatistics {
            cutoff_scores: &own.cutoff_scores,
            cutoff_score_min: own.cutoff_score_min,
            length_adj_per_context: &own.length_adj_per_context,
            eff_searchsp_per_context: &own.eff_searchsp_per_context,
            word_finder_cx: word_finder::WordFinderCx {
                cutoff_scores: &own.cutoff_scores,
                ..self.word_finder_cx
            },
        }
    }
}

/// The cutoffs, length adjustments and search spaces of the subject of `subject_len_nucl` bases
/// (`subject_cutoffs`) where NCBI's preliminary search computes the statistics of each subject
/// (`subject_option_spaces`: the preliminary stage's option spaces when
/// `SubjectSetNumbers::per_subject`), or `None` where the subject reads the batch's. Called once
/// per subject that passed the length check, before its frame loop; the frame loop's word finder
/// saves an extension against the subject's cutoffs (`s_BlastAaWordFinder_TwoHit` reads
/// `word_params->cutoffs`, which `BLAST_OneSubjectUpdateParameters` updated). An overflow of the
/// search space is the error.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1444
pub(super) fn per_subject_cutoffs(
    subject_option_spaces: Option<&[i64]>,
    contexts_ref: &[QueryContext],
    numbers: &SubjectSetNumbers,
    subject_len_nucl: i64,
    evalue_threshold: f64,
) -> Result<Option<SubjectCutoffs>, crate::stats::setup::SearchSpaceOverflow> {
    match subject_option_spaces {
        None => Ok(None),
        Some(spaces) => subject_cutoffs(
            contexts_ref,
            numbers,
            spaces,
            subject_len_nucl,
            evalue_threshold,
        )
        .map(Some),
    }
}
