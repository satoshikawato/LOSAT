//! The traceback of a PSSM search with composition-based statistics (`-comp_based_stats` 1, and
//! 2 and 3, which NCBI runs as 1 with a PSSM): NCBI's `Blast_RedoAlignmentCore_MT` for one
//! position-based query (series PS, item PS4b). Every preliminary list is redone on the PSSM
//! that NCBI's `Blast_CompositionBasedStats` scales for the subject's composition
//! (`core::composition_adjustment::compo_scale`), with the gapped alignment of `gapalign.rs` on
//! the scaled rows. The e-values and bit scores read the gapped Lambda of the PSSM file over the
//! scaling factor and the K of the factor-1 PSSM (`PosMatrixInfo::gapped_k`); the report's
//! footer keeps the file's values (`s_RestoreSearch`).
//!
//! The matches are redone in the order of the preliminary hit list, one after the other, at
//! every thread count: NCBI's one-thread loop, whose early termination test and heap insertions
//! read the matches before (DECISIONS 13 of line PS). The PSSM that the last adjustment left
//! stays for the next match, as NCBI's thread-local copy does.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3058-3066
//! ```c
//!     if (positionBased) {
//!         /* Position based searches can only use traditional
//!          * composition based stats */
//!         if ((int) compo_adjust_mode > 1) {
//!             compo_adjust_mode = eCompositionBasedStats;
//!         }
//!         /* A position-based search can only have one query */
//!         ASSERT(queryInfo->num_queries == 1);
//!         ASSERT(queryBlk->length == (Int4)sbp->psi_matrix->pssm->ncols);
//!     }
//! ```

#[allow(unused_imports)]
use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;

use super::evalue::EvalueStats;
use super::gapalign::{pssm_gapped_alignment_with_traceback, PssmGapScratch};
use super::rows::QueryPssm;
use super::traceback0::final_hsp;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{
    fill_results_from_compo_heaps, BlastCompoHeap, BlastpHitList, BlastpHspList,
    COMPO_HEAP_EVALUE_STRETCH,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::core::composition_adjustment::adjust_scores::{
    blast_get_composition_range, read_aa_composition, BlastAminoAcidComposition,
};
use crate::core::composition_adjustment::compo_scale::{
    pssm_composition_based_stats, FrameRedo, PosMatrixInfo, RedoInstance, XDropRedo,
};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, get_nucl_length, get_seq_frame, get_translated_length,
    test_near_identical, windows_from_translated_aligns, BlastCompoAlignment,
    BlastCompoAlignmentContext, BlastCompoSequenceData, EMatrixAdjustRule,
    MINIMUM_LENGTH_NEAR_IDENTICAL, WINDOW_BORDER,
};
use crate::core::link_hsps::{link_in_place, EvenGapLinkParams, LinkContext, LinkHsp};
use crate::stats::KarlinParams;
use crate::utils::matrix::ncbistdaa;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:681-684
// ```c
// #define SCALING_FACTOR 32
// ```
/// `SCALING_FACTOR`: the redo's scores, Lambda and gap costs are scaled by it.
pub(crate) const SCALING_FACTOR: f64 = 32.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2396-2398
// ```c
// #define NEAR_IDENTICAL_BITS_PER_POSITION (1.74)
// ```
const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/compo_mode_condition.h
// and redo_alignment.c:41-58
// ```c
// #define LOCAL_LN2 0.69314718055994530941723212145818
// /** by what factor might initially reported E-value exceed true E-value */
// #define EVALUE_STRETCH 5
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:43
// ```c
// #define KAPPA_BIT_TOL 2.0
// ```
const KAPPA_BIT_TOL: f64 = 2.0;
const LOCAL_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

/// What the redo of every match reads: the PSSM, its matrix information, the query, the
/// subjects, the scoring and hit-saving options.
pub(crate) struct PssmRedo1<'a> {
    /// The search's PSSM (the scores of the file): the matrix of NCBI's copy of the score block
    /// until the first adjustment overwrites it.
    pub(crate) pssm: &'a QueryPssm,
    /// The position-based `Blast_MatrixInfo` (`PosMatrixInfo::new`).
    pub(crate) info: &'a PosMatrixInfo,
    /// The query's context (its letters, its length).
    pub(crate) context: &'a QueryContext,
    pub(crate) subjects: XRedoSubjects<'a>,
    /// The gap costs of the search, before the scaling.
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// The matrix of the search (`sbp->matrix`), for the positives.
    pub(crate) matrix: crate::config::ScoringMatrix,
    /// `kbp_gap_psi[0]`, the PSSM file's gapped block (`sbp->kbp_gap` before the scaling).
    pub(crate) karlin: &'a KarlinParams,
    /// `extendParams->gap_x_dropoff_final`, the final gapped X-drop of the search (raw score).
    pub(crate) x_drop_final: i32,
    /// `options->gap_x_dropoff_final`, the final gapped X-drop in bits (`-xdrop_gap_final` or the
    /// protein default), which the redo converts at its own Lambda (`RedoRun::new`).
    pub(crate) x_drop_final_bits: f64,
    /// `sbp->gbp` (the redo's copies of the score block have it when the search has it:
    /// blast_stat.c:921-923 and blast_kappa.c:2865-2867) and the context's effective search
    /// space.
    pub(crate) stats: EvalueStats<'a>,
    pub(crate) expect_value: f64,
    /// `hitParams->options->hitlist_size` (psiblast's preliminary size).
    pub(crate) hitlist_size: usize,
    /// `psiOptions->inclusion_ethresh`, the heaps' e-value cutoff.
    pub(crate) inclusion_ethresh: f64,
    /// psiblast's `-sum_stats true` (`do_link_hsps`), `None` without sum statistics.
    pub(crate) link: Option<RedoLink>,
}

/// psiblast's `-sum_stats true` in the redo (`do_link_hsps`): the redone alignments' cutoff and
/// the linking of `s_HitlistEvaluateAndPurge` in the place of the e-values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2453-2461
/// ```c
///     if (do_link_hsps) {
///         ASSERT(hitParams->link_hsp_params != NULL);
///         cutoff_s =
///             (int) (hitParams->cutoff_score_min * context->localScalingFactor);
///     } else {
///         /* There is no cutoff score; we consider e-values instead */
///         cutoff_s = 1;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:410-414
/// ```c
///     if (hitParams->do_sum_stats) {
///         status = BLAST_LinkHsps(program_number, hsp_list, queryInfo,
///                                 subject_length, sbp,
///                                 hitParams->link_hsp_params, TRUE);
///     } else {
/// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct RedoLink {
    /// `params->cutoff_s`: the traceback's `cutoff_score_min` times the scaling factor,
    /// truncated.
    pub(crate) cutoff_s: i32,
    /// The traceback's linking parameters (`EvenGapLinkParams::traceback`).
    pub(crate) params: EvenGapLinkParams,
    /// The query's context; the linking reads the Lambda and logK of the redo's scaled block
    /// (`sbp->kbp_gap` of NCBI's copy, `RedoRun`) in the place of this one's.
    pub(crate) context: LinkContext,
}

impl RedoLink {
    /// The parameters of the redo of the search whose traceback cutoff is `cutoff_score_min`
    /// (one context) and whose query context is `context`.
    pub(crate) fn new(cutoff_score_min: i32, context: LinkContext) -> Self {
        Self {
            cutoff_s: (cutoff_score_min as f64 * SCALING_FACTOR) as i32,
            params: EvenGapLinkParams::traceback(),
            context,
        }
    }
}

/// The values of the search that the redo scales (`s_RescaleSearch`) and derives
/// (`s_GetAlignParams`), and its working memory.
pub(crate) struct RedoRun<'a> {
    /// The query (NCBIstdaa, O as X) and its unmasked letters.
    query: &'a [u8],
    query_nomask: &'a [u8],
    /// The query range's letters with U as C (`s_SequenceGetRange`), for the near-identical test.
    query_data: BlastCompoSequenceData,
    /// `query_info->words`.
    query_words: Vec<u64>,
    /// `query_composition->numTrueAminoAcids`.
    query_true_amino_acids: i32,
    /// The scaled gapped block of the redo's copy: Lambda over the scaling factor, K of the
    /// factor-1 PSSM.
    karlin: KarlinParams,
    log_k: f64,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    near_identical_cutoff: f64,
    /// The PSSM of NCBI's thread-local copy: what the last adjustment left.
    matrix: QueryPssm,
    scratch: PssmGapScratch,
    /// RPS-BLAST: whether each adjusted matrix is checked against `INT4_MAX` (the sum of its
    /// rows' best scores), and the first sum that was above it.
    bound_check: bool,
    bound_exceeded: Option<i64>,
}

/// `BlastCompo_Alignment`: an incoming (preliminary) alignment, scaled, or a redone one with its
/// edit script.
#[derive(Debug, Clone)]
struct Align {
    score: i32,
    rule: EMatrixAdjustRule,
    query_start: i32,
    query_end: i32,
    match_start: i32,
    match_end: i32,
    frame: i32,
    edit_script: Vec<GapEditOp>,
}

impl Align {
    /// The incoming alignment of a preliminary HSP (`s_ResultHspToDistinctAlign`): its score
    /// scaled, its subject frame.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:788-797
    /// ```c
    ///         new_align =
    ///             BlastCompo_AlignmentNew((int) (hsp->score * localScalingFactor),
    ///                                     eDontAdjustMatrix,
    ///                                     hsp->query.offset, hsp->query.end, hsp->context,
    ///                                     hsp->subject.offset, hsp->subject.end,
    ///                                     hsp->subject.frame, hsp);
    /// ```
    /// `C`, the redo's instance, keeps each instance's own copy of the code that collects them
    /// (psiblast's code as before the frames' instance; DECISIONS 63 S3).
    fn incoming<C: RedoInstance>(hsp: &BlastpPreliminaryHsp) -> Self {
        Align {
            score: (hsp.raw_score as f64 * SCALING_FACTOR) as i32,
            rule: EMatrixAdjustRule::DontAdjustMatrix,
            query_start: hsp.query_start,
            query_end: hsp.query_end,
            match_start: hsp.subject_start,
            match_end: hsp.subject_end,
            frame: hsp.subject_frame,
            edit_script: Vec::new(),
        }
    }

    /// The alignment as the composition library holds it (`BlastCompo_Alignment` of query 0),
    /// with `context`.
    #[inline(always)]
    fn compo(&self, context: Option<BlastCompoAlignmentContext>) -> BlastCompoAlignment {
        BlastCompoAlignment {
            score: self.score,
            matrix_adjust_rule: self.rule,
            query_index: 0,
            query_start: self.query_start,
            query_end: self.query_end,
            match_start: self.match_start,
            match_end: self.match_end,
            frame: self.frame,
            context,
            next: None,
        }
    }
}

/// A redone HSP of a match's list (`BlastHSP`), before its identities are counted: its e-value and
/// the size of its linked set (`num`, 1 without sum statistics).
struct Redone {
    align: Align,
    e_value: f64,
    num: i32,
}

impl<'a> RedoRun<'a> {
    /// The values of the search that the redo scales (`s_RescaleSearch` on `kbp_gap` and the gap
    /// costs, `s_GetAlignParams`) and the query's data, from the scaled gapped block `karlin`
    /// (Lambda over the scaling factor, the K of the redo's copy) with its `log_k`, the search's
    /// unscaled gap costs, the final gapped X-drop of the search in bits (`-xdrop_gap_final` or the
    /// protein default) and in raw score, and the PSSM the redo starts from (the matrix of NCBI's
    /// thread-local copy, which each adjustment overwrites). psiblast's run starts from the file's
    /// PSSM; RPS-BLAST's from zeros, one run for each profile.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
    /// ```c
    ///     for (i = 0;  i < num_queries;  i++) {
    ///         if (sbp->kbp_gap[i] != NULL) {
    ///             Blast_KarlinBlk * kbp = sbp->kbp_gap[i];
    ///             kbp->Lambda /= scale_factor;
    ///             kbp->logK = log(kbp->K);
    ///         }
    ///     }
    ///
    ///     sp->gap_open = (Int4)BLAST_Nint(sp->gap_open  * scale_factor);
    ///     sp->gap_extend = (Int4)BLAST_Nint(sp->gap_extend * scale_factor);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2437-2450
    /// ```c
    ///     for (index = queryInfo->first_context;
    ///                      index <= queryInfo->last_context; ++index) {
    ///
    ///         if ((queryInfo->contexts[index].is_valid)) {
    ///     		near_identical_cutoff =
    ///        		 (near_identical_cutoff_bits * NCBIMATH_LN2)
    ///         		/ context->sbp->kbp_gap[index]->Lambda;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2375-2383
    /// ```c
    ///     for (i = 0;  i < num_queries;  i++) {
    ///         if (context->sbp->kbp_gap[i] != NULL &&
    ///             context->sbp->kbp_gap[i]->Lambda < min_lambda) {
    ///             min_lambda = context->sbp->kbp_gap[i]->Lambda;
    ///         }
    ///     }
    ///     gapping_params->x_dropoff = (Int4)
    ///         MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    ///             extendParams->gap_x_dropoff_final);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2331 (`s_GetQueryInfo`:
    /// `Blast_ReadAaComposition` and `s_CreateWordArray` of the query)
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        query: &'a [u8],
        query_nomask: &'a [u8],
        karlin: KarlinParams,
        log_k: f64,
        gap_costs: (i32, i32),
        (x_drop_final_bits, x_drop_final): (f64, i32),
        matrix: QueryPssm,
        scratch: PssmGapScratch,
    ) -> Self {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1682-1691
        // ```c
        //     Uint1 *origData = query->data + q_range->begin;
        //     /* Copy the query sequence (necessary for SEG filtering.) */
        //     queryData->length = q_range->end - q_range->begin;
        //     ...
        //         queryData->data[idx] = (origData[idx] != 24) ? origData[idx] : 3;
        // ```
        let query_u_as_c: Vec<u8> = query
            .iter()
            .map(|&residue| {
                if residue == ncbistdaa::U {
                    ncbistdaa::C
                } else {
                    residue
                }
            })
            .collect();
        let lambda = karlin.lambda;
        RedoRun {
            query,
            query_nomask,
            query_data: BlastCompoSequenceData::from_ncbistdaa(&query_u_as_c),
            query_words: build_query_word_hashes(query),
            query_true_amino_acids: read_aa_composition(query).num_true_amino_acids,
            karlin,
            log_k,
            gap_open: blast_nint(gap_costs.0 as f64 * SCALING_FACTOR),
            gap_extend: blast_nint(gap_costs.1 as f64 * SCALING_FACTOR),
            x_dropoff: crate::core::x_drop::redo_x_dropoff(x_drop_final_bits, lambda, x_drop_final),
            near_identical_cutoff: (NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2)
                / lambda,
            matrix,
            scratch,
            bound_check: false,
            bound_exceeded: None,
        }
    }

    /// `MAX(gap_x_dropoff_final * NCBIMATH_LN2 / min_lambda, extendParams->gap_x_dropoff_final)`
    /// before NCBI's `(Int4)` conversion (blast_kappa.c:2381-2383, quoted at `new`): `new`'s
    /// value with the final X-drop `bits` and the redo's Lambda `lambda`.
    pub(crate) fn x_dropoff_of(bits: f64, lambda: f64, x_drop_final: i32) -> f64 {
        (bits * std::f64::consts::LN_2 / lambda).max(x_drop_final as f64)
    }

    /// Checks every matrix the composition adjustment builds: an alignment uses each row once,
    /// so the sum of the rows' best scores bounds its score, and above `INT4_MAX` NCBI's 32-bit
    /// scores overflow (undefined behavior in C).
    pub(crate) fn with_bound_check(mut self) -> Self {
        self.bound_check = true;
        self
    }

    /// The run of a redo whose query information has no valid context: the cutoff of the
    /// near-identical test is 0, NCBI's old test (`old_test_near_identical`). RPS-BLAST's redo
    /// keeps context 0 of the user query's information only, with its validity: rpstblastn's
    /// query whose frame +1 is invalid (`-strand minus`, a frame +1 without letters) gets it.
    /// A protein query has one context, valid where it has hits.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1244
    /// ```c
    ///           one_query_info->first_context = one_query_info->last_context = 0;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2436-2446
    /// ```c
    ///     /* score block is already scaled by context->localScalingFactor */
    ///     double near_identical_cutoff=0;
    ///     Int4 index;
    ///     for (index = queryInfo->first_context;
    ///                      index <= queryInfo->last_context; ++index) {
    ///
    ///         if ((queryInfo->contexts[index].is_valid)) {
    ///     		near_identical_cutoff =
    ///        		 (near_identical_cutoff_bits * NCBIMATH_LN2)
    ///         		/ context->sbp->kbp_gap[index]->Lambda;
    /// ```
    /// (`Blast_GetOneQueryStructs` copies the query's contexts with `is_valid`, and
    /// `OffsetArrayToContextOffsets` sets only the offsets, lengths, frames and query indices.)
    pub(crate) fn without_valid_context(mut self) -> Self {
        self.near_identical_cutoff = 0.0;
        self
    }

    /// The first sum of best row scores above `INT4_MAX` that an adjusted matrix had.
    pub(crate) fn bound_exceeded(&self) -> Option<i64> {
        self.bound_exceeded
    }

    /// The matrix and the alignment scratch, for the next run (RPS-BLAST redoes one profile
    /// after the other on the same memory).
    pub(crate) fn into_buffers(self) -> (QueryPssm, PssmGapScratch) {
        (self.matrix, self.scratch)
    }
}

impl PssmRedo1<'_> {
    /// The final hit lists of the queries from their preliminary lists, which it takes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3093-3110
    /// ```c
    ///     status_code =
    ///         s_RecordInitialSearch(savedParams, sbp, scoringParams,
    ///                               queryInfo->max_length, compo_adjust_mode,
    ///                               positionBased);
    ///     ...
    ///     if (compo_adjust_mode != eNoCompositionBasedStats) {
    ///         if((0 == strcmp(scoringParams->options->matrix, "BLOSUM62_20"))) {
    ///             localScalingFactor = SCALING_FACTOR / 10;
    ///         } else {
    ///             localScalingFactor = SCALING_FACTOR;
    ///         }
    ///     } else {
    ///         localScalingFactor = 1.0;
    ///     }
    ///     s_RescaleSearch(sbp, scoringParams, numContexts, localScalingFactor);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3770-3789
    /// ```c
    ///     for (i = 0; i < actual_num_threads; ++i) {
    ///         if (seqSrc  &&  status_code == 0) {
    ///             s_FillResultsFromCompoHeaps(
    ///                     results_tld[i],
    ///                     redoneMatches_tld[i],
    ///                     hitParams->options->hitlist_size
    ///             );
    /// ```
    pub(crate) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let mut run = self.redo_run();
        let mut heaps: Vec<BlastCompoHeap> = preliminary_hit_lists
            .iter()
            .map(|_| BlastCompoHeap::new(self.hitlist_size, self.inclusion_ethresh))
            .collect();
        for q_idx in 0..preliminary_hit_lists.len() {
            // A list's HSPs as the redo reads them (`s_ResultHspToDistinctAlign`), with its
            // subject and best e-value, in the order of the sorted preliminary hit list.
            let matches: Vec<(u32, f64, Vec<BlastpPreliminaryHsp>)> =
                if let Some(hit_list) = compact_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list.hsps.iter().map(|hsp| hsp.hsp).collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else if let Some(hit_list) = preliminary_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list
                                .hsps
                                .iter()
                                .map(super::traceback0::preliminary_hsp)
                                .collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else {
                    Vec::new()
                };
            for (oid, best_evalue, hsps) in matches {
                if hsps.is_empty() {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3530
                // ```c
                //                 if (BlastCompo_EarlyTermination(
                //                         localMatch->best_evalue,
                //                         redoneMatches,
                //                         numQueries
                //                 )) {
                //                     Blast_HSPListFree(localMatch);
                //                     if (seqSrc) {
                //                         continue;
                //                     }
                // ```
                if early_termination(best_evalue, &heaps) {
                    continue;
                }
                let Some((hsp_list, best_score, best_evalue)) =
                    self.redo_match(&mut run, oid, &hsps)?
                else {
                    continue;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3709-3720
                // ```c
                //                     if (BlastCompo_HeapWouldInsert(
                //                             &redoneMatches[query_index],
                //                             best_evalue,
                //                             best_score,
                //                             localMatch->oid
                //                     )) {
                //                         *pStatusCode =
                //                                 BlastCompo_HeapInsert(
                // ```
                let heap = &mut heaps[q_idx];
                if heap.would_insert(best_evalue, best_score, oid) {
                    heap.insert(hsp_list, best_evalue, best_score, oid);
                }
            }
        }
        Ok(fill_results_from_compo_heaps(self.hitlist_size, &mut heaps))
    }

    /// The scaled values of the search (`s_RescaleSearch` on `kbp_gap` and the gap costs, the
    /// copy's K from `_PSIUpdateLambdaK`, `s_GetAlignParams`) and the query's data.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
    /// ```c
    ///     for (i = 0;  i < num_queries;  i++) {
    ///         if (sbp->kbp_gap[i] != NULL) {
    ///             Blast_KarlinBlk * kbp = sbp->kbp_gap[i];
    ///             kbp->Lambda /= scale_factor;
    ///             kbp->logK = log(kbp->K);
    ///         }
    ///     }
    ///
    ///     sp->gap_open = (Int4)BLAST_Nint(sp->gap_open  * scale_factor);
    ///     sp->gap_extend = (Int4)BLAST_Nint(sp->gap_extend * scale_factor);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2437-2450
    /// ```c
    ///     for (index = queryInfo->first_context;
    ///                      index <= queryInfo->last_context; ++index) {
    ///
    ///         if ((queryInfo->contexts[index].is_valid)) {
    ///     		near_identical_cutoff =
    ///        		 (near_identical_cutoff_bits * NCBIMATH_LN2)
    ///         		/ context->sbp->kbp_gap[index]->Lambda;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2375-2383
    /// ```c
    ///     for (i = 0;  i < num_queries;  i++) {
    ///         if (context->sbp->kbp_gap[i] != NULL &&
    ///             context->sbp->kbp_gap[i]->Lambda < min_lambda) {
    ///             min_lambda = context->sbp->kbp_gap[i]->Lambda;
    ///         }
    ///     }
    ///     gapping_params->x_dropoff = (Int4)
    ///         MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    ///             extendParams->gap_x_dropoff_final);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2331 (`s_GetQueryInfo`:
    /// `Blast_ReadAaComposition` and `s_CreateWordArray` of the query)
    fn redo_run(&self) -> RedoRun<'_> {
        let ctx = self.context;
        let query = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
        let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
        let query_nomask = &query_nomask_full[1..query_nomask_full.len() - 1];
        let lambda = self.karlin.lambda / SCALING_FACTOR;
        RedoRun::new(
            query,
            query_nomask,
            KarlinParams {
                lambda,
                k: self.info.gapped_k,
                ..*self.karlin
            },
            self.info.gapped_k.glibc_ln(),
            (self.gap_open, self.gap_extend),
            (self.x_drop_final_bits, self.x_drop_final),
            self.pssm.clone(),
            PssmGapScratch::new(),
        )
    }

    /// The redo of one match of the PSSM query: the subject read from the database, then
    /// [`RedoMatch::redo_match`].
    fn redo_match(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
    ) -> Result<Option<(BlastpHspList, i32, f64)>> {
        let subject_read = self.subjects.read(oid)?;
        Ok(RedoMatch {
            info: self.info,
            stats: self.stats,
            expect_value: self.expect_value,
            query_length: self.context.aa_len as i32,
            rps_roles: false,
            matrix: self.matrix,
            link: self.link,
        }
        .redo_match(run, oid, hsps, subject_read.residues()))
    }
}

/// The subject of the X-drop redo (NCBI's matching sequence), and the letters its HSPs are
/// counted on: a protein sequence (psiblast's subject, rpsblast's query: [`ProteinSubject`]) or
/// a translated sequence's six frames (rpstblastn's query: [`FrameSubject`]), whose letters are
/// NCBI's translation buffer (`one_query->sequence`, `one_query->length` bytes: the frames in
/// the context order, one sentinel between two frames with letters).
pub(crate) trait RedoSubject<'a>: Copy {
    /// The instance of the steps the redos share.
    type Instance: RedoInstance;
    /// The letters an HSP of `frame` aligns to, as the redo's subject `subject` holds them.
    fn frame_letters(self, subject: &'a [u8], frame: i32) -> &'a [u8];
    /// The redo of one match of this subject: `RedoMatch::redo_match` (a protein subject) or
    /// `RedoMatch::redo_windows` (frames).
    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)>;
}

/// A protein subject ([`RedoSubject`]): the redo's subject is its letters.
#[derive(Clone, Copy)]
pub(crate) struct ProteinSubject;

impl<'a> RedoSubject<'a> for ProteinSubject {
    type Instance = XDropRedo;

    #[inline(always)]
    fn frame_letters(self, subject: &'a [u8], _frame: i32) -> &'a [u8] {
        subject
    }

    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        matching.redo_match(run, oid, hsps, subject)
    }
}

/// A translated subject ([`RedoSubject`]): the letters of each frame, without sentinels, by
/// context (`BLAST_FrameToContext`).
#[derive(Clone, Copy)]
pub(crate) struct FrameSubject<'a> {
    pub(crate) frames: &'a [&'a [u8]],
}

impl<'a> RedoSubject<'a> for FrameSubject<'a> {
    type Instance = FrameRedo;

    fn frame_letters(self, _subject: &'a [u8], frame: i32) -> &'a [u8] {
        self.frames[crate::core::blast_util::frame_to_context(frame)]
    }

    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        matching.redo_windows(run, oid, hsps, subject, self)
    }
}

/// A window of the subject (`s_WindowInfo`'s `subject_range`): `[begin, end)` of the frame
/// `frame`'s letters. `end` may lie before `begin` (`end - begin` is NCBI's length, then below 0:
/// `rpstblastn -strand`, whose translation buffer is shorter than `GET_NUCL_LENGTH` assumes). A
/// protein subject's window (the whole sequence) is never read.
#[derive(Clone, Copy, Default)]
struct SubjectWindow {
    begin: i32,
    end: i32,
    frame: i32,
}

/// What the redo of one match reads besides its run: the matrix information, the Gumbel block
/// and the search space (`EvalueStats`), `-evalue`, the query's length and the roles of the
/// lengths in the e-values (RPS-BLAST's query is the profile and its subject the user's query:
/// `Blast_HSPListGetEvalues` swaps them).
pub(crate) struct RedoMatch<'a> {
    pub(crate) info: &'a PosMatrixInfo,
    pub(crate) stats: EvalueStats<'a>,
    pub(crate) expect_value: f64,
    pub(crate) query_length: i32,
    pub(crate) rps_roles: bool,
    /// The matrix of the search (`sbp->matrix`), for the positives.
    pub(crate) matrix: crate::config::ScoringMatrix,
    /// psiblast's `-sum_stats true` (`RedoLink`); `None` without sum statistics and for
    /// RPS-BLAST, which refuses it.
    pub(crate) link: Option<RedoLink>,
}

impl RedoMatch<'_> {
    /// The redo of one match (the HSPs of the query and one subject): `Blast_RedoOneMatch`, then
    /// `s_HSPListFromDistinctAlignments`, `s_HitlistReapContained`, `s_HitlistEvaluateAndPurge`,
    /// and for a list whose best e-value passes `-evalue`, `s_HSPListNormalizeScores` and
    /// `s_ComputeNumIdentities`. `None` when no HSP passes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3672-3707
    /// ```c
    ///                 if (hsp_list->hspcnt > 1) {
    ///                     s_HitlistReapContained(hsp_list->hsp_array,
    ///                             &hsp_list->hspcnt);
    ///                 }
    ///                 *pStatusCode =
    ///                         s_HitlistEvaluateAndPurge(&best_score, &best_evalue,
    ///                                 hsp_list,
    ///                                 ...
    ///                 if (best_evalue <= hitParams->options->expect_value) {
    ///                     /* The best alignment is significant */
    ///                     s_HSPListNormalizeScores(hsp_list, kbp->Lambda, kbp->logK,
    ///                             localScalingFactor);
    ///                     s_ComputeNumIdentities(
    /// ```
    pub(crate) fn redo_match(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        let alignments = self.redo_one_match(run, hsps, subject);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-354
        // ```c
        //     for (align = *alignments;  NULL != align;  align = align->next) {
        //         ...
        //         status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
        //     ...
        //         Blast_HSPListSortByScore(hsp_list);
        // ```
        // The sort is glibc's `qsort`, a stable merge sort (see `blastp::hsp::sort_hsplist_by_score`).
        let mut redone: Vec<Redone> = alignments
            .into_iter()
            .map(|align| Redone {
                align,
                e_value: 0.0,
                num: 1,
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained::<XDropRedo>(&mut redone);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:410-446
        // ```c
        //     *pbestEvalue = DBL_MAX;
        //     *pbestScore  = 0;
        //     ...
        //         status =
        //             Blast_HSPListGetEvalues(program_number, queryInfo,
        //             		                s_GetSubjectLength(subject_length, program_number),
        //                                     hsp_list, TRUE, FALSE, sbp,
        //                                     0.0, /* use a non-zero gap decay
        //                                             only when linking HSPs */
        //                                     1.0); /* Use scaling factor equal to
        //                                              1, because both scores and
        //                                              Lambda are scaled, so they
        //                                              will cancel each other. */
        //     ...
        //         Blast_HSPListReapByEvalue(hsp_list, hitParams->options);
        //         if (hsp_list->hspcnt > 0) {
        //             *pbestEvalue = hsp_list->best_evalue;
        //             *pbestScore  = hsp_list->hsp_array[0]->score;
        //         }
        // ```
        // psiblast's e-values are not adjusted by a composition p-value (that is blastp's and
        // blastx's, with `-comp_based_stats 2u`/`3u`).
        let query_length = self.query_length;
        let subject_length = subject.len() as i32;
        let best_evalue = match &self.link {
            Some(link) => link_redone(&mut redone, link, run, subject_length)?,
            None => {
                let mut best_evalue = i32::MAX as f64;
                for hsp in &mut redone {
                    // (RPS-BLAST's roles are swapped, `Blast_HSPListGetEvalues`: its `isRPS`.)
                    let (m, n) = if self.rps_roles {
                        (subject_length, query_length)
                    } else {
                        (query_length, subject_length)
                    };
                    hsp.e_value = self.stats.e_value(hsp.align.score, &run.karlin, m, n);
                    best_evalue = best_evalue.min(hsp.e_value);
                }
                best_evalue
            }
        };
        redone.retain(|hsp| hsp.e_value <= self.expect_value);
        let first = redone.first()?;
        let best_score = first.align.score;
        // NCBI's `best_evalue <= hitParams->options->expect_value` holds: the list keeps an HSP,
        // whose e-value is at most `-evalue`, and `best_evalue` is the smallest of all.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-116
        // ```c
        //         hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);
        //         /* Compute the bit score using the newly computed scaled score. */
        //         hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
        // ```
        let template = &hsps[0];
        let mut hsp_list = BlastpHspList {
            oid,
            query_index: template.q_idx,
            hsps: Vec::with_capacity(redone.len()),
            best_evalue,
        };
        for Redone {
            align,
            e_value,
            num,
        } in redone
        {
            let score = blast_nint(align.score as f64 / SCALING_FACTOR);
            let traced = BlastpPreliminaryHsp {
                query_start: align.query_start,
                query_end: align.query_end,
                subject_start: align.match_start,
                subject_end: align.match_end,
                gapped_query_start: align.query_start,
                gapped_subject_start: align.match_start,
                raw_score: score,
                query_context: template.q_idx as i32,
                query_frame: template.query_frame,
                subject_frame: 0,
                query_length: template.query_length,
                q_idx: template.q_idx,
                s_idx: template.s_idx,
            };
            let mut hsp = final_hsp(
                &traced,
                align.edit_script,
                run.query_nomask,
                subject,
                self.matrix,
            );
            hsp.e_value = e_value;
            hsp.num = num;
            hsp.bit_score = (score as f64 * run.karlin.lambda * SCALING_FACTOR - run.log_k)
                / std::f64::consts::LN_2;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:332-342
            // ```c
            //         switch (align->matrix_adjust_rule) {
            //         case eDontAdjustMatrix:
            //             new_hsp->comp_adjustment_method = eNoCompositionBasedStats;
            //             break;
            //         case eCompoScaleOldMatrix:
            //             new_hsp->comp_adjustment_method = eCompositionBasedStats;
            //             break;
            // ```
            hsp.comp_adjustment_method = match align.rule {
                EMatrixAdjustRule::DontAdjustMatrix => 0,
                EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                _ => 2,
            };
            hsp_list.hsps.push(hsp);
        }
        Some((hsp_list, best_score, best_evalue))
    }

    /// `Blast_RedoOneMatch` for the one window of a protein query and subject (the whole query,
    /// the whole subject, `s_WindowsFromProteinAligns`): the distinct redone alignments, the
    /// list's head first.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1174-1286
    /// ```c
    ///         for (in_align = window->align, hsp_index = 0;
    ///                 in_align != NULL;
    ///                 in_align = in_align->next, hsp_index++) {
    ///             if (hsp_index == 0 || subject_maybe_biased) {
    ///                 nearIdenticalStatus = s_preliminaryTestNearIdentical(query_info,
    ///                         window,
    ///                         in_align,
    ///                         params->near_identical_cutoff);
    ///             }
    ///             if (hsp_index == 0 ||
    ///                     (subject_maybe_biased &&
    ///                             (nearIdenticalStatus != oldNearIdenticalStatus))) {
    ///                 ...
    ///                 status = callbacks->get_range(
    ///                         matchingSeq,
    ///                         &window->subject_range,
    ///                         &subject,
    ///                         &query_info[query_index].seq,
    ///                         &window->query_range,
    ///                         &query,
    ///                         query_info[query_index].words,
    ///                         window->align,
    ///                         nearIdenticalStatus,
    ///                         compo_adjust_mode,
    ///                         FALSE,
    ///                         &subject_maybe_biased
    ///                 );
    ///             ...
    ///             if ( !s_IsContained(in_align, alignments[query_index], Lambda) ) {
    ///                 ...
    ///                 int adjust_search_failed = 0;
    ///                 if (compo_adjust_mode != eNoCompositionBasedStats &&
    ///                         (subject_is_translated || hsp_index == 0
    ///                                 || (nearIdenticalStatus != oldNearIdenticalStatus))) {
    ///                     Blast_AminoAcidComposition subject_composition;
    ///                     s_GetComposition(&subject_composition,
    ///                             alphsize, &subject,
    ///                             &window->subject_range,
    ///                             in_align, FALSE, subject_is_translated);
    ///                     adjust_search_failed =
    ///                             Blast_AdjustScores(matrix, query_composition,
    ///                                     ...
    ///                     num_adjustments++;
    ///                 }
    ///
    ///                 if ( !adjust_search_failed ) {
    ///                     newAlign = callbacks->redo_one_alignment(
    ///                             in_align,
    ///                             matrix_adjust_rule,
    ///                             ...
    ///                     if (newAlign && newAlign->score >= params->cutoff_s) {
    ///                         s_WithDistinctEnds(&newAlign, &alignments[query_index],
    ///                                 callbacks->free_align_traceback,
    ///                                 num_adjustments == 1);
    ///                     } else {
    ///                         ...
    ///             oldNearIdenticalStatus = nearIdenticalStatus;
    ///         } /* end for all alignments in this window */
    /// ```
    /// `params->cutoff_s` is 1 without sum statistics, else `RedoLink::cutoff_s`
    /// (blast_kappa.c:2453-2461).
    fn redo_one_match(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Vec<Align> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:788-797
        // ```c
        //         new_align =
        //             BlastCompo_AlignmentNew((int) (hsp->score * localScalingFactor),
        //                                     eDontAdjustMatrix,
        //                                     hsp->query.offset, hsp->query.end, hsp->context,
        //                                     hsp->subject.offset, hsp->subject.end,
        //                                     hsp->subject.frame, hsp);
        // ```
        let incoming: Vec<Align> = hsps
            .iter()
            .map(|hsp| Align {
                score: (hsp.raw_score as f64 * SCALING_FACTOR) as i32,
                rule: EMatrixAdjustRule::DontAdjustMatrix,
                query_start: hsp.query_start,
                query_end: hsp.query_end,
                match_start: hsp.subject_start,
                match_end: hsp.subject_end,
                frame: hsp.subject_frame,
                edit_script: Vec::new(),
            })
            .collect();
        let lambda = run.karlin.lambda;
        let query_length = run.query.len() as i32;
        let mut alignments: Vec<Align> = Vec::new();
        let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
        let mut subject_data: Vec<u8> = Vec::new();
        let mut near_identical = true;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        for (hsp_index, in_align) in incoming.iter().enumerate() {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical = preliminary_test_near_identical(
                    query_length,
                    in_align,
                    run.near_identical_cutoff,
                );
            }
            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                subject_data = self.subject_range::<XDropRedo>(
                    run,
                    subject,
                    &incoming[0],
                    near_identical,
                    &mut subject_maybe_biased,
                );
            }
            if !is_contained::<XDropRedo>(in_align, &alignments, lambda) {
                let mut adjust_search_failed = false;
                if hsp_index == 0 || near_identical != old_near_identical {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1436-1442
                    // ```c
                    //     if (query_composition->numTrueAminoAcids == 0 ||
                    //         subject_composition->numTrueAminoAcids == 0) {
                    //         /* Either the query or subject contains only ambiguity
                    //            characters, most likely because the entire subject has been
                    //            SEGed.  Compositional adjustment is meaningless. */
                    //         return 1;
                    //     }
                    // ```
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1484-1529
                    // ```c
                    //     if (matrixInfo->positionBased ||
                    //         composition_adjust_mode == eCompositionBasedStats) {
                    //         /* Use old-style composition-based statistics unconditionally. */
                    //         *matrix_adjust_rule =  eCompoScaleOldMatrix;
                    //     ...
                    //     *matrix_adjust_rule = eCompoScaleOldMatrix;
                    //     return Blast_CompositionBasedStats(matrix, ratioToPassBack, matrixInfo,
                    //                                        query_composition->prob,
                    //                                        subject_composition->prob,
                    //                                        calc_lambda,
                    //                                        (compositionTestIndex > 0));
                    // ```
                    let subject_composition = read_aa_composition(&subject_data);
                    if run.query_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
                        pssm_composition_based_stats::<XDropRedo>(
                            &mut run.matrix.rows,
                            self.info,
                            &subject_composition.prob,
                        );
                        if run.bound_check && run.bound_exceeded.is_none() {
                            let bound: i64 = run
                                .matrix
                                .rows
                                .iter()
                                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                                .sum();
                            if bound > i64::from(i32::MAX) {
                                run.bound_exceeded = Some(bound);
                            }
                        }
                    }
                    num_adjustments += 1;
                }
                if !adjust_search_failed {
                    let new_align = self.redo_one_alignment(
                        run,
                        &hsps[hsp_index],
                        matrix_adjust_rule,
                        &subject_data,
                    );
                    if new_align.score >= self.link.map_or(1, |link| link.cutoff_s) {
                        with_distinct_ends::<XDropRedo>(
                            new_align,
                            &mut alignments,
                            num_adjustments == 1,
                        );
                    }
                }
            }
            old_near_identical = near_identical;
        }
        alignments
    }

    /// The redo of one match of a translated subject (rpstblastn): `redo_match`'s steps (quoted
    /// there) over the windows of the frames (`redo_frame_windows`), with the frames compared in
    /// the reap, the subject's length of `s_GetSubjectLength` in the e-values and the frames kept
    /// in the HSPs. A copy (DECISIONS 63 S3): psiblast's and rpsblast's `redo_match` stays dev's
    /// code, which one generic form shared with the frames did not keep instruction for
    /// instruction; its instance of a protein subject is `redo_match`'s for the equivalence test
    /// (`the_window_copy_is_the_x_drop_redo`).
    pub(crate) fn redo_windows<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
        kind: S,
    ) -> Option<(BlastpHspList, i32, f64)> {
        // A protein subject is one window, the whole subject (`s_WindowsFromProteinAligns`).
        let alignments = if S::Instance::SUBJECT_IS_TRANSLATED {
            self.redo_frame_windows::<S>(run, hsps, subject)
        } else {
            self.redo_window::<S>(run, hsps, subject, SubjectWindow::default(), None)
        };
        let mut redone: Vec<Redone> = alignments
            .into_iter()
            .map(|align| Redone {
                align,
                e_value: 0.0,
                num: 1,
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained::<S::Instance>(&mut redone);
        }
        let query_length = self.query_length;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:364-368
        // ```c
        // Int4 s_GetSubjectLength(Int4 total_subj_length, EBlastProgramType program_number)
        // {
        // 	return ((program_number == eBlastTypeRpsTblastn) ?
        // 				(GET_NUCL_LENGTH(total_subj_length) - 1 ) /3 : total_subj_length);
        // }
        // ```
        let subject_length = if S::Instance::SUBJECT_IS_TRANSLATED {
            (get_nucl_length(subject.len() as i32) - 1) / 3
        } else {
            subject.len() as i32
        };
        let best_evalue = match &self.link {
            Some(link) => {
                // The RPS redo has no sum statistics (`rps_search::redo`): only psiblast links.
                assert!(
                    !S::Instance::SUBJECT_IS_TRANSLATED,
                    "the redo of frames has no sum statistics"
                );
                link_redone(&mut redone, link, run, subject_length)?
            }
            None => {
                let mut best_evalue = i32::MAX as f64;
                for hsp in &mut redone {
                    let (m, n) = if self.rps_roles {
                        (subject_length, query_length)
                    } else {
                        (query_length, subject_length)
                    };
                    hsp.e_value = self.stats.e_value(hsp.align.score, &run.karlin, m, n);
                    best_evalue = best_evalue.min(hsp.e_value);
                }
                best_evalue
            }
        };
        redone.retain(|hsp| hsp.e_value <= self.expect_value);
        let first = redone.first()?;
        let best_score = first.align.score;
        let template = &hsps[0];
        let mut hsp_list = BlastpHspList {
            oid,
            query_index: template.q_idx,
            hsps: Vec::with_capacity(redone.len()),
            best_evalue,
        };
        for Redone {
            align,
            e_value,
            num,
        } in redone
        {
            let score = blast_nint(align.score as f64 / SCALING_FACTOR);
            let traced = BlastpPreliminaryHsp {
                query_start: align.query_start,
                query_end: align.query_end,
                subject_start: align.match_start,
                subject_end: align.match_end,
                gapped_query_start: align.query_start,
                gapped_subject_start: align.match_start,
                raw_score: score,
                query_context: template.q_idx as i32,
                query_frame: template.query_frame,
                // (`Blast_HSPInit(.., (Int2) align->frame, ..)`, blast_kappa.c:326-330.)
                subject_frame: if S::Instance::SUBJECT_IS_TRANSLATED {
                    align.frame
                } else {
                    0
                },
                query_length: template.query_length,
                q_idx: template.q_idx,
                s_idx: template.s_idx,
            };
            let mut hsp = final_hsp(
                &traced,
                align.edit_script,
                run.query_nomask,
                kind.frame_letters(subject, align.frame),
                self.matrix,
            );
            hsp.e_value = e_value;
            hsp.num = num;
            hsp.bit_score = (score as f64 * run.karlin.lambda * SCALING_FACTOR - run.log_k)
                / std::f64::consts::LN_2;
            hsp.comp_adjustment_method = match align.rule {
                EMatrixAdjustRule::DontAdjustMatrix => 0,
                EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                _ => 2,
            };
            hsp_list.hsps.push(hsp);
        }
        Some((hsp_list, best_score, best_evalue))
    }

    /// `Blast_RedoOneMatch` of a translated subject: the distinct redone alignments, the list's
    /// head first. (A protein subject is one window, the whole subject with the alignments in
    /// their order, `s_WindowsFromProteinAligns`: [`Self::redo_window`].) The windows are the
    /// alignments' frames around them (`s_WindowsFromTranslatedAligns` with a position-based
    /// matrix, `frame_windows`), each with its alignments in `s_DistinctAlignmentsSort`'s order,
    /// and all windows add to the one list (`alignments[query_index]`). NCBI keeps
    /// `matrix_adjust_rule` from one window to the next; a translated subject's matrix is
    /// adjusted for each alignment that is redone (`redo_window`), so the rule a window starts
    /// with is never read.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1121-1141
    /// ```c
    ///     EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;
    ///     ...
    ///     status =
    ///             s_WindowsFromAligns(incoming_aligns, query_info, hspcnt, numQueries,
    ///                     kWindowBorder, matchingSeq->length, &windows,
    ///                     &nWindows, query_is_translated, subject_is_translated, params->positionBased);
    ///     if (status != 0) {
    ///         goto function_level_cleanup;
    ///     }
    ///     /* for all windows */
    ///     for (window_index = 0;  window_index < nWindows;  window_index++) {
    /// ```
    fn redo_frame_windows<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Vec<Align> {
        let mut alignments: Vec<Align> = Vec::new();
        for (window, order) in frame_windows(hsps, run.query.len() as i32, subject.len() as i32) {
            let window_hsps: Vec<BlastpPreliminaryHsp> = order.iter().map(|&k| hsps[k]).collect();
            alignments =
                self.redo_window::<S>(run, &window_hsps, subject, window, Some(alignments));
        }
        alignments
    }

    /// The redo of one window of `Blast_RedoOneMatch` (`redo_one_match`'s steps, quoted there;
    /// `window`, whose alignments are `hsps` in their order): the list of the windows before
    /// (`before`, none for the first) with the window's distinct redone alignments at its head. A
    /// translated subject's range is the window of its frame (`frame_range`), and the matrix is
    /// adjusted for each alignment that is not contained, to the composition around it in the
    /// window (`s_GetComposition`, `window_composition`).
    fn redo_window<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
        window: SubjectWindow,
        before: Option<Vec<Align>>,
    ) -> Vec<Align> {
        // The incoming alignments (`Align::incoming`) in the window's order.
        let incoming: Vec<Align> = hsps.iter().map(Align::incoming::<S::Instance>).collect();
        let lambda = run.karlin.lambda;
        let query_length = run.query.len() as i32;
        let mut alignments: Vec<Align> = before.unwrap_or_default();
        let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
        let mut subject_data: Vec<u8> = Vec::new();
        let mut near_identical = true;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        for (hsp_index, in_align) in incoming.iter().enumerate() {
            if hsp_index == 0 || subject_maybe_biased {
                // Only rpstblastn's redo (the translated instance) can have the cutoff 0
                // (`RedoRun::without_valid_context`); psiblast's and rpsblast's keep the one test.
                near_identical =
                    if S::Instance::SUBJECT_IS_TRANSLATED && run.near_identical_cutoff <= 0.0 {
                        old_test_near_identical(query_length, incoming.len(), in_align)
                    } else {
                        preliminary_test_near_identical(
                            query_length,
                            in_align,
                            run.near_identical_cutoff,
                        )
                    };
            }
            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                subject_data = if S::Instance::SUBJECT_IS_TRANSLATED {
                    frame_range(
                        run,
                        subject,
                        window,
                        &incoming[0],
                        near_identical,
                        &mut subject_maybe_biased,
                    )
                } else {
                    self.subject_range::<S::Instance>(
                        run,
                        subject,
                        &incoming[0],
                        near_identical,
                        &mut subject_maybe_biased,
                    )
                };
            }
            if !is_contained::<S::Instance>(in_align, &alignments, lambda) {
                let mut adjust_search_failed = false;
                if S::Instance::SUBJECT_IS_TRANSLATED
                    || hsp_index == 0
                    || near_identical != old_near_identical
                {
                    let subject_composition = if S::Instance::SUBJECT_IS_TRANSLATED {
                        window_composition(&subject_data, window, in_align)
                    } else {
                        read_aa_composition(&subject_data)
                    };
                    if run.query_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
                        pssm_composition_based_stats::<S::Instance>(
                            &mut run.matrix.rows,
                            self.info,
                            &subject_composition.prob,
                        );
                        if run.bound_check && run.bound_exceeded.is_none() {
                            let bound: i64 = run
                                .matrix
                                .rows
                                .iter()
                                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                                .sum();
                            if bound > i64::from(i32::MAX) {
                                run.bound_exceeded = Some(bound);
                            }
                        }
                    }
                    num_adjustments += 1;
                }
                if !adjust_search_failed {
                    let new_align = if S::Instance::SUBJECT_IS_TRANSLATED {
                        redo_window_alignment(
                            run,
                            &hsps[hsp_index],
                            matrix_adjust_rule,
                            &subject_data,
                            window,
                        )
                    } else {
                        self.redo_one_alignment(
                            run,
                            &hsps[hsp_index],
                            matrix_adjust_rule,
                            &subject_data,
                        )
                    };
                    if new_align.score >= self.link.map_or(1, |link| link.cutoff_s) {
                        with_distinct_ends::<S::Instance>(
                            new_align,
                            &mut alignments,
                            num_adjustments == 1,
                        );
                    }
                }
            }
            old_near_identical = near_identical;
        }
        alignments
    }

    /// `s_SequenceGetProteinRange` for the whole subject: its letters, masked by SEG ("S 10 1.8
    /// 2.1") unless the window's first alignment is near identical; `subject_maybe_biased`
    /// becomes whether SEG masked anything.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1638
    /// ```c
    ///     if ( !(KAPPA_BLASTP_NO_SEG_SEQUENCE) ) {
    ///         if (compo_adjust_mode
    ///             && (!subject_maybe_biased || *subject_maybe_biased)) {
    ///
    ///             if ( (!shouldTestIdentical)
    ///                  || (shouldTestIdentical
    ///                      && (!s_TestNearIdentical(seqData, 0, queryData,
    ///                                               q_range->begin, query_words,
    ///                                               align)))) {
    ///
    ///                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
    ///                                              subject_maybe_biased);
    ///             }
    ///         }
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1001-1004
    /// ```c
    /// /** NCBIstdaa encoding for 'X' character */
    /// #define BLASTP_MASK_RESIDUE 21
    /// /** Default instructions and mask residue for SEG filtering */
    /// #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
    /// ```
    fn subject_range<C: RedoInstance>(
        &self,
        run: &RedoRun<'_>,
        subject: &[u8],
        window_align: &Align,
        near_identical: bool,
        subject_maybe_biased: &mut bool,
    ) -> Vec<u8> {
        let mut data = subject.to_vec();
        if *subject_maybe_biased {
            let test_identical = near_identical && {
                let align = BlastCompoAlignment {
                    score: window_align.score,
                    matrix_adjust_rule: window_align.rule,
                    query_index: 0,
                    query_start: window_align.query_start,
                    query_end: window_align.query_end,
                    match_start: window_align.match_start,
                    match_end: window_align.match_end,
                    frame: window_align.frame,
                    context: None,
                    next: None,
                };
                test_near_identical(
                    &BlastCompoSequenceData::from_ncbistdaa(subject),
                    0,
                    &run.query_data,
                    0,
                    Some(&run.query_words),
                    &align,
                )
            };
            if !test_identical {
                let intervals =
                    crate::utils::seg::SegMasker::new(10, 1.8, 2.1).mask_sequence(&data);
                *subject_maybe_biased = !intervals.is_empty();
                for interval in &intervals {
                    data[interval.start..interval.end].fill(ncbistdaa::X);
                }
            }
        }
        data
    }

    /// `s_RedoOneAlignment`: the gapped alignment with traceback of the preliminary HSP `hsp`
    /// from its gapped start, on the PSSM of the last adjustment, with the scaled gap costs and
    /// X-drop (`s_NewAlignmentFromGapAlign` of the whole query and subject).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1924-1945
    /// ```c
    ///     q_start = hsp->query.gapped_start - query_range->begin;
    ///     s_start = hsp->subject.gapped_start - subject_range->begin;
    ///
    ///     gapAlign->gap_x_dropoff = gapping_params->x_dropoff;
    ///     ...
    ///     status =
    ///         BLAST_GappedAlignmentWithTraceback(context->prog_number,
    ///                                            query_data->data,
    ///                                            subject_data->data, gapAlign,
    ///                                            context->scoringParams,
    ///                                            q_start, s_start,
    ///                                            query_data->length,
    ///                                            subject_data->length,
    ///                                            &fence_hit);
    ///     if (status == 0) {
    ///         return s_NewAlignmentFromGapAlign(gapAlign, &gapAlign->edit_script,
    ///                                           query_range, subject_range,
    ///                                           matrix_adjust_rule);
    /// ```
    fn redo_one_alignment(
        &self,
        run: &mut RedoRun<'_>,
        hsp: &BlastpPreliminaryHsp,
        rule: EMatrixAdjustRule,
        subject: &[u8],
    ) -> Align {
        let traced = pssm_gapped_alignment_with_traceback(
            &run.matrix,
            subject,
            subject.len(),
            hsp.gapped_query_start as usize,
            hsp.gapped_subject_start as usize,
            run.gap_open,
            run.gap_extend,
            run.x_dropoff,
            &mut run.scratch,
        );
        Align {
            score: traced.score,
            rule,
            query_start: traced.query_start,
            query_end: traced.query_stop,
            match_start: traced.subject_start,
            match_end: traced.subject_stop,
            frame: 0,
            edit_script: traced.edit_script,
        }
    }
}

/// `s_DoSegSequenceData` of the redo's subject letters `data` ("S 10 1.8 2.1", the masked
/// letters X): `subject_maybe_biased` becomes whether SEG masked anything.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1435-1446
/// ```c
///     status = BlastFilteringOptionsFromString(program_name,
///                                              BLASTP_MASK_INSTRUCTIONS,
///                                              &filter_options, NULL);
///     if (status == 0) {
///         status = BlastSetUp_Filter(program_name, seqData->data,
///                                    seqData->length, 0, filter_options,
///                                    &mask_seqloc, NULL);
///         filter_options = SBlastFilterOptionsFree(filter_options);
///     }
///     if (is_seq_biased) {
///         *is_seq_biased = (mask_seqloc != NULL);
///     }
/// ```
#[inline(always)]
fn seg_subject(data: &mut [u8], subject_maybe_biased: &mut bool) {
    let intervals = crate::utils::seg::SegMasker::new(10, 1.8, 2.1).mask_sequence(data);
    *subject_maybe_biased = !intervals.is_empty();
    for interval in &intervals {
        data[interval.start..interval.end].fill(ncbistdaa::X);
    }
}

/// `s_WindowsFromTranslatedAligns` of the incoming alignments of `hsps` with a position-based
/// matrix (`windows_from_translated_aligns`, `is_pos_based`) in a translated subject of
/// `subject_length` bytes (the translation buffer), for a query of `query_length` letters:
/// each window and the indices in `hsps` of its alignments, in the order the redo reads them.
fn frame_windows(
    hsps: &[BlastpPreliminaryHsp],
    query_length: i32,
    subject_length: i32,
) -> Vec<(SubjectWindow, Vec<usize>)> {
    let mut list: Option<Box<BlastCompoAlignment>> = None;
    for (k, hsp) in hsps.iter().enumerate().rev() {
        let mut align = Box::new(
            Align::incoming::<FrameRedo>(hsp)
                .compo(Some(BlastCompoAlignmentContext::PreliminaryHspIndex(k))),
        );
        align.next = list;
        list = Some(align);
    }
    let windows =
        windows_from_translated_aligns(&list, &[query_length], subject_length, WINDOW_BORDER, true)
            .expect("the redo's alignments are of query 0");
    windows
        .into_iter()
        .map(|window| {
            let mut order = Vec::with_capacity(window.hspcnt as usize);
            let mut current = window.align.as_deref();
            while let Some(align) = current {
                if let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(k)) = align.context {
                    order.push(k);
                }
                current = align.next.as_deref();
            }
            let range = window.subject_range;
            (
                SubjectWindow {
                    begin: range.begin,
                    end: range.end,
                    frame: range.context,
                },
                order,
            )
        })
        .collect()
}

/// `s_SequenceGetProteinRange` of a translated subject (`self->index < 0`, `align->frame !=
/// 0`) for the window `window` whose first alignment is `window_align`: NCBI's buffer of the
/// translation buffer's length plus 2 bytes (`calloc`), whose second byte on holds the frame's
/// letters as NCBI finds them in `subject` (`GET_TRANSLATED_LENGTH` letters of `GET_NUCL_LENGTH`
/// of the buffer's length, after the earlier frames' lengths plus a sentinel each), masked by SEG
/// unless the window's first alignment is near identical (`s_TestNearIdentical` on that buffer);
/// then the byte before the window is zeroed (the fit: the window's letters start at byte
/// `1 + window.begin`).
///
/// Every read is inside `subject` and the buffer, also where NCBI's frame arithmetic does not
/// match the buffer (`-strand plus|minus`, whose buffer holds three frames): with `nl =
/// GET_NUCL_LENGTH(l)` the three lengths `GET_TRANSLATED_LENGTH(nl, f)` of `f` = 0, 1, 2 add up
/// to `nl - 2`, so the copy of frame index 5 ends at `2 * (nl - 2) + 5 = 2 * ((l - 5) / 2) + 5 <=
/// l`; a frame has letters only where `nl >= 3`, `l >= 7`. The window's alignments end inside
/// their frame, so before `l`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1597-1643
/// ```c
///     seqData->buffer  = calloc((self->length + 2), sizeof(Uint1));
///     ...
///     seqData->data    = seqData->buffer + 1;
///     seqData->length  = self->length;
///
///     origData = (self->index >= 0) ? local_data->seq_arg.seq->sequence
///                             : seq->sequence;
///     if((self->index < 0) && (align->frame != 0)) {
///     	int i=0, offsets =0;
///     	int f = GET_SEQ_FRAME(align->frame);
///     	int nucl_length = GET_NUCL_LENGTH(self->length);
///     	seqData->length = GET_TRANSLATED_LENGTH(nucl_length, f);
///     	for(; i < f; i++) {
///     		offsets = GET_TRANSLATED_LENGTH(nucl_length, i) +1;
///     		origData += offsets;
///     	}
///     }
///     /* Copy the sequence data */
///     for (idx = 0;  idx < seqData->length;  idx++) {
///         seqData->data[idx] = origData[idx];
///     }
///
///     if ( !(KAPPA_BLASTP_NO_SEG_SEQUENCE) ) {
///         if (compo_adjust_mode
///             && (!subject_maybe_biased || *subject_maybe_biased)) {
///
///             if ( (!shouldTestIdentical)
///                  || (shouldTestIdentical
///                      && (!s_TestNearIdentical(seqData, 0, queryData,
///                                               q_range->begin, query_words,
///                                               align)))) {
///
///                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
///                                              subject_maybe_biased);
///             }
///         }
///     }
///     /* Fit the data to the range. */
///     seqData ->data    = &seqData->data[range->begin - 1];
///     *seqData->data++  = '\0';
///     seqData ->length  = range->end - range->begin;
/// ```
fn frame_range(
    run: &RedoRun<'_>,
    subject: &[u8],
    window: SubjectWindow,
    window_align: &Align,
    near_identical: bool,
    subject_maybe_biased: &mut bool,
) -> Vec<u8> {
    let length = subject.len();
    let mut buffer = vec![0u8; length + 2];
    let f = get_seq_frame(window_align.frame);
    let nucl_length = get_nucl_length(length as i32);
    let frame_length = get_translated_length(nucl_length, f) as usize;
    if frame_length > 0 {
        let offset: usize = (0..f)
            .map(|i| get_translated_length(nucl_length, i) as usize + 1)
            .sum();
        buffer[1..1 + frame_length].copy_from_slice(&subject[offset..offset + frame_length]);
    }
    if *subject_maybe_biased {
        // `s_TestNearIdentical` reads the buffer (not only the frame's letters) at the window's
        // first alignment.
        let data = BlastCompoSequenceData {
            buffer,
            data_offset: 1,
            length: length as i32 + 1,
        };
        let test_identical = near_identical
            && test_near_identical(
                &data,
                0,
                &run.query_data,
                0,
                Some(&run.query_words),
                &window_align.compo(None),
            );
        buffer = data.buffer;
        if !test_identical {
            seg_subject(&mut buffer[1..1 + frame_length], subject_maybe_biased);
        }
    }
    buffer[window.begin as usize] = 0;
    buffer
}

/// `s_GetComposition` of a translated subject: the composition of the letters of `window`
/// (in `data`, `frame_range`'s buffer) around the alignment `align` (`Blast_GetCompositionRange`
/// of its range in the window's `end - begin` letters). Where the window holds fewer letters
/// than the alignment reaches (`-strand plus|minus`), NCBI's range ends at the alignment's end
/// (the scan to the right stops at once): the letters read here reach that end, which stays
/// inside the buffer (`frame_range`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:921-936
/// ```c
///     data = seq->data;
///     length = range->end - range->begin;
///     if (query_is_translated || subject_is_translated) {
///         int start;
///         int end;
///         start = ((query_is_translated) ?
///                  align->queryStart : align->matchStart) - range->begin;
///         end   = ((query_is_translated) ?
///                  align->queryEnd   : align->matchEnd  ) - range->begin;
///         Blast_GetCompositionRange(&left, &right, data, length, start, end);
///     } else {
///         /* Use the whole subject to compute the composition */
///         left = 0;
///         right = length;
///     }
///     Blast_ReadAaComposition(composition, alphsize, &data[left], right-left);
/// ```
fn window_composition(
    data: &[u8],
    window: SubjectWindow,
    align: &Align,
) -> BlastAminoAcidComposition {
    let length = (window.end - window.begin).max(0) as usize;
    let start = (align.match_start - window.begin) as usize;
    let finish = (align.match_end - window.begin) as usize;
    let letters = &data[1 + window.begin as usize..][..length.max(finish)];
    let (left, right) = blast_get_composition_range(letters, start, finish)
        .expect("an alignment's range lies in the letters read");
    read_aa_composition(&letters[left..right])
}

/// `s_RedoOneAlignment` in a window of a translated subject: the gapped alignment of the
/// preliminary HSP `hsp` from its gapped start on the window's letters (in `data`,
/// `frame_range`'s buffer; `end - begin` of them, none where the window ends before it begins:
/// then only the left extension runs, NCBI's `s_start < s_length` failing), shifted back into the
/// frame (`s_NewAlignmentFromGapAlign`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1925-1926
/// ```c
///     q_start = hsp->query.gapped_start - query_range->begin;
///     s_start = hsp->subject.gapped_start - subject_range->begin;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1765
/// ```c
///     queryStart = gap_align->query_start + query_range->begin;
///     queryEnd   = gap_align->query_stop + query_range->begin;
///     queryIndex = query_range->context;
///     matchStart = gap_align->subject_start + subject_range->begin;
///     matchEnd   = gap_align->subject_stop  + subject_range->begin;
///     frame      = subject_range->context;
/// ```
fn redo_window_alignment(
    run: &mut RedoRun<'_>,
    hsp: &BlastpPreliminaryHsp,
    rule: EMatrixAdjustRule,
    data: &[u8],
    window: SubjectWindow,
) -> Align {
    let traced = pssm_gapped_alignment_with_traceback(
        &run.matrix,
        &data[1 + window.begin as usize..],
        (window.end - window.begin).max(0) as usize,
        hsp.gapped_query_start as usize,
        (hsp.gapped_subject_start - window.begin) as usize,
        run.gap_open,
        run.gap_extend,
        run.x_dropoff,
        &mut run.scratch,
    );
    Align {
        score: traced.score,
        rule,
        query_start: traced.query_start,
        query_end: traced.query_stop,
        match_start: traced.subject_start + window.begin,
        match_end: traced.subject_stop + window.begin,
        frame: window.frame,
        edit_script: traced.edit_script,
    }
}

/// `BLAST_LinkHsps` on a match's redone list (psiblast's `-sum_stats true`): every HSP takes its
/// linked set's e-value and size, in the list the linking builds, sorted by score; the best
/// e-value of the list, `None` for an empty list. The scores are the redo's scaled ones and the
/// linking reads the redo's scaled Lambda and logK (`RedoLink::context`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1801-1811
/// ```c
///     /* Sort the HSP array by score */
///     Blast_HSPListSortByScore(hsp_list);
///
///     /* Find and fill the best e-value */
///     hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
///     for (index = 1; index < hsp_list->hspcnt; ++index) {
///         if (hsp_list->hsp_array[index]->evalue < hsp_list->best_evalue)
///             hsp_list->best_evalue = hsp_list->hsp_array[index]->evalue;
///     }
/// ```
fn link_redone(
    redone: &mut Vec<Redone>,
    link: &RedoLink,
    run: &RedoRun<'_>,
    subject_length: i32,
) -> Option<f64> {
    if redone.is_empty() {
        return None;
    }
    let context = LinkContext {
        lambda: run.karlin.lambda,
        log_k: run.log_k,
        ..link.context
    };
    link_in_place(
        redone,
        |hsp| LinkHsp {
            context: 0,
            query_offset: hsp.align.query_start,
            query_end: hsp.align.query_end,
            subject_offset: hsp.align.match_start,
            subject_end: hsp.align.match_end,
            score: hsp.align.score,
        },
        |_| context,
        subject_length,
        &link.params,
        |hsp, linked| {
            hsp.e_value = linked.evalue;
            hsp.num = linked.num;
        },
    );
    // `Blast_HSPListSortByScore`: glibc's `qsort`, a stable merge sort, when the list is not
    // sorted (a stable sort of a sorted list keeps it).
    redone.sort_by(|a, b| score_compare(&a.align, &b.align));
    let best_evalue = redone[1..].iter().fold(redone[0].e_value, |best, hsp| {
        if hsp.e_value < best {
            hsp.e_value
        } else {
            best
        }
    });
    Some(best_evalue)
}

/// `BLAST_Nint`.
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

/// `BlastCompo_EarlyTermination`: every query's heap holds its hit list size and its worst
/// e-value is at most its cutoff, and `evalue` is above five times every cutoff.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1561-1576
/// ```c
///     int i;
///     for (i = 0;  i < numQueries;  i++) {
///         if (BlastCompo_HeapFilledToCutoff(&significantMatches[i])) {
///             double ecutoff = significantMatches[i].ecutoff;
///             /* Only matches with evalue <= ethresh will be saved */
///             if (evalue <= EVALUE_STRETCH * ecutoff) {
///                 /* The evalue if this match is sufficiently small
///                  * that we want to redo it to try to obtain an
///                  * alignment with evalue smaller than ecutoff. */
///                 return FALSE;
///             }
///         } else {
///             return FALSE;
///         }
///     }
///     return TRUE;
/// ```
fn early_termination(evalue: f64, heaps: &[BlastCompoHeap]) -> bool {
    heaps.iter().all(|heap| {
        heap.n >= heap.heap_threshold
            && heap.worst_evalue <= heap.ecutoff
            && evalue > COMPO_HEAP_EVALUE_STRETCH * heap.ecutoff
    })
}

/// `s_preliminaryTestNearIdentical` with a cutoff above 0 (the redo's).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1063-1076
/// ```c
///   if (cutoff > 0) {
///       queryIndex = align->queryIndex;
///       queryLength = query_info[queryIndex].seq.length;
///       if ((align->matchEnd - align->matchStart +1) <
///           (MIN(queryLength,  MINIMUM_LENGTH_NEAR_IDENTICAL)))
///           return(FALSE);
///
///       align_len = MIN(align->queryEnd - align->queryStart,
///                       align->matchEnd - align->matchStart);
///
///       if ((double)align->score / (double)align_len < cutoff) {
///           return (FALSE);
///       }
///   }
/// ```
fn preliminary_test_near_identical(query_length: i32, align: &Align, cutoff: f64) -> bool {
    if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL) {
        return false;
    }
    let align_len = (align.query_end - align.query_start).min(align.match_end - align.match_start);
    // A float division as NCBI's (0 gives an infinity or NaN, neither below the cutoff).
    if (align.score as f64) / (align_len as f64) < cutoff {
        return false;
    }
    true
}

/// `s_preliminaryTestNearIdentical` with the cutoff 0 (`RedoRun::without_valid_context`): the
/// window `hspcnt` holds one alignment, without gaps, at least as long as the profile or 50.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1078-1094
/// ```c
///   else {
///       /* Otherwise fallback to the old near identical test. This is left for
///          compatibility with the code in c toolkit's tools/kappa.c */
///
///       if ((window->hspcnt > 1) || (window->hspcnt < 1))
///           return(FALSE);
///
///       queryIndex = align->queryIndex;
///       queryLength = query_info[queryIndex].seq.length;
///       if ((align->queryEnd - align->queryStart) !=
///           (align->matchEnd - align->matchStart))
///           return(FALSE);
///
///       if ((align->matchEnd - align->matchStart +1) <
///           (MIN(queryLength, MINIMUM_LENGTH_NEAR_IDENTICAL)))
///           return(FALSE);
///   }
/// ```
fn old_test_near_identical(query_length: i32, hspcnt: usize, align: &Align) -> bool {
    hspcnt == 1
        && align.query_end - align.query_start == align.match_end - align.match_start
        && align.match_end - align.match_start + 1
            >= query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
}

/// `KAPPA_CONTAINED_IN_HSP` (redo_alignment.c:72-73) and `CONTAINED_IN_HSP` (blast_kappa.c).
fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    (a <= c && b >= c) && (d <= f && e >= f)
}

/// `s_IsContained`: whether an alignment of the list holds both ends of `in_align` and scores
/// at least `in_align`'s score plus 2 bits.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:966-995
/// ```c
///     double score        = in_align->score;
///     double scoreThresh = score + KAPPA_BIT_TOL * LOCAL_LN2/lambda;
///
///     for (align = alignments;  align != NULL;  align = align->next ) {
///         /* for all elements of alignments */
///         if (KAPPA_SIGN(in_align->frame) == KAPPA_SIGN(align->frame)) {
///             /* hsp1 and hsp2 are in the same query/subject frame */
///             if (KAPPA_CONTAINED_IN_HSP
///                 (align->queryStart, align->queryEnd, query_offset,
///                  align->matchStart, align->matchEnd, subject_offset) &&
///                 KAPPA_CONTAINED_IN_HSP
///                 (align->queryStart, align->queryEnd, query_end,
///                  align->matchStart, align->matchEnd, subject_end) &&
///                 scoreThresh <= align->score) {
///                 return 1;
///             }
///         }
///     }
///     return 0;
/// ```
///
/// `C`, the redo's instance, keeps each instance's own copy (psiblast's code as before the frames'
/// instance; DECISIONS 63 S3).
fn is_contained<C: RedoInstance>(in_align: &Align, alignments: &[Align], lambda: f64) -> bool {
    let score_thresh = in_align.score as f64 + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;
    alignments.iter().any(|align| {
        in_align.frame.signum() == align.frame.signum()
            && contained_in_hsp(
                align.query_start,
                align.query_end,
                in_align.query_start,
                align.match_start,
                align.match_end,
                in_align.match_start,
            )
            && contained_in_hsp(
                align.query_start,
                align.query_end,
                in_align.query_end,
                align.match_start,
                align.match_end,
                in_align.match_end,
            )
            && score_thresh <= align.score as f64
    })
}

/// `s_WithDistinctEnds`: `new_align` joins the list (at its head, without the alignments that
/// share an end with it) unless an alignment that shares an end (with the same adjustment) or a
/// similar end (otherwise) scores at least as much.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:400-449
/// ```c
///     for (align = oldAlignments;  align != NULL;  align = align->next) {
///         ...
///         if (align->frame == newAlign->frame &&
///             ((is_same_adjustment && s_IsSameEndPoint(newAlign, align)) ||
///              (!is_same_adjustment && s_IsSimilarEndPoint(newAlign, align)))) {
///             /* At least one of the endpoints of newAlign matches an endpoint
///                of align. */
///             if (newAlign->score <= align->score) {
///                 /* newAlign cannot be added to the list. */
///                 include_new_align = 0;
///                 break;
///             }
///         }
///     }
///     if (include_new_align) {
///         ...
///         while (align != NULL) {
///             ...
///             if (align->frame == newAlign->frame &&
///                 ((align->queryStart == newAlign->queryStart &&
///                   align->matchStart == newAlign->matchStart)
///                  || (align->queryEnd == newAlign->queryEnd &&
///                      align->matchEnd == newAlign->matchEnd))) {
///                 /* The alignment shares an end with newAlign; */
///                 /* delete it. */
///                 BlastCompo_AlignmentsFree(&align, free_align_context);
///             } else { /* The alignment does not share an end with newAlign; */
///                 /* add it to the output list. */
///                 *tail =  align;
///                 tail  = &align->next;
///             }
///             align = align_next;
///         } /* end while align != NULL */
///         *p_oldAlignments = newAlign;
/// ```
///
/// `C`, the redo's instance, keeps each instance's own copy (psiblast's code as before the frames'
/// instance; DECISIONS 63 S3).
fn with_distinct_ends<C: RedoInstance>(
    new_align: Align,
    alignments: &mut Vec<Align>,
    is_same_adjustment: bool,
) {
    let excluded = alignments.iter().any(|align| {
        align.frame == new_align.frame
            && if is_same_adjustment {
                is_same_end_point(&new_align, align)
            } else {
                is_similar_end_point(&new_align, align)
            }
            && new_align.score <= align.score
    });
    if excluded {
        return;
    }
    alignments
        .retain(|align| !(align.frame == new_align.frame && is_same_end_point(&new_align, align)));
    alignments.insert(0, new_align);
}

/// `s_IsSameEndPoint`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:305-314
/// ```c
///     return ((align->queryStart == newAlign->queryStart &&
///              align->matchStart == newAlign->matchStart)
///             || (align->queryEnd == newAlign->queryEnd &&
///                 align->matchEnd == newAlign->matchEnd));
/// ```
fn is_same_end_point(new_align: &Align, align: &Align) -> bool {
    (align.query_start == new_align.query_start && align.match_start == new_align.match_start)
        || (align.query_end == new_align.query_end && align.match_end == new_align.match_end)
}

/// `s_IsSimilarEndPoint`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:321-348
/// ```c
///     Boolean start_contained =
///         KAPPA_CONTAINED_IN_HSP(align->queryStart, align->queryEnd,
///                                newAlign->queryStart,
///                                align->matchStart, align->matchEnd,
///                                newAlign->matchStart);
///     Boolean end_contained =
///         KAPPA_CONTAINED_IN_HSP(align->queryStart, align->queryEnd,
///                                newAlign->queryEnd,
///                                align->matchStart, align->matchEnd,
///                                newAlign->matchEnd);
///     Boolean result =  (start_contained &&
///                        (newAlign->queryStart - newAlign->matchStart ==
///                         align->queryStart - align->matchStart)) ||
///                       (end_contained &&
///                        (newAlign->queryEnd - newAlign->matchEnd ==
///                         align->queryEnd - align->matchEnd));
/// ```
fn is_similar_end_point(new_align: &Align, align: &Align) -> bool {
    let start_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_start,
        align.match_start,
        align.match_end,
        new_align.match_start,
    );
    let end_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_end,
        align.match_start,
        align.match_end,
        new_align.match_end,
    );
    (start_contained
        && new_align.query_start - new_align.match_start == align.query_start - align.match_start)
        || (end_contained
            && new_align.query_end - new_align.match_end == align.query_end - align.match_end)
}

/// `ScoreCompareHSPs`: the score, highest first, then the subject start, the subject end (last
/// first), the query start and the query end (last first).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
/// ```c
///    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
///        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
///        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
///        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
///        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
///    }
/// ```
fn score_compare(a: &Align, b: &Align) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.match_start.cmp(&b.match_start))
        .then_with(|| b.match_end.cmp(&a.match_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

/// `s_HitlistReapContained`: an HSP whose ends lie in an earlier HSP of the same frames and
/// whose score is at most that HSP's is deleted.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:233-262
/// ```c
///     for (iread = 1;  iread < *hspcnt;  iread++) {
///         ...
///         hsp1 = hsp_array[iread];
///         for (ireadBack = 0;  ireadBack < iread && hsp1 != NULL;  ireadBack++) {
///             ...
///             hsp2 = hsp_array[ireadBack];
///
///             if( hsp2 == NULL ) {  /* hsp2 was deleted in a prior iteration. */
///                 continue;
///             }
///             if (hsp2->query.frame == hsp1->query.frame &&
///                 hsp2->subject.frame == hsp1->subject.frame) {
///                 /* hsp1 and hsp2 are in the same query/subject frame. */
///                 if (CONTAINED_IN_HSP
///                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
///                      hsp2->subject.offset, hsp2->subject.end,
///                      hsp1->subject.offset) &&
///                     CONTAINED_IN_HSP
///                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
///                      hsp2->subject.offset, hsp2->subject.end,
///                      hsp1->subject.end)    &&
///                     hsp1->score <= hsp2->score) {
///                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
///                 }
/// ```
/// Every HSP of a PSSM search has the query frame 0, and the subject frame 0 unless the subject
/// is translated (`RedoInstance::SUBJECT_IS_TRANSLATED`).
fn reap_contained<C: RedoInstance>(redone: &mut Vec<Redone>) {
    let mut keep = vec![true; redone.len()];
    for iread in 1..redone.len() {
        let hsp1 = &redone[iread].align;
        for ireadback in 0..iread {
            if !keep[ireadback] {
                continue;
            }
            let hsp2 = &redone[ireadback].align;
            if (!C::SUBJECT_IS_TRANSLATED || hsp2.frame == hsp1.frame)
                && contained_in_hsp(
                    hsp2.query_start,
                    hsp2.query_end,
                    hsp1.query_start,
                    hsp2.match_start,
                    hsp2.match_end,
                    hsp1.match_start,
                )
                && contained_in_hsp(
                    hsp2.query_start,
                    hsp2.query_end,
                    hsp1.query_end,
                    hsp2.match_start,
                    hsp2.match_end,
                    hsp1.match_end,
                )
                && hsp1.score <= hsp2.score
            {
                keep[iread] = false;
                break;
            }
        }
    }
    let mut index = 0;
    redone.retain(|_| {
        index += 1;
        keep[index - 1]
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039 (quoted in the
// module): `-use_sw_tback` (`eSmithWatermanTbck`) redoes each match with
// `Blast_RedoOneMatchSmithWaterman` in the place of `Blast_RedoOneMatch` (series PR, item
// PS11c-3: RPS-BLAST's redo).
pub(crate) mod smith_waterman;

#[cfg(test)]
mod frame_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn align(score: i32, q: (i32, i32), s: (i32, i32)) -> Align {
        Align {
            score,
            rule: EMatrixAdjustRule::CompoScaleOldMatrix,
            query_start: q.0,
            query_end: q.1,
            match_start: s.0,
            match_end: s.1,
            frame: 0,
            edit_script: Vec::new(),
        }
    }

    // `s_WithDistinctEnds`: a new alignment that shares an end with a better one is dropped; one
    // that beats it replaces every alignment that shares an end and goes to the head.
    #[test]
    fn distinct_ends_keep_the_best_of_a_shared_end() {
        let mut list = vec![align(100, (0, 50), (0, 50)), align(80, (60, 90), (60, 90))];
        with_distinct_ends::<XDropRedo>(align(90, (0, 40), (0, 40)), &mut list, true);
        assert_eq!(list.len(), 2);
        with_distinct_ends::<XDropRedo>(align(120, (10, 90), (10, 90)), &mut list, true);
        assert_eq!(
            list.iter().map(|a| a.score).collect::<Vec<_>>(),
            vec![120, 100]
        );
        // With different adjustments a similar end (contained, same diagonal) counts.
        with_distinct_ends::<XDropRedo>(align(110, (5, 45), (5, 45)), &mut list, false);
        assert_eq!(list.len(), 2);
    }

    // `s_HitlistReapContained` keeps the first of two nested HSPs unless the later scores more.
    #[test]
    fn contained_hsps_are_reaped() {
        let mut redone = vec![
            Redone {
                align: align(100, (0, 50), (0, 50)),
                e_value: 0.0,
                num: 1,
            },
            Redone {
                align: align(90, (10, 40), (10, 40)),
                e_value: 0.0,
                num: 1,
            },
            Redone {
                align: align(80, (60, 70), (60, 70)),
                e_value: 0.0,
                num: 1,
            },
        ];
        reap_contained::<XDropRedo>(&mut redone);
        assert_eq!(
            redone.iter().map(|r| r.align.score).collect::<Vec<_>>(),
            vec![100, 80]
        );
    }

    // `s_preliminaryTestNearIdentical`: the subject span must reach min(query length, 50) and
    // the score per aligned column the cutoff.
    #[test]
    fn the_preliminary_near_identical_test_is_ncbis() {
        let a = align(1000, (0, 60), (0, 60));
        assert!(preliminary_test_near_identical(100, &a, 16.0));
        assert!(!preliminary_test_near_identical(100, &a, 17.0));
        let short = align(1000, (0, 20), (0, 20));
        assert!(!preliminary_test_near_identical(100, &short, 1.0));
        assert!(preliminary_test_near_identical(20, &short, 1.0));
    }
}
