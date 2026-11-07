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

use anyhow::Result;

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
use crate::algorithm::tblastx::constants::X_DROP_GAPPED_FINAL;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::core::composition_adjustment::adjust_scores::read_aa_composition;
use crate::core::composition_adjustment::compo_scale::{
    pssm_composition_based_stats, PosMatrixInfo,
};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, test_near_identical, BlastCompoAlignment, BlastCompoSequenceData,
    EMatrixAdjustRule, MINIMUM_LENGTH_NEAR_IDENTICAL,
};
use crate::stats::{blast_spouge_stoe, BlastGumbelBlk, KarlinParams};
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
    /// `kbp_gap_psi[0]`, the PSSM file's gapped block (`sbp->kbp_gap` before the scaling).
    pub(crate) karlin: &'a KarlinParams,
    /// `extendParams->gap_x_dropoff_final`, the final gapped X-drop of the search (raw score).
    pub(crate) x_drop_final: i32,
    /// `sbp->gbp`.
    pub(crate) gumbel: &'a BlastGumbelBlk,
    pub(crate) expect_value: f64,
    /// `hitParams->options->hitlist_size` (psiblast's preliminary size).
    pub(crate) hitlist_size: usize,
    /// `psiOptions->inclusion_ethresh`, the heaps' e-value cutoff.
    pub(crate) inclusion_ethresh: f64,
}

/// The values of the search that the redo scales (`s_RescaleSearch`) and derives
/// (`s_GetAlignParams`), and its working memory.
struct RedoRun<'a> {
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

/// A redone HSP of a match's list (`BlastHSP`), before its identities are counted.
struct Redone {
    align: Align,
    e_value: f64,
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
        let lambda = self.karlin.lambda / SCALING_FACTOR;
        RedoRun {
            query,
            query_nomask,
            query_data: BlastCompoSequenceData::from_ncbistdaa(&query_u_as_c),
            query_words: build_query_word_hashes(query),
            query_true_amino_acids: read_aa_composition(query).num_true_amino_acids,
            karlin: KarlinParams {
                lambda,
                k: self.info.gapped_k,
                ..*self.karlin
            },
            log_k: self.info.gapped_k.ln(),
            gap_open: blast_nint(self.gap_open as f64 * SCALING_FACTOR),
            gap_extend: blast_nint(self.gap_extend as f64 * SCALING_FACTOR),
            x_dropoff: (X_DROP_GAPPED_FINAL as f64 * std::f64::consts::LN_2 / lambda)
                .max(self.x_drop_final as f64) as i32,
            near_identical_cutoff: (NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2)
                / lambda,
            matrix: self.pssm.clone(),
            scratch: PssmGapScratch::new(),
        }
    }

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
    fn redo_match(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
    ) -> Result<Option<(BlastpHspList, i32, f64)>> {
        let subject_read = self.subjects.read(oid)?;
        let subject = subject_read.residues();
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
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained(&mut redone);
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
        let query_length = self.context.aa_len as i32;
        let subject_length = subject.len() as i32;
        let mut best_evalue = i32::MAX as f64;
        for hsp in &mut redone {
            hsp.e_value = blast_spouge_stoe(
                hsp.align.score,
                &run.karlin,
                self.gumbel,
                query_length,
                subject_length,
            );
            best_evalue = best_evalue.min(hsp.e_value);
        }
        redone.retain(|hsp| hsp.e_value <= self.expect_value);
        let Some(first) = redone.first() else {
            return Ok(None);
        };
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
        for Redone { align, e_value } in redone {
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
            let mut hsp = final_hsp(&traced, align.edit_script, run.query_nomask, subject);
            hsp.e_value = e_value;
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
        Ok(Some((hsp_list, best_score, best_evalue)))
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
    /// `params->cutoff_s` is 1 (no sum statistics, blast_kappa.c:2453-2461).
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
                subject_data = self.subject_range(
                    run,
                    subject,
                    &incoming[0],
                    near_identical,
                    &mut subject_maybe_biased,
                );
            }
            if !is_contained(in_align, &alignments, lambda) {
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
                        pssm_composition_based_stats(
                            &mut run.matrix.rows,
                            self.info,
                            &subject_composition.prob,
                        );
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
                    if new_align.score >= 1 {
                        with_distinct_ends(new_align, &mut alignments, num_adjustments == 1);
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
    fn subject_range(
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
fn is_contained(in_align: &Align, alignments: &[Align], lambda: f64) -> bool {
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
fn with_distinct_ends(new_align: Align, alignments: &mut Vec<Align>, is_same_adjustment: bool) {
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
/// Every HSP of a PSSM search has the query frame and subject frame 0.
fn reap_contained(redone: &mut Vec<Redone>) {
    let mut keep = vec![true; redone.len()];
    for iread in 1..redone.len() {
        let hsp1 = &redone[iread].align;
        for ireadback in 0..iread {
            if !keep[ireadback] {
                continue;
            }
            let hsp2 = &redone[ireadback].align;
            if contained_in_hsp(
                hsp2.query_start,
                hsp2.query_end,
                hsp1.query_start,
                hsp2.match_start,
                hsp2.match_end,
                hsp1.match_start,
            ) && contained_in_hsp(
                hsp2.query_start,
                hsp2.query_end,
                hsp1.query_end,
                hsp2.match_start,
                hsp2.match_end,
                hsp1.match_end,
            ) && hsp1.score <= hsp2.score
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
        with_distinct_ends(align(90, (0, 40), (0, 40)), &mut list, true);
        assert_eq!(list.len(), 2);
        with_distinct_ends(align(120, (10, 90), (10, 90)), &mut list, true);
        assert_eq!(
            list.iter().map(|a| a.score).collect::<Vec<_>>(),
            vec![120, 100]
        );
        // With different adjustments a similar end (contained, same diagonal) counts.
        with_distinct_ends(align(110, (5, 45), (5, 45)), &mut list, false);
        assert_eq!(list.len(), 2);
    }

    // `s_HitlistReapContained` keeps the first of two nested HSPs unless the later scores more.
    #[test]
    fn contained_hsps_are_reaped() {
        let mut redone = vec![
            Redone {
                align: align(100, (0, 50), (0, 50)),
                e_value: 0.0,
            },
            Redone {
                align: align(90, (10, 40), (10, 40)),
                e_value: 0.0,
            },
            Redone {
                align: align(80, (60, 70), (60, 70)),
                e_value: 0.0,
            },
        ];
        reap_contained(&mut redone);
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
