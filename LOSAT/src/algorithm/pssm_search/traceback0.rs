//! The traceback of a PSSM search without composition-based statistics (`-comp_based_stats 0`):
//! NCBI's `Blast_TracebackFromHSPList` on each preliminary HSP list with `gap_align->positionBased`,
//! a PSSM copy of BLASTP's ordinary traceback (`blastp/blast_engine/traceback.rs`), whose kernels
//! read the 28×28 matrix. The starts and the alignments read the PSSM rows (`gapalign.rs`); the
//! identities and positives read the search's matrix (`sbp->matrix`, blast_hits.c:770-772); the
//! e-values and bit scores read `kbp_gap_psi` (`sbp->kbp_gap`, `Blast_QueryIsPssm`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1500
//! ```c
//!     if (ext_params->options->compositionBasedStats > 0 ||
//!         ext_params->options->eTbackExt == eSmithWatermanTbck) {
//!         status =
//!             Blast_RedoAlignmentCore_MT(program_number, ...);
//!     } else {
//!         BlastSeqSrcGetSeqArg seq_arg;
//! ```
//!
//! As for BLASTP, a list's traceback reads nothing another list writes and the hit lists' order
//! is total, so every list is inserted into the one hit list of the query in the lists' order, at
//! every thread count.

use anyhow::Result;

use super::evalue::EvalueStats;
use super::gapalign::{
    pssm_gapped_alignment_with_traceback, pssm_get_offsets_for_gapped_alignment, PssmGapScratch,
};
use super::prelim::{purge_hsps_with_common_endpoints, score_compare};
use super::rows::PssmRows;
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::gapalign::{adjust_subject_range, BlastpPreliminaryHsp};
use crate::algorithm::blastp::hsp::{
    reap_hsplist_by_evalue, BlastpHitList, BlastpHsp, BlastpHspList,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::stats::KarlinParams;
use crate::utils::matrix::standard_score;

/// What the traceback of every list reads: the PSSM, the query's one context, the subjects, the
/// scoring and the hit-saving options.
pub(crate) struct PssmTraceback0<'a, R: PssmRows> {
    pub(crate) rows: &'a R,
    /// The query's context (its letters for the identities, its length).
    pub(crate) context: &'a QueryContext,
    /// `query_blk->length + 1`, the query bound of NCBI's containment tree.
    pub(crate) total_query_span: i32,
    pub(crate) subjects: XRedoSubjects<'a>,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// The matrix of the search (`sbp->matrix`), whose positive scores are the positives.
    pub(crate) matrix: ScoringMatrix,
    /// `ext_params->gap_x_dropoff_final` (`kbp_gap_psi`'s Lambda).
    pub(crate) x_drop_final: i32,
    /// `kbp_gap_psi[0]`.
    pub(crate) karlin: &'a KarlinParams,
    /// `sbp->gbp` and the context's effective search space, for the e-values.
    pub(crate) stats: EvalueStats<'a>,
    pub(crate) expect_value: f64,
    pub(crate) hitlist_size: usize,
}

/// A preliminary list of the query and one subject, full or compact.
enum PrelimList {
    Full(BlastpHspList),
    Compact(BlastpHspList<XPrelimHsp>),
}

/// An HSP that the traceback aligned, with its edit script (`BlastHSP::gap_info`).
struct Traced {
    hsp: BlastpPreliminaryHsp,
    edit_script: Vec<GapEditOp>,
}

impl<R: PssmRows> PssmTraceback0<'_, R> {
    /// The final hit list of the query from its preliminary lists, which it takes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1690-1699
    /// ```c
    ///                 /* Free HSP list if all HSPs have been deleted. */
    ///                 batch->hsplist_array[hsplist_itr] = NULL;
    ///                 if (hsp_list->hspcnt == 0) {
    ///                     hsp_list = Blast_HSPListFree(hsp_list);
    ///                 }
    ///                 else {
    ///                     Blast_HSPResultsInsertHSPList(thread_data->tld[tid]->results, hsp_list,
    ///                                   hit_params->options->hitlist_size);
    ///                 }
    /// ```
    pub(crate) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let mut work: Vec<PrelimList> = Vec::new();
        for (full, compact) in preliminary_hit_lists
            .iter_mut()
            .zip(compact_hit_lists.iter_mut())
        {
            if let Some(hit_list) = compact.take() {
                let count = hit_list.hsplist_count;
                work.extend(
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(PrelimList::Compact),
                );
            } else if let Some(hit_list) = full.take() {
                let count = hit_list.hsplist_count;
                work.extend(
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(PrelimList::Full),
                );
            }
        }
        let traced = self.trace_lists(work, pool)?;
        let mut hit_lists: Vec<Option<BlastpHitList>> = std::iter::repeat_with(|| None)
            .take(preliminary_hit_lists.len())
            .collect();
        for hsp_list in traced.into_iter().flatten() {
            hit_lists[hsp_list.query_index as usize]
                .get_or_insert_with(|| BlastpHitList::new(self.hitlist_size))
                .update(hsp_list);
        }
        Ok(hit_lists)
    }

    /// Each list's traceback, in the lists' order.
    fn trace_lists(
        &self,
        work: Vec<PrelimList>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<Option<BlastpHspList>>> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() && pool.threads() > 1 {
            use rayon::prelude::*;
            return pool.install(|| {
                work.into_par_iter()
                    .map_init(PssmGapScratch::new, |scratch, list| {
                        self.trace(list, scratch)
                    })
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = PssmGapScratch::new();
        work.into_iter()
            .map(|list| self.trace(list, &mut scratch))
            .collect()
    }

    fn trace(
        &self,
        list: PrelimList,
        scratch: &mut PssmGapScratch,
    ) -> Result<Option<BlastpHspList>> {
        match list {
            PrelimList::Full(list) => self.trace_list(
                list.oid,
                list.query_index,
                list.hsps.iter().map(preliminary_hsp),
                scratch,
            ),
            PrelimList::Compact(list) => self.trace_list(
                list.oid,
                list.query_index,
                list.hsps.iter().map(|hsp| hsp.hsp),
                scratch,
            ),
        }
    }

    /// `Blast_TracebackFromHSPList` and `s_HSPListPostTracebackUpdate` on the list of one subject;
    /// `None` when no HSP is left.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-373
    /// ```c
    ///    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
    ///                                  0, (subject_length > 0 ? subject_length :
    ///                                  subject_blk->length / CODON_LENGTH) + 1);
    /// ```
    fn trace_list(
        &self,
        oid: u32,
        query_index: u32,
        hsps: impl Iterator<Item = BlastpPreliminaryHsp>,
        scratch: &mut PssmGapScratch,
    ) -> Result<Option<BlastpHspList>> {
        let ctx = self.context;
        let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
        let query_nomask = &query_nomask_full[1..query_nomask_full.len() - 1];
        let subject_read = self.subjects.read(oid)?;
        let subject = subject_read.residues();
        let mut hsp_list = BlastpHspList {
            oid,
            query_index,
            hsps: trace_hsps(
                &TraceHsps {
                    rows: self.rows,
                    query_nomask,
                    query_length: ctx.aa_len,
                    tree_query_bound: self.total_query_span,
                    subject,
                    gap_open: self.gap_open,
                    gap_extend: self.gap_extend,
                    x_drop_final: self.x_drop_final,
                    contained_before: true,
                    matrix: self.matrix,
                },
                hsps,
                scratch,
            ),
            best_evalue: f64::MAX,
        };
        self.post_traceback_update(&mut hsp_list, subject.len());
        Ok((!hsp_list.hsps.is_empty()).then_some(hsp_list))
    }

    /// `s_HSPListPostTracebackUpdate` (no linking, no RPS, a scale factor of 1): the e-values from
    /// `kbp_gap_psi` (Spouge's finite-size correction with the query's and the subject's
    /// lengths, or without the Gumbel block the context's search space, `EvalueStats::e_value`),
    /// the e-value reap, the bit scores.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:229-250
    /// ```c
    ///       Blast_HSPListGetEvalues(program_number, query_info, subject_length,
    ///                               hsp_list, kGapped, FALSE, sbp, 0,
    ///                               scale_factor);
    ///    }
    ///
    ///    Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
    ///    ...
    ///    s_HSPListRescaleScores(hsp_list, score_params->scale_factor);
    ///    ...
    ///    Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
    /// ```
    fn post_traceback_update(&self, hsp_list: &mut BlastpHspList, subject_length: usize) {
        let query_length = self.context.aa_len as i32;
        for hsp in &mut hsp_list.hsps {
            hsp.e_value = self.stats.e_value(
                hsp.raw_score,
                self.karlin,
                query_length,
                subject_length as i32,
            );
        }
        reap_hsplist_by_evalue(hsp_list, self.expect_value);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1926-1928
        // ```c
        //       hsp->bit_score =
        //          (hsp->score*kbp[hsp->context]->Lambda -
        //           kbp[hsp->context]->logK) / NCBIMATH_LN2;
        // ```
        let log_k = self.karlin.k.ln();
        for hsp in &mut hsp_list.hsps {
            hsp.bit_score =
                (hsp.raw_score as f64 * self.karlin.lambda - log_k) / std::f64::consts::LN_2;
        }
    }
}

/// What the traceback of one list reads: the PSSM rows of the list's query, the query's letters
/// for the identities, the query bound of the containment tree, the subject, the gap costs, the
/// final X-drop, and whether each HSP is first tested for containment in the HSPs traced before
/// it (every program but RPS-BLAST, whose `Blast_TracebackFromHSPList` aligns every HSP).
pub(crate) struct TraceHsps<'a, R: PssmRows> {
    pub(crate) rows: &'a R,
    /// The query's letters without the mask (`query_blk->sequence_nomask`).
    pub(crate) query_nomask: &'a [u8],
    /// `query_info->contexts[hsp->context].query_length`.
    pub(crate) query_length: usize,
    /// `query_blk->length + 1`, the query bound of NCBI's containment tree.
    pub(crate) tree_query_bound: i32,
    /// The subject's letters.
    pub(crate) subject: &'a [u8],
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) x_drop_final: i32,
    pub(crate) contained_before: bool,
    /// The matrix of the search (`sbp->matrix`), for the positives.
    pub(crate) matrix: ScoringMatrix,
}

/// The loop of `Blast_TracebackFromHSPList` over one list's HSPs with a PSSM, then the purge of
/// common endpoints, the sort by score and the containment pass: the traced HSPs with their
/// counts, before `s_HSPListPostTracebackUpdate`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-373
/// ```c
///    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
///                                  0, (subject_length > 0 ? subject_length :
///                                  subject_blk->length / CODON_LENGTH) + 1);
/// ```
pub(crate) fn trace_hsps<R: PssmRows>(
    inputs: &TraceHsps<'_, R>,
    hsps: impl Iterator<Item = BlastpPreliminaryHsp>,
    scratch: &mut PssmGapScratch,
) -> Vec<BlastpHsp> {
    let query_nomask = inputs.query_nomask;
    let subject = inputs.subject;
    let mut tree = BlastIntervalTree::new(0, inputs.tree_query_bound, 0, subject.len() as i32 + 1);
    let mut traced: Vec<Traced> = Vec::new();
    for mut hsp in hsps {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-405
        // ```c
        //       if (program_number == eBlastTypeRpsBlast ||
        //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
        //                                hit_options->min_diag_separation)) {
        // ```
        if inputs.contained_before && tree.contains_hsp(&tree_hsp(&hsp), 0, 0) {
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-447
        // ```c
        //          if (hsp->query.gapped_start == 0 &&
        //                               hsp->subject.gapped_start == 0) {
        //             Boolean retval =
        //                BlastGetOffsetsForGappedAlignment(query, subject, sbp,
        //                    hsp, &q_start, &s_start);
        //             if (!retval)
        //             {  /* Unable to find start for this HSP */
        //                hsp_array[index] = Blast_HSPFree(hsp);
        //                continue;
        //             }
        //             hsp->query.gapped_start = q_start;
        //             hsp->subject.gapped_start = s_start;
        // ```
        if hsp.gapped_query_start == 0 && hsp.gapped_subject_start == 0 {
            let Some((q_start, s_start)) = pssm_get_offsets_for_gapped_alignment(
                inputs.rows,
                subject,
                hsp.query_start,
                hsp.query_end,
                hsp.subject_start,
                hsp.subject_end,
            ) else {
                continue;
            };
            hsp.gapped_query_start = q_start;
            hsp.gapped_subject_start = s_start;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:463-473
        // ```c
        //          if (!kTranslateSubject && !kSmithWaterman) {
        //              AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
        //                                 query_length, &start_shift);
        //              adjusted_subject = subject + start_shift;
        //              /* Shift the gapped start in HSP structure, to compensate for
        //                 a shift in the other direction later. */
        //              hsp->subject.gapped_start = s_start;
        //          }
        // ```
        let (s_start, s_length, shift) = adjust_subject_range(
            hsp.gapped_subject_start as usize,
            subject.len(),
            hsp.gapped_query_start as usize,
            inputs.query_length,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
        // ```c
        //            BLAST_GappedAlignmentWithTraceback(program_number, query,
        //                  adjusted_subject, gap_align, score_params, q_start, s_start,
        //                  query_length, adjusted_s_length,
        //                  fence_hit);
        // ```
        let aligned = pssm_gapped_alignment_with_traceback(
            inputs.rows,
            &subject[shift..],
            s_length,
            hsp.gapped_query_start as usize,
            s_start,
            inputs.gap_open,
            inputs.gap_extend,
            inputs.x_drop_final,
            scratch,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-601
        // ```c
        //            Blast_HSPUpdateWithTraceback(gap_align, hsp);
        //            ...
        //                delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
        //            ...
        //            if (!delete_hsp) {
        //               Blast_HSPAdjustSubjectOffset(hsp, start_shift);
        //               status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
        //                                          eQueryAndSubject);
        // ```
        // `Blast_HSPTest` deletes nothing (psiblast has no percent identity or minimum hit
        // length).
        let shift = shift as i32;
        hsp.raw_score = aligned.score;
        hsp.query_start = aligned.query_start;
        hsp.query_end = aligned.query_stop;
        hsp.subject_start = aligned.subject_start + shift;
        hsp.subject_end = aligned.subject_stop + shift;
        hsp.gapped_subject_start = s_start as i32 + shift;
        tree.add_hsp(tree_hsp(&hsp), 0, IndexMethod::QueryAndSubject);
        traced.push(Traced {
            hsp,
            edit_script: aligned.edit_script,
        });
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:634-692
    // ```c
    //        /* Remove any HSPs that share a starting or ending diagonal
    //           with a higher-scoring HSP. */
    //        Int4 extra_start =
    //            Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    //        ...
    //        /* Sort HSPs by score again, as the scores might have changed. */
    //        Blast_HSPListSortByScore(hsp_list);
    //
    //        /* Remove any HSPs that are contained within other HSPs.
    //           Since the list is sorted by score already, any HSP
    //           contained by a previous HSP is guaranteed to have a
    //           lower score, and may be purged. */
    //        Blast_IntervalTreeReset(tree);
    //        for (index = 0; index < hsp_list->hspcnt; index++) {
    //            hsp = hsp_array[index];
    //
    //            if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                      hit_options->min_diag_separation)) {
    //                hsp_array[index] = Blast_HSPFree(hsp);
    //            }
    //            else {
    //                status = BlastIntervalTreeAddHSP(hsp, tree, query_info,
    //                                        eQueryAndSubject);
    // ```
    // The purge deletes (`purge |= (program != eBlastTypeBlastn)`, blast_hits.c:2464) and
    // returns the HSP count, so the re-evaluation loop from `extra_start`
    // (blast_traceback.c:647-665) has nothing to do.
    purge_hsps_with_common_endpoints(&mut traced, |traced| &traced.hsp);
    traced.sort_by(|a, b| score_compare(&a.hsp, &b.hsp));
    tree.reset();
    let mut final_hsps = Vec::with_capacity(traced.len());
    for Traced { hsp, edit_script } in traced {
        let node = tree_hsp(&hsp);
        if tree.contains_hsp(&node, 0, 0) {
            continue;
        }
        tree.add_hsp(node, 0, IndexMethod::QueryAndSubject);
        final_hsps.push(final_hsp(
            &hsp,
            edit_script,
            query_nomask,
            subject,
            inputs.matrix,
        ));
    }
    final_hsps
}

/// The interval-tree node of an HSP of the context at offset 0.
fn tree_hsp(hsp: &BlastpPreliminaryHsp) -> TreeHsp {
    TreeHsp {
        query_offset: hsp.query_start,
        query_end: hsp.query_end,
        subject_offset: hsp.subject_start,
        subject_end: hsp.subject_end,
        score: hsp.raw_score,
        query_frame: 1,
        query_length: hsp.query_length as i32,
        query_context_offset: 0,
        subject_frame_sign: hsp.subject_frame.signum(),
    }
}

/// A preliminary list's HSP as the traceback reads it (0-based starts).
pub(super) fn preliminary_hsp(local: &BlastpHsp) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: local.q_start as i32 - 1,
        query_end: local.q_end as i32,
        subject_start: local.s_start as i32 - 1,
        subject_end: local.s_end as i32,
        gapped_query_start: local.gapped_q_start,
        gapped_subject_start: local.gapped_s_start,
        raw_score: local.raw_score,
        query_context: local.query_context,
        query_frame: local.query_frame,
        subject_frame: local.subject_frame,
        query_length: local.query_length,
        q_idx: local.q_idx,
        s_idx: local.s_idx,
    }
}

/// The final HSP of a traced one: 1-based starts, no composition adjustment, and the counts of
/// `Blast_HSPGetNumIdentitiesAndPositives` on the unmasked query with the search's matrix
/// (`sbp->matrix`; rpsblast's is BLOSUM62), with the gaps and mismatches that the reports read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:806-833
/// ```c
///       for (index=0; index<esp->size; index++)
///       {
///          switch (esp->op_type[index]) {
///          case eGapAlignSub:
///             for (i=0; i<esp->num[index]; i++) {
///                if (*q == *s) {
///                   num_ident++;
///                }
///                else if (NULL != matrix) {
///             	   if (matrix[*q][*s] > 0)
///             		   num_pos ++;
///                }
///                q++;
///                s++;
///             }
///             break;
///          case eGapAlignDel:
///             s += esp->num[index];
///             break;
///          case eGapAlignIns:
///             q += esp->num[index];
///             break;
/// ```
pub(super) fn final_hsp(
    hsp: &BlastpPreliminaryHsp,
    edit_script: Vec<GapEditOp>,
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
) -> BlastpHsp {
    let (mut q, mut s) = (hsp.query_start as usize, hsp.subject_start as usize);
    let (mut length, mut num_ident, mut num_positives, mut gap_opens, mut gap_letters) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    for op in &edit_script {
        match *op {
            GapEditOp::Sub(n) => {
                for (&ql, &sl) in query_nomask[q..q + n as usize]
                    .iter()
                    .zip(&subject[s..s + n as usize])
                {
                    if ql == sl {
                        num_ident += 1;
                        num_positives += 1;
                    } else if standard_score(matrix, ql, sl) > 0 {
                        num_positives += 1;
                    }
                }
                q += n as usize;
                s += n as usize;
                length += n as usize;
            }
            GapEditOp::Del(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                s += n as usize;
                length += n as usize;
            }
            GapEditOp::Ins(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                q += n as usize;
                length += n as usize;
            }
        }
    }
    BlastpHsp {
        identity: if length > 0 {
            (num_ident as f64 * 100.0) / (length as f64)
        } else {
            0.0
        },
        length,
        mismatch: length as i64 - (num_ident + gap_letters) as i64,
        gapopen: gap_opens,
        q_start: hsp.query_start as usize + 1,
        q_end: hsp.query_end as usize,
        s_start: hsp.subject_start as usize + 1,
        s_end: hsp.subject_end as usize,
        gapped_q_start: hsp.gapped_query_start,
        gapped_s_start: hsp.gapped_subject_start,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident,
        query_context: hsp.query_context,
        query_frame: hsp.query_frame,
        subject_frame: hsp.subject_frame,
        query_length: hsp.query_length,
        q_idx: hsp.q_idx,
        s_idx: hsp.s_idx,
        raw_score: hsp.raw_score,
        gap_info: Some(edit_script),
        num_positives,
        comp_adjustment_method: 0,
    }
}
