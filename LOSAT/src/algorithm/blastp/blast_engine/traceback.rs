//! BLASTP's ordinary traceback (`-comp_based_stats 0`): NCBI's `Blast_TracebackFromHSPList` on
//! each preliminary HSP list, in place of the composition-based redo (`kappa.rs`).
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
//! NCBI runs the lists of an OID batch in parallel threads (`blast_traceback.c:1505-1518`), each
//! with its own results that are merged after the loop (`blast_traceback_mt_priv.c:139-185`).
//! A list's traceback reads nothing another list writes, and the hit lists' order
//! (`s_EvalueCompareHSPLists`: e-value, score, OID) is total, so LOSAT inserts every list into one
//! hit list per query as NCBI's one thread does, at every thread count.
use anyhow::Result;

use super::super::gapalign::{
    adjust_subject_range, blast_gapped_alignment_with_traceback_with_scratch, BlastpPreliminaryHsp,
    GapAlignScratch,
};
use super::super::hsp::{reap_hsplist_by_evalue, BlastpHitList, BlastpHsp, BlastpHspList};
use super::super::x_redo_subjects::XRedoSubjects;
// The preliminary stage's helpers that this module calls are inlined everywhere or have other
// callers; it does not call `preliminary_hit_from_local_match_hsp`,
// `purge_preliminary_hits_with_common_endpoints` or `kappa::recompute_num_identities`, whose one
// caller each inlines them (a second caller would change the default path's machine code).
use super::{
    blastp_standard_score, preliminary_hits_share_query_end_endpoint,
    preliminary_hits_share_query_offset_endpoint, preliminary_query_end_compare,
    preliminary_query_offset_compare, preliminary_score_compare, preliminary_tree_hsp,
    query_nomask_sequence,
};
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::stats::{blast_spouge_stoe, BlastGumbelBlk, KarlinParams};

/// What the traceback of every list reads: the search's queries, subjects, scoring and
/// hit-saving options.
pub(super) struct OrdinaryTraceback<'a> {
    pub(super) contexts: &'a [QueryContext],
    /// `query_blk->length + 1`, the query bound of NCBI's containment tree.
    pub(super) total_query_span: i32,
    pub(super) subjects: XRedoSubjects<'a>,
    pub(super) matrix: ScoringMatrix,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    /// `ext_params->gap_x_dropoff_final` (`blast_parameters.c:457-463`).
    pub(super) x_drop_final: i32,
    /// `sbp->kbp_gap`, unscaled (no composition adjustment).
    pub(super) karlin: &'a KarlinParams,
    pub(super) gumbel: &'a BlastGumbelBlk,
    pub(super) expect_value: f64,
    pub(super) hitlist_size: usize,
}

impl OrdinaryTraceback<'_> {
    /// The final hit list of each query from the preliminary lists, which it takes.
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
    #[inline(never)]
    pub(super) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<Option<BlastpHitList>>> {
        // (query index, the list of one subject)
        let work: Vec<(usize, BlastpHspList)> = preliminary_hit_lists
            .iter_mut()
            .enumerate()
            .filter_map(|(q_idx, hit_list)| Some((q_idx, hit_list.take()?)))
            .flat_map(|(q_idx, hit_list)| {
                let count = hit_list.hsplist_count;
                hit_list
                    .hsplist_array
                    .into_iter()
                    .take(count)
                    .map(move |hsp_list| (q_idx, hsp_list))
            })
            .collect();
        let traced = self.trace_lists(work, pool)?;
        let mut hit_lists: Vec<Option<BlastpHitList>> = std::iter::repeat_with(|| None)
            .take(preliminary_hit_lists.len())
            .collect();
        for (q_idx, hsp_list) in traced {
            if let Some(hsp_list) = hsp_list {
                hit_lists[q_idx]
                    .get_or_insert_with(|| BlastpHitList::new(self.hitlist_size))
                    .update(hsp_list);
            }
        }
        Ok(hit_lists)
    }

    /// Each list's traceback, in the lists' order; a preliminary list is freed when its
    /// traceback ends, as NCBI frees each one it has read.
    fn trace_lists(
        &self,
        work: Vec<(usize, BlastpHspList)>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<(usize, Option<BlastpHspList>)>> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() && pool.threads() > 1 {
            use rayon::prelude::*;
            return pool.install(|| {
                work.into_par_iter()
                    .map_init(GapAlignScratch::new, |scratch, (q_idx, list)| {
                        Ok((q_idx, self.trace_list(q_idx, &list, scratch)?))
                    })
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = GapAlignScratch::new();
        work.into_iter()
            .map(|(q_idx, list)| Ok((q_idx, self.trace_list(q_idx, &list, &mut scratch)?)))
            .collect()
    }

    /// `Blast_TracebackFromHSPList` and `s_HSPListPostTracebackUpdate` on the list of one query
    /// and one subject; `None` when no HSP is left.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-373
    /// ```c
    ///    /* set up the tree for HSP containment tests. subject_length
    ///       is zero only for translated subject sequences, whose maximum
    ///       length is bounded by the length of the first frame */
    ///
    ///    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
    ///                                  0, (subject_length > 0 ? subject_length :
    ///                                  subject_blk->length / CODON_LENGTH) + 1);
    /// ```
    fn trace_list(
        &self,
        q_idx: usize,
        list: &BlastpHspList,
        scratch: &mut GapAlignScratch,
    ) -> Result<Option<BlastpHspList>> {
        if list.hsps.is_empty() {
            return Ok(None);
        }
        let ctx = &self.contexts[q_idx];
        let query = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
        let query_nomask = query_nomask_sequence(ctx);
        let subject_read = self.subjects.read(list.oid)?;
        let subject = subject_read.residues();
        let mut tree =
            BlastIntervalTree::new(0, self.total_query_span, 0, subject.len() as i32 + 1);
        let mut traced: Vec<Traced> = Vec::new();
        for local_hsp in &list.hsps {
            let mut hsp = preliminary_hsp(local_hsp);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-405
            // ```c
            //       if (program_number == eBlastTypeRpsBlast ||
            //           !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
            //                                hit_options->min_diag_separation)) {
            // ```
            if tree.contains_hsp(
                &preliminary_tree_hsp(&hsp, ctx.frame_base),
                ctx.frame_base,
                0,
            ) {
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
                let Some((q_start, s_start)) =
                    get_offsets_for_gapped_alignment(query, subject, &hsp, self.matrix)
                else {
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
                ctx.aa_len,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
            // ```c
            //            BLAST_GappedAlignmentWithTraceback(program_number, query,
            //                  adjusted_subject, gap_align, score_params, q_start, s_start,
            //                  query_length, adjusted_s_length,
            //                  fence_hit);
            // ```
            let Some(aligned) = blast_gapped_alignment_with_traceback_with_scratch(
                query,
                &subject[shift..shift + s_length],
                hsp.gapped_query_start as usize,
                s_start,
                self.matrix,
                None,
                self.gap_open,
                self.gap_extend,
                self.x_drop_final,
                scratch,
                None,
            ) else {
                continue;
            };
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
            // `Blast_HSPTest` deletes nothing: BLASTP has no percent identity or minimum hit
            // length (`-perc_identity` is BLASTN's).
            let shift = shift as i32;
            hsp.raw_score = aligned.score;
            hsp.query_start = aligned.query_start;
            hsp.query_end = aligned.query_stop;
            hsp.subject_start = aligned.subject_start + shift;
            hsp.subject_end = aligned.subject_stop + shift;
            hsp.gapped_subject_start = s_start as i32 + shift;
            tree.add_hsp(
                preliminary_tree_hsp(&hsp, ctx.frame_base),
                ctx.frame_base,
                IndexMethod::QueryAndSubject,
            );
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
        // For BLASTP the purge deletes (`purge |= (program != eBlastTypeBlastn)`,
        // blast_hits.c:2464) and returns the HSP count, so the re-evaluation loop from
        // `extra_start` (blast_traceback.c:647-665) has nothing to do.
        purge_hsps_with_common_endpoints(&mut traced);
        traced.sort_by(|a, b| preliminary_score_compare(&a.hsp, &b.hsp));
        tree.reset();
        let mut hsp_list = BlastpHspList {
            oid: list.oid,
            query_index: list.query_index,
            hsps: Vec::with_capacity(traced.len()),
            best_evalue: f64::MAX,
        };
        for Traced { hsp, edit_script } in traced {
            let node = preliminary_tree_hsp(&hsp, ctx.frame_base);
            if tree.contains_hsp(&node, ctx.frame_base, 0) {
                continue;
            }
            tree.add_hsp(node, ctx.frame_base, IndexMethod::QueryAndSubject);
            hsp_list.hsps.push(final_hsp(
                &hsp,
                edit_script,
                query_nomask,
                subject,
                self.matrix,
            ));
        }
        self.post_traceback_update(&mut hsp_list, ctx, subject.len());
        Ok((!hsp_list.hsps.is_empty()).then_some(hsp_list))
    }

    /// `s_HSPListPostTracebackUpdate` for BLASTP (no linking, no RPS, a scale factor of 1):
    /// the e-values from the unscaled `kbp_gap` (Spouge's finite-size correction with the
    /// query's and the subject's lengths), the e-value reap, the bit scores.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:229-250
    /// ```c
    ///       double scale_factor =
    ///          (Blast_ProgramIsRpsBlast(program_number) ?
    ///          score_params->scale_factor : 1.0);
    ///       ...
    ///       Blast_HSPListGetEvalues(program_number, query_info, subject_length,
    ///                               hsp_list, kGapped, FALSE, sbp, 0,
    ///                               scale_factor);
    ///    }
    ///
    ///    Blast_HSPListReapByEvalue(hsp_list, hit_params->options);
    ///
    ///    /* Rescale the scores by scaling factor, if necessary. This rescaling
    ///     * should be done for all programs where scaling factor is not 1.
    ///     */
    ///    s_HSPListRescaleScores(hsp_list, score_params->scale_factor);
    ///
    ///    /** Calculate and fill the bit scores. @todo: This is not the only time
    ///     * when they are calculated, s_HSPListRescaleScores also does this in
    ///     * blast_kappa.c.
    ///     */
    ///    Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);
    /// ```
    fn post_traceback_update(
        &self,
        hsp_list: &mut BlastpHspList,
        ctx: &QueryContext,
        subject_length: usize,
    ) {
        for hsp in &mut hsp_list.hsps {
            hsp.e_value = blast_spouge_stoe(
                hsp.raw_score,
                self.karlin,
                self.gumbel,
                ctx.aa_len as i32,
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

/// An HSP that the traceback aligned, with its edit script (`BlastHSP::gap_info`).
struct Traced {
    hsp: BlastpPreliminaryHsp,
    edit_script: Vec<GapEditOp>,
}

/// A preliminary list's HSP as the traceback reads it (0-based starts).
fn preliminary_hsp(local: &BlastpHsp) -> BlastpPreliminaryHsp {
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

/// `Blast_HSPListPurgeHSPsWithCommonEndpoints` with `purge` (BLASTP): of the HSPs that share a
/// start, then of those that share an end, the first in NCBI's order stays.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2480-2527
/// ```c
///    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
///    i = 0;
///    while (i < hsp_count) {
///       j = 1;
///       while (i+j < hsp_count &&
///              hsp_array[i] && hsp_array[i+j] &&
///              hsp_array[i]->context == hsp_array[i+j]->context &&
///              hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
///              hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
///              hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
///          ...
///    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
/// ```
fn purge_hsps_with_common_endpoints(hsps: &mut Vec<Traced>) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| preliminary_query_offset_compare(&a.hsp, &b.hsp));
    hsps.dedup_by(|later, kept| {
        preliminary_hits_share_query_offset_endpoint(&kept.hsp, &later.hsp)
    });
    hsps.sort_by(|a, b| preliminary_query_end_compare(&a.hsp, &b.hsp));
    hsps.dedup_by(|later, kept| preliminary_hits_share_query_end_endpoint(&kept.hsp, &later.hsp));
}

/// The final HSP of a traced one: 1-based starts, no composition adjustment, and the counts of
/// `Blast_HSPGetNumIdentitiesAndPositives` on the unmasked query with the search matrix, which
/// the traceback makes for each HSP (blast_traceback.c:586-594), with the gaps and mismatches
/// that the reports read.
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
fn final_hsp(
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
                    } else if blastp_standard_score(matrix, ql, sl) > 0 {
                        num_positives += 1;
                    }
                }
                q += n as usize;
                s += n as usize;
            }
            GapEditOp::Del(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                s += n as usize;
            }
            GapEditOp::Ins(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                q += n as usize;
            }
        }
        length += op_length(op);
    }
    BlastpHsp {
        identity: if length > 0 {
            (num_ident as f64 * 100.0) / (length as f64)
        } else {
            0.0
        },
        length,
        mismatch: length - (num_ident + gap_letters),
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

fn op_length(op: &GapEditOp) -> usize {
    match *op {
        GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => n as usize,
    }
}

/// NCBI's `BlastGetOffsetsForGappedAlignment` (not position-based): the middle of the HSP's
/// best-scoring 11-letter window, or the end window, or none when both score 0 or less.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3259-3320
/// ```c
///     if (q_length <= HSP_MAX_WINDOW) {
///         *q_retval = q_start + q_length/2;
///         *s_retval = s_start + q_length/2;
///         return TRUE;
///     }
///     ...
///     if (max_score > 0)
///     {
///         *q_retval = max_offset;
///         *s_retval = (max_offset - q_start) + s_start;
///         return TRUE;
///     }
///     else  /* Test the window around the ends of the HSP. */
///     {
///         ...
///         if (score > 0)
///         {
///             *q_retval = hsp->query.end - HSP_MAX_WINDOW/2;
///             *s_retval = hsp->subject.end - HSP_MAX_WINDOW/2;
///             return TRUE;
///         }
///     }
///     return FALSE;
/// ```
fn get_offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    hsp: &BlastpPreliminaryHsp,
    matrix: ScoringMatrix,
) -> Option<(i32, i32)> {
    const HSP_MAX_WINDOW: i32 = 11;
    let q_start = hsp.query_start;
    let s_start = hsp.subject_start;
    let q_length = hsp.query_end - q_start;
    let s_length = hsp.subject_end - s_start;
    if q_length <= HSP_MAX_WINDOW {
        return Some((q_start + q_length / 2, s_start + q_length / 2));
    }
    let score_at =
        |q: i32, s: i32| blastp_standard_score(matrix, query[q as usize], subject[s as usize]);
    let mut score: i32 = (0..HSP_MAX_WINDOW)
        .map(|i| score_at(q_start + i, s_start + i))
        .sum();
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    for index in q_start + HSP_MAX_WINDOW..q_start + q_length.min(s_length) {
        let s_index = index - q_start + s_start;
        score -= score_at(index - HSP_MAX_WINDOW, s_index - HSP_MAX_WINDOW);
        score += score_at(index, s_index);
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
    }
    if max_score > 0 {
        return Some((max_offset, max_offset - q_start + s_start));
    }
    let end_score: i32 = (0..HSP_MAX_WINDOW)
        .map(|i| {
            score_at(
                q_start + q_length - HSP_MAX_WINDOW + i,
                s_start + s_length - HSP_MAX_WINDOW + i,
            )
        })
        .sum();
    (end_score > 0).then_some((
        hsp.query_end - HSP_MAX_WINDOW / 2,
        hsp.subject_end - HSP_MAX_WINDOW / 2,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBIstdaa: A 1, C 3, W 20 (BLOSUM62: W/W 11, A/A 4, W/C -2, A/W -3).
    const A: u8 = 1;
    const C: u8 = 3;
    const W: u8 = 20;

    fn hsp(q: (i32, i32), s: (i32, i32)) -> BlastpPreliminaryHsp {
        BlastpPreliminaryHsp {
            query_start: q.0,
            query_end: q.1,
            subject_start: s.0,
            subject_end: s.1,
            gapped_query_start: 0,
            gapped_subject_start: 0,
            raw_score: 0,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
        }
    }

    #[test]
    fn short_hsp_starts_at_its_middle() {
        let seq = [W; 20];
        assert_eq!(
            get_offsets_for_gapped_alignment(
                &seq,
                &seq,
                &hsp((3, 14), (5, 16)),
                ScoringMatrix::Blosum62
            ),
            Some((8, 10))
        );
    }

    #[test]
    fn long_hsp_starts_at_the_end_of_its_first_best_window() {
        // Every window scores 121: the first one's last letter.
        let seq = [W; 30];
        assert_eq!(
            get_offsets_for_gapped_alignment(
                &seq,
                &seq,
                &hsp((2, 22), (4, 24)),
                ScoringMatrix::Blosum62
            ),
            Some((12, 14))
        );
        // A run of A then W: the windows grow until the first all-W one (ending at 21); the
        // later ones only tie.
        let mut seq = [A; 30];
        seq[11..].fill(W);
        assert_eq!(
            get_offsets_for_gapped_alignment(
                &seq,
                &seq,
                &hsp((0, 25), (0, 25)),
                ScoringMatrix::Blosum62
            ),
            Some((21, 21))
        );
    }

    #[test]
    fn hsp_without_a_positive_window_has_no_start() {
        let query = [W; 30];
        let subject = [C; 30];
        assert_eq!(
            get_offsets_for_gapped_alignment(
                &query,
                &subject,
                &hsp((0, 20), (0, 20)),
                ScoringMatrix::Blosum62
            ),
            None
        );
    }
}
