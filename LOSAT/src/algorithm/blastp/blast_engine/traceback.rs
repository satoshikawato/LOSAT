//! BLASTP's ordinary traceback (`-comp_based_stats 0`): NCBI's `Blast_TracebackFromHSPList` on
//! each preliminary HSP list, in place of the composition-based redo (`kappa.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1500
//!
//! NCBI's threads take the OID batches of lists in turn (`c++/src/algo/blast/core/blast_traceback.c:1505-1518`), each
//! thread with its own results that are merged after the loop (`blast_traceback_mt_priv.c:139-185`).
//! A list's traceback reads nothing another list writes, and the hit lists' order
//! (`s_EvalueCompareHSPLists`: e-value, score, OID) is total, so LOSAT inserts every list into one
//! hit list per query as NCBI's one thread does, at every thread count (NCBI's merge of the
//! threads' capped results can keep more lists for `-culling_limit` and `-qcov_hsp_perc` with
//! several threads; LOSAT's output is that of NCBI's `-num_threads 1`).
//!
//! The steps after any traceback on each query's final hit list (`finish_hit_lists`) are here too.
use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;

use super::super::args::ResolvedBlastpArgs;
use super::super::gapalign::{
    adjust_subject_range, blast_gapped_alignment_with_traceback_with_scratch, BlastpPreliminaryHsp,
    GapAlignScratch,
};
use super::super::hsp::{
    reap_hsplist_by_evalue, trim_by_max_hsps, update_best_evalue, BlastpHitList, BlastpHsp,
    BlastpHspList,
};
use super::super::x_prelim_compact::XPrelimHsp;
use super::super::x_redo_subjects::XRedoSubjects;
use super::{besthit, culling};
// The preliminary stage's helpers that this module calls are inlined everywhere or have other
// callers; it does not call `preliminary_hit_from_local_match_hsp`,
// `purge_preliminary_hits_with_common_endpoints` or `kappa::recompute_num_identities`, whose one
// caller each inlines them (a second caller would change the default path's machine code).
use super::{
    preliminary_hits_share_query_end_endpoint, preliminary_hits_share_query_offset_endpoint,
    preliminary_query_end_compare, preliminary_query_offset_compare, preliminary_score_compare,
    preliminary_tree_hsp, query_nomask_sequence,
};
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::core::blast_hits::reap_by_query_coverage;
use crate::stats::karlin::evalue_from_raw_score;
use crate::stats::{blast_spouge_stoe, BlastGumbelBlk, KarlinParams};
use crate::utils::matrix::standard_score;

/// What the traceback of every list reads: the search's queries, subjects, scoring and
/// hit-saving options.
pub(super) struct OrdinaryTraceback<'a> {
    pub(super) contexts: &'a [QueryContext],
    /// One more than `query_blk->length`, the query bound of NCBI's containment tree.
    pub(super) total_query_span: i32,
    pub(super) subjects: XRedoSubjects<'a>,
    pub(super) matrix: ScoringMatrix,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    /// `ext_params->gap_x_dropoff_final` (`blast_parameters.c:457-463`).
    pub(super) x_drop_final: i32,
    /// `sbp->kbp_gap`, unscaled (no composition adjustment).
    pub(super) karlin: &'a KarlinParams,
    /// None with `-searchsp`: the E-values then come from the context's search space.
    pub(super) gumbel: Option<&'a BlastGumbelBlk>,
    pub(super) search_spaces: &'a [crate::stats::search_space::SearchSpace],
    pub(super) expect_value: f64,
    pub(super) hitlist_size: usize,
    /// Whether NCBI sorts the results at the end of the traceback (`DbStatistics::sorts_results`):
    /// otherwise the hit lists keep the order in which the lists come, by OID (`run`).
    pub(super) sorts_results: bool,
    /// psiblast's `-sum_stats true`: the traceback's linking parameters
    /// (`EvenGapLinkParams::traceback`), `None` without sum statistics.
    pub(super) link: Option<crate::core::link_hsps::EvenGapLinkParams>,
}

impl OrdinaryTraceback<'_> {
    /// The final hit list of each query from the preliminary lists, which it takes: a query's
    /// lists are full (`preliminary_hit_lists`) or compact (`compact_hit_lists`, which the
    /// traceback reads without expanding them).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1690-1706
    #[inline(never)]
    pub(super) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let work = prelim_lists_of_queries(preliminary_hit_lists, compact_hit_lists);
        let traced = self.trace_lists(work, pool)?;
        Ok(hit_lists_of_traced(
            traced,
            preliminary_hit_lists.len(),
            self.hitlist_size,
            self.sorts_results,
        ))
    }

    /// Each list's traceback, in the lists' order; a preliminary list is freed when its
    /// traceback ends, as NCBI frees each one it has read.
    fn trace_lists(
        &self,
        work: Vec<(usize, PrelimList)>,
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
                        Ok((q_idx, list.trace(self, q_idx, scratch)?))
                    })
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = GapAlignScratch::new();
        work.into_iter()
            .map(|(q_idx, list)| Ok((q_idx, list.trace(self, q_idx, &mut scratch)?)))
            .collect()
    }

    /// `Blast_TracebackFromHSPList` and `s_HSPListPostTracebackUpdate` on the list of one query
    /// and one subject; `None` when no HSP is left.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-373
    fn trace_list(
        &self,
        q_idx: usize,
        oid: u32,
        query_index: u32,
        hsps: impl Iterator<Item = BlastpPreliminaryHsp>,
        scratch: &mut GapAlignScratch,
    ) -> Result<Option<BlastpHspList>> {
        let ctx = &self.contexts[q_idx];
        let query = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
        let query_nomask = query_nomask_sequence(ctx);
        let subject_read = self.subjects.read(oid)?;
        let subject = subject_read.residues();
        let mut tree =
            BlastIntervalTree::new(0, self.total_query_span, 0, subject.len() as i32 + 1);
        let mut traced: Vec<Traced> = Vec::new();
        for mut hsp in hsps {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-405
            if tree.contains_hsp(
                &preliminary_tree_hsp(&hsp, ctx.frame_base),
                ctx.frame_base,
                0,
            ) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-447
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
            let (s_start, s_length, shift) = adjust_subject_range(
                hsp.gapped_subject_start as usize,
                subject.len(),
                hsp.gapped_query_start as usize,
                ctx.aa_len,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
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
        // For BLASTP the purge deletes (`purge` is set when `program` is not `eBlastTypeBlastn`,
        // blast_hits.c:2464) and returns the HSP count, so the re-evaluation loop from
        // `extra_start` (blast_traceback.c:647-665) has nothing to do.
        purge_hsps_with_common_endpoints(&mut traced);
        traced.sort_by(|a, b| preliminary_score_compare(&a.hsp, &b.hsp));
        tree.reset();
        let mut hsp_list = BlastpHspList {
            oid,
            query_index,
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

    /// `s_HSPListPostTracebackUpdate` for BLASTP (no RPS, a scale factor of 1): the e-values
    /// from the unscaled `kbp_gap` (Spouge's finite-size correction with the query's and the
    /// subject's lengths), or with psiblast's `-sum_stats true` the linked sets' e-values (large
    /// gaps only), the e-value reap, the bit scores.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:224-250
    fn post_traceback_update(
        &self,
        hsp_list: &mut BlastpHspList,
        ctx: &QueryContext,
        subject_length: usize,
    ) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:213-216
        // psiblast's `-sum_stats true`: the linked sets' e-values (large gaps only,
        // `EvenGapLinkParams::traceback`).
        if let Some(link) = self.link.as_ref() {
            crate::algorithm::blastp::hsp::link_hsp_list(
                hsp_list,
                super::blastp_link_context(
                    self.karlin,
                    ctx.aa_len,
                    &self.search_spaces[hsp_list.query_index as usize],
                ),
                subject_length as i32,
                link,
            );
        } else {
            self.get_evalues(hsp_list, ctx, subject_length);
        }
        reap_hsplist_by_evalue(hsp_list, self.expect_value);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1926-1928
        let log_k = self.karlin.k.glibc_ln();
        for hsp in &mut hsp_list.hsps {
            hsp.bit_score =
                (hsp.raw_score as f64 * self.karlin.lambda - log_k) / std::f64::consts::LN_2;
        }
    }

    /// `Blast_HSPListGetEvalues` without linking: the e-values from the unscaled `kbp_gap`
    /// (Spouge's finite-size correction with the query's and the subject's lengths).
    fn get_evalues(&self, hsp_list: &mut BlastpHspList, ctx: &QueryContext, subject_length: usize) {
        match self.gumbel {
            Some(gumbel) => {
                for hsp in &mut hsp_list.hsps {
                    hsp.e_value = blast_spouge_stoe(
                        hsp.raw_score,
                        self.karlin,
                        gumbel,
                        ctx.aa_len as i32,
                        subject_length as i32,
                    );
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1886-1891
            None => {
                let search_space =
                    self.search_spaces[hsp_list.query_index as usize].effective_space;
                for hsp in &mut hsp_list.hsps {
                    hsp.e_value = evalue_from_raw_score(hsp.raw_score, self.karlin, search_space);
                }
            }
        }
    }
}

/// The preliminary lists of every query as the final stage reads them: a query's lists are full
/// (`preliminary_hit_lists`) or compact (`compact_hit_lists`, which are taken without being
/// expanded); each with its query index.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1636-1639
fn prelim_lists_of_queries(
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> Vec<(usize, PrelimList)> {
    // (query index, the list of one subject)
    let mut work: Vec<(usize, PrelimList)> = Vec::new();
    for (q_idx, (full, compact)) in preliminary_hit_lists
        .iter_mut()
        .zip(compact_hit_lists.iter_mut())
        .enumerate()
    {
        if let Some(hit_list) = compact.take() {
            let count = hit_list.hsplist_count;
            work.extend(
                hit_list
                    .hsplist_array
                    .into_iter()
                    .take(count)
                    .map(|list| (q_idx, PrelimList::Compact(list))),
            );
        } else if let Some(hit_list) = full.take() {
            let count = hit_list.hsplist_count;
            work.extend(
                hit_list
                    .hsplist_array
                    .into_iter()
                    .take(count)
                    .map(|list| (q_idx, PrelimList::Full(list))),
            );
        }
    }
    work
}

/// The final hit list of each of `queries` queries from the final lists of `traced`, inserted as
/// NCBI's results are (`Blast_HSPResultsInsertHSPList` with `hitlist_size`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:591-608 (cited
/// at `HitList::put_in_traceback_order` of BLASTN): the traceback reads the preliminary
/// lists by OID. Where the results are sorted at the end, the order of insertion changes
/// nothing (the comparison is total); where they are not, it is the report's order.
pub(super) fn hit_lists_of_traced(
    mut traced: Vec<(usize, Option<BlastpHspList>)>,
    queries: usize,
    hitlist_size: usize,
    sorts_results: bool,
) -> Vec<Option<BlastpHitList>> {
    if !sorts_results {
        traced.sort_by_key(|(_, list)| list.as_ref().map_or(u32::MAX, |list| list.oid));
    }
    let mut hit_lists: Vec<Option<BlastpHitList>> =
        std::iter::repeat_with(|| None).take(queries).collect();
    for (q_idx, hsp_list) in traced {
        if let Some(hsp_list) = hsp_list {
            hit_lists[q_idx]
                .get_or_insert_with(|| BlastpHitList::new(hitlist_size))
                .update(hsp_list);
        }
    }
    hit_lists
}

/// What the final stage of an ungapped search reads (`-ungapped`: no traceback).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1685-1704
pub(super) struct UngappedFinish<'a> {
    pub(super) contexts: &'a [QueryContext],
    pub(super) subjects: XRedoSubjects<'a>,
    /// The search matrix, for the positives of each HSP.
    pub(super) matrix: ScoringMatrix,
    pub(super) hitlist_size: usize,
    /// Whether NCBI sorts the results at the end of the traceback (`DbStatistics::sorts_results`).
    pub(super) sorts_results: bool,
}

impl UngappedFinish<'_> {
    /// The final hit list of each query from the preliminary lists, which it takes: each HSP is
    /// the ungapped hit of the preliminary stage (one run of aligned letters, its E-value kept),
    /// with its identities and positives counted on the unmasked query and its bit score from the
    /// context's own block (`kbp_std`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1918-1928
    /// An HSP of length 0 has no letter to align: it goes on to the lists of the report, which
    /// write it as NCBI does (`-ungapped` with a low `-threshold`, `-window_size 0`).
    pub(super) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let work = prelim_lists_of_queries(preliminary_hit_lists, compact_hit_lists);
        let mut finished = Vec::with_capacity(work.len());
        for (q_idx, list) in work {
            let (oid, query_index, hsps): (u32, u32, Vec<(BlastpPreliminaryHsp, f64)>) = match list
            {
                PrelimList::Full(list) => (
                    list.oid,
                    list.query_index,
                    list.hsps
                        .iter()
                        .map(|hsp| (preliminary_hsp(hsp), hsp.e_value))
                        .collect(),
                ),
                PrelimList::Compact(list) => (
                    list.oid,
                    list.query_index,
                    list.hsps.iter().map(|hsp| (hsp.hsp, hsp.e_value)).collect(),
                ),
            };
            finished.push((q_idx, self.finish_list(q_idx, oid, query_index, &hsps)?));
        }
        Ok(hit_lists_of_traced(
            finished,
            preliminary_hit_lists.len(),
            self.hitlist_size,
            self.sorts_results,
        ))
    }

    /// The final list of one query and one subject; `None` when it has no HSP.
    fn finish_list(
        &self,
        q_idx: usize,
        oid: u32,
        query_index: u32,
        hsps: &[(BlastpPreliminaryHsp, f64)],
    ) -> Result<Option<BlastpHspList>> {
        let ctx = &self.contexts[q_idx];
        let query_nomask = query_nomask_sequence(ctx);
        let subject_read = self.subjects.read(oid)?;
        let subject = subject_read.residues();
        let log_k = ctx.karlin_params.k.glibc_ln();
        let mut hsp_list = BlastpHspList {
            oid,
            query_index,
            hsps: Vec::with_capacity(hsps.len()),
            best_evalue: f64::MAX,
        };
        for &(hsp, e_value) in hsps {
            // The kernel deleted the HSPs of a negative length (`Blast_HSPTest`); one of length 0
            // is `Sub(0)`.
            let length = hsp.query_end - hsp.query_start;
            debug_assert!(length >= 0, "ungapped_subject deletes a negative length");
            let mut final_hsp = final_hsp(
                &hsp,
                vec![GapEditOp::Sub(length as u32)],
                query_nomask,
                subject,
                self.matrix,
            );
            final_hsp.e_value = e_value;
            final_hsp.bit_score =
                (hsp.raw_score as f64 * ctx.karlin_params.lambda - log_k) / std::f64::consts::LN_2;
            hsp_list.hsps.push(final_hsp);
        }
        Ok((!hsp_list.hsps.is_empty()).then_some(hsp_list))
    }
}

/// A preliminary list of one query and one subject, full or compact.
enum PrelimList {
    Full(BlastpHspList),
    Compact(BlastpHspList<XPrelimHsp>),
}

impl PrelimList {
    fn trace(
        self,
        traceback: &OrdinaryTraceback<'_>,
        q_idx: usize,
        scratch: &mut GapAlignScratch,
    ) -> Result<Option<BlastpHspList>> {
        match self {
            Self::Full(list) => traceback.trace_list(
                q_idx,
                list.oid,
                list.query_index,
                list.hsps.iter().map(preliminary_hsp),
                scratch,
            ),
            Self::Compact(list) => traceback.trace_list(
                q_idx,
                list.oid,
                list.query_index,
                list.hsps.iter().map(|hsp| hsp.hsp),
                scratch,
            ),
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2527
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:796-833
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
                    } else if standard_score(matrix, ql, sl) > 0 {
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
        num: 1,
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

fn op_length(op: &GapEditOp) -> usize {
    match *op {
        GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => n as usize,
    }
}

/// NCBI's `BlastGetOffsetsForGappedAlignment` (not position-based): the middle of the HSP's
/// best-scoring 11-letter window, or the end window, or none when both score 0 or less.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3259-3320
///
/// The windows read 11 letters from the HSP's starts and ends. An HSP shorter than that in the
/// subject (only longer in the query, through gaps) would make NCBI read outside the HSP and
/// possibly outside the subject; BLASTP's preliminary gapped starts are never both 0 (the start
/// is the middle of a window), so this is not reached, and a window outside a sequence gives no
/// start instead of reading there.
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
    let score_at = |q: i32, s: i32| -> Option<i32> {
        let q = query.get(usize::try_from(q).ok()?)?;
        let s = subject.get(usize::try_from(s).ok()?)?;
        Some(standard_score(matrix, *q, *s))
    };
    let mut score = 0;
    for i in 0..HSP_MAX_WINDOW {
        score += score_at(q_start + i, s_start + i)?;
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    for index in q_start + HSP_MAX_WINDOW..q_start + q_length.min(s_length) {
        let s_index = index - q_start + s_start;
        score -= score_at(index - HSP_MAX_WINDOW, s_index - HSP_MAX_WINDOW)?;
        score += score_at(index, s_index)?;
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
    }
    if max_score > 0 {
        return Some((max_offset, max_offset - q_start + s_start));
    }
    let mut end_score = 0;
    for i in 0..HSP_MAX_WINDOW {
        end_score += score_at(
            q_start + q_length - HSP_MAX_WINDOW + i,
            s_start + s_length - HSP_MAX_WINDOW + i,
        )?;
    }
    (end_score > 0).then_some((
        hsp.query_end - HSP_MAX_WINDOW / 2,
        hsp.subject_end - HSP_MAX_WINDOW / 2,
    ))
}

/// The steps after the traceback on each query's final hit list (`run_call`, after the redo or
/// the ordinary traceback), in the order that the comments below cite: `-culling_limit`'s or
/// `-best_hit_*`'s pipe, `-max_hsps`, `-qcov_hsp_perc` and `-subject_besthit`, then the e-value
/// sort (not when `sorts_results` is false) and the cut to the traceback's hit list size (neither
/// for rpsblast over several volumes, `rps_multi_volume`).
pub(super) fn finish_hit_lists(
    args: &ResolvedBlastpArgs,
    hit_lists: &mut [Option<BlastpHitList>],
    prelim_hitlist_size: usize,
    rps_multi_volume: bool,
    sorts_results: bool,
) {
    for hit_list_opt in hit_lists.iter_mut() {
        if let Some(hit_list) = hit_list_opt.as_mut() {
            // `-culling_limit`: the culling pipe is the first step after the traceback
            // (`BlastHSPStreamTBackClose`, blast_traceback.c:1721), before `-max_hsps` and
            // `-qcov_hsp_perc` (blast_traceback.c:1763-1766), before the e-value sort and the
            // prune, and before the `-query_loc` shift of the coordinates (`run_call`). The best
            // hit pipe of `-best_hit_*` is in the same place (the options exclude each other).
            if args.culling_limit > 0 {
                culling::cull_traceback(hit_list, args.culling_limit);
            } else if let Some(params) = args.best_hit {
                besthit::best_hit_traceback(hit_list, params, prelim_hitlist_size);
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:845-865
            // The e-value sort below reads the best e-value and the first HSP of each list,
            // which `-max_hsps` keeps, so it follows the cut as NCBI's does.
            if args.max_hsps_per_subject > 0
                || args.query_cov_hsp_perc > 0.0
                || args.subject_besthit
            {
                for hsp_list in hit_list
                    .hsplist_array
                    .iter_mut()
                    .take(hit_list.hsplist_count)
                {
                    trim_by_max_hsps(hsp_list, args.max_hsps_per_subject);
                    if reap_by_query_coverage(&mut hsp_list.hsps, args.query_cov_hsp_perc, |hsp| {
                        (
                            hsp.q_start as i64 - 1,
                            hsp.q_end as i64,
                            hsp.query_length as i64,
                        )
                    }) {
                        update_best_evalue(hsp_list);
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:858-861
                    // `-subject_besthit` (one context per query; the best e-value stays).
                    if args.subject_besthit {
                        crate::core::blast_hits::subject_besthit(
                            &mut hsp_list.hsps,
                            false,
                            |hsp| crate::core::blast_hits::BestHitExtent {
                                context: 0,
                                query_frame: 0,
                                offset: hsp.q_start as i64 - 1,
                                end: hsp.q_end as i64,
                                query_length: 0,
                            },
                        );
                    }
                }
                if args.query_cov_hsp_perc > 0.0 {
                    hit_list.purge_empty_hsp_lists();
                }
            }
            // (rpsblast over several volumes: the merged lists keep `s_MergeAlignSet`'s order
            // and the formatter cuts them, `psiblast_hits_found` in `run_call`.)
            if !rps_multi_volume {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1772
                // (cited at `DbStatistics::sorts_results`): not for a database length of 0 or of
                // 2^63 or more, whose lists keep the traceback's order (by OID with
                // `-comp_based_stats 0`, the composition heap's otherwise; rpsblast's engine
                // sorted by the same gate, `RpsUnit::sorts_results`).
                if sorts_results {
                    hit_list.sort_by_evalue();
                }
                hit_list.prune_by_size(args.traceback_hitlist_size());
            }
        }
    }
}

#[cfg(test)]
#[path = "traceback_tests.rs"]
mod tests;
