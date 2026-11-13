//! The preliminary stage of one user query against one profile volume
//! (`s_BlastSearchEngineCore` with the concatenated profiles as the query): the word finder,
//! `BLAST_GetGappedScore` over the profiles (one context each), the purge of common endpoints,
//! the sort by score, the e-values and the reap by the preliminary e-value, then the collector's
//! split into one list per profile (`s_BlastHSPCollectorRun_RPS`) and the query's hit list.

use anyhow::{bail, Result};

use super::wordfinder::{word_finder, RpsInitHsp, RpsWordParams, RpsWordScratch};
use super::ProfileSlice;
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHspList};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::gapalign::{
    pssm_get_start_for_gapped_alignment, pssm_score_only_gapped_alignment, PssmGapScratch,
};
use crate::algorithm::pssm_search::prelim::{purge_hsps_with_common_endpoints, score_compare};
use crate::stats::KarlinParams;

/// What the preliminary stage of every query of a volume reads besides the query.
pub(crate) struct RpsPrelim<'a> {
    pub(crate) word: RpsWordParams<'a>,
    /// The concatenated rows (`lookup->rps_pssm`).
    pub(crate) rows: &'a [[i32; 28]],
    /// The profiles' row offsets (`lookup->rps_seq_offsets`, `num_profiles + 1` of them).
    pub(crate) offsets: &'a [i32],
    /// `score_params->gap_open` and `gap_extend` (the options' costs times the scale).
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// `ext_params->gap_x_dropoff` (times the scale when it is above 1).
    pub(crate) x_drop_gapped: i32,
    /// `hit_params->cutoffs[query].cutoff_score` (times the scale).
    pub(crate) hit_cutoff: i32,
    /// The first gapped Karlin block of the batch (`kbp_gap_std`, BLOSUM62 11/1's table values)
    /// with its Lambda divided by the scale, as `Blast_HSPListGetEvalues` divides it.
    pub(crate) karlin: KarlinParams,
    /// The Gumbel block, if any, and the query's search space (`EvalueStats`): the e-values and
    /// the reap by `prelim_evalue` run with the block only.
    pub(crate) stats: EvalueStats<'a>,
    /// `hit_params->prelim_evalue`.
    pub(crate) prelim_evalue: f64,
    /// `prelim_hitlist_size` of the collector.
    pub(crate) prelim_hitlist_size: usize,
    /// The query's index (the lists' `query_index`).
    pub(crate) query_index: u32,
    /// The program of LOSAT's refusals.
    pub(crate) program: &'static str,
}

/// The memory of one worker's preliminary stage, kept from query to query.
pub(crate) struct RpsPrelimScratch {
    pub(crate) word: RpsWordScratch,
    init_hsps: Vec<RpsInitHsp>,
    tree: BlastIntervalTree,
    hsps: Vec<BlastpPreliminaryHsp>,
    gap: PssmGapScratch,
    /// The query's HSPs of its contexts so far (`finish`).
    combined: Vec<BlastpPreliminaryHsp>,
}

impl RpsPrelimScratch {
    pub(crate) fn new(word: RpsWordScratch) -> Self {
        Self {
            word,
            init_hsps: Vec::new(),
            tree: BlastIntervalTree::new(0, 1, 0, 1),
            hsps: Vec::new(),
            gap: PssmGapScratch::new(),
            combined: Vec::new(),
        }
    }
}

/// `BSearchContextInfo` on the profiles' offsets (no length bounds are set for the
/// concatenation's query information, so the search covers every context): the last profile
/// that starts at or before `n`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
/// ```c
/// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
/// {
///     Int4 m=0, b=0, e=0, size=0;
///
///     size = A->last_context+1;
///
///     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
///         ...
///     }
///     else {
///         b = 0;
///         e = size;
///     }
///
///     while (b < e - 1) {
/// 	m = (b + e) / 2;
/// 	if (A->contexts[m].query_offset > n)
/// 	    e = m;
/// 	else
/// 	    b = m;
///     }
///     return b;
/// }
/// ```
pub(crate) fn context_of(offsets: &[i32], n: i32) -> usize {
    let size = offsets.len() - 1;
    let (mut b, mut e) = (0usize, size);
    while b + 1 < e {
        let m = (b + e) / 2;
        if offsets[m] > n {
            e = m;
        } else {
            b = m;
        }
    }
    b
}

/// `score_compare_match` (`Blast_InitHitListSortByScore`), stable as glibc's qsort.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-292
/// ```c
///    if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
///                                 h1->ungapped_data->score)) &&
///        0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
///                                 h2->ungapped_data->s_start)) &&
///        0 == (result = BLAST_CMP(h2->ungapped_data->length,
///                                 h1->ungapped_data->length)) &&
///        0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
///                                 h2->ungapped_data->q_start))) {
///        result = BLAST_CMP(h2->ungapped_data->length,
///                           h1->ungapped_data->length);
///    }
/// ```
fn sort_init_hsps_by_score(init_hsps: &mut [RpsInitHsp]) {
    init_hsps.sort_by(|a, b| {
        b.ungapped
            .score
            .cmp(&a.ungapped.score)
            .then_with(|| a.ungapped.s_start.cmp(&b.ungapped.s_start))
            .then_with(|| b.ungapped.length.cmp(&a.ungapped.length))
            .then_with(|| a.ungapped.q_start.cmp(&b.ungapped.q_start))
    });
}

impl RpsPrelim<'_> {
    /// The preliminary HSPs of one context of the query on this volume (`subject_seq`: its
    /// letters and the trailing sentinel, `subject->sequence`; `frame`: its frame, 0 for a protein
    /// query), sorted by score, after those of the query's earlier contexts in `scratch`
    /// (`finish`), as `s_BlastSearchEngineOneContext` gives them to `Blast_HSPListAppend`.
    pub(crate) fn search_frame(
        &self,
        scratch: &mut RpsPrelimScratch,
        subject_seq: &[u8],
        frame: i8,
    ) -> Result<()> {
        let subject_len = subject_seq.len() - 1;
        let subject = &subject_seq[..subject_len];
        let RpsPrelimScratch {
            word,
            init_hsps,
            tree,
            hsps,
            gap,
            combined,
        } = scratch;
        init_hsps.clear();
        word_finder(&self.word, word, subject, init_hsps);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:497
        // ```c
        //             if (init_hitlist->total == 0) continue;
        // ```
        if init_hsps.is_empty() {
            return Ok(());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:282
        // ```c
        //     Blast_InitHitListSortByScore(init_hitlist);
        // ```
        sort_init_hsps_by_score(init_hsps);
        hsps.clear();
        self.gapped_score(init_hsps, tree, hsps, gap, subject_seq, subject_len)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:541-555
        // ```c
        //             /* Removes redundant HSPs. */
        //             Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
        //         ...
        //         Blast_HSPListSortByScore(hsp_list);
        // ```
        purge_hsps_with_common_endpoints(hsps, |hsp| hsp);
        if hsps
            .windows(2)
            .any(|pair| score_compare(&pair[0], &pair[1]) == std::cmp::Ordering::Greater)
        {
            hsps.sort_by(score_compare);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
        // ```c
        //                 Blast_HSPListAppend(&hsp_list, &hsp_list_out, kHspNumMax);
        // ```
        // Without a cap (`hsp_num_max` 0: `INT4_MAX`, blast_hits.c:213-226) the lists of the
        // frames are concatenated in the frame order (`s_BlastHSPListsCombineByScore`).
        combined.extend(hsps.iter().map(|hsp| BlastpPreliminaryHsp {
            subject_frame: frame.into(),
            ..*hsp
        }));
        Ok(())
    }

    /// The query's preliminary hit list on this volume from the HSPs of its contexts
    /// (`search_frame`; emptied here), with the e-values of the length `stat_length`
    /// (`RpsQuery::stat_length`), or `None` when it has no HSP. The lists' oids are the profiles'
    /// local ordinal ids.
    pub(crate) fn finish(
        &self,
        scratch: &mut RpsPrelimScratch,
        stat_length: usize,
    ) -> Option<BlastpHitList<XPrelimHsp>> {
        let hsps = &mut scratch.combined;
        if hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:882-899
        // ```c
        //         double scale_factor = 1.0;
        //         if (isRPS) {
        //             scale_factor = score_params->scale_factor;
        //         }
        //         Blast_HSPListGetEvalues(program_number, query_info,
        //                                          stat_length, hsp_list_out,
        //                                          score_options->gapped_calculation,
        //                                          isRPS, gap_align->sbp, 0, scale_factor);
        //     ...
        //        status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1876-1883
        // ```c
        //           } else {
        //               /* for RPS blast, query and subject is swapped */
        //               hsp->evalue =
        //                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
        //                                subject_length,
        //                                query_info->contexts[hsp->context].query_length);
        //           }
        // ```
        // `subject_length` is the user query's length (`stat_length`), the context's length the
        // profile's.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:874-880
        // ```c
        //     } else if (!Blast_ProgramIsPhiBlast(program_number)
        //            && !(isRPS && !sbp->gbp)
        //            /* do not calculate E-values for mapping */
        //            && program_number != eBlastTypeMapping ) {
        //         /* Calculate e-values for all HSPs. Skip this step
        //            for PHI or RPS with old FSC, since calculating the E values
        //            requires precomputation that has not been done yet */
        // ```
        // Without the Gumbel block (`-searchsp` at any value, `OLD_FSC`) the HSPs keep their
        // e-value of 0; the reap by `prelim_evalue` (blast_engine.c:896-899, outside that branch)
        // then keeps every one, and the hit list orders its lists by score, then ordinal id
        // (every best e-value is equal).
        let mut evalued: Vec<XPrelimHsp> = Vec::with_capacity(hsps.len());
        for hsp in hsps.iter() {
            let e_value = if self.stats.gumbel.is_some() {
                self.stats.e_value(
                    hsp.raw_score,
                    &self.karlin,
                    stat_length as i32,
                    hsp.query_length as i32,
                )
            } else {
                0.0
            };
            if e_value <= self.prelim_evalue {
                evalued.push(XPrelimHsp { hsp: *hsp, e_value });
            }
        }
        hsps.clear();
        if evalued.is_empty() {
            return None;
        }
        Some(self.collect(evalued))
    }

    /// `s_BlastHSPCollectorRun_RPS`: the HSPs sorted by context, then by score, one list per
    /// profile (`oid` the context, every HSP's context 0), each inserted into the query's hit
    /// list in that order.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:228-262
    /// ```c
    ///    qsort(hsplist_in->hsp_array, hsplist_in->hspcnt, sizeof(BlastHSP*),
    ///          s_ScoreCompareHSPWithContext);
    ///
    ///    /* Sequentially extract HSPs corresponding to one subject into a new
    ///       HSPList, and save these new HSPLists in a normal way, as in all other
    ///       BLAST programs. */
    ///    next_index = 0;
    ///
    ///    for (index = 0; index < hsplist_in->hspcnt; index = next_index) {
    ///        BlastHSPList* hsp_list;
    ///        Int4 oid = hsplist_in->hsp_array[index]->context;
    ///        Int4 hspcnt;
    ///        /* Find the first HSP that corresponds to a different subject.
    ///           At the same time, set all HSP contexts to 0, since this is what
    ///           traceback code expects. */
    ///        for (next_index = index; next_index < hsplist_in->hspcnt;
    ///             ++next_index) {
    ///            if (hsplist_in->hsp_array[next_index]->context != oid)
    ///                break;
    ///            hsplist_in->hsp_array[next_index]->context = 0;
    ///        }
    ///        ...
    ///        /* Insert this HSPList into this query's hit list. */
    ///        Blast_HitListUpdate(hit_list, hsp_list);
    ///    }
    /// ```
    /// `s_ScoreCompareHSPWithContext` (hspfilter_collector.c:191-193): the context, then
    /// `ScoreCompareHSPs`.
    fn collect(&self, mut hsps: Vec<XPrelimHsp>) -> BlastpHitList<XPrelimHsp> {
        hsps.sort_by(|a, b| {
            a.hsp
                .query_context
                .cmp(&b.hsp.query_context)
                .then_with(|| score_compare(&a.hsp, &b.hsp))
        });
        let mut hit_list = BlastpHitList::new(self.prelim_hitlist_size);
        let mut start = 0usize;
        while start < hsps.len() {
            let oid = hsps[start].hsp.query_context;
            let end = start
                + hsps[start..]
                    .iter()
                    .take_while(|hsp| hsp.hsp.query_context == oid)
                    .count();
            let list_hsps: Vec<XPrelimHsp> = hsps[start..end]
                .iter()
                .map(|hsp| {
                    let mut hsp = *hsp;
                    hsp.hsp.query_context = 0;
                    hsp.hsp.s_idx = oid as u32;
                    hsp
                })
                .collect();
            hit_list.update(BlastpHspList {
                oid: oid as u32,
                query_index: self.query_index,
                hsps: list_hsps,
                best_evalue: f64::MAX,
            });
            start = end;
        }
        hit_list
    }

    /// `BLAST_GetGappedScore` over the profiles of the concatenation: each initial HSP's profile
    /// by its word's offset, the PSSM moved to the profile's start, the query's cutoff, the
    /// containment tree over the concatenated coordinates. No restricted alignment (blastp's
    /// only) and no chaining (an option of blastp-fast).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3858-4085
    /// ```c
    ///    for (index=0; index<init_hitlist->total; index++)
    ///    {
    ///       ...
    ///       s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
    ///                                         &query_tmp, &context);
    ///       ...
    ///       if (rpsblast_pssms)
    ///          gap_align->sbp->psi_matrix->pssm->data = rpsblast_pssms +
    ///              query_info->contexts[context].query_offset;
    ///       ...
    ///          q_start = init_hsp->ungapped_data->q_start;
    ///          q_end = q_start + init_hsp->ungapped_data->length;
    ///          s_start = init_hsp->ungapped_data->s_start;
    ///          s_end = s_start + init_hsp->ungapped_data->length;
    ///          score = init_hsp->ungapped_data->score;
    ///       ...
    ///       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
    ///                                         hit_options->min_diag_separation))
    ///       {
    ///          ...
    ///          if (is_rpsblast)
    ///             cutoff = hit_params->cutoffs[rps_cutoff_index].cutoff_score;
    ///          ...
    ///             max_offset =
    ///                BlastGetStartForGappedAlignment(query_tmp.sequence,
    ///                   subject->sequence, gap_align->sbp,
    ///                   init_hsp->ungapped_data->q_start,
    ///                   init_hsp->ungapped_data->length,
    ///                   init_hsp->ungapped_data->s_start,
    ///                   init_hsp->ungapped_data->length);
    ///             init_hsp->offsets.qs_offsets.s_off +=
    ///                 max_offset - init_hsp->offsets.qs_offsets.q_off;
    ///             init_hsp->offsets.qs_offsets.q_off = max_offset;
    ///          ...
    ///             status = s_BlastProtGappedAlignment(program_number, &query_tmp,
    ///                                                 subject, gap_align,
    ///                                                 score_params, init_hsp,
    ///                                                 restricted_alignment,
    ///                                                 fence_hit);
    ///          ...
    ///          if (gap_align->score >= cutoff) {
    ///              ...
    ///              status = Blast_HSPInit(gap_align->query_start,
    ///                            gap_align->query_stop, gap_align->subject_start,
    ///                            gap_align->subject_stop,
    ///                            init_hsp->offsets.qs_offsets.q_off,
    ///                            init_hsp->offsets.qs_offsets.s_off, context,
    ///                            query_frame, subject->frame, gap_align->score,
    ///                            &(gap_align->edit_script), &new_hsp);
    ///              ...
    ///              status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
    ///              ...
    ///              status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
    ///                                      eQueryAndSubject);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3824-3832
    /// ```c
    ///    else {
    ///       tree = Blast_IntervalTreeInit(0, query->length+1,
    ///                                     0, subject->length+1);
    ///    }
    /// ```
    fn gapped_score(
        &self,
        init_hsps: &[RpsInitHsp],
        tree: &mut BlastIntervalTree,
        hsps: &mut Vec<BlastpPreliminaryHsp>,
        gap: &mut PssmGapScratch,
        subject_seq: &[u8],
        subject_len: usize,
    ) -> Result<()> {
        let concat_length = self.rows.len() as i32;
        tree.reset_with_bounds(0, concat_length + 1, 0, subject_len as i32 + 1);
        for init in init_hsps {
            let mut init = *init;
            let context = context_of(self.offsets, init.q_off);
            let query_start = self.offsets[context];
            let profile_length = (self.offsets[context + 1] - query_start - 1) as usize;
            // `s_AdjustInitialHSPOffsets`: the offsets of the profile.
            init.q_off -= query_start;
            init.ungapped.q_start -= query_start;
            let ungapped = init.ungapped;
            if ungapped.q_start < 0
                || (ungapped.q_start + ungapped.length) as usize > profile_length
            {
                bail!(
                    "a profile database whose scores extend an ungapped alignment across two profiles (NCBI BLAST+ asserts against it) is not supported by LOSAT's {}",
                    self.program
                );
            }
            let candidate = TreeHsp {
                query_offset: ungapped.q_start,
                query_end: ungapped.q_start + ungapped.length,
                subject_offset: ungapped.s_start,
                subject_end: ungapped.s_start + ungapped.length,
                score: ungapped.score,
                query_frame: 0,
                query_length: profile_length as i32,
                query_context_offset: query_start,
                subject_frame_sign: 0,
            };
            if tree.contains_hsp(&candidate, query_start, 0) {
                continue;
            }
            let profile = ProfileSlice {
                rows: self.rows,
                start: query_start as usize,
                len: profile_length,
            };
            let max_offset = pssm_get_start_for_gapped_alignment(
                &profile,
                subject_seq,
                ungapped.q_start,
                ungapped.length,
                ungapped.s_start,
            );
            init.s_off += max_offset - init.q_off;
            init.q_off = max_offset;
            let aligned = pssm_score_only_gapped_alignment(
                &profile,
                subject_seq,
                subject_len,
                init.q_off as usize,
                init.s_off as usize,
                self.gap_open,
                self.gap_extend,
                self.x_drop_gapped,
                gap,
            );
            if aligned.score < self.hit_cutoff {
                continue;
            }
            let hsp = BlastpPreliminaryHsp {
                query_start: aligned.query_start,
                query_end: aligned.query_stop,
                subject_start: aligned.subject_start,
                subject_end: aligned.subject_stop,
                gapped_query_start: init.q_off,
                gapped_subject_start: init.s_off,
                raw_score: aligned.score,
                query_context: context as i32,
                query_frame: 0,
                subject_frame: 0,
                query_length: profile_length,
                q_idx: self.query_index,
                s_idx: context as u32,
            };
            tree.add_hsp(
                TreeHsp {
                    query_offset: hsp.query_start,
                    query_end: hsp.query_end,
                    subject_offset: hsp.subject_start,
                    subject_end: hsp.subject_end,
                    score: hsp.raw_score,
                    query_frame: 0,
                    query_length: profile_length as i32,
                    query_context_offset: query_start,
                    subject_frame_sign: 0,
                },
                query_start,
                IndexMethod::QueryAndSubject,
            );
            hsps.push(hsp);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;

    fn prelim_list(oid: u32, score: i32, e_value: f64) -> BlastpHspList<XPrelimHsp> {
        BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![XPrelimHsp {
                hsp: BlastpPreliminaryHsp {
                    query_start: 0,
                    query_end: 10,
                    subject_start: 0,
                    subject_end: 10,
                    gapped_query_start: 0,
                    gapped_subject_start: 0,
                    raw_score: score,
                    query_context: 0,
                    query_frame: 0,
                    subject_frame: 0,
                    query_length: 100,
                    q_idx: 0,
                    s_idx: oid,
                },
                e_value,
            }],
            best_evalue: f64::MAX,
        }
    }

    // blast_engine.c:874-899, blast_hits.c (`Blast_HitListUpdate`, `s_EvalueCompareHSPLists`):
    // without the Gumbel block every preliminary HSP has the e-value 0, so the lists are ordered
    // by the first HSP's score (down), then the ordinal id (down), which is the order of the
    // composition-based redo's reads: the last list read is the lowest, and its profile's K is
    // the batch's (`RpsTraceback::trace_query_redo`).
    #[test]
    fn lists_of_e_value_zero_are_ordered_by_score_then_ordinal_id() {
        let mut hit_list = BlastpHitList::new(10);
        for (oid, score) in [(0, 50), (1, 80), (2, 80), (3, 65)] {
            hit_list.update(prelim_list(oid, score, 0.0));
        }
        hit_list.sort_by_evalue();
        let oids: Vec<u32> = hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|list| list.oid)
            .collect();
        assert_eq!(oids, vec![2, 1, 3, 0]);
        // With e-values the best e-value decides first.
        let mut hit_list = BlastpHitList::new(10);
        for (oid, score, e_value) in [(0, 50, 1e-5), (1, 80, 1e-2), (2, 80, 1e-9)] {
            hit_list.update(prelim_list(oid, score, e_value));
        }
        hit_list.sort_by_evalue();
        let oids: Vec<u32> = hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|list| list.oid)
            .collect();
        assert_eq!(oids, vec![2, 0, 1]);
    }

    fn kept_oids(lists: &[(u32, i32, f64)]) -> Vec<u32> {
        let size = crate::algorithm::blastn::hsp::get_prelim_hitlist_size(1, false, true);
        assert_eq!(size, 10);
        let mut hit_list = BlastpHitList::new(size as usize);
        for &(oid, score, e_value) in lists {
            hit_list.update(prelim_list(oid, score, e_value));
        }
        hit_list.sort_by_evalue();
        hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|list| list.oid)
            .collect()
    }

    // The prelim hit list holds `min(max(2*size, 10), size + 50)` lists, 10 for `-max_target_seqs`
    // 1 to 5 (blast_hits.c:68), and a list past that capacity goes through
    // `Blast_HitListUpdate` (blast_hits.c:3243-3297) into the heap whose root is the worst list:
    // ```c
    //       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    //       if (evalue_order < 0) { Blast_HSPListFree(hsp_list); }  /* new is worse: discard */
    //       else { s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list); } /* replaces the root */
    // ```
    // with `s_EvalueCompareHSPLists` (blast_hits.c:3078-3106) ordering by best E-value (up), then
    // the first HSP's score (down), then `BLAST_CMP(h2->oid, h1->oid)`: of two lists equal in
    // E-value and score the one of the higher ordinal id is the better, so the root is the lowest
    // id of the worst tie and a new list of a higher id replaces it ("newer hits with score and
    // e-value both equal to the current worst will be saved"), one of a lower id is discarded.
    // The kept lists are therefore the best 10 of the total order whatever the insertion order,
    // and sorted they read best first.
    #[test]
    fn a_full_prelim_list_of_equal_e_values_keeps_the_best_ten_of_the_total_order() {
        let ascending: Vec<u32> = (0..18).collect();
        let descending: Vec<u32> = (0..18).rev().collect();
        let shuffled = [
            9u32, 3, 17, 0, 12, 5, 16, 1, 8, 14, 2, 11, 6, 15, 4, 13, 7, 10,
        ];
        for order in [ascending.as_slice(), descending.as_slice(), &shuffled] {
            // Equal e-value 0 and equal score: the ten highest ordinal ids, highest first.
            let equal: Vec<(u32, i32, f64)> = order.iter().map(|&oid| (oid, 100, 0.0)).collect();
            assert_eq!(
                kept_oids(&equal),
                vec![17, 16, 15, 14, 13, 12, 11, 10, 9, 8],
                "order {order:?}"
            );
            // Even ids score 90 and odd ids 70: the nine of score 90 (high id first), then
            // the highest id of score 70.
            let mixed: Vec<(u32, i32, f64)> = order
                .iter()
                .map(|&oid| (oid, if oid % 2 == 0 { 90 } else { 70 }, 0.0))
                .collect();
            assert_eq!(
                kept_oids(&mixed),
                vec![16, 14, 12, 10, 8, 6, 4, 2, 0, 17],
                "order {order:?}"
            );
            // The E-value decides before the score and the id: the ten lists of E-value 0 and
            // score 10 (ids 8 to 17) beat the eight of E-value 1e-30 and score 100 (ids 0 to 7).
            let with_e: Vec<(u32, i32, f64)> = order
                .iter()
                .map(|&oid| {
                    if oid < 8 {
                        (oid, 100, 1e-30)
                    } else {
                        (oid, 10, 0.0)
                    }
                })
                .collect();
            let expected = vec![17, 16, 15, 14, 13, 12, 11, 10, 9, 8];
            assert_eq!(kept_oids(&with_e), expected, "order {order:?}");
        }
    }

    #[test]
    fn the_context_of_an_offset_is_the_profile_that_holds_it() {
        let offsets = [0, 399, 1130, 3482, 3772, 4621, 5035];
        assert_eq!(context_of(&offsets, 0), 0);
        assert_eq!(context_of(&offsets, 398), 0);
        assert_eq!(context_of(&offsets, 399), 1);
        assert_eq!(context_of(&offsets, 5033), 5);
        assert_eq!(context_of(&offsets, 5034), 5);
    }
}
