//! The preliminary stage of one user query against one profile volume
//! (`s_BlastSearchEngineCore` with the concatenated profiles as the query): the word finder,
//! `BLAST_GetGappedScore` over the profiles (one context each), the purge of common endpoints,
//! the sort by score, the e-values and the reap by the preliminary e-value, then the collector's
//! split into one list per profile (`s_BlastHSPCollectorRun_RPS`) and the query's hit list.

use anyhow::{bail, Result};

use super::wordfinder::{word_finder, RpsInitHsp, RpsWordParams, RpsWordScratch};
use super::ProfileSlice;
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, TreeHsp};
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
    /// The unit's concatenated length `Lc` (`RpsRun::concat_length`, `query->length` of
    /// `BLAST_GetGappedScore`): the interval tree covers `[0, Lc + 1]`.
    pub(crate) concat_length: i32,
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
    /// A query's HSPs of its contexts, concatenated for `finish` (kept from query to query).
    pub(crate) combined: Vec<BlastpPreliminaryHsp>,
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-300
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
    /// query), sorted by score, at the end of `out`, as `s_BlastSearchEngineOneContext` gives them
    /// to `Blast_HSPListAppend`; the query's contexts' HSPs are concatenated in the context order
    /// for `finish`.
    pub(crate) fn search_frame(
        &self,
        scratch: &mut RpsPrelimScratch,
        subject_seq: &[u8],
        frame: i8,
        out: &mut Vec<BlastpPreliminaryHsp>,
    ) -> Result<()> {
        let subject_len = subject_seq.len() - 1;
        let subject = &subject_seq[..subject_len];
        let RpsPrelimScratch {
            word,
            init_hsps,
            tree,
            hsps,
            gap,
            combined: _,
        } = scratch;
        init_hsps.clear();
        word_finder(&self.word, word, subject, init_hsps);
        // An extension that reads past the rows what NCBI cannot (`RpsRows::past`; only where
        // `Lc` is above them).
        if let Some(past) = self.word.rows.past {
            if past
                .overran
                .swap(false, std::sync::atomic::Ordering::Relaxed)
            {
                bail!(
                    "an RPS-BLAST volume searched as {} concatenated positions, above the {} rows of its profiles (an alias file's LENGTH and NSEQ over one volume, or a volume reached twice), whose ungapped extension of a query runs past the row after the last profile, or reads that row past the .rps file's last memory page (NCBI BLAST+ 2.17.0 then reads past its score matrix), is not supported by LOSAT's {}",
                    self.concat_length,
                    self.rows.len(),
                    self.program
                );
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:497
        if init_hsps.is_empty() {
            return Ok(());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:282
        sort_init_hsps_by_score(init_hsps);
        hsps.clear();
        self.gapped_score(init_hsps, tree, hsps, gap, subject_seq, subject_len)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:541-555
        purge_hsps_with_common_endpoints(hsps, |hsp| hsp);
        if hsps
            .windows(2)
            .any(|pair| score_compare(&pair[0], &pair[1]) == std::cmp::Ordering::Greater)
        {
            hsps.sort_by(score_compare);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
        // Without a cap (`hsp_num_max` 0: `INT4_MAX`, blast_hits.c:213-226) the lists of the
        // frames are concatenated in the frame order (`s_BlastHSPListsCombineByScore`).
        out.extend(hsps.iter().map(|hsp| BlastpPreliminaryHsp {
            subject_frame: frame.into(),
            ..*hsp
        }));
        Ok(())
    }

    /// The query's preliminary hit list on this volume from the HSPs of its contexts `hsps`
    /// (`search_frame`'s, in the context order; emptied here), with the e-values of the length
    /// `stat_length` (`RpsQuery::stat_length`), or `None` when it has no HSP. The lists' oids are
    /// the profiles' local ordinal ids.
    pub(crate) fn finish(
        &self,
        hsps: &mut Vec<BlastpPreliminaryHsp>,
        stat_length: usize,
    ) -> Option<BlastpHitList<XPrelimHsp>> {
        if hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:881-899
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1876-1884
        // `subject_length` is the user query's length (`stat_length`), the context's length the
        // profile's.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:874-880
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-670
            // NCBI drops an HSP whose e-value is above `prelim_evalue`: a NaN e-value stays, and a
            // `prelim_evalue` of NaN (an `-evalue` of NaN, which NCBI's option check lets pass)
            // keeps every HSP.
            if crate::common::evalue_not_above(e_value, self.prelim_evalue) {
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850-4087
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3824-3832
    fn gapped_score(
        &self,
        init_hsps: &[RpsInitHsp],
        tree: &mut BlastIntervalTree,
        hsps: &mut Vec<BlastpPreliminaryHsp>,
        gap: &mut PssmGapScratch,
        subject_seq: &[u8],
        subject_len: usize,
    ) -> Result<()> {
        // `Lc + 1` is at most 2^30 + 1 (`RpsUnit::check_engine_length`); an `Lc` below the
        // rows leaves the regions of later profiles past the tree's end.
        tree.reset_with_bounds(0, self.concat_length + 1, 0, subject_len as i32 + 1);
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603-657,703-712
            // and :102-109 (s_IntervalNodeInit)
            // The tree's root is `[0, R]`, where `R` is one more than `Lc`. A region (the HSP's
            // profile offset plus its query range) that starts past `R` goes right at every node
            // `[l, R]` of the right spine (its middle is at most `R`) and becomes the leaf at the
            // spine's end; a second one in the same tree (a query's, a frame's for RPSTBLASTN)
            // meets that leaf, and the split makes `[R + 1, R]`, whose middle is `R` again: both go
            // right, and `BlastIntervalTreeAddHSP` splits that node without end while its node pool
            // doubles (memory runs out). A region that starts at or before `R` and ends past it
            // straddles a middle of the spine and ends in a subject tree. Below `R = 0` (an `Lc` of
            // -2 or less) every middle of the spine is below the regions' starts. `try_add_hsp`
            // finds the insertion that makes such a node, after NCBI's end-point checks (which may
            // drop the HSP or the leaf it would meet): only an `Lc` below the rows (an alias file's
            // `LENGTH` and `NSEQ` over one volume) leaves regions past the end, and the search
            // stops there (LOSAT's refusal of the query batch; at most one region past the end in
            // a tree is NCBI's search).
            if !tree.try_add_hsp(
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
            ) {
                bail!(
                    "an RPS-BLAST volume searched as {} concatenated positions, below the {} rows of its profiles (an alias file's LENGTH and NSEQ, or NCBI's totals of the OIDs its OIDLIST keeps, over one volume, or profiles whose .pin lengths are below their .rps rows), where two gapped alignments of a query start past position {} (NCBI BLAST+ 2.17.0's interval tree then splits a node without end), is not supported by LOSAT's {}",
                    self.concat_length,
                    self.rows.len(),
                    i64::from(self.concat_length) + 1,
                    self.program
                );
            }
            hsps.push(hsp);
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "prelim_tests.rs"]
mod tests;
