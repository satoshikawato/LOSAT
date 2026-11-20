//! The preliminary search of one subject with a PSSM: the subject scan of blastp's lookup table
//! (built from the PSSM, `lookup.rs`), the PSSM ungapped extensions (`ungapped.rs`), and NCBI's
//! `BLAST_GetGappedScore` for the one context of a PSSM query, without restricted alignment and
//! without chaining, then the purge of common endpoints, the sort by score, the e-values of
//! `kbp_gap_psi` and the reap by the preliminary e-value (`s_BlastSearchEngineCore`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:490-555
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:445
//!
//! The PSSM query has one context at query offset 0 (`s_AdjustHspOffsetsAndGetQueryData`
//! subtracts nothing), so `query_info->contexts[0]` stands for every context lookup of the
//! reference code. The chaining of `BLAST_GetGappedScore` is blastp-fast's only (psiblast's task
//! is blastp's), and `hit_params->low_score` is set by no psiblast option.

use super::evalue::EvalueStats;
use super::gapalign::{
    pssm_get_start_for_gapped_alignment, pssm_score_only_gapped_alignment, PssmGapScratch,
};
use super::rows::PssmRows;
use super::ungapped::{pssm_extend_one_hit, pssm_extend_two_hit};
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::extension::BlastpUngappedData;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{
    link_hsp_list, reap_hsplist_by_evalue, update_best_evalue, BlastpHspList,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::tblastx::blast_aascan::{
    s_blast_aa_scan_subject_one_range, BlastOffsetPair as OffsetPair,
};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::core::link_hsps::{EvenGapLinkParams, LinkContext};
use crate::stats::KarlinParams;

/// NCBI's `BlastInitHSP` with its `BlastUngappedData`.
#[derive(Clone, Copy, Debug)]
struct PssmInitHsp {
    q_off: i32,
    s_off: i32,
    ungapped: BlastpUngappedData,
}

/// What the preliminary search of every subject reads: the PSSM, its lookup table, and the
/// search's cutoffs (computed with `kbp_psi`, `kbp_std` and `kbp_gap_psi` as NCBI does).
pub(crate) struct PssmPrelim<'a, R: PssmRows> {
    pub(crate) rows: &'a R,
    pub(crate) lookup: &'a BlastAaLookupTable,
    pub(crate) word_size: i32,
    /// `-window_size`; 0 for the one-hit search.
    pub(crate) window: i32,
    /// `cutoffs->x_dropoff` of the context (`ceil(x_dropoff ln 2 / kbp_psi->Lambda)` of the
    /// search's `-xdrop_ungap`, 7 bits by default).
    pub(crate) x_dropoff: i32,
    /// `word_params->cutoffs[0].cutoff_score`.
    pub(crate) word_cutoff: i32,
    /// `hit_params->cutoffs[0].cutoff_score`.
    pub(crate) hit_cutoff: i32,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// `ext_params->gap_x_dropoff`.
    pub(crate) x_drop_gapped: i32,
    /// `kbp_gap_psi[0]`, for the preliminary e-values.
    pub(crate) karlin: &'a KarlinParams,
    /// `sbp->gbp` and the context's effective search space, for the preliminary e-values.
    pub(crate) stats: EvalueStats<'a>,
    /// `hit_params->prelim_evalue`.
    pub(crate) prelim_evalue: f64,
    /// psiblast's `-sum_stats true`: the preliminary linking parameters (`hit_params->
    /// link_hsp_params`) and the query's context as the linking reads it; `None` without sum
    /// statistics.
    pub(crate) link: Option<(EvenGapLinkParams, LinkContext)>,
    /// `-subject_besthit`: the subject's list is filtered after its sort by score.
    pub(crate) subject_besthit: bool,
}

/// The memory of one subject's preliminary search, kept from subject to subject as NCBI keeps
/// its `BlastInitHitList`, offset pairs, interval tree and `BlastGapAlignStruct`.
pub(crate) struct PssmSubjectScratch {
    init_hsps: Vec<PssmInitHsp>,
    tree: BlastIntervalTree,
    hsps: Vec<BlastpPreliminaryHsp>,
    gap: PssmGapScratch,
}

impl PssmSubjectScratch {
    pub(crate) fn new() -> Self {
        Self {
            init_hsps: Vec::new(),
            tree: BlastIntervalTree::new(0, 1, 0, 1),
            hsps: Vec::new(),
            gap: PssmGapScratch::new(),
        }
    }
}

impl Default for PssmSubjectScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// `BlastSaveInitHsp`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:358-372
#[inline]
fn save_init_hsp(
    init_hsps: &mut Vec<PssmInitHsp>,
    ungapped: BlastpUngappedData,
    q_off: i32,
    s_off: i32,
) {
    init_hsps.push(PssmInitHsp {
        q_off,
        s_off,
        ungapped,
    });
}

/// `score_compare_match` (the sort of `Blast_InitHitListSortByScore`): score down, subject start
/// up, length down, query start up. glibc's qsort under NCBI BLAST+ is a stable merge sort, which
/// `sort_by` is.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-309
fn sort_init_hsps_by_score(init_hsps: &mut [PssmInitHsp]) {
    init_hsps.sort_by(|a, b| {
        b.ungapped
            .score
            .cmp(&a.ungapped.score)
            .then_with(|| a.ungapped.s_start.cmp(&b.ungapped.s_start))
            .then_with(|| b.ungapped.length.cmp(&a.ungapped.length))
            .then_with(|| a.ungapped.q_start.cmp(&b.ungapped.q_start))
    });
}

/// The interval-tree node of a preliminary HSP of the context at offset 0.
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
        subject_frame_sign: 0,
    }
}

/// `ScoreCompareHSPs`: score down, subject offset up, subject end down, query offset up, query
/// end down.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1350
pub(crate) fn score_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

/// `Blast_HSPListPurgeHSPsWithCommonEndpoints` with `purge` (one context, subject frame 0): of
/// the HSPs that share a start, then of those that share an end, the first in NCBI's order stays.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2527
/// `s_QueryOffsetCompareHSPs` (blast_hits.c:2268-2315): context, query offset, subject offset up,
/// score, query end, subject end down; `s_QueryEndCompareHSPs` (2333-2379): context, query end,
/// subject end up, score, query offset, subject offset down.
pub(crate) fn purge_hsps_with_common_endpoints<T>(
    hsps: &mut Vec<T>,
    hsp: impl Fn(&T) -> &BlastpPreliminaryHsp,
) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        let (a, b) = (hsp(a), hsp(b));
        a.query_context
            .cmp(&b.query_context)
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| a.subject_start.cmp(&b.subject_start))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| b.query_end.cmp(&a.query_end))
            .then_with(|| b.subject_end.cmp(&a.subject_end))
    });
    hsps.dedup_by(|later, kept| {
        let (later, kept) = (hsp(later), hsp(kept));
        kept.query_context == later.query_context
            && kept.query_start == later.query_start
            && kept.subject_start == later.subject_start
            && kept.subject_frame == later.subject_frame
    });
    hsps.sort_by(|a, b| {
        let (a, b) = (hsp(a), hsp(b));
        a.query_context
            .cmp(&b.query_context)
            .then_with(|| a.query_end.cmp(&b.query_end))
            .then_with(|| a.subject_end.cmp(&b.subject_end))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| b.query_start.cmp(&a.query_start))
            .then_with(|| b.subject_start.cmp(&a.subject_start))
    });
    hsps.dedup_by(|later, kept| {
        let (later, kept) = (hsp(later), hsp(kept));
        kept.query_context == later.query_context
            && kept.query_end == later.query_end
            && kept.subject_end == later.subject_end
            && kept.subject_frame == later.subject_frame
    });
}

impl<R: PssmRows> PssmPrelim<'_, R> {
    /// The preliminary HSP list of subject `s_idx` (`aa_seq`: its letters between two sentinels,
    /// `aa_len` of them), or `None` when it has none. `diag_array`, `diag_mask` and
    /// `diag_offset` are the diagonal table of the subjects searched before it, which the caller
    /// moves on after the subject (`BlastExtendWordExit`); `offset_pairs` is the scan's buffer
    /// (`GetOffsetArraySize` pairs) of the one-hit search.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn search_subject(
        &self,
        scratch: &mut PssmSubjectScratch,
        offset_pairs: &mut [OffsetPair],
        diag_array: &mut [DiagStruct],
        diag_mask: i32,
        diag_offset: i32,
        s_idx: u32,
        aa_seq: &[u8],
        aa_len: usize,
    ) -> Option<BlastpHspList<XPrelimHsp>> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
        let subject_sequence = &aa_seq[1..];
        let subject_raw = &aa_seq[1..aa_seq.len() - 1];
        let PssmSubjectScratch {
            init_hsps,
            tree,
            hsps,
            gap,
        } = scratch;
        init_hsps.clear();
        self.word_finder(
            init_hsps,
            offset_pairs,
            diag_array,
            diag_mask,
            diag_offset,
            subject_sequence,
            subject_raw,
            aa_len,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:497,500
        if init_hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:228-231
        sort_init_hsps_by_score(init_hsps);
        hsps.clear();
        self.gapped_score(init_hsps, tree, hsps, gap, s_idx, subject_sequence, aa_len);
        purge_hsps_with_common_endpoints(hsps, |hsp| hsp);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1373-1381
        if hsps
            .windows(2)
            .any(|pair| score_compare(&pair[0], &pair[1]) == std::cmp::Ordering::Greater)
        {
            hsps.sort_by(score_compare);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:586-591
        // `-subject_besthit` on the subject's list (one context), sorted by score, before the
        // e-values and the reap by the preliminary e-value.
        if self.subject_besthit {
            crate::core::blast_hits::subject_besthit(hsps, false, |hsp| {
                crate::core::blast_hits::BestHitExtent {
                    context: 0,
                    query_frame: 0,
                    offset: i64::from(hsp.query_start),
                    end: i64::from(hsp.query_end),
                    query_length: 0,
                }
            });
        }
        if hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1522-1539
        // `sbp->kbp_gap` is `kbp_gap_psi` (`Blast_QueryIsPssm`); Spouge's finite-size correction
        // reads the query's and the subject's lengths, the e-value without the Gumbel block the
        // context's search space (`EvalueStats::e_value`, blast_hits.c:1872-1890). psiblast
        // computes them without the block too (only RPS-BLAST skips them then).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:874-880
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
        // With psiblast's `-sum_stats true` the HSPs take the e-values of their linked sets
        // (`kbp_gap_psi`), in the list the linking sorts by score.
        let query_length = self.rows.len() as i32;
        let mut list = BlastpHspList {
            oid: s_idx,
            query_index: 0,
            hsps: hsps
                .iter()
                .map(|hsp| XPrelimHsp {
                    hsp: *hsp,
                    e_value: match self.link {
                        Some(_) => f64::MAX,
                        None => self.stats.e_value(
                            hsp.raw_score,
                            self.karlin,
                            query_length,
                            aa_len as i32,
                        ),
                    },
                })
                .collect(),
            best_evalue: f64::MAX,
        };
        match &self.link {
            Some((params, context)) => link_hsp_list(&mut list, *context, aa_len as i32, params),
            None => update_best_evalue(&mut list),
        }
        reap_hsplist_by_evalue(&mut list, self.prelim_evalue);
        (!list.hsps.is_empty()).then_some(list)
    }

    /// `s_BlastAaWordFinder_TwoHit` (with a window) or `s_BlastAaWordFinder_OneHit` with
    /// `use_pssm`.
    #[allow(clippy::too_many_arguments)]
    fn word_finder(
        &self,
        init_hsps: &mut Vec<PssmInitHsp>,
        offset_pairs: &mut [OffsetPair],
        diag_array: &mut [DiagStruct],
        diag_mask: i32,
        diag_offset: i32,
        subject_sequence: &[u8],
        subject_raw: &[u8],
        aa_len: usize,
    ) {
        let rows = self.rows;
        let word_size = self.word_size;
        let x_dropoff = self.x_dropoff;
        let cutoff = self.word_cutoff;
        if self.window != 0 {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:551-603
            // The scan and the diagonal test are blastp's (`x_scan_two_hit`, whose loop gives
            // the reference's outcomes); this closure is the part from the context check on.
            let candidate = |query_offset: i32,
                             subject_offset: i32,
                             last_hit: i32,
                             diff: i32,
                             diag_entry: &mut DiagStruct| {
                if query_offset - diff < 0 {
                    diag_entry.set_last_hit(subject_offset + diag_offset);
                    return;
                }
                let hit = pssm_extend_two_hit(
                    rows,
                    subject_raw,
                    last_hit + word_size,
                    subject_offset,
                    query_offset,
                    x_dropoff,
                    word_size,
                );
                if hit.ungapped_data.score >= cutoff {
                    save_init_hsp(init_hsps, hit.ungapped_data, query_offset, subject_offset);
                }
                if hit.right_extend {
                    diag_entry.set_flag(1);
                    diag_entry.set_last_hit(hit.s_last_off - (word_size - 1) + diag_offset);
                } else {
                    diag_entry.set_last_hit(subject_offset + diag_offset);
                }
            };
            crate::algorithm::blastp::x_seed_fused::x_scan_two_hit(
                self.lookup,
                subject_sequence,
                aa_len,
                diag_array,
                diag_mask,
                diag_offset,
                self.window,
                word_size,
                candidate,
            );
            return;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:681-802
        let array_size = offset_pairs.len() as i32;
        let mut scan_range = [0i32, 0, aa_len as i32 - word_size];
        while scan_range[1] <= scan_range[2] {
            let hits = s_blast_aa_scan_subject_one_range(
                self.lookup,
                subject_sequence,
                offset_pairs,
                array_size,
                &mut scan_range,
            );
            for pair in &offset_pairs[..hits as usize] {
                let query_offset = pair.q_off as i32;
                let subject_offset = pair.s_off as i32;
                let diag_coord = ((subject_offset - query_offset) & diag_mask) as usize;
                let diag_entry = &mut diag_array[diag_coord];
                let diff = subject_offset - (diag_entry.last_hit() - diag_offset);
                if diff < 0 {
                    continue;
                }
                let hit = pssm_extend_one_hit(
                    rows,
                    subject_raw,
                    subject_offset,
                    query_offset,
                    x_dropoff,
                    word_size,
                );
                if hit.ungapped_data.score >= cutoff {
                    save_init_hsp(init_hsps, hit.ungapped_data, query_offset, subject_offset);
                }
                diag_entry.set_last_hit(hit.s_last_off - (word_size - 1) + diag_offset);
            }
        }
    }

    /// `BLAST_GetGappedScore` for the one context of a PSSM query, without chaining and without
    /// restricted alignment (`hit_params->restricted_align` is blastp's,
    /// `blast_parameters.c:851-855`), so there is no recomputation of a query's alignments.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850-4092
    #[allow(clippy::too_many_arguments)]
    fn gapped_score(
        &self,
        init_hsps: &[PssmInitHsp],
        tree: &mut BlastIntervalTree,
        hsps: &mut Vec<BlastpPreliminaryHsp>,
        gap: &mut PssmGapScratch,
        s_idx: u32,
        subject_sequence: &[u8],
        aa_len: usize,
    ) {
        let query_length = self.rows.len();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3826-3832
        tree.reset_with_bounds(0, query_length as i32 + 1, 0, aa_len as i32 + 1);
        for init in init_hsps {
            let mut init = *init;
            let ungapped = init.ungapped;
            let candidate = TreeHsp {
                query_offset: ungapped.q_start,
                query_end: ungapped.q_start + ungapped.length,
                subject_offset: ungapped.s_start,
                subject_end: ungapped.s_start + ungapped.length,
                score: ungapped.score,
                query_frame: 1,
                query_length: query_length as i32,
                query_context_offset: 0,
                subject_frame_sign: 0,
            };
            if tree.contains_hsp(&candidate, 0, 0) {
                continue;
            }
            let max_offset = pssm_get_start_for_gapped_alignment(
                self.rows,
                subject_sequence,
                ungapped.q_start,
                ungapped.length,
                ungapped.s_start,
            );
            init.s_off += max_offset - init.q_off;
            init.q_off = max_offset;
            let aligned = pssm_score_only_gapped_alignment(
                self.rows,
                subject_sequence,
                aa_len,
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
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length,
                q_idx: 0,
                s_idx,
            };
            tree.add_hsp(tree_hsp(&hsp), 0, IndexMethod::QueryAndSubject);
            hsps.push(hsp);
        }
    }
}
