//! The preliminary search of one subject with a PSSM: the subject scan of blastp's lookup table
//! (built from the PSSM, `lookup.rs`), the PSSM ungapped extensions (`ungapped.rs`), and NCBI's
//! `BLAST_GetGappedScore` for the one context of a PSSM query, without restricted alignment and
//! without chaining, then the purge of common endpoints, the sort by score, the e-values of
//! `kbp_gap_psi` and the reap by the preliminary e-value (`s_BlastSearchEngineCore`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:490-547
//! ```c
//!         if (aux_struct->WordFinder) {
//!             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
//!                                    word_params, aux_struct->ewp,
//!                                    aux_struct->offset_pairs,
//!                                    kScanSubjectOffsetArraySize,
//!                                    init_hitlist, ungapped_stats);
//!
//!             if (init_hitlist->total == 0) continue;
//!         }
//!         ...
//!         if (aux_struct->GetGappedScore) {
//!             status = aux_struct->GetGappedScore(program_number, query,
//!                     query_info,
//!                     subject, gap_align, score_params, ext_params, hit_params,
//!                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
//!         }
//!         ...
//!         /* No need to do this for short reads */
//!         if (aux_struct->GetGappedScore) {
//!
//!             /* Removes redundant HSPs. */
//!             Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
//!         ...
//!         }
//!
//!         Blast_HSPListSortByScore(hsp_list);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:445
//! ```c
//!     Int4** matrix = (gap_align->positionBased) ? gap_align->sbp->psi_matrix->pssm->data : gap_align->sbp->matrix->data;
//! ```
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
/// ```c
///     ungapped_data->q_start = q_start;
///     ungapped_data->s_start = s_start;
///     ungapped_data->length = len;
///     ungapped_data->score = score;
///
///     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
/// ```
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
///    ...
///     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
///           sizeof(BlastInitHSP), score_compare_match);
/// ```
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
/// ```c
///     if (0 == (result = BLAST_CMP(h2->score, h1->score)) &&
///         0 == (result = BLAST_CMP(h1->subject.offset, h2->subject.offset)) &&
///         0 == (result = BLAST_CMP(h2->subject.end, h1->subject.end)) &&
///         0 == (result = BLAST_CMP(h1->query.offset, h2->query.offset))) {
///         result = BLAST_CMP(h2->query.end, h1->query.end);
///     }
/// ```
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
///          hsp_count--;
///          hsp = hsp_array[i+j];
///          ...
///    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
///    i = 0;
///    while (i < hsp_count) {
///       j = 1;
///       while (i+j < hsp_count &&
///              hsp_array[i] && hsp_array[i+j] &&
///              hsp_array[i]->context == hsp_array[i+j]->context &&
///              hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
///              hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
///              hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
/// ```
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-121
        // ```c
        //     (*seq_blk)->sequence_start = (Uint1 *) buffer;
        //     /* The first byte is a sentinel byte. */
        //     (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
        //     (*seq_blk)->length = length;
        // ```
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:497
        // ```c
        //             if (init_hitlist->total == 0) continue;
        // ```
        if init_hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:228-231
        // ```c
        //     Blast_InitHitListSortByScore(init_hitlist);
        // ```
        sort_init_hsps_by_score(init_hsps);
        hsps.clear();
        self.gapped_score(init_hsps, tree, hsps, gap, s_idx, subject_sequence, aa_len);
        purge_hsps_with_common_endpoints(hsps, |hsp| hsp);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1373-1380
        // ```c
        //     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
        //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        //               ScoreCompareHSPs);
        //     }
        // ```
        if hsps
            .windows(2)
            .any(|pair| score_compare(&pair[0], &pair[1]) == std::cmp::Ordering::Greater)
        {
            hsps.sort_by(score_compare);
        }
        if hsps.is_empty() {
            return None;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1523-1539
        // ```c
        //             Blast_HSPListGetEvalues(program_number, query_info,
        //                                     stat_length, hsp_list,
        //                                     gapped_calculation, FALSE,
        //                                     sbp, 0, 1.0);
        //             status = s_Blast_HSPListReapByPrelimEvalue(hsp_list, hit_params);
        // ```
        // `sbp->kbp_gap` is `kbp_gap_psi` (`Blast_QueryIsPssm`); Spouge's finite-size correction
        // reads the query's and the subject's lengths, the e-value without the Gumbel block the
        // context's search space (`EvalueStats::e_value`, blast_hits.c:1872-1890). psiblast
        // computes them without the block too (only RPS-BLAST skips them then).
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
        // ```c
        //     if (hit_params->link_hsp_params) {
        //         status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
        //                   subject->length, gap_align->sbp, hit_params->link_hsp_params,
        //                   score_options->gapped_calculation);
        // ```
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:551-596
            // ```c
            //                 curr_context = BSearchContextInfo(query_offset, query_info);
            //
            //                 /* Check if the last hit hits current query. Because last_hit
            //                    is never reset, it may contain a hit to the previous
            //                    concatenated query */
            //
            //                 if (query_offset - diff <
            //                     query_info->contexts[curr_context].query_offset) {
            //
            //                     /* there was no last hit for this diagnol; start a new hit */
            //                     diag_array[diag_coord].last_hit =
            //                         subject_offset + diag_offset;
            //                     continue;
            //                 }
            //
            //                 cutoffs = word_params->cutoffs + curr_context;
            //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
            //                                               last_hit + wordsize,
            //                                               subject_offset, query_offset,
            //                                               cutoffs->x_dropoff,
            //                                               &hsp_q, &hsp_s,
            //                                               &hsp_len, use_pssm,
            //                                               wordsize, &right_extend,
            //                                               &s_last_off);
            //
            //                 ++hits_extended;
            //
            //                 /* if the hsp meets the score threshold, report it */
            //                 if (score >= cutoffs->cutoff_score)
            //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
            //                                      query_offset, subject_offset, hsp_len,
            //                                      score);
            //
            //                 /* If an extension to the right happened, reset the last hit
            //                    so that future hits to this diagonal must start over. */
            //
            //                 if (right_extend) {
            //                     diag_array[diag_coord].flag = 1;
            //                     diag_array[diag_coord].last_hit =
            //                         s_last_off - (wordsize - 1) + diag_offset;
            //                 }
            //
            //                 /* Otherwise, make the present hit into the previous hit for
            //                    this diagonal */
            //
            //                 else {
            //                     diag_array[diag_coord].last_hit =
            //                         subject_offset + diag_offset;
            //                 }
            // ```
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:763-802
        // ```c
        //     scan_range[0] = 0;
        //     scan_range[1] = subject->seq_ranges[0].left;
        //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
        //
        //     while (scan_range[1] <= scan_range[2]) {
        //         hits = scansub(lookup_wrap, subject,
        //                                  offset_pairs, array_size, scan_range);
        //         ...
        //         for (i = 0; i < hits; ++i) {
        //             Int4 query_offset = offset_pairs[i].qs_offsets.q_off;
        //             Int4 subject_offset = offset_pairs[i].qs_offsets.s_off;
        //             diag_coord = (subject_offset - query_offset) & diag_mask;
        //             diff = subject_offset -
        //                 (diag_array[diag_coord].last_hit - diag_offset);
        //
        //             /* do an extension, but only if we have not already extended
        //                this far */
        //             if (diff >= 0) {
        //                 Int4 curr_context = BSearchContextInfo(query_offset,
        //                                                        query_info);
        //                 BlastUngappedCutoffs *cutoffs = word_params->cutoffs +
        //                                                 curr_context;
        //                 score = s_BlastAaExtendOneHit(matrix, subject, query,
        //                                               subject_offset, query_offset,
        //                                               cutoffs->x_dropoff,
        //                                               &hsp_q, &hsp_s, &hsp_len,
        //                                               wordsize, use_pssm, &s_last_off);
        //                 ++hits_extended;
        //
        //                 if (score >= cutoffs->cutoff_score) {
        //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
        //                                      query_offset, subject_offset, hsp_len,
        //                                      score);
        //                 }
        //                 diag_array[diag_coord].last_hit =
        //                     s_last_off - (wordsize - 1) + diag_offset;
        //             }
        //         }
        //     }
        // ```
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3858-4092
    /// ```c
    ///    for (index=0; index<init_hitlist->total; index++)
    ///    {
    ///       ...
    ///       tmp_init_hsp = init_hsp_array[index];
    ///       ...
    ///       s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
    ///                                         &query_tmp, &context);
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
    ///             cutoff = hit_params->cutoffs[context].cutoff_score;
    ///          ...
    ///          if(is_prot && !score_params->options->is_ooframe) {
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
    ///          }
    ///          if (is_prot) {
    ///             status = s_BlastProtGappedAlignment(program_number, &query_tmp,
    ///                                                 subject, gap_align,
    ///                                                 score_params, init_hsp,
    ///                                                 restricted_alignment,
    ///                                                 fence_hit);
    ///          ...
    ///          if (gap_align->score >= cutoff) {
    ///             ...
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
        // ```c
        //       tree = Blast_IntervalTreeInit(0, query->length+1,
        //                                     0, subject->length+1);
        // ```
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
