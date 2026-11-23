//! Internal TBLASTN Kappa result hitlist order and replacement boundary.

use std::cmp::Ordering;

use anyhow::{ensure, Result};

use crate::common::GapEditOp;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;

use super::stage_d_linking::{
    compare_preliminary_lists_for_kappa, score_compare, LinkedHsp, LinkedHspList,
};

// `-culling_limit`: NCBI's culling writer and pipe over these hit lists (a child module for
// the private fields of `KappaResultHitList`).
#[path = "culling.rs"]
mod culling;
pub(super) use culling::PrelimCulling;

// `-best_hit_overhang` / `-best_hit_score_edge`: NCBI's best hit writer and pipe over these hit
// lists (a child module for the private fields of `KappaResultHitList`).
#[path = "besthit.rs"]
mod besthit;
pub(super) use besthit::PrelimBestHit;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:152-181;
// core/blast_kappa.c:2494-2515:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2510-2510
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-358,
// 3690-3706; core/blast_hits.h:126-143:
// Initializes `num_ident` to zero and fills it after post-redo evaluation.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct KappaHspPayload {
    pub bit_score: f64,
    pub num_ident: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:515-526;
    // c++/src/algo/blast/core/blast_hits.c:767-811,966-989:
    // Alignment length and gap counts follow the same gap_info operation lengths.
    pub num_positives: usize,
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:971-1021,1090-1095:
    // Matching query and subject display characters (`m_QuerySeq` and `m_SubjectSeq`)
    // increment both `num_ident` and `num_positives`.
    // Report-time counts use sequence strings, after the core HSP counts above.
    pub report_num_ident: usize,
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:478
    // The columns less the identities and the gap letters: fewer than none where a PSSM query's
    // gap residue faces a gap, a column that is an identity and a gap
    // (`stage_e_report::recount_shown_counts`).
    pub report_mismatches: i64,
    pub align_length: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub gap_letters: usize,
    pub edit_script: Vec<GapEditOp>,
    pub matrix_adjust_rule: EMatrixAdjustRule,
}

#[derive(Clone, Debug)]
pub(super) struct KappaResultList {
    pub oid: i32,
    pub hsps: LinkedHspList,
    // One payload per HSP, in the same order as `hsps` after each sort.
    // Empty only in the older numeric-only comparison fixtures.
    pub payloads: Vec<KappaHspPayload>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3125-3135,
// 3243-3297,3420-3437:
#[allow(dead_code)] // The public TBLASTN result path remains gated until C/D/E pass.
pub(super) struct KappaResultHitList {
    lists: Vec<KappaResultList>,
    max: usize,
    worst_evalue: f64,
    low_score: i32,
    heapified: bool,
}

#[allow(dead_code)] // The public TBLASTN result path remains gated until C/D/E pass.
impl KappaResultHitList {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3125-3135:
    pub(super) fn new(hitlist_size: usize) -> Result<Self> {
        ensure!(hitlist_size > 0, "TBLASTN hitlist size must be positive");
        Ok(Self {
            lists: Vec::new(),
            max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3297:
    pub(super) fn update(
        &mut self,
        mut incoming: KappaResultList,
    ) -> Result<Option<KappaResultList>> {
        ensure!(
            !incoming.hsps.hsps.is_empty(),
            "Kappa result list must contain an HSP"
        );
        ensure!(
            incoming.payloads.is_empty() || incoming.payloads.len() == incoming.hsps.hsps.len(),
            "Kappa result payload count must match HSP count"
        );
        refresh_best_evalue(&mut incoming.hsps);
        if self.lists.len() < self.max {
            self.worst_evalue = self.worst_evalue.max(incoming.hsps.best_evalue);
            self.low_score = self.low_score.min(incoming.hsps.hsps[0].hsp.score);
            self.lists.push(incoming);
            return Ok(None);
        }
        if !self.heapified {
            for list in &mut self.lists {
                sort_hsps_by_evalue(list);
                refresh_best_evalue(&mut list.hsps);
            }
            for start in (0..self.lists.len() / 2).rev() {
                heapify_down(&mut self.lists, start);
            }
            self.heapified = true;
        }
        sort_hsps_by_evalue(&mut incoming);
        refresh_best_evalue(&mut incoming.hsps);
        let discarded = if compare_result_lists(&self.lists[0], &incoming) == Ordering::Less {
            incoming
        } else {
            let old = std::mem::replace(&mut self.lists[0], incoming);
            heapify_down(&mut self.lists, 0);
            old
        };
        self.worst_evalue = self.lists[0].hsps.best_evalue;
        self.low_score = self.lists[0].hsps.hsps[0].hsp.score;
        Ok(Some(discarded))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3420-3437:
    pub(super) fn reverse_order(&mut self) {
        self.lists.reverse();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1776;
    // c++/src/algo/blast/core/blast_hits.c:3383-3400:
    pub(super) fn sort_by_evalue_for_positive_totlen(&mut self) {
        self.lists.sort_by(compare_result_lists);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:827-857;
    // c++/src/algo/blast/core/blast_hits.c:3331-3337 (`core::blast_hspstream`)
    /// `BlastHSPCBSStreamClose`'s trim of this query's preliminary hit list (`ADAPTIVE_CBS`):
    /// sorted by E-value (`Blast_HitListSortByEvalue`: the lists left without HSPs, last,
    /// removed) and cut where the E-values pass the limit, when it holds enough lists.
    pub(super) fn cbs_close_trim(&mut self, hitlist_size: usize) {
        use crate::core::blast_hspstream::{cbs_trim_count, cbs_trim_sorts};
        let count = self.lists.len();
        if !cbs_trim_sorts(count, hitlist_size) {
            return;
        }
        self.lists.sort_unstable_by(compare_result_lists);
        let sorted = self
            .lists
            .iter()
            .position(|list| list.hsps.hsps.is_empty())
            .unwrap_or(count);
        self.lists.truncate(sorted);
        self.heapified = false;
        let lists = &self.lists;
        if let Some(keep) =
            cbs_trim_count(count, sorted, hitlist_size, |i| lists[i].hsps.best_evalue)
        {
            self.lists.truncate(keep);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577:
    // Sort at report conversion, after Kappa has finished using the list.
    pub(super) fn sort_hsps_for_report(&mut self) {
        for list in &mut self.lists {
            sort_hsps_by_evalue(list);
            seqalign_scores_of(list);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1485-1490
    /// Leaves out the HSPs of score 0 (kept by a very large -evalue), which have no
    /// Seq-align and so are in no report, and the subjects left without an HSP.
    pub(super) fn drop_zero_score_hsps(&mut self) {
        for list in &mut self.lists {
            if list.payloads.is_empty() {
                list.hsps.hsps.retain(|linked| linked.hsp.score != 0);
            } else {
                let (hsps, payloads): (Vec<_>, Vec<_>) = std::mem::take(&mut list.hsps.hsps)
                    .into_iter()
                    .zip(std::mem::take(&mut list.payloads))
                    .filter(|(linked, _)| linked.hsp.score != 0)
                    .unzip();
                list.hsps.hsps = hsps;
                list.payloads = payloads;
            }
        }
        self.lists.retain(|list| !list.hsps.hsps.is_empty());
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:169-181:
    // Stores the bounded array of HSP lists (`hsplist_array`) with its current count and capacity.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:845-865
    /// `s_FilterBlastResults` of `-max_hsps` (`max_hsps`, 0 keeps every HSP) and
    /// `-qcov_hsp_perc` (`query_cov_hsp_perc`, 0 removes none) on the final hit list of a
    /// query of `query_length` residues: each list cut to its first HSPs, then its HSPs below
    /// the coverage removed (its best e-value set again when one is), and the lists left
    /// empty removed. The payloads follow their HSPs.
    pub(super) fn filter_results(
        &mut self,
        max_hsps: usize,
        query_cov_hsp_perc: f64,
        query_length: usize,
    ) {
        use crate::core::blast_hits::{reap_by_query_coverage, trim_by_max_hsps};
        let extent = |linked: &LinkedHsp| {
            (
                i64::from(linked.hsp.q_start),
                i64::from(linked.hsp.q_end),
                query_length as i64,
            )
        };
        for list in &mut self.lists {
            trim_by_max_hsps(&mut list.hsps.hsps, max_hsps);
            trim_by_max_hsps(&mut list.payloads, max_hsps);
            let removed = if query_cov_hsp_perc <= 0.0 {
                false
            } else if list.payloads.is_empty() {
                reap_by_query_coverage(&mut list.hsps.hsps, query_cov_hsp_perc, extent)
            } else {
                let mut paired: Vec<_> = std::mem::take(&mut list.hsps.hsps)
                    .into_iter()
                    .zip(std::mem::take(&mut list.payloads))
                    .collect();
                let removed =
                    reap_by_query_coverage(&mut paired, query_cov_hsp_perc, |(linked, _)| {
                        extent(linked)
                    });
                (list.hsps.hsps, list.payloads) = paired.into_iter().unzip();
                removed
            };
            if removed {
                refresh_best_evalue(&mut list.hsps);
            }
        }
        if query_cov_hsp_perc > 0.0 {
            self.lists.retain(|list| !list.hsps.hsps.is_empty());
        }
    }

    /// `-subject_besthit` on each final list of the query (`s_FilterBlastResults`, after
    /// `-qcov_hsp_perc`; a query is one context): the later HSPs within a kept HSP's query range
    /// widened by 3 go, with their payloads. The best e-value is not set again, as NCBI's is not.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:858-861
    pub(super) fn subject_besthit(&mut self) {
        use crate::core::blast_hits::{subject_besthit, BestHitExtent};
        let extent = |linked: &LinkedHsp| BestHitExtent {
            context: 0,
            query_frame: 0,
            offset: i64::from(linked.hsp.q_start),
            end: i64::from(linked.hsp.q_end),
            query_length: 0,
        };
        for list in &mut self.lists {
            if list.payloads.is_empty() {
                subject_besthit(&mut list.hsps.hsps, false, extent);
            } else {
                let mut paired: Vec<_> = std::mem::take(&mut list.hsps.hsps)
                    .into_iter()
                    .zip(std::mem::take(&mut list.payloads))
                    .collect();
                subject_besthit(&mut paired, false, |(linked, _)| extent(linked));
                (list.hsps.hsps, list.payloads) = paired.into_iter().unzip();
            }
        }
    }

    pub(super) fn lists(&self) -> &[KappaResultList] {
        &self.lists
    }

    /// The lists, for a caller that renames their subjects (`database::reported_subjects`:
    /// ordinal ids to ranks, in the same order).
    pub(super) fn lists_mut(&mut self) -> &mut [KappaResultList] {
        &mut self.lists
    }

    /// The hit list after `cut` removes lists from it (`database::drop_unaligned`).
    pub(super) fn with_lists_mut(&mut self, cut: impl FnOnce(&mut Vec<KappaResultList>)) {
        cut(&mut self.lists);
    }

    pub(super) fn worst_evalue(&self) -> f64 {
        self.worst_evalue
    }

    pub(super) fn low_score(&self) -> i32 {
        self.low_score
    }

    pub(super) fn heapified(&self) -> bool {
        self.heapified
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1734-1750:
fn refresh_best_evalue(list: &mut LinkedHspList) {
    list.best_evalue = list
        .hsps
        .iter()
        .fold(i32::MAX as f64, |best, hsp| best.min(hsp.evalue));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1385-1455:
/// The order of a list's HSPs in the reports (`sort_hsps_by_evalue`, `evalue_order`).
fn compare_report_hsps(a: &LinkedHsp, b: &LinkedHsp) -> Ordering {
    let evalue = if a.evalue < 1.0e-180 && b.evalue < 1.0e-180 {
        Ordering::Equal
    } else if a.evalue < b.evalue {
        Ordering::Less
    } else if a.evalue > b.evalue {
        Ordering::Greater
    } else {
        Ordering::Equal
    };
    evalue.then_with(|| score_compare(&a.hsp, &b.hsp))
}

/// The indices of `list`'s HSPs (and payloads) in the order `sort_hsps_by_evalue` puts them
/// (`KappaResultHitList::sort_hsps_for_report`, before the reports), without moving them: the
/// same stable sort by `compare_report_hsps`.
pub(super) fn evalue_order(list: &KappaResultList) -> Vec<usize> {
    let hsps = &list.hsps.hsps;
    let mut order: Vec<usize> = (0..hsps.len()).collect();
    order.sort_by(|&a, &b| compare_report_hsps(&hsps[a], &hsps[b]));
    order
}

fn sort_hsps_by_evalue(list: &mut KappaResultList) {
    let compare = compare_report_hsps;
    if list.payloads.is_empty() {
        list.hsps.hsps.sort_by(compare);
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1385-1455:
        // Sorts HSP pointers via `ScoreCompareHSPs` so each bit score, identity count,
        // and edit script moves with its HSP.
        let mut paired: Vec<_> = std::mem::take(&mut list.hsps.hsps)
            .into_iter()
            .zip(std::mem::take(&mut list.payloads))
            .collect();
        paired.sort_by(|a, b| compare(&a.0, &b.0));
        let (hsps, payloads) = paired.into_iter().unzip();
        list.hsps.hsps = hsps;
        list.payloads = payloads;
    }
}

/// The E-values and bit scores of a list's HSPs as the Seq-align that NCBI makes of them carries
/// them (`seqalign_scores`: a negative E-value is 0, a NaN one stays as the E-value that the
/// Seq-align lacks, and a bit score below 0, which a PSSM whose gapped Karlin-Altschul block is
/// not valid gives, is the readers' -1).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1185-1197
fn seqalign_scores_of(list: &mut KappaResultList) {
    for linked in &mut list.hsps.hsps {
        linked.evalue = crate::report::serial::seqalign::seqalign_scores(linked.evalue, 0.0).0;
    }
    for payload in &mut list.payloads {
        payload.bit_score =
            crate::report::serial::seqalign::seqalign_scores(0.0, payload.bit_score).1;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3071-3111:
fn compare_result_lists(a: &KappaResultList, b: &KappaResultList) -> Ordering {
    compare_preliminary_lists_for_kappa(a.oid, &a.hsps, b.oid, &b.hsps)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1627-1676,
// 3194-3210:
fn heapify_down(lists: &mut [KappaResultList], start: usize) {
    let mut root = start;
    loop {
        let left = root * 2 + 1;
        if left >= lists.len() {
            break;
        }
        let right = left + 1;
        let large = if right < lists.len()
            && compare_result_lists(&lists[left], &lists[right]) == Ordering::Less
        {
            right
        } else {
            left
        };
        if compare_result_lists(&lists[root], &lists[large]) != Ordering::Less {
            break;
        }
        lists.swap(root, large);
        root = large;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:476-534
/// NCBI's `BlastHSPStreamMerge` of a query chunk (RP-4): the HSPs of each part's preliminary
/// hit list move to the part's query (a protein query has one context and frame 0), and the
/// list merges with the query's hit list of the chunks before. `chunk_lists[i]` is the hit
/// list of the chunk's part `parts[i] = (query, offset)`. The split blocks of a gapped search
/// allow gaps (split_query_cxx.cpp:864-865, `GetGappedMode`).
pub(super) fn merge_query_chunk(
    merged: &mut [Option<KappaResultHitList>],
    chunk_lists: Vec<Option<KappaResultHitList>>,
    parts: &[(usize, i32)],
    chunk_overlap_size: i32,
) -> Result<()> {
    ensure!(
        chunk_lists.len() == parts.len(),
        "one preliminary hit list is required per query part of a chunk"
    );
    for (hitlist, &(query, offset)) in chunk_lists.into_iter().zip(parts) {
        let Some(mut hitlist) = hitlist else {
            continue;
        };
        for list in &mut hitlist.lists {
            for linked in &mut list.hsps.hsps {
                linked.context = query;
                linked.hsp.q_start += offset;
                linked.hsp.q_end += offset;
                linked.hsp.q_gapped_start += offset;
            }
        }
        hit_list_merge(hitlist, &mut merged[query], offset, chunk_overlap_size)?;
    }
    for hitlist in merged.iter_mut().flatten() {
        for list in &mut hitlist.lists {
            sort_hsps_by_score(list);
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2132-2217
// `Blast_HitListUpdate` keeps the best `prelim_hitlist_size` subjects; the subjects of a hit
// list are distinct, so the OID sort is total.
fn hit_list_merge(
    mut hitlist1: KappaResultHitList,
    combined: &mut Option<KappaResultHitList>,
    split_offset: i32,
    chunk_overlap_size: i32,
) -> Result<()> {
    let Some(mut hitlist2) = combined.take() else {
        *combined = Some(hitlist1);
        return Ok(());
    };
    let mut new_hitlist = KappaResultHitList::new(hitlist1.max)?;
    hitlist1.lists.sort_by_key(|list| list.oid);
    hitlist2.lists.sort_by_key(|list| list.oid);
    let query_is_split = split_offset > 0;
    let mut lists1 = std::mem::take(&mut hitlist1.lists).into_iter().peekable();
    let mut lists2 = std::mem::take(&mut hitlist2.lists).into_iter().peekable();
    while let (Some(list1), Some(list2)) = (lists1.peek(), lists2.peek()) {
        match list1.oid.cmp(&list2.oid) {
            Ordering::Less => {
                new_hitlist.update(lists1.next().expect("peeked"))?;
            }
            Ordering::Greater => {
                new_hitlist.update(lists2.next().expect("peeked"))?;
            }
            Ordering::Equal => {
                let list1 = lists1.next().expect("peeked");
                let mut list2 = lists2.next().expect("peeked");
                if query_is_split {
                    query_hsp_lists_merge(list1, &mut list2, split_offset, chunk_overlap_size);
                } else {
                    hsp_list_append(list1, &mut list2);
                }
                new_hitlist.update(list2)?;
            }
        }
    }
    for list in lists1.chain(lists2) {
        new_hitlist.update(list)?;
    }
    *combined = Some(new_hitlist);
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2921-3029
// The query split branch of `Blast_HSPListsMerge` for a protein query (one context, frame 0),
// with gaps allowed; TBLASTN has no `-max_hsps`, so `hsp_num_max` keeps every HSP.
fn query_hsp_lists_merge(
    new: KappaResultList,
    combined: &mut KappaResultList,
    split_offset: i32,
    chunk_overlap_size: i32,
) {
    if new.hsps.hsps.is_empty() {
        return;
    }
    let mut incoming = new.hsps.hsps;
    let old = &mut combined.hsps.hsps;
    let mut hspcnt1 = 0;
    if split_offset >= 0 {
        for index1 in 0..old.len() {
            if old[index1].hsp.q_end > split_offset {
                old.swap(hspcnt1, index1);
                hspcnt1 += 1;
            }
        }
    }
    let mut hspcnt2 = 0;
    if split_offset >= 0 {
        for index2 in 0..incoming.len() {
            if incoming[index2].hsp.q_start < split_offset + chunk_overlap_size {
                incoming.swap(hspcnt2, index2);
                hspcnt2 += 1;
            }
        }
    }
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut slots: Vec<Option<LinkedHsp>> = incoming.into_iter().map(Some).collect();
        for index1 in 0..hspcnt1 {
            for slot in slots.iter_mut().take(hspcnt2) {
                let Some(candidate) = slot.as_ref() else {
                    continue;
                };
                if old[index1].context != candidate.context {
                    continue;
                }
                let end_diag = old[index1].hsp.q_end - old[index1].hsp.s_end;
                let start_diag = candidate.hsp.q_start - candidate.hsp.s_start;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1537
                if (end_diag - start_diag).abs() < 10
                    && super::search_gapped::merge_two_chunk_hsps(
                        &mut old[index1].hsp,
                        &candidate.hsp,
                        true,
                    )
                {
                    *slot = None;
                }
            }
        }
        incoming = slots.into_iter().flatten().collect();
    }
    combine_by_score(incoming, combined);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2819-2829,2849
fn hsp_list_append(new: KappaResultList, combined: &mut KappaResultList) {
    combine_by_score(new.hsps.hsps, combined);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
fn combine_by_score(incoming: Vec<LinkedHsp>, combined: &mut KappaResultList) {
    combined.hsps.hsps.extend(incoming);
    sort_hsps_by_score(combined);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
// Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5).
fn sort_hsps_by_score(list: &mut KappaResultList) {
    let hsps = &mut list.hsps.hsps;
    if hsps
        .windows(2)
        .any(|pair| score_compare(&pair[0].hsp, &pair[1].hsp) == Ordering::Greater)
    {
        hsps.sort_by(|a, b| score_compare(&a.hsp, &b.hsp));
    }
}

#[cfg(test)]
#[path = "stage_d_results_tests.rs"]
mod tests;
