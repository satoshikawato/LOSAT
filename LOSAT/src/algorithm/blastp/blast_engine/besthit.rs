//! BLASTP `-best_hit_overhang` and `-best_hit_score_edge` (psiblast's too): NCBI's best hit
//! writer of the preliminary stage and its best hit pipe after the traceback
//! (`core::hspfilter_besthit` holds the lists both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//!
//! This module is a child of `blast_engine` for the subject result type that it keeps private.
//! Nothing here runs without one of the two options: the choice is made once per search call,
//! where `call_setup::prelim_sinks` picks the writer of the preliminary hit lists
//! (`PrelimFilter`), and once per query where `blast_engine` finishes a final hit list.

use super::BlastpSubjectResult;
use crate::algorithm::blastp::hsp::{
    score_compare_hsps, update_best_evalue, BlastpHitList, BlastpHsp, BlastpHspList, HitListHsp,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::core::hspfilter_besthit::{BestHitKey, BestHitLists, BestHitParams};

/// The key of an HSP for the lists.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-243
/// A blastp context has frame 0, so the start is the HSP's own `query.offset` (`q_start - 1`,
/// within the searched range of the query or the query part of a chunk) and the query length is
/// not read.
fn best_hit_key<H: HitListHsp>(hsp: &H) -> BestHitKey {
    let (q_start, q_end, _, _) = hsp.offsets();
    let begin = q_start as i64 - 1;
    BestHitKey {
        begin: begin as i32,
        len: (q_end as i64 - begin) as i32,
        score: hsp.raw_score(),
        evalue: hsp.e_value(),
    }
}

/// One subject's HSP list of one query into the lists, in the order of its HSPs.
fn save_list<H: HitListHsp>(lists: &mut BestHitLists<H>, query: usize, list: BlastpHspList<H>) {
    let oid = list.oid;
    for hsp in list.hsps {
        let key = best_hit_key(&hsp);
        lists.save(query, oid, key, hsp);
    }
}

/// `s_BlastHSPBestHitFinal`: the hit list of each query with kept HSPs goes into `out`: the
/// kept lists of the query (one per subject in the order of its first HSP in the query's list,
/// each sorted by score), with their best e-values, sorted by e-value and saved into a hit list
/// of `prelim_hitlist_size` lists, which keeps the lists of the best e-values when there are
/// more (`Blast_HitListUpdate`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:190-221
/// `s_ExportToHitlist` saves the lists into a hit list as large as the HSPs (it is never full,
/// so `Blast_HitListUpdate` only sets each list's best e-value) and `Blast_HitListSortByEvalue`
/// sorts it.
fn finish_lists<H: HitListHsp>(
    lists: BestHitLists<H>,
    prelim_hitlist_size: usize,
    mut out: impl FnMut(usize, BlastpHitList<H>),
) {
    for (query, kept) in lists.finalize(score_compare_hsps::<H>) {
        let mut exported = BlastpHitList::new(kept.len());
        for (oid, hsps) in kept {
            exported.update(BlastpHspList {
                oid,
                query_index: query as u32,
                hsps,
                best_evalue: 0.0,
            });
        }
        exported.sort_by_evalue();
        let mut hit_list = BlastpHitList::new(prelim_hitlist_size);
        for list in exported.hsplist_array {
            hit_list.update(list);
        }
        out(query, hit_list);
    }
}

/// NCBI's best hit writer in place of the HSP collector (`update_blastp_preliminary_hit_lists`):
/// each subject's preliminary HSPs go into the list of their query, and `finish` makes the hit
/// lists of the queries. The subjects must come in OID order and the HSPs of a subject in the
/// order of its score-sorted HSP list (the order in which the collector receives them), because
/// the result depends on the order of the writes.
///
/// The preliminary lists are the compact lists (`x_prelim_compact.rs`), so the lists hold
/// `XPrelimHsp`s.
pub(super) struct PrelimBestHit {
    lists: BestHitLists<XPrelimHsp>,
    prelim_hitlist_size: usize,
}

impl PrelimBestHit {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(super) fn new(prelim_hitlist_size: usize, params: BestHitParams) -> Self {
        Self {
            lists: BestHitLists::writer(params, prelim_hitlist_size),
            prelim_hitlist_size,
        }
    }

    /// `s_BlastHSPBestHitRun` for the HSPs of one subject (all queries of the call; in a chunk
    /// call the parts of the chunk are its queries). The subject result is emptied, as the
    /// collector empties it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:228-333
    pub(super) fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        for (query, list) in subject_result.x_compact_hsp_lists.drain(..) {
            save_list(&mut self.lists, query, list);
        }
    }

    /// `s_BlastHSPBestHitFinal`: the hit list of each query with kept HSPs goes into
    /// `compact_lists`, whose lists the collector would have filled.
    pub(super) fn finish(self, compact_lists: &mut [Option<BlastpHitList<XPrelimHsp>>]) {
        finish_lists(self.lists, self.prelim_hitlist_size, |query, hit_list| {
            compact_lists[query] = Some(hit_list);
        });
    }
}

/// NCBI's best hit pipe on the final hit list of one query, the first step after the traceback
/// (`BlastHSPStreamTBackClose`, before the mask level, `-max_hsps`, `-qcov_hsp_perc` and
/// `-subject_besthit`): the HSP lists sorted by their best e-values (the HSPs of each list keep
/// their order), every HSP in that order put into a new list with the overhang and the score
/// edge of the options, and the kept HSPs made into the query's hit list again.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539
/// The lists are per query, so the pipe of a batch is this function for each of its queries.
/// The HSPs are those of the traceback, the redo of `-comp_based_stats` or the ungapped stage:
/// the offsets are those within the searched range of the query.
pub(super) fn best_hit_traceback(
    hit_list: &mut BlastpHitList,
    params: BestHitParams,
    prelim_hitlist_size: usize,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
    // `Blast_HSPResultsSortByEvalue`: the lists by their best e-values (the smallest e-value of
    // each list, as the traceback's `Blast_HitListUpdate` set it).
    hit_list.hsplist_array.truncate(hit_list.hsplist_count);
    for list in &mut hit_list.hsplist_array {
        update_best_evalue(list);
    }
    hit_list.sort_by_evalue();
    let query_index = hit_list
        .hsplist_array
        .first()
        .map_or(0, |list| list.query_index);

    let mut lists = BestHitLists::<BlastpHsp>::pipe(params, prelim_hitlist_size);
    for list in std::mem::take(&mut hit_list.hsplist_array) {
        save_list(&mut lists, 0, list);
    }
    let mut new_list = BlastpHitList::new(prelim_hitlist_size);
    finish_lists(lists, prelim_hitlist_size, |_, finished| {
        new_list = finished;
    });
    for list in &mut new_list.hsplist_array {
        list.query_index = query_index;
    }
    *hit_list = new_list;
}

#[cfg(test)]
#[path = "besthit_tests.rs"]
mod tests;
