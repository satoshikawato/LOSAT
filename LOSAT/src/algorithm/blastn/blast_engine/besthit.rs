//! BLASTN `-best_hit_overhang` and `-best_hit_score_edge`: NCBI's best hit writer of the
//! preliminary stage and its best hit pipe after the traceback (`core::hspfilter_besthit` holds
//! the lists both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//!
//! This module is a child of `run` for the preliminary HSP types that `run` keeps private.
//! Nothing here runs without one of the two options: the choice is made once per search, where
//! `run` picks the sink of the preliminary hit lists, and once per query where it finishes a
//! final hit list.

use super::{score_compare_prelim_hits, PrelimHit, PrelimHitLists, PrelimHspList};
use crate::algorithm::blastn::hsp::{
    score_compare_hsps, update_best_evalue, BlastnHitList, BlastnHspList, HitList,
};
use crate::core::hspfilter_besthit::{BestHitKey, BestHitLists, BestHitParams, BestHitQuery};

/// An `Int4` length as NCBI's `query_length` holds it (the lengths of LOSAT's queries fit).
fn int4(length: usize) -> i32 {
    i32::try_from(length).unwrap_or(i32::MAX)
}

/// NCBI's best hit writer in place of the HSP collector (`PrelimHitCollector`): each subject's
/// preliminary HSPs go into the list of their query, and `finish` makes the hit lists of the
/// queries. The subjects must come in OID order and the HSPs of a subject in the order of its
/// score-sorted HSP list (the order in which the collector receives them), because the result
/// depends on the order of the writes.
pub(super) struct PrelimBestHit {
    lists: BestHitLists<PrelimHit>,
    /// The length of each query (of each query part in the chunk stage of a split batch).
    query_lengths: Vec<i32>,
    prelim_hitlist_size: usize,
}

impl PrelimBestHit {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(super) fn new(
        query_lengths: &[usize],
        prelim_hitlist_size: usize,
        params: BestHitParams,
    ) -> Self {
        Self {
            lists: BestHitLists::writer(params, prelim_hitlist_size),
            query_lengths: query_lengths.iter().map(|&length| int4(length)).collect(),
            prelim_hitlist_size,
        }
    }

    /// `s_BlastHSPBestHitRun` for the HSPs of subject `oid` (all queries, both strands).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-352
    pub(super) fn add(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        for hit in hits {
            let query = hit.query_idx as usize;
            let key = BestHitKey::new(
                hit.context_idx % 2 == 1,
                self.query_lengths[query],
                hit.prelim_qs,
                hit.prelim_qe,
                hit.prelim_score,
                hit.prelim_evalue,
            );
            self.lists.save(query, oid as u32, key, hit);
        }
    }

    /// `s_BlastHSPBestHitFinal`: the hit list of each query with kept HSPs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231
    pub(super) fn finish(self) -> PrelimHitLists {
        let mut hit_lists: PrelimHitLists = Vec::with_capacity(self.query_lengths.len());
        hit_lists.resize_with(self.query_lengths.len(), || None);
        let kept = self.lists.finalize(score_compare_prelim_hits);
        for (query, lists) in kept {
            let lists: Vec<PrelimHspList> = lists
                .into_iter()
                .map(|(oid, hsps)| PrelimHspList {
                    oid,
                    hsps,
                    best_evalue: 0.0,
                })
                .collect();
            hit_lists[query] = Some(HitList::best_hit_final(lists, self.prelim_hitlist_size));
        }
        hit_lists
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
/// `query_len` is the length of each context of the query (`finish_query_hit_list`). The HSP
/// offsets are those within the HSP's context and the context is minus for a negative frame.
pub(super) fn best_hit_traceback(
    hit_list: &mut BlastnHitList,
    query_len: usize,
    params: BestHitParams,
    prelim_hitlist_size: usize,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
    // `Blast_HSPResultsSortByEvalue`: the lists by their best e-values (the smallest e-value of
    // each list, as the traceback's `Blast_HitListUpdate` set it).
    for list in &mut hit_list.hsplist_array {
        update_best_evalue(list);
    }
    hit_list.sort_by_evalue();
    let query_index = hit_list
        .hsplist_array
        .first()
        .map_or(0, |list| list.query_index);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:525-532
    let qlen = int4(query_len);
    let mut lists = BestHitLists::pipe(params, prelim_hitlist_size);
    for list in std::mem::take(&mut hit_list.hsplist_array) {
        let oid = list.oid;
        for hsp in list.hsps {
            let key = BestHitKey::new(
                hsp.query_frame < 0,
                qlen,
                hsp.internal_q_offset_0,
                hsp.internal_q_end_0,
                hsp.raw_score,
                hsp.e_value,
            );
            lists.save(0, oid, key, hsp);
        }
    }

    // `s_BlastHSPBestHitFinal`: the new hit list (as `PrelimBestHit::finish`, with
    // `ScoreCompareHSPs` on the final HSPs).
    let kept: BestHitQuery<_> = lists
        .finalize(score_compare_hsps)
        .into_iter()
        .flat_map(|(_, lists)| lists)
        .collect();
    let lists: Vec<BlastnHspList> = kept
        .into_iter()
        .map(|(oid, hsps)| BlastnHspList {
            oid,
            query_index,
            hsps,
            best_evalue: 0.0,
        })
        .collect();
    *hit_list = HitList::best_hit_final(lists, prelim_hitlist_size);
}

#[cfg(test)]
#[path = "besthit_tests.rs"]
mod tests;
