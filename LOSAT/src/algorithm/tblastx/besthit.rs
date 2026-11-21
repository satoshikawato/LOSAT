//! TBLASTX `-best_hit_overhang` and `-best_hit_score_edge`: NCBI's best hit writer of the
//! preliminary stage and its best hit pipe after the (skipped) traceback
//! (`core::hspfilter_besthit` holds the lists both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329,377-385
//!
//! TBLASTX is ungapped and has no composition-based statistics, so the preliminary hit list has
//! the final size (`GetPrelimHitlistSize`: `-max_target_seqs`) in the writer, the pipe and the
//! traceback alike. Nothing here runs without one of the two options: the hit collector of a
//! query batch picks this writer once, in place of NCBI's HSP collector
//! (`HitCollector::new`).

use crate::algorithm::blastn::hsp::{HitList, HitListEntry};
use crate::core::blast_hits::trim_by_max_hsps;
use crate::core::hspfilter_besthit::{
    translated_query_length, BestHitKey, BestHitLists, BestHitParams, BestHitQuery,
};

use super::blast_engine::TblastxHsp;
use super::hit_collector::TblastxHspList;
use super::hsp_culling::tblastx_context_lengths;
use super::report::HitFilters;

/// The key of an HSP for the lists: its range in its frame (`query.offset`, `query.end` of the
/// HSP: the frame offsets of `Blast_HSPInit`), folded over `qlen` for a minus frame.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266
/// `qlen` is the query's `BlastQueryInfoGetQueryLength`, 2 plus the lengths of its three frames
/// (`translated_query_length`): a DNA length, while the offsets are those of the frame's
/// residues, as in NCBI.
fn best_hit_key(hsp: &TblastxHsp, qlen: i32) -> BestHitKey {
    let hit = &hsp.hit;
    BestHitKey::new(
        hit.query_frame < 0,
        qlen,
        hit.sort_query_offset as i32,
        hit.sort_query_end as i32,
        hit.raw_score,
        hit.e_value,
    )
}

/// The kept lists of a query as `TblastxHspList`s, each sorted by score
/// (`Blast_HSPListSortByScore` in `s_BlastHSPBestHitFinal`).
fn kept_lists(kept: BestHitQuery<TblastxHsp>) -> Vec<TblastxHspList> {
    kept.into_iter()
        .map(|(oid, hsps)| TblastxHspList {
            oid,
            hsps,
            best_evalue: 0.0,
        })
        .collect()
}

fn score_order(a: &TblastxHsp, b: &TblastxHsp) -> std::cmp::Ordering {
    crate::common::score_compare_hsps(&a.hit, &b.hit)
}

/// NCBI's best hit writer in place of the HSP collector of a query batch: each subject's HSPs of
/// a query (in the order of the subject's list, after the e-value and coverage reaps) go into the
/// list of the query, and `finish` runs the rest of NCBI's search on the kept HSPs. The subjects
/// must come in OID order (`HitCollector::push`), because the result depends on the order of the
/// writes.
pub(super) struct BestHitWriter {
    lists: BestHitLists<TblastxHsp>,
    params: BestHitParams,
    /// `BlastQueryInfoGetQueryLength` of each query of the batch.
    query_lengths: Vec<i32>,
    hitlist_size: usize,
}

impl BestHitWriter {
    /// The writer of a batch whose queries have `nucleotide_lengths` (the searched lengths).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
    pub(super) fn new(
        params: BestHitParams,
        hitlist_size: usize,
        nucleotide_lengths: &[usize],
    ) -> Self {
        Self {
            lists: BestHitLists::writer(params, hitlist_size),
            params,
            query_lengths: nucleotide_lengths
                .iter()
                .map(|&length| translated_query_length(&tblastx_context_lengths(length)))
                .collect(),
            hitlist_size,
        }
    }

    /// `s_BlastHSPBestHitRun` for subject `oid`'s HSPs of query `q_idx`, in their order.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-352
    pub(super) fn add(&mut self, q_idx: u32, oid: u32, hsps: Vec<TblastxHsp>) {
        let query = q_idx as usize;
        let qlen = self.query_lengths[query];
        for hsp in hsps {
            let key = best_hit_key(&hsp, qlen);
            self.lists.save(query, oid, key, hsp);
        }
    }

    /// The batch's HSPs in report order: per query, the writer's hit list
    /// (`s_BlastHSPBestHitFinal`), the traceback's hit list of its lists in OID order, the pipe on
    /// it, then `-max_hsps`, `-subject_besthit`, the e-value sort (`sorts_results`) and the prune
    /// to `-max_target_seqs`, as `report::culled_hit_order` ends.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1700-1706 (the
    /// ungapped traceback inserts each list of the stream, read in OID order, into the query's
    /// hit list of `hitlist_size`), 1720-1777 (`BlastHSPStreamTBackClose`, then the filters, the
    /// sort and the prune)
    pub(super) fn finish(self, filters: &HitFilters<'_>, sorts_results: bool) -> Vec<TblastxHsp> {
        let Self {
            lists,
            params,
            query_lengths,
            hitlist_size,
        } = self;
        let mut ordered = Vec::new();
        for (query, kept) in lists.finalize(score_order) {
            let prelim = HitList::best_hit_final(kept_lists(kept), hitlist_size);
            let mut lists: Vec<TblastxHspList> = prelim
                .hsplist_array
                .into_iter()
                .take(prelim.hsplist_count)
                .collect();
            lists.sort_by_key(|list| list.oid);
            let mut traceback = HitList::new(hitlist_size);
            for list in lists {
                traceback.update(list);
            }
            let mut final_list =
                best_hit_pipe(traceback, params, hitlist_size, query_lengths[query]);
            for list in &mut final_list.hsplist_array {
                trim_by_max_hsps(&mut list.hsps, filters.max_hsps);
                filters.subject_besthit(&mut list.hsps);
            }
            if sorts_results {
                final_list.sort_by_evalue();
            }
            final_list.prune_by_size(hitlist_size);
            for mut list in final_list.hsplist_array {
                HitListEntry::sort_by_evalue(&mut list);
                ordered.extend(list.hsps);
            }
        }
        ordered
    }
}

/// NCBI's best hit pipe on the hit list of one query (`BlastHSPStreamTBackClose`, before
/// `-max_hsps` and `-subject_besthit`): the lists sorted by their best e-values (the HSPs of each
/// list keep their order: by score, or by e-value where the writer's hit list became a heap),
/// every HSP in that order put into a new list with the overhang and the score edge, and the
/// kept HSPs made into the query's hit list again.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
/// (`Blast_HSPResultsSortByEvalue`: the best e-values that the traceback's `Blast_HitListUpdate`
/// set)
fn best_hit_pipe(
    mut hit_list: HitList<TblastxHspList>,
    params: BestHitParams,
    hitlist_size: usize,
    qlen: i32,
) -> HitList<TblastxHspList> {
    hit_list.sort_by_evalue();
    let mut lists = BestHitLists::pipe(params, hitlist_size);
    for list in hit_list
        .hsplist_array
        .into_iter()
        .take(hit_list.hsplist_count)
    {
        let oid = list.oid;
        for hsp in list.hsps {
            let key = best_hit_key(&hsp, qlen);
            lists.save(0, oid, key, hsp);
        }
    }
    let kept = lists
        .finalize(score_order)
        .into_iter()
        .flat_map(|(_, kept)| kept)
        .collect();
    HitList::best_hit_final(kept_lists(kept), hitlist_size)
}

#[cfg(test)]
#[path = "besthit_tests.rs"]
mod tests;
