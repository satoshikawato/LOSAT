//! BLASTP `-culling_limit`: NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback (`core::hspfilter_culling` holds the interval tree both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c
//!
//! This module is a child of `blast_engine` for the subject result type that it keeps private.
//! Nothing here runs without `-culling_limit` (a limit of 0 is NCBI's `SetCullingLimit(0)`,
//! which is ignored): the choice is made once per search call, where `blast_engine` picks the
//! sink of the preliminary hit lists, and once per query where it finishes a final hit list.

use super::BlastpSubjectResult;
use crate::algorithm::blastp::hsp::{
    score_compare_hsps, sort_hsplist_by_evalue, update_best_evalue, BlastpHitList, BlastpHsp,
    BlastpHspList, HitListHsp,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::core::hspfilter_culling::{CullKey, CullingForest};

/// An `Int4` length as NCBI's `query_length` holds it (the lengths of LOSAT's queries fit).
fn int4(length: usize) -> i32 {
    i32::try_from(length).unwrap_or(i32::MAX)
}

/// The key of an HSP for the tree and the length of the tree's query.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:612-633
/// Not blastn: the range is the HSP's own, `query.offset` = `q_start - 1`, `query.end` =
/// `q_end` (`HitListHsp::offsets`), and `subject.offset` = `s_start - 1`, which
/// `s_DominateTest` reads for the tie of two HSPs of one subject (hspfilter_culling.c:100-112).
/// `qlen` is the length of the HSP's context: the searched range of the query, or the part of
/// it in a query chunk (`HitListHsp::query_length`).
fn cull_key<H: HitListHsp>(hsp: &H) -> (i32, CullKey) {
    let (q_start, q_end, s_start, _) = hsp.offsets();
    let key = CullKey {
        begin: int4(q_start.saturating_sub(1)),
        end: int4(q_end),
        score: hsp.raw_score(),
        subject_offset: s_start.saturating_sub(1) as i64,
    };
    (int4(hsp.query_length()), key)
}

/// NCBI's culling writer in place of the HSP collector (`update_blastp_preliminary_hit_lists`):
/// each subject's preliminary HSPs go into the tree of their query, and `finish` makes the hit
/// lists of the queries. The subjects must come in OID order and the HSPs of a subject in the
/// order of its score-sorted HSP list (the order in which the collector receives them), because
/// the result depends on the order of the writes.
///
/// The preliminary lists are the compact lists (`x_prelim_compact.rs`), so the forest holds
/// `XPrelimHsp`s.
pub(super) struct PrelimCulling {
    compact: CullingForest<XPrelimHsp>,
    prelim_hitlist_size: usize,
}

impl PrelimCulling {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:318-342
    pub(super) fn new(prelim_hitlist_size: usize, culling_limit: i32) -> Self {
        let culling_max = CullingForest::<BlastpHsp>::prelim_max(culling_limit);
        Self {
            compact: CullingForest::new(culling_max),
            prelim_hitlist_size,
        }
    }

    /// `s_BlastHSPCullingRun` for the HSPs of one subject (all queries of the call).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// The tree id is the query of the call (blastp has one context per query; in a chunk call
    /// the parts of the chunk are its queries). An HSP that is not saved is freed. The subject
    /// result is emptied, as the collector empties it.
    pub(super) fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        for (query, list) in subject_result.x_compact_hsp_lists.drain(..) {
            save_list(&mut self.compact, query, list);
        }
    }

    /// `s_BlastHSPCullingFinal`: the hit list of each query with a tree goes into
    /// `compact_lists`. It has one HSP list for each subject in the
    /// order of its first HSP in the rip of the tree, each sorted by score, and as many lists as
    /// the tree keeps: the hit list is made with `prelim_hitlist_size` as its maximum but is not
    /// cut to it, and it is not heapified.
    ///
    /// The merge of a split query (`query_split.rs`, `Blast_HitListMerge`) reads `hsplist_max`
    /// and the lists; the first chunk's hit list is moved as it is, with more lists than
    /// `hsplist_max`. The stream close (`sort_by_evalue`) sorts and purges the first
    /// `hsplist_count` lists.
    pub(super) fn finish(self, compact_lists: &mut [Option<BlastpHitList<XPrelimHsp>>]) {
        finish_forest(self.compact, self.prelim_hitlist_size, compact_lists);
    }
}

/// One subject's HSP list of one query into the forest.
fn save_list<H: HitListHsp>(forest: &mut CullingForest<H>, query: usize, list: BlastpHspList<H>) {
    let oid = list.oid;
    for hsp in list.hsps {
        let (qlen, key) = cull_key(&hsp);
        forest.save(query, qlen, oid, key, hsp);
    }
}

/// `s_BlastHSPCullingFinal` for one forest.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-587
fn finish_forest<H: HitListHsp>(
    forest: CullingForest<H>,
    prelim_hitlist_size: usize,
    out: &mut [Option<BlastpHitList<H>>],
) {
    for (query, culled) in forest.finalize(|query| query, score_compare_hsps::<H>) {
        out[query] = Some(hit_list_of(query, culled, prelim_hitlist_size));
    }
}

/// The hit list of a query's culled lists (the Final of the writer and of the pipe): not cut to
/// `hsplist_max`, not heapified, the best e-value of each list, the worst e-value and the low
/// score of the hit list as `s_BlastHSPCullingFinal` sets them.
fn hit_list_of<H: HitListHsp>(
    query: usize,
    culled: Vec<(u32, Vec<H>)>,
    hsplist_max: usize,
) -> BlastpHitList<H> {
    let mut hit_list = BlastpHitList::new(hsplist_max);
    let mut worst_evalue = 0.0f64;
    let mut low_score = i32::MAX;
    for (oid, hsps) in culled {
        let mut list = BlastpHspList {
            oid,
            query_index: query as u32,
            hsps,
            best_evalue: 0.0,
        };
        update_best_evalue(&mut list);
        worst_evalue = worst_evalue.max(list.best_evalue);
        if let Some(hsp) = list.hsps.first() {
            low_score = low_score.min(hsp.raw_score());
        }
        hit_list.hsplist_array.push(list);
    }
    hit_list.hsplist_count = hit_list.hsplist_array.len();
    hit_list.hsplist_current = hit_list.hsplist_count;
    hit_list.worst_evalue = worst_evalue;
    hit_list.low_score = low_score;
    hit_list
}

/// NCBI's culling pipe on the final hit list of one query, the first step after the traceback
/// (`BlastHSPStreamTBackClose`, before the mask level, `-max_hsps`, `-qcov_hsp_perc` and
/// `-subject_besthit`): the HSP lists sorted by e-value, then every HSP in that order culled
/// into a new tree with the merit `culling_limit` (without the three of the preliminary
/// writer), and the kept HSPs made into the query's hit list again.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:699-733
/// The trees are per query, so the pipe of a batch is this function for each of its queries.
/// The HSPs are those of the redo (or, in NCBI, of the traceback): the offsets are those within
/// the searched range of the query, and `qlen` is the range's length (`query_length`).
/// `hsplist_max` of the new hit list is the old one's: nothing reads it afterwards
/// (`prune_by_size` and the sorts use `hsplist_count`).
pub(super) fn cull_traceback(hit_list: &mut BlastpHitList, culling_limit: i32) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:707-713
    hit_list.hsplist_array.truncate(hit_list.hsplist_count);
    for list in &mut hit_list.hsplist_array {
        sort_hsplist_by_evalue(list);
        if let Some(first) = list.hsps.first() {
            list.best_evalue = first.e_value;
        }
    }
    hit_list.sort_by_evalue();
    let query_index = hit_list
        .hsplist_array
        .first()
        .map_or(0, |list| list.query_index);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:715-725
    let mut forest = CullingForest::new(i64::from(culling_limit));
    for list in std::mem::take(&mut hit_list.hsplist_array) {
        save_list(&mut forest, 0, list);
    }
    let hsplist_max = hit_list.hsplist_max;
    let culled = forest.finalize(|_| 0, score_compare_hsps::<BlastpHsp>);
    let mut new_list = BlastpHitList::new(hsplist_max);
    for (_, lists) in culled {
        new_list = hit_list_of(query_index as usize, lists, hsplist_max);
    }
    *hit_list = new_list;
}

#[cfg(test)]
#[path = "culling_tests.rs"]
mod tests;
