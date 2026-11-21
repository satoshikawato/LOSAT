//! BLASTX `-best_hit_overhang` and `-best_hit_score_edge`: NCBI's best hit writer of the
//! preliminary stage, in place of the HSP collector of each stream (the batch's, or each query
//! chunk's), and its best hit pipe after the traceback (`core::hspfilter_besthit` holds the lists
//! both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//!
//! Nothing here runs without one of the two options: the writer is chosen once per stream
//! (`Runtime::batch_start`, `chunk_start`), the pipe once per batch (`Runtime::batch_end`).

use super::{
    preliminary::compare_score,
    query_setup::QueryContext,
    results::{compare_list, HitList, Hsp, HspList},
};
use crate::core::hspfilter_besthit::{
    translated_query_length, BestHitKey, BestHitLists, BestHitParams,
};

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:87-89
/// The contexts of a BLASTX query: its six frames (`NUM_FRAMES`).
const CONTEXTS_PER_QUERY: usize = 6;

/// What the key of an HSP reads of the query information of its stream (NCBI's `query_info` of
/// the writer or the pipe): each context's length and whether its frame is on the minus strand.
/// A query chunk's stream reads the chunk's own contexts, as NCBI's writer of a chunk reads the
/// chunk's query information.
pub(crate) struct QueryFrames {
    lengths: Vec<i32>,
    minus: Vec<bool>,
}

impl QueryFrames {
    pub(crate) fn of(contexts: &[QueryContext]) -> Self {
        Self {
            lengths: contexts.iter().map(|c| c.length as i32).collect(),
            minus: contexts.iter().map(|c| c.frame < 0).collect(),
        }
    }

    /// The query of an HSP and its key: the start on the plus strand of the query, which for a
    /// minus frame is the query's length less the HSP's end; the length of the HSP's query range;
    /// its raw score and its e-value.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266
    /// The query's length is `BlastQueryInfoGetQueryLength` of a translated query (a length of
    /// nucleotides) and the HSP's offsets are those of its frame (residues): NCBI mixes the two
    /// units, and so does the key.
    fn key(&self, hsp: &Hsp) -> (usize, BestHitKey) {
        let context = hsp.hsp.context;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:40-49
        let query = context / CONTEXTS_PER_QUERY;
        let first = query * CONTEXTS_PER_QUERY;
        let qlen = translated_query_length(&self.lengths[first..first + CONTEXTS_PER_QUERY]);
        (
            query,
            BestHitKey::new(
                self.minus[context],
                qlen,
                hsp.hsp.q_start,
                hsp.hsp.q_end,
                hsp.hsp.score,
                hsp.evalue,
            ),
        )
    }

    /// One subject's HSP list into the lists, in the order of its HSPs.
    fn save(&self, lists: &mut BestHitLists<Hsp>, oid: usize, hsps: Vec<Hsp>) {
        let sid = u32::try_from(oid).expect("a BLASTX subject's OID is an Int4");
        for hsp in hsps {
            let (query, key) = self.key(&hsp);
            lists.save(query, sid, key, hsp);
        }
    }
}

/// `s_BlastHSPBestHitFinal`: the hit list of each of the `num_queries` queries, of `max` lists
/// (`prelim_hitlist_size`): its kept HSPs, one list per subject in the order of the subject's
/// first HSP in the query's list, each sorted by score (`ScoreCompareHSPs`), with its best e-value,
/// the lists sorted by e-value and saved into the hit list in that order (`Blast_HitListUpdate`
/// keeps those of the best e-values when there are more than `max`). A query without kept HSPs has
/// an empty hit list.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231
/// `s_ExportToHitlist` saves the lists into a hit list as large as the HSPs (never full, so
/// `Blast_HitListUpdate` only sets each list's best e-value) and `Blast_HitListSortByEvalue`
/// sorts it (`s_EvalueCompareHSPLists`, `compare_list`: total for the lists of distinct OIDs).
fn final_hit_lists(lists: BestHitLists<Hsp>, num_queries: usize, max: usize) -> Vec<HitList> {
    let mut queries: Vec<HitList> = (0..num_queries).map(|_| HitList::new(max)).collect();
    for (query, kept) in lists.finalize(|a, b| compare_score(&a.hsp, &b.hsp)) {
        let mut exported: Vec<HspList> = kept
            .into_iter()
            .map(|(oid, hsps)| {
                let mut list = HspList {
                    query_index: query,
                    oid: oid as usize,
                    hsps,
                    best_evalue: 0.0,
                };
                list.refresh();
                list
            })
            .collect();
        exported.sort_by(compare_list);
        for list in exported {
            queries[query].update(list);
        }
    }
    queries
}

/// NCBI's best hit writer in place of the HSP collector of a stream (`HspStream::write_observed`):
/// each subject's HSPs, as the stream receives them (the subjects in OID order, each subject's
/// HSPs in the order of its list), go into the list of their query, and `finalize` makes the
/// hit lists of the queries that the collector would have made.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
pub(crate) struct BestHitWriter {
    lists: BestHitLists<Hsp>,
    frames: QueryFrames,
}

impl BestHitWriter {
    /// The writer of a stream whose query information is `frames`, with the preliminary hit list
    /// size `prelim_hitlist_size` (`results::prelim_hitlist_size`, the collector's).
    pub(crate) fn new(
        params: BestHitParams,
        frames: QueryFrames,
        prelim_hitlist_size: usize,
    ) -> Self {
        Self {
            lists: BestHitLists::writer(params, prelim_hitlist_size),
            frames,
        }
    }

    /// `s_BlastHSPBestHitRun` for the HSPs of subject `oid`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-352
    pub(crate) fn write(&mut self, oid: usize, hsps: Vec<Hsp>) {
        self.frames.save(&mut self.lists, oid, hsps);
    }

    /// `s_BlastHSPBestHitFinal`: the hit lists of the stream's `num_queries` queries, of `max`
    /// lists each (`final_hit_lists`).
    pub(crate) fn finalize(self, num_queries: usize, max: usize) -> Vec<HitList> {
        final_hit_lists(self.lists, num_queries, max)
    }
}

/// NCBI's best hit pipe on the queries' final hit lists, the first step after the traceback
/// (`BlastHSPStreamTBackClose`, before `-max_hsps`, `-qcov_hsp_perc`, `-subject_besthit`, the
/// e-value sort and the cut to the hit list size): each query's lists sorted by their stored best
/// e-values (the HSPs of each list keep their order), every HSP in that order put into a new list
/// with the overhang and the score edge of the options, and the kept HSPs made into the queries'
/// hit lists again, of `prelim_hitlist_size` lists. `frames` is the batch's query information.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539
/// The HSPs are those of the traceback, of the redo of `-comp_based_stats` or of the ungapped
/// stage, in the order of their lists in the hit lists (`Blast_HitListUpdate` keeps a list's
/// order, or sorts it by e-value once the hit list is full), with the best e-value that it set.
pub(crate) fn traceback_pipe(
    queries: &mut Vec<HitList>,
    frames: &QueryFrames,
    params: BestHitParams,
    prelim_hitlist_size: usize,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
    // `Blast_HSPResultsSortByEvalue`: each query's lists by their best e-values; a list without
    // HSPs (there is none after the traceback) would go last and be removed.
    for query in queries.iter_mut() {
        query.lists.sort_by(compare_list);
        query.lists.retain(|list| !list.hsps.is_empty());
    }
    let mut lists = BestHitLists::pipe(params, prelim_hitlist_size);
    for query in queries.iter_mut() {
        for list in std::mem::take(&mut query.lists) {
            frames.save(&mut lists, list.oid, list.hsps);
        }
    }
    *queries = final_hit_lists(lists, queries.len(), prelim_hitlist_size);
}

#[cfg(test)]
#[path = "besthit_tests.rs"]
mod tests;
