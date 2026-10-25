//! BLASTN `-culling_limit`: NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback (`core::hspfilter_culling` holds the interval tree both share).
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
//!
//! This module is a child of `run` for the preliminary HSP types that `run` keeps private.
//! Nothing here runs without `-culling_limit` (a limit of 0 is NCBI's `SetCullingLimit(0)`,
//! which is ignored): the choice is made once per search, where `run` picks the sink of the
//! preliminary hit lists, and once per query where it finishes a final hit list.

use super::{score_compare_prelim_hits, PrelimHit, PrelimHitLists, PrelimHspList};
use crate::algorithm::blastn::hsp::{
    score_compare_hsps, sort_hsplist_by_evalue, update_best_evalue, BlastnHitList, BlastnHspList,
    HitList, HitListEntry,
};
use crate::core::hspfilter_culling::{CullKey, CullingForest};

/// An `Int4` length as NCBI's `query_length` holds it (the lengths of LOSAT's queries fit).
fn int4(length: usize) -> i32 {
    i32::try_from(length).unwrap_or(i32::MAX)
}

/// The key of an HSP for the tree: the range on the plus strand of the query. A minus
/// context folds its range over the query length.
///
/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:617-633
/// ```c
///       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
///       A.sid   = hsp_list->oid;
///       A.merit = params->culling_max;
///       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
///       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
///     	  A.begin = qlen - A.hsp->query.end;
///     	  A.end   = qlen -  A.hsp->query.offset;
///       }
///       else {
///     	  A.begin = A.hsp->query.offset;
///     	  A.end   = A.hsp->query.end;
///       }
/// ```
/// `query_offset` and `query_end` are the HSP's offsets within its context.
fn cull_key(
    minus: bool,
    qlen: i32,
    query_offset: usize,
    query_end: usize,
    score: i32,
    subject_offset: usize,
) -> CullKey {
    let (query_offset, query_end) = (query_offset as i32, query_end as i32);
    let (begin, end) = if minus {
        (qlen - query_end, qlen - query_offset)
    } else {
        (query_offset, query_end)
    };
    CullKey {
        begin,
        end,
        score,
        subject_offset: subject_offset as i64,
    }
}

/// NCBI's culling writer in place of the HSP collector (`PrelimHitCollector`): each subject's
/// preliminary HSPs go into the tree of their query, and `finish` makes the hit lists of the
/// queries. The subjects must come in OID order and the HSPs of a subject in the order of its
/// score-sorted HSP list (the order in which the collector receives them), because the
/// result depends on the order of the writes.
pub(super) struct PrelimCulling {
    forest: CullingForest<PrelimHit>,
    /// The length of each query (of each query part in the chunk stage of a split batch),
    /// `contexts[context].query_length` of both its contexts.
    query_lengths: Vec<i32>,
    prelim_hitlist_size: usize,
}

impl PrelimCulling {
    /// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:318-342
    /// ```c
    ///         else if (filt_opts->culling_opts &&
    ///                  (filt_opts->culling_stage & ePrelimSearch))
    ///         {
    ///             BlastHSPCullingParams* params =
    ///                 BlastHSPCullingParamsNew(opts_memento->m_HitSaveOpts,
    ///                      filt_opts->culling_opts,
    ///                      opts_memento->m_ExtnOpts->compositionBasedStats,
    ///                      opts_memento->m_ScoringOpts->gapped_calculation);
    ///             if(params->culling_max > 1){
    ///             	params->culling_max += 3;
    ///             }
    ///             writer_info = BlastHSPCullingInfoNew(params);
    ///         }
    /// ```
    pub(super) fn new(
        query_lengths: &[usize],
        prelim_hitlist_size: usize,
        culling_limit: i32,
    ) -> Self {
        Self {
            forest: CullingForest::new(CullingForest::<PrelimHit>::prelim_max(culling_limit)),
            query_lengths: query_lengths.iter().map(|&length| int4(length)).collect(),
            prelim_hitlist_size,
        }
    }

    /// `s_BlastHSPCullingRun` for the HSPs of subject `oid` (all queries, both strands).
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// ```c
    ///    for (i=0; i<hsp_list->hspcnt; ++i) {
    ///       /* wrap the hsp with a LinkedHSP structure */
    ///       A.hsp   = hsp_list->hsp_array[i];
    ///       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
    ///       A.sid   = hsp_list->oid;
    ///       A.merit = params->culling_max;
    /// ...
    ///       if (! c_tree[A.cid]) {
    ///          c_tree[A.cid] = s_CTreeNew(qlen);
    ///       }
    ///
    ///       if(s_SaveHSP(c_tree[A.cid], &A)){
    ///     	 hsp_list->hsp_array[i] = NULL;
    ///       }
    ///    }
    /// ```
    /// The tree id is the query (`cid` is its first context, and the trees are ripped in the
    /// order of their `cid`). An HSP that is not saved is freed.
    pub(super) fn add(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        for hit in hits {
            let query = hit.query_idx as usize;
            let qlen = self.query_lengths[query];
            let key = cull_key(
                hit.context_idx % 2 == 1,
                qlen,
                hit.prelim_qs,
                hit.prelim_qe,
                hit.prelim_score,
                hit.prelim_ss,
            );
            self.forest.save(query, qlen, oid as u32, key, hit);
        }
    }

    /// `s_BlastHSPCullingFinal`: the hit list of each query with a tree. It has one HSP
    /// list for each subject in the order of its first HSP in the rip of the tree, each
    /// sorted by score, and as many lists as the tree keeps: the hit list is made with
    /// `prelim_hitlist_size` as its maximum but is not cut to it, and it is not heapified.
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:516-581
    /// ```c
    ///          if (!results->hitlist_array[qid]) {
    ///             results->hitlist_array[qid] = Blast_HitListNew(params->prelim_hitlist_size);
    ///          }
    /// ...
    ///          /* sort hsplist */
    ///          worst_evalue = 0.0;
    ///          low_score = INT4_MAX;
    ///          for (sid=0; sid < hitlist->hsplist_count; ++sid) {
    ///             list = hitlist->hsplist_array[sid];
    ///             best_evalue = (double) INT4_MAX;
    ///             for (id=0; id < list->hspcnt; ++id) {
    ///                 best_evalue = MIN(list->hsp_array[id]->evalue, best_evalue);
    ///             }
    ///             Blast_HSPListSortByScore(list);
    ///             list->best_evalue = best_evalue;
    ///             worst_evalue = MAX(worst_evalue, best_evalue);
    ///             low_score = MIN(list->hsp_array[0]->score, low_score);
    ///          }
    ///          hitlist->worst_evalue = worst_evalue;
    ///          hitlist->low_score = low_score;
    /// ```
    /// The merge of a split query (`merge_prelim_hit_list`, `Blast_HitListMerge`) and the
    /// traceback input (`prelim_hits_by_subject`, `listed_prelim_hits`) read a hit list
    /// through its lists only; the first chunk's list is moved as it is, with more lists
    /// than `hsplist_max`.
    pub(super) fn finish(self) -> PrelimHitLists {
        let mut hit_lists: PrelimHitLists = Vec::with_capacity(self.query_lengths.len());
        hit_lists.resize_with(self.query_lengths.len(), || None);
        let culled = self
            .forest
            .finalize(|query| query, score_compare_prelim_hits);
        for (query, lists) in culled {
            let mut hit_list = HitList::new(self.prelim_hitlist_size);
            let mut worst_evalue = 0.0f64;
            let mut low_score = i32::MAX;
            for (oid, hsps) in lists {
                let mut list = PrelimHspList {
                    oid,
                    hsps,
                    best_evalue: 0.0,
                };
                list.update_best_evalue();
                worst_evalue = worst_evalue.max(list.best_evalue);
                if let Some(score) = list.first_score() {
                    low_score = low_score.min(score);
                }
                hit_list.hsplist_array.push(list);
            }
            hit_list.hsplist_count = hit_list.hsplist_array.len();
            hit_list.hsplist_current = hit_list.hsplist_count;
            hit_list.worst_evalue = worst_evalue;
            hit_list.low_score = low_score;
            hit_lists[query] = Some(hit_list);
        }
        hit_lists
    }
}

/// NCBI's culling pipe on the final hit list of one query, the first step after the traceback
/// (`BlastHSPStreamTBackClose`, before the mask level, `-max_hsps`, `-qcov_hsp_perc` and
/// `-subject_besthit`): the HSP lists sorted by e-value, then every HSP in that order culled
/// into a new tree with the merit `culling_limit` (without the three of the preliminary
/// writer), and the kept HSPs made into the query's hit list again.
///
/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:699-733
/// ```c
/// static int
/// s_BlastHSPCullingPipeRun(void* data, BlastHSPResults* results)
/// {
///    int qid, sid, num_list;
///    BlastHSPList* hsp_list = NULL;
///    s_BlastHSPCullingInit(data, results);
///    for (qid = 0; qid < results->num_queries; ++qid) {
///          if (!(results->hitlist_array[qid])) continue;
///          num_list = results->hitlist_array[qid]->hsplist_count;
///          for (sid = 0; sid < num_list; ++sid) {
///         	 hsp_list = results->hitlist_array[qid]->hsplist_array[sid];
///         	 Blast_HSPListSortByEvalue(hsp_list);
///         	 hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
///          }
///          Blast_HitListSortByEvalue(results->hitlist_array[qid]);
///    }
///
///    for (qid = 0; qid < results->num_queries; ++qid) {
///       if (!(results->hitlist_array[qid])) continue;
///       num_list = results->hitlist_array[qid]->hsplist_count;
///       for (sid = 0; sid < num_list; ++sid) {
///          s_BlastHSPCullingRun(data,
///                    results->hitlist_array[qid]->hsplist_array[sid]);
///          results->hitlist_array[qid]->hsplist_array[sid] = NULL;
///       }
///       results->hitlist_array[qid]->hsplist_count = 0;
///       Blast_HitListFree(results->hitlist_array[qid]);
///       results->hitlist_array[qid] = NULL;
///    }
///    s_BlastHSPCullingFinal(data, results);
///    return 0;
/// }
/// ```
/// The trees are per query, so the pipe of a batch is this function for each of its queries.
/// `query_len` is the length of each context of the query (`finish_query_hit_list`). The HSP
/// offsets are those within the HSP's context and the context is minus for a negative frame
/// (`BLAST_ContextToFrame`).
pub(super) fn cull_traceback(
    hit_list: &mut BlastnHitList,
    query_len: usize,
    culling_limit: i32,
    prelim_hitlist_size: usize,
) {
    // NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:707-713
    // ```c
    //          for (sid = 0; sid < num_list; ++sid) {
    //         	 hsp_list = results->hitlist_array[qid]->hsplist_array[sid];
    //         	 Blast_HSPListSortByEvalue(hsp_list);
    //         	 hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
    //          }
    //          Blast_HitListSortByEvalue(results->hitlist_array[qid]);
    // ```
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

    // NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:715-725
    // ```c
    //       for (sid = 0; sid < num_list; ++sid) {
    //          s_BlastHSPCullingRun(data,
    //                    results->hitlist_array[qid]->hsplist_array[sid]);
    // ```
    let qlen = int4(query_len);
    let mut forest = CullingForest::new(i64::from(culling_limit));
    for list in std::mem::take(&mut hit_list.hsplist_array) {
        let oid = list.oid;
        for hsp in list.hsps {
            let key = cull_key(
                hsp.query_frame < 0,
                qlen,
                hsp.internal_q_offset_0,
                hsp.internal_q_end_0,
                hsp.raw_score,
                hsp.internal_s_offset_0,
            );
            forest.save(0, qlen, oid, key, hsp);
        }
    }

    // `s_BlastHSPCullingFinal`: the new hit list (as `PrelimCulling::finish`, with
    // `ScoreCompareHSPs` on the final HSPs).
    // `hsplist_max` is set from `prelim_hitlist_size` as NCBI's pipe Final sets it
    // (`Blast_HitListNew(params->prelim_hitlist_size)`); nothing reads it afterwards: the hit
    // list is finished by `finish_query_hit_list`, which trims by `hitlist_size`, and neither
    // `HitList::update` nor the array growth runs on it.
    let mut new_list = BlastnHitList::new(prelim_hitlist_size);
    let mut worst_evalue = 0.0f64;
    let mut low_score = i32::MAX;
    let culled = forest.finalize(|_| 0, score_compare_hsps);
    for (_, lists) in culled {
        for (oid, hsps) in lists {
            let mut list = BlastnHspList {
                oid,
                query_index,
                hsps,
                best_evalue: 0.0,
            };
            update_best_evalue(&mut list);
            worst_evalue = worst_evalue.max(list.best_evalue);
            if let Some(first) = list.hsps.first() {
                low_score = low_score.min(first.raw_score);
            }
            new_list.hsplist_array.push(list);
        }
    }
    new_list.hsplist_count = new_list.hsplist_array.len();
    new_list.hsplist_current = new_list.hsplist_count;
    new_list.worst_evalue = worst_evalue;
    new_list.low_score = low_score;
    *hit_list = new_list;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prelim(
        query: u32,
        minus: bool,
        qs: usize,
        qe: usize,
        ss: usize,
        score: i32,
        evalue: f64,
    ) -> PrelimHit {
        PrelimHit {
            context_idx: 2 * query + u32::from(minus),
            query_idx: query,
            query_frame: if minus { -1 } else { 1 },
            query_context_offset: 0,
            prelim_qs: qs,
            prelim_qe: qe,
            prelim_ss: ss,
            prelim_se: ss + (qe - qs),
            prelim_score: score,
            seed_qs: qs,
            seed_ss: ss,
            prelim_evalue: evalue,
        }
    }

    fn scores(lists: &PrelimHitLists, query: usize) -> Vec<(u32, Vec<i32>)> {
        lists[query]
            .as_ref()
            .map(|hit_list| {
                hit_list
                    .hsplist_array
                    .iter()
                    .map(|list| (list.oid, list.hsps.iter().map(|h| h.prelim_score).collect()))
                    .collect()
            })
            .unwrap_or_default()
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:617-633: a minus HSP is folded to
    // the plus strand, so a minus HSP over the same plus range as a better plus HSP is
    // dominated by it, and one elsewhere is not.
    #[test]
    fn a_minus_hsp_is_folded_to_the_plus_strand() {
        let mut culling = PrelimCulling::new(&[1000], 10, 1);
        // Minus (700, 800) is plus (200, 300) of a query of 1000.
        culling.add(
            0,
            vec![
                prelim(0, false, 200, 300, 5, 100, 1e-20),
                prelim(0, true, 700, 800, 5, 60, 1e-10),
                prelim(0, true, 200, 300, 9, 50, 1e-9),
            ],
        );
        let lists = culling.finish();
        assert_eq!(scores(&lists, 0), vec![(0, vec![100, 50])]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:516-581, blast_hits.c:3125-3133: the
    // hit list is made with `prelim_hitlist_size` and holds more lists than that, in the
    // order of the first HSP of each subject, with the best e-value, worst e-value and low
    // score set.
    #[test]
    fn the_hit_list_is_not_cut_and_keeps_the_rip_order() {
        let mut culling = PrelimCulling::new(&[1000, 500], 2, 5);
        for oid in 0..4usize {
            culling.add(
                oid,
                vec![
                    prelim(
                        0,
                        false,
                        100 + 50 * oid,
                        120 + 50 * oid,
                        1,
                        40 - oid as i32,
                        1e-5,
                    ),
                    prelim(1, false, 10, 30, 1, 30, 1e-3),
                    prelim(0, false, 100 + 50 * oid, 110 + 50 * oid, 2, 20, 1e-2),
                ],
            );
        }
        let lists = culling.finish();
        let q0 = lists[0].as_ref().unwrap();
        assert_eq!(q0.hsplist_max, 2);
        assert_eq!(q0.hsplist_count, 4);
        assert!(!q0.heapified);
        assert_eq!(q0.low_score, 37);
        assert!((q0.worst_evalue - 1e-5).abs() < 1e-30);
        // Limit 5 (+3 = 8) keeps everything. The root list holds the newest HSP first, so
        // the rip meets the last subject first; every list is score sorted.
        assert_eq!(
            scores(&lists, 0),
            vec![
                (3, vec![37, 20]),
                (2, vec![38, 20]),
                (1, vec![39, 20]),
                (0, vec![40, 20])
            ]
        );
        assert_eq!(lists[1].as_ref().unwrap().hsplist_count, 4);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:100-112: the lower OID wins a tie,
    // and the subject must come in OID order for the later one to lose.
    #[test]
    fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
        let mut culling = PrelimCulling::new(&[1000], 10, 1);
        culling.add(0, vec![prelim(0, false, 100, 200, 5, 80, 1e-10)]);
        culling.add(1, vec![prelim(0, false, 100, 200, 5, 80, 1e-10)]);
        let lists = culling.finish();
        assert_eq!(scores(&lists, 0), vec![(0, vec![80])]);
    }
}
