//! BLASTP `-culling_limit`: NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback (`core::hspfilter_culling` holds the interval tree both share).
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
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
/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:612-633
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
/// There are two forests, because the preliminary lists have two representations
/// (EXPERIMENT LOSAT_X_BPPCOMPACT, `x_prelim_compact.rs`): the compact lists (the default),
/// the full ones (`LOSAT_X_DEFAULT=0`), or both (the shadow mode). Each forest receives what the
/// subject result holds of its kind and gives the same culled lists as the other.
pub(super) struct PrelimCulling {
    full: CullingForest<BlastpHsp>,
    compact: CullingForest<XPrelimHsp>,
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
    pub(super) fn new(prelim_hitlist_size: usize, culling_limit: i32) -> Self {
        let culling_max = CullingForest::<BlastpHsp>::prelim_max(culling_limit);
        Self {
            full: CullingForest::new(culling_max),
            compact: CullingForest::new(culling_max),
            prelim_hitlist_size,
        }
    }

    /// `s_BlastHSPCullingRun` for the HSPs of one subject (all queries of the call).
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
    /// The tree id is the query of the call (blastp has one context per query; in a chunk call
    /// the parts of the chunk are its queries). An HSP that is not saved is freed. The subject
    /// result is emptied, as the collector empties it.
    pub(super) fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        for (query, list) in subject_result.query_hsp_lists.drain(..) {
            save_list(&mut self.full, query, list);
        }
        for (query, list) in subject_result.x_compact_hsp_lists.drain(..) {
            save_list(&mut self.compact, query, list);
        }
    }

    /// `s_BlastHSPCullingFinal` for both forests: the hit list of each query with a tree goes
    /// into `lists` (full HSPs) or `compact_lists`. It has one HSP list for each subject in the
    /// order of its first HSP in the rip of the tree, each sorted by score, and as many lists as
    /// the tree keeps: the hit list is made with `prelim_hitlist_size` as its maximum but is not
    /// cut to it, and it is not heapified.
    ///
    /// The merge of a split query (`query_split.rs`, `Blast_HitListMerge`) reads `hsplist_max`
    /// and the lists; the first chunk's hit list is moved as it is, with more lists than
    /// `hsplist_max`. The stream close (`sort_by_evalue`) sorts and purges the first
    /// `hsplist_count` lists.
    pub(super) fn finish(
        self,
        lists: &mut [Option<BlastpHitList>],
        compact_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) {
        finish_forest(self.full, self.prelim_hitlist_size, lists);
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
/// The HSPs are those of the redo (or, in NCBI, of the traceback): the offsets are those within
/// the searched range of the query, and `qlen` is the range's length (`query_length`).
/// `hsplist_max` of the new hit list is the old one's: nothing reads it afterwards
/// (`prune_by_size` and the sorts use `hsplist_count`).
pub(super) fn cull_traceback(hit_list: &mut BlastpHitList, culling_limit: i32) {
    // NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:707-713
    // ```c
    //          for (sid = 0; sid < num_list; ++sid) {
    //         	 hsp_list = results->hitlist_array[qid]->hsplist_array[sid];
    //         	 Blast_HSPListSortByEvalue(hsp_list);
    //         	 hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
    //          }
    //          Blast_HitListSortByEvalue(results->hitlist_array[qid]);
    // ```
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

    // NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:715-725
    // ```c
    //       for (sid = 0; sid < num_list; ++sid) {
    //          s_BlastHSPCullingRun(data,
    //                    results->hitlist_array[qid]->hsplist_array[sid]);
    // ```
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
mod tests {
    use super::*;
    use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;

    fn prelim(
        query: u32,
        qlen: usize,
        qs: i32,
        qe: i32,
        ss: i32,
        score: i32,
    ) -> BlastpPreliminaryHsp {
        BlastpPreliminaryHsp {
            query_start: qs,
            query_end: qe,
            subject_start: ss,
            subject_end: ss + (qe - qs),
            gapped_query_start: qs,
            gapped_subject_start: ss,
            raw_score: score,
            query_context: query as i32,
            query_frame: 0,
            subject_frame: 0,
            query_length: qlen,
            q_idx: query,
            s_idx: 0,
        }
    }

    fn full_hsp(hsp: &BlastpPreliminaryHsp, oid: u32, e_value: f64) -> BlastpHsp {
        let mut full = super::super::blastp_hsp_from_preliminary(hsp, e_value);
        full.s_idx = oid;
        full
    }

    /// The subject result of `oid` with the HSPs of one query in both representations.
    fn subject(oid: usize, query: usize, hsps: &[(BlastpPreliminaryHsp, f64)]) -> BlastpSubjectResult {
        let mut result = BlastpSubjectResult::new(oid);
        result.query_hsp_lists.push((
            query,
            BlastpHspList {
                oid: oid as u32,
                query_index: query as u32,
                hsps: hsps.iter().map(|(h, e)| full_hsp(h, oid as u32, *e)).collect(),
                best_evalue: 0.0,
            },
        ));
        result.x_compact_hsp_lists.push((
            query,
            BlastpHspList {
                oid: oid as u32,
                query_index: query as u32,
                hsps: hsps
                    .iter()
                    .map(|(h, e)| XPrelimHsp {
                        hsp: BlastpPreliminaryHsp {
                            s_idx: oid as u32,
                            ..*h
                        },
                        e_value: *e,
                    })
                    .collect(),
                best_evalue: 0.0,
            },
        ));
        result
    }

    fn finished(
        culling: PrelimCulling,
        queries: usize,
    ) -> (Vec<Option<BlastpHitList>>, Vec<Option<BlastpHitList<XPrelimHsp>>>) {
        let mut lists: Vec<Option<BlastpHitList>> = (0..queries).map(|_| None).collect();
        let mut compact: Vec<Option<BlastpHitList<XPrelimHsp>>> =
            (0..queries).map(|_| None).collect();
        culling.finish(&mut lists, &mut compact);
        (lists, compact)
    }

    fn scores(lists: &[Option<BlastpHitList>], query: usize) -> Vec<(u32, Vec<i32>)> {
        lists[query]
            .as_ref()
            .map(|hit_list| {
                hit_list
                    .hsplist_array
                    .iter()
                    .map(|l| (l.oid, l.hsps.iter().map(|h| h.raw_score).collect()))
                    .collect()
            })
            .unwrap_or_default()
    }

    fn compact_scores(lists: &[Option<BlastpHitList<XPrelimHsp>>], query: usize) -> Vec<(u32, Vec<i32>)> {
        lists[query]
            .as_ref()
            .map(|hit_list| {
                hit_list
                    .hsplist_array
                    .iter()
                    .map(|l| (l.oid, l.hsps.iter().map(|h| h.hsp.raw_score).collect()))
                    .collect()
            })
            .unwrap_or_default()
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:612-633: the key of a protein HSP is its
    // own query range (0-based start, exclusive end), its score and its 0-based subject offset;
    // `qlen` is the HSP's query length.
    #[test]
    fn the_key_is_the_context_relative_range_score_and_subject_offset() {
        let hsp = full_hsp(&prelim(0, 321, 10, 60, 7, 55), 3, 1e-5);
        let (qlen, key) = cull_key(&hsp);
        assert_eq!(qlen, 321);
        assert_eq!(
            key,
            CullKey {
                begin: 10,
                end: 60,
                score: 55,
                subject_offset: 7
            }
        );
        let compact = XPrelimHsp {
            hsp: prelim(0, 321, 10, 60, 7, 55),
            e_value: 1e-5,
        };
        assert_eq!(cull_key(&compact), (qlen, key));
    }

    // NCBI hspfilter_culling.c:516-581 and blast_hits.c:3125-3133: the hit list is made with
    // `prelim_hitlist_size` and holds more lists than that, in the order of the first HSP of
    // each subject in the rip, with the best e-value, worst e-value and low score set, and not
    // heapified; and the compact forest gives the same lists as the full one.
    #[test]
    fn the_hit_list_is_not_cut_and_keeps_the_rip_order_in_both_forests() {
        let mut culling = PrelimCulling::new(2, 5);
        for oid in 0..4usize {
            let a = prelim(0, 1000, 100 + 50 * oid as i32, 120 + 50 * oid as i32, 1, 40 - oid as i32);
            let b = prelim(0, 1000, 100 + 50 * oid as i32, 110 + 50 * oid as i32, 2, 20);
            let c = prelim(1, 500, 10, 30, 1, 30);
            let mut first = subject(oid, 0, &[(a, 1e-5), (b, 1e-2)]);
            let second = subject(oid, 1, &[(c, 1e-3)]);
            first.query_hsp_lists.extend(second.query_hsp_lists);
            first.x_compact_hsp_lists.extend(second.x_compact_hsp_lists);
            culling.add(&mut first);
            assert!(first.query_hsp_lists.is_empty() && first.x_compact_hsp_lists.is_empty());
        }
        let (lists, compact) = finished(culling, 2);
        let q0 = lists[0].as_ref().unwrap();
        assert_eq!(q0.hsplist_max, 2);
        assert_eq!(q0.hsplist_count, 4);
        assert!(!q0.heapified);
        assert_eq!(q0.low_score, 37);
        assert!((q0.worst_evalue - 1e-5).abs() < 1e-30);
        // Limit 5 (+3 = 8) keeps everything. The root list holds the newest HSP first, so the
        // rip meets the last subject first; every list is score sorted.
        let expected = vec![
            (3, vec![37, 20]),
            (2, vec![38, 20]),
            (1, vec![39, 20]),
            (0, vec![40, 20]),
        ];
        assert_eq!(scores(&lists, 0), expected);
        assert_eq!(compact_scores(&compact, 0), expected);
        assert_eq!(lists[1].as_ref().unwrap().hsplist_count, 4);
        assert_eq!(compact[1].as_ref().unwrap().hsplist_count, 4);
        let q0c = compact[0].as_ref().unwrap();
        assert_eq!(q0c.low_score, q0.low_score);
        assert_eq!(q0c.worst_evalue.to_bits(), q0.worst_evalue.to_bits());
    }

    // NCBI hspfilter_culling.c:100-112: the lower OID wins a tie, and the subjects must come in
    // OID order for the later one to lose.
    #[test]
    fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
        let mut culling = PrelimCulling::new(10, 1);
        let h = prelim(0, 1000, 100, 200, 5, 80);
        culling.add(&mut subject(0, 0, &[(h, 1e-10)]));
        culling.add(&mut subject(1, 0, &[(h, 1e-10)]));
        let (lists, compact) = finished(culling, 1);
        assert_eq!(scores(&lists, 0), vec![(0, vec![80])]);
        assert_eq!(compact_scores(&compact, 0), vec![(0, vec![80])]);
    }

    // NCBI hspfilter_culling.c:602-644: the tree of a query is made with the length of its HSPs,
    // and a chunk call's parts are queries of their own: a better HSP of query 1 does not
    // dominate an HSP of query 0 over the same range.
    #[test]
    fn each_query_has_its_own_tree() {
        let mut culling = PrelimCulling::new(10, 1);
        culling.add(&mut subject(0, 0, &[(prelim(0, 400, 10, 90, 0, 50), 1e-3)]));
        culling.add(&mut subject(1, 1, &[(prelim(1, 300, 10, 90, 0, 90), 1e-9)]));
        culling.add(&mut subject(2, 0, &[(prelim(0, 400, 10, 90, 0, 60), 1e-5)]));
        let (lists, _) = finished(culling, 2);
        // Limit 1: score 60 dominates score 50 in query 0.
        assert_eq!(scores(&lists, 0), vec![(2, vec![60])]);
        assert_eq!(scores(&lists, 1), vec![(1, vec![90])]);
    }

    fn pipe_hsp(oid: u32, qs: usize, qe: usize, score: i32, e_value: f64) -> BlastpHsp {
        let hsp = prelim(0, 1000, qs as i32 - 1, qe as i32, 0, score);
        full_hsp(&hsp, oid, e_value)
    }

    // NCBI hspfilter_culling.c:699-733: the lists are sorted by e-value (best first, then the
    // first score, the higher OID first), the HSPs are saved in that order into a tree with the
    // limit itself (no +3), and the kept HSPs are made into lists again, each sorted by score.
    #[test]
    fn the_pipe_saves_in_e_value_order_with_the_plain_limit() {
        let mut hit_list = BlastpHitList::new(10);
        // OID 1 has the better e-value, so it is saved first; OID 0's weaker HSP over the same
        // range is dominated by it (limit 1).
        for (oid, score, e_value) in [(0u32, 60, 1e-6), (1, 90, 1e-12)] {
            hit_list.hsplist_array.push(BlastpHspList {
                oid,
                query_index: 0,
                hsps: vec![pipe_hsp(oid, 11, 90, score, e_value)],
                best_evalue: e_value,
            });
        }
        hit_list.hsplist_count = 2;
        cull_traceback(&mut hit_list, 1);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 1);
        assert_eq!(hit_list.hsplist_max, 10);
        assert!(!hit_list.heapified);
        assert_eq!(hit_list.low_score, 90);
    }

    // With the limit 2 the second dominated HSP stays; the lists keep the rip order and each is
    // sorted by score.
    #[test]
    fn the_pipe_keeps_the_hsps_that_fewer_than_the_limit_dominate() {
        let mut hit_list = BlastpHitList::new(10);
        for (oid, score, e_value) in [(0u32, 60, 1e-6), (1, 90, 1e-12), (2, 70, 1e-8)] {
            hit_list.hsplist_array.push(BlastpHspList {
                oid,
                query_index: 0,
                hsps: vec![pipe_hsp(oid, 11, 90, score, e_value)],
                best_evalue: e_value,
            });
        }
        hit_list.hsplist_count = 3;
        cull_traceback(&mut hit_list, 2);
        let mut kept: Vec<u32> = hit_list.hsplist_array.iter().map(|l| l.oid).collect();
        kept.sort_unstable();
        // Saved in the order 1 (1e-12), 2 (1e-8), 0 (1e-6): OID 0 is dominated by both.
        assert_eq!(kept, vec![1, 2]);
        assert_eq!(hit_list.hsplist_count, 2);
    }
}
