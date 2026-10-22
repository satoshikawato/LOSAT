//! TBLASTN `-culling_limit`: NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback (`core::hspfilter_culling` holds the interval tree both share).
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
//!
//! This module is a child of `stage_d_results` for the private fields of
//! `KappaResultHitList`. Nothing here runs without `-culling_limit` (a limit of 0 is NCBI's
//! `SetCullingLimit(0)`, which is ignored): the choice is made once per call of
//! `run_preliminary_subjects`, where it picks the sink of the preliminary hit lists, and once
//! per query batch where the search finishes the final hit lists.
//!
//! TBLASTN has one context for each protein query, and the six frames of a subject are not
//! contexts: `s_DominateTest` reads the score, the query range, the OID and
//! `hsp->subject.offset` (the offset in the frame's own translation) and never the frame,
//! so the HSPs of all six frames of a subject go into the tree of the query.

use anyhow::{ensure, Result};

use super::super::stage_d_linking::{score_compare, LinkedHsp, LinkedHspList};
use super::{
    compare_result_lists, refresh_best_evalue, sort_hsps_by_evalue, KappaHspPayload,
    KappaResultHitList, KappaResultList,
};
use crate::core::hspfilter_culling::{CullKey, CullingForest};

/// The key of an HSP for the tree: the protein range of the query and the frame-local offset
/// of the subject. NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:612-633 (the
/// non-blastn branch: `begin = query.offset`, `end = query.end`, no strand fold) and
/// c++/src/algo/blast/core/hspfilter_culling.c:79-120 (`hsp->score`, `hsp->subject.offset`).
fn cull_key(linked: &LinkedHsp) -> CullKey {
    CullKey {
        begin: linked.hsp.q_start,
        end: linked.hsp.q_end,
        score: linked.hsp.score,
        subject_offset: i64::from(linked.hsp.s_start),
    }
}

/// The hit list that `s_BlastHSPCullingFinal` makes of the lists the trees keep for one query.
///
/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:516-581
/// ```c
///          if (!results->hitlist_array[qid]) {
///             results->hitlist_array[qid] = Blast_HitListNew(params->prelim_hitlist_size);
///          }
/// ...
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
/// The hit list is made with `max` as its maximum but is not cut to it and is not
/// heapified (`Blast_HitListNew` does not heapify); the merge of a split query's chunks
/// (`hit_list_merge`) and the lists read through `lists()` use only the lists.
fn hit_list_of_culled(
    lists: Vec<(u32, Vec<LinkedHsp>)>,
    payloads: Vec<Vec<Option<KappaHspPayload>>>,
    max: usize,
) -> KappaResultHitList {
    debug_assert_eq!(lists.len(), payloads.len());
    let mut worst_evalue = 0.0f64;
    let mut low_score = i32::MAX;
    let mut result = Vec::with_capacity(lists.len());
    for ((oid, hsps), payloads) in lists.into_iter().zip(payloads) {
        let payloads: Vec<KappaHspPayload> = payloads.into_iter().flatten().collect();
        debug_assert!(payloads.is_empty() || payloads.len() == hsps.len());
        let mut list = KappaResultList {
            oid: oid as i32,
            hsps: LinkedHspList {
                hsps,
                best_evalue: 0.0,
            },
            payloads,
        };
        refresh_best_evalue(&mut list.hsps);
        worst_evalue = worst_evalue.max(list.hsps.best_evalue);
        low_score = low_score.min(list.hsps.hsps[0].hsp.score);
        result.push(list);
    }
    KappaResultHitList {
        lists: result,
        max,
        worst_evalue,
        low_score,
        heapified: false,
    }
}

/// NCBI's culling writer in place of the hit-list collector: each subject's preliminary HSPs
/// go into the tree of their query, and `finish` makes the hit lists of the queries. The
/// subjects must come in OID order and the HSPs of a subject in the order of its HSP list
/// (score order), because the result depends on the order of the writes.
pub(in crate::algorithm::tblastn) struct PrelimCulling {
    forest: CullingForest<LinkedHsp>,
    /// The length of each query (of each query part in the chunk stage of a split batch),
    /// `contexts[context].query_length`.
    query_lengths: Vec<i32>,
    prelim_size: usize,
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
    pub(in crate::algorithm::tblastn) fn new(
        query_lengths: &[i32],
        prelim_size: usize,
        culling_limit: i32,
    ) -> Self {
        Self {
            forest: CullingForest::new(CullingForest::<LinkedHsp>::prelim_max(culling_limit)),
            query_lengths: query_lengths.to_vec(),
            prelim_size,
        }
    }

    /// `s_BlastHSPCullingRun` for the HSPs of subject `oid`: all queries and all frames, in
    /// the order of the subject's list.
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// ```c
    ///    for (i=0; i<hsp_list->hspcnt; ++i) {
    ///       /* wrap the hsp with a LinkedHSP structure */
    ///       A.hsp   = hsp_list->hsp_array[i];
    ///       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
    ///       A.sid   = hsp_list->oid;
    ///       A.merit = params->culling_max;
    ///       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
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
    /// An HSP that is not saved is freed.
    pub(in crate::algorithm::tblastn) fn add(
        &mut self,
        oid: usize,
        hsps: &[LinkedHsp],
    ) -> Result<()> {
        let sid = u32::try_from(oid)?;
        for linked in hsps {
            let qlen = self.query_lengths[linked.context];
            self.forest
                .save(linked.context, qlen, sid, cull_key(linked), *linked);
        }
        Ok(())
    }

    /// `s_BlastHSPCullingFinal`: the hit list of each query with a tree (see
    /// `hit_list_of_culled`); the HSP lists of a query are sorted by `ScoreCompareHSPs`
    /// after each tree.
    pub(in crate::algorithm::tblastn) fn finish(self) -> Result<Vec<Option<KappaResultHitList>>> {
        let mut hit_lists: Vec<Option<KappaResultHitList>> = Vec::new();
        hit_lists.resize_with(self.query_lengths.len(), || None);
        let culled = self
            .forest
            .finalize(|query| query, |a, b| score_compare(&a.hsp, &b.hsp));
        for (query, lists) in culled {
            let payloads = lists
                .iter()
                .map(|(_, hsps)| vec![None; hsps.len()])
                .collect();
            hit_lists[query] = Some(hit_list_of_culled(lists, payloads, self.prelim_size));
        }
        Ok(hit_lists)
    }
}

impl KappaResultHitList {
    /// NCBI's culling pipe on the final hit list of one query, the first step after the
    /// traceback (`BlastHSPStreamTBackClose`, before the mask level, `-max_hsps`,
    /// `-qcov_hsp_perc` and `-subject_besthit`): the HSP lists sorted by e-value, then every
    /// HSP in that order culled into a new tree with the merit `culling_limit` (without the
    /// three of the preliminary writer), and the kept HSPs made into the query's hit list
    /// again. The lists leave in the rip order of the tree, so the caller sorts them by
    /// e-value afterwards (NCBI's `Blast_HSPResultsSortByEvalue` after the pipes does).
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
    /// The trees are per query, so the pipe of a batch is this function for each of its
    /// queries. `query_length` is the length of the query (`contexts[0].query_length`).
    pub(in crate::algorithm::tblastn) fn cull(
        &mut self,
        culling_limit: i32,
        query_length: usize,
    ) -> Result<()> {
        let qlen = i32::try_from(query_length)?;
        for list in &mut self.lists {
            sort_hsps_by_evalue(list);
            if let Some(first) = list.hsps.hsps.first() {
                list.hsps.best_evalue = first.evalue;
            }
        }
        self.lists.sort_by(compare_result_lists);

        let mut forest = CullingForest::new(i64::from(culling_limit));
        for list in std::mem::take(&mut self.lists) {
            let sid = u32::try_from(list.oid)?;
            let KappaResultList { hsps, payloads, .. } = list;
            ensure!(
                payloads.is_empty() || payloads.len() == hsps.hsps.len(),
                "Kappa result payload count must match HSP count"
            );
            let mut payloads = payloads.into_iter();
            for linked in hsps.hsps {
                let key = cull_key(&linked);
                forest.save(0, qlen, sid, key, (linked, payloads.next()));
            }
        }

        let culled = forest.finalize(
            |_| 0,
            |a: &(LinkedHsp, Option<KappaHspPayload>), b| score_compare(&a.0.hsp, &b.0.hsp),
        );
        let mut lists = Vec::new();
        let mut payloads = Vec::new();
        for (_, query_lists) in culled {
            for (oid, pairs) in query_lists {
                let (hsps, list_payloads): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
                lists.push((oid, hsps));
                payloads.push(list_payloads);
            }
        }
        *self = hit_list_of_culled(lists, payloads, self.max);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::search_gapped::GappedHsp;

    fn hsp(context: usize, frame: i8, qs: i32, qe: i32, ss: i32, score: i32) -> LinkedHsp {
        LinkedHsp {
            context,
            hsp: GappedHsp {
                frame,
                score,
                q_start: qs,
                q_end: qe,
                q_gapped_start: qs,
                s_start: ss,
                s_end: ss + (qe - qs),
                s_gapped_start: ss,
            },
            num: 1,
            evalue: 1e-5 / f64::from(score),
            source_index: 0,
        }
    }

    fn summary(list: &KappaResultHitList) -> Vec<(i32, Vec<i32>)> {
        list.lists()
            .iter()
            .map(|l| (l.oid, l.hsps.hsps.iter().map(|h| h.hsp.score).collect()))
            .collect()
    }

    #[test]
    fn the_key_is_the_query_range_and_the_frame_local_subject_offset() {
        let key = cull_key(&hsp(0, -2, 10, 50, 7, 33));
        assert_eq!(
            key,
            CullKey {
                begin: 10,
                end: 50,
                score: 33,
                subject_offset: 7
            }
        );
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:79-120, 612-633: the frame is no key,
    // so an HSP on another frame over the same query range is dominated by a better one,
    // and the tree is the query's (`cid` = context), with the length of that query.
    #[test]
    fn the_frames_of_a_subject_share_the_tree_of_the_query() {
        let mut culling = PrelimCulling::new(&[1000, 500], 10, 1);
        culling
            .add(
                0,
                &[
                    hsp(0, 1, 100, 200, 50, 90),
                    hsp(0, -1, 100, 200, 60, 40),
                    hsp(0, 2, 400, 500, 10, 30),
                    // Another query: its own tree.
                    hsp(1, -3, 100, 200, 50, 20),
                ],
            )
            .unwrap();
        let lists = culling.finish().unwrap();
        assert_eq!(summary(lists[0].as_ref().unwrap()), vec![(0, vec![90, 30])]);
        assert_eq!(summary(lists[1].as_ref().unwrap()), vec![(0, vec![20])]);
    }

    // Two frames tie on all five keys when their frame-local offsets are equal; the one
    // written later loses (`s_DominateTest`, the later insert is the dominated one).
    #[test]
    fn frames_that_tie_on_the_five_keys_keep_the_first_written() {
        let mut culling = PrelimCulling::new(&[1000], 10, 1);
        let mut first = hsp(0, 1, 100, 200, 50, 90);
        first.num = 11;
        let mut second = hsp(0, -1, 100, 200, 50, 90);
        second.num = 22;
        culling.add(0, &[first, second]).unwrap();
        let lists = culling.finish().unwrap();
        let kept = lists[0].as_ref().unwrap();
        assert_eq!(kept.lists().len(), 1);
        assert_eq!(kept.lists()[0].hsps.hsps.len(), 1);
        assert_eq!(kept.lists()[0].hsps.hsps[0].num, 11);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:516-581: the hit list is made with
    // `prelim_hitlist_size` and holds more lists than that, in the rip order, with the best
    // e-value, worst e-value and low score set, and is not heapified.
    #[test]
    fn the_hit_list_is_not_cut_and_keeps_the_rip_order() {
        let mut culling = PrelimCulling::new(&[1000], 2, 5);
        for oid in 0..4usize {
            culling
                .add(
                    oid,
                    &[
                        hsp(
                            0,
                            1,
                            100 + 50 * oid as i32,
                            120 + 50 * oid as i32,
                            1,
                            40 - oid as i32,
                        ),
                        hsp(0, 1, 100 + 50 * oid as i32, 110 + 50 * oid as i32, 2, 20),
                    ],
                )
                .unwrap();
        }
        let lists = culling.finish().unwrap();
        let q0 = lists[0].as_ref().unwrap();
        assert_eq!(q0.max, 2);
        assert_eq!(q0.lists().len(), 4);
        assert!(!q0.heapified());
        assert_eq!(q0.low_score(), 37);
        // The best e-value of subject 3's list: the score 37 HSP.
        assert!((q0.worst_evalue() - 1e-5 / 37.0).abs() < 1e-18);
        // Limit 5 (+3 = 8) keeps everything: the newest subject first, scores descending.
        assert_eq!(
            summary(q0),
            vec![
                (3, vec![37, 20]),
                (2, vec![38, 20]),
                (1, vec![39, 20]),
                (0, vec![40, 20])
            ]
        );
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:100-112: the lower OID wins a tie, and
    // the subjects must come in OID order for the later one to lose.
    #[test]
    fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
        let mut culling = PrelimCulling::new(&[1000], 10, 1);
        culling.add(0, &[hsp(0, 1, 100, 200, 5, 80)]).unwrap();
        culling.add(1, &[hsp(0, 1, 100, 200, 5, 80)]).unwrap();
        let lists = culling.finish().unwrap();
        assert_eq!(summary(lists[0].as_ref().unwrap()), vec![(0, vec![80])]);
    }

    fn payload(bit_score: f64) -> KappaHspPayload {
        use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
        KappaHspPayload {
            bit_score,
            num_ident: 0,
            num_positives: 0,
            report_num_ident: 0,
            report_num_positives: 0,
            report_mismatches: 0,
            align_length: 0,
            mismatches: 0,
            gap_opens: 0,
            gap_letters: 0,
            edit_script: Vec::new(),
            matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
        }
    }

    // The pipe: lists sorted by e-value (payloads stay with their HSPs), the culled HSPs
    // dropped with their payloads, the survivors score sorted in the new lists.
    #[test]
    fn the_pipe_culls_in_evalue_order_and_keeps_the_payloads_paired() {
        let mut hit_list = KappaResultHitList::new(500).unwrap();
        // Subject 0: a weak HSP (e-value order puts it after the strong one of subject 1).
        let weak = hsp(0, 1, 100, 200, 5, 40);
        let strong = hsp(0, 1, 100, 200, 5, 90);
        let other = hsp(0, 1, 400, 480, 5, 30);
        let list0 = KappaResultList {
            oid: 0,
            hsps: LinkedHspList {
                hsps: vec![weak, other],
                best_evalue: weak.evalue,
            },
            payloads: vec![payload(40.0), payload(30.0)],
        };
        let list1 = KappaResultList {
            oid: 1,
            hsps: LinkedHspList {
                hsps: vec![strong],
                best_evalue: strong.evalue,
            },
            payloads: vec![payload(90.0)],
        };
        hit_list.update(list0).unwrap();
        hit_list.update(list1).unwrap();
        hit_list.cull(1, 1000).unwrap();
        // The weak HSP is dominated by the strong one (overlap 100 %, better score). The
        // lists leave in the rip order of the tree: the node of the 400-480 HSP holds it
        // before the left subtree that has the strong HSP.
        assert_eq!(summary(&hit_list), vec![(0, vec![30]), (1, vec![90])]);
        let bits: Vec<Vec<f64>> = hit_list
            .lists()
            .iter()
            .map(|l| l.payloads.iter().map(|p| p.bit_score).collect())
            .collect();
        assert_eq!(bits, vec![vec![30.0], vec![90.0]]);
        assert_eq!(hit_list.max, 500);
        assert_eq!(hit_list.low_score(), 30);
        assert!(!hit_list.heapified());
    }
}
