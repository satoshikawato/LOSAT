//! TBLASTN `-culling_limit`: NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback (`core::hspfilter_culling` holds the interval tree both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c
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
/// of the subject. NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:612-633 (the
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-587
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:318-342
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:602-644
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
    /// again. The lists leave in the rip order of the tree; the caller sorts them by e-value
    /// afterwards where NCBI's `Blast_HSPResultsSortByEvalue` after the pipes does (a positive
    /// `BlastSeqSrcGetTotLen`, blast_traceback.c:1768-1772).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:699-733
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
#[path = "culling_tests.rs"]
mod tests;
