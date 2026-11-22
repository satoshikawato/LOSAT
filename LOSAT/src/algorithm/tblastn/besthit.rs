//! TBLASTN `-best_hit_overhang` and `-best_hit_score_edge`: NCBI's best hit writer of the
//! preliminary stage and its best hit pipe after the traceback (`core::hspfilter_besthit` holds
//! the lists both share).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c
//!
//! This module is a child of `stage_d_results` for the private fields of
//! `KappaResultHitList`. Nothing here runs without one of the two options: the choice is made
//! once per call of `run_preliminary_subjects`, where it picks the sink of the preliminary hit
//! lists, and once per query batch where the search finishes the final hit lists.
//!
//! TBLASTN has one context for each protein query, of frame 0: the key of an HSP is its own
//! query range and the query length is never read. The six frames of a subject are not
//! contexts, so the HSPs of all six frames of a subject go into the list of their query.

use anyhow::{ensure, Result};

use super::super::stage_d_linking::{score_compare, LinkedHsp, LinkedHspList};
use super::{
    compare_result_lists, refresh_best_evalue, KappaHspPayload, KappaResultHitList, KappaResultList,
};
use crate::core::hspfilter_besthit::{BestHitKey, BestHitLists, BestHitParams};

/// The key of an HSP for the lists.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266
/// The context's frame is 0, so the start is the HSP's own `query.offset` (`q_start`, 0-based,
/// within the searched range of the query or the query part of a chunk) and the extent is
/// `query.end - query.offset`.
fn best_hit_key(linked: &LinkedHsp) -> BestHitKey {
    BestHitKey {
        begin: linked.hsp.q_start,
        len: linked.hsp.q_end - linked.hsp.q_start,
        score: linked.hsp.score,
        evalue: linked.evalue,
    }
}

/// One kept list of a query: the subject's OID, its HSPs, and their payloads (empty in the
/// preliminary stage, which has none).
type KeptList = (u32, Vec<LinkedHsp>, Vec<KappaHspPayload>);

/// The hit list that `s_BlastHSPBestHitFinal` makes of the lists kept for one query (one per
/// subject in the order of its first HSP in the query's list, each sorted by score): their best
/// e-values set, the lists sorted by e-value and saved into a hit list of `prelim_size` lists,
/// which keeps the lists of the best e-values when there are more (`Blast_HitListUpdate`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231
/// `s_ExportToHitlist` saves the lists into a hit list as large as the HSPs (it is never full,
/// so `Blast_HitListUpdate` only sets each list's best e-value) and `Blast_HitListSortByEvalue`
/// sorts it (`compare_result_lists`, a total order: the OIDs differ).
fn hit_list_of_kept(kept: Vec<KeptList>, prelim_size: usize) -> Result<KappaResultHitList> {
    let mut exported = Vec::with_capacity(kept.len());
    for (oid, hsps, payloads) in kept {
        ensure!(
            payloads.is_empty() || payloads.len() == hsps.len(),
            "Kappa result payload count must match HSP count"
        );
        let mut list = KappaResultList {
            oid: i32::try_from(oid)?,
            hsps: LinkedHspList {
                hsps,
                best_evalue: 0.0,
            },
            payloads,
        };
        refresh_best_evalue(&mut list.hsps);
        exported.push(list);
    }
    exported.sort_by(compare_result_lists);
    let mut hit_list = KappaResultHitList::new(prelim_size)?;
    for list in exported {
        hit_list.update(list)?;
    }
    Ok(hit_list)
}

/// NCBI's best hit writer in place of the hit-list collector: each subject's preliminary HSPs
/// go into the list of their query, and `finish` makes the hit lists of the queries. The
/// subjects must come in OID order and the HSPs of a subject in the order of its HSP list (the
/// order in which the collector receives them), because the result depends on the order of the
/// writes.
pub(in crate::algorithm::tblastn) struct PrelimBestHit {
    lists: BestHitLists<LinkedHsp>,
    /// The queries of the call (the query parts of a chunk in the chunk stage of a split batch).
    queries: usize,
    prelim_size: usize,
}

impl PrelimBestHit {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:318-329
    pub(in crate::algorithm::tblastn) fn new(
        queries: usize,
        prelim_size: usize,
        params: BestHitParams,
    ) -> Self {
        Self {
            lists: BestHitLists::writer(params, prelim_size),
            queries,
            prelim_size,
        }
    }

    /// `s_BlastHSPBestHitRun` for the HSPs of subject `oid`: all queries and all frames, in the
    /// order of the subject's list; the query of an HSP is its context.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:237-352
    pub(in crate::algorithm::tblastn) fn add(
        &mut self,
        oid: usize,
        hsps: &[LinkedHsp],
    ) -> Result<()> {
        let sid = u32::try_from(oid)?;
        for linked in hsps {
            self.lists
                .save(linked.context, sid, best_hit_key(linked), *linked);
        }
        Ok(())
    }

    /// `s_BlastHSPBestHitFinal`: the hit list of each query with kept HSPs (see
    /// `hit_list_of_kept`), `None` for the others.
    pub(in crate::algorithm::tblastn) fn finish(self) -> Result<Vec<Option<KappaResultHitList>>> {
        let mut hit_lists: Vec<Option<KappaResultHitList>> = Vec::new();
        hit_lists.resize_with(self.queries, || None);
        let prelim_size = self.prelim_size;
        let kept = self
            .lists
            .finalize(|a: &LinkedHsp, b: &LinkedHsp| score_compare(&a.hsp, &b.hsp));
        for (query, lists) in kept {
            let lists = lists
                .into_iter()
                .map(|(oid, hsps)| (oid, hsps, Vec::new()))
                .collect();
            hit_lists[query] = Some(hit_list_of_kept(lists, prelim_size)?);
        }
        Ok(hit_lists)
    }
}

impl KappaResultHitList {
    /// NCBI's best hit pipe on the final hit list of one query, the first step after the
    /// traceback (`BlastHSPStreamTBackClose`, before the mask level, `-max_hsps`,
    /// `-qcov_hsp_perc` and `-subject_besthit`): the HSP lists sorted by their best e-values (the
    /// HSPs of each list keep their order), every HSP in that order put into a new list with the
    /// overhang and the score edge of the options, and the kept HSPs made into the query's hit
    /// list again (of `prelim_size` lists, the size of the preliminary stage's). The payloads
    /// follow their HSPs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539
    /// The lists are per query, so the pipe of a batch is this function for each of its
    /// queries. The HSPs of each list are in the order the traceback left them (the order of
    /// NCBI's lists at `BlastHSPStreamTBackClose`: by score, or by e-value where the insert into
    /// the hit list heapified it, `KappaResultHitList::update`), which is the order of the
    /// writes; unlike the culling pipe nothing sorts them here.
    pub(in crate::algorithm::tblastn) fn best_hit(
        &mut self,
        params: BestHitParams,
        prelim_size: usize,
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
        // `Blast_HSPResultsSortByEvalue`: the lists by their best e-values (the smallest e-value
        // of each list, as the traceback's `Blast_HitListUpdate` set it).
        for list in &mut self.lists {
            refresh_best_evalue(&mut list.hsps);
        }
        self.lists.sort_by(compare_result_lists);

        let mut lists = BestHitLists::pipe(params, prelim_size);
        for list in std::mem::take(&mut self.lists) {
            let sid = u32::try_from(list.oid)?;
            let KappaResultList { hsps, payloads, .. } = list;
            ensure!(
                payloads.is_empty() || payloads.len() == hsps.hsps.len(),
                "Kappa result payload count must match HSP count"
            );
            let mut payloads = payloads.into_iter();
            for linked in hsps.hsps {
                let key = best_hit_key(&linked);
                lists.save(0, sid, key, (linked, payloads.next()));
            }
        }

        let kept = lists.finalize(|a: &(LinkedHsp, Option<KappaHspPayload>), b| {
            score_compare(&a.0.hsp, &b.0.hsp)
        });
        let mut kept_lists = Vec::new();
        for (_, query_lists) in kept {
            for (oid, pairs) in query_lists {
                let (hsps, payloads): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
                kept_lists.push((oid, hsps, payloads.into_iter().flatten().collect()));
            }
        }
        *self = hit_list_of_kept(kept_lists, prelim_size)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "besthit_tests.rs"]
mod tests;
