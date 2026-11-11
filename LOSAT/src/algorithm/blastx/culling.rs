//! BLASTX culling (`-culling_limit`): NCBI's culling writer of the preliminary stage and its
//! pipe after the traceback, on the interval tree that every program shares
//! (`core::hspfilter_culling::CullingForest`).
//!
//! BLASTX is not BLASTN in `s_BlastHSPCullingRun`: the tree of an HSP is its context (one of
//! the six frames of a query), its range is the context's own protein range (no fold onto the
//! plus strand) and the tree's length is the context's length. A query's hit list is filled
//! from the trees of its six contexts (`Blast_GetQueryIndexFromContext`: context / 6).
use anyhow::Result;

use super::{
    preliminary::compare_score,
    results::{HitList, Hsp, HspList, ResultsObserver},
};
use crate::core::hspfilter_culling::{CullKey, CullingForest};

/// The keys of `LinkedHSP` and `s_DominateTest` for a BLASTX HSP.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:621-628
/// ```c++
///       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
///     	  A.begin = qlen - A.hsp->query.end;
///     	  A.end   = qlen -  A.hsp->query.offset;
///       }
///       else {
///     	  A.begin = A.hsp->query.offset;
///     	  A.end   = A.hsp->query.end;
///       }
/// ```
/// and `s_DominateTest` (hspfilter_culling.c:79-120) reads `hsp->score` and
/// `hsp->subject.offset`.
fn cull_key(hsp: &Hsp) -> CullKey {
    CullKey {
        begin: hsp.hsp.q_start,
        end: hsp.hsp.q_end,
        score: hsp.hsp.score,
        subject_offset: i64::from(hsp.hsp.s_start),
    }
}

/// NCBI's culling writer (or pipe) over the HSPs of one batch or chunk: `write` is
/// `s_BlastHSPCullingRun` for the HSP list of one subject, `finalize` is
/// `s_BlastHSPCullingFinal`. The result depends on the order of the writes, which must be
/// NCBI's (subjects in OID order, the HSPs of a subject in the order of its list).
pub(crate) struct Culling {
    forest: CullingForest<Hsp>,
    /// The length of each context, `query_info->contexts[context].query_length`.
    lengths: Vec<i32>,
}

impl Culling {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:487-493
    /// ```c++
    /// static int
    /// s_BlastHSPCullingInit(void* data, void* results)
    /// {
    ///     BlastHSPCullingData * cull_data = data;
    ///     cull_data->c_tree = calloc(cull_data->num_contexts, sizeof(CTreeNode *));
    ///     return 0;
    /// }
    /// ```
    /// `culling_max` is the limit of the pipe, or `CullingForest::prelim_max` of it for the
    /// preliminary writer.
    pub fn new(lengths: Vec<i32>, culling_max: i64) -> Self {
        Self {
            forest: CullingForest::new(culling_max),
            lengths,
        }
    }

    /// `s_BlastHSPCullingRun` for the HSP list of subject `oid`; an HSP that is not saved is
    /// freed.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:614-638
    /// ```c++
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
    pub fn write(&mut self, oid: usize, hsps: Vec<Hsp>) -> Result<()> {
        let sid = u32::try_from(oid)?;
        for hsp in hsps {
            let context = hsp.hsp.context;
            let key = cull_key(&hsp);
            self.forest
                .save(context, self.lengths[context], sid, key, hsp);
        }
        Ok(())
    }

    /// `s_BlastHSPCullingFinal`: a hit list of at most `max` HSP lists for each of the
    /// `num_queries` queries (empty for a query without a tree, as the collector leaves it);
    /// the HSP lists of a query in the order of their subject's first HSP in the rip, each
    /// sorted by `ScoreCompareHSPs` and with its best e-value.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:517-582
    /// ```c++
    ///    for (cid=0; cid < cull_data->num_contexts; ++cid) {
    ///       if (c_tree[cid]) {
    ///          qid = Blast_GetQueryIndexFromContext(cid, params->program);
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
    /// ```
    /// The best e-value is the minimum over the list, so taking it once after the last sort
    /// gives NCBI's value.
    pub fn finalize(self, num_queries: usize, max: usize) -> Vec<HitList> {
        let mut queries: Vec<HitList> = (0..num_queries).map(|_| HitList::new(max)).collect();
        let culled = self
            .forest
            .finalize(|context| context / 6, |a, b| compare_score(&a.hsp, &b.hsp));
        for (query_index, lists) in culled {
            queries[query_index].lists = lists
                .into_iter()
                .map(|(sid, hsps)| {
                    let mut list = HspList {
                        query_index,
                        oid: sid as usize,
                        hsps,
                        best_evalue: 0.0,
                    };
                    list.refresh();
                    list
                })
                .collect();
        }
        queries
    }
}

/// NCBI's culling pipe after the traceback: the HSP lists of every query, each sorted by
/// e-value and the lists of a query sorted as the hit list, go through one culling forest
/// with the option's limit (no `+ 3`), whose final replaces the hit lists.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:705-728
/// ```c++
///    for (qid = 0; qid < results->num_queries; ++qid) {
///       if (!(results->hitlist_array[qid])) continue;
///       num_list = results->hitlist_array[qid]->hsplist_count;
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
/// ```
pub(crate) fn traceback_pipe(
    queries: &mut Vec<HitList>,
    lengths: Vec<i32>,
    limit: i32,
    max: usize,
    trace: &mut dyn ResultsObserver,
) -> Result<()> {
    let mut culling = Culling::new(lengths, i64::from(limit));
    for query in queries.iter_mut() {
        for list in &mut query.lists {
            list.sort_evalue();
            if let Some(h) = list.hsps.first() {
                list.best_evalue = h.evalue;
            }
        }
        query.lists.sort_by(super::results::compare_list);
    }
    for query in queries.iter_mut() {
        for list in std::mem::take(&mut query.lists) {
            culling.write(list.oid, list.hsps)?;
        }
    }
    *queries = culling.finalize(queries.len(), max);
    for q in queries.iter() {
        for l in &q.lists {
            trace.hsps("CULL_FINAL", l.oid, &l.hsps);
        }
    }
    Ok(())
}
