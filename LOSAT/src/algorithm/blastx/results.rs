//! Internal BLASTX HSP ownership, collector, stream and result reduction.
use super::{
    args::ResolvedOptions,
    linking::{LinkedHsp, LinkedHspList},
    preliminary::{compare_score, PreliminaryHsp},
    split::{QueryChunk, OVERLAP_AA},
};
use crate::common::GapEditOp;
use anyhow::{ensure, Result};
use std::cmp::Ordering;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
// ```c++
// Int2 Blast_HitListUpdate(BlastHitList* hit_list,
//                          BlastHSPList* hsp_list)
// {
//    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
//
// #ifndef NDEBUG
//    ASSERT(s_BlastCheckBestEvalue(hsp_list) == TRUE); /* NCBI_FAKE_WARNING */
// #endif /* _DEBUG */
//
//    if (hit_list->hsplist_count < hit_list->hsplist_max) {
//       /* If the array of HSP lists for this query is not yet allocated,
//          do it here */
//       if (hit_list->hsplist_current == hit_list->hsplist_count)
//       {
//          Int2 status = s_Blast_HitListGrowHSPListArray(hit_list);
//          if (status)
//            return status;
//       }
//       /* Just add to the end; sort later */
//       hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
//       hit_list->worst_evalue =
//          MAX(hsp_list->best_evalue, hit_list->worst_evalue);
//       hit_list->low_score =
//          MIN(hsp_list->hsp_array[0]->score, hit_list->low_score);
//    } else {
//       int evalue_order = 0;
//       if (!hit_list->heapified) {
//     	  /* make sure all hsp_list is sorted */
//           int index;
//           for (index =0; index < hit_list->hsplist_count; index++) {
//               Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
//               hit_list->hsplist_array[index]->best_evalue = s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
//           }
//           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
//                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
//           hit_list->heapified = TRUE;
//       }
//
//       /* make sure the hsp_list is sorted.  We actually do not need to sort
//          the full list: all that we need is the best score.   However, the
//          following code assumes hsp_list->hsp_array[0] has the best score. */
//       Blast_HSPListSortByEvalue(hsp_list);
//       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
//       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
//       if (evalue_order < 0) {
//          /* This hit list is less significant than any of those already saved;
//             discard it. Note that newer hits with score and e-value both equal
//             to the current worst will be saved, at the expense of some older
//             hit.
//          */
//          Blast_HSPListFree(hsp_list);
//       } else {
//          s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
//       }
//    }
//    return 0;
// ```
pub trait ResultsObserver {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:878-894
    // ```c++
    // s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // {
    //    Int4 query_index, subject_index;
    //    BlastHitList* hit_list;
    //
    //    for (query_index = 0; query_index < results->num_queries; ++query_index) {
    //       if (!(hit_list = results->hitlist_array[query_index]))
    //          continue;
    //       for (subject_index = hitlist_size;
    //            subject_index < hit_list->hsplist_count; ++subject_index) {
    //          hit_list->hsplist_array[subject_index] =
    //          Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
    //       }
    //       hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
    //    }
    // }
    //
    // ```
    fn targets(&mut self, _stage: &str, _q: usize, _lists: &[HspList]) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676
    // ```c++
    // s_Blast_HSPListReapByPrelimEvalue(BlastHSPList* hsp_list, const BlastHitSavingParameters* hit_params)
    // {
    //    BlastHSP* hsp;
    //    BlastHSP** hsp_array;
    //    Int4 hsp_cnt = 0;
    //    Int4 index;
    //    double cutoff;
    //
    //    if (hsp_list == NULL)
    //       return 0;
    //
    //    cutoff = hit_params->prelim_evalue;
    //
    //    hsp_array = hsp_list->hsp_array;
    //    for (index = 0; index < hsp_list->hspcnt; index++) {
    //       hsp = hsp_array[index];
    //
    //       ASSERT(hsp != NULL);
    //
    //       if (hsp->evalue > cutoff) {
    //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    //       } else {
    //          if (index > hsp_cnt)
    //             hsp_array[hsp_cnt] = hsp_array[index];
    //          hsp_cnt++;
    //       }
    //    }
    //
    //    hsp_list->hspcnt = hsp_cnt;
    //
    //    return 0;
    // }
    //
    // /** The core of the BLAST search: comparison between the (concatenated)
    // ```
    fn reap(&mut self, _stage: &str, _oid: usize, _list: &LinkedHspList, _cutoff: f64) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-433
    // ```c++
    // int BlastHSPStreamMerge(SSplitQueryBlk *squery_blk,
    //                              Uint4 chunk_num,
    //                              BlastHSPStream* stream1,
    //                              BlastHSPStream* stream2)
    // {
    //    Int4 i, j, k;
    //    BlastHSPResults *results1 = NULL;
    //    BlastHSPResults *results2 = NULL;
    //    Int4 contexts_per_query = 0;
    // #ifdef _DEBUG
    //    Int4 num_queries = 0, num_ctx = 0, num_ctx_offsets = 0;
    //    Int4 max_ctx;
    // #endif
    //
    //    Uint4 *query_list = NULL, *offset_list = NULL, num_contexts = 0;
    //    Int4 *context_list = NULL;
    //
    //
    //    if (!stream1 || !stream2)
    //        return kBlastHSPStream_Error;
    //
    //    s_FinalizeWriter(stream1);
    //    s_FinalizeWriter(stream2);
    //
    //    results1 = stream1->results;
    //    results2 = stream2->results;
    //
    //    contexts_per_query = BLAST_GetNumberOfContexts(stream2->program);
    //
    //    SplitQueryBlk_GetQueryIndicesForChunk(squery_blk, chunk_num, &query_list);
    //    SplitQueryBlk_GetQueryContextsForChunk(squery_blk, chunk_num,
    //                                           &context_list, &num_contexts);
    //    SplitQueryBlk_GetContextOffsetsForChunk(squery_blk, chunk_num, &offset_list);
    //
    // #if defined(_DEBUG_VERBOSE)
    // ```
    fn merge(&mut self, _chunk: &QueryChunk, _enter: bool) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:268-281
    // ```c++
    //  * @param hsp_list_out The read HSP list. [out]
    //  * @return Success, error, or end of reading, when nothing left to read.
    //  */
    // int BlastHSPStreamRead(BlastHSPStream* hsp_stream, BlastHSPList** hsp_list_out)
    // {
    //    *hsp_list_out = NULL;
    //
    //    if (!hsp_stream)
    //       return kBlastHSPStream_Error;
    //
    //    if (!hsp_stream->results)
    //       return kBlastHSPStream_Eof;
    //
    //    /* If this stream is not yet closed for writing, close it. In particular,
    // ```
    fn read(&mut self, _list: Option<&HspList>) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:591-615
    // ```c++
    //    /* return all the HSPlists with the same subject OID as the
    //       last HSPList in the collection stored. We assume there is
    //       at most one HSPList per query sequence */
    //
    //    num_hsplists = hsp_stream->num_hsplists;
    //    if (num_hsplists == 0)
    //       return kBlastHSPStream_Eof;
    //
    //    hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1];
    //    target_oid = hsplist->oid;
    //
    //    for (i = 0; i < num_hsplists; i++) {
    //        hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1 - i];
    //        if (hsplist->oid != target_oid)
    //            break;
    //
    //        batch->hsplist_array[i] = hsplist;
    //    }
    //
    //    hsp_stream->num_hsplists = num_hsplists - i;
    //    batch->num_hsplists = i;
    //
    //    return kBlastHSPStream_Success;
    // }
    //
    // ```
    fn batch_read(&mut self, _lists: &[HspList]) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-940
    // ```c++
    // BlastIntervalTreeContainsHSP(const BlastIntervalTree *tree,
    //                              const BlastHSP *hsp,
    //                              const BlastQueryInfo *query_info,
    //                              Int4 min_diag_separation)
    // {
    //     SIntervalNode *node = tree->nodes;
    //     Int4 query_start = s_GetQueryStrandOffset(query_info, hsp->context);
    //     Int4 region_start = query_start + hsp->query.offset;
    //     Int4 region_end = query_start + hsp->query.end;
    //     Int8 middle;
    // ```
    fn containment(&mut self, _h: &PreliminaryHsp, _contained: bool) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-568
    // ```c++
    //             Blast_HSPListAdjustOddBlastnScores(hsp_list, score_options->gapped_calculation, gap_align->sbp);
    // #endif
    //
    //         }
    //
    //         Blast_HSPListSortByScore(hsp_list);
    //
    //
    //         if (score_options->is_ooframe && kTranslatedSubject)
    //             subject->length = prot_length;
    //         } else {
    //             BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
    //                     hit_params->options, &hsp_list);
    //         }
    //
    //         if (hsp_list->hspcnt == 0) continue;
    //
    //         /* The subject ordinal id is not yet filled in this HSP list */
    //         hsp_list->oid = subject->oid;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-150
    // ```c++
    // static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
    //     Int8 b1 = p->begin;
    //     Int8 b2 = y->begin;
    //     Int8 e1 = p->end;
    //     Int8 e2 = y->end;
    //     Int8 s1 = p->hsp->score;
    //     Int8 s2 = y->hsp->score;
    //     Int8 l1 = e1 - b1;
    //     Int8 l2 = e2 - b2;
    //     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
    //     Int8 d = 0;
    //
    //     // If not overlap by more than 50%
    //     if(2 *overlap < l2) {
    //     	return FALSE;
    //     }
    //
    //     /* the main criterion:
    //        2 * (%diff in score) + 1 * (%diff in length) */
    //     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
    //     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
    //     // If identical, use oid as tie breaker
    //     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
    //     	if(s1 != s2) {
    //     		return (s1>s2);
    //     	}
    //     	if(p->sid != y->sid) {
    //     		return (p->sid < y->sid);
    //     	}
    //
    //     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
    //     		return FALSE;
    //     	}
    //     	return TRUE;
    //     }
    //
    //    	if (d < 0) {
    //    		return FALSE;
    //     }
    //
    //     return TRUE;
    // }
    //
    // /** check how many hsps in list dominates y, and update merit of y accordingly */
    // static Boolean s_FullPass(LinkedHSP *list, LinkedHSP *y) {
    //     LinkedHSP *p = list;
    //     while (p) {
    //        if (s_DominateTest(p, y)) {
    //           (y->merit)--;
    //           if (y->merit <= 0) return FALSE;
    //        }
    //        p = p->next;
    //     }
    //     return TRUE;
    // }
    //
    // /** update merit for hsps in list; also returns the number of hsps in list */
    // static Int4 s_ProcessHSPList(LinkedHSP **list, LinkedHSP *y) {
    //     Int4 num = 0;
    //     LinkedHSP *p = *list, *q, *r;
    //     q = p;
    //     while (p) {
    //        ++num;
    //        r = p;
    //        p = p->next;
    //        if (r != y && s_DominateTest(y, r)) {
    //           (r->merit)--;
    //           if (r->merit <= 0) {
    //              if (r == *list) {
    //                  *list = p;
    //                  q = p;
    //              } else {
    // ```
    fn cull_compare(&mut self, _ao: usize, _a: &Hsp, _bo: usize, _b: &Hsp, _result: bool) {}
    fn cull_delete(&mut self, _oid: usize, _hsp: &Hsp, _merit: i32) {}
    fn cull_save(&mut self, _oid: usize, _hsp: &Hsp) {}
    fn cull_saved(&mut self, _result: bool) {}
    fn cull_fork(&mut self, _begin: i32, _end: i32) {}
    fn cull_final(&mut self, _num_queries: usize) {}
    fn preliminary_call(&mut self, _oid: usize) {}
    fn numeric(&mut self, _stage: &str, _oid: usize, _hsps: &[PreliminaryHsp]) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2008
    // ```c++
    // Int2 Blast_HSPListReapByEvalue(BlastHSPList* hsp_list,
    //         const BlastHitSavingOptions* hit_options)
    // {
    //    BlastHSP* hsp;
    //    BlastHSP** hsp_array;
    //    Int4 hsp_cnt = 0;
    //    Int4 index;
    //    double cutoff;
    //
    //    if (hsp_list == NULL)
    //       return 0;
    //
    //    cutoff = hit_options->expect_value;
    //
    //    hsp_array = hsp_list->hsp_array;
    //    for (index = 0; index < hsp_list->hspcnt; index++) {
    //       hsp = hsp_array[index];
    //
    //       ASSERT(hsp != NULL);
    //
    //       if (hsp->evalue > cutoff) {
    //          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    //       } else {
    //          if (index > hsp_cnt)
    //             hsp_array[hsp_cnt] = hsp_array[index];
    //          hsp_cnt++;
    //       }
    //    }
    //
    //    hsp_list->hspcnt = hsp_cnt;
    //
    //    return 0;
    // }
    // ```
    fn linked(&mut self, _stage: &str, _oid: usize, _hsps: &LinkedHspList) {}
    fn hsps(&mut self, _stage: &str, _oid: usize, _hsps: &[Hsp]) {}
    fn hitlist(
        &mut self,
        _stage: &str,
        _oid: usize,
        _max: usize,
        _lists: &[HspList],
        _candidate: Option<&HspList>,
    ) {
    }
}
struct NoopResults;
impl ResultsObserver for NoopResults {}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1917-1934
// ```c++
//
//    kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
//
//    for (index=0; index<hsp_list->hspcnt; index++) {
//       hsp = hsp_list->hsp_array[index];
//       ASSERT(hsp != NULL);
// #if 0
//       ASSERT(sbp->round_down == FALSE || (hsp->score & 1) == 0);
// #endif
//       hsp->bit_score =
//          (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
//          NCBIMATH_LN2;
//    }
//
//    return 0;
// }
//
// void Blast_HSPListPHIGetBitScores(BlastHSPList* hsp_list, BlastScoreBlk* sbp)
// ```
#[derive(Clone, Debug)]
pub struct Hsp {
    pub hsp: PreliminaryHsp,
    pub num: i32,
    pub evalue: f64,
    pub bit_score: f64,
    pub identity: usize,
    pub positive: usize,
    pub edit_script: Vec<GapEditOp>,
    pub composition_method: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1390-1435
// ```c++
// s_EvalueComp(double evalue1, double evalue2)
// {
//     const double epsilon = 1.0e-180;
//     if (evalue1 < epsilon && evalue2 < epsilon) {
//         return 0;
//     }
//
//     if (evalue1 < evalue2) {
//         return -1;
//     } else if (evalue1 > evalue2) {
//         return 1;
//     } else {
//         return 0;
//     }
// }
//
// /** Comparison callback function for sorting HSPs by e-value and score, before
//  * saving BlastHSPList in a BlastHitList. E-value has priority over score,
//  * because lower scoring HSPs might have lower e-values, if they are linked
//  * with sum statistics.
//  * E-values are compared only up to a certain precision.
//  * @param v1 Pointer to first HSP [in]
//  * @param v2 Pointer to second HSP [in]
//  */
// static int
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    BlastHSP* h1,* h2;
//    int retval = 0;
//
//    h1 = *((BlastHSP**) v1);
//    h2 = *((BlastHSP**) v2);
//
//    /* Check if one or both of these are null. Those HSPs should go to the end */
//    if (!h1 && !h2)
//       return 0;
//    else if (!h1)
//       return 1;
//    else if (!h2)
//       return -1;
//
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
pub(crate) fn compare_evalue(a: f64, b: f64) -> Ordering {
    if a < 1.0e-180 && b < 1.0e-180 {
        Ordering::Equal
    } else {
        a.partial_cmp(&b).expect("finite BLAST E-value")
    }
}
pub(crate) fn compare_hsp(a: &Hsp, b: &Hsp) -> Ordering {
    compare_evalue(a.evalue, b.evalue).then_with(|| compare_score(&a.hsp, &b.hsp))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3078-3106
// ```c++
// s_EvalueCompareHSPLists(const void* v1, const void* v2)
// {
//    BlastHSPList* h1,* h2;
//    int retval = 0;
//
//    h1 = *(BlastHSPList**) v1;
//    h2 = *(BlastHSPList**) v2;
//
//    /* If any of the HSP lists is empty, it is considered "worse" than the
//       other, unless the other is also empty. */
//    if (h1->hspcnt == 0 && h2->hspcnt == 0)
//       return 0;
//    else if (h1->hspcnt == 0)
//       return 1;
//    else if (h2->hspcnt == 0)
//       return -1;
//
//    if ((retval = s_EvalueComp(h1->best_evalue,
//                                    h2->best_evalue)) != 0)
//       return retval;
//
//    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score)
//       return -1;
//    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score)
//       return 1;
//
//    /* In case of equal best E-values and scores, order will be determined
//       by ordinal ids of the subject sequences */
//    return BLAST_CMP(h2->oid, h1->oid);
// ```
#[derive(Clone, Debug)]
pub struct HspList {
    pub query_index: usize,
    pub oid: usize,
    pub hsps: Vec<Hsp>,
    pub best_evalue: f64,
}
impl HspList {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1737-1749
    // ```c++
    //  * @return TRUE if OK, FALSE otherwise.
    //  */
    // static double
    // s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
    // {
    //     int index = 0;
    //     double best_evalue = (double) INT4_MAX;
    //
    //     for (index=0; index<hsp_list->hspcnt; index++)
    //        best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
    //
    //     return best_evalue;
    // }
    // ```
    pub(crate) fn refresh(&mut self) {
        self.best_evalue = self
            .hsps
            .iter()
            .fold(i32::MAX as f64, |best, h| h.evalue.min(best));
    }
    pub(crate) fn sort_evalue(&mut self) {
        self.hsps.sort_by(compare_hsp);
    }
    pub(crate) fn linked(&self) -> LinkedHspList {
        LinkedHspList {
            hsps: self
                .hsps
                .iter()
                .enumerate()
                .map(|(source_index, h)| LinkedHsp {
                    hsp: h.hsp.clone(),
                    num: h.num,
                    evalue: h.evalue,
                    source_index,
                })
                .collect(),
            best_evalue: self.best_evalue,
        }
    }
}
pub(crate) fn compare_list(a: &HspList, b: &HspList) -> Ordering {
    match (a.hsps.first(), b.hsps.first()) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Greater,
        (Some(_), None) => Ordering::Less,
        (Some(ah), Some(bh)) => compare_evalue(a.best_evalue, b.best_evalue)
            .then_with(|| bh.hsp.score.cmp(&ah.hsp.score))
            .then_with(|| b.oid.cmp(&a.oid)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
// ```c++
// Int2 Blast_HitListUpdate(BlastHitList* hit_list,
//                          BlastHSPList* hsp_list)
// {
//    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
//
// #ifndef NDEBUG
//    ASSERT(s_BlastCheckBestEvalue(hsp_list) == TRUE); /* NCBI_FAKE_WARNING */
// #endif /* _DEBUG */
//
//    if (hit_list->hsplist_count < hit_list->hsplist_max) {
//       /* If the array of HSP lists for this query is not yet allocated,
//          do it here */
//       if (hit_list->hsplist_current == hit_list->hsplist_count)
//       {
//          Int2 status = s_Blast_HitListGrowHSPListArray(hit_list);
//          if (status)
//            return status;
//       }
//       /* Just add to the end; sort later */
//       hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
//       hit_list->worst_evalue =
//          MAX(hsp_list->best_evalue, hit_list->worst_evalue);
//       hit_list->low_score =
//          MIN(hsp_list->hsp_array[0]->score, hit_list->low_score);
//    } else {
//       int evalue_order = 0;
//       if (!hit_list->heapified) {
//     	  /* make sure all hsp_list is sorted */
//           int index;
//           for (index =0; index < hit_list->hsplist_count; index++) {
//               Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
//               hit_list->hsplist_array[index]->best_evalue = s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
//           }
//           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
//                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
//           hit_list->heapified = TRUE;
//       }
//
//       /* make sure the hsp_list is sorted.  We actually do not need to sort
//          the full list: all that we need is the best score.   However, the
//          following code assumes hsp_list->hsp_array[0] has the best score. */
//       Blast_HSPListSortByEvalue(hsp_list);
//       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
//       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
//       if (evalue_order < 0) {
//          /* This hit list is less significant than any of those already saved;
//             discard it. Note that newer hits with score and e-value both equal
//             to the current worst will be saved, at the expense of some older
//             hit.
//          */
//          Blast_HSPListFree(hsp_list);
//       } else {
//          s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
//       }
//    }
//    return 0;
// ```
pub(crate) struct HitList {
    pub lists: Vec<HspList>,
    max: usize,
    heapified: bool,
}
impl HitList {
    pub fn new(max: usize) -> Self {
        Self {
            lists: Vec::new(),
            max,
            heapified: false,
        }
    }
    fn down(&mut self, mut i: usize) {
        while i * 2 + 1 < self.lists.len() {
            let l = i * 2 + 1;
            let r = l + 1;
            let worse = if r == self.lists.len()
                || compare_list(&self.lists[l], &self.lists[r]) != Ordering::Less
            {
                l
            } else {
                r
            };
            if compare_list(&self.lists[i], &self.lists[worse]) != Ordering::Less {
                break;
            }
            self.lists.swap(i, worse);
            i = worse;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
    // ```c++
    // Int2 Blast_HitListUpdate(BlastHitList* hit_list,
    //                          BlastHSPList* hsp_list)
    // {
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //
    // #ifndef NDEBUG
    //    ASSERT(s_BlastCheckBestEvalue(hsp_list) == TRUE); /* NCBI_FAKE_WARNING */
    // #endif /* _DEBUG */
    //
    //    if (hit_list->hsplist_count < hit_list->hsplist_max) {
    //       /* If the array of HSP lists for this query is not yet allocated,
    //          do it here */
    //       if (hit_list->hsplist_current == hit_list->hsplist_count)
    //       {
    //          Int2 status = s_Blast_HitListGrowHSPListArray(hit_list);
    //          if (status)
    //            return status;
    //       }
    //       /* Just add to the end; sort later */
    //       hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
    //       hit_list->worst_evalue =
    //          MAX(hsp_list->best_evalue, hit_list->worst_evalue);
    //       hit_list->low_score =
    //          MIN(hsp_list->hsp_array[0]->score, hit_list->low_score);
    //    } else {
    //       int evalue_order = 0;
    //       if (!hit_list->heapified) {
    //     	  /* make sure all hsp_list is sorted */
    //           int index;
    //           for (index =0; index < hit_list->hsplist_count; index++) {
    //               Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
    //               hit_list->hsplist_array[index]->best_evalue = s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
    //           }
    //           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
    //                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //           hit_list->heapified = TRUE;
    //       }
    //
    //       /* make sure the hsp_list is sorted.  We actually do not need to sort
    //          the full list: all that we need is the best score.   However, the
    //          following code assumes hsp_list->hsp_array[0] has the best score. */
    //       Blast_HSPListSortByEvalue(hsp_list);
    //       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    //       if (evalue_order < 0) {
    //          /* This hit list is less significant than any of those already saved;
    //             discard it. Note that newer hits with score and e-value both equal
    //             to the current worst will be saved, at the expense of some older
    //             hit.
    //          */
    //          Blast_HSPListFree(hsp_list);
    //       } else {
    //          s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
    //       }
    //    }
    //    return 0;
    // ```
    pub fn update_observed(
        &mut self,
        candidate: HspList,
        trace: &mut dyn ResultsObserver,
    ) -> Option<HspList> {
        let oid = candidate.oid;
        trace.hitlist("IN", oid, self.max, &self.lists, Some(&candidate));
        let discarded = self.update(candidate);
        trace.hitlist("OUT", oid, self.max, &self.lists, None);
        discarded
    }
    pub fn update(&mut self, mut candidate: HspList) -> Option<HspList> {
        candidate.refresh();
        if self.lists.len() < self.max {
            self.lists.push(candidate);
            return None;
        }
        if !self.heapified {
            for l in &mut self.lists {
                l.sort_evalue();
            }
            for i in (0..self.lists.len() / 2).rev() {
                self.down(i);
            }
            self.heapified = true;
        }
        candidate.sort_evalue();
        if compare_list(&self.lists[0], &candidate) == Ordering::Less {
            return Some(candidate);
        }
        let discarded = std::mem::replace(&mut self.lists[0], candidate);
        self.down(0);
        Some(discarded)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:44-70
// ```c++
// GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats, Boolean gapped_calculation)
// {
//     Int4 prelim_hitlist_size = hitlist_size;
//     char * ADAPTIVE_CBS_ENV = getenv("ADAPTIVE_CBS");
//     if (compositionBasedStats) {
//     	if(ADAPTIVE_CBS_ENV != NULL) {
//     		if(hitlist_size < 1000) {
//     			prelim_hitlist_size = MAX(prelim_hitlist_size + 1000, 1500);
//     		}
//     		else {
//     			prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//     		}
//     	}
//     	else {
//     		if(hitlist_size <= 500) {
//     			prelim_hitlist_size = 1050;
//     		}
//     		else {
//     			prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//     		}
//
//     	}
//     }
//     else if (gapped_calculation) {
//          prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10), prelim_hitlist_size + 50);
//     }
//     return prelim_hitlist_size;
// ```
pub(crate) fn prelim_hitlist_size(o: &ResolvedOptions) -> usize {
    let n = o.hitlist_size as usize;
    if o.composition > 0 {
        if n <= 500 {
            1050
        } else {
            n * 2 + 50
        }
    } else if o.gapped {
        (2 * n).max(10).min(n + 50)
    } else {
        n
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2538-2588
// ```c++
// Blast_HSPListSubjectBestHit(EBlastProgramType program,
// 		                    const BlastHSPSubjectBestHitOptions* subject_besthit_opts,
// 		                    const BlastQueryInfo *query_info,
//                             BlastHSPList* hsp_list)
// {
//    BlastHSP** hsp_array;  /* hsp_array to purge. */
//    const int range_diff = subject_besthit_opts->max_range_diff;
//    Boolean isBlastn = (program == eBlastTypeBlastn);
//    unsigned int i, j;
//    int o, e;
//    int curr_context, target_context;
//    int qlen;
//
//    /* If HSP list is empty, return immediately. */
//    if (hsp_list == NULL || hsp_list->hspcnt == 0)
//        return 0;
//
//    if (Blast_ProgramIsPhiBlast(program))
//        return hsp_list->hspcnt;
//
//    hsp_array = hsp_list->hsp_array;
//
//    // The hsp list is sorted by score
//    for(i=0; i < hsp_list->hspcnt -1; i++) {
// 	  if(hsp_array[i] == NULL){
// 		  continue;
// 	  }
//       j = 1;
//       o = hsp_array[i]->query.offset - range_diff;
//       e = hsp_array[i]->query.end + range_diff;
//       if (o < 0) o = 0;
//       if (e < 0) e = hsp_array[i]->query.end;
//       while (i+j < hsp_list->hspcnt) {
//           if (hsp_array[i+j] && hsp_array[i]->context == hsp_array[i+j]->context &&
//               ((hsp_array[i+j]->query.offset >= o) &&
//                (hsp_array[i+j]->query.end <= e))){
//        	      hsp_array[i+j] = Blast_HSPFree(hsp_array[i+j]);
//           }
//           j++;
//       }
//    }
//
//    Blast_HSPListPurgeNullHSPs(hsp_list);
//
//    if(isBlastn) {
// 	   for(i=0; i < hsp_list->hspcnt -1; i++) {
// 	   	  if(hsp_array[i] == NULL){
// 	   		  continue;
// 	   	  }
// 	   	  // Flip query offsets of current hsp to target context frame
// 	   	  j = 1;
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:1211-1211
// ```c++
// #define DEFAULT_SUBJECT_BESTHIT_PROT_MAX_RANGE_DIFF 3
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1960-1973
// ```c++
// }
//
// BlastHSPSubjectBestHitOptions*
// BlastHSPSubjectBestHitOptionsNew(Boolean isProtein)
// {
//     BlastHSPSubjectBestHitOptions* retval =
//         (BlastHSPSubjectBestHitOptions*) calloc(1, sizeof(BlastHSPSubjectBestHitOptions));
//     if(isProtein){
//         retval->max_range_diff = DEFAULT_SUBJECT_BESTHIT_PROT_MAX_RANGE_DIFF;
//     }
//     else {
//         retval->max_range_diff = DEFAULT_SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF;
//     }
//     return retval;
// ```
pub(crate) fn subject_besthit(hsps: &mut Vec<Hsp>) -> Vec<Hsp> {
    let mut removed = vec![false; hsps.len()];
    for i in 0..hsps.len().saturating_sub(1) {
        if removed[i] {
            continue;
        }
        let a = &hsps[i].hsp;
        let begin = (a.q_start - 3).max(0);
        let end = a.q_end + 3;
        for j in i + 1..hsps.len() {
            let b = &hsps[j].hsp;
            if !removed[j] && a.context == b.context && b.q_start >= begin && b.q_end <= end {
                removed[j] = true;
            }
        }
    }
    let mut discarded = Vec::new();
    let mut i = 0;
    hsps.retain(|h| {
        let keep = !removed[i];
        i += 1;
        if !keep {
            discarded.push(h.clone());
        }
        keep
    });
    discarded
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-421
// ```c++
// int BlastHSPStreamMerge(SSplitQueryBlk *squery_blk,
//                              Uint4 chunk_num,
//                              BlastHSPStream* stream1,
//                              BlastHSPStream* stream2)
// {
//    Int4 i, j, k;
//    BlastHSPResults *results1 = NULL;
//    BlastHSPResults *results2 = NULL;
//    Int4 contexts_per_query = 0;
// #ifdef _DEBUG
//    Int4 num_queries = 0, num_ctx = 0, num_ctx_offsets = 0;
//    Int4 max_ctx;
// #endif
//
//    Uint4 *query_list = NULL, *offset_list = NULL, num_contexts = 0;
//    Int4 *context_list = NULL;
//
//
//    if (!stream1 || !stream2)
//        return kBlastHSPStream_Error;
//
//    s_FinalizeWriter(stream1);
//    s_FinalizeWriter(stream2);
// ```
pub(crate) struct HspStream {
    pub queries: Vec<HitList>,
    sorted: Option<Vec<HspList>>,
    finalized: bool,
    culling: Option<super::culling::Culling>,
    max: usize,
}
impl HspStream {
    pub fn new(num_queries: usize, max: usize) -> Self {
        Self {
            queries: (0..num_queries).map(|_| HitList::new(max)).collect(),
            sorted: None,
            finalized: false,
            culling: None,
            max,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-161
    // ```c++
    //    /* Rearrange HSPs into multiple hit lists if more than one query */
    //    if (results->num_queries > 1) {
    //       BlastHSP* hsp;
    //       BlastHSPList** hsp_list_array;
    //       BlastHSPList* tmp_hsp_list;
    //       Int4 index;
    //
    //       hsp_list_array = calloc(results->num_queries, sizeof(BlastHSPList*));
    //       if (hsp_list_array == NULL)
    //          return -1;
    //
    //       for (index = 0; index < hsp_list->hspcnt; index++) {
    //          Int4 query_index;
    //          hsp = hsp_list->hsp_array[index];
    //          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
    //
    //          if (!(tmp_hsp_list = hsp_list_array[query_index])) {
    //             hsp_list_array[query_index] = tmp_hsp_list =
    //                Blast_HSPListNew(params->hsp_num_max);
    //             if (tmp_hsp_list == NULL)
    //             {
    //                  sfree(hsp_list_array);
    //                  return -1;
    //             }
    //             tmp_hsp_list->oid = hsp_list->oid;
    //          }
    //
    //          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
    //          hsp_list->hsp_array[index] = NULL;
    //       }
    //
    //       /* All HSPs from the hsp_list structure are now moved to the results
    //          structure, so set the HSP count back to 0 */
    //       hsp_list->hspcnt = 0;
    //       Blast_HSPListFree(hsp_list);
    //
    //       /* Insert the hit list(s) into the appropriate places in the results
    //          structure */
    //       for (index = 0; index < results->num_queries; index++) {
    //          if (hsp_list_array[index]) {
    //             if (!results->hitlist_array[index]) {
    //                results->hitlist_array[index] =
    //                   Blast_HitListNew(params->prelim_hitlist_size);
    //             }
    //             Blast_HitListUpdate(results->hitlist_array[index],
    //                                 hsp_list_array[index]);
    //          }
    //       }
    //       sfree(hsp_list_array);
    //    } else if (hsp_list->hspcnt > 0) {
    //       /* Single query; save the HSP list directly into the results
    //          structure */
    //       if (!results->hitlist_array[0]) {
    //          results->hitlist_array[0] =
    //             Blast_HitListNew(params->prelim_hitlist_size);
    //       }
    //       Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:330-342
    // ```c++
    //         else if (filt_opts->culling_opts &&
    //                  (filt_opts->culling_stage & ePrelimSearch))
    //         {
    //             BlastHSPCullingParams* params =
    //                 BlastHSPCullingParamsNew(opts_memento->m_HitSaveOpts,
    //                      filt_opts->culling_opts,
    //                      opts_memento->m_ExtnOpts->compositionBasedStats,
    //                      opts_memento->m_ScoringOpts->gapped_calculation);
    //             if(params->culling_max > 1){
    //             	params->culling_max += 3;
    //             }
    //             writer_info = BlastHSPCullingInfoNew(params);
    //         }
    // ```
    pub fn with_culling(lengths: Vec<i32>, limit: i32, max: usize) -> Self {
        let mut stream = Self::new(lengths.len() / 6, max);
        stream.culling = Some(super::culling::Culling::new(
            lengths,
            if limit > 1 { limit + 3 } else { limit },
        ));
        stream
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:604-647
    // ```c++
    // {
    //    Int4 i, qlen;
    //    LinkedHSP A;
    //
    //    BlastHSPCullingData * cull_data = data;
    //    BlastHSPCullingParams* params = cull_data->params;
    //    CTreeNode **c_tree = cull_data->c_tree;
    //    Boolean isBlastn = (params->program == eBlastTypeBlastn);
    //    if (!hsp_list) return 0;
    //
    //    for (i=0; i<hsp_list->hspcnt; ++i) {
    //       /* wrap the hsp with a LinkedHSP structure */
    //       A.hsp   = hsp_list->hsp_array[i];
    //       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
    //       A.sid   = hsp_list->oid;
    //       A.merit = params->culling_max;
    //       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
    //       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
    //     	  A.begin = qlen - A.hsp->query.end;
    //     	  A.end   = qlen -  A.hsp->query.offset;
    //       }
    //       else {
    //     	  A.begin = A.hsp->query.offset;
    //     	  A.end   = A.hsp->query.end;
    //       }
    //       A.next  = NULL;
    //
    //       if (! c_tree[A.cid]) {
    //          c_tree[A.cid] = s_CTreeNew(qlen);
    //       }
    //
    //       if(s_SaveHSP(c_tree[A.cid], &A)){
    //     	 hsp_list->hsp_array[i] = NULL;
    //       }
    //    }
    //
    //    /* now all good hits have moved to tree, we can remove hsp_list */
    //    Blast_HSPListFree(hsp_list);
    //
    //    return 0;
    // }
    //
    // /** Free the writer
    //  * @param writer The writer to free [in]
    // ```
    pub fn write(&mut self, oid: usize, hsps: Vec<Hsp>) -> Result<Vec<HspList>> {
        self.write_observed(oid, hsps, &mut NoopResults)
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:105-161
    // ```c++
    //    /* Rearrange HSPs into multiple hit lists if more than one query */
    //    if (results->num_queries > 1) {
    //       BlastHSP* hsp;
    //       BlastHSPList** hsp_list_array;
    //       BlastHSPList* tmp_hsp_list;
    //       Int4 index;
    //
    //       hsp_list_array = calloc(results->num_queries, sizeof(BlastHSPList*));
    //       if (hsp_list_array == NULL)
    //          return -1;
    //
    //       for (index = 0; index < hsp_list->hspcnt; index++) {
    //          Int4 query_index;
    //          hsp = hsp_list->hsp_array[index];
    //          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
    //
    //          if (!(tmp_hsp_list = hsp_list_array[query_index])) {
    //             hsp_list_array[query_index] = tmp_hsp_list =
    //                Blast_HSPListNew(params->hsp_num_max);
    //             if (tmp_hsp_list == NULL)
    //             {
    //                  sfree(hsp_list_array);
    //                  return -1;
    //             }
    //             tmp_hsp_list->oid = hsp_list->oid;
    //          }
    //
    //          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
    //          hsp_list->hsp_array[index] = NULL;
    //       }
    //
    //       /* All HSPs from the hsp_list structure are now moved to the results
    //          structure, so set the HSP count back to 0 */
    //       hsp_list->hspcnt = 0;
    //       Blast_HSPListFree(hsp_list);
    //
    //       /* Insert the hit list(s) into the appropriate places in the results
    //          structure */
    //       for (index = 0; index < results->num_queries; index++) {
    //          if (hsp_list_array[index]) {
    //             if (!results->hitlist_array[index]) {
    //                results->hitlist_array[index] =
    //                   Blast_HitListNew(params->prelim_hitlist_size);
    //             }
    //             Blast_HitListUpdate(results->hitlist_array[index],
    //                                 hsp_list_array[index]);
    //          }
    //       }
    //       sfree(hsp_list_array);
    //    } else if (hsp_list->hspcnt > 0) {
    //       /* Single query; save the HSP list directly into the results
    //          structure */
    //       if (!results->hitlist_array[0]) {
    //          results->hitlist_array[0] =
    //             Blast_HitListNew(params->prelim_hitlist_size);
    //       }
    //       Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
    // ```
    pub fn write_observed(
        &mut self,
        oid: usize,
        hsps: Vec<Hsp>,
        trace: &mut dyn ResultsObserver,
    ) -> Result<Vec<HspList>> {
        ensure!(
            !self.finalized && self.sorted.is_none(),
            "BLASTX HSP stream is finalized"
        );
        if let Some(culling) = &mut self.culling {
            culling.write(oid, hsps, trace);
            return Ok(Vec::new());
        }
        let mut lists: Vec<Vec<Hsp>> = (0..self.queries.len()).map(|_| Vec::new()).collect();
        for h in hsps {
            lists[h.hsp.context / 6].push(h);
        }
        let mut discarded = Vec::new();
        for (query_index, hsps) in lists.into_iter().enumerate() {
            if hsps.is_empty() {
                continue;
            }
            let list = HspList {
                query_index,
                oid,
                hsps,
                best_evalue: 0.0,
            };
            if let Some(l) = self.queries[query_index].update_observed(list, trace) {
                discarded.push(l);
            }
        }
        Ok(discarded)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:99-130
    // ```c++
    //  * and therefore must be finalized before reading/merging
    //  */
    // static void s_FinalizeWriter(BlastHSPStream* hsp_stream)
    // {
    //    BlastHSPPipe *pipe;
    //    if (!hsp_stream || !hsp_stream->results || hsp_stream->writer_finalized)
    //       return;
    //
    //    /* perform post-writer clean ups */
    //    if (hsp_stream->writer) {
    //        if (!hsp_stream->writer_initialized) {
    //            /* some filter (e.g. hsp_queue) always needs finalization */
    //            (hsp_stream->writer->InitFnPtr)
    //                 (hsp_stream->writer->data, hsp_stream->results);
    //        }
    //        (hsp_stream->writer->FinalFnPtr)
    //             (hsp_stream->writer->data, hsp_stream->results);
    //    }
    //
    //    /* apply preliminary stage pipes */
    //    while (hsp_stream->pre_pipe) {
    //        pipe = hsp_stream->pre_pipe;
    //        hsp_stream->pre_pipe = pipe->next;
    //        (pipe->RunFnPtr) (pipe->data, hsp_stream->results);
    //        (pipe->FreeFnPtr) (pipe);
    //    }
    //
    //    hsp_stream->writer_finalized = TRUE;
    // }
    //
    // /** Prohibit any future writing to the HSP stream when all results are written.
    //  * Also perform sorting of results here to prepare them for reading.
    // ```
    pub fn finalize(&mut self) {
        self.finalize_observed(&mut NoopResults);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:105-130
    // ```c++
    //       return;
    //
    //    /* perform post-writer clean ups */
    //    if (hsp_stream->writer) {
    //        if (!hsp_stream->writer_initialized) {
    //            /* some filter (e.g. hsp_queue) always needs finalization */
    //            (hsp_stream->writer->InitFnPtr)
    //                 (hsp_stream->writer->data, hsp_stream->results);
    //        }
    //        (hsp_stream->writer->FinalFnPtr)
    //             (hsp_stream->writer->data, hsp_stream->results);
    //    }
    //
    //    /* apply preliminary stage pipes */
    //    while (hsp_stream->pre_pipe) {
    //        pipe = hsp_stream->pre_pipe;
    //        hsp_stream->pre_pipe = pipe->next;
    //        (pipe->RunFnPtr) (pipe->data, hsp_stream->results);
    //        (pipe->FreeFnPtr) (pipe);
    //    }
    //
    //    hsp_stream->writer_finalized = TRUE;
    // }
    //
    // /** Prohibit any future writing to the HSP stream when all results are written.
    //  * Also perform sorting of results here to prepare them for reading.
    // ```
    pub fn finalize_observed(&mut self, trace: &mut dyn ResultsObserver) {
        if self.finalized {
            return;
        }
        if let Some(culling) = &mut self.culling {
            self.queries = culling.finalize(self.queries.len(), self.max);
            trace.cull_final(self.queries.len());
            for q in &self.queries {
                for l in &q.lists {
                    trace.hsps("CULL_FINAL", l.oid, &l.hsps);
                }
            }
        }
        self.finalized = true;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:456-533
    // ```c++
    //    for (i = 0; i < results1->num_queries; i++) {
    //        BlastHitList *hitlist = results1->hitlist_array[i];
    //        Int4 global_query = query_list[i];
    //        Int4 split_points[NUM_FRAMES];
    // #ifdef _DEBUG
    //        ASSERT(i < num_queries);
    // #endif
    //
    //        if (hitlist == NULL) {
    // #if defined(_DEBUG_VERBOSE)
    // fprintf(stderr, "No hits to query %d\n", global_query);
    // #endif
    //            continue;
    //        }
    //
    //        /* we will be mapping HSPs from the local context to
    //           their place on the unsplit concatenated query. Once
    //           that's done, overlapping HSPs need to get merged, and
    //           to do that we must know the offset within each context
    //           where the last chunk ended and the current chunk begins */
    //        for (j = 0; j < contexts_per_query; j++) {
    //            split_points[j] = -1;
    //        }
    //
    //        for (j = 0; j < contexts_per_query; j++) {
    //            Int4 local_context = i * contexts_per_query + j;
    //            if (context_list[local_context] >= 0) {
    //                split_points[context_list[local_context] % contexts_per_query] =
    //                                 offset_list[local_context];
    //            }
    //        }
    //
    // #if defined(_DEBUG_VERBOSE)
    //        fprintf(stderr, "query %d split points: ", i);
    //        for (j = 0; j < contexts_per_query; j++) {
    //            fprintf(stderr, "%d ", split_points[j]);
    //        }
    //        fprintf(stderr, "\n");
    // #endif
    //
    //        for (j = 0; j < hitlist->hsplist_count; j++) {
    //            BlastHSPList *hsplist = hitlist->hsplist_array[j];
    //
    //            for (k = 0; k < hsplist->hspcnt; k++) {
    //                BlastHSP *hsp = hsplist->hsp_array[k];
    //                Int4 local_context = hsp->context;
    // #ifdef _DEBUG
    //                ASSERT(local_context <= max_ctx);
    //                ASSERT(local_context < num_ctx);
    //                ASSERT(local_context < num_ctx_offsets);
    // #endif
    //
    //                hsp->context = context_list[local_context];
    //                hsp->query.offset += offset_list[local_context];
    //                hsp->query.end += offset_list[local_context];
    //                hsp->query.gapped_start += offset_list[local_context];
    //                hsp->query.frame = BLAST_ContextToFrame(stream2->program,
    //                                                        hsp->context);
    //            }
    //
    //            hsplist->query_index = global_query;
    //        }
    //
    //        Blast_HitListMerge(results1->hitlist_array + i,
    //                           results2->hitlist_array + global_query,
    //                           contexts_per_query, split_points,
    //                           (Int4)SplitQueryBlk_GetChunkOverlapSize(squery_blk),
    //                           SplitQueryBlk_AllowGap(squery_blk));
    //    }
    //
    //    /* Sort to the canonical order, which the merge may not have done. */
    //    for (i = 0; i < results2->num_queries; i++) {
    //        BlastHitList *hitlist = results2->hitlist_array[i];
    //        if (hitlist == NULL)
    //            continue;
    //
    //        for (j = 0; j < hitlist->hsplist_count; j++)
    //            Blast_HSPListSortByScore(hitlist->hsplist_array[j]);
    // ```
    pub fn merge(&mut self, chunk_stream: Self, chunk: &QueryChunk) -> Result<Vec<HspList>> {
        self.merge_observed(chunk_stream, chunk, &mut NoopResults)
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:419-421
    // ```c++
    //
    //    s_FinalizeWriter(stream1);
    //    s_FinalizeWriter(stream2);
    // ```
    pub fn merge_observed(
        &mut self,
        mut chunk_stream: Self,
        chunk: &QueryChunk,
        trace: &mut dyn ResultsObserver,
    ) -> Result<Vec<HspList>> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:399-421
        // ```c++
        // int BlastHSPStreamMerge(SSplitQueryBlk *squery_blk,
        //                              Uint4 chunk_num,
        //                              BlastHSPStream* stream1,
        //                              BlastHSPStream* stream2)
        // {
        //    Int4 i, j, k;
        //    BlastHSPResults *results1 = NULL;
        //    BlastHSPResults *results2 = NULL;
        //    Int4 contexts_per_query = 0;
        // #ifdef _DEBUG
        //    Int4 num_queries = 0, num_ctx = 0, num_ctx_offsets = 0;
        //    Int4 max_ctx;
        // #endif
        //
        //    Uint4 *query_list = NULL, *offset_list = NULL, num_contexts = 0;
        //    Int4 *context_list = NULL;
        //
        //
        //    if (!stream1 || !stream2)
        //        return kBlastHSPStream_Error;
        //
        //    s_FinalizeWriter(stream1);
        //    s_FinalizeWriter(stream2);
        // ```
        trace.merge(chunk, true);
        ensure!(self.sorted.is_none(), "cannot merge a read BLASTX stream");
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:419-420
        // ```c++
        //
        //    s_FinalizeWriter(stream1);
        // ```
        chunk_stream.finalize_observed(trace);
        self.finalize_observed(trace);
        let mut discarded = Vec::new();
        for (local_query, mut incoming_hitlist) in chunk_stream.queries.into_iter().enumerate() {
            if incoming_hitlist.lists.is_empty() {
                continue;
            }
            let global_query = chunk.absolute_contexts[local_query * 6..local_query * 6 + 6]
                .iter()
                .copied()
                .find(|c| *c >= 0)
                .expect("query with HSP has context") as usize
                / 6;
            let mut points = [-1; 6];
            for local in local_query * 6..local_query * 6 + 6 {
                if chunk.absolute_contexts[local] >= 0 {
                    points[chunk.absolute_contexts[local] as usize % 6] =
                        chunk.corrections[local] as i32;
                }
            }
            for list in &mut incoming_hitlist.lists {
                list.query_index = global_query;
                for h in &mut list.hsps {
                    let local = h.hsp.context;
                    ensure!(
                        chunk.absolute_contexts[local] >= 0,
                        "HSP in invalid split context"
                    );
                    h.hsp.context = chunk.absolute_contexts[local] as usize;
                    let offset = chunk.corrections[local] as i32;
                    h.hsp.q_start += offset;
                    h.hsp.q_end += offset;
                    h.hsp.q_gapped_start += offset;
                    h.hsp.frame = [1, 2, 3, -1, -2, -3][h.hsp.context % 6];
                }
            }
            let target = &mut self.queries[global_query];
            if target.lists.is_empty() {
                *target = incoming_hitlist;
                continue;
            }
            incoming_hitlist.lists.sort_by_key(|l| l.oid);
            target.lists.sort_by_key(|l| l.oid);
            let mut combined = HitList::new(incoming_hitlist.max);
            let mut old = std::mem::take(&mut target.lists).into_iter().peekable();
            let mut new = incoming_hitlist.lists.into_iter().peekable();
            while old.peek().is_some() || new.peek().is_some() {
                let list = match (old.peek(), new.peek()) {
                    (Some(a), Some(b)) if a.oid == b.oid => {
                        let mut a = old.next().unwrap();
                        let b = new.next().unwrap();
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2884
                        // ```c++
                        // Int2 Blast_HSPListsMerge(BlastHSPList** hsp_list_ptr,
                        //                    BlastHSPList** combined_hsp_list_ptr,
                        //                    Int4 hsp_num_max, Int4 *split_offsets,
                        //                    Int4 contexts_per_query, Int4 chunk_overlap_size,
                        //                    Boolean allow_gap, Boolean short_reads)
                        // {
                        //    BlastHSPList* combined_hsp_list = *combined_hsp_list_ptr;
                        //    BlastHSPList* hsp_list = *hsp_list_ptr;
                        //    BlastHSP* hsp1, *hsp2, *hsp_var;
                        //    BlastHSP** hspp1,** hspp2;
                        //    Int4 index1, index2;
                        //    Int4 hspcnt1, hspcnt2, new_hspcnt = 0;
                        //    Int4 start_diag, end_diag;
                        //    Int4 offset_idx;
                        //    BlastHSP** new_hsp_array;
                        //
                        //    if (!hsp_list || hsp_list->hspcnt == 0)
                        //       return 0;
                        //
                        //    /* If no previous HSP list, just return a copy of the new one. */
                        //    if (!combined_hsp_list) {
                        //       *combined_hsp_list_ptr = hsp_list;
                        //       *hsp_list_ptr = NULL;
                        //       return 0;
                        //    }
                        //
                        //    /* Merge the two HSP lists for successive chunks of the subject sequence.
                        //       First put all HSPs that intersect the overlap region at the front of
                        // ```
                        if points.iter().any(|p| *p > 0) {
                            trace.hsps("OVERLAP_OLD", a.oid, &a.hsps);
                            trace.hsps("OVERLAP_NEW", b.oid, &b.hsps);
                            merge_hsps(&mut a.hsps, b.hsps, &points);
                            trace.hsps("OVERLAP_OUT", a.oid, &a.hsps);
                        } else {
                            a.hsps.extend(b.hsps);
                            a.hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
                        }
                        a
                    }
                    (Some(a), Some(b)) if a.oid > b.oid => new.next().unwrap(),
                    (Some(_), _) => old.next().unwrap(),
                    _ => new.next().unwrap(),
                };
                if let Some(l) = combined.update_observed(list, trace) {
                    discarded.push(l);
                }
            }
            *target = combined;
        }
        for q in &mut self.queries {
            for l in &mut q.lists {
                l.hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:521-537
        // ```c++
        //                           contexts_per_query, split_points,
        //                           (Int4)SplitQueryBlk_GetChunkOverlapSize(squery_blk),
        //                           SplitQueryBlk_AllowGap(squery_blk));
        //    }
        //
        //    /* Sort to the canonical order, which the merge may not have done. */
        //    for (i = 0; i < results2->num_queries; i++) {
        //        BlastHitList *hitlist = results2->hitlist_array[i];
        //        if (hitlist == NULL)
        //            continue;
        //
        //        for (j = 0; j < hitlist->hsplist_count; j++)
        //            Blast_HSPListSortByScore(hitlist->hsplist_array[j]);
        //    }
        //
        //    stream2->results_sorted = FALSE;
        //
        // ```
        trace.merge(chunk, false);
        Ok(discarded)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:173-202
    // ```c++
    //            Int4 alloc = MAX(num_hsplists + hitlist->hsplist_count + 100,
    //                             2 * hsp_stream->num_hsplists_alloc);
    //            hsp_stream->num_hsplists_alloc = alloc;
    //            hsp_stream->sorted_hsplists = (BlastHSPList **)realloc(
    //                                          hsp_stream->sorted_hsplists,
    //                                          alloc * sizeof(BlastHSPList *));
    //        }
    //
    //        for (j = k = 0; j < hitlist->hsplist_count; j++) {
    //
    //            BlastHSPList *hsplist = hitlist->hsplist_array[j];
    //            if (hsplist == NULL)
    //                continue;
    //
    //            hsplist->query_index = i;
    //            hsp_stream->sorted_hsplists[num_hsplists + k] = hsplist;
    //            k++;
    //        }
    //
    //        hitlist->hsplist_count = 0;
    //        num_hsplists += k;
    //    }
    //
    //    /* sort in order of decreasing subject OID. HSPLists will be
    //       read out from the end of hsplist_array later */
    //
    //    hsp_stream->num_hsplists = num_hsplists;
    //    if (num_hsplists > 1) {
    //       qsort(hsp_stream->sorted_hsplists, num_hsplists,
    //                     sizeof(BlastHSPList *), s_SortHSPListByOid);
    // ```
    pub fn close(&mut self) {
        self.close_observed(&mut NoopResults);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    // ```c++
    //        if (hsp_stream->sort_by_score->sort_on_read) {
    //            Blast_HSPResultsReverseSort(hsp_stream->results);
    //        } else {
    //            /* Reverse the order of HSP lists, because they will be returned
    //               starting from end, for the sake of convenience */
    //            Blast_HSPResultsReverseOrder(hsp_stream->results);
    //        }
    //        hsp_stream->results_sorted = TRUE;
    //        hsp_stream->x_lock = MT_LOCK_Delete(hsp_stream->x_lock);
    //        return;
    //    }
    //
    //    results = hsp_stream->results;
    //    num_hsplists = hsp_stream->num_hsplists;
    //
    //    /* concatenate all the HSPLists from 'results' */
    //
    //    for (i = 0; i < results->num_queries; i++) {
    //
    //        BlastHitList *hitlist = results->hitlist_array[i];
    // ```
    pub fn close_observed(&mut self, trace: &mut dyn ResultsObserver) {
        if self.sorted.is_some() {
            return;
        }
        self.finalize_observed(trace);
        let mut lists = Vec::new();
        for (query_index, q) in self.queries.iter_mut().enumerate() {
            for mut l in std::mem::take(&mut q.lists) {
                l.query_index = query_index;
                lists.push(l);
            }
        }
        lists.sort_by(|a, b| b.oid.cmp(&a.oid));
        self.sorted = Some(lists);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:320-329
    // ```c++
    //    } else {
    //        /* return the next HSPlist out of the collection stored */
    //
    //        if (!hsp_stream->num_hsplists)
    //           return kBlastHSPStream_Eof;
    //
    //        *hsp_list_out =
    //            hsp_stream->sorted_hsplists[--hsp_stream->num_hsplists];
    //
    //    }
    // ```
    pub fn read(&mut self) -> Option<HspList> {
        self.read_observed(&mut NoopResults)
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    // ```c++
    //        if (hsp_stream->sort_by_score->sort_on_read) {
    //            Blast_HSPResultsReverseSort(hsp_stream->results);
    //        } else {
    //            /* Reverse the order of HSP lists, because they will be returned
    //               starting from end, for the sake of convenience */
    //            Blast_HSPResultsReverseOrder(hsp_stream->results);
    //        }
    //        hsp_stream->results_sorted = TRUE;
    //        hsp_stream->x_lock = MT_LOCK_Delete(hsp_stream->x_lock);
    //        return;
    //    }
    //
    //    results = hsp_stream->results;
    //    num_hsplists = hsp_stream->num_hsplists;
    //
    //    /* concatenate all the HSPLists from 'results' */
    //
    //    for (i = 0; i < results->num_queries; i++) {
    //
    //        BlastHitList *hitlist = results->hitlist_array[i];
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:268-281
    // ```c++
    //  * @param hsp_list_out The read HSP list. [out]
    //  * @return Success, error, or end of reading, when nothing left to read.
    //  */
    // int BlastHSPStreamRead(BlastHSPStream* hsp_stream, BlastHSPList** hsp_list_out)
    // {
    //    *hsp_list_out = NULL;
    //
    //    if (!hsp_stream)
    //       return kBlastHSPStream_Error;
    //
    //    if (!hsp_stream->results)
    //       return kBlastHSPStream_Eof;
    //
    //    /* If this stream is not yet closed for writing, close it. In particular,
    // ```
    pub fn read_observed(&mut self, trace: &mut dyn ResultsObserver) -> Option<HspList> {
        self.close_observed(trace);
        let list = self.sorted.as_mut().unwrap().pop();
        trace.read(list.as_ref());
        list
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:599-613
    // ```c++
    //    hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1];
    //    target_oid = hsplist->oid;
    //
    //    for (i = 0; i < num_hsplists; i++) {
    //        hsplist = hsp_stream->sorted_hsplists[num_hsplists - 1 - i];
    //        if (hsplist->oid != target_oid)
    //            break;
    //
    //        batch->hsplist_array[i] = hsplist;
    //    }
    //
    //    hsp_stream->num_hsplists = num_hsplists - i;
    //    batch->num_hsplists = i;
    //
    //    return kBlastHSPStream_Success;
    // ```
    pub fn batch_read(&mut self) -> Vec<HspList> {
        self.batch_read_observed(&mut NoopResults)
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:145-164
    // ```c++
    //        if (hsp_stream->sort_by_score->sort_on_read) {
    //            Blast_HSPResultsReverseSort(hsp_stream->results);
    //        } else {
    //            /* Reverse the order of HSP lists, because they will be returned
    //               starting from end, for the sake of convenience */
    //            Blast_HSPResultsReverseOrder(hsp_stream->results);
    //        }
    //        hsp_stream->results_sorted = TRUE;
    //        hsp_stream->x_lock = MT_LOCK_Delete(hsp_stream->x_lock);
    //        return;
    //    }
    //
    //    results = hsp_stream->results;
    //    num_hsplists = hsp_stream->num_hsplists;
    //
    //    /* concatenate all the HSPLists from 'results' */
    //
    //    for (i = 0; i < results->num_queries; i++) {
    //
    //        BlastHitList *hitlist = results->hitlist_array[i];
    // ```
    pub fn batch_read_observed(&mut self, trace: &mut dyn ResultsObserver) -> Vec<HspList> {
        self.close_observed(trace);
        let lists = self.sorted.as_mut().unwrap();
        let Some(oid) = lists.last().map(|l| l.oid) else {
            trace.batch_read(&[]);
            return Vec::new();
        };
        let mut batch = Vec::new();
        while lists.last().is_some_and(|l| l.oid == oid) {
            batch.push(lists.pop().unwrap());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:606-615
        // ```c++
        //
        //        batch->hsplist_array[i] = hsplist;
        //    }
        //
        //    hsp_stream->num_hsplists = num_hsplists - i;
        //    batch->num_hsplists = i;
        //
        //    return kBlastHSPStream_Success;
        // }
        //
        // ```
        trace.batch_read(&batch);
        batch
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2911-3003
// ```c++
//    }
//    else {            /* query seq is split */
//
//       /* An HSP can be a candidate for merging if it lies in the
//          overlap region. Whether this is true depends on whether the
//          HSP starts to the left of the split point, or ends to the
//          right of the overlap region. A complication is that 'left'
//          and 'right' have opposite meaning when the HSP is on the
//          minus strand of the query sequence */
//
//       for (index1 = 0; index1 < combined_hsp_list->hspcnt; index1++) {
//          hsp1 = combined_hsp_list->hsp_array[index1];
//          offset_idx = hsp1->context % contexts_per_query;
//          if (split_offsets[offset_idx] < 0) continue;
//          if ((hsp1->query.frame >= 0 && hsp1->query.end >
//                          split_offsets[offset_idx]) ||
//              (hsp1->query.frame < 0 && hsp1->query.offset <
//                          split_offsets[offset_idx] + chunk_overlap_size)) {
//             /* At least part of this HSP lies in the overlap strip. */
//             hsp_var = combined_hsp_list->hsp_array[hspcnt1];
//             combined_hsp_list->hsp_array[hspcnt1] = hsp1;
//             combined_hsp_list->hsp_array[index1] = hsp_var;
//             ++hspcnt1;
//          }
//       }
//       for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
//          hsp2 = hsp_list->hsp_array[index2];
//          offset_idx = hsp2->context % contexts_per_query;
//          if (split_offsets[offset_idx] < 0) continue;
//          if ((hsp2->query.frame < 0 && hsp2->query.end >
//                          split_offsets[offset_idx]) ||
//              (hsp2->query.frame >= 0 && hsp2->query.offset <
//                          split_offsets[offset_idx] + chunk_overlap_size)) {
//             /* At least part of this HSP lies in the overlap strip. */
//             hsp_var = hsp_list->hsp_array[hspcnt2];
//             hsp_list->hsp_array[hspcnt2] = hsp2;
//             hsp_list->hsp_array[index2] = hsp_var;
//             ++hspcnt2;
//          }
//       }
//    }
//
//    /* the merge process is independent of whether merging happens
//       between query chunks or subject chunks */
//
//    if (hspcnt1 > 0 && hspcnt2 > 0) {
//       hspp1 = combined_hsp_list->hsp_array;
//       hspp2 = hsp_list->hsp_array;
//
//       for (index1 = 0; index1 < hspcnt1; index1++) {
//
//          hsp1 = hspp1[index1];
//
//          for (index2 = 0; index2 < hspcnt2; index2++) {
//
//             hsp2 = hspp2[index2];
//
//             /* Skip already deleted HSPs, or HSPs from different contexts */
//             if (!hsp2 || hsp1->context != hsp2->context)
//                continue;
//
//             /* Short read qureies are shorter than the overlap region and may
//                already have a traceback */
//             if (short_reads) {
//                 hspp2[index2] = Blast_HSPFree(hsp2);
//                 continue;
//             }
//
//             /* we have to determine the starting diagonal of one HSP
//                and the ending diagonal of the other */
//
//             if (contexts_per_query < 0 || hsp1->query.frame >= 0) {
//                end_diag = s_HSPEndDiag(hsp1);
//                start_diag = s_HSPStartDiag(hsp2);
//             }
//             else {
//                end_diag = s_HSPEndDiag(hsp2);
//                start_diag = s_HSPStartDiag(hsp1);
//             }
//
//             if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
//                if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) {
//                   /* Free the second HSP. */
//                   hspp2[index2] = Blast_HSPFree(hsp2);
//                }
//             }
//          }
//       }
//
//       /* Purge the nulled out HSPs from the new HSP list */
//       Blast_HSPListPurgeNullHSPs(hsp_list);
//    }
//
// ```
fn merge_hsps(combined: &mut Vec<Hsp>, mut incoming: Vec<Hsp>, points: &[i32; 6]) {
    let mut nold = 0;
    for i in 0..combined.len() {
        let h = &combined[i].hsp;
        let p = points[h.context % 6];
        if p >= 0
            && if h.frame >= 0 {
                h.q_end > p
            } else {
                h.q_start < p + OVERLAP_AA
            }
        {
            combined.swap(nold, i);
            nold += 1;
        }
    }
    let mut nnew = 0;
    for i in 0..incoming.len() {
        let h = &incoming[i].hsp;
        let p = points[h.context % 6];
        if p >= 0
            && if h.frame < 0 {
                h.q_end > p
            } else {
                h.q_start < p + OVERLAP_AA
            }
        {
            incoming.swap(nnew, i);
            nnew += 1;
        }
    }
    let mut removed = vec![false; incoming.len()];
    for a in &mut combined[..nold] {
        for (i, b) in incoming[..nnew].iter().enumerate() {
            if removed[i] || a.hsp.context != b.hsp.context {
                continue;
            }
            let ah = &a.hsp;
            let bh = &b.hsp;
            let (end, start) = if ah.frame >= 0 {
                (ah.q_end - ah.s_end, bh.q_start - bh.s_start)
            } else {
                (bh.q_end - bh.s_end, ah.q_start - ah.s_start)
            };
            if (end - start).abs() < 10 && super::split::merge_two(&mut a.hsp, &b.hsp) {
                removed[i] = true;
            }
        }
    }
    combined.extend(
        incoming
            .into_iter()
            .enumerate()
            .filter_map(|(i, h)| (!removed[i]).then_some(h)),
    );
    combined.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101
// ```c++
// s_CompoHeapRecordCompare(BlastCompo_HeapRecord * place1,
//                          BlastCompo_HeapRecord * place2)
// {
//     int result;
//     if (0 == (result = CMP(place1->bestEvalue, place2->bestEvalue)) &&
//         0 == (result = CMP(place2->bestScore, place1->bestScore))) {
//         result = CMP(place2->subject_index, place1->subject_index);
//     }
//     return result > 0;
// }
//
//
// /** Swap two records in the heap. */
// ```
pub struct HeapEntry {
    pub list: HspList,
    pub best_score: i32,
}
fn worse(a: &HeapEntry, b: &HeapEntry) -> bool {
    a.list
        .best_evalue
        .partial_cmp(&b.list.best_evalue)
        .unwrap()
        .then_with(|| b.best_score.cmp(&a.best_score))
        .then_with(|| b.list.oid.cmp(&a.list.oid))
        == Ordering::Greater
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:402-419
// ```c++
//
// /* Documented in compo_heap.h. */
// int
// BlastCompo_HeapFilledToCutoff(const BlastCompo_Heap * self)
// {
//     return self->n >= self->heapThreshold &&
//         self->worstEvalue <= self->ecutoff;
// }
//
//
// /* Documented in compo_heap.h. */
// int
// BlastCompo_HeapInitialize(BlastCompo_Heap * self, int heapThreshold,
//                           double ecutoff)
// {
//     self->n             = 0;
//     self->heapThreshold = heapThreshold;
//     self->ecutoff       = ecutoff;
// ```
pub struct CompoHeap {
    pub entries: Vec<HeapEntry>,
    threshold: usize,
    cutoff: f64,
    pub worst: f64,
    pub heapified: bool,
}
impl CompoHeap {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:163-163
    // ```c++
    // #define PSI_INCLUSION_ETHRESH 0.002 /**< Inclusion threshold for PSI BLAST */
    // ```
    pub fn new(threshold: usize) -> Self {
        Self {
            entries: Vec::new(),
            threshold,
            cutoff: 0.002,
            worst: 0.0,
            heapified: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:163-245
    // ```c++
    //  * @param n            the size of the entire heap.
    //  */
    // static void
    // s_CompoHeapifyDown(BlastCompo_HeapRecord * heapArray,
    //                        int top, int n)
    // {
    //     int i, left, right, largest;    /* placeholders for indices in swapping */
    //
    //     largest = top;
    //     do {
    //         i = largest;
    //         left  = 2 * i;
    //         right = 2 * i + 1;
    //         if (left <= n &&
    //             s_CompoHeapRecordCompare(&heapArray[left],
    //                                      &heapArray[i])) {
    //             largest = left;
    //         } else {
    //             largest = i;
    //         }
    //         if (right <= n &&
    //             s_CompoHeapRecordCompare(&heapArray[right],
    //                                      &heapArray[largest])) {
    //             largest = right;
    //         }
    //         if (largest != i) {
    //             s_CompoHeapRecordSwap(&heapArray[i], &heapArray[largest]);
    //         }
    //     } while (largest != i);
    //     if (COMPO_INTENSE_DEBUG) {
    //         assert(s_CompoHeapIsValid(heapArray, top, n));
    //     }
    // }
    //
    //
    // /**
    //  * Relocate a leaf in the heap so that the entire heap is in valid
    //  * heap order.  On entry, all elements but the leaf must be in valid
    //  * heap order.
    //  *
    //  * @param heapArray      array representing the heap as a binary tree
    //  * @param i              element in heap array that may be out of order [in]
    //  */
    // static void
    // s_CompoHeapifyUp(BlastCompo_HeapRecord * heapArray, int i)
    // {
    //     int parent = i / 2;          /* index to the node that is the
    //                                     parent of node i */
    //     while (parent >= 1 && s_CompoHeapRecordCompare(&heapArray[i],
    //                                                    &heapArray[parent]))
    //     {
    //         s_CompoHeapRecordSwap(&heapArray[i], &heapArray[parent]);
    //
    //         i       = parent;
    //         parent /= 2;
    //     }
    //     if (COMPO_INTENSE_DEBUG) {
    //         assert(s_CompoHeapIsValid(heapArray, 1, i));
    //     }
    // }
    //
    //
    // /** Convert a BlastCompo_Heap from a representation as an unordered array to
    //  *  a representation as a heap-ordered array.
    //  *
    //  *  @param self         the BlastCompo_Heap to convert
    //  */
    // static void
    // s_ConvertToHeap(BlastCompo_Heap * self)
    // {
    //     if (NULL != self->array) {    /* If we aren't already a heap */
    //         int i;                     /* heap node index */
    //         int n;                     /* number of elements in the heap */
    //         self->heapArray = self->array;
    //         self->array     = NULL;
    //
    //         n = self->n;
    //         for (i = n / 2;  i >= 1;  --i) {
    //             s_CompoHeapifyDown(self->heapArray, i, n);
    //         }
    //     }
    //     if (COMPO_INTENSE_DEBUG) {
    //         assert(s_CompoHeapIsValid(self->heapArray, 1, self->n));
    // ```
    fn down(&mut self, mut i: usize) {
        loop {
            let l = i * 2 + 1;
            let r = l + 1;
            let mut large = i;
            if l < self.entries.len() && worse(&self.entries[l], &self.entries[i]) {
                large = l;
            }
            if r < self.entries.len() && worse(&self.entries[r], &self.entries[large]) {
                large = r;
            }
            if large == i {
                break;
            }
            self.entries.swap(i, large);
            i = large;
        }
    }
    fn up(&mut self, mut i: usize) {
        while i > 0 {
            let p = (i - 1) / 2;
            if !worse(&self.entries[i], &self.entries[p]) {
                break;
            }
            self.entries.swap(i, p);
            i = p;
        }
    }
    fn heapify(&mut self) {
        if !self.heapified {
            self.heapified = true;
            for i in (0..self.entries.len() / 2).rev() {
                self.down(i);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    // ```c++
    // BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
    //                            double eValue,
    //                            int score,
    //                            int subject_index)
    // {
    //     if (self->n < self->heapThreshold ||
    //         eValue <= self->ecutoff ||
    //         eValue <  self->worstEvalue) {
    //         return TRUE;
    //     } else {
    //         /* self is either currently a heap, or must be converted to
    //          * one; use s_CompoHeapRecordCompare to compare against
    //          * the worst element in the heap */
    //         BlastCompo_HeapRecord heapRecord; /* temporary record to
    //                                              compare against */
    //         if (self->heapArray == NULL) s_ConvertToHeap(self);
    //
    //         heapRecord.bestEvalue       = eValue;
    //         heapRecord.bestScore        = score;
    //         heapRecord.subject_index    = subject_index;
    //         heapRecord.theseAlignments  = NULL;
    //
    //         return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
    //     }
    // ```
    pub fn would_insert(&mut self, e: f64, score: i32, oid: usize) -> bool {
        if self.entries.len() < self.threshold || e <= self.cutoff || e < self.worst {
            return true;
        }
        self.heapify();
        let root = &self.entries[0];
        root.list
            .best_evalue
            .partial_cmp(&e)
            .unwrap()
            .then_with(|| score.cmp(&root.best_score))
            .then_with(|| oid.cmp(&root.list.oid))
            == Ordering::Greater
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391
    // ```c++
    // BlastCompo_HeapInsert(BlastCompo_Heap * self,
    //                       void * alignments,
    //                       double eValue,
    //                       int score,
    //                       int subject_index,
    //                       void ** discardedAlignments)
    // {
    //     *discardedAlignments = NULL;
    //     if (self->array && self->n >= self->heapThreshold) {
    //         s_ConvertToHeap(self);
    //     }
    //     if (self->array != NULL) {
    //         /* "self" is currently a list. Add the new alignments to the end */
    //         int status =
    //             s_CompHeapRecordInsertAtEnd(&self->array, &self->n,
    //                                         &self->capacity, alignments,
    //                                         eValue, score,
    //                                         subject_index);
    //         if (status != 0) { /* out of memory */
    //             return -1;
    //         }
    //         if (self->worstEvalue < eValue) {
    //             self->worstEvalue = eValue;
    //         }
    //     } else {                      /* "self" is currently a heap */
    //         if (self->n < self->heapThreshold ||
    //             (eValue <= self->ecutoff &&
    //              self->worstEvalue <= self->ecutoff)) {
    //             /* The new alignments must be inserted into the heap, and all old
    //              * alignments retained */
    //             int status =
    //                 s_CompHeapRecordInsertAtEnd(&self->heapArray,
    //                                             &self->n,
    //                                             &self->capacity,
    //                                             alignments, eValue,
    //                                             score, subject_index);
    //             if (status != 0) { /* out of memory */
    //                 return -1;
    //             }
    //             s_CompoHeapifyUp(self->heapArray, self->n);
    //         } else {
    //             /* Some set of alignments must be discarded; discardedAlignments
    //              * will hold a pointer to these alignments. */
    //             BlastCompo_HeapRecord heapRecord;   /* Candidate record
    //                                                    for insertion */
    //             heapRecord.bestEvalue      = eValue;
    //             heapRecord.bestScore       = score;
    //             heapRecord.theseAlignments = alignments;
    //             heapRecord.subject_index   = subject_index;
    //
    //             if (s_CompoHeapRecordCompare(&self->heapArray[1],
    //                                              &heapRecord)) {
    //                 /* The new record should be inserted, and the largest
    //                  * element currently in the heap may be discarded */
    //                 *discardedAlignments = self->heapArray[1].theseAlignments;
    //                 memcpy(&self->heapArray[1], &heapRecord,
    //                        sizeof(BlastCompo_HeapRecord));
    //             } else {
    //                 *discardedAlignments = heapRecord.theseAlignments;
    //             }
    //             s_CompoHeapifyDown(self->heapArray, 1, self->n);
    //         }
    // ```
    pub fn insert(&mut self, entry: HeapEntry) -> Option<HeapEntry> {
        if !self.heapified && self.entries.len() >= self.threshold {
            self.heapify();
        }
        if !self.heapified {
            self.worst = self.worst.max(entry.list.best_evalue);
            self.entries.push(entry);
            return None;
        }
        let discarded = if self.entries.len() < self.threshold
            || (entry.list.best_evalue <= self.cutoff && self.worst <= self.cutoff)
        {
            self.entries.push(entry);
            self.up(self.entries.len() - 1);
            None
        } else if worse(&self.entries[0], &entry) {
            let d = std::mem::replace(&mut self.entries[0], entry);
            self.down(0);
            Some(d)
        } else {
            self.down(0);
            Some(entry)
        };
        self.worst = self.entries[0].list.best_evalue;
        discarded
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:439-461
    // ```c++
    // }
    //
    //
    // /* Documented in compo_heap.h. */
    // void *
    // BlastCompo_HeapPop(BlastCompo_Heap * self)
    // {
    //     void * results = NULL;   /* the list of SeqAligns to be returned */
    //
    //     s_ConvertToHeap(self);
    //     if (self->n > 0) { /* The heap is not empty */
    //         BlastCompo_HeapRecord *first, *last; /* The first and last
    //                                                 elements of the array
    //                                                 that represents the
    //                                                 heap.  */
    //         first = &self->heapArray[1];
    //         last  = &self->heapArray[self->n];
    //
    //         results = first->theseAlignments;
    //         if (--self->n > 0) {
    //             /* The heap is still not empty */
    //             memcpy(first, last, sizeof(BlastCompo_HeapRecord));
    //             s_CompoHeapifyDown(self->heapArray, 1, self->n);
    // ```
    pub fn pop(&mut self) -> Option<HeapEntry> {
        self.heapify();
        if self.entries.is_empty() {
            return None;
        }
        let entry = self.entries.swap_remove(0);
        if !self.entries.is_empty() {
            self.down(0);
        }
        Some(entry)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1582
    // ```c++
    // BlastCompo_EarlyTermination(double evalue,
    //                             BlastCompo_Heap significantMatches[],
    //                             int numQueries)
    // {
    //     int i;
    //     for (i = 0;  i < numQueries;  i++) {
    //         if (BlastCompo_HeapFilledToCutoff(&significantMatches[i])) {
    //             double ecutoff = significantMatches[i].ecutoff;
    //             /* Only matches with evalue <= ethresh will be saved. */
    //             if (evalue <= EVALUE_STRETCH * ecutoff) {
    //                 /* The evalue if this match is sufficiently small
    //                  * that we want to redo it to try to obtain an
    //                  * alignment with evalue smaller than ecutoff. */
    //                 return FALSE;
    //             }
    //         } else {
    //             return FALSE;
    //         }
    //     }
    //     return TRUE;
    // }
    // ```
    pub fn early(e: f64, heaps: &[Self]) -> bool {
        heaps
            .iter()
            .all(|h| h.entries.len() >= h.threshold && h.worst <= h.cutoff && e > 5.0 * h.cutoff)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
// ```c++
// BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
//                            double eValue,
//                            int score,
//                            int subject_index)
// {
//     if (self->n < self->heapThreshold ||
//         eValue <= self->ecutoff ||
//         eValue <  self->worstEvalue) {
//         return TRUE;
//     } else {
//         /* self is either currently a heap, or must be converted to
//          * one; use s_CompoHeapRecordCompare to compare against
//          * the worst element in the heap */
//         BlastCompo_HeapRecord heapRecord; /* temporary record to
//                                              compare against */
//         if (self->heapArray == NULL) s_ConvertToHeap(self);
//
//         heapRecord.bestEvalue       = eValue;
//         heapRecord.bestScore        = score;
//         heapRecord.subject_index    = subject_index;
//         heapRecord.theseAlignments  = NULL;
//
//         return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
//     }
// ```
#[cfg(test)]
mod boundary_tests {
    use super::*;
    use std::fmt::Write;
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:75-85
    // ```c++
    // typedef struct BlastCompo_HeapRecord {
    //     double        bestEvalue;     /**< best (smallest) evalue of all
    //                                        alignments in the record */
    //     int           bestScore;      /**< best (largest) score; used to
    //                                        break ties between records with
    //                                        the same e-value */
    //     int           subject_index;  /**< index of the subject sequence in
    //                                        the database */
    //     void *        theseAlignments;  /**< a collection of alignments */
    // } BlastCompo_HeapRecord;
    //
    // ```
    fn dump(out: &mut String, stage: &str, oid: i32, result: i32, h: &CompoHeap) {
        write!(
            out,
            "H\t{stage}\t{oid}\t{result}\t{}\t{:016x}\t{}",
            h.entries.len(),
            h.worst.to_bits(),
            i32::from(h.heapified)
        )
        .unwrap();
        for e in &h.entries {
            write!(
                out,
                "\t{}:{}:{:016x}",
                e.list.oid,
                e.best_score,
                e.list.best_evalue.to_bits()
            )
            .unwrap();
        }
        out.push('\n');
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    // ```c++
    // BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
    //                            double eValue,
    //                            int score,
    //                            int subject_index)
    // {
    //     if (self->n < self->heapThreshold ||
    //         eValue <= self->ecutoff ||
    //         eValue <  self->worstEvalue) {
    //         return TRUE;
    //     } else {
    //         /* self is either currently a heap, or must be converted to
    //          * one; use s_CompoHeapRecordCompare to compare against
    //          * the worst element in the heap */
    //         BlastCompo_HeapRecord heapRecord; /* temporary record to
    //                                              compare against */
    //         if (self->heapArray == NULL) s_ConvertToHeap(self);
    //
    //         heapRecord.bestEvalue       = eValue;
    //         heapRecord.bestScore        = score;
    //         heapRecord.subject_index    = subject_index;
    //         heapRecord.theseAlignments  = NULL;
    //
    //         return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391
    // ```c++
    // BlastCompo_HeapInsert(BlastCompo_Heap * self,
    //                       void * alignments,
    //                       double eValue,
    //                       int score,
    //                       int subject_index,
    //                       void ** discardedAlignments)
    // {
    //     *discardedAlignments = NULL;
    //     if (self->array && self->n >= self->heapThreshold) {
    //         s_ConvertToHeap(self);
    //     }
    //     if (self->array != NULL) {
    //         /* "self" is currently a list. Add the new alignments to the end */
    //         int status =
    //             s_CompHeapRecordInsertAtEnd(&self->array, &self->n,
    //                                         &self->capacity, alignments,
    //                                         eValue, score,
    //                                         subject_index);
    //         if (status != 0) { /* out of memory */
    //             return -1;
    //         }
    //         if (self->worstEvalue < eValue) {
    //             self->worstEvalue = eValue;
    //         }
    //     } else {                      /* "self" is currently a heap */
    //         if (self->n < self->heapThreshold ||
    //             (eValue <= self->ecutoff &&
    //              self->worstEvalue <= self->ecutoff)) {
    //             /* The new alignments must be inserted into the heap, and all old
    //              * alignments retained */
    //             int status =
    //                 s_CompHeapRecordInsertAtEnd(&self->heapArray,
    //                                             &self->n,
    //                                             &self->capacity,
    //                                             alignments, eValue,
    //                                             score, subject_index);
    //             if (status != 0) { /* out of memory */
    //                 return -1;
    //             }
    //             s_CompoHeapifyUp(self->heapArray, self->n);
    //         } else {
    //             /* Some set of alignments must be discarded; discardedAlignments
    //              * will hold a pointer to these alignments. */
    //             BlastCompo_HeapRecord heapRecord;   /* Candidate record
    //                                                    for insertion */
    //             heapRecord.bestEvalue      = eValue;
    //             heapRecord.bestScore       = score;
    //             heapRecord.theseAlignments = alignments;
    //             heapRecord.subject_index   = subject_index;
    //
    //             if (s_CompoHeapRecordCompare(&self->heapArray[1],
    //                                              &heapRecord)) {
    //                 /* The new record should be inserted, and the largest
    //                  * element currently in the heap may be discarded */
    //                 *discardedAlignments = self->heapArray[1].theseAlignments;
    //                 memcpy(&self->heapArray[1], &heapRecord,
    //                        sizeof(BlastCompo_HeapRecord));
    //             } else {
    //                 *discardedAlignments = heapRecord.theseAlignments;
    //             }
    //             s_CompoHeapifyDown(self->heapArray, 1, self->n);
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:439-461
    // ```c++
    // }
    //
    //
    // /* Documented in compo_heap.h. */
    // void *
    // BlastCompo_HeapPop(BlastCompo_Heap * self)
    // {
    //     void * results = NULL;   /* the list of SeqAligns to be returned */
    //
    //     s_ConvertToHeap(self);
    //     if (self->n > 0) { /* The heap is not empty */
    //         BlastCompo_HeapRecord *first, *last; /* The first and last
    //                                                 elements of the array
    //                                                 that represents the
    //                                                 heap.  */
    //         first = &self->heapArray[1];
    //         last  = &self->heapArray[self->n];
    //
    //         results = first->theseAlignments;
    //         if (--self->n > 0) {
    //             /* The heap is still not empty */
    //             memcpy(first, last, sizeof(BlastCompo_HeapRecord));
    //             s_CompoHeapifyDown(self->heapArray, 1, self->n);
    // ```
    #[test]
    fn pinned_c_heap_threshold_cutoff_tie_and_empty_pop() {
        let expected = include_str!("../../../tests/unit/blastx_d_boundary_expected.tsv")
            .lines()
            .filter(|l| !l.starts_with("C\t"))
            .map(|l| format!("{l}\n"))
            .collect::<String>();
        let mut actual = String::new();
        for threshold in 1..=2 {
            let mut heap = CompoHeap::new(threshold);
            writeln!(actual, "T\t{threshold}").unwrap();
            dump(&mut actual, "INIT", -1, 0, &heap);
            let es = [0.1, 0.1, 0.2, 0.09, 0.0021, 0.002, 0.001, 0.002, 0.003];
            let scores = [100, 100, 101, 99, 99, 98, 97, 98, 1000];
            for (i, (&e, &score)) in es.iter().zip(&scores).enumerate() {
                let yes = heap.would_insert(e, score, i);
                dump(&mut actual, "WOULD", i as i32, i32::from(yes), &heap);
                if yes {
                    let discarded = heap.insert(HeapEntry {
                        list: HspList {
                            query_index: 0,
                            oid: i,
                            hsps: Vec::new(),
                            best_evalue: e,
                        },
                        best_score: score,
                    });
                    dump(
                        &mut actual,
                        "INSERT",
                        i as i32,
                        discarded.map_or(-1, |e| e.list.oid as i32),
                        &heap,
                    );
                }
            }
            loop {
                let entry = heap.pop();
                dump(
                    &mut actual,
                    "POP",
                    entry.as_ref().map_or(-1, |e| e.list.oid as i32),
                    0,
                    &heap,
                );
                if entry.is_none() {
                    break;
                }
            }
        }
        assert_eq!(actual, expected);
    }
}
