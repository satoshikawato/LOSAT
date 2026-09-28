//! Comparison-only observation of internal runtime; no public BLASTX execution/report success.
use anyhow::{bail, Result};
use LOSAT::{
    algorithm::blastx::{
        input::read_fasta,
        runtime::{search_internal_observed, RuntimeObserver},
    },
    cli::{Cli, Commands},
    common::GapEditOp,
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1785
// ```c++
//     if(hit_params->options->query_cov_hsp_perc > 0 || hit_params->options->max_hsps_per_subject > 0 ||
//        (hit_params->options->hsp_filt_opt != NULL && hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
//     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
//     }
//
//     /* Re-sort the hit lists according to their best e-values, because they
//        could have changed. Only do this for a database search. */
//     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
//         Blast_HSPResultsSortByEvalue(results);
//     }
//
//
//     /* Eliminate extra hits from results, if preliminary hit list size is
//        larger than the final hit list size */
//     s_BlastPruneExtraHits(results, hit_params->options->hitlist_size);
//
//     if (retval == BLASTERR_INTERRUPTED) {
//         results = Blast_HSPResultsFree(results);
//     }
//
//     *results_out = results;
//
//     return retval;
// ```

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1239-1274
// ```c++
//                             in_align, FALSE, subject_is_translated);
//                     adjust_search_failed =
//                             Blast_AdjustScores(matrix, query_composition,
//                                     query.length,
//                                     &subject_composition,
//                                     subject.length,
//                                     scaledMatrixInfo, compo_adjust_mode,
//                                     RE_pseudocounts, NRrecord,
//                                     &matrix_adjust_rule,
//                                     callbacks->calc_lambda,
//                                     pvalueForThisPair,
//                                     compositionTestIndex,
//                                     LambdaRatio);
//                     if (adjust_search_failed < 0) { /* fatal error */
//                         status = adjust_search_failed;
//                         goto window_index_loop_cleanup;
//                     }
//                     num_adjustments++;
//                 }
//
//                 if ( !adjust_search_failed ) {
//                     newAlign = callbacks->redo_one_alignment(
//                             in_align,
//                             matrix_adjust_rule,
//                             &query,
//                             &window->query_range,
//                             ccat_query_length,
//                             &subject,
//                             &window->subject_range,
//                             matchingSeq->length,
//                             gapping_params
//                     );
//                     if (newAlign && newAlign->score >= params->cutoff_s) {
//                         s_WithDistinctEnds(&newAlign, &alignments[query_index],
//                                 callbacks->free_align_traceback,
//                                 num_adjustments == 1);
// ```
struct Observer {
    context: usize,
    before: u64,
    adjust_count: usize,
}
fn matrix_hash(
    m: Option<&LOSAT::core::composition_adjustment::adjust_scores::AdjustedProteinMatrix>,
) -> u64 {
    let mut h = 14695981039346656037u64;
    for q in 0..28 {
        for s in 0..28 {
            let v = m.map_or_else(
                || LOSAT::utils::matrix::blosum62_score_ncbistdaa_direct(q as u8, s as u8),
                |m| m.scores[q][s],
            );
            h ^= v as u32 as u64;
            h = h.wrapping_mul(1099511628211);
        }
    }
    h
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:53-62
// ```c++
//     BlastHSP * hsp;
//     Int4 cid;    /* context id for hsp */
//     Int4 sid;    /* OID for hsp*/
//     Int4 begin;  /* query offset in plus strand */
//     Int4 end;    /* query end in plus strand */
//     Int4 merit;  /* how many other hsps in the tree dominates me? */
//     struct LinkedHSP *next;
// } LinkedHSP;
//
// /** functions to manipulate LinkedHSPs */
// ```
fn inline_hsp(oid: usize, h: &LOSAT::algorithm::blastx::preliminary::PreliminaryHsp) {
    print!(
        "\t{oid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        h.context,
        h.frame,
        h.score,
        h.q_start,
        h.q_end,
        h.q_gapped_start,
        h.s_start,
        h.s_end,
        h.s_gapped_start
    );
}
fn dump_set(
    stage: &str,
    oid: usize,
    hsps: &[LOSAT::algorithm::blastx::preliminary::PreliminaryHsp],
) {
    println!("B_SET\t{stage}\t{oid}\t{}", hsps.len());
    for (i, h) in hsps.iter().enumerate() {
        println!(
            "B_HSP\t{stage}\t{oid}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            h.context,
            h.frame,
            h.score,
            h.q_start,
            h.q_end,
            h.q_gapped_start,
            h.s_start,
            h.s_end,
            h.s_gapped_start
        );
    }
}
impl LOSAT::algorithm::blastx::results::ResultsObserver for Observer {
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
    fn targets(
        &mut self,
        stage: &str,
        q: usize,
        lists: &[LOSAT::algorithm::blastx::results::HspList],
    ) {
        print!("B_TARGETS\t{stage}\t{q}\t{}", lists.len());
        for l in lists {
            print!(
                "\t{}:{}:{:016x}",
                l.oid,
                l.hsps.len(),
                l.best_evalue.to_bits()
            );
        }
        println!();
    }

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
    fn merge(&mut self, c: &LOSAT::algorithm::blastx::split::QueryChunk, enter: bool) {
        print!(
            "B_MERGE\t{}\t{}",
            i32::from(enter),
            c.absolute_contexts.len()
        );
        for (context, offset) in c.absolute_contexts.iter().zip(&c.corrections) {
            print!("\t{context}:{offset}");
        }
        println!();
    }
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
    fn read(&mut self, l: Option<&LOSAT::algorithm::blastx::results::HspList>) {
        if let Some(l) = l {
            println!("B_READ\t0\t{}\t{}", l.query_index, l.oid);
        } else {
            println!("B_READ\t1\t-1\t-1");
        }
    }
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
    fn batch_read(&mut self, ls: &[LOSAT::algorithm::blastx::results::HspList]) {
        print!("B_BATCHREAD\t{}\t{}", i32::from(ls.is_empty()), ls.len());
        for l in ls {
            print!("\t{}:{}", l.query_index, l.oid);
        }
        println!();
    }

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
    fn cull_compare(
        &mut self,
        ao: usize,
        a: &LOSAT::algorithm::blastx::results::Hsp,
        bo: usize,
        b: &LOSAT::algorithm::blastx::results::Hsp,
        result: bool,
    ) {
        print!("B_CULL_COMPARE");
        inline_hsp(ao, &a.hsp);
        inline_hsp(bo, &b.hsp);
        println!("\t{}", u8::from(result));
    }
    fn cull_delete(&mut self, oid: usize, h: &LOSAT::algorithm::blastx::results::Hsp, merit: i32) {
        print!("B_CULL_DELETE");
        inline_hsp(oid, &h.hsp);
        println!("\t{merit}");
    }
    fn cull_save(&mut self, oid: usize, h: &LOSAT::algorithm::blastx::results::Hsp) {
        print!("B_CULL_SAVE_IN");
        inline_hsp(oid, &h.hsp);
        println!();
    }
    fn cull_saved(&mut self, result: bool) {
        println!("B_CULL_SAVE_OUT\t{}", u8::from(result));
    }
    fn cull_fork(&mut self, begin: i32, end: i32) {
        println!("B_CULL_FORK\t{begin}\t{end}");
    }
    fn cull_final(&mut self, n: usize) {
        println!("B_CULL_FINAL\t{n}");
    }
    fn preliminary_call(&mut self, oid: usize) {
        println!("B_PRE_SOURCE_OID\t{oid}");
    }
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
    fn reap(
        &mut self,
        stage: &str,
        oid: usize,
        l: &LOSAT::algorithm::blastx::linking::LinkedHspList,
        cutoff: f64,
    ) {
        self.linked(stage, oid, l);
        println!(
            "B_REAP_STATS\t{stage}\t{oid}\t{:016x}\t{:016x}",
            cutoff.to_bits(),
            l.best_evalue.to_bits()
        );
        for (i, h) in l.hsps.iter().enumerate() {
            println!(
                "B_REAP_HSP\t{stage}\t{oid}\t{i}\t{}\t{:016x}",
                h.num,
                h.evalue.to_bits()
            );
        }
    }
    fn linked(
        &mut self,
        stage: &str,
        oid: usize,
        h: &LOSAT::algorithm::blastx::linking::LinkedHspList,
    ) {
        dump_set(
            stage,
            oid,
            &h.hsps.iter().map(|h| h.hsp.clone()).collect::<Vec<_>>(),
        );
    }
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
    fn containment(
        &mut self,
        h: &LOSAT::algorithm::blastx::preliminary::PreliminaryHsp,
        contained: bool,
    ) {
        println!(
            "B_CONTAINS\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            h.context,
            h.frame,
            h.score,
            h.q_start,
            h.q_end,
            h.s_start,
            h.s_end,
            i32::from(contained)
        );
    }
    fn numeric(
        &mut self,
        stage: &str,
        oid: usize,
        h: &[LOSAT::algorithm::blastx::preliminary::PreliminaryHsp],
    ) {
        dump_set(stage, oid, h);
    }
    fn hsps(&mut self, stage: &str, oid: usize, h: &[LOSAT::algorithm::blastx::results::Hsp]) {
        dump_set(
            stage,
            oid,
            &h.iter().map(|h| h.hsp.clone()).collect::<Vec<_>>(),
        );
    }
    fn hitlist(
        &mut self,
        stage: &str,
        oid: usize,
        max: usize,
        lists: &[LOSAT::algorithm::blastx::results::HspList],
        candidate: Option<&LOSAT::algorithm::blastx::results::HspList>,
    ) {
        println!("B_HITLIST\t{stage}\t{oid}\t{max}\t{}", lists.len());
        for list in lists {
            self.hsps(stage, list.oid, &list.hsps);
        }
        if let Some(c) = candidate {
            self.hsps("CANDIDATE", c.oid, &c.hsps);
        }
    }
}
impl RuntimeObserver for Observer {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    // ```c++
    //         // Restore the full query sequence for the traceback stage!
    //         if (m_InternalData->m_Queries == NULL) {
    //             CRef<ILocalQueryData> query_data
    //                 (m_QueryFactory->MakeLocalQueryData(&*m_Options));
    //             // Query masking info is calculated as a side-effect
    //             CBlastScoreBlk sbp
    //                 (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
    //                                         NULL, m_Messages, NULL, NULL));
    //             m_InternalData->m_Queries = query_data->GetSequenceBlk();
    //         }
    // ```
    fn restored(
        &mut self,
        b: &LOSAT::algorithm::blastx::query_setup::PreparedQueryBatch,
        p: &[LOSAT::algorithm::blastx::parameters::ContextParameters],
    ) {
        println!(
            "B_RESTORED\t{}\t{}\t{}",
            b.sequence_start.len() - 2,
            b.first_context,
            b.contexts.len() - 1
        );
        for (i, (c, p)) in b.contexts.iter().zip(p).enumerate() {
            println!(
                "B_RESTORED_CONTEXT\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.query_index,
                c.frame,
                c.offset,
                c.length,
                i32::from(p.valid),
                p.search_space,
                p.length_adjustment
            );
        }
        print!("B_RESTORED_BUFFER\tMASKED\t");
        for x in &b.sequence_start {
            print!("{x:02x}");
        }
        println!();
        print!("B_RESTORED_BUFFER\tNOMASK\t");
        for x in &b.sequence_start_nomask {
            print!("{x:02x}");
        }
        println!();
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
    fn heap(
        &mut self,
        q: usize,
        stage: &str,
        e: f64,
        score: i32,
        oid: i64,
        decision: i64,
        h: &LOSAT::algorithm::blastx::results::CompoHeap,
    ) {
        print!(
            "B_HEAP_{stage}\t{q}\t{:016x}\t{score}\t{oid}\t{decision}\t{}\t{:016x}\t{}",
            e.to_bits(),
            h.entries.len(),
            h.worst.to_bits(),
            u8::from(h.heapified)
        );
        for r in &h.entries {
            print!(
                "\t{}:{}:{:016x}",
                r.list.oid,
                r.best_score,
                r.list.best_evalue.to_bits()
            );
        }
        println!();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
    // ```c++
    // {
    //     Int4 index;
    //
    //     if (!hsp_list || hsp_list->hspcnt == 0)
    //         return 0;
    //
    //     ASSERT(link_hsp_params);
    //
    //     /* Remove any information on number of linked HSPs from previous
    //        linking. */
    //     for (index = 0; index < hsp_list->hspcnt; ++index)
    //         hsp_list->hsp_array[index]->num = 1;
    //
    //     /* Link up the HSP's for this hsp_list. */
    //     if (link_hsp_params->longest_intron <= 0) {
    //         s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
    //                               subject_length, sbp, link_hsp_params,
    //                               gapped_calculation);
    //         /* The HSP's may be in a different order than they were before,
    //            but hsp contains the first one. */
    //     } else {
    //         Blast_HSPListAdjustOddBlastnScores(hsp_list, gapped_calculation, sbp);
    //         /* Calculate individual HSP e-values first - they'll be needed to
    //            compare with sum e-values. Use decay rate to compensate for
    //            multiple tests. */
    // ```
    fn link_setup(
        &mut self,
        length: i32,
        db: i64,
        state: ([i32; 2], f64),
        link: &LOSAT::algorithm::blastx::linking::LinkParameters,
    ) {
        println!(
            "B_LINK_SETUP\t{length}\t{db}\t{}\t{}\t{:016x}\t{}\t{}\t{}\t{:016x}",
            state.0[0],
            state.0[1],
            state.1.to_bits(),
            link.longest_intron,
            link.gap_size,
            link.overlap_size,
            link.gap_decay_rate.to_bits()
        );
    }
    fn link_run(
        &mut self,
        length: i32,
        gapped: bool,
        link: &LOSAT::algorithm::blastx::linking::LinkParameters,
    ) {
        println!(
            "B_LINK_RUN\t{length}\t{}\t{}\t{}\t{}\t{:016x}",
            u8::from(gapped),
            link.longest_intron,
            link.gap_size,
            link.overlap_size,
            link.gap_decay_rate.to_bits()
        );
    }
    fn early(&mut self, e: f64, n: usize, decision: bool) {
        println!("B_EARLY\t{:016x}\t{n}\t{}", e.to_bits(), u8::from(decision));
    }

    fn matrix(
        &mut self,
        oid: usize,
        event: LOSAT::core::composition_adjustment::redo_alignment::BlastRedoTraceEvent<'_>,
    ) {
        use LOSAT::core::composition_adjustment::redo_alignment::BlastRedoTraceEvent::*;
        match event {
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:966-1000
// ```c++
// s_IsContained(BlastCompo_Alignment * in_align,
//               BlastCompo_Alignment * alignments,
//               double lambda)
// {
//     BlastCompo_Alignment * align;     /* represents the current alignment
//                                             in the main loop */
//     /* Endpoints of the alignment */
//     int query_offset    = in_align->queryStart;
//     int query_end       = in_align->queryEnd;
//     int subject_offset  = in_align->matchStart;
//     int subject_end     = in_align->matchEnd;
//     double score        = in_align->score;
//     double scoreThresh = score + KAPPA_BIT_TOL * LOCAL_LN2/lambda;
// 
//     for (align = alignments;  align != NULL;  align = align->next ) {
//         /* for all elements of alignments */
//         if (KAPPA_SIGN(in_align->frame) == KAPPA_SIGN(align->frame)) {
//             /* hsp1 and hsp2 are in the same query/subject frame */
//             if (KAPPA_CONTAINED_IN_HSP
//                 (align->queryStart, align->queryEnd, query_offset,
//                  align->matchStart, align->matchEnd, subject_offset) &&
//                 KAPPA_CONTAINED_IN_HSP
//                 (align->queryStart, align->queryEnd, query_end,
//                  align->matchStart, align->matchEnd, subject_end) &&
//                 scoreThresh <= align->score) {
//                 return 1;
//             }
//         }
//     }
//     return 0;
// }
// 
// 
// /* Documented in redo_alignment.h. */
// void
// ```
   Contained{incoming,result}=>println!("K_CONTAINED\t{oid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",incoming.query_index,incoming.score,incoming.query_start,incoming.query_end,incoming.match_start,incoming.match_end,i32::from(result)),
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:232-261
// ```c++
//     for (iread = 1;  iread < *hspcnt;  iread++) {
//         /* for all HSPs in the hitlist */
//         Int4      ireadBack;  /* iterator over indices less than iread */
//         BlastHSP *hsp1;       /* an HSP that is a candidate for deletion */
// 
//         hsp1 = hsp_array[iread];
//         for (ireadBack = 0;  ireadBack < iread && hsp1 != NULL;  ireadBack++) {
//             /* for all HSPs before hsp1 in the hitlist and while hsp1
//              * has not been deleted */
//             BlastHSP *hsp2;    /* an HSP that occurs earlier in hsp_array
//                                 * than hsp1 */
//             hsp2 = hsp_array[ireadBack];
// 
//             if( hsp2 == NULL ) {  /* hsp2 was deleted in a prior iteration. */
//                 continue;
//             }
//             if (hsp2->query.frame == hsp1->query.frame &&
//                 hsp2->subject.frame == hsp1->subject.frame) {
//                 /* hsp1 and hsp2 are in the same query/subject frame. */
//                 if (CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.offset) &&
//                     CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.end)    &&
//                     hsp1->score <= hsp2->score) {
//                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
//                 }
// ```
   ReapContained{frame,score,q_start,q_end,s_start,s_end,result}=>println!("K_REAP_CONTAINED\t{oid}\t{frame}\t{score}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t{}",i32::from(result)),
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:266-273
// ```c
//     /* Condense the hsp_array, removing any NULL items. */
//     iwrite = 0;
//     for (iread = 0;  iread < *hspcnt;  iread++) {
//         if (hsp_array[iread] != NULL) {
//             hsp_array[iwrite++] = hsp_array[iread];
//         }
//     }
//     *hspcnt = iwrite;
// ```
   ReapList{phase,index,context,frame,score,q_start,q_end,s_start,s_end}=>println!("X_REAP_LIST\t{oid}\t{phase}\t{index}\t{context}\t{frame}\t{score}\t{q_start}\t{q_end}\t{s_start}\t{s_end}"),
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:238-260
// ```c
//         for (ireadBack = 0;  ireadBack < iread && hsp1 != NULL;  ireadBack++) {
//             /* for all HSPs before hsp1 in the hitlist and while hsp1
//              * has not been deleted */
//             BlastHSP *hsp2;    /* an HSP that occurs earlier in hsp_array
//                                 * than hsp1 */
//             hsp2 = hsp_array[ireadBack];
// 
//             if( hsp2 == NULL ) {  /* hsp2 was deleted in a prior iteration. */
//                 continue;
//             }
//             if (hsp2->query.frame == hsp1->query.frame &&
//                 hsp2->subject.frame == hsp1->subject.frame) {
//                 /* hsp1 and hsp2 are in the same query/subject frame. */
//                 if (CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.offset) &&
//                     CONTAINED_IN_HSP
//                     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
//                      hsp2->subject.offset, hsp2->subject.end,
//                      hsp1->subject.end)    &&
//                     hsp1->score <= hsp2->score) {
//                     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
// ```
   ReapCompare{index,previous_index,result}=>println!("X_REAP_COMPARE\t{oid}\t{index}\t{previous_index}\t{}",i32::from(result)),
   MatchStart{context,count,matrix}=>{self.context=context;self.adjust_count=0;println!("K_MATCH\t{oid}\t{context}\t{count}\t{:016x}",matrix_hash(matrix));},
   MatchEnd{context,matrix}=>println!("K_MATCH_END\t{oid}\t{context}\t0\t{:016x}",matrix_hash(matrix)),
   Adjustment{before,status,query_count,subject_count,rule,matrix}=>{if before{self.before=matrix_hash(matrix);}else{self.adjust_count+=1;println!("K_ADJUST\t{oid}\t{}\t{status}\t{query_count}\t{subject_count}\t{}\t{:016x}\t{:016x}",self.context,rule as i32,self.before,matrix_hash(matrix));}},
   Redo{incoming,rule,matrix}=>{println!("K_REDO\t{oid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}",incoming.query_index,incoming.query_start,incoming.query_end,incoming.match_start,incoming.match_end,self.adjust_count,rule as i32,matrix_hash(matrix));self.adjust_count=0;},
  }
    }
}
fn main() -> Result<()> {
    let cli: Cli = LOSAT::cli::try_parse_from(std::env::args_os())?;
    let Commands::Blastx(args) = cli.command else {
        bail!("BLASTX internal runtime comparison only")
    };
    let options = args.resolve()?;
    let subjects = read_fasta(&args.subject, true, args.lcase_masking)?;
    let queries = read_fasta(&args.query, false, args.lcase_masking)?;
    for (b, batch) in search_internal_observed(
        &queries,
        &subjects,
        &options,
        &mut Observer {
            context: 0,
            before: 0,
            adjust_count: 0,
        },
    )?
    .into_iter()
    .enumerate()
    {
        println!("R_BATCH\t{b}\t{}", batch.queries.len());
        for (q, lists) in batch.queries.iter().enumerate() {
            println!("R_Q\t{b}\t{q}\t{}", lists.len());
            for (s, list) in lists.iter().enumerate() {
                println!("R_LIST\t{b}\t{q}\t{s}\t{}\t{}", list.oid, list.hsps.len());
                for (i, h) in list.hsps.iter().enumerate() {
                    let p = &h.hsp;
                    println!("R_HSP\t{b}\t{q}\t{s}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}\t{:016x}\t{}\t{}\t{}",p.context,p.frame,p.q_start,p.q_end,p.q_gapped_start,p.s_start,p.s_end,p.s_gapped_start,p.score,h.num,h.evalue.to_bits(),h.bit_score.to_bits(),h.identity,h.positive,h.composition_method);
                    print!("R_EDIT\t{b}\t{q}\t{s}\t{i}\t{}", h.edit_script.len());
                    for op in &h.edit_script {
                        let (kind, len) = match op {
                            GapEditOp::Sub(n) => (3, n),
                            GapEditOp::Ins(n) => (6, n),
                            GapEditOp::Del(n) => (0, n),
                        };
                        print!("\t{kind}:{len}");
                    }
                    println!();
                }
            }
        }
    }
    Ok(())
}
