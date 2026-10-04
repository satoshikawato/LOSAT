//! NCBI's merge of the preliminary hit lists of BLASTP query chunks (`BlastHSPStreamMerge`):
//! the chunks themselves are `common/protein_query_split.rs`, and each chunk's preliminary
//! stage is a search of its query parts (`blast_engine.rs`).

use std::cmp::Ordering;

use super::hsp::{score_compare_hsps, BlastpHitList, BlastpHspList};

/// One query part of a chunk: the batch query it belongs to (the index of its preliminary
/// hit list among the batch's), its offset in that query and that query's length.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BlastpChunkPart {
    pub batch_query: usize,
    pub q_idx: u32,
    pub offset: i32,
    pub query_length: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hspstream.c:476-534
// ```c
//        for (j = 0; j < contexts_per_query; j++) {
//            Int4 local_context = i * contexts_per_query + j;
//            if (context_list[local_context] >= 0) {
//                split_points[context_list[local_context] % contexts_per_query] =
//                                 offset_list[local_context];
//            }
//        }
//        ...
//                hsp->context = context_list[local_context];
//                hsp->query.offset += offset_list[local_context];
//                hsp->query.end += offset_list[local_context];
//                hsp->query.gapped_start += offset_list[local_context];
//                hsp->query.frame = BLAST_ContextToFrame(stream2->program,
//                                                        hsp->context);
//        ...
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
//    }
// ```
/// The merge of one query chunk (RP-4): `chunk_lists[i]` is the preliminary hit list of the
/// chunk's part `parts[i]`; its HSPs move to the part's query (a blastp query has one
/// context, `query_context` its index, and frame 0) and the list merges with that query's hit
/// list of the chunks before, `merged[parts[i].batch_query]`. The split blocks of a gapped
/// search allow gaps (split_query_cxx.cpp:864-865, `GetGappedMode`).
pub(crate) fn merge_query_chunk(
    merged: &mut [Option<BlastpHitList>],
    chunk_lists: Vec<Option<BlastpHitList>>,
    parts: &[BlastpChunkPart],
    chunk_overlap_size: i32,
) {
    for (hitlist, part) in chunk_lists.into_iter().zip(parts) {
        let Some(mut hitlist) = hitlist else {
            continue;
        };
        let offset = part.offset as usize;
        for list in hitlist.hsplist_array.iter_mut().take(hitlist.hsplist_count) {
            for hsp in &mut list.hsps {
                hsp.query_context = part.q_idx as i32;
                hsp.q_idx = part.q_idx;
                hsp.query_length = part.query_length;
                hsp.q_start += offset;
                hsp.q_end += offset;
                hsp.gapped_q_start += part.offset;
            }
            list.query_index = part.q_idx;
        }
        hit_list_merge(
            hitlist,
            &mut merged[part.batch_query],
            part.offset,
            chunk_overlap_size,
        );
    }
    for hitlist in merged.iter_mut().flatten() {
        let count = hitlist.hsplist_count;
        for list in hitlist.hsplist_array.iter_mut().take(count) {
            sort_by_score(list);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2132-2217
// ```c
//     if (hitlist1 == NULL)
//         return 0;
//     if (hitlist2 == NULL) {
//         *combined_hit_list_ptr = hitlist1;
//         *old_hit_list_ptr = NULL;
//         return 0;
//     }
//     num_hsplists1 = hitlist1->hsplist_count;
//     num_hsplists2 = hitlist2->hsplist_count;
//     new_hitlist = Blast_HitListNew(hitlist1->hsplist_max);
//
//     /* sort the lists of HSPs by oid */
//
//     if (num_hsplists1 > 1) {
//         qsort(hitlist1->hsplist_array, num_hsplists1,
//               sizeof(BlastHSPList*), s_SortHSPListByOid);
//     }
//     ...
//     query_is_split = FALSE;
//     for (i = 0; i < contexts_per_query; i++) {
//         if (split_offsets[i] > 0) {
//             query_is_split = TRUE;
//     ...
//         if (hsplist1->oid < hsplist2->oid) {
//             Blast_HitListUpdate(new_hitlist, hsplist1);
//             i++;
//         }
//         else if (hsplist1->oid > hsplist2->oid) {
//             Blast_HitListUpdate(new_hitlist, hsplist2);
//             j++;
//         }
//         else {
//             ...
//             if (query_is_split) {
//                 Blast_HSPListsMerge(hitlist1->hsplist_array + i,
//                                     hitlist2->hsplist_array + j,
//                                     hsplist2->hsp_max, split_offsets,
//                                     contexts_per_query,
//                                     chunk_overlap_size,
//                                     allow_gap, FALSE);
//             }
//             else {
//                 Blast_HSPListAppend(hitlist1->hsplist_array + i,
//                                     hitlist2->hsplist_array + j,
//                                     hsplist2->hsp_max);
//             }
//             Blast_HitListUpdate(new_hitlist, hitlist2->hsplist_array[j]);
//     ...
//     *old_hit_list_ptr = NULL;
//     *combined_hit_list_ptr = new_hitlist;
// ```
// The subjects of a hit list are distinct, so the OID sort is total. BLASTP has no
// `-max_hsps` merge limit before the traceback (`hsp_max` keeps every HSP).
fn hit_list_merge(
    mut hitlist1: BlastpHitList,
    combined: &mut Option<BlastpHitList>,
    split_offset: i32,
    chunk_overlap_size: i32,
) {
    let Some(mut hitlist2) = combined.take() else {
        *combined = Some(hitlist1);
        return;
    };
    let mut new_hitlist = BlastpHitList::new(hitlist1.hsplist_max);
    let mut lists1: Vec<BlastpHspList> = hitlist1
        .hsplist_array
        .drain(..)
        .take(hitlist1.hsplist_count)
        .collect();
    let mut lists2: Vec<BlastpHspList> = hitlist2
        .hsplist_array
        .drain(..)
        .take(hitlist2.hsplist_count)
        .collect();
    lists1.sort_by_key(|list| list.oid);
    lists2.sort_by_key(|list| list.oid);
    let query_is_split = split_offset > 0;
    let mut lists1 = lists1.into_iter().peekable();
    let mut lists2 = lists2.into_iter().peekable();
    while let (Some(list1), Some(list2)) = (lists1.peek(), lists2.peek()) {
        match list1.oid.cmp(&list2.oid) {
            Ordering::Less => new_hitlist.update(lists1.next().expect("peeked")),
            Ordering::Greater => new_hitlist.update(lists2.next().expect("peeked")),
            Ordering::Equal => {
                let list1 = lists1.next().expect("peeked");
                let mut list2 = lists2.next().expect("peeked");
                if query_is_split {
                    query_hsp_lists_merge(list1, &mut list2, split_offset, chunk_overlap_size);
                } else {
                    combine_by_score(list1, &mut list2);
                }
                new_hitlist.update(list2);
            }
        }
    }
    for list in lists1.chain(lists2) {
        new_hitlist.update(list);
    }
    *combined = Some(new_hitlist);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2921-3029
// ```c
//       for (index1 = 0; index1 < combined_hsp_list->hspcnt; index1++) {
//          hsp1 = combined_hsp_list->hsp_array[index1];
//          offset_idx = hsp1->context % contexts_per_query;
//          if (split_offsets[offset_idx] < 0) continue;
//          if ((hsp1->query.frame >= 0 && hsp1->query.end >
//                          split_offsets[offset_idx]) ||
//          ...
//             hsp_var = combined_hsp_list->hsp_array[hspcnt1];
//             combined_hsp_list->hsp_array[hspcnt1] = hsp1;
//             combined_hsp_list->hsp_array[index1] = hsp_var;
//             ++hspcnt1;
//       ...
//          if ((hsp2->query.frame < 0 && hsp2->query.end >
//                          split_offsets[offset_idx]) ||
//              (hsp2->query.frame >= 0 && hsp2->query.offset <
//                          split_offsets[offset_idx] + chunk_overlap_size)) {
//       ...
//             /* Skip already deleted HSPs, or HSPs from different contexts */
//             if (!hsp2 || hsp1->context != hsp2->context)
//                continue;
//       ...
//             if (contexts_per_query < 0 || hsp1->query.frame >= 0) {
//                end_diag = s_HSPEndDiag(hsp1);
//                start_diag = s_HSPStartDiag(hsp2);
//             }
//       ...
//             if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
//                if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) {
//                   /* Free the second HSP. */
//                   hspp2[index2] = Blast_HSPFree(hsp2);
//       ...
//    new_hspcnt =
//       MIN(hsp_list->hspcnt + combined_hsp_list->hspcnt, hsp_num_max);
//    ...
//    s_BlastHSPListsCombineByScore(hsp_list, combined_hsp_list, new_hspcnt);
// ```
// The query split branch of `Blast_HSPListsMerge` for a protein query (one context, frame
// 0), with gaps allowed. `BlastpHsp` keeps `query.offset + 1` in `q_start` and `query.end`
// in `q_end` (likewise for the subject).
fn query_hsp_lists_merge(
    new: BlastpHspList,
    combined: &mut BlastpHspList,
    split_offset: i32,
    chunk_overlap_size: i32,
) {
    if new.hsps.is_empty() {
        return;
    }
    let mut incoming = new.hsps;
    let old = &mut combined.hsps;
    let mut hspcnt1 = 0;
    let mut hspcnt2 = 0;
    if split_offset >= 0 {
        for index1 in 0..old.len() {
            if old[index1].q_end as i64 > i64::from(split_offset) {
                old.swap(hspcnt1, index1);
                hspcnt1 += 1;
            }
        }
        for index2 in 0..incoming.len() {
            if (incoming[index2].q_start as i64 - 1)
                < i64::from(split_offset) + i64::from(chunk_overlap_size)
            {
                incoming.swap(hspcnt2, index2);
                hspcnt2 += 1;
            }
        }
    }
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut slots: Vec<Option<super::hsp::BlastpHsp>> =
            incoming.into_iter().map(Some).collect();
        for hsp1 in old.iter_mut().take(hspcnt1) {
            for slot in slots.iter_mut().take(hspcnt2) {
                let Some(hsp2) = slot.as_ref() else {
                    continue;
                };
                if hsp1.query_context != hsp2.query_context {
                    continue;
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1463-1478
                // ```c
                // s_HSPStartDiag(const BlastHSP *hsp)
                // {
                //     return hsp->query.offset - hsp->subject.offset;
                // }
                // ...
                // s_HSPEndDiag(const BlastHSP *hsp)
                // {
                //     return hsp->query.end - hsp->subject.end;
                // }
                // ```
                let end_diag = hsp1.q_end as i64 - hsp1.s_end as i64;
                let start_diag = (hsp2.q_start as i64 - 1) - (hsp2.s_start as i64 - 1);
                // NCBI blast_hits.c:1537: #define OVERLAP_DIAG_CLOSE 10
                if (end_diag - start_diag).abs() < 10 && merge_two_hsps(hsp1, hsp2) {
                    *slot = None;
                }
            }
        }
        incoming = slots.into_iter().flatten().collect();
    }
    combined.hsps.extend(incoming);
    sort_by_score(combined);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1488-1530
// ```c
//    /* do not merge off-diagonal hsps for ungapped search */
//    if (!allow_gap &&
//        hsp1->subject.offset - hsp2->subject.offset -hsp1->query.offset + hsp2->query.offset)
//    {
//        return FALSE;
//    }
//
//    if(hsp1->subject.frame != hsp2->subject.frame)
// 	   return FALSE;
//
//    /* combine the boundaries of the two HSPs,
//       assuming they intersect at all */
//    if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end,
//                         hsp2->query.offset,
//                         hsp1->subject.offset, hsp1->subject.end,
//                         hsp2->subject.offset) ||
//        CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end,
//                         hsp2->query.end,
//                         hsp1->subject.offset, hsp1->subject.end,
//                         hsp2->subject.end)) {
//
// 	  double score_density =  (hsp1->score + hsp2->score) *(1.0) /
// 			                  ((hsp1->query.end - hsp1->query.offset) +
// 			                   (hsp2->query.end - hsp2->query.offset));
//       hsp1->query.offset = MIN(hsp1->query.offset, hsp2->query.offset);
//       hsp1->subject.offset = MIN(hsp1->subject.offset, hsp2->subject.offset);
//       hsp1->query.end = MAX(hsp1->query.end, hsp2->query.end);
//       hsp1->subject.end = MAX(hsp1->subject.end, hsp2->subject.end);
//       if (hsp2->score > hsp1->score) {
//           hsp1->query.gapped_start = hsp2->query.gapped_start;
//           hsp1->subject.gapped_start = hsp2->subject.gapped_start;
// 	  hsp1->score = hsp2->score;
//       }
//
//       hsp1->score = MAX((int) (score_density *(hsp1->query.end - hsp1->query.offset)), hsp1->score);
//       return TRUE;
//    }
// ```
// `CONTAINED_IN_HSP(a,b,c,d,e,f)` is `a <= c && b >= c && d <= f && e >= f`
// (blast_hits.h). The gaps are allowed, and `length` follows the merged query span.
fn merge_two_hsps(hsp1: &mut super::hsp::BlastpHsp, hsp2: &super::hsp::BlastpHsp) -> bool {
    if hsp1.subject_frame != hsp2.subject_frame {
        return false;
    }
    let (q1_offset, q1_end) = (hsp1.q_start as i64 - 1, hsp1.q_end as i64);
    let (s1_offset, s1_end) = (hsp1.s_start as i64 - 1, hsp1.s_end as i64);
    let (q2_offset, q2_end) = (hsp2.q_start as i64 - 1, hsp2.q_end as i64);
    let (s2_offset, s2_end) = (hsp2.s_start as i64 - 1, hsp2.s_end as i64);
    let contained = |q: i64, s: i64| q1_offset <= q && q1_end >= q && s1_offset <= s && s1_end >= s;
    if !contained(q2_offset, s2_offset) && !contained(q2_end, s2_end) {
        return false;
    }
    let score_density = f64::from(hsp1.raw_score + hsp2.raw_score)
        / ((q1_end - q1_offset) + (q2_end - q2_offset)) as f64;
    let q_offset = q1_offset.min(q2_offset);
    let s_offset = s1_offset.min(s2_offset);
    let q_end = q1_end.max(q2_end);
    let s_end = s1_end.max(s2_end);
    hsp1.q_start = (q_offset + 1) as usize;
    hsp1.s_start = (s_offset + 1) as usize;
    hsp1.q_end = q_end as usize;
    hsp1.s_end = s_end as usize;
    hsp1.length = (q_end - q_offset) as usize;
    if hsp2.raw_score > hsp1.raw_score {
        hsp1.gapped_q_start = hsp2.gapped_q_start;
        hsp1.gapped_s_start = hsp2.gapped_s_start;
        hsp1.raw_score = hsp2.raw_score;
    }
    hsp1.raw_score = ((score_density * (q_end - q_offset) as f64) as i32).max(hsp1.raw_score);
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2758-2766,2819-2829,2849
// ```c
//    if (new_hspcnt >= hsp_list->hspcnt + combined_hsp_list->hspcnt) {
//       /* All HSPs from both arrays are saved */
//       for (index=combined_hsp_list->hspcnt, index1=0;
//            index1<hsp_list->hspcnt; index1++) {
//          if (hsp_list->hsp_array[index1] != NULL)
//             combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
//       }
//       combined_hsp_list->hspcnt = new_hspcnt;
//       Blast_HSPListSortByScore(combined_hsp_list);
//    ...
//    /* If no previous HSP list, return a pointer to the old one */
//    if (!combined_hsp_list) {
//       *combined_hsp_list_ptr = hsp_list;
//       *old_hsp_list_ptr = NULL;
//       return 0;
//    }
//    ...
//    s_BlastHSPListsCombineByScore(hsp_list, combined_hsp_list, new_hspcnt);
// ```
// `Blast_HSPListAppend` of a query that is not split here.
fn combine_by_score(new: BlastpHspList, combined: &mut BlastpHspList) {
    combined.hsps.extend(new.hsps);
    sort_by_score(combined);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1383
// ```c
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!hsp_list || hsp_list->hspcnt <= 1)
//         return;
//
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
//     }
// }
// ```
// Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5): HSPs that two chunks found
// with other gapped starts keep their order.
fn sort_by_score(list: &mut BlastpHspList) {
    if list
        .hsps
        .windows(2)
        .any(|pair| score_compare_hsps(&pair[0], &pair[1]) == Ordering::Greater)
    {
        list.hsps.sort_by(score_compare_hsps);
    }
}
