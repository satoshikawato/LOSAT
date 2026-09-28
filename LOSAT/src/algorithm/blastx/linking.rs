//! BLASTX uneven-gap linking in internal context coordinates.
use super::{
    args::ResolvedOptions,
    parameters::ContextParameters,
    preliminary::{compare_score, PreliminaryHsp},
};
use crate::stats::{
    spouge::{blast_spouge_stoe, BlastGumbelBlk},
    sum_statistics::{gap_decay_divisor, uneven_gap_sum_e},
    tables::KarlinParams,
};
use anyhow::{ensure, Result};
use std::cmp::Ordering;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:598-617
// ```c++
// Int2 BlastLinkHSPParametersNew(EBlastProgramType program_number,
//                                Boolean gapped_calculation,
//                                BlastLinkHSPParameters** link_hsp_params)
// {
//    BlastLinkHSPParameters* params;
//
//    if (!link_hsp_params)
//       return -1;
//
//    params = (BlastLinkHSPParameters*)
//       calloc(1, sizeof(BlastLinkHSPParameters));
//
//    if (program_number == eBlastTypeBlastn || !gapped_calculation) {
//       params->gap_prob = BLAST_GAP_PROB;
//       params->gap_decay_rate = BLAST_GAP_DECAY_RATE;
//    } else {
//       params->gap_prob = BLAST_GAP_PROB_GAPPED;
//       params->gap_decay_rate = BLAST_GAP_DECAY_RATE_GAPPED;
//    }
//    params->gap_size = BLAST_GAP_SIZE;
// ```
#[derive(Clone, Debug)]
pub struct LinkParameters {
    pub gap_size: i32,
    pub overlap_size: i32,
    pub longest_intron: i32,
    pub gap_decay_rate: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:781-815
// ```c++
//    if (params->do_sum_stats) {
//       BlastLinkHSPParametersNew(program_number, gapped_calculation,
//                                 &params->link_hsp_params);
//
//       if((Blast_QueryIsTranslated(program_number) ||
// 	  Blast_SubjectIsTranslated(program_number)) &&
// 	 program_number != eBlastTypeTblastx) {
//           /* The program may use Blast_UnevenGapLinkHSPs find significant
//              collections of distinct alignments */
//           Int4 max_protein_gap; /* the largest gap permitted in the
//                                * translated sequence */
//
//           max_protein_gap = (options->longest_intron - 2)/3;
//           if(gapped_calculation) {
//               if(options->longest_intron == 0) {
//                   /* a zero value of longest_intron invokes the
//                    * default behavior, which for gapped calculation is
//                    * to set longest_intron to a predefined value. */
//                   params->link_hsp_params->longest_intron =
//                       (DEFAULT_LONGEST_INTRON - 2) / 3;
//               } else if(max_protein_gap <= 0) {
//                   /* A nonpositive value of max_protein_gap disables linking */
//                   params->link_hsp_params =
//                       BlastLinkHSPParametersFree(params->link_hsp_params);
//                   params->do_sum_stats = FALSE;
//               } else { /* the value of max_protein_gap is positive */
//                   params->link_hsp_params->longest_intron = max_protein_gap;
//               }
//           } else { /* This is an ungapped calculation. */
//               /* For ungapped calculations, we preserve the old behavior
//                * of the longest_intron parameter to maintain
//                * backward-compatibility with older versions of BLAST. */
//               params->link_hsp_params->longest_intron =
//                 MAX(max_protein_gap, 0);
//           }
// ```
pub fn link_parameters(options: &ResolvedOptions) -> Option<LinkParameters> {
    if !options.sum_stats {
        return None;
    }
    let max_gap = (options.max_intron_length - 2) / 3;
    let intron = if options.gapped {
        if options.max_intron_length == 0 {
            40
        } else if max_gap <= 0 {
            return None;
        } else {
            max_gap
        }
    } else {
        max_gap.max(0)
    };
    Some(LinkParameters {
        gap_size: 40,
        overlap_size: 9,
        longest_intron: intron,
        gap_decay_rate: if options.gapped { 0.1 } else { 0.5 },
    })
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
#[derive(Clone, Debug)]
pub struct LinkedHsp {
    pub hsp: PreliminaryHsp,
    pub num: i32,
    pub evalue: f64,
    pub source_index: usize,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1810
// ```c++
//     /* Sort the HSP array by score */
//     Blast_HSPListSortByScore(hsp_list);
//
//     /* Find and fill the best e-value */
//     hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
//     for (index = 1; index < hsp_list->hspcnt; ++index) {
//         if (hsp_list->hsp_array[index]->evalue < hsp_list->best_evalue)
//             hsp_list->best_evalue = hsp_list->hsp_array[index]->evalue;
//     }
// ```
#[derive(Clone, Debug)]
pub struct LinkedHspList {
    pub hsps: Vec<LinkedHsp>,
    pub best_evalue: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1093-1101
// ```c++
// /******************************************************************************
//  * Structures and functions used only in uneven gap linking method.           *
//  ******************************************************************************/
//
// /** Simple doubly linked list of HSPs, used for calculating sum statistics. */
// typedef struct BlastLinkedHSPSet {
//     BlastHSP* hsp;                 /**< HSP for the current link in the chain. */
//     Uint4  queryId;                /**< Used for support of OOF linking */
//     struct BlastLinkedHSPSet* next;/**< Next link in the chain. */
// ```
#[derive(Clone)]
struct WorkHsp {
    value: LinkedHsp,
    prev: Option<usize>,
    next: Option<usize>,
    sum_score: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1567-1594
// ```c++
// s_LinkedHSPSetArrayIndexQueryEnds(BlastLinkedHSPSet** hsp_array, Int4 hspcnt,
//                                   Int4** qend_index_ptr)
// {
//     Int4 index;
//     Int4* qend_index_array = NULL;
//     BlastLinkedHSPSet* link;
//     Int4 current_end = 0;
//     Int4 current_index = 0;
//
//     /* Allocate the array. */
//     *qend_index_ptr = qend_index_array = (Int4*) calloc(hspcnt, sizeof(Int4));
//     if (!qend_index_array)
//         return -1;
//
//     current_end = hsp_array[0]->hsp->query.end;
//
//     for (index = 1; index < hspcnt; ++index) {
//         link = hsp_array[index];
//         if (link->queryId > hsp_array[current_index]->queryId ||
//             link->hsp->query.end > current_end) {
//             current_index = index;
//             current_end = link->hsp->query.end;
//         }
//         qend_index_array[index] = current_index;
//     }
//     return 0;
// }
//
// ```
fn query_end_prefix_index(work: &[WorkHsp], offset_order: &[usize]) -> Vec<usize> {
    let mut index = vec![0; offset_order.len()];
    let mut current = 0;
    let mut current_end = work[offset_order[0]].value.hsp.q_end;
    for i in 1..offset_order.len() {
        let candidate = offset_order[i];
        let current_id = offset_order[current];
        if work[candidate].value.hsp.context / 3 > work[current_id].value.hsp.context / 3
            || work[candidate].value.hsp.q_end > current_end
        {
            current = i;
            current_end = work[candidate].value.hsp.q_end;
        }
        index[i] = current;
    }
    index
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1234-1263
// ```c++
//  *                  BlastLinkedHSPSet structures [in]
//  * @param size Number of elements in the array [in]
//  * @param queryId Context of the target HSP [in]
//  * @param offset The target offset to search for [in]
//  * @return The index in the array of the HSP whose start/end offset
//  *         is closest to but >= the value 'offset'
//  */
// static Int4
// s_HSPOffsetBinarySearch(BlastLinkedHSPSet** hsp_array, Int4 size,
//                         Uint4 queryId, Int4 offset)
// {
//    Int4 index, begin, end;
//
//    begin = 0;
//    end = size;
//    while (begin < end) {
//       index = (begin + end) / 2;
//
//       if (hsp_array[index]->queryId < queryId)
//           begin = index + 1;
//       else if (hsp_array[index]->queryId > queryId)
//           end = index;
//       else {
//           if (hsp_array[index]->hsp->query.offset >= offset)
//               end = index;
//           else
//               begin = index + 1;
//       }
//    }
//
// ```
fn first_start_at_or_after(
    work: &[WorkHsp],
    order: &[usize],
    context: usize,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, order.len());
    while begin < end {
        let middle = (begin + end) / 2;
        let h = &work[order[middle]].value;
        if h.hsp.context / 3 < context || (h.hsp.context / 3 == context && h.hsp.q_start < offset) {
            begin = middle + 1;
        } else {
            end = middle;
        }
    }
    end
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1282-1306
// ```c++
// {
//    Int4 begin, end;
//
//    begin = 0;
//    end = size;
//    while (begin < end) {
//        Int4 right_index = (begin + end) / 2;
//        Int4 left_index = qend_index_array[right_index];
//
//        if (hsp_array[right_index]->queryId < queryId)
//            begin = right_index + 1;
//        else if (hsp_array[right_index]->queryId > queryId)
//            end = left_index;
//        else {
//            if (hsp_array[left_index]->hsp->query.end >= offset)
//                end = left_index;
//            else
//                begin = right_index + 1;
//        }
//    }
//
//    return end;
// }
//
// /** Merges HSPs from two linked HSP sets into an array of HSPs, sorted in
// ```
fn first_end_at_or_after(
    work: &[WorkHsp],
    order: &[usize],
    end_index: &[usize],
    context: usize,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, order.len());
    while begin < end {
        let right = (begin + end) / 2;
        let left = end_index[right];
        let right_context = work[order[right]].value.hsp.context / 3;
        if right_context < context {
            begin = right + 1;
        } else if right_context > context {
            end = left;
        } else if work[order[left]].value.hsp.q_end >= offset {
            end = left;
        } else {
            begin = right + 1;
        }
    }
    end
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1325
// ```c++
//  * @param merged_size The total number of HSPs in two sets. [out]
//  * @return The array of pointers to HSPs representing a merged set.
//  */
// static BlastLinkedHSPSet**
// s_MergeLinkedHSPSets(BlastLinkedHSPSet* hsp_set1, BlastLinkedHSPSet* hsp_set2,
//                      Int4* merged_size)
// {
//     Int4 index;
//     Int4 length;
//     BlastLinkedHSPSet** merged_hsps;
//
//     /* Find the first link of the old HSP chain. */
//     while (hsp_set1->prev)
//         hsp_set1 = hsp_set1->prev;
//     /* Find first and last link in the new HSP chain. */
//     while (hsp_set2->prev)
// ```
fn chain_head(work: &[WorkHsp], mut id: usize) -> usize {
    while let Some(previous) = work[id].prev {
        id = previous;
    }
    id
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1354
// ```c++
//  * @param merged_size The total number of HSPs in two sets. [out]
//  * @return The array of pointers to HSPs representing a merged set.
//  */
// static BlastLinkedHSPSet**
// s_MergeLinkedHSPSets(BlastLinkedHSPSet* hsp_set1, BlastLinkedHSPSet* hsp_set2,
//                      Int4* merged_size)
// {
//     Int4 index;
//     Int4 length;
//     BlastLinkedHSPSet** merged_hsps;
//
//     /* Find the first link of the old HSP chain. */
//     while (hsp_set1->prev)
//         hsp_set1 = hsp_set1->prev;
//     /* Find first and last link in the new HSP chain. */
//     while (hsp_set2->prev)
//         hsp_set2 = hsp_set2->prev;
//
//     *merged_size = length = hsp_set1->hsp->num + hsp_set2->hsp->num;
//
//     merged_hsps = (BlastLinkedHSPSet**)
//         malloc(length*sizeof(BlastLinkedHSPSet*));
//
//     index = 0;
//     while (hsp_set1 || hsp_set2) {
//         /* NB: HSP sets for which some HSPs have identical query offsets cannot
//            possibly be admissible, so it doesn't matter how to deal with equal
//            offsets. */
//         if (!hsp_set2 || (hsp_set1 &&
//             hsp_set1->hsp->query.offset < hsp_set2->hsp->query.offset)) {
//             merged_hsps[index] = hsp_set1;
//             hsp_set1 = hsp_set1->next;
//         } else {
//             merged_hsps[index] = hsp_set2;
//             hsp_set2 = hsp_set2->next;
//         }
//         ++index;
//     }
//     return merged_hsps;
// }
//
// /** Combines two linked sets of HSPs into a single set.
//  * @param hsp_set1 First set of HSPs [in]
//  * @param hsp_set2 Second set of HSPs [in]
//  * @param sum_score The sum score of the combined linked set
// ```
fn merged_chain(work: &[WorkHsp], first: usize, second: usize) -> Vec<usize> {
    let (mut a, mut b) = (
        Some(chain_head(work, first)),
        Some(chain_head(work, second)),
    );
    let mut merged = Vec::new();
    while a.is_some() || b.is_some() {
        if b.is_none()
            || a.is_some_and(|id| work[id].value.hsp.q_start < work[b.unwrap()].value.hsp.q_start)
        {
            let id = a.unwrap();
            merged.push(id);
            a = work[id].next;
        } else {
            let id = b.unwrap();
            merged.push(id);
            b = work[id].next;
        }
    }
    merged
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1416-1492
// ```c++
//
//     if (!hsp_set1 || !hsp_set2 || !link_hsp_params)
//         return FALSE;
//
//     /* The first input HSP must be the head of its set. */
//     if (hsp_set1->prev)
//         return FALSE;
//
//     /* The second input HSP may not be the head of its set. Hence follow the
//        previous pointers to get to the head. */
//     for ( ; hsp_set2->prev; hsp_set2 = hsp_set2->prev);
//
//     /* If left and right HSP are the same, return inadmissible status. */
//     if (hsp_set1 == hsp_set2)
//         return FALSE;
//
//     /* Check if these HSPs are for the same protein sequence (same queryId) */
//     if (hsp_set1->queryId != hsp_set2->queryId)
//         return FALSE;
//
//     /* Check if new HSP and hsp_set2 are on the same nucleotide sequence strand.
//        (same sign of subject frame) */
//     if (SIGN(hsp_set1->hsp->subject.frame) !=
//         SIGN(hsp_set2->hsp->subject.frame))
//         return FALSE;
//
//     /* Merge the two sets into an array with increasing order of query
//        offsets. */
//     merged_hsps = s_MergeLinkedHSPSets(hsp_set1, hsp_set2, &combined_size);
//
//     gap_s = link_hsp_params->longest_intron; /* Maximal gap size in
//                                                          subject */
//     gap_q = link_hsp_params->gap_size; /* Maximal gap size in query */
//
//     overlap = link_hsp_params->overlap_size; /* Maximal overlap size in
//                                                      query or subject */
//
//     /* swap gap_s and gap_q if blastx */
//     if (program == eBlastTypeBlastx) {
//         gap_s = link_hsp_params->gap_size;
//         gap_q = link_hsp_params->longest_intron;
//     }
//
//     for (index = 0; index < combined_size - 1; ++index) {
//         BlastLinkedHSPSet* left_hsp = merged_hsps[index];
//         BlastLinkedHSPSet* right_hsp = merged_hsps[index+1];
//
//
//         /* If the new HSP is too far to the left from the right_hsp, indicate this by
//            setting the boolean output value to TRUE. */
//         if (left_hsp->hsp->query.end < right_hsp->hsp->query.offset - gap_q)
//             break;
//
//         /* Check if the left HSP's query offset is to the right of the right HSP's
//            offset, i.e. they came in wrong order. */
//         if (left_hsp->hsp->query.offset >= right_hsp->hsp->query.offset)
//             break;
//
//         /* Check the remaining condition for query offsets: left HSP cannot end
//            further than the maximal allowed overlap from the right HSP's offset;
//            and left HSP must end before the right HSP. */
//         if (left_hsp->hsp->query.end > right_hsp->hsp->query.offset + overlap ||
//             left_hsp->hsp->query.end >= right_hsp->hsp->query.end)
//             break;
//
//         /* Check the subject offsets conditions. */
//         if (left_hsp->hsp->subject.end >
//             right_hsp->hsp->subject.offset + overlap ||
//             left_hsp->hsp->subject.end <
//             right_hsp->hsp->subject.offset - gap_s ||
//             left_hsp->hsp->subject.offset >= right_hsp->hsp->subject.offset ||
//             left_hsp->hsp->subject.end >= right_hsp->hsp->subject.end)
//             break;
//     }
//
//     sfree(merged_hsps);
//
// ```
fn sets_admissible(work: &[WorkHsp], first: usize, second: usize, params: &LinkParameters) -> bool {
    if work[first].prev.is_some() {
        return false;
    }
    let second = chain_head(work, second);
    if first == second {
        return false;
    }
    if work[first].value.hsp.context / 3 != work[second].value.hsp.context / 3 {
        return false;
    }
    let merged = merged_chain(work, first, second);
    for pair in merged.windows(2) {
        let left = &work[pair[0]].value.hsp;
        let right = &work[pair[1]].value.hsp;
        if left.q_end < right.q_start - params.longest_intron
            || left.q_start >= right.q_start
            || left.q_end > right.q_start + params.overlap_size
            || left.q_end >= right.q_end
            || left.s_end > right.s_start + params.overlap_size
            || left.s_end < right.s_start - params.gap_size
            || left.s_start >= right.s_start
            || left.s_end >= right.s_end
        {
            return false;
        }
    }
    true
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1114-1155
// ```c++
//  * @return E-value of all the HSPs together
//  */
// static double
// s_SumHSPEvalue(EBlastProgramType program_number,
//    const BlastQueryInfo* query_info, Int4 subject_length,
//    const BlastLinkHSPParameters* link_hsp_params,
//    BlastLinkedHSPSet* head_hsp, BlastLinkedHSPSet* new_hsp, double* sum_score)
// {
//    double gap_decay_rate, sum_evalue;
//    Int2 num;
//    Int4 subject_eff_length, query_eff_length, len_adj;
//    Int4 context = head_hsp->hsp->context;
//    Int4 query_window_size;
//    Int4 subject_window_size;
//
//    ASSERT(program_number != eBlastTypeTblastx);
//
//    subject_eff_length = (Blast_SubjectIsTranslated(program_number)) ?
//        subject_length/3 : subject_length;
//
//    gap_decay_rate = link_hsp_params->gap_decay_rate;
//
//    num = head_hsp->hsp->num + new_hsp->hsp->num;
//
//    len_adj = query_info->contexts[context].length_adjustment;
//
//    query_eff_length = MAX(query_info->contexts[context].query_length - len_adj, 1);
//
//    subject_eff_length = MAX(subject_eff_length - len_adj, 1);
//
//    *sum_score = new_hsp->sum_score + head_hsp->sum_score;
//
//    query_window_size =
//       link_hsp_params->overlap_size + link_hsp_params->gap_size + 1;
//    subject_window_size =
//       link_hsp_params->overlap_size + link_hsp_params->longest_intron + 1;
//
//    sum_evalue =
//        BLAST_UnevenGapSumE(query_window_size, subject_window_size,
//           num, *sum_score, query_eff_length, subject_eff_length,
//           query_info->contexts[context].eff_searchsp,
//           BLAST_GapDecayDivisor(gap_decay_rate, num));
// ```
fn linked_sum_evalue(
    work: &[WorkHsp],
    head: usize,
    candidate: usize,
    query_lengths: &[i32],
    lengths: &[ContextParameters],
    subject_length: i32,
    params: &LinkParameters,
) -> (f64, f64) {
    let context = work[head].value.hsp.context;
    let adjustment = lengths[context].length_adjustment as i32;
    let query_eff_length = (query_lengths[context] - adjustment).max(1);
    let subject_eff_length = (subject_length - adjustment).max(1);
    let num = (work[head].value.num + work[candidate].value.num) as i16;
    let sum_score = work[candidate].sum_score + work[head].sum_score;
    let evalue = uneven_gap_sum_e(
        params.overlap_size + params.gap_size + 1,
        params.overlap_size + params.longest_intron + 1,
        num,
        sum_score,
        query_eff_length,
        subject_eff_length,
        lengths[context].search_space,
        gap_decay_divisor(params.gap_decay_rate, num as usize),
    );
    (evalue, sum_score)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1361-1404
// ```c++
// {
//     BlastLinkedHSPSet** merged_hsps;
//     BlastLinkedHSPSet* head_hsp;
//     Int4 index, new_num;
//
//     if (!hsp_set2)
//         return hsp_set1;
//     else if (!hsp_set1)
//         return hsp_set2;
//
//     merged_hsps = s_MergeLinkedHSPSets(hsp_set1, hsp_set2, &new_num);
//
//     head_hsp = merged_hsps[0];
//     head_hsp->prev = NULL;
//     for (index = 0; index < new_num; ++index) {
//         BlastLinkedHSPSet* link = merged_hsps[index];
//         if (index < new_num - 1) {
//             BlastLinkedHSPSet* next_link = merged_hsps[index+1];
//             link->next = next_link;
//             next_link->prev = link;
//         } else {
//             link->next = NULL;
//         }
//         link->sum_score = sum_score;
//         link->hsp->evalue = evalue;
//         link->hsp->num = new_num;
//     }
//
//     sfree(merged_hsps);
//     return head_hsp;
// }
//
// /** Checks if new candidate HSP is admissible to be linked to a set of HSPs on
//  * the left. The new HSP must start strictly before the parent HSP in both query
//  * and subject, and its end must lie within an interval from the parent HSP's
//  * start, determined by the allowed gap and overlap sizes in query and subject.
//  * This function also indicates whether parent is already too far to the right
//  * of the candidate HSP, via a boolean pointer.
//  * @param hsp_set1 First linked set of HSPs. [in]
//  * @param hsp_set2 Second linked set of HSPs. [in]
//  * @param link_hsp_params Parameters for linking HSPs. [in]
//  * @param program Type of BLAST program (blastx or tblastn) [in]
//  * @return Do the two sets satisfy the admissibility criteria to form a
//  *         combined set?
// ```
fn combine_sets(
    work: &mut [WorkHsp],
    first: usize,
    second: usize,
    sum_score: f64,
    evalue: f64,
) -> usize {
    let merged = merged_chain(work, first, second);
    let num = merged.len() as i32;
    for (position, &id) in merged.iter().enumerate() {
        work[id].prev = position.checked_sub(1).map(|i| merged[i]);
        work[id].next = merged.get(position + 1).copied();
        work[id].sum_score = sum_score;
        work[id].value.num = num;
        work[id].value.evalue = evalue;
    }
    merged[0]
}
// The public BLASTX path stays gated until C/D/E pass.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1627-1666
// ```c++
//    Int4* qend_index_array = NULL;
//
//    /* Check input arguments. */
//    if (!link_hsp_params || !sbp || !query_info)
//        return -1;
//
//    /* If HSP list is not available or has <= 1 HSPs, there is nothing to do. */
//    if (!hsp_list || hsp_list->hspcnt <= 1)
//        return 0;
//
//    if(gapped_calculation) {
//        kbp_array = sbp->kbp_gap;
//    } else {
//        kbp_array = sbp->kbp;
//    }
//
//    /* max gap size in query */
//    gap_size = (program == eBlastTypeBlastx) ?
//               link_hsp_params->longest_intron :
//               link_hsp_params->gap_size;
//
//    hspcnt = hsp_list->hspcnt;
//    hsp_array = hsp_list->hsp_array;
//
//    /* Set up an array of HSP structure wrappers. */
//    link_hsp_array =
//        s_LinkedHSPSetArraySetUp(hsp_array, hspcnt, kbp_array, program);
//
//    /* Allocate, fill and sort the auxiliary arrays. */
//    score_hsp_array =
//        (BlastLinkedHSPSet**) malloc(hspcnt*sizeof(BlastLinkedHSPSet*));
//    memcpy(score_hsp_array, link_hsp_array, hspcnt*sizeof(BlastLinkedHSPSet*));
//    qsort(score_hsp_array, hspcnt, sizeof(BlastLinkedHSPSet*),
//          s_SumScoreCompareLinkedHSPSets);
//    offset_hsp_array =
//        (BlastLinkedHSPSet**) malloc(hspcnt*sizeof(BlastLinkedHSPSet*));
//    memcpy(offset_hsp_array, link_hsp_array, hspcnt*sizeof(BlastLinkedHSPSet*));
//    qsort(offset_hsp_array, hspcnt, sizeof(BlastLinkedHSPSet*),
//          s_FwdCompareLinkedHSPSets);
//
// ```
pub fn link_uneven(
    input: &[PreliminaryHsp],
    query_lengths: &[i32],
    lengths: &[ContextParameters],
    subject_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    link: &LinkParameters,
) -> Result<LinkedHspList> {
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    ensure!(
        link.longest_intron > 0,
        "BLASTX uneven-gap link parameters required"
    );
    ensure!(
        query_lengths.len() == lengths.len() && lengths.len() == gapped_params.len(),
        "one statistical context per query"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1524-1540
    // ```c++
    //             kbp_array[hsp->context]->logK;
    //         link_hsp_array[index]->queryId =
    //             (program == eBlastTypeBlastx) ?
    //             hsp->context / 3 : hsp->context;
    //
    //         hsp_array[index]->num = 1;
    //     }
    //
    //     return link_hsp_array;
    // }
    //
    // /** Frees the array of special structures, used for linking HSPs and restores
    //  * the original contexts and subject/query order in BlastHSP structures, when
    //  * necessary.
    //  * @param link_hsp_array Array of wrapper HSP structures, used for linking. [in]
    //  * @param hspcnt Size of the array. [in]
    //  * @return NULL.
    // ```
    let mut work: Vec<WorkHsp> = input
        .iter()
        .enumerate()
        .map(|(source_index, hsp)| {
            let context = hsp.context;
            let hsp = hsp.clone();
            let kbp = &gapped_params[context];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1851-1896
            // ```c++
            //       kbp_context = hsp->context;
            //       if (RPS_prelim) {
            //           /* All kbp in preliminary stage are equivalent.  However, some
            //              may be invalid.  Search for the first populated kbp */
            //           int i;
            //           for (i=0; i < sbp->number_of_contexts; ++i) {
            //               if (kbp[i]) break;
            //           }
            //           kbp_context = i;
            //       }
            //       ASSERT(kbp[kbp_context]);
            //       kbp[kbp_context]->Lambda /= scaling_factor;
            //
            //       /* Round score down to even number for E-value calculations only. */
            //       /* Added 2018/5/16 by rackerst, SB-2303 */
            //       score = hsp->score;
            //       if (hsp_list  &&  hsp_list->hspcnt != 0
            //               &&  gapped_calculation  &&  sbp->round_down) {
            //           score &= ~1;
            //       }
            //
            //       if (sbp->gbp) {
            //           /* Only try Spouge's method if gumbel parameters are available */
            //           if (!isRPS) {
            //               hsp->evalue =
            //                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
            //                                query_info->contexts[hsp->context].query_length,
            //                                subject_length);
            //           } else {
            //               /* for RPS blast, query and subject is swapped */
            //               hsp->evalue =
            //                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
            //                                subject_length,
            //                                query_info->contexts[hsp->context].query_length);
            //           }
            //       } else {
            //           /* Get effective search space from the query information block */
            //           hsp->evalue =
            //               BLAST_KarlinStoE_simple(score, kbp[kbp_context],
            //                                query_info->contexts[hsp->context].eff_searchsp);
            //       }
            //
            //       hsp->evalue /= gap_decay_divisor;
            //       /* Put back the unscaled value of Lambda. */
            //       kbp[kbp_context]->Lambda *= scaling_factor;
            //    }
            // ```
            let evalue = if let Some(gumbel) = gumbel {
                blast_spouge_stoe(
                    hsp.score,
                    kbp,
                    gumbel,
                    query_lengths[context],
                    subject_length,
                )
            } else {
                lengths[context].search_space as f64
                    * (-kbp.lambda * hsp.score as f64 + kbp.k.ln()).exp()
            } / gap_decay_divisor(link.gap_decay_rate, 1);
            let sum_score = kbp.lambda * hsp.score as f64 - kbp.k.ln();
            WorkHsp {
                value: LinkedHsp {
                    hsp,
                    num: 1,
                    evalue,
                    source_index,
                },
                prev: None,
                next: None,
                sum_score,
            }
        })
        .collect();
    if work.len() > 1 {
        let mut score_order: Vec<usize> = (0..work.len()).collect();
        score_order.sort_by(|&a, &b| {
            if work[a].sum_score < work[b].sum_score {
                Ordering::Greater
            } else if work[a].sum_score > work[b].sum_score {
                Ordering::Less
            } else {
                compare_score(&work[a].value.hsp, &work[b].value.hsp)
            }
        });
        let mut offset_order: Vec<usize> = (0..work.len()).collect();
        offset_order.sort_by(|&a, &b| {
            work[a]
                .value
                .hsp
                .context
                .div_euclid(3)
                .cmp(&work[b].value.hsp.context.div_euclid(3))
                .then(work[a].value.hsp.q_start.cmp(&work[b].value.hsp.q_start))
                .then(work[a].value.hsp.s_start.cmp(&work[b].value.hsp.s_start))
        });
        let end_index = query_end_prefix_index(&work, &offset_order);
        let mut head = None;
        let mut index = 0;
        while index < score_order.len() {
            if head.is_none() {
                while index < score_order.len() {
                    let id = score_order[index];
                    if work[id].prev.is_none() && work[id].next.is_none() {
                        break;
                    }
                    index += 1;
                }
                if index == score_order.len() {
                    break;
                }
                head = Some(score_order[index]);
            }
            let current = head.unwrap();
            let mut tail = current;
            while let Some(next) = work[tail].next {
                tail = next;
            }
            let left_offset = work[current].value.hsp.q_start - link.longest_intron;
            let left = first_end_at_or_after(
                &work,
                &offset_order,
                &end_index,
                work[current].value.hsp.context / 3,
                left_offset,
            );
            let right = first_start_at_or_after(
                &work,
                &offset_order,
                work[tail].value.hsp.context / 3,
                work[tail].value.hsp.q_end + link.longest_intron,
            );
            let mut best_evalue = work[current].value.evalue;
            let mut best = None;
            let mut best_sum_score = 0.0;
            for &candidate in &offset_order[left..right] {
                if work[candidate]
                    .prev
                    .is_some_and(|prev| work[prev].value.hsp.q_end >= left_offset)
                {
                    continue;
                }
                if sets_admissible(&work, current, candidate, link) {
                    let (evalue, sum_score) = linked_sum_evalue(
                        &work,
                        current,
                        candidate,
                        query_lengths,
                        lengths,
                        subject_length,
                        link,
                    );
                    if evalue < best_evalue.min(work[candidate].value.evalue) {
                        best = Some(candidate);
                        best_evalue = evalue;
                        best_sum_score = sum_score;
                    }
                }
            }
            if let Some(candidate) = best {
                head = Some(combine_sets(
                    &mut work,
                    current,
                    candidate,
                    best_sum_score,
                    best_evalue,
                ));
            } else {
                head = None;
                index += 1;
            }
        }
    }
    work.sort_by(|a, b| compare_score(&a.value.hsp, &b.value.hsp));
    let hsps: Vec<_> = work.into_iter().map(|item| item.value).collect();
    let best_evalue = hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .reduce(f64::min)
        .unwrap_or(f64::MAX);
    Ok(LinkedHspList { hsps, best_evalue })
}
