//! Context-local BLASTX preliminary HSPs before linking/traceback.
use super::{
    args::ResolvedOptions, parameters::ContextParameters, query_setup::PreparedQueryBatch,
    seed::InitHsp,
};
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    BlastpGappedAlignmentMode, GapAlignScratch,
};
use crate::config::ScoringMatrix;
use anyhow::Result;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-103
// ```c++
// typedef struct BlastSeg {
//    Int2 frame;  /**< Translation frame */
//    Int4 offset; /**< Start of hsp */
//    Int4 end;    /**< End of hsp */
//    Int4 gapped_start;/**< Where the gapped extension started. */
// } BlastSeg;
//
// /** In PHI BLAST: information about pattern match in a given HSP. */
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-143
// ```c++
// typedef struct BlastHSP {
//    Int4 score;           /**< This HSP's raw score */
//    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
//    double bit_score;     /**< Bit score, calculated from score */
//    double evalue;        /**< This HSP's e-value */
//    BlastSeg query;       /**< Query sequence info. */
//    BlastSeg subject;     /**< Subject sequence info. */
//    Int4     context;     /**< Context number of query */
//    GapEditScript* gap_info;/**< ALL gapped alignment is here */
//    Int4 num;             /**< How many HSP's are linked together for sum
//                               statistics evaluation? If unset (0), this HSP is
//                               not part of a linked set, i.e. value 0 is treated
//                               the same way as 1. */
//    Int2		comp_adjustment_method;  /**< which mode of composition
//                                               adjustment was used; relevant
//                                               only for blastp and tblastn */
//    SPHIHspInfo* pat_info; /**< In PHI BLAST, information about this pattern
//                                  match. */
// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreliminaryHsp {
    pub context: usize,
    pub frame: i8,
    pub score: i32,
    pub q_start: i32,
    pub q_end: i32,
    pub q_gapped_start: i32,
    pub s_start: i32,
    pub s_end: i32,
    pub s_gapped_start: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
// ```c++
// s_GetQueryStrandOffset(const BlastQueryInfo *query_info,
//                        Int4 context)
// {
//     Int4 c = context;
//
//     while (c) {
//         Int4 frame = query_info->contexts[c].frame;
//         if (frame == 0 || SIGN(frame) !=
//             SIGN(query_info->contexts[c-1].frame)) {
//             break;
//         }
//         c--;
//     }
//
//     return query_info->contexts[c].query_offset;
// ```
fn tree_hsp(h: &PreliminaryHsp, batch: &PreparedQueryBatch) -> TreeHsp {
    let strand_context = h.context / 6 * 6 + if h.context % 6 < 3 { 0 } else { 3 };
    TreeHsp {
        query_offset: h.q_start,
        query_end: h.q_end,
        subject_offset: h.s_start,
        subject_end: h.s_end,
        score: h.score,
        query_frame: h.frame as i32,
        query_length: batch.contexts[h.context].length as i32,
        query_context_offset: batch.contexts[strand_context].offset as i32,
        subject_frame_sign: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1635-1675
// ```c++
//       if (left_son == last)
//          large_son = left_son;
//       else
//          large_son = (*compar)(left_son, left_son+width) >= 0 ?
//             left_son : left_son+width;
//       if ((*compar)(base, large_son) < 0) {
//          for (i=0; i<width; ++i) {
//             ch = base[i];
//             base[i] = large_son[i];
//             large_son[i] = ch;
//          }
//          base = large_son;
//          left_son = base0 + 2*(base-base0) + width;
//       } else
//          break;
//    }
// }
//
// /** Creates a heap of elements based on a comparison function.
//  * @param b An array [in] [out]
//  * @param nel Number of elements in b [in]
//  * @param width The size of each element [in]
//  * @param compar Callback to compare two heap elements [in]
//  */
// static void
// s_CreateHeap (void* b, size_t nel, size_t width,
//    int (*compar )(const void*, const void* ))
// {
//    char*    base = (char*)b;
//    size_t i;
//    char*    base0 = (char*)base,* lim,* basef;
//
//    if (nel < 2)
//       return;
//
//    lim = &base[((nel-2)/2)*width];
//    basef = &base[(nel-1)*width];
//    i = nel/2;
//    for (base = &base0[(i - 1)*width]; i > 0; base = base - width) {
//       s_Heapify(base0, base, lim, basef, width, compar);
//       i--;
// ```
fn heapify(hits: &mut [PreliminaryHsp], mut base: usize) {
    while hits.len() >= 2 && base <= (hits.len() - 2) / 2 {
        let left = 2 * base + 1;
        let worst = if left == hits.len() - 1 || compare_score(&hits[left], &hits[left + 1]).is_ge()
        {
            left
        } else {
            left + 1
        };
        if compare_score(&hits[base], &hits[worst]).is_lt() {
            hits.swap(base, worst);
            base = worst;
        } else {
            break;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1687-1706
// ```c++
// s_BlastHSPListInsertHSPInHeap(BlastHSPList* hsp_list,
//                              BlastHSP** hsp)
// {
//     BlastHSP** hsp_array = hsp_list->hsp_array;
//     if (ScoreCompareHSPs(hsp, &hsp_array[0]) > 0)
//     {
//          Blast_HSPFree(*hsp);
//          return;
//     }
//     else
//          Blast_HSPFree(hsp_array[0]);
//
//     hsp_array[0] = *hsp;
//     if (hsp_list->hspcnt >= 2) {
//         s_Heapify((char*)hsp_array, (char*)hsp_array,
//                 (char*)&hsp_array[hsp_list->hspcnt/2 - 1],
//                  (char*)&hsp_array[hsp_list->hspcnt-1],
//                  sizeof(BlastHSP*), ScoreCompareHSPs);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1751-1809
// ```c++
// /* Comments in blast_hits.h
//  */
// Int2
// Blast_HSPListSaveHSP(BlastHSPList* hsp_list, BlastHSP* new_hsp)
// {
//    BlastHSP** hsp_array;
//    Int4 hspcnt;
//    Int4 hsp_allocated; /* how many hsps are in the array. */
//    Int2 status = 0;
//
//    hspcnt = hsp_list->hspcnt;
//    hsp_allocated = hsp_list->allocated;
//    hsp_array = hsp_list->hsp_array;
//
//
//    /* Check if list is already full, then reallocate. */
//    if (hspcnt >= hsp_allocated && hsp_list->do_not_reallocate == FALSE)
//    {
//       Int4 new_allocated = MIN(2*hsp_list->allocated, hsp_list->hsp_max);
//       if (new_allocated > hsp_list->allocated) {
//          hsp_array = (BlastHSP**)
//             realloc(hsp_array, new_allocated*sizeof(BlastHSP*));
//          if (hsp_array == NULL)
//          {
//             hsp_list->do_not_reallocate = TRUE;
//             hsp_array = hsp_list->hsp_array;
//             /** Return a non-zero status, because restriction on number
//                 of HSPs here is a result of memory allocation failure. */
//             status = -1;
//          } else {
//             hsp_list->hsp_array = hsp_array;
//             hsp_list->allocated = new_allocated;
//             hsp_allocated = new_allocated;
//          }
//       } else {
//          hsp_list->do_not_reallocate = TRUE;
//       }
//       /* If it is the first time when the HSP array is filled to capacity,
//          create a heap now. */
//       if (hsp_list->do_not_reallocate) {
//           s_CreateHeap(hsp_array, hspcnt, sizeof(BlastHSP*), ScoreCompareHSPs);
//       }
//    }
//
//    /* If there is space in the allocated HSP array, simply save the new HSP.
//       Othewise, if the new HSP has lower score than the worst HSP in the heap,
//       then delete it, else insert it in the heap. */
//    if (hspcnt < hsp_allocated)
//    {
//       hsp_array[hsp_list->hspcnt] = new_hsp;
//       (hsp_list->hspcnt)++;
//       return status;
//    } else {
//        /* Insert the new HSP in heap. */
//        s_BlastHSPListInsertHSPInHeap(hsp_list, &new_hsp);
//    }
//
//    return status;
// }
// ```
fn save_hsp(hits: &mut Vec<PreliminaryHsp>, h: PreliminaryHsp, max: usize, heap: &mut bool) {
    if hits.len() < max {
        hits.push(h);
        return;
    }
    if !*heap {
        for i in (0..hits.len() / 2).rev() {
            heapify(hits, i);
        }
        *heap = true;
    }
    if compare_score(&h, &hits[0]).is_gt() {
        return;
    }
    hits[0] = h;
    heapify(hits, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3860-3869
// ```c++
//       tmp_init_hsp = init_hsp_array[index];
//       if (tmp_init_hsp.ungapped_data) {
//           tmp_ungapped_data = *(init_hsp_array[index].ungapped_data);
//           tmp_init_hsp.ungapped_data = &tmp_ungapped_data;
//       }
//       init_hsp = &tmp_init_hsp;
//
//       s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
//                                         &query_tmp, &context);
//
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3953
// ```c++
//       } else {
//          q_start = init_hsp->ungapped_data->q_start;
//          q_end = q_start + init_hsp->ungapped_data->length;
//          s_start = init_hsp->ungapped_data->s_start;
//          s_end = s_start + init_hsp->ungapped_data->length;
//          score = init_hsp->ungapped_data->score;
//       }
//
//       tmp_hsp.score = score;
//       tmp_hsp.context = context;
//       tmp_hsp.query.offset = q_start;
//       tmp_hsp.query.end = q_end;
//       tmp_hsp.query.frame = query_info->contexts[context].frame;
//       tmp_hsp.subject.offset = s_start;
//       tmp_hsp.subject.end = s_end;
//       tmp_hsp.subject.frame = subject->frame;
//
//       /* use priate interval tree when recomputing alignments */
//       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
//                                         hit_options->min_diag_separation))
//       {
//          BlastHSP* new_hsp;
//          Int4 cutoff, restricted_cutoff = 0;
//
//          if (is_rpsblast)
//             cutoff = hit_params->cutoffs[rps_cutoff_index].cutoff_score;
//          else
//             cutoff = hit_params->cutoffs[context].cutoff_score;
//
//          if (restricted_alignment)
//              restricted_cutoff = (Int4)(kRestrictedMult * cutoff);
//
//          if (gapped_stats)
//             ++gapped_stats->extensions;
//
//          if(is_prot && !score_params->options->is_ooframe) {
//             max_offset =
//                BlastGetStartForGappedAlignment(query_tmp.sequence,
//                   subject->sequence, gap_align->sbp,
//                   init_hsp->ungapped_data->q_start,
//                   init_hsp->ungapped_data->length,
//                   init_hsp->ungapped_data->s_start,
//                   init_hsp->ungapped_data->length);
//             init_hsp->offsets.qs_offsets.s_off +=
//                 max_offset - init_hsp->offsets.qs_offsets.q_off;
//             init_hsp->offsets.qs_offsets.q_off = max_offset;
//          }
//
//          if (is_prot) {
//             status = s_BlastProtGappedAlignment(program_number, &query_tmp,
//                                                 subject, gap_align,
//                                                 score_params, init_hsp,
//                                                 restricted_alignment,
//                                                 fence_hit);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4087
// ```c++
//          if (gap_align->score >= cutoff) {
//              Int2 query_frame = 0;
//              /* For mixed-frame search, the query frame is determined
//                 from the offset, not only from context. */
//              if (score_params->options->is_ooframe &&
//                  program_number == eBlastTypeBlastx) {
//                  query_frame = gap_align->query_start % CODON_LENGTH + 1;
//                  if ((context % NUM_FRAMES) >= CODON_LENGTH)
//                      query_frame = -query_frame;
//              } else {
//                  query_frame = query_info->contexts[context].frame;
//              }
//
//              status = Blast_HSPInit(gap_align->query_start,
//                            gap_align->query_stop, gap_align->subject_start,
//                            gap_align->subject_stop,
//                            init_hsp->offsets.qs_offsets.q_off,
//                            init_hsp->offsets.qs_offsets.s_off, context,
//                            query_frame, subject->frame, gap_align->score,
//                            &(gap_align->edit_script), &new_hsp);
//              if (status)
// 	     {
//    		sfree(found_high_score);
//    		tree = Blast_IntervalTreeFree(tree);
//                 return status;
// 	     }
//              status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
//              if (status)
//                  break;
//              status = BlastIntervalTreeAddHSP(new_hsp, tree, query_info,
// ```
pub fn gapped(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    initial: &[InitHsp],
    xdrop: i32,
) -> Result<Vec<PreliminaryHsp>> {
    gapped_observed(
        batch,
        parameters,
        options,
        subject,
        initial,
        xdrop,
        &mut |_, _| {},
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3909-3921
// ```c++
//       tmp_hsp.context = context;
//       tmp_hsp.query.offset = q_start;
//       tmp_hsp.query.end = q_end;
//       tmp_hsp.query.frame = query_info->contexts[context].frame;
//       tmp_hsp.subject.offset = s_start;
//       tmp_hsp.subject.end = s_end;
//       tmp_hsp.subject.frame = subject->frame;
//
//       /* use priate interval tree when recomputing alignments */
//       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
//                                         hit_options->min_diag_separation))
//       {
//          BlastHSP* new_hsp;
// ```
pub(crate) fn gapped_observed(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    // EXPERIMENT (LOSAT_X_BXLEAN): most (chunk, subject) pairs have no initial
    // HSP; the loop below then does nothing, so the tree and the DP scratch
    // need not be built.
    if initial.is_empty() && super::runtime::x_bx_lean() {
        return Ok(Vec::new());
    }
    let last = batch.contexts.last().expect("BLASTX contexts");
    let mut tree = BlastIntervalTree::new(
        0,
        (last.offset + last.length + 1) as i32,
        0,
        subject.len() as i32,
    );
    let mut scratch = GapAlignScratch::new();
    let mut hits = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:213-227
    // ```c++
    // Int4 BlastHspNumMax(Boolean gapped_calculation, const BlastHitSavingOptions* options)
    // {
    //    Int4 retval=0;
    //
    //    /* per-subject HSP limits do not apply to gapped searches; JIRA SB-616 */
    //    if (options->hsp_num_max <= 0)
    //    {
    //       retval = INT4_MAX;
    //    }
    //    else
    //    {
    //       retval = options->hsp_num_max;
    //    }
    //
    //    return retval;
    // ```
    // hsp_num_max defaults to zero; max_hsps_per_subject is a later collector filter.
    let max = i32::MAX as usize;
    let mut heap = false;
    for init in initial {
        let context = batch
            .contexts
            .partition_point(|c| c.offset as i32 <= init.q_seed)
            .saturating_sub(1);
        let c = &batch.contexts[context];
        let mut h = PreliminaryHsp {
            context,
            frame: c.frame,
            score: init.score,
            q_start: init.q_start - c.offset as i32,
            q_end: init.q_start - c.offset as i32 + init.length,
            q_gapped_start: init.q_seed - c.offset as i32,
            s_start: init.s_start,
            s_end: init.s_start + init.length,
            s_gapped_start: init.s_seed,
        };
        let candidate = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3909-3921
        // ```c++
        //       tmp_hsp.context = context;
        //       tmp_hsp.query.offset = q_start;
        //       tmp_hsp.query.end = q_end;
        //       tmp_hsp.query.frame = query_info->contexts[context].frame;
        //       tmp_hsp.subject.offset = s_start;
        //       tmp_hsp.subject.end = s_end;
        //       tmp_hsp.subject.frame = subject->frame;
        //
        //       /* use priate interval tree when recomputing alignments */
        //       if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
        //                                         hit_options->min_diag_separation))
        //       {
        //          BlastHSP* new_hsp;
        // ```
        let contained = tree.contains_hsp(&candidate, candidate.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
        let start = blastp_get_start_for_gapped_alignment(
            query,
            &subject[..subject.len() - 1],
            h.q_start as usize,
            init.length as usize,
            h.s_start as usize,
            init.length as usize,
            ScoringMatrix::Blosum62,
        );
        h.s_gapped_start += start as i32 - h.q_gapped_start;
        h.q_gapped_start = start as i32;
        let a = blastp_score_only_gapped_alignment_with_scratch(
            query,
            &subject[..subject.len() - 1],
            start,
            h.s_gapped_start as usize,
            ScoringMatrix::Blosum62,
            options.gap_open,
            options.gap_extend,
            xdrop,
            BlastpGappedAlignmentMode::Exact,
            &mut scratch,
        );
        if a.score < parameters[context].hit_cutoff {
            continue;
        }
        h.score = a.score;
        h.q_start = a.query_start;
        h.q_end = a.query_stop;
        h.s_start = a.subject_start;
        h.s_end = a.subject_stop;
        let node = tree_hsp(&h, batch);
        save_hsp(&mut hits, h, max, &mut heap);
        tree.add_hsp(
            node,
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
    }
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4745-4777
// ```c++
//       Int4 context = 0;
//       init_hsp = &init_hitlist->init_hsp_array[index];
//       if (!init_hsp->ungapped_data)
//          continue;
//
//       if (!hsp_list) {
//          hsp_list = Blast_HSPListNew(kHspNumMax);
//          *hsp_list_ptr = hsp_list;
//       }
//       /* Adjust the initial HSP's coordinates in case of concatenated
//          multiple queries/strands/frames */
//       context = s_GetUngappedHSPContext(query_info, init_hsp);
//       s_AdjustInitialHSPOffsets(init_hsp,
//                                 query_info->contexts[context].query_offset);
//       ungapped_data = init_hsp->ungapped_data;
//       Blast_HSPInit(ungapped_data->q_start,
//                     ungapped_data->length+ungapped_data->q_start,
//                     ungapped_data->s_start,
//                     ungapped_data->length+ungapped_data->s_start,
//                     init_hsp->offsets.qs_offsets.q_off,
//                     init_hsp->offsets.qs_offsets.s_off,
//                     context, query_info->contexts[context].frame,
//                     subject->frame, ungapped_data->score, NULL, &new_hsp);
//       Blast_HSPListSaveHSP(hsp_list, new_hsp);
//    }
//
//    /* Sort the HSP array by score */
//    Blast_HSPListSortByScore(hsp_list);
//
//    return 0;
// }
//
// ```
pub fn ungapped(batch: &PreparedQueryBatch, initial: &[InitHsp]) -> Vec<PreliminaryHsp> {
    let mut hits = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:213-227
    // ```c++
    // Int4 BlastHspNumMax(Boolean gapped_calculation, const BlastHitSavingOptions* options)
    // {
    //    Int4 retval=0;
    //
    //    /* per-subject HSP limits do not apply to gapped searches; JIRA SB-616 */
    //    if (options->hsp_num_max <= 0)
    //    {
    //       retval = INT4_MAX;
    //    }
    //    else
    //    {
    //       retval = options->hsp_num_max;
    //    }
    //
    //    return retval;
    // ```
    // hsp_num_max defaults to zero; max_hsps_per_subject is a later collector filter.
    let max = i32::MAX as usize;
    let mut heap = false;
    for init in initial {
        let context = batch
            .contexts
            .partition_point(|c| c.offset as i32 <= init.q_seed)
            .saturating_sub(1);
        let c = &batch.contexts[context];
        save_hsp(
            &mut hits,
            PreliminaryHsp {
                context,
                frame: c.frame,
                score: init.score,
                q_start: init.q_start - c.offset as i32,
                q_end: init.q_start - c.offset as i32 + init.length,
                q_gapped_start: init.q_seed - c.offset as i32,
                s_start: init.s_start,
                s_end: init.s_start + init.length,
                s_gapped_start: init.s_seed,
            },
            max,
            &mut heap,
        );
    }
    sort_by_score(&mut hits);
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1358
// ```c++
// ScoreCompareHSPs(const void* h1, const void* h2)
// {
//    BlastHSP* hsp1,* hsp2;   /* the HSPs to be compared */
//    int result = 0;      /* the result of the comparison */
//
//    hsp1 = *((BlastHSP**) h1);
//    hsp2 = *((BlastHSP**) h2);
//
//    /* Null HSPs are "greater" than any non-null ones, so they go to the end
//       of a sorted list. */
//    if (!hsp1 && !hsp2)
//        return 0;
//    else if (!hsp1)
//        return 1;
//    else if (!hsp2)
//        return -1;
//
//    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
//        /* if all other test can't distinguish the HSPs, then the final
//           test is the result */
//        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
//    }
//    return result;
// }
//
// Boolean Blast_HSPListIsSortedByScore(const BlastHSPList* hsp_list)
// ```
pub fn compare_score(a: &PreliminaryHsp, b: &PreliminaryHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
// ```c++
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
pub fn sort_by_score(hits: &mut [PreliminaryHsp]) {
    hits.sort_by(compare_score);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2295-2321
// ```c++
//
//    if (h1->subject.offset < h2->subject.offset)
//       return -1;
//    if (h1->subject.offset > h2->subject.offset)
//       return 1;
//
//    /* tie breakers: sort by decreasing score, then
//       by increasing size of query range, then by
//       increasing subject range. */
//
//    if (h1->score < h2->score)
//       return 1;
//    if (h1->score > h2->score)
//       return -1;
//
//    if (h1->query.end < h2->query.end)
//       return 1;
//    if (h1->query.end > h2->query.end)
//       return -1;
//
//    if (h1->subject.end < h2->subject.end)
//       return 1;
//    if (h1->subject.end > h2->subject.end)
//       return -1;
//
//    return 0;
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2362-2387
// ```c++
//       return -1;
//    if (h1->subject.end > h2->subject.end)
//       return 1;
//
//    /* tie breakers: sort by decreasing score, then
//       by increasing size of query range, then by
//       increasing size of subject range. The shortest range
//       means the *largest* sequence offset must come
//       first */
//    if (h1->score < h2->score)
//       return 1;
//    if (h1->score > h2->score)
//       return -1;
//
//    if (h1->query.offset < h2->query.offset)
//       return 1;
//    if (h1->query.offset > h2->query.offset)
//       return -1;
//
//    if (h1->subject.offset < h2->subject.offset)
//       return 1;
//    if (h1->subject.offset > h2->subject.offset)
//       return -1;
//
//    return 0;
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2505
// ```c++
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
//    i = 0;
//    while (i < hsp_count) {
//       j = 1;
//       while (i+j < hsp_count &&
//              hsp_array[i] && hsp_array[i+j] &&
//              hsp_array[i]->context == hsp_array[i+j]->context &&
//              hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
//              hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
//              hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
//          hsp_count--;
//          hsp = hsp_array[i+j];
//          if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
//              s_CutOffGapEditScript(hsp, hsp_array[i]->query.end,
//                                         hsp_array[i]->subject.end, TRUE);
//          } else {
//              hsp = Blast_HSPFree(hsp);
//          }
//          for (k=i+j; k<hsp_count; k++) {
//              hsp_array[k] = hsp_array[k+1];
//          }
//          hsp_array[hsp_count] = hsp;
//       }
//       i += j;
//    }
//
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
//    i = 0;
// ```
pub fn purge_endpoints(hits: &mut Vec<PreliminaryHsp>) {
    hits.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then(a.q_start.cmp(&b.q_start))
            .then(a.s_start.cmp(&b.s_start))
            .then(b.score.cmp(&a.score))
            .then(b.q_end.cmp(&a.q_end))
            .then(b.s_end.cmp(&a.s_end))
    });
    hits.dedup_by(|b, a| {
        a.context == b.context && a.q_start == b.q_start && a.s_start == b.s_start
    });
    hits.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then(a.q_end.cmp(&b.q_end))
            .then(a.s_end.cmp(&b.s_end))
            .then(b.score.cmp(&a.score))
            .then(b.q_start.cmp(&a.q_start))
            .then(b.s_start.cmp(&a.s_start))
    });
    hits.dedup_by(|b, a| a.context == b.context && a.q_end == b.q_end && a.s_end == b.s_end);
    sort_by_score(hits);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
// ```c++
// s_GetQueryStrandOffset(const BlastQueryInfo *query_info,
//                        Int4 context)
// {
//     Int4 c = context;
//
//     while (c) {
//         Int4 frame = query_info->contexts[c].frame;
//         if (frame == 0 || SIGN(frame) !=
//             SIGN(query_info->contexts[c-1].frame)) {
//             break;
//         }
//         c--;
//     }
//
//     return query_info->contexts[c].query_offset;
// ```
#[cfg(test)]
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
    // ```c++
    // s_GetQueryStrandOffset(const BlastQueryInfo *query_info,
    //                        Int4 context)
    // {
    //     Int4 c = context;
    //
    //     while (c) {
    //         Int4 frame = query_info->contexts[c].frame;
    //         if (frame == 0 || SIGN(frame) !=
    //             SIGN(query_info->contexts[c-1].frame)) {
    //             break;
    //         }
    //         c--;
    //     }
    //
    //     return query_info->contexts[c].query_offset;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819-834
    // ```c++
    //     if (in_q_start != tree_q_start)
    //         return FALSE;
    //
    //     if (in_hsp->score <= tree_hsp->score &&
    //         SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
    //         CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end,
    //                               in_hsp->query.offset,
    //                               tree_hsp->subject.offset, tree_hsp->subject.end,
    //                               in_hsp->subject.offset) &&
    //         CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end,
    //                              in_hsp->query.end,
    //                              tree_hsp->subject.offset, tree_hsp->subject.end,
    //                              in_hsp->subject.end)) {
    //
    //         if (min_diag_separation == 0)
    //             return TRUE;
    // ```
    #[test]
    fn containment_groups_three_frames_by_query_strand() {
        let (_, batch, _) = super::super::stage_c_tests::fixture("S06_seed.fna");
        let h = PreliminaryHsp {
            context: 0,
            frame: 1,
            score: 100,
            q_start: 0,
            q_end: 70,
            q_gapped_start: 2,
            s_start: 0,
            s_end: 70,
            s_gapped_start: 2,
        };
        let mut tree = BlastIntervalTree::new(0, 600, 0, 100);
        let node = tree_hsp(&h, &batch);
        tree.add_hsp(
            node,
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
        for c in 0..6 {
            let mut candidate = h.clone();
            candidate.context = c;
            candidate.frame = batch.contexts[c].frame;
            candidate.score = 90;
            candidate.q_start = 5;
            candidate.q_end = 60;
            candidate.s_start = 5;
            candidate.s_end = 60;
            let node = tree_hsp(&candidate, &batch);
            assert_eq!(
                tree.contains_hsp(&node, node.query_context_offset, 0),
                c < 3
            );
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1358
    // ```c++
    // ScoreCompareHSPs(const void* h1, const void* h2)
    // {
    //    BlastHSP* hsp1,* hsp2;   /* the HSPs to be compared */
    //    int result = 0;      /* the result of the comparison */
    //
    //    hsp1 = *((BlastHSP**) h1);
    //    hsp2 = *((BlastHSP**) h2);
    //
    //    /* Null HSPs are "greater" than any non-null ones, so they go to the end
    //       of a sorted list. */
    //    if (!hsp1 && !hsp2)
    //        return 0;
    //    else if (!hsp1)
    //        return 1;
    //    else if (!hsp2)
    //        return -1;
    //
    //    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
    //        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
    //        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
    //        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
    //        /* if all other test can't distinguish the HSPs, then the final
    //           test is the result */
    //        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    //    }
    //    return result;
    // }
    //
    // Boolean Blast_HSPListIsSortedByScore(const BlastHSPList* hsp_list)
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1687-1706
    // ```c++
    // s_BlastHSPListInsertHSPInHeap(BlastHSPList* hsp_list,
    //                              BlastHSP** hsp)
    // {
    //     BlastHSP** hsp_array = hsp_list->hsp_array;
    //     if (ScoreCompareHSPs(hsp, &hsp_array[0]) > 0)
    //     {
    //          Blast_HSPFree(*hsp);
    //          return;
    //     }
    //     else
    //          Blast_HSPFree(hsp_array[0]);
    //
    //     hsp_array[0] = *hsp;
    //     if (hsp_list->hspcnt >= 2) {
    //         s_Heapify((char*)hsp_array, (char*)hsp_array,
    //                 (char*)&hsp_array[hsp_list->hspcnt/2 - 1],
    //                  (char*)&hsp_array[hsp_list->hspcnt-1],
    //                  sizeof(BlastHSP*), ScoreCompareHSPs);
    //     }
    // }
    // ```
    #[test]
    fn capped_save_heap_ties_preserve_ncbi_comparator() {
        let h = |score, s_start| PreliminaryHsp {
            context: 0,
            frame: 1,
            score,
            q_start: 0,
            q_end: 5,
            q_gapped_start: 1,
            s_start,
            s_end: s_start + 5,
            s_gapped_start: s_start + 1,
        };
        let mut hits = Vec::new();
        let mut heap = false;
        for x in [h(30, 5), h(20, 0), h(30, 1), h(10, 0)] {
            save_hsp(&mut hits, x, 2, &mut heap);
        }
        sort_by_score(&mut hits);
        assert_eq!(hits, vec![h(30, 1), h(30, 5)]);
    }
}
