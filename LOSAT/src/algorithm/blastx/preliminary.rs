//! Context-local BLASTX preliminary HSPs before linking/traceback.
use super::{
    args::ResolvedOptions, parameters::ContextParameters, query_setup::PreparedQueryBatch,
    seed::InitHsp,
};
use crate::algorithm::blastn::interval_tree::{
    BlastIntervalTree, IndexMethod, TreeHsp, X_MAX_CELLS,
};
use crate::algorithm::blastp::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    gapped_start_int4, gapped_start_window, BlastpGappedAlignmentMode, GapAlignScratch,
    GAPPED_START_WINDOW,
};
use crate::config::ScoringMatrix;
use anyhow::Result;

/// The bytes that NCBI's buffer of the searched subject holds after the sentinel byte that follows
/// its residues, and the name of that buffer (`subjects::Subjects::buffer_after`), asked only where
/// a one-hit gapped start's window passes the subject's sentinel.
pub(crate) type SubjectAfter<'a> = &'a dyn Fn() -> Result<(usize, &'static str)>;

/// LOSAT's refusal of a search whose one-hit gapped start reads past NCBI's buffers
/// (`refuse_window_past_buffers`), not an error of NCBI's engine (`native::database_search_error`,
/// `local::run`).
#[derive(Debug)]
pub(crate) struct OneHitWindowRefusal(String);

impl std::fmt::Display for OneHitWindowRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for OneHitWindowRefusal {}

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

// No NCBI counterpart: grid cells per HSP of a subject's tree (`x_tree_cells`); memory of the
// grid index only.
const X_TREE_CELLS_PER_HSP: usize = 64;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
// ```c
//       tree = Blast_IntervalTreeInit(0, query->length+1,
//                                     0, subject->length+1);
// ```
/// No NCBI counterpart: the grid limit (`reset_with_bounds_cells`) of the tree NCBI makes for each
/// subject. The tree holds at most one HSP per initial HSP, so its grid index gets
/// `X_TREE_CELLS_PER_HSP` cells per initial HSP, at most the grid's own limit (`X_MAX_CELLS`, four
/// bytes per cell, zeroed on the first HSP): a subject with a few initial HSPs no longer takes a
/// grid of up to 256 KiB for them. The count only sizes the grid: a smaller grid lists more leaves
/// per cell and answers the same, whatever the tree holds.
fn x_tree_cells(initial: usize) -> usize {
    initial
        .saturating_mul(X_TREE_CELLS_PER_HSP)
        .clamp(X_TREE_CELLS_PER_HSP, X_MAX_CELLS)
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
    subject_after: SubjectAfter<'_>,
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    gapped_observed_with(
        super::seed::Blosum62Bx,
        batch,
        parameters,
        options,
        subject,
        subject_after,
        initial,
        xdrop,
        observe,
    )
}

/// [`gapped_observed`] with another matrix than BLOSUM62 (item SC6 of line NP). Not inlined, so that
/// BLOSUM62's keeps its code.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn gapped_observed_of<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    subject_after: SubjectAfter<'_>,
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    gapped_observed_with(
        scoring,
        batch,
        parameters,
        options,
        subject,
        subject_after,
        initial,
        xdrop,
        observe,
    )
}

/// The body of [`gapped_observed`], with the matrix of `scoring`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn gapped_observed_with<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    subject_after: SubjectAfter<'_>,
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3707-3708
    // ```c
    //    if (init_hitlist->total == 0)
    //       return 0;
    // ```
    // NCBI returns at once when the initial hit list is empty (and the engine skips the gapped
    // stage, blast_engine.c:500). This port returns the empty list before it builds the interval
    // tree and the DP scratch; the loop below would have done nothing. Most (chunk, subject)
    // pairs have no initial HSP.
    if initial.is_empty() {
        return Ok(Vec::new());
    }
    let last = batch.contexts.last().expect("BLASTX contexts");
    let (q_max, s_max) = ((last.offset + last.length + 1) as i32, subject.len() as i32);
    let mut tree = BlastIntervalTree::new(0, q_max, 0, s_max);
    tree.reset_with_bounds_cells(0, q_max, 0, s_max, x_tree_cells(initial.len()));
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3935-3942
        // ```c
        //          if(is_prot && !score_params->options->is_ooframe) {
        //             max_offset =
        //                BlastGetStartForGappedAlignment(query_tmp.sequence,
        //                   subject->sequence, gap_align->sbp,
        //                   init_hsp->ungapped_data->q_start,
        //                   init_hsp->ungapped_data->length,
        //                   init_hsp->ungapped_data->s_start,
        //                   init_hsp->ungapped_data->length);
        // ```
        // A negative length comes only from the one-hit extension (`-window_size 0`,
        // aa_ungapped.c:1054,1083), which NCBI passes as `Uint4` (`gapped_start_negative`).
        let start = if init.length < 0 {
            gapped_start_negative(
                batch,
                c.offset,
                query,
                &subject[..subject.len() - 1],
                h.q_start as usize,
                h.s_start as usize,
                scoring.matrix(),
                subject_after,
            )?
        } else {
            blastp_get_start_for_gapped_alignment(
                query,
                &subject[..subject.len() - 1],
                h.q_start as usize,
                init.length as usize,
                h.s_start as usize,
                init.length as usize,
                scoring.matrix(),
            )
        };
        h.s_gapped_start += start as i32 - h.q_gapped_start;
        h.q_gapped_start = start as i32;
        let a = blastp_score_only_gapped_alignment_with_scratch(
            query,
            &subject[..subject.len() - 1],
            start,
            h.s_gapped_start as usize,
            scoring.matrix(),
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083 and
// c++/src/algo/blast/core/blast_gapalign.c:3394,3401,3405-3413,3419-3420
// ```c
//     init_hit_width = q_right_off - q_left_off + 1;
//     *hsp_len = left_disp + right_disp + init_hit_width;
// ...
//    const BlastScoreBlk* sbp, Uint4 q_start, Uint4 q_length,
//     if (q_length <= HSP_MAX_WINDOW) {
//     hsp_end = q_start + HSP_MAX_WINDOW;
//     query_var = query + q_start;
//     subject_var = subject + s_start;
// ...
//     hsp_end = q_start + MIN(q_length, s_length);
//     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
// ```
/// The gapped start of a one-hit HSP of negative length at `q_start` of the context at `offset`
/// of the query batch (`query` its letters) and at `s_start` of `subject` (its letters). NCBI
/// reads the length as `Uint4`, so it scores only the first window of `HSP_MAX_WINDOW` letters
/// (`gapped_start_int4`). A window that passes the end of its context or of the subject holds the
/// NULLB after it, so its sum is below 0 and the start is `q_start` whatever bytes follow; the
/// window reads those bytes from NCBI's buffers, and past them it reads memory that is not NCBI's
/// (`refuse_window_past_buffers`). Out of line, so that the gapped stage of the other HSPs keeps
/// its code.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn gapped_start_negative(
    batch: &PreparedQueryBatch,
    offset: usize,
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
    subject_after: SubjectAfter<'_>,
) -> Result<usize> {
    refuse_window_past_buffers(
        batch,
        offset + q_start,
        subject.len(),
        s_start,
        subject_after,
    )?;
    Ok(gapped_start_int4(
        &gapped_start_window(query, q_start),
        &gapped_start_window(subject, s_start),
        q_start,
        matrix,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,885-893 and
// c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1802
// ```c++
//     int buflen = QueryInfo_GetSeqBufLen(qinfo);
//     TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
// ...
//     buflen = CalculateSeqBufferLength(sv.size(), eBlastEncodingProtein);
//     _ASSERT(buflen != 0);
//     buf = buf_var = (Uint1*) malloc(sizeof(Uint1)*buflen);
// ...
//     *buf_var++ = GetSentinelByte(eBlastEncodingProtein);
// ...
//         *buffer = m_Seq->GetFileDataPtr(start_offset-1) + 1;
// ```
/// The refusal (PD-LOSAT-NCBI-DEFECTS, rule 3) of the window of `HSP_MAX_WINDOW` letters from
/// `q_abs` of the query batch (`query->sequence`) and from `s_start` of a subject of `s_length`
/// residues, where it reads past NCBI's buffers and NCBI's result is not defined: past the query
/// batch's buffer (`calloc(buflen + 1)`, `PreparedQueryBatch::batch_buffer_length`), or past the
/// subject's (`-subject`: `malloc(len + 2)`, its residues between two sentinel bytes; `-db`: the
/// volume's sequence file, mapped whole, `subject_after`). Inside them the search gives NCBI's
/// result. This is the search's one check for it: reading the bytes past NCBI's buffers as
/// sentinels too, as BLASTP's exception 3 does (PD-LOSAT-NCBI-DEFECTS 1.3), would remove it.
fn refuse_window_past_buffers(
    batch: &PreparedQueryBatch,
    q_abs: usize,
    s_length: usize,
    s_start: usize,
    subject_after: SubjectAfter<'_>,
) -> Result<()> {
    let past = if q_abs + GAPPED_START_WINDOW > batch.batch_buffer_length() {
        "query batch"
    } else if s_start + GAPPED_START_WINDOW > s_length + 1 {
        let (after, buffer) = subject_after()?;
        if s_start + GAPPED_START_WINDOW <= s_length + 1 + after {
            return Ok(());
        }
        buffer
    } else {
        return Ok(());
    };
    Err(OneHitWindowRefusal(format!(
        "BLASTX -window_size 0: the gapped start of a one-hit alignment of negative length reads past the end of the {past} (NCBI BLAST+ reads past its buffer there), which LOSAT does not support"
    ))
    .into())
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
    // No NCBI counterpart: the grid limit of a subject's tree (`x_tree_cells`).
    #[test]
    fn tree_grid_follows_the_initial_hsps() {
        assert_eq!(x_tree_cells(0), X_TREE_CELLS_PER_HSP);
        assert_eq!(x_tree_cells(1), X_TREE_CELLS_PER_HSP);
        assert_eq!(x_tree_cells(10), 10 * X_TREE_CELLS_PER_HSP);
        assert_eq!(x_tree_cells(1 << 20), X_MAX_CELLS);
        assert_eq!(x_tree_cells(usize::MAX), X_MAX_CELLS);
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:245-250,
    // c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,1236-1239 (quoted at
    // `PreparedQueryBatch::batch_buffer_length`, `Subjects::buffer_after` and
    // `refuse_window_past_buffers`)
    /// A one-hit gapped start of negative length is refused exactly where its window of 11 letters
    /// reads past NCBI's buffers: past calloc's byte after the query batch's last NULLB, past the
    /// sentinel after a `-subject` record, past what a database sequence's volume holds after its
    /// sentinel; the subject's bytes after its sentinel are asked for only when the window passes
    /// it.
    #[test]
    fn a_one_hit_window_is_refused_only_past_ncbis_buffers() {
        let (_, batch, _) = super::super::stage_c_tests::fixture("S06_seed.fna");
        let last = batch.contexts.last().unwrap();
        assert!(last.length > 0);
        let buflen = batch.batch_buffer_length();
        // `query->sequence` is `calloc(buflen + 1)` from its first NULLB: the last context, its
        // NULLB and calloc's byte.
        assert_eq!(buflen, last.offset + last.length + 2);
        assert_eq!(buflen, batch.sequence_start.len());
        let fasta = || -> Result<(usize, &'static str)> { Ok((0, "subject")) };
        let unasked = || -> Result<(usize, &'static str)> { panic!("asked for the bytes after") };
        let refused = |result: Result<()>, buffer: &str| {
            let error = result.unwrap_err();
            assert!(error.is::<OneHitWindowRefusal>());
            assert!(
                error
                    .to_string()
                    .contains(&format!("reads past the end of the {buffer} ")),
                "{error}"
            );
        };
        // Query: the window's last letter on calloc's byte, or one past it.
        refuse_window_past_buffers(&batch, buflen - 11, 30, 0, &unasked).unwrap();
        refused(
            refuse_window_past_buffers(&batch, buflen - 10, 30, 0, &unasked),
            "query batch",
        );
        // A subject of 30 letters: its sentinel at 30 is the window's last byte from 20.
        refuse_window_past_buffers(&batch, 0, 30, 20, &unasked).unwrap();
        refused(
            refuse_window_past_buffers(&batch, 0, 30, 21, &fasta),
            "subject",
        );
        // A database sequence followed by 3 bytes in its volume.
        let volume = || -> Result<(usize, &'static str)> { Ok((3, "database volume")) };
        for s_start in 21..=23 {
            refuse_window_past_buffers(&batch, 0, 30, s_start, &volume).unwrap();
        }
        refused(
            refuse_window_past_buffers(&batch, 0, 30, 24, &volume),
            "database volume",
        );
        // The query's check comes first.
        refused(
            refuse_window_past_buffers(&batch, buflen, 30, 24, &volume),
            "query batch",
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3405-3437 (quoted at
    // `gapped_start_negative`)
    /// The gapped start of a one-hit HSP of negative length inside NCBI's buffers: `q_start + 5`
    /// for a window inside both sequences that sums above 0, `q_start` for one that sums to 0 or
    /// less or holds a sequence's NULLB (whatever follows it).
    #[test]
    fn a_one_hit_gapped_start_of_negative_length_scores_its_first_window() {
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        let (_, batch, _) = super::super::stage_c_tests::fixture("S06_seed.fna");
        assert!(batch.batch_buffer_length() > 40);
        let fasta = || -> Result<(usize, &'static str)> { Ok((0, "subject")) };
        let start = |query: &[u8], subject: &[u8], q_start: usize, s_start: usize| {
            gapped_start_negative(
                &batch,
                0,
                query,
                subject,
                q_start,
                s_start,
                ScoringMatrix::Blosum62,
                &fasta,
            )
            .unwrap()
        };
        let code = |letters: &[u8]| -> Vec<u8> {
            letters.iter().map(|&l| aa_char_to_ncbistdaa(l)).collect()
        };
        // BLOSUM62: W/W 11, W/G -2.
        let query = code(b"MKWWWWWWWWWWWK");
        let subject = code(b"WWWWWWWWWWWWW");
        assert_eq!(start(&query, &subject, 2, 0), 2 + 5);
        // The subject ends one letter before the window's end: its NULLB is in the window.
        assert_eq!(start(&query, &subject[..10], 2, 0), 2);
        // The window passes the query's end (its context's NULLB).
        assert_eq!(start(&query, &subject, 4, 0), 4);
        // Inside both, sums of -22, -9 (one W/W) and 4 (two).
        assert_eq!(start(&query, &code(b"GGGGGGGGGGGG"), 2, 0), 2);
        assert_eq!(start(&query, &code(b"WGGGGGGGGGGG"), 2, 0), 2);
        assert_eq!(start(&query, &code(b"WWGGGGGGGGGG"), 2, 0), 2 + 5);
        // From s_start on the subject side.
        assert_eq!(start(&query, &code(b"GGWWGGGGGGGGGG"), 2, 2), 2 + 5);
    }
}
