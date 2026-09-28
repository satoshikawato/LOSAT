//! BLASTX ungapped even-gap linking, with protein-subject lengths in AA.
use super::{
    linking::{LinkedHsp, LinkedHspList},
    parameters::ContextParameters,
    preliminary::{compare_score, PreliminaryHsp},
    query_setup::PreparedQueryBatch,
};
use crate::stats::sum_statistics::{gap_decay_divisor, large_gap_sum_e, small_gap_sum_e};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:998-1097
// ```c++
// CalculateLinkHSPCutoffs(EBlastProgramType program, BlastQueryInfo* query_info,
//    const BlastScoreBlk* sbp, BlastLinkHSPParameters* link_hsp_params,
//    const BlastInitialWordParameters* word_params,
//    Int8 db_length, Int4 subject_length)
// {
//     Blast_KarlinBlk* kbp = NULL;
//     double gap_prob, gap_decay_rate, x_variable, y_variable;
//     Int4 expected_length, window_size, query_length;
//     Int8 search_sp;
//     const double kEpsilon = 1.0e-9;
//
//     if (!link_hsp_params)
//         return;
//
//     /* Get KarlinBlk for context with smallest lambda (still greater than zero) */
//     s_BlastFindSmallestLambda(sbp->kbp, query_info, &kbp);
//     if (!kbp)
//         return;
//
//     window_size
//         = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
//     gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
//     gap_decay_rate = link_hsp_params->gap_decay_rate;
//     /* Use average query length */
//
//     query_length =
//         (query_info->contexts[query_info->last_context].query_offset +
//         query_info->contexts[query_info->last_context].query_length - 1)
//         / (query_info->last_context + 1);
//
//     if (Blast_SubjectIsTranslated(program) || program == eBlastTypeRpsTblastn) {
//         /* Lengths in subsequent calculations should be on the protein scale */
//         subject_length /= CODON_LENGTH;
//         db_length /= CODON_LENGTH;
//     }
//
//
//     /* Subtract off the expected score. */
//    expected_length = (Int4)BLAST_Nint(log(kbp->K*((double) query_length)*
//                                     ((double) subject_length))/(kbp->H));
//    query_length = query_length - expected_length;
//
//    subject_length = subject_length - expected_length;
//    query_length = MAX(query_length, 1);
//    subject_length = MAX(subject_length, 1);
//
//    /* If this is a database search, use database length, else the single
//       subject sequence length */
//    if (db_length > subject_length) {
//       y_variable = log((double) (db_length)/(double) subject_length)*(kbp->K)/
//          (gap_decay_rate);
//    } else {
//       y_variable = log((double) (subject_length + expected_length)/
//                        (double) subject_length)*(kbp->K)/(gap_decay_rate);
//    }
//
//    search_sp = ((Int8) query_length)* ((Int8) subject_length);
//    x_variable = 0.25*y_variable*((double) search_sp);
//
//    /* To use "small" gaps the query and subject must be "large" compared to
//       the gap size. If small gaps may be used, then the cutoff values must be
//       adjusted for the "bayesian" possibility that both large and small gaps
//       are being checked for. */
//
//    if (search_sp > 8*window_size*window_size) {
//       x_variable /= (1.0 - gap_prob + kEpsilon);
//       link_hsp_params->cutoff_big_gap =
//          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
//       x_variable = y_variable*(window_size*window_size);
//       x_variable /= (gap_prob + kEpsilon);
//       link_hsp_params->cutoff_small_gap =
//          MAX(word_params->cutoff_score_min,
//              (Int4) floor((log(x_variable)/kbp->Lambda)) + 1);
//    } else {
//       link_hsp_params->cutoff_big_gap =
//          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
//       /* The following is equivalent to forcing small gap rule to be ignored
//          when linking HSPs. */
//       link_hsp_params->gap_prob = 0;
//       link_hsp_params->cutoff_small_gap = 0;
//    }
//
//    link_hsp_params->cutoff_big_gap *= (Int4)sbp->scale_factor;
//    link_hsp_params->cutoff_small_gap *= (Int4)sbp->scale_factor;
// }
//
// /* For debugging within the C core -RMH- */
// void printBlastInitialWordParamters ( BlastInitialWordParameters *word_params, BlastQueryInfo *query_info )
// {
//   int context;
//   printf("BlastInitialWordParamters:\n");
//   printf("  x_dropoff_max = %d\n", word_params->x_dropoff_max );
//   printf("  cutoff_score_min = %d\n", word_params->cutoff_score_min );
//   printf("  cutoffs:\n");
//   for (context = query_info->first_context;
//        context <= query_info->last_context; ++context)
//   {
//     if (!(query_info->contexts[context].is_valid))
//       continue;
//     printf("    %d x_dropoff_init = %d\n", context, word_params->cutoffs[context].x_dropoff_init );
// ```
pub(crate) fn cutoffs(
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    subject_length: i32,
    db_length: i64,
) -> ([i32; 2], f64) {
    let ka = params
        .iter()
        .filter(|p| p.valid && p.ungapped.lambda > 0.0)
        .min_by(|a, b| a.ungapped.lambda.partial_cmp(&b.ungapped.lambda).unwrap())
        .expect("valid context")
        .ungapped;
    let last = batch.contexts.last().unwrap();
    let mut query_length =
        (last.offset as i32 + last.length as i32 - 1) / batch.contexts.len() as i32;
    let value = (ka.k * query_length as f64 * subject_length as f64).ln() / ka.h;
    let expected = (value + if value >= 0.0 { 0.5 } else { -0.5 }) as i32;
    query_length = (query_length - expected).max(1);
    let subject = (subject_length - expected).max(1);
    let y = if db_length > subject as i64 {
        (db_length as f64 / subject as f64).ln() * ka.k / 0.5
    } else {
        ((subject + expected) as f64 / subject as f64).ln() * ka.k / 0.5
    };
    let search = query_length as i64 * subject as i64;
    let mut x = 0.25 * y * search as f64;
    if search > 8 * 50 * 50 {
        x /= 1.0 - 0.5 + 1.0e-9;
        let big = (x.ln() / ka.lambda).floor() as i32 + 1;
        x = y * (50 * 50) as f64;
        x /= 0.5 + 1.0e-9;
        let small = ((x.ln() / ka.lambda).floor() as i32 + 1).max(
            params
                .iter()
                .filter(|p| p.valid)
                .map(|p| p.word_cutoff)
                .min()
                .unwrap(),
        );
        ([small, big], 0.5)
    } else {
        ([0, (x.ln() / ka.lambda).floor() as i32 + 1], 0.0)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:53-107
// ```c++
// /* Forward declaration */
// struct LinkHSPStruct;
//
// /** The following structure is used in "link_hsps" to decide between
//  * two different "gapping" models.  Here link is used to hook up
//  * a chain of HSP's, num is the number of links, and sum is the sum score.
//  * Once the best gapping model has been found, this information is
//  * transferred up to the LinkHSPStruct.  This structure should not be
//  * used outside of the function Blast_EvenGapLinkHSPs.
//  */
// typedef struct BlastHSPLink {
//    struct LinkHSPStruct* link[eOrderingMethods]; /**< Best
//                                                choice of HSP to link with */
//    Int2 num[eOrderingMethods]; /**< number of HSP in the ordering. */
//    Int4 sum[eOrderingMethods]; /**< Sum-Score of HSP. */
//    double xsum[eOrderingMethods]; /**< Sum-Score of HSP,
//                                      multiplied by the appropriate Lambda. */
//    Int4 changed; /**< Has the link been changed since previous access? */
// } BlastHSPLink;
//
// /** Structure containing all information internal to the process of linking
//  * HSPs.
//  */
// typedef struct LinkHSPStruct {
//    BlastHSP* hsp;      /**< Specific HSP this structure corresponds to */
//    struct LinkHSPStruct* prev;         /**< Previous HSP in a set, if any */
//    struct LinkHSPStruct* next;         /**< Next HSP in a set, if any */
//    BlastHSPLink  hsp_link; /**< Auxiliary structure for keeping track of sum
//                               scores, etc. */
//    Boolean linked_set;     /**< Is this HSp part of a linked set? */
//    Boolean start_of_chain; /**< If TRUE, this HSP starts a chain along the
//                               "link" pointer. */
//    Int4 linked_to;         /**< Where this HSP is linked to? */
//    double xsum;              /**< Normalized score of a set of HSPs */
//    ELinkOrderingMethod ordering_method;   /**< Which method (max or
//                                             no max for gaps) was
//                                             used for linking HSPs? */
//    Int4 q_offset_trim;     /**< Start of trimmed hsp in query */
//    Int4 q_end_trim;        /**< End of trimmed HSP in query */
//    Int4 s_offset_trim;     /**< Start of trimmed hsp in subject */
//    Int4 s_end_trim;        /**< End of trimmed HSP in subject */
// } LinkHSPStruct;
//
// /** The helper array contains the info used frequently in the inner
//  * for loops of the HSP linking algorithm.
//  * One array of helpers will be allocated for each thread.
//  */
// typedef struct LinkHelpStruct {
//   LinkHSPStruct* ptr;         /**< The HSP to which the info belongs */
//   Int4 q_off_trim;            /**< query start of trimmed HSP */
//   Int4 s_off_trim;            /**< subject start of trimmed HSP */
//   Int4 sum[eOrderingMethods]; /**< raw score of linked set containing HSP(?) */
//   Int4 maxsum1;               /**< threshold for stopping link attempts (?) */
//   Int4 next_larger;           /**< offset into array of HelpStructs containing
//                                    HSP with higher score, used in bailout
// ```
#[derive(Clone, Default)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-68
// ```c++
// typedef struct BlastHSPLink {
//    struct LinkHSPStruct* link[eOrderingMethods]; /**< Best
//                                                choice of HSP to link with */
//    Int2 num[eOrderingMethods]; /**< number of HSP in the ordering. */
//    Int4 sum[eOrderingMethods]; /**< Sum-Score of HSP. */
//    double xsum[eOrderingMethods]; /**< Sum-Score of HSP,
// ```
struct Path {
    sum: [i32; 2],
    num: [i16; 2],
    xsum: [f64; 2],
    next: [Option<usize>; 2],
    changed: bool,
    linked_to: i32,
    method: usize,
    linked: bool,
    head: bool,
}
#[derive(Clone, Copy, Default)]
struct Helper {
    id: Option<usize>,
    sum: [i32; 2],
    qtrim: i32,
    strim: i32,
    next_larger: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:415-588
// ```c++
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
//    const BlastQueryInfo* query_info, Int4 subject_length,
//    const BlastScoreBlk* sbp, const BlastLinkHSPParameters* link_hsp_params,
//    Boolean gapped_calculation)
// {
// 	LinkHSPStruct* H,* H2,* best[2],* first_hsp,* last_hsp,** hp_frame_start;
// 	LinkHSPStruct* hp_start = NULL;
//    BlastHSP* hsp;
//    BlastHSP** hsp_array;
// 	Blast_KarlinBlk** kbp;
// 	Int4 maxscore, cutoff[2];
// 	Boolean linked_set, ignore_small_gaps;
// 	double gap_decay_rate, gap_prob, prob[2];
// 	Int4 index, index1, num_links, frame_index;
//    ELinkOrderingMethod ordering_method;
//    Int4 num_query_frames, num_subject_frames;
// 	Int4 *hp_frame_number;
// 	Int4 window_size, trim_size;
//    Int4 number_of_hsps, total_number_of_hsps;
//    Int4 query_length, length_adjustment;
//    Int4 subject_length_orig = subject_length;
// 	LinkHSPStruct* link;
// 	Int4 H2_index,H_index;
// 	Int4 i;
//  	Int4 path_changed;  /* will be set if an element is removed that may change an existing path */
//  	Int4 first_pass, use_current_max;
// 	LinkHelpStruct *lh_helper=0;
//    Int4 lh_helper_size;
// 	Int4 query_context; /* AM: to support query concatenation. */
//    const Boolean kTranslatedQuery = Blast_QueryIsTranslated(program_number);
//    LinkHSPStruct** link_hsp_array;
//
// 	if (hsp_list == NULL)
// 		return -1;
//
//    hsp_array = hsp_list->hsp_array;
//
//    lh_helper_size = MAX(1024,hsp_list->hspcnt+5);
//    lh_helper = (LinkHelpStruct *)
//       calloc(lh_helper_size, sizeof(LinkHelpStruct));
//
// 	if (gapped_calculation)
// 		kbp = sbp->kbp_gap;
// 	else
// 		kbp = sbp->kbp;
//
// 	total_number_of_hsps = hsp_list->hspcnt;
//
// 	number_of_hsps = total_number_of_hsps;
//
//    /* For convenience, include overlap size into the gap size */
// 	window_size = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
//    trim_size = (link_hsp_params->overlap_size + 1) / 2;
// 	gap_prob = link_hsp_params->gap_prob;
// 	gap_decay_rate = link_hsp_params->gap_decay_rate;
//
//    if (Blast_SubjectIsTranslated(program_number))
//       num_subject_frames = NUM_STRANDS;
//    else
//       num_subject_frames = 1;
//
//    link_hsp_array =
//       (LinkHSPStruct**) malloc(total_number_of_hsps*sizeof(LinkHSPStruct*));
//    for (index = 0; index < total_number_of_hsps; ++index) {
//       link_hsp_array[index] = (LinkHSPStruct*) calloc(1, sizeof(LinkHSPStruct));
//       link_hsp_array[index]->hsp = hsp_array[index];
//    }
//
//    /* Sort by (reverse) position. */
//    if (kTranslatedQuery) {
//       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//             s_RevCompareHSPsTbx);
//    } else {
//       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//             s_RevCompareHSPsTbn);
//    }
//
//    cutoff[0] = link_hsp_params->cutoff_small_gap;
//    cutoff[1] = link_hsp_params->cutoff_big_gap;
//
//    ignore_small_gaps = (cutoff[0] == 0);
//
//    /* If query is nucleotide, it has 2 strands that should be separated. */
//    if (Blast_QueryIsNucleotide(program_number))
//       num_query_frames = NUM_STRANDS*query_info->num_queries;
//    else
//       num_query_frames = query_info->num_queries;
//
//    hp_frame_start =
//        calloc(num_subject_frames*num_query_frames, sizeof(LinkHSPStruct*));
//    hp_frame_number = calloc(num_subject_frames*num_query_frames, sizeof(Int4));
//
//    /* hook up the HSP's */
//    hp_frame_start[0] = link_hsp_array[0];
//
//    /* Put entries from different strands into separate 'query_frame's. */
//    {
//       Int4 cur_frame=0;
//       Int4 strand_factor = (kTranslatedQuery ? 3 : 1);
//       for (index=0;index<number_of_hsps;index++)
//       {
//         H=link_hsp_array[index];
//         H->start_of_chain = FALSE;
//         hp_frame_number[cur_frame]++;
//
//         H->prev= index ? link_hsp_array[index-1] : NULL;
//         H->next= index<(number_of_hsps-1) ? link_hsp_array[index+1] : NULL;
//         if (H->prev != NULL &&
//             ((H->hsp->context/strand_factor) !=
// 	     (H->prev->hsp->context/strand_factor) ||
//              (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
//         { /* If frame switches, then start new list. */
//            hp_frame_number[cur_frame]--;
//            hp_frame_number[++cur_frame]++;
//            hp_frame_start[cur_frame] = H;
//            H->prev->next = NULL;
//            H->prev = NULL;
//         }
//       }
//       num_query_frames = cur_frame+1;
//    }
//
//    /* trim_size is the maximum amount q.offset can differ from
//       q.offset_trim */
//    /* This is used to break out of H2 loop early */
//    for (index=0;index<number_of_hsps;index++)
//    {
//        Int4 q_length, s_length;
//        H = link_hsp_array[index];
//        hsp = H->hsp;
//        q_length = (hsp->query.end - hsp->query.offset) / 4;
//        s_length = (hsp->subject.end - hsp->subject.offset) / 4;
//        H->q_offset_trim = hsp->query.offset + MIN(q_length, trim_size);
//        H->q_end_trim = hsp->query.end - MIN(q_length, trim_size);
//        H->s_offset_trim = hsp->subject.offset + MIN(s_length, trim_size);
//        H->s_end_trim = hsp->subject.end - MIN(s_length, trim_size);
//    }
//
// 	for (frame_index=0; frame_index<num_query_frames; frame_index++)
// 	{
//       hp_start = s_LinkHSPStructReset(hp_start);
//       hp_start->next = hp_frame_start[frame_index];
//       hp_frame_start[frame_index]->prev = hp_start;
//       number_of_hsps = hp_frame_number[frame_index];
//       query_context = hp_start->next->hsp->context;
//       length_adjustment = query_info->contexts[query_context].length_adjustment;
//       query_length = query_info->contexts[query_context].query_length;
//       query_length = MAX(query_length - length_adjustment, 1);
//       subject_length = subject_length_orig; /* in nucleotides even for tblast[nx] */
//       /* If subject is translated, length adjustment is given in nucleotide
//          scale. */
//       if (Blast_SubjectIsTranslated(program_number))
//       {
//          length_adjustment /= CODON_LENGTH;
//          subject_length /= CODON_LENGTH;
//       }
//       subject_length = MAX(subject_length - length_adjustment, 1);
//
//       lh_helper[0].ptr = hp_start;
//       lh_helper[0].q_off_trim = 0;
//       lh_helper[0].s_off_trim = 0;
//       lh_helper[0].maxsum1  = -10000;
//       lh_helper[0].next_larger  = 0;
//
//       /* lh_helper[0]  = empty     = additional end marker
//        * lh_helper[1]  = hsp_start = empty entry used in original code
//        * lh_helper[2]  = hsp_array->next = hsp_array[0]
//        * lh_helper[i]  = ... = hsp_array[i-2] (for i>=2)
//        */
//       first_pass=1;    /* do full search */
//       path_changed=1;
//       for (H=hp_start->next; H!=NULL; H=H->next)
//          H->hsp_link.changed=1;
//
// ```
pub(crate) fn link_even(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    subject_length: i32,
    prepared: ([i32; 2], f64),
) -> LinkedHspList {
    if input.is_empty() {
        return LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        };
    }
    let (cutoff, gap_prob) = prepared;
    let ignore = cutoff[0] == 0;
    let mut order: Vec<_> = (0..input.len()).collect();
    order.sort_by(|&i, &j| {
        let a = &input[i];
        let b = &input[j];
        (a.context / 3)
            .cmp(&(b.context / 3))
            .then_with(|| b.q_start.cmp(&a.q_start))
            .then_with(|| b.q_end.cmp(&a.q_end))
            .then_with(|| b.s_start.cmp(&a.s_start))
            .then_with(|| b.s_end.cmp(&a.s_end))
    });
    let mut paths = vec![Path::default(); input.len()];
    let mut evalues = vec![0.0; input.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && input[order[start]].context / 3 == input[order[end]].context / 3
        {
            end += 1;
        }
        let mut active = order[start..end].to_vec();
        let ctx = input[active[0]].context;
        let adjustment = params[ctx].length_adjustment;
        let qlen = (batch.contexts[ctx].length as i32 - adjustment).max(1);
        let slen = (subject_length - adjustment).max(1);
        for &id in &active {
            paths[id].changed = true;
        }
        let mut first = true;
        let mut path_changed = true;
        while !active.is_empty() {
            let mut best = [None, None];
            let mut use_current = false;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:598-652
            // ```c++
            //           *  - Find the max paths (based on old scores).
            //           *  - If no paths were changed by removal of nodes (ie research==0)
            //           *    then these max paths are still the best.
            //           *  - else if these max paths were unchanged, then they are still the best.
            //           */
            //          use_current_max=0;
            //          if (!first_pass){
            //             Int4 max0,max1;
            //             /* Find the current max sums */
            //             if(!ignore_small_gaps){
            //                max0 = -cutoff[0];
            //                max1 = -cutoff[1];
            //                for (H=hp_start->next; H!=NULL; H=H->next) {
            //                   Int4 sum0=H->hsp_link.sum[0];
            //                   Int4 sum1=H->hsp_link.sum[1];
            //                   if(sum0>=max0)
            //                   {
            //                      max0=sum0;
            //                      best[0]=H;
            //                   }
            //                   if(sum1>=max1)
            //                   {
            //                      max1=sum1;
            //                      best[1]=H;
            //                   }
            //                }
            //             } else {
            //                maxscore = -cutoff[1];
            //                for (H=hp_start->next; H!=NULL; H=H->next) {
            //                   Int4  sum=H->hsp_link.sum[1];
            //                   if(sum>=maxscore)
            //                   {
            //                      maxscore=sum;
            //                      best[1]=H;
            //                   }
            //                }
            //             }
            //             if(path_changed==0){
            //                /* No path was changed, use these max sums. */
            //                use_current_max=1;
            //             }
            //             else{
            //                /* If max path hasn't chaged, we can use it */
            //                /* Walk down best, give up if we find a removed item in path */
            //                use_current_max=1;
            //                if(!ignore_small_gaps){
            //                   for (H=best[0]; H!=NULL; H=H->hsp_link.link[0])
            //                      if (H->linked_to==-1000) {use_current_max=0; break;}
            //                }
            //                if(use_current_max)
            //                   for (H=best[1]; H!=NULL; H=H->hsp_link.link[1])
            //                      if (H->linked_to==-1000) {use_current_max=0; break;}
            //
            //             }
            //          }
            // ```
            if !first {
                for method in if ignore { 1..2 } else { 0..2 } {
                    let mut max = -cutoff[method];
                    for &id in &active {
                        if paths[id].sum[method] >= max {
                            max = paths[id].sum[method];
                            best[method] = Some(id);
                        }
                    }
                }
                use_current = !path_changed;
                if path_changed {
                    use_current = true;
                    for method in if ignore { 1..2 } else { 0..2 } {
                        let mut node = best[method];
                        while let Some(id) = node {
                            if paths[id].linked_to == -1000 {
                                use_current = false;
                                break;
                            }
                            node = paths[id].next[method];
                        }
                        if !use_current {
                            break;
                        }
                    }
                }
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:660-697
            // ```c++
            //             for (H=hp_start,H_index=1; H!=NULL; H=H->next,H_index++) {
            //                Int4 s_frame = H->hsp->subject.frame;
            //                Int4 s_off_t = H->s_offset_trim;
            //                Int4 q_off_t = H->q_offset_trim;
            //                lh_helper[H_index].ptr = H;
            //                lh_helper[H_index].q_off_trim = q_off_t;
            //                lh_helper[H_index].s_off_trim = s_off_t;
            //                for(i=0;i<eOrderingMethods;i++)
            //                   lh_helper[H_index].sum[i] = H->hsp_link.sum[i];
            //                max[SIGN(s_frame)+1]=
            //                   MAX(max[SIGN(s_frame)+1],H->hsp_link.sum[1]);
            //                lh_helper[H_index].maxsum1 =max[SIGN(s_frame)+1];
            //
            //                /* set next_larger to link back to closest entry with a sum1
            //                   larger than this */
            //                {
            //                   Int4 cur_sum=lh_helper[H_index].sum[1];
            //                   Int4 prev = H_index-1;
            //                   Int4 prev_sum = lh_helper[prev].sum[1];
            //                   while((cur_sum>=prev_sum) && (prev>0)){
            //                      prev=lh_helper[prev].next_larger;
            //                      prev_sum = lh_helper[prev].sum[1];
            //                   }
            //                   lh_helper[H_index].next_larger = prev;
            //                }
            //                H->linked_to = 0;
            //             }
            //
            //             lh_helper[1].maxsum1 = -10000;
            //
            //             /****** loop iter for index = 0  **************************/
            //             if(!ignore_small_gaps)
            //             {
            //                index=0;
            //                maxscore = -cutoff[index];
            //                H_index = 2;
            //                for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
            //                {
            // ```
            if !use_current {
                let mut help = vec![Helper::default(); active.len() + 2];
                help[0].sum[1] = -10000;
                for (pos, &id) in active.iter().enumerate() {
                    let h = &input[id];
                    let index = pos + 2;
                    help[index] = Helper {
                        id: Some(id),
                        sum: paths[id].sum,
                        qtrim: h.q_start + ((h.q_end - h.q_start) / 4).min(5),
                        strim: h.s_start + ((h.s_end - h.s_start) / 4).min(5),
                        next_larger: 0,
                    };
                    let mut prev = index - 1;
                    while prev > 0 && help[index].sum[1] >= help[prev].sum[1] {
                        prev = help[prev].next_larger;
                    }
                    help[index].next_larger = prev;
                    paths[id].linked_to = 0;
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:699-772
                // ```c++
                //                   Int4 H_hsp_sum=0;
                //                   double H_hsp_xsum=0.0;
                //                   LinkHSPStruct* H_hsp_link=NULL;
                //                   if (H->hsp->score > cutoff[index]) {
                //                      Int4 H_query_etrim = H->q_end_trim;
                //                      Int4 H_sub_etrim = H->s_end_trim;
                //                      Int4 H_q_et_gap = H_query_etrim+window_size;
                //                      Int4 H_s_et_gap = H_sub_etrim+window_size;
                //
                //                      /* We only walk down hits with the same frame sign */
                //                      /* for (H2=H->prev; H2!=NULL; H2=H2->prev,H2_index--) */
                //                      for (H2_index=H_index-1; H2_index>1; H2_index=H2_index-1)
                //                      {
                //                         Int4 b1,b2,b4,b5;
                //                         Int4 q_off_t,s_off_t,sum;
                //
                //                         /* s_frame = lh_helper[H2_index].s_frame; */
                //                         q_off_t = lh_helper[H2_index].q_off_trim;
                //                         s_off_t = lh_helper[H2_index].s_off_trim;
                //
                //                         /* combine tests to reduce mispredicts -cfj */
                //                         b1 = q_off_t <= H_query_etrim;
                //                         b2 = s_off_t <= H_sub_etrim;
                //                         sum = lh_helper[H2_index].sum[index];
                //
                //
                //                         b4 = ( q_off_t > H_q_et_gap ) ;
                //                         b5 = ( s_off_t > H_s_et_gap ) ;
                //
                //                         /* list is sorted by q_off, so q_off should only increase.
                //                          * q_off_t can only differ from q_off by trim_size
                //                          * So once q_off_t is large enough (ie it exceeds limit
                //                          * by trim_size), we can stop.  -cfj
                //                          */
                //                         if(q_off_t > (H_q_et_gap+trim_size))
                //                            break;
                //
                //                         if (b1|b2|b5|b4) continue;
                //
                //                         if (sum>H_hsp_sum)
                //                         {
                //                            H2=lh_helper[H2_index].ptr;
                //                            H_hsp_num=H2->hsp_link.num[index];
                //                            H_hsp_sum=H2->hsp_link.sum[index];
                //                            H_hsp_xsum=H2->hsp_link.xsum[index];
                //                            H_hsp_link=H2;
                //                         }
                //                      } /* end for H2... */
                //                   }
                //                   {
                //                      Int4 score=H->hsp->score;
                //                      double new_xsum =
                //                        H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
                //                        kbp[H->hsp->context]->logK;
                //                      Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
                //
                //                      H->hsp_link.sum[index] = new_sum;
                //                      H->hsp_link.num[index] = H_hsp_num+1;
                //                      H->hsp_link.link[index] = H_hsp_link;
                //                      lh_helper[H_index].sum[index] = new_sum;
                //                      if (new_sum >= maxscore)
                //                      {
                //                         maxscore=new_sum;
                //                         best[index]=H;
                //                      }
                //                      H->hsp_link.xsum[index] = new_xsum;
                //                      if(H_hsp_link)
                //                         ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                //                   }
                //                } /* end for H=... */
                //             }
                //             /****** loop iter for index = 1  **************************/
                //             index=1;
                //             maxscore = -cutoff[index];
                // ```
                if !ignore {
                    let mut max = -cutoff[0];
                    for (pos, &id) in active.iter().enumerate() {
                        let h = &input[id];
                        let index = pos + 2;
                        let qend = h.q_end - ((h.q_end - h.q_start) / 4).min(5);
                        let send = h.s_end - ((h.s_end - h.s_start) / 4).min(5);
                        let mut link = None;
                        let mut sum = 0;
                        if h.score > cutoff[0] {
                            for prev in (2..index).rev() {
                                let p = &help[prev];
                                if p.qtrim > qend + 50 + 5 {
                                    break;
                                }
                                if p.qtrim <= qend
                                    || p.strim <= send
                                    || p.qtrim > qend + 50
                                    || p.strim > send + 50
                                {
                                    continue;
                                }
                                if p.sum[0] > sum {
                                    sum = p.sum[0];
                                    link = p.id;
                                }
                            }
                        }
                        let num = link.map_or(0, |i| i32::from(paths[i].num[0]));
                        let xsum = link.map_or(0.0, |i| paths[i].xsum[0]);
                        let ka = params[h.context].ungapped;
                        paths[id].sum[0] = sum + (h.score - cutoff[0]);
                        paths[id].num[0] = (num + 1) as i16;
                        paths[id].next[0] = link;
                        paths[id].xsum[0] = xsum + h.score as f64 * ka.lambda - ka.k.ln();
                        help[index].sum[0] = paths[id].sum[0];
                        if paths[id].sum[0] >= max {
                            max = paths[id].sum[0];
                            best[0] = Some(id);
                        }
                        if let Some(i) = link {
                            paths[i].linked_to += 1;
                        }
                    }
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:774-895
                // ```c++
                //             for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
                //             {
                //                Int4 H_hsp_num=0;
                //                Int4 H_hsp_sum=0;
                //                double H_hsp_xsum=0.0;
                //                LinkHSPStruct* H_hsp_link=NULL;
                //
                //                H->hsp_link.changed=1;
                //                H2 = H->hsp_link.link[index];
                //                if ((!first_pass) && ((H2==0) || (H2->hsp_link.changed==0)))
                //                {
                //                   /* If the best choice last time has not been changed, then
                //                      it is still the best choice, so no need to walk down list.
                //                   */
                //                   if(H2){
                //                      H_hsp_num=H2->hsp_link.num[index];
                //                      H_hsp_sum=H2->hsp_link.sum[index];
                //                      H_hsp_xsum=H2->hsp_link.xsum[index];
                //                   }
                //                   H_hsp_link=H2;
                //                   H->hsp_link.changed=0;
                //                } else if (H->hsp->score > cutoff[index]) {
                //                   Int4 H_query_etrim = H->q_end_trim;
                //                   Int4 H_sub_etrim = H->s_end_trim;
                //
                //
                //                   /* Here we look at what was the best choice last time, if it's
                //                    * still around, and set this to the initial choice. By
                //                    * setting the best score to a (potentially) large value
                //                    * initially, we can reduce the number of hsps checked.
                //                    */
                //
                //                   /* Currently we set the best score to a value just less than
                //                    * the real value. This is not really necessary, but doing
                //                    * this ensures that in the case of a tie, we make the same
                //                    * selection the original code did.
                //                    */
                //
                //                   if(!first_pass&&H2&&H2->linked_to>=0){
                //                      if(1){
                //                         /* We set this to less than the real value to keep the
                //                            original ordering in case of ties. */
                //                         H_hsp_sum=H2->hsp_link.sum[index]-1;
                //                      }else{
                //                         H_hsp_num=H2->hsp_link.num[index];
                //                         H_hsp_sum=H2->hsp_link.sum[index];
                //                         H_hsp_xsum=H2->hsp_link.xsum[index];
                //                         H_hsp_link=H2;
                //                      }
                //                   }
                //
                //                   /* We now only walk down hits with the same frame sign */
                //                   /* for (H2=H->prev; H2!=NULL; H2=H2->prev,H2_index--) */
                //                   for (H2_index=H_index-1; H2_index>1;)
                //                   {
                //                      Int4 b0,b1,b2;
                //                      Int4 q_off_t,s_off_t,sum,next_larger;
                //                      LinkHelpStruct * H2_helper=&lh_helper[H2_index];
                //                      sum = H2_helper->sum[index];
                //                      next_larger = H2_helper->next_larger;
                //
                //                      s_off_t = H2_helper->s_off_trim;
                //                      q_off_t = H2_helper->q_off_trim;
                //
                //                      b0 = sum <= H_hsp_sum;
                //
                //                      /* Compute the next H2_index */
                //                      H2_index--;
                //                      if(b0){	 /* If this sum is too small to beat H_hsp_sum, advance to a larger sum */
                //                         H2_index=next_larger;
                //                      }
                //
                //                      /* combine tests to reduce mispredicts -cfj */
                //                      b1 = q_off_t <= H_query_etrim;
                //                      b2 = s_off_t <= H_sub_etrim;
                //
                //                      if(0) if(H2_helper->maxsum1<=H_hsp_sum)break;
                //
                //                      if (!(b0|b1|b2) )
                //                      {
                //                         H2 = H2_helper->ptr;
                //
                //                         H_hsp_num=H2->hsp_link.num[index];
                //                         H_hsp_sum=H2->hsp_link.sum[index];
                //                         H_hsp_xsum=H2->hsp_link.xsum[index];
                //                         H_hsp_link=H2;
                //                      }
                //                   } /* end for H2_index... */
                //                } /* end if(H->score>cuttof[]) */
                //                {
                //                   Int4 score=H->hsp->score;
                //                   double new_xsum =
                //                      H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
                //                      kbp[H->hsp->context]->logK;
                //                   Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
                //
                //                   H->hsp_link.sum[index] = new_sum;
                //                   H->hsp_link.num[index] = H_hsp_num+1;
                //                   H->hsp_link.link[index] = H_hsp_link;
                //                   lh_helper[H_index].sum[index] = new_sum;
                //                   lh_helper[H_index].maxsum1 = MAX(lh_helper[H_index-1].maxsum1, new_sum);
                //                   /* Update this entry's 'next_larger' field */
                //                   {
                //                      Int4 cur_sum=lh_helper[H_index].sum[1];
                //                      Int4 prev = H_index-1;
                //                      Int4 prev_sum = lh_helper[prev].sum[1];
                //                      while((cur_sum>=prev_sum) && (prev>0)){
                //                         prev=lh_helper[prev].next_larger;
                //                         prev_sum = lh_helper[prev].sum[1];
                //                      }
                //                      lh_helper[H_index].next_larger = prev;
                //                   }
                //
                //                   if (new_sum >= maxscore)
                //                   {
                //                      maxscore=new_sum;
                //                      best[index]=H;
                //                   }
                //                   H->hsp_link.xsum[index] = new_xsum;
                //                   if(H_hsp_link)
                //                      ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                //                }
                // ```
                let mut max = -cutoff[1];
                for (pos, &id) in active.iter().enumerate() {
                    let h = &input[id];
                    let index = pos + 2;
                    let qend = h.q_end - ((h.q_end - h.q_start) / 4).min(5);
                    let send = h.s_end - ((h.s_end - h.s_start) / 4).min(5);
                    paths[id].changed = true;
                    let previous = paths[id].next[1];
                    let mut link = None;
                    let mut sum = 0;
                    if !first && previous.is_none_or(|i| !paths[i].changed) {
                        link = previous;
                        sum = link.map_or(0, |i| paths[i].sum[1]);
                        paths[id].changed = false;
                    } else if h.score > cutoff[1] {
                        if !first && previous.is_some_and(|i| paths[i].linked_to >= 0) {
                            sum = paths[previous.unwrap()].sum[1] - 1;
                        }
                        let mut prev = index - 1;
                        while prev > 1 {
                            let p = help[prev];
                            let too_small = p.sum[1] <= sum;
                            prev = if too_small { p.next_larger } else { prev - 1 };
                            if !too_small && p.qtrim > qend && p.strim > send {
                                sum = p.sum[1];
                                link = p.id;
                            }
                        }
                    }
                    let num = link.map_or(0, |i| i32::from(paths[i].num[1]));
                    let xsum = link.map_or(0.0, |i| paths[i].xsum[1]);
                    let ka = params[h.context].ungapped;
                    paths[id].sum[1] = sum + (h.score - cutoff[1]);
                    paths[id].num[1] = (num + 1) as i16;
                    paths[id].next[1] = link;
                    paths[id].xsum[1] = xsum + h.score as f64 * ka.lambda - ka.k.ln();
                    help[index].sum[1] = paths[id].sum[1];
                    let mut prev = index - 1;
                    while prev > 0 && help[index].sum[1] >= help[prev].sum[1] {
                        prev = help[prev].next_larger;
                    }
                    help[index].next_larger = prev;
                    if paths[id].sum[1] >= max {
                        max = paths[id].sum[1];
                        best[1] = Some(id);
                    }
                    if let Some(i) = link {
                        paths[i].linked_to += 1;
                    }
                }
                path_changed = false;
                first = false;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:898-984
            // ```c++
            //             first_pass=0;
            //          }
            //          /********************************/
            //          if (!ignore_small_gaps)
            //          {
            //             /* Select the best ordering method.
            //                First we add back in the value cutoff[index] * the number
            //                of links, as this was subtracted out for purposes of the
            //                comparison above. */
            //             best[0]->hsp_link.sum[0] +=
            //                (best[0]->hsp_link.num[0])*cutoff[0];
            //
            //             prob[0] = BLAST_SmallGapSumE(window_size,
            //                          best[0]->hsp_link.num[0], best[0]->hsp_link.xsum[0],
            //                          query_length, subject_length,
            //                          query_info->contexts[query_context].eff_searchsp,
            //                          BLAST_GapDecayDivisor(gap_decay_rate,
            //                                               best[0]->hsp_link.num[0]) );
            //
            //             /* Adjust the e-value because we are performing multiple tests */
            //             if( best[0]->hsp_link.num[0] > 1 ) {
            //               if( gap_prob == 0 || (prob[0] /= gap_prob) > INT4_MAX ) {
            //                 prob[0] = INT4_MAX;
            //               }
            //             }
            //
            //             prob[1] = BLAST_LargeGapSumE(best[1]->hsp_link.num[1],
            //                          best[1]->hsp_link.xsum[1],
            //                          query_length, subject_length,
            //                          query_info->contexts[query_context].eff_searchsp,
            //                          BLAST_GapDecayDivisor(gap_decay_rate,
            //                                               best[1]->hsp_link.num[1]));
            //
            //             if( best[1]->hsp_link.num[1] > 1 ) {
            //               if( 1 - gap_prob == 0 || (prob[1] /= 1 - gap_prob) > INT4_MAX ) {
            //                 prob[1] = INT4_MAX;
            //               }
            //             }
            //             ordering_method =
            //                prob[0]<=prob[1] ? eLinkSmallGaps : eLinkLargeGaps;
            //          }
            //          else
            //          {
            //             /* We only consider the case of big gaps. */
            //             best[1]->hsp_link.sum[1] +=
            //                (best[1]->hsp_link.num[1])*cutoff[1];
            //
            //             prob[1] = BLAST_LargeGapSumE(
            //                          best[1]->hsp_link.num[1],
            //                          best[1]->hsp_link.xsum[1],
            //                          query_length, subject_length,
            //                          query_info->contexts[query_context].eff_searchsp,
            //                          BLAST_GapDecayDivisor(gap_decay_rate,
            //                                               best[1]->hsp_link.num[1]));
            //             ordering_method = eLinkLargeGaps;
            //          }
            //
            //          best[ordering_method]->start_of_chain = TRUE;
            //          best[ordering_method]->hsp->evalue    = prob[ordering_method];
            //
            //          /* remove the links that have been ordered already. */
            //          if (best[ordering_method]->hsp_link.link[ordering_method])
            //             linked_set = TRUE;
            //          else
            //             linked_set = FALSE;
            //
            //          if (best[ordering_method]->linked_to>0) path_changed=1;
            //          for (H=best[ordering_method]; H!=NULL;
            //               H=H->hsp_link.link[ordering_method])
            //          {
            //             if (H->linked_to>1) path_changed=1;
            //             H->linked_to=-1000;
            //             H->hsp_link.changed=1;
            //             /* record whether this is part of a linked set. */
            //             H->linked_set = linked_set;
            //             H->ordering_method = ordering_method;
            //             H->hsp->evalue = prob[ordering_method];
            //             if (H->next)
            //                (H->next)->prev=H->prev;
            //             if (H->prev)
            //                (H->prev)->next=H->next;
            //             number_of_hsps--;
            //          }
            //
            //       } /* end while num_hsps... */
            // 	} /* end for frame_index ... */
            //
            // ```
            let mut probs = [0.0, 0.0];
            let method;
            let b = best[1].unwrap();
            let n = i32::from(paths[b].num[1]);
            if !ignore {
                let a = best[0].unwrap();
                let n0 = i32::from(paths[a].num[0]);
                paths[a].sum[0] += n0 * cutoff[0];
                probs[0] = small_gap_sum_e(
                    50,
                    n0 as i16,
                    paths[a].xsum[0],
                    qlen,
                    slen,
                    params[ctx].search_space,
                    gap_decay_divisor(0.5, n0 as usize),
                );
                if n0 > 1 {
                    probs[0] = if gap_prob == 0.0 {
                        i32::MAX as f64
                    } else {
                        (probs[0] / gap_prob).min(i32::MAX as f64)
                    };
                }
                probs[1] = large_gap_sum_e(
                    n as i16,
                    paths[b].xsum[1],
                    qlen,
                    slen,
                    params[ctx].search_space,
                    gap_decay_divisor(0.5, n as usize),
                );
                if n > 1 {
                    probs[1] = if 1.0 - gap_prob == 0.0 {
                        i32::MAX as f64
                    } else {
                        (probs[1] / (1.0 - gap_prob)).min(i32::MAX as f64)
                    };
                }
                method = if probs[0] <= probs[1] { 0 } else { 1 };
            } else {
                paths[b].sum[1] += n * cutoff[1];
                probs[1] = large_gap_sum_e(
                    n as i16,
                    paths[b].xsum[1],
                    qlen,
                    slen,
                    params[ctx].search_space,
                    gap_decay_divisor(0.5, n as usize),
                );
                method = 1;
            }
            let head = best[method].unwrap();
            paths[head].head = true;
            let linked = paths[head].next[method].is_some();
            if paths[head].linked_to > 0 {
                path_changed = true;
            }
            let mut node = Some(head);
            while let Some(id) = node {
                if paths[id].linked_to > 1 {
                    path_changed = true;
                }
                paths[id].linked_to = -1000;
                paths[id].changed = true;
                paths[id].linked = linked;
                paths[id].method = method;
                evalues[id] = probs[method];
                node = paths[id].next[method];
            }
            active.retain(|id| paths[*id].linked_to != -1000);
        }
        start = end;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-1088
    // ```c++
    //    if (kTranslatedQuery) {
    //       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_RevCompareHSPsTransl);
    //       qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_FwdCompareHSPsTransl);
    //    } else {
    //       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_RevCompareHSPs);
    //       qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_FwdCompareHSPs);
    //    }
    //
    //    /* Sort by starting position. */
    //
    //
    // 	for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
    // 	{
    // 		H = link_hsp_array[index];
    // 		H->prev = NULL;
    // 		H->next = NULL;
    // 	}
    //
    //    /* hook up the HSP's. */
    // 	first_hsp = NULL;
    // 	for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
    //    {
    // 		H = link_hsp_array[index];
    //
    //       /* If this is not a single piece or the start of a chain, then Skip it. */
    //       if (H->linked_set == TRUE && H->start_of_chain == FALSE)
    // 			continue;
    //
    //       /* If the HSP has no "link" connect the "next", otherwise follow the "link"
    //          chain down, connecting them with "next" and "prev". */
    // 		if (last_hsp == NULL)
    // 			first_hsp = H;
    // 		H->prev = last_hsp;
    // 		ordering_method = H->ordering_method;
    // 		if (H->hsp_link.link[ordering_method] == NULL)
    // 		{
    //          /* Grab the next HSP that is not part of a chain or the start of a chain */
    //          /* The "next" pointers are not hooked up yet in HSP's further down array. */
    //          index1=index;
    //          H2 = index1<(total_number_of_hsps-1) ? link_hsp_array[index1+1] : NULL;
    //          while (H2 && H2->linked_set == TRUE &&
    //                 H2->start_of_chain == FALSE)
    //          {
    //             index1++;
    // 		     	H2 = index1<(total_number_of_hsps-1) ? link_hsp_array[index1+1] : NULL;
    //          }
    //          H->next= H2;
    // 		}
    // 		else
    // 		{
    // 			/* The first one has the number of links correct. */
    // 			num_links = H->hsp_link.num[ordering_method];
    // 			link = H->hsp_link.link[ordering_method];
    // 			while (link)
    // 			{
    // 				H->hsp->num = num_links;
    //                 H->xsum = H->hsp_link.xsum[ordering_method];
    // 				H->next = (LinkHSPStruct*) link;
    // 				H->prev = last_hsp;
    // 				last_hsp = H;
    // 				H = H->next;
    // 				if (H != NULL)
    // 				    link = H->hsp_link.link[ordering_method];
    // 				else
    // 				    break;
    // 			}
    // 			/* Set these for last link in chain. */
    // 			H->hsp->num = num_links;
    //             H->xsum = H->hsp_link.xsum[ordering_method];
    //          /* Grab the next HSP that is not part of a chain or the start of a chain */
    //          index1=index;
    //          H2 = index1<(total_number_of_hsps-1) ? link_hsp_array[index1+1] : NULL;
    //          while (H2 && H2->linked_set == TRUE &&
    //                 H2->start_of_chain == FALSE)
    // 		   {
    //             index1++;
    //             H2 = index1<(total_number_of_hsps-1) ? link_hsp_array[index1+1] : NULL;
    // 			}
    //          H->next= H2;
    // 			H->prev = last_hsp;
    // 		}
    // 		last_hsp = H;
    // 	}
    //
    //    /* The HSP's may be in a different order than they were before,
    //       but first_hsp contains the first one. */
    //    for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
    //       hsp_list->hsp_array[index] = H->hsp;
    //       /* Free the wrapper structure */
    //       H2 = H->next;
    //       sfree(H);
    //       H = H2;
    //    }
    //    sfree(link_hsp_array);
    //    sfree(lh_helper);
    // ```
    order.sort_by(|&i, &j| {
        let a = &input[i];
        let b = &input[j];
        (a.context / 3)
            .cmp(&(b.context / 3))
            .then_with(|| b.q_start.cmp(&a.q_start))
            .then_with(|| b.s_start.cmp(&a.s_start))
    });
    order.sort_by(|&i, &j| {
        let a = &input[i];
        let b = &input[j];
        (a.context / 3)
            .cmp(&(b.context / 3))
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| a.s_start.cmp(&b.s_start))
    });
    let mut output = Vec::new();
    for id in order {
        if paths[id].linked && !paths[id].head {
            continue;
        }
        let method = paths[id].method;
        let n = i32::from(paths[id].num[method]);
        let mut node = Some(id);
        while let Some(i) = node {
            output.push(LinkedHsp {
                hsp: input[i].clone(),
                num: if paths[id].linked { n } else { 1 },
                evalue: evalues[i],
                source_index: i,
            });
            node = paths[i].next[method];
        }
    }
    output.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
    let best_evalue = output
        .iter()
        .map(|h| h.evalue)
        .reduce(f64::min)
        .unwrap_or(0.0);
    LinkedHspList {
        hsps: output,
        best_evalue,
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:998-1097
// ```c++
// CalculateLinkHSPCutoffs(EBlastProgramType program, BlastQueryInfo* query_info,
//    const BlastScoreBlk* sbp, BlastLinkHSPParameters* link_hsp_params,
//    const BlastInitialWordParameters* word_params,
//    Int8 db_length, Int4 subject_length)
// {
//     Blast_KarlinBlk* kbp = NULL;
//     double gap_prob, gap_decay_rate, x_variable, y_variable;
//     Int4 expected_length, window_size, query_length;
//     Int8 search_sp;
//     const double kEpsilon = 1.0e-9;
//
//     if (!link_hsp_params)
//         return;
//
//     /* Get KarlinBlk for context with smallest lambda (still greater than zero) */
//     s_BlastFindSmallestLambda(sbp->kbp, query_info, &kbp);
//     if (!kbp)
//         return;
//
//     window_size
//         = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
//     gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
//     gap_decay_rate = link_hsp_params->gap_decay_rate;
//     /* Use average query length */
//
//     query_length =
//         (query_info->contexts[query_info->last_context].query_offset +
//         query_info->contexts[query_info->last_context].query_length - 1)
//         / (query_info->last_context + 1);
//
//     if (Blast_SubjectIsTranslated(program) || program == eBlastTypeRpsTblastn) {
//         /* Lengths in subsequent calculations should be on the protein scale */
//         subject_length /= CODON_LENGTH;
//         db_length /= CODON_LENGTH;
//     }
//
//
//     /* Subtract off the expected score. */
//    expected_length = (Int4)BLAST_Nint(log(kbp->K*((double) query_length)*
//                                     ((double) subject_length))/(kbp->H));
//    query_length = query_length - expected_length;
//
//    subject_length = subject_length - expected_length;
//    query_length = MAX(query_length, 1);
//    subject_length = MAX(subject_length, 1);
//
//    /* If this is a database search, use database length, else the single
//       subject sequence length */
//    if (db_length > subject_length) {
//       y_variable = log((double) (db_length)/(double) subject_length)*(kbp->K)/
//          (gap_decay_rate);
//    } else {
//       y_variable = log((double) (subject_length + expected_length)/
//                        (double) subject_length)*(kbp->K)/(gap_decay_rate);
//    }
//
//    search_sp = ((Int8) query_length)* ((Int8) subject_length);
//    x_variable = 0.25*y_variable*((double) search_sp);
//
//    /* To use "small" gaps the query and subject must be "large" compared to
//       the gap size. If small gaps may be used, then the cutoff values must be
//       adjusted for the "bayesian" possibility that both large and small gaps
//       are being checked for. */
//
//    if (search_sp > 8*window_size*window_size) {
//       x_variable /= (1.0 - gap_prob + kEpsilon);
//       link_hsp_params->cutoff_big_gap =
//          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
//       x_variable = y_variable*(window_size*window_size);
//       x_variable /= (gap_prob + kEpsilon);
//       link_hsp_params->cutoff_small_gap =
//          MAX(word_params->cutoff_score_min,
//              (Int4) floor((log(x_variable)/kbp->Lambda)) + 1);
//    } else {
//       link_hsp_params->cutoff_big_gap =
//          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
//       /* The following is equivalent to forcing small gap rule to be ignored
//          when linking HSPs. */
//       link_hsp_params->gap_prob = 0;
//       link_hsp_params->cutoff_small_gap = 0;
//    }
//
//    link_hsp_params->cutoff_big_gap *= (Int4)sbp->scale_factor;
//    link_hsp_params->cutoff_small_gap *= (Int4)sbp->scale_factor;
// }
//
// /* For debugging within the C core -RMH- */
// void printBlastInitialWordParamters ( BlastInitialWordParameters *word_params, BlastQueryInfo *query_info )
// {
//   int context;
//   printf("BlastInitialWordParamters:\n");
//   printf("  x_dropoff_max = %d\n", word_params->x_dropoff_max );
//   printf("  cutoff_score_min = %d\n", word_params->cutoff_score_min );
//   printf("  cutoffs:\n");
//   for (context = query_info->first_context;
//        context <= query_info->last_context; ++context)
//   {
//     if (!(query_info->contexts[context].is_valid))
//       continue;
//     printf("    %d x_dropoff_init = %d\n", context, word_params->cutoffs[context].x_dropoff_init );
// ```
#[cfg(test)]
mod cutoff_boundary_tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
    // ```c++
    //          searches. */
    //       if (hit_params->link_hsp_params && !kNucleotide &&
    //           !gapped_calculation) {
    //           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
    //             hit_params->link_hsp_params, word_params, db_length,
    //             seq_arg.seq->length);
    //       }
    //
    //       if (Blast_SubjectIsTranslated(program_number)) {
    // ```
    #[test]
    fn pinned_official_subject_aa_cutoff_boundaries() {
        let golden = include_str!("../../../tests/unit/blastx_even_gap_cutoffs_expected.tsv");
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../docs/evidence/losatx_stage_a/inputs/generated");
        let mut options = super::super::stage_c_tests::options();
        options.gapped = false;
        for line in golden.lines() {
            let f: Vec<_> = line.split('\t').collect();
            let gap = f[0];
            assert_eq!(f[1], "B_LINK_SETUP");
            let records = super::super::input::read_fasta(
                &root.join(format!("S08_gap{gap}.fna")),
                false,
                false,
            )
            .unwrap();
            let batch = super::super::query_setup::prepare_queries(&records, &options).unwrap();
            let mut parameters = super::super::parameters::score_block(&batch);
            let length: i32 = f[2].parse().unwrap();
            let db: i64 = f[3].parse().unwrap();
            super::super::parameters::subject_parameters(
                &batch,
                &mut parameters,
                &options,
                db,
                1,
                length as usize,
                None,
            )
            .unwrap();
            let (scores, prob) = cutoffs(&batch, &parameters, length, db);
            assert_eq!(
                scores,
                [f[4].parse().unwrap(), f[5].parse().unwrap()],
                "{line}"
            );
            assert_eq!(
                prob.to_bits(),
                u64::from_str_radix(f[6], 16).unwrap(),
                "{line}"
            );
        }
    }
}
