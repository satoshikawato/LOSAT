/* Comparison-only extension of the pinned gapped probe. The original source
 * and its retained fixture hashes remain unchanged.
 * NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
 * cutoff = hit_params->cutoffs[context].cutoff_score;
 */
#define BLAST_GetGappedScore BLAST_GetGappedScore_base
#include "ncbi_gapped_trace.c"
#undef BLAST_GetGappedScore

short BLAST_GetGappedScore(int program, void *query, void *query_info,
                          void *subject, void *gap_align, void *score_params,
                          void *ext_params, void *hit_params,
                          void *word_params, InitHitList *init_hitlist,
                          HspList **hsp_list, void *gapped_stats,
                          int *fence_hit)
{
    static unsigned long call_index = 0;
    QueryInfo *qinfo = (QueryInfo *)query_info;
    HitParams *hit = (HitParams *)hit_params;
    unsigned long call = call_index++;
    for (int32_t context = qinfo->first_context;
         context <= qinfo->last_context; ++context) {
        fprintf(stderr, "GAPPED_CONTEXT_CUTOFF\t%lu\t%d\t%d\t%d\t%d\n",
                call, context, qinfo->contexts[context].query_index,
                hit->cutoffs[context].cutoff_score,
                hit->cutoffs[context].cutoff_score_max);
    }
    return BLAST_GetGappedScore_base(program, query, query_info,
                                    subject, gap_align, score_params,
                                    ext_params, hit_params, word_params,
                                    init_hitlist, hsp_list, gapped_stats,
                                    fence_hit);
}
