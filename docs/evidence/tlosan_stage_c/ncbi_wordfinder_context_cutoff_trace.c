/* Comparison-only NCBI per-context WordFinder cutoff probe.
 * Pinned aa_ungapped.c:547-583 selects word_params->cutoffs after
 * BSearchContextInfo(query_offset, query_info), then uses its x-drop and
 * cutoff_score for extension and HSP saving. The base probe logs every
 * scansub pair at aa_ungapped.c:478-505 without changing NCBI output.
 */
#define BlastAaWordFinder BlastAaWordFinder_base
#include "ncbi_candidate_trace.c"
#undef BlastAaWordFinder

typedef struct { int32_t first_context,last_context,num_queries; void *contexts; } QueryInfoPrefix;
short BlastAaWordFinder(void* subject, void* query, void* query_info,
        void* lookup, void* matrix, void* word_params, void* ewp,
        void* offset_pairs, int32_t offset_array_size,
        void* init_hitlist, void* ungapped_stats)
{
    QueryInfoPrefix *qi=(QueryInfoPrefix*)query_info;
    WordParamsPrefix *params=(WordParamsPrefix*)word_params;
    unsigned long call=call_index;
    for (int32_t context=qi->first_context; context<=qi->last_context; ++context) {
        WordCutoffs c=params->cutoffs[context];
        fprintf(stderr,"WORD_CONTEXT_CUTOFF\t%lu\t%d\t%d\t%d\t%d\n",
                call,context,c.x_dropoff_init,c.x_dropoff,c.cutoff_score);
    }
    return BlastAaWordFinder_base(subject,query,query_info,lookup,matrix,
            word_params,ewp,offset_pairs,offset_array_size,init_hitlist,
            ungapped_stats);
}
