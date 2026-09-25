/* Comparison-only NCBI score-setup probe. Never linked into LOSAT.
 * Pinned c++/src/algo/blast/core/blast_setup.c:456-464,652-654 calls
 * BlastSetup_ScoreBlkInit(query_blk, query_info, scoring_options, ...).
 * Context layout: c++/include/algo/blast/core/blast_query_info.h:60-92.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    int32_t query_offset, query_length;
    int64_t eff_searchsp;
    int32_t length_adjustment, query_index;
    int8_t frame;
    uint8_t is_valid;
    int32_t segment_flags;
} QueryContext;
typedef struct {
    int32_t first_context, last_context, num_queries;
    QueryContext *contexts;
    uint32_t max_length, min_length;
    void *pattern_info;
} QueryInfo;
typedef short (*real_fn)(void *, const QueryInfo *, const void *, int,
                         void **, double, void **, void *);
short BlastSetup_ScoreBlkInit(void *query_blk, const QueryInfo *query_info,
                             const void *scoring_options, int program_number,
                             void **sbpp, double scale_factor,
                             void **blast_message, void *get_path) {
    real_fn real = (real_fn)dlsym(RTLD_NEXT, "BlastSetup_ScoreBlkInit");
    short status = real(query_blk, query_info, scoring_options, program_number,
                        sbpp, scale_factor, blast_message, get_path);
    fprintf(stderr, "G_CONTEXT status=%d program=%d first=%d last=%d\n",
            status, program_number, query_info->first_context,
            query_info->last_context);
    for (int i = query_info->first_context; i <= query_info->last_context; ++i)
        fprintf(stderr, "G_CONTEXT index=%d valid=%d length=%d space=%lld\n",
                i, query_info->contexts[i].is_valid,
                query_info->contexts[i].query_length,
                (long long)query_info->contexts[i].eff_searchsp);
    return status;
}
