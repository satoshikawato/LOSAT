/* Comparison-only pinned NCBI query state probe. Never linked by LOSAT.
 * c++/src/algo/blast/core/blast_setup.c:614-625 masks the working query
 * when mask_at_hash is false. blast_engine.c:484-525 passes the same
 * BLAST_SequenceBlk to WordFinder and GetGappedScore. blast_gapalign.c:
 * 2410-2442 points single_query->sequence into that same working block.
 * BLAST_SequenceBlk prefix is in c++/include/algo/blast/core/blast_def.h:242-246.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { uint8_t *sequence, *sequence_start; int32_t length; } QueryPrefix;
static void trace_query(const char *label, unsigned long call, void *query)
{
    const QueryPrefix *q = (const QueryPrefix *)query;
    fprintf(stderr, "%s\t%lu\t%d\t", label, call, q->length);
    for (int32_t i = 0; i < q->length; ++i) fprintf(stderr, "%02x", q->sequence[i]);
    fputc('\n', stderr);
}
typedef short (*word_fn)(void*,void*,void*,void*,void*,void*,void*,void*,int32_t,void*,void*);
short BlastAaWordFinder(void *subject, void *query, void *query_info,
    void *lookup, void *matrix, void *word_params, void *ewp,
    void *offset_pairs, int32_t offset_array_size, void *init_hitlist,
    void *ungapped_stats)
{
    static unsigned long call = 0;
    word_fn real = (word_fn)dlsym(RTLD_NEXT, "BlastAaWordFinder");
    if (!real) return -1;
    trace_query("WORD_QUERY_BYTES", call++, query);
    return real(subject, query, query_info, lookup, matrix, word_params,
        ewp, offset_pairs, offset_array_size, init_hitlist, ungapped_stats);
}
typedef short (*gap_fn)(int,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*);
short BLAST_GetGappedScore(int program, void *query, void *query_info,
    void *subject, void *gap_align, void *score_params, void *ext_params,
    void *hit_params, void *word_params, void *init_hitlist,
    void *hsp_list, void *gapped_stats, void *fence_hit)
{
    static unsigned long call = 0;
    gap_fn real = (gap_fn)dlsym(RTLD_NEXT, "BLAST_GetGappedScore");
    if (!real) return -1;
    trace_query("GAPPED_QUERY_BYTES", call++, query);
    return real(program, query, query_info, subject, gap_align, score_params,
        ext_params, hit_params, word_params, init_hitlist, hsp_list,
        gapped_stats, fence_hit);
}
/* Pinned blast_traceback.c:583-596 passes query_nomask here; blast_hits.h:
 * 339-345 and blast_hits.c:966-991 define the identity result fields. */
typedef struct { int16_t frame; int32_t offset, end, gapped_start; } Segment;
typedef struct {
    int32_t score, num_ident;
    double bit_score, evalue;
    Segment query, subject;
    int32_t context;
    void *gap_info;
} Hsp;
typedef short (*identity_fn)(const uint8_t*,const uint8_t*,Hsp*,const void*,int32_t*,const void*);
short Blast_HSPGetNumIdentitiesAndPositives(const uint8_t *query,
    const uint8_t *subject, Hsp *hsp, const void *score_options,
    int32_t *align_length, const void *sbp)
{
    static unsigned long call = 0;
    identity_fn real = (identity_fn)dlsym(RTLD_NEXT,
        "Blast_HSPGetNumIdentitiesAndPositives");
    if (!real) return -1;
    short status = real(query, subject, hsp, score_options, align_length, sbp);
    fprintf(stderr, "IDENTITY_STATE\t%lu\t%d\t%d\t%d\t%d\t%d\t",
        call++, hsp->score, hsp->num_ident, *align_length,
        hsp->query.offset, hsp->query.end);
    for (int32_t i = 0; i < hsp->query.end; ++i)
        fprintf(stderr, "%02x", query[i]);
    fputc('\n', stderr);
    return status;
}
