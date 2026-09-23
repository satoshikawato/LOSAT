/* Comparison-only NCBI function probe. Never linked or called by LOSAT.
 * Pinned NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
 *   status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
 *   Blast_InitHitListSortByScore(init_hitlist);
 * Pinned headers:
 *   c++/include/algo/blast/core/blast_extend.h:142-163:
 *   typedef struct BlastInitHSP { BlastOffsetPair offsets;
 *       BlastUngappedData* ungapped_data; } BlastInitHSP;
 *   typedef struct BlastInitHitList { Int4 total; Int4 allocated;
 *       BlastInitHSP* init_hsp_array; Boolean do_not_reallocate; } ...
 *   c++/include/algo/blast/core/blast_def.h:141-150:
 *   struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { int32_t q_start, s_start, length, score; } Ungapped;
typedef struct { uint32_t q_off, s_off; Ungapped *ungapped_data; } InitHsp;
typedef struct { int32_t total, allocated; InitHsp *init_hsp_array; int32_t do_not_reallocate; } InitHitList;
typedef short (*wordfinder_fn)(void*, void*, void*, void*, void*, void*,
        void*, void*, int32_t, InitHitList*, void*);

short BlastAaWordFinder(void* subject, void* query, void* query_info,
        void* lookup, void* matrix, void* word_params, void* ewp,
        void* offset_pairs, int32_t offset_array_size,
        InitHitList* init_hitlist, void* ungapped_stats)
{
    static unsigned long call_index = 0;
    wordfinder_fn real = (wordfinder_fn)dlsym(RTLD_NEXT, "BlastAaWordFinder");
    if (!real) return -1;
    short status = real(subject, query, query_info, lookup, matrix, word_params,
            ewp, offset_pairs, offset_array_size, init_hitlist, ungapped_stats);
    unsigned long call = call_index++;
    fprintf(stderr, "WORD_FINDER\t%lu\t%d\n", call, init_hitlist->total);
    for (int32_t i = 0; i < init_hitlist->total; ++i) {
        InitHsp* hit = &init_hitlist->init_hsp_array[i];
        Ungapped* data = hit->ungapped_data;
        fprintf(stderr, "INIT\t%lu\t%d\t%u\t%u\t%d\t%d\t%d\t%d\n",
                call, i, hit->q_off, hit->s_off,
                data ? data->q_start : -1, data ? data->s_start : -1,
                data ? data->length : -1, data ? data->score : -1);
    }
    return status;
}
