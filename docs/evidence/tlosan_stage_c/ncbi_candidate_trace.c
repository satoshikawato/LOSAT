/* Comparison-only NCBI scan-callback probe; never link this into LOSAT.
 * Pinned c++/src/algo/blast/core/aa_ungapped.c:478-505:
 *   scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
 *   hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
 * Pinned c++/include/algo/blast/core/lookup_wrap.h:50-58 and
 * blast_aalookup.h:99-132 describe the callback holder layout.
 * Pinned c++/include/algo/blast/core/blast_aascan.h:45-49 gives its ABI.
 * Pinned c++/include/algo/blast/core/blast_parameters.h:95-114
 * gives the per-query-context initial word cutoff fields.
 * Pinned c++/include/algo/blast/core/blast_def.h:141-150 gives pair layout.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { uint32_t q_off, s_off; } OffsetPair;
typedef struct { int lut_type; void* lut; } LookupWrapPrefix;
typedef struct {
    int32_t threshold, mask, charsize, word_length, lut_word_length;
    int32_t alphabet_size, backbone_size, longest_chain;
    int32_t **thin_backbone;
    int32_t bone_type;
    void *thick_backbone, *overflow;
    int32_t overflow_size;
    void *pv;
    int32_t use_pssm;
    void *scansub_callback;
} AaLookupPrefix;
typedef struct {
    int32_t x_dropoff_init, x_dropoff, cutoff_score, reduced_nucl_cutoff_score;
} WordCutoffs;
typedef struct {
    void* options;
    int32_t x_dropoff_max, cutoff_score_min;
    WordCutoffs* cutoffs;
} WordParamsPrefix;
typedef int32_t (*scan_fn)(const void*, const void*, OffsetPair*, int32_t, int32_t*);
typedef short (*wordfinder_fn)(void*, void*, void*, void*, void*, void*,
        void*, void*, int32_t, void*, void*);

static scan_fn original_scan;
static unsigned long current_call;
static unsigned long call_index;
static unsigned long ordinal;

static int32_t trace_scan(const void* lookup, const void* subject,
        OffsetPair* pairs, int32_t capacity, int32_t* range)
{
    int32_t count = original_scan(lookup, subject, pairs, capacity, range);
    for (int32_t i = 0; i < count; ++i) {
        fprintf(stderr, "CAND\t%lu\t%lu\t%u\t%u\n",
                current_call, ordinal++, pairs[i].q_off, pairs[i].s_off);
    }
    return count;
}

short BlastAaWordFinder(void* subject, void* query, void* query_info,
        void* lookup, void* matrix, void* word_params, void* ewp,
        void* offset_pairs, int32_t offset_array_size,
        void* init_hitlist, void* ungapped_stats)
{
    wordfinder_fn real = (wordfinder_fn)dlsym(RTLD_NEXT, "BlastAaWordFinder");
    if (!real) return -1;
    LookupWrapPrefix* wrap = (LookupWrapPrefix*)lookup;
    if (wrap->lut_type != 3) return -1; /* eAaLookupTable */
    AaLookupPrefix* aa = (AaLookupPrefix*)wrap->lut;
    original_scan = (scan_fn)aa->scansub_callback;
    current_call = call_index++;
    ordinal = 0;
    WordCutoffs* cutoffs = ((WordParamsPrefix*)word_params)->cutoffs;
    fprintf(stderr, "PARAM\t%lu\t%d\t%d\t%d\n", current_call,
            cutoffs[0].x_dropoff_init, cutoffs[0].x_dropoff,
            cutoffs[0].cutoff_score);
    aa->scansub_callback = (void*)trace_scan;
    short status = real(subject, query, query_info, lookup, matrix, word_params,
            ewp, offset_pairs, offset_array_size, init_hitlist, ungapped_stats);
    aa->scansub_callback = (void*)original_scan;
    fprintf(stderr, "END\t%lu\t%lu\n", current_call, ordinal);
    return status;
}
