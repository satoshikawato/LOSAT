/* Comparison-only NCBI GetGappedScore probe. Never linked by LOSAT.
 * Pinned c++/src/algo/blast/core/blast_engine.c:522-529 calls
 *   aux_struct->GetGappedScore(..., init_hitlist, &hsp_list, ...);
 * Pinned c++/src/algo/blast/core/blast_gapalign.c:3827-3830 requires
 *   Blast_InitHitListIsSortedByScore(init_hitlist);
 * Pinned c++/include/algo/blast/core/blast_hits.h:96-160 defines the
 *   BlastSeg, BlastHSP, and BlastHSPList fields read below.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { int32_t q_start, s_start, length, score; } Ungapped;
typedef struct { uint32_t q_off, s_off; Ungapped *ungapped_data; } InitHsp;
typedef struct { int32_t total, allocated; InitHsp *init_hsp_array; int32_t do_not_reallocate; } InitHitList;
typedef struct { int16_t frame; int32_t offset, end, gapped_start; } Seg;
typedef struct {
    int32_t score, num_ident;
    double bit_score, evalue;
    Seg query, subject;
    int32_t context;
    void *gap_info;
} Hsp;
typedef struct {
    int32_t oid, query_index;
    Hsp **hsp_array;
    int32_t hspcnt, allocated, hsp_max, do_not_reallocate;
    double best_evalue;
} HspList;
/* Pinned c++/include/algo/blast/core/blast_query_info.h:60-99:
 *   query_offset, query_length, eff_searchsp, length_adjustment,
 *   query_index, frame, is_valid, segment_flags;
 *   first_context, last_context, num_queries, contexts.
 * Pinned blast_query_info.c:68-96 assigns one protein context per query.
 */
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

/* Pinned c++/src/algo/blast/core/blast_engine.c:840-850:
 *   Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
 * Pinned blast_hits.c:2809-2864 appends the new list, score-sorts the
 * combined list, and enforces hsp_num_max. Log both input lists and the
 * returned list so per-frame merging and any cap are observable.
 */
static void trace_append_list(const char *event, unsigned long call,
                              const HspList *list)
{
    if (!list) return;
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        const Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "%s\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                event, call, i, hit->context, hit->score, hit->subject.frame,
                hit->query.offset, hit->query.end,
                hit->subject.offset, hit->subject.end);
    }
}

typedef short (*append_fn)(HspList **, HspList **, int32_t);
short Blast_HSPListAppend(HspList **incoming, HspList **combined,
                          int32_t hsp_num_max)
{
    static unsigned long call_index = 0;
    append_fn real = (append_fn)dlsym(RTLD_NEXT, "Blast_HSPListAppend");
    if (!real) return -1;
    unsigned long call = call_index++;
    fprintf(stderr, "APPEND_INPUT\t%lu\t%d\t%d\t%d\n", call,
            hsp_num_max, *incoming ? (*incoming)->hspcnt : 0,
            *combined ? (*combined)->hspcnt : 0);
    trace_append_list("APPEND_IN_HSP", call, *incoming);
    trace_append_list("APPEND_OLD_HSP", call, *combined);
    short status = real(incoming, combined, hsp_num_max);
    fprintf(stderr, "APPEND_OUTPUT\t%lu\t%d\t%d\t%d\n", call, status,
            *incoming ? (*incoming)->hspcnt : 0,
            *combined ? (*combined)->hspcnt : 0);
    trace_append_list("APPEND_OUT_HSP", call, *combined);
    return status;
}
/* Pinned c++/src/algo/blast/core/blast_engine.c:572-586:
 *   Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
 *   Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, kHspNumMax,
 *                       &(backup.offset), INT4_MIN, overlap, ...);
 * Pinned blast_hits.c:2857-3035 merges intersecting overlap HSPs before
 * the frame-level Blast_HSPListAppend. Capture adjusted input and output.
 */
static void trace_merge_list(const char *event, unsigned long call,
                             const HspList *list)
{
    if (!list) return;
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        const Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "%s\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                event, call, i, hit->context, hit->score, hit->subject.frame,
                hit->query.offset, hit->query.end, hit->query.gapped_start,
                hit->subject.offset, hit->subject.end, hit->subject.gapped_start);
    }
}

typedef short (*merge_fn)(HspList **, HspList **, int32_t, int32_t *,
                          int32_t, int32_t, uint8_t, uint8_t);
short Blast_HSPListsMerge(HspList **incoming, HspList **combined,
                          int32_t hsp_num_max, int32_t *split_offsets,
                          int32_t contexts_per_query, int32_t overlap,
                          uint8_t allow_gap, uint8_t short_reads)
{
    static unsigned long call_index = 0;
    merge_fn real = (merge_fn)dlsym(RTLD_NEXT, "Blast_HSPListsMerge");
    if (!real) return -1;
    unsigned long call = call_index++;
    fprintf(stderr, "MERGE_INPUT\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            call, hsp_num_max, split_offsets ? split_offsets[0] : -1,
            contexts_per_query, overlap, allow_gap, short_reads,
            *incoming ? (*incoming)->hspcnt : 0,
            *combined ? (*combined)->hspcnt : 0);
    trace_merge_list("MERGE_IN_HSP", call, *incoming);
    trace_merge_list("MERGE_OLD_HSP", call, *combined);
    short status = real(incoming, combined, hsp_num_max, split_offsets,
                        contexts_per_query, overlap, allow_gap, short_reads);
    fprintf(stderr, "MERGE_OUTPUT\t%lu\t%d\t%d\t%d\n", call, status,
            *incoming ? (*incoming)->hspcnt : 0,
            *combined ? (*combined)->hspcnt : 0);
    trace_merge_list("MERGE_OUT_HSP", call, *combined);
    return status;
}

/* Pinned c++/include/algo/blast/core/blast_def.h:311-319 and
 * c++/include/algo/blast/core/ncbi_std.h:94:
 *   Boolean partial; Int4 num_frames; Int4* range;
 *   typedef Uint1 Boolean;
 * Pinned c++/src/algo/blast/core/blast_hits.c:1147-1229:
 *   start = target_t->range[2*context];
 *   stop = target_t->range[2*context+1];
 *   return target_t->translations[context] - target_t->range[2*context] + 1;
 * Capture each HSP's actual translation window in the comparison oracle.
 */
typedef struct {
    int32_t program_number;
    const uint8_t *gen_code_string;
    uint8_t **translations;
    uint8_t partial;
    int32_t num_frames;
    int32_t *range;
    void *subject_blk;
} TargetTranslation;
typedef const uint8_t *(*target_translation_fn)(TargetTranslation *, const Hsp *, int32_t *);
const uint8_t *Blast_HSPGetTargetTranslation(TargetTranslation *target,
                                             const Hsp *hsp,
                                             int32_t *translated_length)
{
    target_translation_fn real = (target_translation_fn)dlsym(
        RTLD_NEXT, "Blast_HSPGetTargetTranslation");
    if (!real) return NULL;
    const uint8_t *result = real(target, hsp, translated_length);
    int context = hsp->subject.frame > 0 ? hsp->subject.frame - 1 :
                  2 - hsp->subject.frame;
    fprintf(stderr, "TARGET_TRANSLATION\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            hsp->subject.frame, hsp->subject.offset, hsp->subject.end,
            target->range[2 * context], target->range[2 * context + 1],
            translated_length ? *translated_length : -1, target->partial);
    return result;
}

/* Pinned c++/include/algo/blast/core/blast_parameters.h:129-136,153-209. */
typedef struct { void *options; int32_t gap_x_dropoff, gap_x_dropoff_final; } ExtParams;
typedef struct { int32_t cutoff_score, cutoff_score_max; } GappedCutoffs;
typedef struct {
    void *options;
    int32_t cutoff_score_min;
    GappedCutoffs *cutoffs;
    void *link_hsp_params;
    uint8_t restricted_align, do_sum_stats;
    int32_t mask_level;
    int32_t *low_score;
    double prelim_evalue;
} HitParams;
typedef struct {
    void *options;
    int16_t reward, penalty;
    int32_t gap_open, gap_extend, shift_pen;
    double scale_factor;
} ScoreParams;
typedef short (*gapped_fn)(int, void *, void *, void *, void *, void *,
                           void *, void *, void *, InitHitList *,
                           HspList **, void *, int *);

short BLAST_GetGappedScore(int program, void *query, void *query_info,
                          void *subject, void *gap_align, void *score_params,
                          void *ext_params, void *hit_params,
                          void *word_params, InitHitList *init_hitlist,
                          HspList **hsp_list, void *gapped_stats,
                          int *fence_hit)
{
    static unsigned long call_index = 0;
    gapped_fn real = (gapped_fn)dlsym(RTLD_NEXT, "BLAST_GetGappedScore");
    if (!real) return -1;
    unsigned long call = call_index++;
    QueryInfo *qinfo = (QueryInfo *)query_info;
    fprintf(stderr, "QUERY_INFO\t%lu\t%d\t%d\t%d\n", call,
            qinfo->first_context, qinfo->last_context, qinfo->num_queries);
    for (int32_t context = qinfo->first_context;
         context <= qinfo->last_context; ++context) {
        QueryContext *q = &qinfo->contexts[context];
        fprintf(stderr, "QUERY_CONTEXT\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\n",
                call, context, q->query_index, q->query_offset,
                q->query_length, q->frame, q->is_valid);
    }
    ExtParams *ext = (ExtParams *)ext_params;
    HitParams *hit = (HitParams *)hit_params;
    ScoreParams *scoring = (ScoreParams *)score_params;
    fprintf(stderr, "GAPPED_PARAMS\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            call, scoring->gap_open, scoring->gap_extend,
            ext->gap_x_dropoff, ext->gap_x_dropoff_final,
            hit->cutoff_score_min, hit->cutoffs[0].cutoff_score,
            hit->cutoffs[0].cutoff_score_max, hit->restricted_align,
            hit->do_sum_stats);
    fprintf(stderr, "GAPPED_INPUT\t%lu\t%d\n", call, init_hitlist->total);
    for (int32_t i = 0; i < init_hitlist->total; ++i) {
        InitHsp *hit = &init_hitlist->init_hsp_array[i];
        Ungapped *data = hit->ungapped_data;
        fprintf(stderr, "GAPPED_INIT\t%lu\t%d\t%u\t%u\t%d\t%d\t%d\t%d\n",
                call, i, hit->q_off, hit->s_off,
                data ? data->q_start : -1, data ? data->s_start : -1,
                data ? data->length : -1, data ? data->score : -1);
    }
    short status = real(program, query, query_info, subject, gap_align,
                        score_params, ext_params, hit_params, word_params,
                        init_hitlist, hsp_list, gapped_stats, fence_hit);
    int32_t count = *hsp_list ? (*hsp_list)->hspcnt : 0;
    fprintf(stderr, "GAPPED_OUTPUT\t%lu\t%d\t%d\n", call, status, count);
    for (int32_t i = 0; i < count; ++i) {
        Hsp *hit = (*hsp_list)->hsp_array[i];
        fprintf(stderr, "GAPPED_HSP\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                call, i, hit->score, hit->context,
                hit->query.frame, hit->query.offset, hit->query.end,
                hit->query.gapped_start, hit->subject.frame,
                hit->subject.offset, hit->subject.end,
                hit->subject.gapped_start);
    }
    return status;
}

/* Pinned c++/src/algo/blast/core/blast_traceback.c:635-666:
 *   purge common endpoints; then for each surviving HSP
 *   delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
 *   if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp);
 * Pinned c++/src/algo/blast/core/blast_hits.c:479-660 updates raw score,
 * offsets, edit script, and returns the deletion decision.
 */
typedef uint8_t (*reevaluate_fn)(Hsp *, const uint8_t *, int32_t,
                                 const uint8_t *, int32_t, const HitParams *,
                                 const ScoreParams *, const void *);
uint8_t Blast_HSPReevaluateWithAmbiguitiesGapped(
        Hsp *hsp, const uint8_t *query, int32_t query_length,
        const uint8_t *subject, int32_t subject_length,
        const HitParams *hit_params, const ScoreParams *score_params,
        const void *score_block)
{
    static unsigned long call_index = 0;
    reevaluate_fn real = (reevaluate_fn)dlsym(
        RTLD_NEXT, "Blast_HSPReevaluateWithAmbiguitiesGapped");
    if (!real) return 1;
    unsigned long call = call_index++;
    fprintf(stderr, "REEVAL_INPUT\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            call, hsp->score, hsp->subject.frame, hsp->query.offset,
            hsp->query.end, hsp->subject.offset, hsp->subject.end,
            query_length, subject_length,
            hit_params->cutoffs[hsp->context].cutoff_score);
    uint8_t deleted = real(hsp, query, query_length, subject, subject_length,
                           hit_params, score_params, score_block);
    fprintf(stderr, "REEVAL_OUTPUT\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            call, deleted, hsp->score, hsp->query.offset, hsp->query.end,
            hsp->subject.offset, hsp->subject.end, hsp->subject.frame);
    return deleted;
}

/* Pinned c++/src/algo/blast/core/blast_traceback.c:259-288,383-435,
 * 503-589,613-721:
 *   Blast_TracebackFromHSPList(..., hsp_list, ...);
 *   BLAST_GappedAlignmentWithTraceback(...);
 *   Blast_HSPUpdateWithTraceback(gap_align, hsp);
 *   Blast_HSPListSortByScore(hsp_list);
 * The list is traced before and after the actual NCBI call to locate changes
 * caused by target translation, traceback, reevaluation, and deletion.
 */
typedef short (*traceback_fn)(int, HspList *, void *, void *, void *,
                              void *, void *, void *, void *, void *,
                              const uint8_t *, uint8_t *);
short Blast_TracebackFromHSPList(
        int program, HspList *list, void *query, void *subject,
        void *query_info, void *gap_align, void *score_block,
        void *score_params, void *ext_options, void *hit_params,
        const uint8_t *gen_code, uint8_t *fence_hit)
{
    traceback_fn real = (traceback_fn)dlsym(RTLD_NEXT, "Blast_TracebackFromHSPList");
    if (!real) return -1;
    fprintf(stderr, "TRACEBACK_INPUT\t%d\t%d\n", list->oid, list->hspcnt);
    /* Pinned blast_traceback.c:1242-1261,1644-1684 passes the query-indexed
     * HSP list through traceback and retries the same list after a fence. */
    fprintf(stderr, "TRACEBACK_CONTEXT\t%d\t%d\n",
            list->query_index, list->hspcnt);
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "TRACEBACK_IN_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                list->oid, i, hit->score, hit->subject.frame,
                hit->query.offset, hit->query.end,
                hit->subject.offset, hit->subject.end);
        fprintf(stderr, "TRACEBACK_IN_CONTEXT_HSP\t%d\t%d\t%d\n",
                list->query_index, i, hit->context);
    }
    short status = real(program, list, query, subject, query_info, gap_align,
                        score_block, score_params, ext_options, hit_params,
                        gen_code, fence_hit);
    fprintf(stderr, "TRACEBACK_OUTPUT\t%d\t%d\t%d\t%d\n",
            list->oid, status, list->hspcnt, (int)*fence_hit);
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "TRACEBACK_OUT_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                list->oid, i, hit->score, hit->subject.frame,
                hit->query.offset, hit->query.end,
                hit->subject.offset, hit->subject.end);
    }
    return status;
}

/* Pinned c++/src/algo/blast/core/blast_traceback.c:401-409,585-603,
 * 677-696 calls containment before traceback, HSPTest after alignment,
 * and containment again after the score sort. Pinned blast_itree.c:931-953
 * defines the containment predicate; blast_hits.h:356-359 defines HSPTest.
 * Record the exact deletion predicate for the six-frame comparison case.
 */
typedef uint8_t (*contains_fn)(const void *, const Hsp *, const void *, int32_t);
uint8_t BlastIntervalTreeContainsHSP(const void *tree, const Hsp *hsp,
                                     const void *query_info, int32_t separation)
{
    contains_fn real = (contains_fn)dlsym(RTLD_NEXT, "BlastIntervalTreeContainsHSP");
    if (!real) return 0;
    uint8_t result = real(tree, hsp, query_info, separation);
    if (hsp->score >= 0)
        fprintf(stderr, "CONTAINS\t%d\t%d\t%d\t%d\t%d\t%d\n",
                result, hsp->subject.frame, hsp->query.offset,
                hsp->query.end, hsp->subject.offset, hsp->subject.end);
    return result;
}
typedef uint8_t (*hsp_test_fn)(Hsp *, const void *, int32_t);
uint8_t Blast_HSPTest(Hsp *hsp, const void *hit_options, int32_t align_length)
{
    hsp_test_fn real = (hsp_test_fn)dlsym(RTLD_NEXT, "Blast_HSPTest");
    if (!real) return 0;
    uint8_t result = real(hsp, hit_options, align_length);
    if (hsp->score >= 0)
        fprintf(stderr, "HSP_TEST\t%d\t%d\t%d\t%d\t%d\t%d\n",
                result, align_length, hsp->query.offset, hsp->query.end,
                hsp->subject.offset, hsp->subject.end);
    return result;
}

/* Pinned c++/src/algo/blast/core/blast_traceback.c:635-666:
 *   extra_start = Blast_HSPListPurgeHSPsWithCommonEndpoints(
 *       program_number, hsp_list, FALSE);
 * Pinned c++/src/algo/blast/core/blast_hits.c:2455-2537:
 *   purge |= (program != eBlastTypeBlastn);
 *   qsort(..., s_QueryOffsetCompareHSPs); ...
 *   qsort(..., s_QueryEndCompareHSPs);
 * Record list inputs, outputs, and returned count without changing NCBI.
 */
typedef int32_t (*purge_endpoints_fn)(int32_t, HspList *, uint8_t);
int32_t Blast_HSPListPurgeHSPsWithCommonEndpoints(int32_t program,
                                                   HspList *list,
                                                   uint8_t purge)
{
    purge_endpoints_fn real = (purge_endpoints_fn)dlsym(
        RTLD_NEXT, "Blast_HSPListPurgeHSPsWithCommonEndpoints");
    if (!real) return -1;
    fprintf(stderr, "ENDPOINT_INPUT\t%d\t%d\n", purge, list->hspcnt);
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "ENDPOINT_IN_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                i, hit->score, hit->subject.frame, hit->query.offset,
                hit->query.end, hit->subject.offset, hit->subject.end);
    }
    int32_t retained = real(program, list, purge);
    fprintf(stderr, "ENDPOINT_OUTPUT\t%d\t%d\n", retained, list->hspcnt);
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "ENDPOINT_OUT_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                i, hit->score, hit->subject.frame, hit->query.offset,
                hit->query.end, hit->subject.offset, hit->subject.end);
    }
    return retained;
}
