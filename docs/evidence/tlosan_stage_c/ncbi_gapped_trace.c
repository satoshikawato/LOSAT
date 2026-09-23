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
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        Hsp *hit = list->hsp_array[i];
        fprintf(stderr, "TRACEBACK_IN_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                list->oid, i, hit->score, hit->subject.frame,
                hit->query.offset, hit->query.end,
                hit->subject.offset, hit->subject.end);
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
