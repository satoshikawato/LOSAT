/* Comparison-only probe for pinned NCBI Stage D calls. Never linked by LOSAT.
 * Source: c++/src/algo/blast/core/blast_engine.c:870-905;
 * blast_traceback.c:1481-1499; blast_kappa.c:407-441,2981-2996;
 * link_hsps.c:1765-1810; blast_hits.c:1811-1926.
 * Layout: c++/include/algo/blast/core/blast_hits.h:111-165,
 * blast_query_info.h:60-92, blast_parameters.h:137-150,168-191.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { int16_t frame; int32_t offset, end, gapped_start; } Seg;
typedef struct {
    int32_t score, num_ident;
    double bit_score, evalue;
    Seg query, subject;
    int32_t context;
    void *gap_info;
    int32_t num;
    int16_t comp_adjustment_method;
} Hsp;
typedef struct {
    int32_t oid, query_index;
    Hsp **hsp_array;
    int32_t hspcnt, allocated, hsp_max;
    uint8_t do_not_reallocate;
    double best_evalue;
} HspList;
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
typedef struct {
    double gap_prob;
    int32_t gap_size, overlap_size;
    double gap_decay_rate;
    int32_t cutoff_small_gap, cutoff_big_gap, longest_intron;
} LinkParams;
typedef struct { int32_t cutoff_score, cutoff_score_max; } GappedCutoffs;
typedef struct {
    void *options;
    int32_t cutoff_score_min;
    GappedCutoffs *cutoffs;
    LinkParams *link_hsp_params;
    uint8_t restricted_align, do_sum_stats;
    int32_t mask_level;
    int32_t *low_score;
    double prelim_evalue;
} HitParams;

static unsigned long event_number;
/* Pinned blast_setup.c:1006-1024 calls CalcEffLengths, hit-cutoff update,
 * then word-cutoff and link update only when word_params is non-NULL.
 * blast_engine.c:1434-1442 calls this before local-subject search;
 * blast_traceback.c:1600-1630 calls it in the ordinary traceback branch.
 */
typedef short (*update_fn)(int32_t, uint32_t, const void *, QueryInfo *,
                           const void *, HitParams *, void *, void *);
short BLAST_OneSubjectUpdateParameters(int32_t program, uint32_t subject_length,
                                       const void *scoring, QueryInfo *query,
                                       const void *sbp, HitParams *hit,
                                       void *word, void *eff)
{
    update_fn real = (update_fn)dlsym(RTLD_NEXT,
                                      "BLAST_OneSubjectUpdateParameters");
    unsigned long n = event_number++;
    fprintf(stderr, "D_CALL\t%lu\tupdate\t%d\t%u\t%d\t%d\n",
            n, program, subject_length, word != NULL,
            query ? query->num_queries : -1);
    short status = real(program, subject_length, scoring, query, sbp,
                        hit, word, eff);
    fprintf(stderr, "D_RETURN\t%lu\tupdate\t%d\t%d\t%d\t%.17g\n",
            n, status, hit ? hit->do_sum_stats : -1,
            hit && hit->link_hsp_params ?
                hit->link_hsp_params->longest_intron : -1,
            hit ? hit->prelim_evalue : -1.);
    if (query) {
        for (int32_t i = query->first_context; i <= query->last_context; ++i) {
            const QueryContext *c = &query->contexts[i];
            fprintf(stderr,
                    "D_UPDATE_CONTEXT\t%lu\t%d\t%d\t%d\t%lld\t%d\t%d\t%d\n",
                    n, i, c->query_index, c->query_length,
                    (long long)c->eff_searchsp, c->length_adjustment,
                    hit && hit->cutoffs ? hit->cutoffs[i].cutoff_score : -1,
                    hit && hit->cutoffs ? hit->cutoffs[i].cutoff_score_max : -1);
        }
    }
    return status;
}
static void dump_list(const char *stage, unsigned long event, const HspList *list)
{
    if (!list) {
        fprintf(stderr, "D_LIST\t%lu\t%s\tNULL\n", event, stage);
        return;
    }
    fprintf(stderr, "D_LIST\t%lu\t%s\t%d\t%d\t%d\t%.17g\n",
            event, stage, list->oid, list->query_index,
            list->hspcnt, list->best_evalue);
    for (int32_t i = 0; i < list->hspcnt; ++i) {
        const Hsp *h = list->hsp_array[i];
        fprintf(stderr,
                "D_HSP\t%lu\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%d\n",
                event, stage, i, h->context, h->subject.frame,
                h->query.offset, h->query.end, h->subject.offset,
                h->subject.end, h->score, h->num_ident, h->num,
                h->evalue, h->bit_score, h->comp_adjustment_method,
                h->subject.gapped_start);
    }
}

typedef short (*link_fn)(int32_t, HspList *, const QueryInfo *, int32_t,
                         const void *, const LinkParams *, uint8_t);
short BLAST_LinkHsps(int32_t program, HspList *list, const QueryInfo *query,
                     int32_t subject_length, const void *sbp,
                     const LinkParams *params, uint8_t gapped)
{
    link_fn real = (link_fn)dlsym(RTLD_NEXT, "BLAST_LinkHsps");
    unsigned long n = event_number++;
    fprintf(stderr,
            "D_CALL\t%lu\tlink\t%d\t%d\t%d\t%d\t%.17g\t%d\t%d\t%d\n",
            n, program, subject_length, gapped,
            params ? params->longest_intron : -1,
            params ? params->gap_decay_rate : -1.,
            params ? params->gap_size : -1,
            params ? params->cutoff_small_gap : -1,
            params ? params->cutoff_big_gap : -1);
    if (query && list && list->hspcnt) {
        int context = list->hsp_array[0]->context;
        fprintf(stderr, "D_CONTEXT\t%lu\t%d\t%d\t%lld\t%d\n",
                n, context, query->contexts[context].query_length,
                (long long)query->contexts[context].eff_searchsp,
                query->contexts[context].length_adjustment);
    }
    dump_list("link_before", n, list);
    short status = real(program, list, query, subject_length, sbp, params, gapped);
    fprintf(stderr, "D_RETURN\t%lu\tlink\t%d\n", n, status);
    dump_list("link_after", n, list);
    return status;
}

typedef short (*evalue_fn)(int32_t, const QueryInfo *, int32_t, HspList *,
                           uint8_t, uint8_t, const void *, double, double);
short Blast_HSPListGetEvalues(int32_t program, const QueryInfo *query,
                              int32_t subject_length, HspList *list,
                              uint8_t gapped, uint8_t rps, const void *sbp,
                              double gap_decay, double scaling)
{
    evalue_fn real = (evalue_fn)dlsym(RTLD_NEXT, "Blast_HSPListGetEvalues");
    unsigned long n = event_number++;
    fprintf(stderr, "D_CALL\t%lu\tevalue\t%d\t%d\t%d\t%d\t%.17g\t%.17g\n",
            n, program, subject_length, gapped, rps, gap_decay, scaling);
    if (query && list && list->hspcnt) {
        int context = list->hsp_array[0]->context;
        fprintf(stderr, "D_CONTEXT\t%lu\t%d\t%d\t%lld\t%d\n",
                n, context, query->contexts[context].query_length,
                (long long)query->contexts[context].eff_searchsp,
                query->contexts[context].length_adjustment);
    }
    dump_list("evalue_before", n, list);
    short status = real(program, query, subject_length, list, gapped,
                        rps, sbp, gap_decay, scaling);
    fprintf(stderr, "D_RETURN\t%lu\tevalue\t%d\n", n, status);
    dump_list("evalue_after", n, list);
    return status;
}

typedef short (*reap_fn)(HspList *, const void *);
short Blast_HSPListReapByEvalue(HspList *list, const void *options)
{
    reap_fn real = (reap_fn)dlsym(RTLD_NEXT, "Blast_HSPListReapByEvalue");
    unsigned long n = event_number++;
    fprintf(stderr, "D_CALL\t%lu\treap\n", n);
    dump_list("reap_before", n, list);
    short status = real(list, options);
    fprintf(stderr, "D_RETURN\t%lu\treap\t%d\n", n, status);
    dump_list("reap_after", n, list);
    return status;
}

typedef short (*bit_fn)(HspList *, uint8_t, const void *);
short Blast_HSPListGetBitScores(HspList *list, uint8_t gapped, const void *sbp)
{
    bit_fn real = (bit_fn)dlsym(RTLD_NEXT, "Blast_HSPListGetBitScores");
    unsigned long n = event_number++;
    fprintf(stderr, "D_CALL\t%lu\tbits\t%d\n", n, gapped);
    dump_list("bits_before", n, list);
    short status = real(list, gapped, sbp);
    fprintf(stderr, "D_RETURN\t%lu\tbits\t%d\n", n, status);
    dump_list("bits_after", n, list);
    return status;
}

typedef short (*redo_fn)(int32_t, uint32_t, void *, const QueryInfo *, void *,
                         void *, const void *, int32_t, HspList *, void *,
                         void *, const void *, const void *, const void *, void *);
short Blast_RedoAlignmentCore_MT(int32_t program, uint32_t threads,
                                 void *query_blk, const QueryInfo *query,
                                 void *sbp, void *subject_blk,
                                 const void *seqsrc, int32_t code,
                                 HspList *match, void *stream,
                                 void *scoring, const void *extension,
                                 const void *hit, const void *psi, void *results)
{
    redo_fn real = (redo_fn)dlsym(RTLD_NEXT, "Blast_RedoAlignmentCore_MT");
    unsigned long n = event_number++;
    fprintf(stderr, "D_CALL\t%lu\tredo\t%d\t%u\t%d\t%d\t%d\n",
            n, program, threads, code, seqsrc != NULL,
            query ? query->num_queries : -1);
    dump_list("redo_before", n, match);
    short status = real(program, threads, query_blk, query, sbp, subject_blk,
                        seqsrc, code, match, stream, scoring, extension,
                        hit, psi, results);
    fprintf(stderr, "D_RETURN\t%lu\tredo\t%d\n", n, status);
    dump_list("redo_after", n, match);
    return status;
}
