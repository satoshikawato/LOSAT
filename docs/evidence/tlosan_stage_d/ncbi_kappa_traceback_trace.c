/* Comparison-only TBLASTN Kappa traceback probe. Never linked into LOSAT.
 * Pinned NCBI: composition_adjustment/redo_alignment.c:1101-1111;
 * core/blast_kappa.c:1896-1957;
 * core/blast_gapalign.h:70-100,167-173.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    int32_t *op_type;
    int32_t *num;
    int32_t size;
} EditScript;

typedef struct {
    uint8_t position_based;
    void *state_struct, *edit_script, *fwd_prelim_tback, *rev_prelim_tback;
    void *greedy_align_mem, *dp_mem;
    int32_t dp_mem_alloc;
    void *sbp;
    int32_t gap_x_dropoff, max_mismatches, mismatch_window;
    int32_t query_start, query_stop, subject_start, subject_stop;
    int32_t greedy_query_seed_start, greedy_subject_seed_start, score;
} GapAlign;

typedef struct {
    void *options;
    int16_t reward, penalty;
    int32_t gap_open, gap_extend, shift_pen;
    double scale_factor;
} Scoring;

typedef struct {
    void *matrix_info, *gapping_params;
    int32_t compo_mode, position_based, pseudocounts;
    int32_t subject_is_translated, query_is_translated;
    int32_t ccat_query_length, cutoff_s;
    double cutoff_e;
    int32_t do_link_hsps;
    void *callbacks;
    double near_identical_cutoff;
} RedoParams;
typedef struct { int32_t length, index; void *local_data; } MatchingSeq;
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
typedef struct Alignment {
    int32_t score, rule, query_index, query_start, query_end;
    int32_t match_start, match_end, frame;
    Hsp *context;
    struct Alignment *next;
} Alignment;
/* NCBI comparison ABI: core/blast_hits.h:153-166;
 * composition_adjustment/compo_heap.h:82-101. */
typedef struct {
    int32_t oid, query_index;
    Hsp **hsp_array;
    int32_t hspcnt, allocated, hsp_max;
    uint8_t do_not_reallocate;
    double best_evalue;
} HspList;
typedef struct {
    int n, capacity, heap_threshold;
    double ecutoff, worst_evalue;
    void *array, *heap_array;
} CompoHeap;
static unsigned long redo_id;
static unsigned long traceback_id;
static int inside_redo;
static unsigned long active_redo;

static void dump_bytes(const char *kind, unsigned long event, const uint8_t *data, int32_t length)
{
    fprintf(stderr, "K_TRACE_%s\t%lu\t", kind, event);
    for (int32_t i = 0; i < length; ++i) fprintf(stderr, "%02x", data[i]);
    fputc('\n', stderr);
}

typedef int (*redo_fn)(void **, RedoParams *, void *, int, double,
                       MatchingSeq *, int, void *, int, int **, int,
                       void *, double *, int, double *);
int Blast_RedoOneMatch(void **out, RedoParams *params, void *incoming,
                       int hspcnt, double lambda, MatchingSeq *matching,
                       int ccat_query_length, void *query_info, int num_queries,
                       int **matrix, int alphsize, void *workspace,
                       double *pair_pvalue, int test_index, double *lambda_ratio)
{
    redo_fn real = (redo_fn)dlsym(RTLD_NEXT, "Blast_RedoOneMatch");
    unsigned long event = redo_id++;
    unsigned long previous = active_redo;
    active_redo = event;
    int index = 0;
    for (const Alignment *align = (const Alignment *)incoming;
         align && index < hspcnt; align = align->next, ++index) {
        const Hsp *hsp = align->context;
        fprintf(stderr,
                "K_TRACE_PRELIM\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                event, index, hsp ? hsp->score : -1,
                hsp ? hsp->context : -1,
                hsp ? hsp->subject.frame : -1,
                hsp ? hsp->query.offset : -1,
                hsp ? hsp->query.end : -1,
                hsp ? hsp->query.gapped_start : -1,
                hsp ? hsp->subject.offset : -1,
                hsp ? hsp->subject.end : -1,
                hsp ? hsp->subject.gapped_start : -1);
    }
    ++inside_redo;
    int status = real(out, params, incoming, hspcnt, lambda, matching,
                      ccat_query_length, query_info, num_queries, matrix,
                      alphsize, workspace, pair_pvalue, test_index, lambda_ratio);
    --inside_redo;
    active_redo = previous;
    return status;
}

typedef int16_t (*trace_fn)(int32_t, const uint8_t *, const uint8_t *,
                            GapAlign *, const void *, int32_t, int32_t,
                            int32_t, int32_t, uint8_t *);
int16_t BLAST_GappedAlignmentWithTraceback(int32_t program,
        const uint8_t *query, const uint8_t *subject, GapAlign *gap_align,
        const void *score_params, int32_t q_start, int32_t s_start,
        int32_t query_length, int32_t subject_length, uint8_t *fence_hit)
{
    trace_fn real = (trace_fn)dlsym(RTLD_NEXT, "BLAST_GappedAlignmentWithTraceback");
    unsigned long event = traceback_id++;
    if (inside_redo) {
        fprintf(stderr,
                "K_TRACE_ENTER\t%lu\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\n",
                active_redo, event, program, q_start, s_start,
                query_length, subject_length,
                gap_align ? gap_align->gap_x_dropoff : -1,
                ((const Scoring *)score_params)->gap_open,
                ((const Scoring *)score_params)->gap_extend,
                fence_hit ? 1 : 0,
                ((const Scoring *)score_params)->scale_factor);
        if (query && query_length > 0) dump_bytes("QUERY", event, query, query_length);
        if (subject && subject_length > 0) dump_bytes("SUBJECT", event, subject, subject_length);
    }
    int16_t status = real(program, query, subject, gap_align, score_params,
                          q_start, s_start, query_length, subject_length, fence_hit);
    if (inside_redo) {
        fprintf(stderr,
                "K_TRACE_RETURN\t%lu\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                active_redo, event, status,
                gap_align ? gap_align->score : -1,
                gap_align ? gap_align->query_start : -1,
                gap_align ? gap_align->query_stop : -1,
                gap_align ? gap_align->subject_start : -1,
                gap_align ? gap_align->subject_stop : -1,
                fence_hit ? *fence_hit : -1);
        const EditScript *edit = gap_align ? (const EditScript *)gap_align->edit_script : NULL;
        fprintf(stderr, "K_TRACE_EDIT\t%lu\t%d", event, edit ? edit->size : -1);
        if (edit) {
            for (int32_t i = 0; i < edit->size; ++i)
                fprintf(stderr, "\t%d:%d", edit->op_type[i], edit->num[i]);
        }
        fputc('\n', stderr);
    }
    return status;
}

/* Pinned NCBI core/blast_kappa.c:3689-3736 calls these after normalized
 * score and identity, with the local -subject seqSrc non-NULL. */
typedef int (*heap_would_fn)(CompoHeap *, double, int, int);
int BlastCompo_HeapWouldInsert(CompoHeap *heap, double evalue, int score, int oid)
{
    heap_would_fn real = (heap_would_fn)dlsym(RTLD_NEXT, "BlastCompo_HeapWouldInsert");
    int accepted = real(heap, evalue, score, oid);
    fprintf(stderr, "K_TRACE_HEAP_WOULD\t%d\t%.17g\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\n",
            oid, evalue, score, heap->n, heap->heap_threshold,
            heap->capacity, heap->ecutoff, heap->worst_evalue, accepted);
    return accepted;
}

typedef int (*heap_insert_fn)(CompoHeap *, void *, double, int, int, void **);
int BlastCompo_HeapInsert(CompoHeap *heap, void *alignments, double evalue,
                          int score, int oid, void **discarded)
{
    HspList *list = (HspList *)alignments;
    fprintf(stderr, "K_TRACE_HEAP_INSERT\t%d\t%.17g\t%d\t%d\t%d\t%.17g\n",
            oid, evalue, score, heap->n, list ? list->hspcnt : -1,
            list ? list->best_evalue : -1.);
    if (list) for (int i=0; i<list->hspcnt; ++i) {
        Hsp *h = list->hsp_array[i];
        fprintf(stderr, "K_TRACE_HEAP_HSP\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                oid, i, h->score, h->bit_score, h->evalue,
                h->num_ident, h->context, h->subject.frame,
                h->query.offset, h->query.end, h->subject.offset,
                h->subject.end, h->num);
    }
    heap_insert_fn real = (heap_insert_fn)dlsym(RTLD_NEXT, "BlastCompo_HeapInsert");
    int status = real(heap, alignments, evalue, score, oid, discarded);
    HspList *removed = discarded ? (HspList *)*discarded : NULL;
    fprintf(stderr, "K_TRACE_HEAP_INSERT_RETURN\t%d\t%d\t%d\t%.17g\t%d\n",
            oid, status, heap->n, heap->worst_evalue,
            removed ? removed->oid : -1);
    return status;
}

typedef void *(*heap_pop_fn)(CompoHeap *);
void *BlastCompo_HeapPop(CompoHeap *heap)
{
    heap_pop_fn real = (heap_pop_fn)dlsym(RTLD_NEXT, "BlastCompo_HeapPop");
    HspList *list = (HspList *)real(heap);
    fprintf(stderr, "K_TRACE_HEAP_POP\t%d\t%d\t%.17g\n",
            list ? list->oid : -1, heap->n,
            list ? list->best_evalue : -1.);
    return list;
}
