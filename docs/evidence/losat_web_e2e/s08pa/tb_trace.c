/* Comparison-only probe (S08+a, TN-2). Never linked by LOSAT.
 * NCBI blast_traceback.c:259-721 (Blast_TracebackFromHSPList),
 * blast_hits.c:1147-1232 (Blast_HSPGetTargetTranslation).
 * Layout: blast_def.h:242-254,311-319; blast_hits.h:111-165. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { int16_t frame; int32_t offset, end, gapped_start; } Seg;
typedef struct {
    int32_t score, num_ident; double bit_score, evalue;
    Seg query, subject; int32_t context; void *gap_info; int32_t num;
    int16_t comp_adjustment_method;
} Hsp;
typedef struct {
    int32_t oid, query_index; Hsp **hsp_array;
    int32_t hspcnt, allocated, hsp_max; uint8_t do_not_reallocate; double best_evalue;
} HspList;
typedef struct { uint8_t *sequence, *sequence_start; int32_t length; } SeqBlk;
typedef struct {
    int32_t program_number; const uint8_t *gen_code_string; uint8_t **translations;
    uint8_t partial; int32_t num_frames; int32_t *range; SeqBlk *subject_blk;
} TargetT;
static unsigned long ev;
typedef short (*tb_fn)(int32_t, HspList *, const void *, SeqBlk *, const void *,
                       void *, const void *, const void *, const void *,
                       const void *, const uint8_t *, uint8_t *);
short Blast_TracebackFromHSPList(int32_t program, HspList *list, const void *qblk,
                                 SeqBlk *sblk, const void *qi, void *ga,
                                 const void *sbp, const void *sp, const void *eo,
                                 const void *hp, const uint8_t *gc, uint8_t *fence)
{
    tb_fn real = (tb_fn)dlsym(RTLD_NEXT, "Blast_TracebackFromHSPList");
    unsigned long n = ev++;
    fprintf(stderr, "T_CALL\t%lu\t%d\t%d\t%d\t%d\n", n, list->oid, list->query_index,
            sblk->length, fence ? *fence : -1);
    for (int i = 0; i < list->hspcnt; ++i) {
        Hsp *h = list->hsp_array[i];
        fprintf(stderr, "T_IN\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n, i,
                h->context, h->subject.frame, h->query.offset, h->query.end,
                h->subject.offset, h->subject.end, h->score,
                h->query.gapped_start, h->subject.gapped_start);
    }
    short st = real(program, list, qblk, sblk, qi, ga, sbp, sp, eo, hp, gc, fence);
    fprintf(stderr, "T_RET\t%lu\t%d\t%d\n", n, st, fence ? *fence : -1);
    return st;
}
typedef const uint8_t *(*tt_fn)(TargetT *, const Hsp *, int32_t *);
const uint8_t *Blast_HSPGetTargetTranslation(TargetT *t, const Hsp *h, int32_t *len)
{
    tt_fn real = (tt_fn)dlsym(RTLD_NEXT, "Blast_HSPGetTargetTranslation");
    const uint8_t *r = real(t, h, len);
    fprintf(stderr, "T_XL\t%d\t%d\t%d\t%d\t%d\n", h ? h->subject.frame : 0,
            h ? h->subject.offset : 0, h ? h->subject.end : 0, len ? *len : -1,
            t->subject_blk ? t->subject_blk->length : -1);
    return r;
}
