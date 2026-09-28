/* Comparison-only NCBI traceback input intervention.
 * Pinned blast_traceback.c:436-449 calls BlastGetOffsetsForGappedAlignment
 * when both gapped starts are zero, then frees the HSP on FALSE.
 * Pinned blast_gapalign.c:3248-3321 computes its score-window predicate.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { int16_t frame; int32_t offset,end,gapped_start; } Seg;
typedef struct { int32_t score,num_ident; double bit_score,evalue; Seg query,subject; int32_t context; void *gap_info; } Hsp;
typedef struct { int32_t oid,query_index; Hsp **hsp_array; int32_t hspcnt,allocated,hsp_max,do_not_reallocate; double best_evalue; } HspList;
typedef short (*traceback_fn)(int,HspList*,void*,void*,void*,void*,void*,void*,void*,void*,const uint8_t*,uint8_t*);
short Blast_TracebackFromHSPList(int program,HspList *list,void *query,void *subject,void *qi,void *ga,void *sb,void *sp,void *eo,void *hp,const uint8_t *gc,uint8_t *fence)
{
    traceback_fn real=(traceback_fn)dlsym(RTLD_NEXT,"Blast_TracebackFromHSPList");
    if (list->query_index==0 && list->hspcnt>0) {
        Hsp *h=list->hsp_array[0];
        const int positive=getenv("TLOSAN_START_POSITIVE") != NULL;
        const int so=positive ? 200 : 0;
        const int length=positive ? 120 : 20;
        h->query.offset=0;h->query.end=length;h->query.gapped_start=0;
        h->subject.frame=1;h->subject.offset=so;h->subject.end=so+length;h->subject.gapped_start=0;
        fprintf(stderr,"INJECT_START\t%d\t%d\t%d\t%d\n",list->query_index,h->query.offset,h->subject.offset,h->subject.end);
    }
    short retval=real(program,list,query,subject,qi,ga,sb,sp,eo,hp,gc,fence);
    if (list->query_index==0 && list->hspcnt>0) {
        Hsp *h=list->hsp_array[0];
        fprintf(stderr,"AFTER_START_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",list->query_index,h->score,h->query.offset,h->query.end,h->query.gapped_start,h->subject.offset,h->subject.end,h->subject.gapped_start);
    }
    return retval;
}
typedef uint8_t (*start_fn)(const uint8_t*,const uint8_t*,const void*,Hsp*,int32_t*,int32_t*);
uint8_t BlastGetOffsetsForGappedAlignment(const uint8_t *query,const uint8_t *subject,const void *sb,Hsp *h,int32_t *q,int32_t *s)
{
    start_fn real=(start_fn)dlsym(RTLD_NEXT,"BlastGetOffsetsForGappedAlignment");
    uint8_t retval=real(query,subject,sb,h,q,s);
    fprintf(stderr,"START_RESULT\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",retval,h->query.offset,h->query.end,h->subject.offset,h->subject.end,retval?*q:-1,retval?*s:-1);
    return retval;
}
