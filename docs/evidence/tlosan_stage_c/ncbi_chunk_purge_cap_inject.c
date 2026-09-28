/* Comparison-only input-state intervention at pinned NCBI real call sites.
 * blast_engine.c:539-552 purges chunk HSPs before score sorting, and
 * blast_engine.c:840-850 appends frame lists with kHspNumMax. The first
 * purge input is changed to share a lower-scoring HSP's query/subject
 * start with the top HSP; the append argument is set to three.
 * No part of this probe is linked or called by LOSAT.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { int16_t frame; int32_t offset,end,gapped_start; } Seg;
typedef struct { int32_t score,num_ident; double bit_score,evalue; Seg query,subject; int32_t context; void *gap_info; } Hsp;
typedef struct { int32_t oid,query_index; Hsp **hsp_array; int32_t hspcnt,allocated,hsp_max,do_not_reallocate; double best_evalue; } HspList;
typedef int32_t (*purge_fn)(int32_t,HspList*,uint8_t);
int32_t Blast_HSPListPurgeHSPsWithCommonEndpoints(int32_t program,HspList *list,uint8_t purge)
{
    static int intervened=0;
    purge_fn real=(purge_fn)dlsym(RTLD_NEXT,"Blast_HSPListPurgeHSPsWithCommonEndpoints");
    if (!intervened && purge && list->hspcnt==3 &&
        list->hsp_array[0]->score==656 && list->hsp_array[2]->score==20 &&
        list->hsp_array[0]->context==list->hsp_array[2]->context) {
        Hsp *top=list->hsp_array[0],*lower=list->hsp_array[2];
        lower->query.offset=top->query.offset;
        lower->subject.offset=top->subject.offset;
        intervened=1;
        fprintf(stderr,"INJECT_ENDPOINT\t%d\t%d\t%d\t%d\t%d\n",
                lower->context,lower->score,lower->query.offset,
                lower->subject.offset,list->hspcnt);
    }
    return real(program,list,purge);
}
typedef short (*append_fn)(HspList**,HspList**,int32_t);
short Blast_HSPListAppend(HspList **incoming,HspList **combined,int32_t cap)
{
    static int call=0;
    append_fn real=(append_fn)dlsym(RTLD_NEXT,"Blast_HSPListAppend");
    fprintf(stderr,"INJECT_APPEND\t%d\t%d\t3\n",call++,cap);
    return real(incoming,combined,3);
}
