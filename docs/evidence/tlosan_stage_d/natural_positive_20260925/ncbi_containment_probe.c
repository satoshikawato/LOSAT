/* Comparison-only probe. Pinned NCBI core/blast_kappa.c:3577-3680:
 * Blast_RedoOneMatch(..., alignments, ...); s_HSPListFromDistinctAlignments(...);
 * if (hsp_list->hspcnt > 1) s_HitlistReapContained(...);
 * s_HitlistEvaluateAndPurge(...) calls BLAST_LinkHsps.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct Alignment {
 int32_t score, rule, query_index, query_start, query_end;
 int32_t match_start, match_end, frame;
 void *context;
 struct Alignment *next;
} Alignment;
typedef struct { int32_t oid, query_index; void **hsp_array; int32_t hspcnt; } HspList;
static unsigned long event;
static int redone;
typedef int (*redo_fn)(void **, void *, void *, int, double, void *, int, void *, int, int **, int, void *, double *, int, double *);
int Blast_RedoOneMatch(void **out, void *params, void *incoming, int hspcnt, double lambda, void *matching,
                       int ccat_query_length, void *query_info, int num_contexts, int **matrix, int alphsize,
                       void *workspace, double *pair_pvalue, int test_index, double *lambda_ratio) {
 redo_fn real=(redo_fn)dlsym(RTLD_NEXT,"Blast_RedoOneMatch");
 int status=real(out,params,incoming,hspcnt,lambda,matching,ccat_query_length,query_info,num_contexts,matrix,alphsize,workspace,pair_pvalue,test_index,lambda_ratio);
 int count=0;
 for(int q=0;q<num_contexts;q++) for(Alignment *a=(Alignment*)out[q];a;a=a->next) {count++; fprintf(stderr,"CP_A_HSP\t%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",event,q,a->score,a->query_start,a->query_end,a->match_start,a->match_end,a->frame,a->rule);}
 fprintf(stderr,"CP_A\t%lu\t%d\t%d\t%d\n",event++,hspcnt,count,status);
 redone=1;
 return status;
}
typedef short (*link_fn)(int32_t, HspList *, const void *, int32_t, const void *, const void *, uint8_t);
short BLAST_LinkHsps(int32_t program,HspList *list,const void *query,int32_t length,const void *sbp,const void *params,uint8_t gapped){
 link_fn real=(link_fn)dlsym(RTLD_NEXT,"BLAST_LinkHsps");
 if(redone) fprintf(stderr,"CP_L\t%d\t%d\n",list?list->oid:-1,list?list->hspcnt:-1);
 redone=0;
 return real(program,list,query,length,sbp,params,gapped);
}
