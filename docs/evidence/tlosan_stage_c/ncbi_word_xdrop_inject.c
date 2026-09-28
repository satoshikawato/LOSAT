/* Comparison-only pinned NCBI WordFinder input-state intervention.
 * c++/src/algo/blast/core/aa_ungapped.c:562-590 selects
 * word_params->cutoffs[curr_context].x_dropoff before two-hit extension.
 * Set context 1 x-drop to 1 for each actual WordFinder call, restore
 * the value afterward, and leave LOSAT runtime/build untouched.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { int32_t x_dropoff_init,x_dropoff,cutoff_score,reduced_nucl_cutoff_score; } WordCutoffs;
typedef struct { void *options; int32_t x_dropoff_max,cutoff_score_min; WordCutoffs *cutoffs; } WordParamsPrefix;
typedef short (*wordfinder_fn)(void*,void*,void*,void*,void*,void*,void*,void*,int32_t,void*,void*);
short BlastAaWordFinder(void *subject,void *query,void *query_info,void *lookup,void *matrix,void *word_params,void *ewp,void *pairs,int32_t size,void *hits,void *stats) {
 static int call=0;
 wordfinder_fn real=(wordfinder_fn)dlsym(RTLD_NEXT,"BlastAaWordFinder");
 WordParamsPrefix *params=(WordParamsPrefix*)word_params;
 int32_t old=params->cutoffs[1].x_dropoff;
 params->cutoffs[1].x_dropoff=1;
 fprintf(stderr,"WORD_XDROP_INJECT\t%d\t%d\t1\n",call++,old);
 short result=real(subject,query,query_info,lookup,matrix,word_params,ewp,pairs,size,hits,stats);
 params->cutoffs[1].x_dropoff=old;
 return result;
}
