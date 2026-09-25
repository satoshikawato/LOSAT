/* Comparison-only pinned NCBI oracle; never linked into LOSAT.
 * c++/src/algo/blast/composition_adjustment/redo_alignment.c:103-124:
 * BlastCompo_AlignmentNew(score, matrix_adjust_rule, queryStart, queryEnd,
 *                         queryIndex, matchStart, matchEnd, frame, context)
 * writes the raw matrix adjustment rule into each owned redo alignment.
 * core/blast_kappa.c:325-342 subsequently maps that rule to HSP metadata. */
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
typedef Alignment *(*new_fn)(int, int, int, int, int, int, int, int, void *);
Alignment *BlastCompo_AlignmentNew(int score, int rule, int q0, int q1,
                                   int query_index, int s0, int s1, int frame,
                                   void *context)
{
    new_fn real = (new_fn)dlsym(RTLD_NEXT, "BlastCompo_AlignmentNew");
    Alignment *align = real(score, rule, q0, q1, query_index, s0, s1, frame, context);
    if (align)
        fprintf(stderr, "K_RULE_NEW\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                align->score, align->rule, align->query_index,
                align->query_start, align->query_end,
                align->match_start, align->match_end, align->frame);
    return align;
}
