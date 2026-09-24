/* Comparison-only pinned NCBI Kappa stream and early-termination probe.
 * core/blast_kappa.c:3383-3427,3525-3535;
 * core/blast_hspstream.c:271-319;
 * composition_adjustment/redo_alignment.c:1560-1582;
 * composition_adjustment/compo_heap.c:405-410.
 * Never linked into LOSAT runtime or build.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    int32_t oid, query_index;
    void *hsp_array;
    int32_t hspcnt, allocated, hsp_max;
    uint8_t do_not_reallocate;
    double best_evalue;
} HspList;
typedef struct {
    int n, capacity, heap_threshold;
    double ecutoff, worst_evalue;
    void *array, *heap_array;
} CompoHeap;

typedef int (*stream_read_fn)(void *, HspList **);
int BlastHSPStreamRead(void *stream, HspList **out)
{
    stream_read_fn real = (stream_read_fn)dlsym(RTLD_NEXT, "BlastHSPStreamRead");
    int status = real(stream, out);
    HspList *list = out ? *out : NULL;
    fprintf(stderr, "K_EARLY_STREAM\t%d\t%d\t%d\t%.17g\t%d\n",
            status, list ? list->oid : -1, list ? list->query_index : -1,
            list ? list->best_evalue : -1., list ? list->hspcnt : -1);
    return status;
}

typedef int (*early_fn)(double, CompoHeap *, int);
int BlastCompo_EarlyTermination(double evalue, CompoHeap *heaps, int num_queries)
{
    early_fn real = (early_fn)dlsym(RTLD_NEXT, "BlastCompo_EarlyTermination");
    int terminated = real(evalue, heaps, num_queries);
    for (int q = 0; q < num_queries; ++q) {
        CompoHeap *heap = &heaps[q];
        fprintf(stderr, "K_EARLY_EVAL\t%d\t%.17g\t%d\t%d\t%.17g\t%.17g\t%d\n",
                q, evalue, heap->n, heap->heap_threshold, heap->worst_evalue,
                heap->ecutoff, terminated);
    }
    return terminated;
}
