/* Comparison-only pinned NCBI C API oracle for Blast_HSPListAppend's cap.
 * Never linked, built, or called by LOSAT.
 * Pinned c++/src/algo/blast/core/blast_hits.c:151-188,1558-1574,
 * 2762-2864: construct preliminary HSPs, save two lists, and combine
 * them by ScoreCompareHSPs under hsp_num_max=3.
 * Pinned c++/include/algo/blast/core/blast_hits.h:96-160,232-253:
 * BlastSeg/BlastHSP/BlastHSPList prefixes and Blast_HSPInit parameters.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

typedef HspList *(*list_new_fn)(int32_t);
typedef short (*hsp_init_fn)(int32_t, int32_t, int32_t, int32_t,
                             int32_t, int32_t, int32_t, int16_t,
                             int16_t, int32_t, void **, Hsp **);
typedef short (*save_fn)(HspList *, Hsp *);
typedef short (*append_fn)(HspList **, HspList **, int32_t);
typedef HspList *(*list_free_fn)(HspList *);

static void add_hsp(HspList *list, hsp_init_fn init, save_fn save,
                    int score, int context, int frame, int subject_start)
{
    Hsp *hsp = NULL;
    if (init(0, 20, subject_start, subject_start + 20, 5,
             subject_start + 5, context, 0, frame, score, NULL, &hsp)
        || !hsp || save(list, hsp)) {
        fprintf(stderr, "NCBI HSP constructor/save failed\n");
        exit(1);
    }
}

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!library) {
        fprintf(stderr, "%s\n", dlerror());
        return 1;
    }
    list_new_fn make = (list_new_fn)dlsym(library, "Blast_HSPListNew");
    hsp_init_fn init = (hsp_init_fn)dlsym(library, "Blast_HSPInit");
    save_fn save = (save_fn)dlsym(library, "Blast_HSPListSaveHSP");
    append_fn append = (append_fn)dlsym(library, "Blast_HSPListAppend");
    list_free_fn free_list = (list_free_fn)dlsym(library, "Blast_HSPListFree");
    if (!make || !init || !save || !append || !free_list) return 1;

    HspList *old = make(0);
    HspList *incoming = make(0);
    add_hsp(old, init, save, 100, 0, 1, 100);
    add_hsp(old, init, save, 80, 0, 1, 300);
    add_hsp(old, init, save, 60, 0, 1, 500);
    add_hsp(incoming, init, save, 90, 1, 2, 200);
    add_hsp(incoming, init, save, 80, 1, 2, 300);
    add_hsp(incoming, init, save, 70, 1, 2, 400);
    printf("CAP_INPUT\t3\t%d\t%d\n", old->hspcnt, incoming->hspcnt);
    short status = append(&incoming, &old, 3);
    printf("CAP_OUTPUT\t%d\t%d\n", status, old->hspcnt);
    for (int32_t i = 0; i < old->hspcnt; ++i) {
        Hsp *hsp = old->hsp_array[i];
        printf("CAP_HSP\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
               i, hsp->context, hsp->score, hsp->subject.frame,
               hsp->query.offset, hsp->query.end,
               hsp->subject.offset, hsp->subject.end);
    }
    if (incoming) free_list(incoming);
    free_list(old);
    dlclose(library);
    return status ? 1 : 0;
}
