/* Comparison-only translated Kappa SEG sequence trace; never linked into LOSAT.
 * Pinned core/blast_kappa.c:1428-1455,1475-1549,2423-2478;
 * core/blast_filter.h:250-265; core/blast_util.c:1141-1205.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef int (*redo_fn)(void **, void *, void *, int, double, void *, int,
                       void *, int, int **, int, void *, double *, int, double *);
typedef int (*partial_fn)(const uint8_t *, int32_t, int16_t, const uint8_t *,
                          uint8_t **, int32_t *, uint8_t **);
typedef void (*mask_fn)(uint8_t *, int32_t, int32_t, const void *, int32_t, int32_t);

static int inside_redo;
static unsigned long redo_event, translation_event;

int Blast_RedoOneMatch(void **out, void *params, void *incoming, int hspcnt,
                       double lambda, void *matching, int ccat_length,
                       void *query_info, int num_queries, int **matrix,
                       int alphsize, void *workspace, double *pair_pvalue,
                       int test_index, double *lambda_ratio)
{
    redo_fn real = (redo_fn)dlsym(RTLD_NEXT, "Blast_RedoOneMatch");
    unsigned long saved = redo_event++;
    inside_redo = 1;
    int status = real(out, params, incoming, hspcnt, lambda, matching,
                      ccat_length, query_info, num_queries, matrix,
                      alphsize, workspace, pair_pvalue, test_index, lambda_ratio);
    inside_redo = 0;
    fprintf(stderr, "K_SEG_REDO\t%lu\t%d\t%d\n", saved, hspcnt, status);
    return status;
}

int Blast_GetPartialTranslation(const uint8_t *nucl, int32_t nucl_length,
                                int16_t frame, const uint8_t *genetic_code,
                                uint8_t **translation, int32_t *protein_length,
                                uint8_t **mixed)
{
    partial_fn real = (partial_fn)dlsym(RTLD_NEXT, "Blast_GetPartialTranslation");
    int status = real(nucl, nucl_length, frame, genetic_code,
                      translation, protein_length, mixed);
    if (inside_redo && status == 0 && translation && *translation && protein_length) {
        unsigned long event = translation_event++;
        fprintf(stderr, "K_SEG_RAW\t%lu\t%d\t%d\t%d\t", event,
                nucl_length, frame, *protein_length);
        for (int32_t i = 0; i < *protein_length; ++i)
            fprintf(stderr, "%02x", (*translation)[i+1]);
        fputc('\n', stderr);
    }
    return status;
}

void Blast_MaskTheResidues(uint8_t *buffer, int32_t length, int32_t is_na,
                            const void *mask_loc, int32_t reverse, int32_t offset)
{
    mask_fn real = (mask_fn)dlsym(RTLD_NEXT, "Blast_MaskTheResidues");
    real(buffer, length, is_na, mask_loc, reverse, offset);
    if (inside_redo && !is_na) {
        fprintf(stderr, "K_SEG_MASKED\t%lu\t%d\t", translation_event - 1, length);
        for (int32_t i = 0; i < length; ++i)
            fprintf(stderr, "%02x", buffer[i]);
        fputc('\n', stderr);
    }
}
