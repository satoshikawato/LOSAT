/* Comparison-only NCBI Blast_AdjustScores input probe; never linked into LOSAT.
 * Pinned composition_adjustment/composition_adjustment.h:51,54-59,76-94,280-294;
 * composition_adjustment/redo_alignment.c:1219-1254.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    double prob[28];
    int32_t num_true;
} Composition;

typedef struct {
    char *name;
    int **start_matrix;
    double **start_freq_ratios;
    int32_t rows, cols, position_based;
    double ungapped_lambda;
} MatrixInfo;

typedef int (*adjust_fn)(int **, const Composition *, int, const Composition *, int,
                         const void *, int, int, void *, int *, void *,
                         double *, int, double *);

static unsigned long event_no;
static void dump_composition(unsigned long event, const char *which,
                             const Composition *comp)
{
    fprintf(stderr, "K_COMP\t%lu\t%s\t%d", event, which,
            comp ? comp->num_true : -1);
    if (comp) {
        for (int i = 0; i < 28; ++i) {
            union { double value; uint64_t bits; } bitcast = { comp->prob[i] };
            fprintf(stderr, "\t%016llx", (unsigned long long)bitcast.bits);
        }
    }
    fputc('\n', stderr);
}

int Blast_AdjustScores(int **matrix, const Composition *query_comp,
                       int query_length, const Composition *subject_comp,
                       int subject_length, const void *matrix_info, int mode,
                       int pseudocounts, void *workspace, int *rule,
                       void *calc_lambda, double *pair_pvalue,
                       int test_index, double *lambda_ratio)
{
    adjust_fn real = (adjust_fn)dlsym(RTLD_NEXT, "Blast_AdjustScores");
    unsigned long event = event_no++;
    fprintf(stderr, "K_COMP_CALL\t%lu\t%d\t%d\t%d\t%d\t%d\n",
            event, query_length, subject_length, mode, pseudocounts, test_index);
    const MatrixInfo *info = (const MatrixInfo *)matrix_info;
    fprintf(stderr, "K_MATRIX\t%lu\t%d\t%d\t%d\t%.17g\t%d\t%d\t%d\n",
            event, info->rows, info->cols, info->position_based,
            info->ungapped_lambda, info->start_matrix[1][1],
            info->start_matrix[1][16], info->start_matrix[22][22]);
    dump_composition(event, "query", query_comp);
    dump_composition(event, "subject", subject_comp);
    int status = real(matrix, query_comp, query_length, subject_comp,
                      subject_length, matrix_info, mode, pseudocounts,
                      workspace, rule, calc_lambda, pair_pvalue,
                      test_index, lambda_ratio);
    fprintf(stderr, "K_COMP_RESULT\t%lu\t%d\t%d\t%.17g\n", event,
            status, rule ? *rule : -1, lambda_ratio ? *lambda_ratio : -1.);
    fprintf(stderr, "K_ADJUSTED\t%lu", event);
    for (int row = 0; row < info->rows; ++row)
        for (int col = 0; col < info->cols; ++col)
            fprintf(stderr, "\t%d", matrix[row][col]);
    fputc('\n', stderr);
    return status;
}
