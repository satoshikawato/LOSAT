/* Comparison-only probe; never linked into LOSAT.
 * Pinned NCBI c++/src/algo/blast/core/blast_setup.c:1011-1024 calls
 * BlastHitSavingParametersUpdate after BLAST_CalcEffLengths.
 * Pinned blast_parameters.c:902-999 chooses Spouge or BLAST_Cutoffs,
 * then reduces the cutoff for sum statistics.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    double lambda, c, g, a, alpha, sigma, a_un, alpha_un;
    double b, beta, tau;
    int64_t db_length;
    uint8_t filled;
} Gumbel;

typedef struct { int32_t cutoff_score, cutoff_score_max; } Cutoffs;
typedef struct {
    void *options;
    int32_t cutoff_score_min;
    Cutoffs *cutoffs;
    void *link_hsp_params;
    uint8_t restricted_align, do_sum_stats;
    int32_t mask_level;
    int32_t *low_score;
    double prelim_evalue;
} HitParams;
typedef struct {
    int32_t first_context, last_context, num_queries;
    void *contexts;
} QueryInfo;

typedef int32_t (*spouge_fn)(double, void *, Gumbel *, int32_t, int32_t);
int32_t BLAST_SpougeEtoS(double evalue, void *karlin, Gumbel *gumbel,
                         int32_t query_length, int32_t subject_length)
{
    spouge_fn real = (spouge_fn)dlsym(RTLD_NEXT, "BLAST_SpougeEtoS");
    int32_t result = real(evalue, karlin, gumbel, query_length, subject_length);
    fprintf(stderr, "D_PARAM_SPOUGE\t%.17g\t%d\t%d\t%lld\t%d\t%d\n",
            evalue, query_length, subject_length,
            (long long)gumbel->db_length, gumbel->filled, result);
    return result;
}

typedef short (*hit_fn)(int32_t, const void *, const QueryInfo *, int32_t,
                        int32_t, HitParams *);
short BlastHitSavingParametersUpdate(int32_t program, const void *score,
                                    const QueryInfo *query, int32_t subject_length,
                                    int32_t composition, HitParams *params)
{
    hit_fn real = (hit_fn)dlsym(RTLD_NEXT, "BlastHitSavingParametersUpdate");
    short status = real(program, score, query, subject_length, composition, params);
    fprintf(stderr, "D_PARAM_HIT\t%d\t%d\t%d\t%d\t%.17g\t%d\n",
            program, subject_length, composition, params->do_sum_stats,
            params->prelim_evalue, params->cutoff_score_min);
    for (int32_t i = query->first_context; i <= query->last_context; ++i)
        fprintf(stderr, "D_PARAM_CONTEXT\t%d\t%d\t%d\n", i,
                params->cutoffs[i].cutoff_score,
                params->cutoffs[i].cutoff_score_max);
    return status;
}
