/* Comparison-only pinned NCBI BLAST_SpougeStoE input-state probe.
 * c++/src/algo/blast/core/blast_hits.c:1855-1918 passes HSP score,
 * query length and translated subject length into blast_stat.c:5176-5231.
 * Layout: c++/include/algo/blast/core/blast_stat.h:65-74,94-112.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { double lambda, k, logk, h, paramc; } Karlin;
typedef struct {
    double lambda, c, g, a, alpha, sigma, a_un, alpha_un;
    double b, beta, tau;
    int64_t db_length;
    uint8_t filled;
} Gumbel;
typedef double (*spouge_fn)(int32_t, Karlin *, Gumbel *, int32_t, int32_t);
double BLAST_SpougeStoE(int32_t score, Karlin *karlin, Gumbel *gumbel,
                       int32_t query_length, int32_t subject_length)
{
    spouge_fn real = (spouge_fn)dlsym(RTLD_NEXT, "BLAST_SpougeStoE");
    double result = real(score, karlin, gumbel, query_length, subject_length);
    fprintf(stderr, "D_SPOUGE\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%lld\t%d\t%.17g\n",
            score, query_length, subject_length,
            karlin->lambda, karlin->k, karlin->logk,
            gumbel->lambda, gumbel->c, gumbel->g, gumbel->a,
            gumbel->alpha, gumbel->sigma, gumbel->a_un,
            gumbel->alpha_un, gumbel->b, gumbel->beta, gumbel->tau,
            (long long)gumbel->db_length, gumbel->filled, result);
    return result;
}
