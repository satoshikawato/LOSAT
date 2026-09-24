/* Comparison-only ErfC probe. Pinned blast_stat.c:5216,5223 calls ErfC;
 * pinned boost_erf.c:130-233 evaluates the polynomial with expl(-z*z).
 * Never linked into LOSAT.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>

double ErfC(double z)
{
    double (*real)(double) = dlsym(RTLD_NEXT, "ErfC");
    double result = real(z);
    fprintf(stderr, "D_ERFC\t%.17g\t%.17g\n", z, result);
    return result;
}
