/* Comparison-only valid API-state intervention; never linked into LOSAT.
 * Pinned c++/src/algo/blast/core/blast_traceback.c:585-605 calls
 * Blast_HSPTest after identity calculation and frees rejected HSPs.
 * Pinned c++/include/algo/blast/core/blast_options.h:391-400 gives
 * percent_identity in BlastHitSavingOptions. Set it to 100.0 only for
 * this call to exercise the real traceback deletion branch.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
typedef struct {
    double expect_value;
    int32_t cutoff_score;
    int32_t cutoff_score_fun[2];
    double percent_identity;
} OptionsPrefix;
typedef uint8_t (*hsp_test_fn)(void *, void *, int32_t);
uint8_t Blast_HSPTest(void *hsp, void *hit_options, int32_t align_length)
{
    hsp_test_fn real = (hsp_test_fn)dlsym(RTLD_NEXT, "Blast_HSPTest");
    OptionsPrefix *opts = (OptionsPrefix *)hit_options;
    double original = opts->percent_identity;
    opts->percent_identity = 100.0;
    uint8_t delete_hsp = real(hsp, hit_options, align_length);
    opts->percent_identity = original;
    fprintf(stderr, "INJECT_HSP_TEST\t%d\t%.1f\t100.0\t%d\n",
            align_length, original, delete_hsp);
    return delete_hsp;
}
