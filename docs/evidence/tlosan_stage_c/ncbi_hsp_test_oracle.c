/* Comparison-only pinned NCBI C API oracle for traceback hit deletion.
 * Never linked, built, or called by LOSAT.
 * Pinned c++/src/algo/blast/core/blast_hits.c:993-1001,1027-1032:
 * Blast_HSPTest compares num_ident*100.0 with align_length*percent_identity,
 * then compares align_length with min_hit_length.
 * Structure prefixes: c++/include/algo/blast/core/blast_hits.h:96-143 and
 * c++/include/algo/blast/core/blast_options.h:391-427.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { int16_t frame; int32_t offset, end, gapped_start; } Seg;
typedef struct {
    int32_t score, num_ident;
    double bit_score, evalue;
    Seg query, subject;
    int32_t context;
    void *gap_info;
} Hsp;
typedef struct {
    double expect_value;
    int32_t cutoff_score;
    int32_t cutoff_score_fun[2];
    double percent_identity;
    int32_t max_edit_distance, hitlist_size, hsp_num_max,
            total_hsp_limit, culling_limit, mask_level, do_sum_stats,
            longest_intron, min_hit_length;
} OptionsPrefix;
typedef int (*hsp_test_fn)(Hsp *, const OptionsPrefix *, int32_t);

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!library) { fprintf(stderr, "%s\n", dlerror()); return 1; }
    hsp_test_fn test = (hsp_test_fn)dlsym(library, "Blast_HSPTest");
    if (!test) return 1;
    const struct { int num_ident, align_length, min_length; double pct; } cases[] = {
        {10, 10, 0, 100.0}, {8, 10, 0, 80.0},
        {8, 10, 0, 80.1}, {8, 10, 11, 0.0},
        {8, 10, 10, 0.0}, {0, 0, 0, 0.0},
    };
    for (unsigned i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        Hsp hsp = {0};
        OptionsPrefix options = {0};
        hsp.num_ident = cases[i].num_ident;
        options.percent_identity = cases[i].pct;
        options.min_hit_length = cases[i].min_length;
        printf("HSP_TEST\t%d\t%d\t%.1f\t%d\t%d\n",
               cases[i].num_ident, cases[i].align_length, cases[i].pct,
               cases[i].min_length, test(&hsp, &options, cases[i].align_length));
    }
    dlclose(library);
    return 0;
}
