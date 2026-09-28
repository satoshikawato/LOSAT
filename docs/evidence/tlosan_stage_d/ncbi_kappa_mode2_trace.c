/* Comparison-only Kappa-mode-2 call probe; never linked into LOSAT.
 * Pinned composition_adjustment/redo_alignment.h:57-72,328-354,466-481;
 * core/blast_kappa.c:2423-2478,3625-3645,1475-1565;
 * composition_adjustment/composition_adjustment.h:280-294.
 */
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
typedef struct {
    void *matrix_info, *gapping_params;
    int32_t compo_mode, position_based, re_pseudocounts;
    int32_t subject_is_translated, query_is_translated;
    int32_t ccat_query_length, cutoff_s;
    double cutoff_e;
    int32_t do_link_hsps;
    void *callbacks;
    double near_identical_cutoff;
} RedoParams;
typedef struct { int32_t length, index; void *local_data; } MatchingSeq;

static unsigned long event_no;
static int inside_redo;

static void dump_alignments(unsigned long event, const char *stage,
                            int query, const Alignment *align)
{
    for (int i = 0; align && i < 1024; ++i, align = align->next) {
        fprintf(stderr,
                "K_ALIGN\t%lu\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                event, stage, query, i, align->score, align->rule,
                align->query_index, align->query_start, align->query_end,
                align->match_start, align->match_end, align->frame);
    }
}

typedef int (*redo_fn)(Alignment **, RedoParams *, Alignment *, int,
                       double, MatchingSeq *, int, void *, int, int **,
                       int, void *, double *, int, double *);
int Blast_RedoOneMatch(Alignment **out, RedoParams *params,
                       Alignment *incoming, int hspcnt, double lambda,
                       MatchingSeq *matching, int ccat_query_length,
                       void *query_info, int num_queries, int **matrix,
                       int alphsize, void *workspace, double *pair_pvalue,
                       int test_index, double *lambda_ratio)
{
    redo_fn real = (redo_fn)dlsym(RTLD_NEXT, "Blast_RedoOneMatch");
    unsigned long event = event_no++;
    fprintf(stderr,
            "K_CALL\t%lu\tredo_enter\t%d\t%.17g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\t%d\t%d\n",
            event, hspcnt, lambda, ccat_query_length, num_queries,
            matching ? matching->length : -1, params ? params->compo_mode : -1,
            params ? params->subject_is_translated : -1,
            params ? params->query_is_translated : -1,
            params ? params->cutoff_s : -1,
            params ? params->do_link_hsps : -1,
            params ? params->cutoff_e : -1., alphsize, test_index);
    dump_alignments(event, "incoming", -1, incoming);
    ++inside_redo;
    int status = real(out, params, incoming, hspcnt, lambda, matching,
                      ccat_query_length, query_info, num_queries, matrix,
                      alphsize, workspace, pair_pvalue, test_index, lambda_ratio);
    --inside_redo;
    fprintf(stderr, "K_CALL\t%lu\tredo_return\t%d\t%.17g\t%.17g\n",
            event, status, pair_pvalue ? *pair_pvalue : -1.,
            lambda_ratio ? *lambda_ratio : -1.);
    for (int query = 0; out && query < num_queries; ++query)
        dump_alignments(event, "redone", query, out[query]);
    return status;
}

typedef int (*partial_fn)(const uint8_t *, int32_t, int16_t, const uint8_t *,
                          uint8_t **, int32_t *, uint8_t **);
int Blast_GetPartialTranslation(const uint8_t *nucl, int32_t nucl_length,
                                int16_t frame, const uint8_t *genetic_code,
                                uint8_t **translation, int32_t *protein_length,
                                uint8_t **mixed)
{
    partial_fn real = (partial_fn)dlsym(RTLD_NEXT, "Blast_GetPartialTranslation");
    int status = real(nucl, nucl_length, frame, genetic_code,
                      translation, protein_length, mixed);
    if (inside_redo) {
        unsigned long event = event_no++;
        fprintf(stderr, "K_CALL\t%lu\tpartial_translation\t%d\t%d\t%d\t%d\n",
                event, nucl_length, frame, status,
                protein_length ? *protein_length : -1);
        if (status == 0 && translation && *translation && protein_length) {
            fprintf(stderr, "K_TRANSLATED\t%lu\t", event);
            for (int32_t i = 0; i < *protein_length; ++i)
                fprintf(stderr, "%02x", (*translation)[i+1]);
            fputc('\n', stderr);
        }
    }
    return status;
}

typedef int (*adjust_fn)(int **, const void *, int, const void *, int,
                         const void *, int, int, void *, int *, void *,
                         double *, int, double *);
int Blast_AdjustScores(int **matrix, const void *query_comp, int query_length,
                       const void *subject_comp, int subject_length,
                       const void *matrix_info, int mode, int pseudocounts,
                       void *workspace, int *rule, void *calc_lambda,
                       double *pair_pvalue, int test_index, double *lambda_ratio)
{
    adjust_fn real = (adjust_fn)dlsym(RTLD_NEXT, "Blast_AdjustScores");
    int status = real(matrix, query_comp, query_length, subject_comp,
                      subject_length, matrix_info, mode, pseudocounts,
                      workspace, rule, calc_lambda, pair_pvalue,
                      test_index, lambda_ratio);
    if (inside_redo) {
        fprintf(stderr,
                "K_CALL\t%lu\tadjust_scores\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\t%.17g\n",
                event_no++, query_length, subject_length, mode,
                pseudocounts, status, rule ? *rule : -1,
                pair_pvalue ? *pair_pvalue : -1.,
                lambda_ratio ? *lambda_ratio : -1.);
    }
    return status;
}
