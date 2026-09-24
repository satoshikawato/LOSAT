/* Comparison-only NCBI call-state probe; never linked into LOSAT.
 * Pinned blast_engine.c:1407,1434-1443 reads BlastSeqSrcGetTotLen and
 * invokes BLAST_OneSubjectUpdateParameters only when its return is zero.
 * Pinned blast_setup.c:1001-1024 passes compositionBasedStats=0 in update.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

int64_t BlastSeqSrcGetTotLen(const void *seq_src)
{
    int64_t (*real)(const void *) = dlsym(RTLD_NEXT, "BlastSeqSrcGetTotLen");
    int64_t result = real(seq_src);
    fprintf(stderr, "D_SEQSRC_TOTLEN\t%lld\n", (long long)result);
    return result;
}

short BLAST_OneSubjectUpdateParameters(int32_t program, uint32_t subject_length,
                                      const void *scoring_options, void *query_info,
                                      const void *sbp, void *hit_params,
                                      void *word_params, void *eff_len_params)
{
    short (*real)(int32_t, uint32_t, const void *, void *, const void *,
                  void *, void *, void *) =
        dlsym(RTLD_NEXT, "BLAST_OneSubjectUpdateParameters");
    fprintf(stderr, "D_ONE_SUBJECT_UPDATE_ENTER\t%d\t%u\n", program,
            subject_length);
    short status = real(program, subject_length, scoring_options, query_info,
                        sbp, hit_params, word_params, eff_len_params);
    fprintf(stderr, "D_ONE_SUBJECT_UPDATE_RETURN\t%d\n", status);
    return status;
}
