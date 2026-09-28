/* Comparison-only probe. Never link this into LOSAT.
 * NCBI c++/src/algo/blast/core/blast_engine.c:804-841:
 *   subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
 *   subject->sequence = translation_buffer + frame_offsets[context] + 1;
 *   subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
 *   status = s_BlastSearchEngineOneContext(...);
 * NCBI c++/include/algo/blast/core/blast_def.h:242-247:
 *   Uint1* sequence; Uint1* sequence_start; Int4 length; Int2 frame;
 * The WordFinder entry receives the translated preliminary-search subject.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    uint8_t *sequence;
    uint8_t *sequence_start;
    int32_t length;
    int16_t frame;
} SubjectPrefix;

typedef short (*wordfinder_fn)(void*, void*, void*, void*, void*, void*,
        void*, void*, int32_t, void*, void*);

short BlastAaWordFinder(void* subject, void* query, void* query_info,
        void* lookup, void* matrix, void* word_params, void* ewp,
        void* offset_pairs, int32_t offset_array_size,
        void* init_hitlist, void* ungapped_stats)
{
    static unsigned long call_index = 0;
    wordfinder_fn real = (wordfinder_fn)dlsym(RTLD_NEXT, "BlastAaWordFinder");
    if (!real) return -1;
    SubjectPrefix* s = (SubjectPrefix*)subject;
    unsigned long call = call_index++;
    fprintf(stderr, "FRAME\t%lu\t%d\t%d\t", call, s->frame, s->length);
    for (int32_t i = 0; i < s->length; ++i) {
        fprintf(stderr, "%02x", s->sequence[i]);
    }
    fputc('\n', stderr);
    return real(subject, query, query_info, lookup, matrix, word_params,
            ewp, offset_pairs, offset_array_size, init_hitlist, ungapped_stats);
}
