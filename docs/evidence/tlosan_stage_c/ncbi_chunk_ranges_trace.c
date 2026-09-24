/* Comparison-only NCBI input-state probe; never part of LOSAT runtime/build.
 * Pinned c++/include/algo/blast/core/blast_def.h:155-158,242-284 defines
 * SSeqRange and BLAST_SequenceBlk through seq_ranges/num_seq_ranges.
 * Pinned c++/src/algo/blast/core/blast_engine.c:259-310,478-487,804-841
 * fills chunk-local seq_ranges before calling BlastAaWordFinder.
 * Pinned c++/include/algo/blast/core/ncbi_std.h:94: Boolean is Uint1.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { int32_t left, right; } SSeqRange;
typedef struct {
    uint8_t *sequence, *sequence_start;
    int32_t length;
    int16_t frame, subject_strand;
    int32_t oid;
    uint8_t sequence_allocated, sequence_start_allocated;
    uint8_t *sequence_start_nomask, *sequence_nomask;
    uint8_t nomask_allocated;
    uint8_t *oof_sequence;
    uint8_t oof_sequence_allocated;
    uint8_t *compressed_nuc_seq, *compressed_nuc_seq_start;
    void *lcase_mask;
    uint8_t lcase_mask_allocated;
    int32_t chunk;
    uint8_t *gen_code_string;
    SSeqRange *seq_ranges;
    uint32_t num_seq_ranges;
} SubjectRanges;
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
    SubjectRanges *frame = (SubjectRanges *)subject;
    unsigned long call = call_index++;
    fprintf(stderr, "RANGE_CALL\t%lu\t%d\t%d\t%u\n",
            call, frame->frame, frame->length, frame->num_seq_ranges);
    for (uint32_t i = 0; i < frame->num_seq_ranges; ++i) {
        fprintf(stderr, "RANGE\t%lu\t%u\t%d\t%d\n", call, i,
                frame->seq_ranges[i].left, frame->seq_ranges[i].right);
    }
    return real(subject, query, query_info, lookup, matrix, word_params,
            ewp, offset_pairs, offset_array_size, init_hitlist, ungapped_stats);
}

/* Pinned c++/src/algo/blast/core/blast_util.c:182-217: local-subject
 * setup copies the mask complement and fills its first/last boundaries. */
typedef short (*set_ranges_fn)(void*, SSeqRange*, uint32_t, uint8_t, int32_t);
short BlastSeqBlkSetSeqRanges(void* subject, SSeqRange* ranges,
        uint32_t count, uint8_t copy, int32_t mask_type)
{
    set_ranges_fn real = (set_ranges_fn)dlsym(RTLD_NEXT, "BlastSeqBlkSetSeqRanges");
    if (!real) return -1;
    short status = real(subject, ranges, count, copy, mask_type);
    SubjectRanges *frame = (SubjectRanges *)subject;
    fprintf(stderr, "SET_RANGES\t%u\n", frame->num_seq_ranges);
    for (uint32_t i = 0; i < frame->num_seq_ranges; ++i) {
        fprintf(stderr, "SET_RANGE\t%u\t%d\t%d\n", i,
                frame->seq_ranges[i].left, frame->seq_ranges[i].right);
    }
    return status;
}
