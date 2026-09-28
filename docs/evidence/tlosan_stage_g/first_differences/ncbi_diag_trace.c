/* Comparison-only NCBI scan-callback probe; never link this into LOSAT.
 * Pinned c++/src/algo/blast/core/aa_ungapped.c:478-505:
 *   scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
 *   hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
 * Pinned c++/include/algo/blast/core/lookup_wrap.h:50-58 and
 * blast_aalookup.h:99-132 describe the callback holder layout.
 * Pinned c++/include/algo/blast/core/blast_aascan.h:45-49 gives its ABI.
 * Pinned c++/include/algo/blast/core/blast_parameters.h:95-114
 * gives the per-query-context initial word cutoff fields.
 * Pinned c++/include/algo/blast/core/blast_def.h:141-150 gives pair layout.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { uint32_t q_off, s_off; } OffsetPair;
typedef struct { int lut_type; void* lut; } LookupWrapPrefix;
typedef struct {
    int32_t threshold, mask, charsize, word_length, lut_word_length;
    int32_t alphabet_size, backbone_size, longest_chain;
    int32_t **thin_backbone;
    int32_t bone_type;
    void *thick_backbone, *overflow;
    int32_t overflow_size;
    void *pv;
    int32_t use_pssm;
    void *scansub_callback;
} AaLookupPrefix;
typedef struct { uint8_t* sequence; uint8_t* sequence_start; int32_t length; int16_t frame; } SeqProbe;
typedef struct {
    int32_t last_hit : 31;
    uint32_t flag : 1;
} DiagStructProbe;
typedef struct {
    DiagStructProbe* array;
    uint8_t* hit_len_array;
    int32_t array_length, mask, offset, window, multiple_hits, actual_window;
} DiagTableProbe;
typedef struct { DiagTableProbe* diag; void* hash; } ExtendWordProbe;
typedef struct {
    int32_t x_dropoff_init, x_dropoff, cutoff_score, reduced_nucl_cutoff_score;
} WordCutoffs;
typedef struct {
    void* options;
    int32_t x_dropoff_max, cutoff_score_min;
    WordCutoffs* cutoffs;
} WordParamsPrefix;
typedef int32_t (*scan_fn)(const void*, const void*, OffsetPair*, int32_t, int32_t*);
typedef short (*wordfinder_fn)(void*, void*, void*, void*, void*, void*,
        void*, void*, int32_t, void*, void*);

static scan_fn original_scan;
static unsigned long current_call;
static unsigned long call_index;
static unsigned long ordinal;
static DiagTableProbe* current_diag;

static int32_t trace_scan(const void* lookup, const void* subject,
        OffsetPair* pairs, int32_t capacity, int32_t* range)
{
    if (current_call == 3 && ordinal > 14800 && ordinal < 14950) fprintf(stderr, "SCAN_BEFORE\t%lu\t%d\t%d\t%d\n", ordinal, current_diag->offset, current_diag->array[52].last_hit, current_diag->array[52].flag);
    int32_t count = original_scan(lookup, subject, pairs, capacity, range);
    for (int32_t i = 0; i < count; ++i) {
        fprintf(stderr, "CAND\t%lu\t%lu\t%u\t%u\n",
                current_call, ordinal++, pairs[i].q_off, pairs[i].s_off);
    }
    return count;
}

static void direct_probe(const SeqProbe* query, const SeqProbe* subject, int32_t** matrix) {
    const uint8_t *q = query->sequence, *s = subject->sequence;
    int qright=63, sright=76811, sleft=76781, dropoff=15;
    int score=0, leftscore=0, rightd=0;
    for (int i=0; i<3; i++) { score += matrix[q[qright+i]][s[sright+i]]; if(score>leftscore) {leftscore=score;rightd=i+1;} }
    qright+=rightd; sright+=rightd;
    int n=qright-1 < sright-1 ? qright-1 : sright-1;
    int best_i=n+1, run=0; leftscore=0;
    for (int i=n; i>=0; i--) {
        run += matrix[q[qright-1-n+i]][s[sright-1-n+i]];
        if (run>leftscore) {leftscore=run;best_i=i;}
        if (leftscore-run >= dropoff) break;
    }
    int leftd=n-best_i+1, rightscore=0, rightlen=0, slast=sright, reached=leftd >= sright-sleft;
    if (reached) {
        int best=-1; run=leftscore;rightscore=leftscore;
        int limit=query->length-qright < subject->length-sright ? query->length-qright : subject->length-sright;
        int i=0;
        for (; i<limit; i++) {
            run += matrix[q[qright+i]][s[sright+i]];
            if(run>rightscore){rightscore=run;best=i;}
            if(run<=0 || rightscore-run>=dropoff)break;
        }
        rightlen=best+1;slast=sright+i;
    }
    fprintf(stderr,"DIRECT_EXT q63/s76811 score=%d left=%d right=%d leftd=%d rightlen=%d qstart=%d sstart=%d slast=%d reached=%d\n",leftscore>rightscore?leftscore:rightscore,leftscore,rightscore,leftd,rightlen,qright-leftd,sright-leftd,slast,reached);
}

short BlastAaWordFinder(void* subject, void* query, void* query_info,
        void* lookup, void* matrix, void* word_params, void* ewp,
        void* offset_pairs, int32_t offset_array_size,
        void* init_hitlist, void* ungapped_stats)
{
    wordfinder_fn real = (wordfinder_fn)dlsym(RTLD_NEXT, "BlastAaWordFinder");
    if (!real) return -1;
    LookupWrapPrefix* wrap = (LookupWrapPrefix*)lookup;
    if (wrap->lut_type != 3) return -1; /* eAaLookupTable */
    AaLookupPrefix* aa = (AaLookupPrefix*)wrap->lut;
    original_scan = (scan_fn)aa->scansub_callback;
    current_call = call_index++;
    ordinal = 0;
    DiagTableProbe *diag = ((ExtendWordProbe*)ewp)->diag;
    current_diag = diag;
    if (current_call == 3) {
        SeqProbe* sq = (SeqProbe*)query; SeqProbe* ss = (SeqProbe*)subject;
        direct_probe(sq, ss, (int32_t**)matrix);
        fprintf(stderr, "NCBI_QUERY_LEN %d SUBJECT_LEN %d\n", sq->length, ss->length);
        fprintf(stderr, "NCBI_QUERY_HEX "); for (int j=25;j<100;j++) fprintf(stderr, "%02x", sq->sequence[j]); fprintf(stderr, "\n");
        fprintf(stderr, "NCBI_SUBJ_HEX "); for (int j=76770;j<76850;j++) fprintf(stderr, "%02x", ss->sequence[j]); fprintf(stderr, "\n");
    }
    if (current_call == 3) fprintf(stderr, "DIAG_BEFORE\t%d\t%d\t%d\t%d\t%d\t%d\n", diag->array_length, diag->mask, diag->offset, diag->window, diag->array[52].last_hit, diag->array[52].flag);
    WordCutoffs* cutoffs = ((WordParamsPrefix*)word_params)->cutoffs;
    fprintf(stderr, "PARAM\t%lu\t%d\t%d\t%d\n", current_call,
            cutoffs[0].x_dropoff_init, cutoffs[0].x_dropoff,
            cutoffs[0].cutoff_score);
    aa->scansub_callback = (void*)trace_scan;
    short status = real(subject, query, query_info, lookup, matrix, word_params,
            ewp, offset_pairs, 16, init_hitlist, ungapped_stats);
    aa->scansub_callback = (void*)original_scan;
    if (current_call == 3) fprintf(stderr, "DIAG_AFTER\t%d\t%d\t%d\t%d\n", diag->offset, diag->mask, diag->array[52].last_hit, diag->array[52].flag);
    fprintf(stderr, "END\t%lu\t%lu\n", current_call, ordinal);
    return status;
}
