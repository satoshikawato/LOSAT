//! Port of the non-Smith-Waterman protein-only control flow from
//! `composition_adjustment/redo_alignment.{h,c}` used by `blastp` postprocess.
//!
//! Reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c

use anyhow::{bail, Result};
use std::cell::Cell;
use std::cmp::Ordering;
use std::ffi::c_void;
use std::ptr::NonNull;

use crate::common::GapEditOp;
use crate::utils::matrix::BLASTAA_SIZE;

use super::adjust_scores::{
    blast_adjust_scores, read_aa_composition, AdjustedProteinMatrix, BlastAminoAcidComposition,
    BlastMatrixInfo,
};

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:41-46
// ```c
// /** The natural log of 2 ... */
// #define LOCAL_LN2 0.69314718055994530941723212145818
// ```
const LOCAL_LN2: f64 = 0.6931471805599453;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:791-797
// ```c
// #define KAPPA_BIT_TOL 2.0
// ```
const KAPPA_BIT_TOL: f64 = 2.0;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:79-84
// ```c
// static const int kWindowBorder = 200;
// static const int kReMatrixAdjustmentPseudocounts = 20;
// ```
pub const WINDOW_BORDER: i32 = 200;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1050
// ```c
// #define MINIMUM_LENGTH_NEAR_IDENTICAL 50
// ```
pub const MINIMUM_LENGTH_NEAR_IDENTICAL: i32 = 50;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1269-1271
// ```c
// const double kMinFractionNearIdentical = 0.95;
// int max_shift = 8;
// ```
const MIN_FRACTION_NEAR_IDENTICAL: f64 = 0.95;
const NEAR_IDENTICAL_MAX_SHIFT: usize = 8;
const NEAR_IDENTICAL_WORD_SIZE: usize = 8;
const NEAR_IDENTICAL_HASH_MASK: u64 = 0xFFFFFFFFFF;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1124-1128
// ```c
// ECompoAdjustModes compo_adjust_mode = params->compo_adjust_mode;
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastCompoAdjustMode {
    NoCompositionBasedStats,
    CompositionBasedStats,
    CompositionMatrixAdjust,
    ForceFullMatrixAdjust,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_constants.h:61-68
// ```c
// typedef enum EMatrixAdjustRule {
//     eDontAdjustMatrix              = (-1),
//     eCompoScaleOldMatrix           = 0,
//     eUnconstrainedRelEntropy       = 1,
//     eRelEntropyOldMatrixNewContext = 2,
//     eRelEntropyOldMatrixOldContext = 3,
//     eUserSpecifiedRelEntropy       = 4
// } EMatrixAdjustRule;
// ```
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EMatrixAdjustRule {
    DontAdjustMatrix = -1,
    CompoScaleOldMatrix = 0,
    UnconstrainedRelEntropy = 1,
    RelEntropyOldMatrixNewContext = 2,
    RelEntropyOldMatrixOldContext = 3,
    UserSpecifiedRelEntropy = 4,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-67
// ```c
// typedef struct BlastCompo_Alignment {
//     int score;
//     EMatrixAdjustRule matrix_adjust_rule;
//     int queryIndex;
//     int queryStart;
//     int queryEnd;
//     int matchStart;
//     int matchEnd;
//     int frame;
//     void * context;
//     struct BlastCompo_Alignment * next;
// } BlastCompo_Alignment;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:316-323
// ```c
// GapEditScript * editScript = align->context;
// ...
// Blast_HSPInit(..., align->score, &editScript, &new_hsp);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1921-1926
// ```c
// BlastHSP * hsp = in_align->context;
// q_start = hsp->query.gapped_start - query_range->begin;
// s_start = hsp->subject.gapped_start - subject_range->begin;
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BlastCompoAlignmentContext {
    PreliminaryHspIndex(usize),
    EditScript(Vec<GapEditOp>),
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-67
// ```c
// typedef struct BlastCompo_Alignment {
//     ...
//     void * context;
//     struct BlastCompo_Alignment * next;
// } BlastCompo_Alignment;
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastCompoAlignment {
    pub score: i32,
    pub matrix_adjust_rule: EMatrixAdjustRule,
    pub query_index: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub match_start: i32,
    pub match_end: i32,
    pub frame: i32,
    pub context: Option<BlastCompoAlignmentContext>,
    pub next: Option<Box<BlastCompoAlignment>>,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-106
// ```c
// typedef struct BlastCompo_GappingParams {
//     int gap_open;
//     int gap_extend;
//     int decline_align;
//     int x_dropoff;
//     void * context;
// } BlastCompo_GappingParams;
// ```
#[derive(Debug)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub decline_align: i32,
    pub x_dropoff: i32,
    pub context: Cell<Option<NonNull<c_void>>>,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:111-118
// ```c
// typedef struct BlastCompo_SequenceRange
// {
//     int begin;
//     int end;
//     int context;
// } BlastCompo_SequenceRange;
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastCompoSequenceRange {
    pub begin: i32,
    pub end: i32,
    pub context: i32,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:123-130
// ```c
// typedef struct BlastCompo_SequenceData {
//     Uint1 * data;
//     int length;
//     Uint1 * buffer;
// } BlastCompo_SequenceData;
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastCompoSequenceData {
    pub buffer: Vec<u8>,
    pub data_offset: usize,
    pub length: i32,
}

impl BlastCompoSequenceData {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence_start = (Uint1 *) buffer;
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    pub fn from_ncbistdaa(seq: &[u8]) -> Self {
        let mut buffer = Vec::with_capacity(seq.len() + 2);
        buffer.push(0);
        buffer.extend_from_slice(seq);
        buffer.push(0);
        Self {
            buffer,
            data_offset: 1,
            length: seq.len() as i32,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
    // ```c
    // const Uint1* query, const Uint1* subject
    // ```
    #[inline]
    pub fn data(&self) -> &[u8] {
        let start = self.data_offset;
        let end = start + self.length as usize;
        &self.buffer[start..end]
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1619-1646
    // ```c
    // queryData->length = q_range->end - q_range->begin;
    // ...
    // for (idx = 0;  idx < queryData->length;  idx++) {
    //     queryData->data[idx] = origData[idx];
    // }
    // ```
    pub fn slice_range(&self, range: &BlastCompoSequenceRange) -> Self {
        let begin = range.begin.max(0) as usize;
        let end = range.end.max(range.begin) as usize;
        Self::from_ncbistdaa(&self.data()[begin..end])
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1677-1704
    // ```c
    // origData = query->data + q_range->begin;
    // queryData->length = q_range->end - q_range->begin;
    // queryData->buffer = calloc((queryData->length + 2), sizeof(Uint1));
    // queryData->data   = queryData->buffer + 1;
    // for (idx = 0;  idx < queryData->length;  idx++) {
    //     queryData->data[idx] = (origData[idx] != 24) ? origData[idx] : 3;
    // }
    // ```
    pub fn copy_query_range_with_selenocysteine_fix(
        orig_query: &BlastCompoSequenceData,
        query_range: &BlastCompoSequenceRange,
    ) -> Self {
        let begin = usize::try_from(query_range.begin)
            .expect("NCBI blast_kappa.c requires non-negative query_range.begin");
        let end = usize::try_from(query_range.end)
            .expect("NCBI blast_kappa.c requires non-negative query_range.end");
        let query_len = end.saturating_sub(begin);
        let mut buffer = vec![0u8; query_len + 2];
        for (dst, &src) in buffer[1..1 + query_len]
            .iter_mut()
            .zip(orig_query.data()[begin..end].iter())
        {
            *dst = if src == 24 { 3 } else { src };
        }
        Self {
            buffer,
            data_offset: 1,
            length: i32::try_from(query_len)
                .expect("NCBI blast_kappa.c query range length must fit in Int4"),
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1638-1646
    // ```c
    // seqData->data    = &seqData->data[range->begin - 1];
    // *seqData->data++ = '\0';
    // seqData->length  = range->end - range->begin;
    // ```
    pub fn fit_protein_range_in_place(&mut self, range: &BlastCompoSequenceRange) {
        let begin = usize::try_from(range.begin)
            .expect("NCBI blast_kappa.c requires non-negative subject_range.begin");
        let end = usize::try_from(range.end)
            .expect("NCBI blast_kappa.c requires non-negative subject_range.end");
        let prefix_index = if begin == 0 {
            self.data_offset
                .checked_sub(1)
                .expect("NCBI blast_kappa.c requires a leading sentinel")
        } else {
            self.data_offset + begin - 1
        };
        self.buffer[prefix_index] = 0;
        self.data_offset = prefix_index + 1;
        self.length = i32::try_from(end.saturating_sub(begin))
            .expect("NCBI blast_kappa.c subject range length must fit in Int4");
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:143-146
// ```c
// typedef struct BlastCompo_MatchingSequence {
//   Int4 length;
//   Int4 index;
//   void * local_data;
// } BlastCompo_MatchingSequence;
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastCompoMatchingSequence<'a> {
    pub length: i32,
    pub index: i32,
    pub data: &'a [u8],
}

impl<'a> BlastCompoMatchingSequence<'a> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1067-1106
    // ```c
    // self->length = BlastSeqSrcGetSeqLen(...);
    // self->index = subject_index;
    // ```
    pub fn new(index: i32, data: &'a [u8]) -> Self {
        Self {
            length: data.len() as i32,
            index,
            data,
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:166-178
// ```c
// typedef struct BlastCompo_QueryInfo {
//     int origin;
//     BlastCompo_SequenceData seq;
//     ...
//     double eff_search_space;
//     Uint8* words;
// } BlastCompo_QueryInfo;
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: BlastAminoAcidComposition,
    pub eff_search_space: f64,
    pub words: Option<Vec<u64>>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:85-97
// ```c
// typedef struct s_WindowInfo
// {
//     BlastCompo_SequenceRange query_range;
//     BlastCompo_SequenceRange subject_range;
//     BlastCompo_Alignment * align;
//     int hspcnt;
// } s_WindowInfo;
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowInfo {
    pub query_range: BlastCompoSequenceRange,
    pub subject_range: BlastCompoSequenceRange,
    pub align: Option<Box<BlastCompoAlignment>>,
    pub hspcnt: i32,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:328-347
// ```c
// typedef struct Blast_RedoAlignParams {
//     Blast_MatrixInfo * matrix_info;
//     BlastCompo_GappingParams * gapping_params;
//     ...
//     double near_identical_cutoff;
// } Blast_RedoAlignParams;
// ```
#[derive(Debug)]
pub struct BlastRedoAlignParams {
    pub matrix_info: BlastMatrixInfo,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: BlastCompoAdjustMode,
    pub alphsize: i32,
    pub composition_test_index: i32,
    pub unified_p: bool,
    pub log_k: f64,
    pub score_divisor: f64,
    pub restricted_alignment: bool,
    pub smith_waterman: bool,
    pub is_same_adjustment: bool,
    pub near_identical_cutoff: f64,
    pub position_based: bool,
    pub re_matrix_adjustment_pseudocounts: i32,
    pub ccat_query_length: i32,
    pub query_is_translated: bool,
    pub subject_is_translated: bool,
    pub cutoff_score: i32,
    pub cutoff_evalue: f64,
    pub do_link_hsps: bool,
}

impl BlastRedoAlignParams {
    #[inline]
    pub fn uses_composition_based_stats(&self) -> bool {
        self.compo_adjust_mode != BlastCompoAdjustMode::NoCompositionBasedStats
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:214-236
// ```c
// get_range_type(..., Boolean shouldTestIdentical, ...);
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastRedoRangeResult {
    pub query: BlastCompoSequenceData,
    pub subject: BlastCompoSequenceData,
    pub subject_maybe_biased: bool,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:312-322
// ```c
// typedef struct Blast_RedoAlignCallbacks {
//     calc_lambda_type * calc_lambda;
//     get_range_type * get_range;
//     redo_one_alignment_type * redo_one_alignment;
//     new_xdrop_align_type * new_xdrop_align;
//     free_align_traceback_type * free_align_traceback;
// } Blast_RedoAlignCallbacks;
// ```
pub type BlastCalcLambdaFn = fn(&[f64], i32, i32, f64) -> Result<f64>;
pub type BlastGetRangeFn = for<'seq> fn(
    &BlastCompoMatchingSequence<'seq>,
    &BlastCompoSequenceRange,
    &BlastCompoSequenceData,
    &BlastCompoSequenceRange,
    Option<&[u64]>,
    &BlastCompoAlignment,
    bool,
    bool,
    &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult>;
pub type BlastRedoOneAlignmentFn = fn(
    &BlastCompoAlignment,
    EMatrixAdjustRule,
    Option<&AdjustedProteinMatrix>,
    &BlastCompoSequenceData,
    &BlastCompoSequenceRange,
    i32,
    &BlastCompoSequenceData,
    &BlastCompoSequenceRange,
    i32,
    &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>>;
pub type BlastNewXdropAlignFn = fn(
    i32,
    i32,
    i32,
    i32,
    i32,
    &BlastCompoSequenceData,
    &BlastCompoSequenceRange,
    i32,
    &BlastCompoSequenceData,
    &BlastCompoSequenceRange,
    i32,
    &BlastRedoAlignParams,
    EMatrixAdjustRule,
) -> Result<Option<Box<BlastCompoAlignment>>>;
pub type BlastFreeAlignTracebackFn = fn(usize);

pub struct BlastRedoAlignCallbacks {
    pub calc_lambda: Option<BlastCalcLambdaFn>,
    pub get_range: BlastGetRangeFn,
    pub redo_one_alignment: BlastRedoOneAlignmentFn,
    pub new_xdrop_align: Option<BlastNewXdropAlignFn>,
    pub free_align_traceback: Option<BlastFreeAlignTracebackFn>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1292
// ```c
// int
// Blast_RedoOneMatch(...,
//                    double *pvalueForThisPair,
//                    int compositionTestIndex,
//                    double *LambdaRatio)
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRedoOneMatchResult {
    pub alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>>,
    pub pvalue_for_this_pair: Option<f64>,
    pub lambda_ratio: Option<f64>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1110-1121
// ```c
// static Uint8 s_GetHash(const Uint1* data, int word_size)
// {
//     Uint8 hash = 0;
//     int k;
//     for (k=0;k < word_size;k++) {
//         hash <<= 5;
//         hash += (Int8)data[k];
//     }
//     return hash;
// }
// ```
fn get_hash(data: &[u8], word_size: usize) -> u64 {
    let mut hash = 0u64;
    for residue in data.iter().take(word_size) {
        hash <<= 5;
        hash += u64::from(*residue);
    }
    hash
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2235-2272
// ```c
// static int
// s_CreateWordArray(const Uint1* seq_data, Int4 seq_len, Uint8** words)
// {
//     ...
//     query_hashes = (Uint8*)calloc((seq_len - word_size + 1),
//                                   sizeof(Uint8));
//     query_hashes[0] = s_GetHash(&seq_data[0], word_size);
//     for (i = 1; i < seq_len - word_size; i++) {
//         query_hashes[i] = query_hashes[i - 1];
//         query_hashes[i] <<= 5;
//         query_hashes[i] &= mask;
//         query_hashes[i] += (Uint8)seq_data[i + word_size - 1];
//     }
// }
// ```
pub fn build_query_word_hashes(query_seq: &[u8]) -> Vec<u64> {
    if query_seq.len() < NEAR_IDENTICAL_WORD_SIZE {
        return Vec::new();
    }

    let mut hashes = vec![0u64; query_seq.len() - NEAR_IDENTICAL_WORD_SIZE + 1];
    hashes[0] = get_hash(query_seq, NEAR_IDENTICAL_WORD_SIZE);
    for pos in 1..(query_seq.len() - NEAR_IDENTICAL_WORD_SIZE) {
        hashes[pos] = hashes[pos - 1];
        hashes[pos] <<= 5;
        hashes[pos] &= NEAR_IDENTICAL_HASH_MASK;
        hashes[pos] += u64::from(query_seq[pos + NEAR_IDENTICAL_WORD_SIZE - 1]);
    }
    hashes
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:103-118
// ```c
// BlastCompo_AlignmentNew(int score,
//                         EMatrixAdjustRule matrix_adjust_rule,
//                         int queryStart, int queryEnd, int queryIndex,
//                         int matchStart, int matchEnd, int frame,
//                         void * context)
// {
//     ...
//     align->queryIndex = queryIndex;
//     align->queryStart = queryStart;
// }
// ```
pub fn blast_compo_alignment_new(
    score: i32,
    which_rule: EMatrixAdjustRule,
    query_start: i32,
    query_end: i32,
    query_index: i32,
    match_start: i32,
    match_end: i32,
    frame: i32,
    context: Option<BlastCompoAlignmentContext>,
) -> Box<BlastCompoAlignment> {
    Box::new(BlastCompoAlignment {
        score,
        matrix_adjust_rule: which_rule,
        query_index,
        query_start,
        query_end,
        match_start,
        match_end,
        frame,
        context,
        next: None,
    })
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:88-95
// ```c
// void BlastCompo_AlignmentsFree(BlastCompo_Alignment ** palign,
//                                void (*free_context)(void*));
// ```
pub fn blast_compo_alignments_free(palign: &mut Option<Box<BlastCompoAlignment>>) {
    *palign = None;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:160-178
// ```c
// static int
// s_AlignmentCmp(const BlastCompo_Alignment * a,
//                const BlastCompo_Alignment * b)
// {
//     if (0 == (result = CMP(b->score, a->score)) &&
//         0 == (result = CMP(a->matchStart, b->matchStart)) &&
//         0 == (result = CMP(b->matchEnd, a->matchEnd)) &&
//         0 == (result = CMP(a->queryStart, b->queryStart))) {
//         result = CMP(b->queryEnd, a->queryEnd);
//     }
//     return result;
// }
// ```
pub fn blast_compo_alignment_cmp(a: &BlastCompoAlignment, b: &BlastCompoAlignment) -> Ordering {
    match b.score.cmp(&a.score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a.match_start.cmp(&b.match_start) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match b.match_end.cmp(&a.match_end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a.query_start.cmp(&b.query_start) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b.query_end.cmp(&a.query_end)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:200-208
// ```c
// static int
// s_DistinctAlignmentsLength(BlastCompo_Alignment * list)
// {
//     int length = 0;
//     for ( ;  list != NULL;  list = list->next) {
//         length++;
//     }
//     return length;
// }
// ```
pub fn distinct_alignments_length(list: &Option<Box<BlastCompoAlignment>>) -> usize {
    let mut length = 0usize;
    let mut current = list.as_deref();
    while let Some(node) = current {
        length += 1;
        current = node.next.as_deref();
    }
    length
}

#[cfg(test)]
fn list_to_vec(mut list: Option<Box<BlastCompoAlignment>>) -> Vec<BlastCompoAlignment> {
    let mut out = Vec::new();
    while let Some(mut node) = list {
        list = node.next.take();
        out.push(*node);
    }
    out
}

#[cfg(test)]
fn vec_to_list(mut values: Vec<BlastCompoAlignment>) -> Option<Box<BlastCompoAlignment>> {
    let mut head = None;
    while let Some(mut value) = values.pop() {
        value.next = head;
        head = Some(Box::new(value));
    }
    head
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:214-280
// ```c
// static void
// s_DistinctAlignmentsSort(BlastCompo_Alignment ** plist, int hspcnt)
// {
//     ...
// }
// ```
pub fn distinct_alignments_sort(plist: &mut Option<Box<BlastCompoAlignment>>) {
    let hspcnt = distinct_alignments_length(plist);
    distinct_alignments_sort_internal(plist, hspcnt);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:214-280
// ```c
// static void
// s_DistinctAlignmentsSort(BlastCompo_Alignment ** plist, int hspcnt)
// {
//     ...
//     leftcnt = hspcnt/2;
//     rightcnt = hspcnt - leftcnt;
//     ...
//     while (leftlist != NULL || rightlist != NULL) {
//         ...
//     }
// }
// ```
fn distinct_alignments_sort_internal(plist: &mut Option<Box<BlastCompoAlignment>>, hspcnt: usize) {
    if hspcnt <= 1 {
        return;
    }

    let leftcnt = hspcnt / 2;
    let rightcnt = hspcnt - leftcnt;
    let mut leftlist = plist.take();
    let mut split = leftlist
        .as_mut()
        .expect("NCBI redo_alignment.c requires non-empty left list when hspcnt > 1");
    for _ in 0..leftcnt.saturating_sub(1) {
        split = split
            .next
            .as_mut()
            .expect("NCBI redo_alignment.c split point must exist inside the left half");
    }
    let mut rightlist = split.next.take();

    if leftcnt > 1 {
        distinct_alignments_sort_internal(&mut leftlist, leftcnt);
    }
    if rightcnt > 1 {
        distinct_alignments_sort_internal(&mut rightlist, rightcnt);
    }

    let mut merged = None;
    let mut tail = &mut merged;
    while leftlist.is_some() || rightlist.is_some() {
        if leftlist.is_none() {
            *tail = rightlist.take();
            break;
        } else if rightlist.is_none() {
            *tail = leftlist.take();
            break;
        } else {
            let take_left = {
                let left = leftlist
                    .as_ref()
                    .expect("left list is present when merging distinct alignments");
                let right = rightlist
                    .as_ref()
                    .expect("right list is present when merging distinct alignments");
                blast_compo_alignment_cmp(left, right) == Ordering::Less
            };
            let mut elt = if take_left {
                let mut next = leftlist
                    .take()
                    .expect("left list node exists when chosen for merge");
                leftlist = next.next.take();
                next
            } else {
                let mut next = rightlist
                    .take()
                    .expect("right list node exists when chosen for merge");
                rightlist = next.next.take();
                next
            };
            elt.next = None;
            *tail = Some(elt);
            tail = &mut tail
                .as_mut()
                .expect("merged list tail exists after append")
                .next;
        }
    }
    *plist = merged;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:138-152
// ```c
// static void
// s_AlignmentsRev(BlastCompo_Alignment ** plist)
// {
//     ...
// }
// ```
fn alignments_rev(plist: &mut Option<Box<BlastCompoAlignment>>) {
    let mut list = plist.take();
    let mut reversed = None;
    while let Some(mut node) = list {
        list = node.next.take();
        node.next = reversed;
        reversed = Some(node);
    }
    *plist = reversed;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:282-294
// ```c
// static BlastCompo_Alignment *
// s_AlignmentCopy(const BlastCompo_Alignment * align)
// {
//     return BlastCompo_AlignmentNew(..., align->context);
// }
// ```
fn alignment_copy(align: &BlastCompoAlignment) -> Box<BlastCompoAlignment> {
    blast_compo_alignment_new(
        align.score,
        align.matrix_adjust_rule,
        align.query_start,
        align.query_end,
        align.query_index,
        align.match_start,
        align.match_end,
        align.frame,
        align.context.clone(),
    )
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:303-316
// ```c
// static Boolean
// s_IsSameEndPoint(const BlastCompo_Alignment* newAlign,
//                  const BlastCompo_Alignment* align)
// {
//     ASSERT(newAlign->frame == align->frame);
//     return ((align->queryStart == newAlign->queryStart &&
//              align->matchStart == newAlign->matchStart)
//             || (align->queryEnd == newAlign->queryEnd &&
//                 align->matchEnd == newAlign->matchEnd));
// }
// ```
fn is_same_endpoint(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    new_align.frame == align.frame
        && ((align.query_start == new_align.query_start
            && align.match_start == new_align.match_start)
            || (align.query_end == new_align.query_end && align.match_end == new_align.match_end))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:322-351
// ```c
// static Boolean
// s_IsSimilarEndPoint(const BlastCompo_Alignment* newAlign,
//                     const BlastCompo_Alignment* align)
// {
//     Boolean start_contained = KAPPA_CONTAINED_IN_HSP(...);
//     Boolean end_contained = KAPPA_CONTAINED_IN_HSP(...);
//     Boolean result = (start_contained &&
//                       (newAlign->queryStart - newAlign->matchStart ==
//                        align->queryStart - align->matchStart)) ||
//                      (end_contained &&
//                       (newAlign->queryEnd - newAlign->matchEnd ==
//                        align->queryEnd - align->matchEnd));
//     return result;
// }
// ```
fn is_similar_endpoint(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    let start_contained = align.query_start <= new_align.query_start
        && align.query_end >= new_align.query_start
        && align.match_start <= new_align.match_start
        && align.match_end >= new_align.match_start;
    let end_contained = align.query_start <= new_align.query_end
        && align.query_end >= new_align.query_end
        && align.match_start <= new_align.match_end
        && align.match_end >= new_align.match_end;
    (start_contained
        && (new_align.query_start - new_align.match_start == align.query_start - align.match_start))
        || (end_contained
            && (new_align.query_end - new_align.match_end == align.query_end - align.match_end))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:386-451
// ```c
// static void
// s_WithDistinctEnds(BlastCompo_Alignment **p_newAlign,
//                    BlastCompo_Alignment **p_oldAlignments,
//                    void free_align_context(void *),
//                    Boolean is_same_adjustment)
// {
//     ...
// }
// ```
pub fn with_distinct_ends(
    p_new_align: &mut Option<Box<BlastCompoAlignment>>,
    p_old_alignments: &mut Option<Box<BlastCompoAlignment>>,
    is_same_adjustment: bool,
) {
    let Some(mut new_align) = p_new_align.take() else {
        return;
    };
    let mut include_new_align = true;

    let mut current = p_old_alignments.as_deref();
    while let Some(align) = current {
        let same_end = align.frame == new_align.frame
            && if is_same_adjustment {
                is_same_endpoint(&new_align, align)
            } else {
                is_similar_endpoint(&new_align, align)
            };
        if same_end && new_align.score <= align.score {
            include_new_align = false;
            break;
        }
        current = align.next.as_deref();
    }

    if !include_new_align {
        return;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:420-441
    // ```c
    // tail  = &newAlign->next;
    // align = oldAlignments;
    // while (align != NULL) {
    //     BlastCompo_Alignment * align_next = align->next;
    //     align->next = NULL;
    //     if (...) {
    //         BlastCompo_AlignmentsFree(&align, free_align_context);
    //     } else {
    //         *tail =  align;
    //         tail  = &align->next;
    //     }
    //     align = align_next;
    // }
    // *p_oldAlignments = newAlign;
    // ```
    let mut tail = &mut new_align.next;
    let mut align = p_old_alignments.take();
    while let Some(mut current_align) = align {
        let align_next = current_align.next.take();
        let same_exact_end = current_align.frame == new_align.frame
            && ((current_align.query_start == new_align.query_start
                && current_align.match_start == new_align.match_start)
                || (current_align.query_end == new_align.query_end
                    && current_align.match_end == new_align.match_end));
        if !same_exact_end {
            *tail = Some(current_align);
            tail = &mut tail
                .as_mut()
                .expect("alignment appended to distinct-end list")
                .next;
        }
        align = align_next;
    }
    *p_old_alignments = Some(new_align);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:468-489
// ```c
// static s_WindowInfo *
// s_WindowInfoNew(int begin, int end, int context,
//                 int queryOrigin, int queryLength, int query_index,
//                 BlastCompo_Alignment * align)
// {
//     ...
// }
// ```
fn window_info_new(
    begin: i32,
    end: i32,
    context: i32,
    query_origin: i32,
    query_length: i32,
    query_index: i32,
    align: Option<Box<BlastCompoAlignment>>,
) -> WindowInfo {
    let hspcnt = distinct_alignments_length(&align) as i32;
    WindowInfo {
        subject_range: BlastCompoSequenceRange {
            begin,
            end,
            context,
        },
        query_range: BlastCompoSequenceRange {
            begin: query_origin,
            end: query_origin + query_length,
            context: query_index,
        },
        align,
        hspcnt,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:538-564
// ```c
// static void
// s_WindowInfoJoin(s_WindowInfo * win1, s_WindowInfo ** pwin2)
// {
//     ...
//     tail = &win1->align;
//     for (align = win1->align;  align != NULL;  align = align->next) {
//         tail = &align->next;
//     }
//     *tail = win2->align;
// }
// ```
fn window_info_join(win1: &mut WindowInfo, win2: WindowInfo) {
    win1.subject_range.begin = win1.subject_range.begin.min(win2.subject_range.begin);
    win1.subject_range.end = win1.subject_range.end.max(win2.subject_range.end);
    win1.hspcnt += win2.hspcnt;

    let mut tail = &mut win1.align;
    while let Some(node) = tail {
        tail = &mut node.next;
    }
    *tail = win2.align;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:497-505
// ```c
// static void
// s_WindowSwapRange(s_WindowInfo * window)
// {
//     BlastCompo_SequenceRange tmp = window->query_range;
//     window->query_range = window->subject_range;
//     window->subject_range = tmp;
// }
// ```
fn window_swap_range(window: &mut WindowInfo) {
    std::mem::swap(&mut window.query_range, &mut window.subject_range);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:569-593
// ```c
// static int
// s_LocationCompareWindows(const void * vp1, const void *vp2)
// {
//     ...
// }
// ```
fn location_compare_windows(w1: &WindowInfo, w2: &WindowInfo) -> Ordering {
    match w1.query_range.context.cmp(&w2.query_range.context) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.subject_range.context.cmp(&w2.subject_range.context) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.subject_range.begin.cmp(&w2.subject_range.begin) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.subject_range.end.cmp(&w2.subject_range.end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.query_range.begin.cmp(&w2.query_range.begin) {
        Ordering::Equal => {}
        ord => return ord,
    }
    w1.query_range.end.cmp(&w2.query_range.end)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:598-621
// ```c
// static int
// s_SubjectCompareWindows(const void * vp1, const void *vp2)
// {
//     ...
// }
// ```
fn subject_compare_windows(w1: &WindowInfo, w2: &WindowInfo) -> Ordering {
    match w1.subject_range.begin.cmp(&w2.subject_range.begin) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.subject_range.end.cmp(&w2.subject_range.end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.subject_range.context.cmp(&w2.subject_range.context) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.query_range.begin.cmp(&w2.query_range.begin) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match w1.query_range.end.cmp(&w2.query_range.end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    w1.query_range.context.cmp(&w2.query_range.context)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:738-787
// ```c
// static int
// s_WindowsFromProteinAligns(BlastCompo_Alignment * alignments,
//                            BlastCompo_QueryInfo * query_info,
//                            int numQueries,
//                            int sequence_length,
//                            s_WindowInfo ***pwindows,
//                            int * nWindows)
// {
//     ...
// }
// ```
fn windows_from_protein_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_info: &BlastCompoQueryInfo,
    sequence_length: i32,
) -> Result<Vec<WindowInfo>> {
    let hspcnt = distinct_alignments_length(alignments);
    if hspcnt == 0 {
        return Ok(Vec::new());
    }

    let query_count = alignments
        .as_deref()
        .map(|head| {
            let mut max_query_index = head.query_index;
            let mut current = Some(head);
            while let Some(align) = current {
                max_query_index = max_query_index.max(align.query_index);
                current = align.next.as_deref();
            }
            max_query_index.saturating_add(1)
        })
        .unwrap_or_default();

    let mut windows: Vec<Option<WindowInfo>> = vec![None; query_count as usize];
    let mut align = alignments.as_deref();
    while let Some(current) = align {
        let query_index = usize::try_from(current.query_index).unwrap_or_default();
        let query_length = query_info.seq.length;

        if windows[query_index].is_none() {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:800-806
            // ```c
            // if (windows[query_index] == NULL) {
            //     windows[query_index] =
            //         s_WindowInfoNew(0, sequence_length, 0,
            //                         0, query_length, query_index, NULL);
            // }
            // ```
            windows[query_index] = Some(window_info_new(
                0,
                sequence_length,
                0,
                0,
                query_length,
                current.query_index,
                None,
            ));
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:807-812
        // ```c
        // copiedAlign = s_AlignmentCopy(align);
        // copiedAlign->next = windows[query_index]->align;
        // windows[query_index]->align = copiedAlign;
        // windows[query_index]->hspcnt++;
        // ```
        let mut align_copy = alignment_copy(current);
        if let Some(window) = windows[query_index].as_mut() {
            align_copy.next = window.align.take();
            window.align = Some(align_copy);
            window.hspcnt += 1;
        }

        align = current.next.as_deref();
    }

    let mut packed = Vec::with_capacity(query_count as usize);
    for mut window in windows.into_iter().flatten() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:815-820
        // ```c
        // if (windows[query_index] != NULL) {
        //     windows[window_index] = windows[query_index];
        //     s_AlignmentsRev(&windows[window_index]->align);
        //     window_index++;
        // }
        // ```
        alignments_rev(&mut window.align);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:821-833
        // ```c
        // /* shrink to fit */
        // ...
        // qsort(windows, *nWindows, sizeof(s_WindowInfo*),
        //       s_SubjectCompareWindows);
        // ```
        //
        // Protein-path windows stop at `s_AlignmentsRev`; NCBI does not call
        // `s_DistinctAlignmentsSort` here.
        packed.push(window);
    }

    packed.sort_unstable_by(subject_compare_windows);
    Ok(packed)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:804-833
// ```c
// static Boolean
// s_IsContained(BlastCompo_Alignment * in_align,
//               BlastCompo_Alignment * alignments,
//               double lambda)
// {
//     ...
//     double scoreThresh = score + KAPPA_BIT_TOL * LOCAL_LN2/lambda;
//     ...
// }
// ```
pub fn is_contained(
    in_align: &BlastCompoAlignment,
    alignments: &Option<Box<BlastCompoAlignment>>,
    lambda: f64,
) -> bool {
    let score_thresh = in_align.score as f64 + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;
    let mut current = alignments.as_deref();
    while let Some(align) = current {
        if in_align.frame.signum() == align.frame.signum()
            && align.query_start <= in_align.query_start
            && align.query_end >= in_align.query_start
            && align.match_start <= in_align.match_start
            && align.match_end >= in_align.match_start
            && align.query_start <= in_align.query_end
            && align.query_end >= in_align.query_end
            && align.match_start <= in_align.match_end
            && align.match_end >= in_align.match_end
            && score_thresh <= align.score as f64
        {
            return true;
        }
        current = align.next.as_deref();
    }
    false
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:936-1013
// ```c
// static int s_ExtendRight(Uint1* query_seq, int query_len,
//                          Uint1* subject_seq, int subject_len,
//                          int max_shift,
//                          int* query_ext_len, int* subject_ext_len,
//                          int* align_len)
// {
//     ...
// }
// ```
fn extend_right(
    query_seq: &[u8],
    subject_seq: &[u8],
    max_shift: usize,
) -> (usize, usize, usize, usize) {
    let mut num_identical = 0usize;
    let mut q_pos = 0usize;
    let mut s_pos = 0usize;
    let mut gaps_in_query = 0usize;
    let mut gaps_in_subject = 0usize;

    while q_pos < query_seq.len() && s_pos < subject_seq.len() {
        let mut matched = false;
        while q_pos < query_seq.len()
            && s_pos < subject_seq.len()
            && query_seq[q_pos] == subject_seq[s_pos]
        {
            num_identical += 1;
            q_pos += 1;
            s_pos += 1;
        }

        for n in 1..max_shift {
            if q_pos + n + 1 >= query_seq.len() || s_pos + n + 1 >= subject_seq.len() || matched {
                break;
            }

            if query_seq[q_pos + n] == subject_seq[s_pos + n]
                && query_seq[q_pos + n + 1] == subject_seq[s_pos + n + 1]
            {
                q_pos += n + 2;
                s_pos += n + 2;
                num_identical += 2;
                matched = true;
            }

            if !matched
                && query_seq[q_pos + n] == subject_seq[s_pos]
                && query_seq[q_pos + n + 1] == subject_seq[s_pos + 1]
            {
                q_pos += n + 2;
                s_pos += 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            if !matched
                && query_seq[q_pos] == subject_seq[s_pos + n]
                && query_seq[q_pos + 1] == subject_seq[s_pos + n + 1]
            {
                q_pos += 2;
                s_pos += n + 2;
                num_identical += 2;
                gaps_in_query += n;
                matched = true;
            }
        }

        if !matched {
            break;
        }
    }

    let align_len = if q_pos > s_pos {
        q_pos + gaps_in_query
    } else {
        s_pos + gaps_in_subject
    };
    (num_identical, q_pos, s_pos, align_len)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1025-1103
// ```c
// static int s_ExtendLeft(Uint1* query_seq, int query_len,
//                         Uint1* subject_seq, int subject_len,
//                         int max_shift,
//                         int* query_ext_len, int* subject_ext_len,
//                         int* align_len)
// {
//     ...
// }
// ```
fn extend_left(
    query_seq: &[u8],
    subject_seq: &[u8],
    max_shift: usize,
) -> (usize, usize, usize, usize) {
    if query_seq.is_empty() || subject_seq.is_empty() {
        return (0, 0, 0, 0);
    }

    let mut q_pos = query_seq.len() - 1;
    let mut s_pos = subject_seq.len() - 1;
    let mut num_identical = 0usize;
    let mut gaps_in_query = 0usize;
    let mut gaps_in_subject = 0usize;

    loop {
        let mut matched = false;

        while q_pos > 0 && s_pos > 0 && query_seq[q_pos] == subject_seq[s_pos] {
            num_identical += 1;
            q_pos -= 1;
            s_pos -= 1;
        }

        for n in 1..max_shift {
            if q_pos <= n + 1 || s_pos <= n + 1 || matched {
                break;
            }

            if query_seq[q_pos - n] == subject_seq[s_pos - n]
                && query_seq[q_pos - n - 1] == subject_seq[s_pos - n - 1]
            {
                q_pos -= n + 2;
                s_pos -= n + 2;
                num_identical += 2;
                matched = true;
            }

            if !matched
                && query_seq[q_pos - n] == subject_seq[s_pos]
                && query_seq[q_pos - n - 1] == subject_seq[s_pos - 1]
            {
                q_pos -= n + 2;
                s_pos -= 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            if !matched
                && query_seq[q_pos] == subject_seq[s_pos - n]
                && query_seq[q_pos - 1] == subject_seq[s_pos - n - 1]
            {
                q_pos -= 2;
                s_pos -= n + 2;
                num_identical += 2;
                gaps_in_query += n;
                matched = true;
            }
        }

        if !matched || q_pos == 0 || s_pos == 0 {
            break;
        }
    }

    let query_ext_len = query_seq.len() - q_pos - 1;
    let subject_ext_len = subject_seq.len() - s_pos - 1;
    let align_len = if query_ext_len > subject_ext_len {
        query_ext_len + gaps_in_query
    } else {
        subject_ext_len + gaps_in_subject
    };
    (num_identical, query_ext_len, subject_ext_len, align_len)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1124-1204
// ```c
// static int s_FindNumIdentical(Uint1* query_seq,
//                               const Uint8* query_hashes,
//                               int query_len,
//                               Uint1* subject_seq,
//                               int subject_len,
//                               int max_shift)
// {
//     ...
// }
// ```
fn find_num_identical(
    query_seq: &[u8],
    query_hashes: &[u64],
    subject_seq: &[u8],
    max_shift: usize,
) -> usize {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1160-1164
    // ```c
    // if (!query_seq || !query_hashes || !subject_seq
    //     || query_len < word_size || subject_len < word_size) {
    //     return 0;
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1166-1188
    // ```c
    // for (s_pos = 0; s_pos < subject_len - word_size; s_pos++) {
    //     ...
    //     for (q_pos = query_from; q_pos < query_len - word_size; q_pos++) {
    //         if (query_hashes[q_pos] == hash) {
    //             break;
    //         }
    //     }
    // }
    // ```
    //
    // `s_FindNumIdentical` indexes `query_hashes` only up to
    // `query_len - word_size - 1`, so a query slice of length `L` needs
    // exactly `L - word_size` usable hashes here.
    let required_query_hashes = query_seq.len().saturating_sub(NEAR_IDENTICAL_WORD_SIZE);
    if query_seq.len() < NEAR_IDENTICAL_WORD_SIZE
        || query_hashes.len() < required_query_hashes
        || subject_seq.len() < NEAR_IDENTICAL_WORD_SIZE
    {
        return 0;
    }

    let mut query_from = 0usize;
    let mut subject_from = 0usize;
    let mut num_identical = 0usize;
    let mut matched_last = false;
    let mut hash = 0u64;

    let last_subject_start = subject_seq.len() - NEAR_IDENTICAL_WORD_SIZE;
    let last_query_start = query_seq.len() - NEAR_IDENTICAL_WORD_SIZE;
    let mut s_pos = 0usize;
    while s_pos < last_subject_start {
        if s_pos == 0 || matched_last {
            hash = get_hash(&subject_seq[s_pos..], NEAR_IDENTICAL_WORD_SIZE);
        } else {
            hash <<= 5;
            hash &= NEAR_IDENTICAL_HASH_MASK;
            hash += u64::from(subject_seq[s_pos + NEAR_IDENTICAL_WORD_SIZE - 1]);
        }

        let mut matched_q = None;
        for q_pos in query_from..last_query_start {
            if query_hashes[q_pos] == hash {
                matched_q = Some(q_pos);
                break;
            }
        }

        if let Some(query_start) = matched_q {
            let subject_start = s_pos;
            let (left_identical, _query_left_len, _subject_left_len, _align_len_left) = extend_left(
                &query_seq[query_from..query_start],
                &subject_seq[subject_from..subject_start],
                max_shift,
            );
            let (right_identical, query_right_len, subject_right_len, _align_len_right) =
                extend_right(
                    &query_seq[query_start + NEAR_IDENTICAL_WORD_SIZE..],
                    &subject_seq[subject_start + NEAR_IDENTICAL_WORD_SIZE..],
                    max_shift,
                );

            matched_last = true;
            num_identical += NEAR_IDENTICAL_WORD_SIZE + left_identical + right_identical;
            query_from = query_start + NEAR_IDENTICAL_WORD_SIZE + query_right_len;
            subject_from = subject_start + NEAR_IDENTICAL_WORD_SIZE + subject_right_len;
            s_pos = subject_from.saturating_sub(1);
            if s_pos >= last_subject_start {
                break;
            }
        } else {
            matched_last = false;
        }
        s_pos += 1;
    }

    num_identical
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1259-1312
// ```c
// static Boolean
// s_TestNearIdentical(const BlastCompo_SequenceData* seqData,
//                     const int seqOffset,
//                     const BlastCompo_SequenceData* queryData,
//                     const int queryOffset,
//                     const Uint8* query_words,
//                     const BlastCompo_Alignment* align)
// {
//     ...
// }
// ```
pub fn test_near_identical(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    query_data: &BlastCompoSequenceData,
    query_offset: i32,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
) -> bool {
    let q_start_i32 = align.query_start - query_offset;
    let q_end_i32 = align.query_end - query_offset - 1;
    let s_start_i32 = align.match_start - seq_offset;
    let s_end_i32 = align.match_end - seq_offset - 1;
    if q_start_i32 < 0 || q_end_i32 < q_start_i32 || s_start_i32 < 0 || s_end_i32 < s_start_i32 {
        return false;
    }

    let q_start = q_start_i32 as usize;
    let q_end = q_end_i32 as usize;
    let s_start = s_start_i32 as usize;
    let s_end = s_end_i32 as usize;

    if q_end >= query_data.data().len()
        || s_end >= seq_data.data().len()
        || q_start > q_end
        || s_start > s_end
    {
        return false;
    }

    let query_seq = &query_data.data()[q_start..=q_end];
    let subject_seq = &seq_data.data()[s_start..=s_end];
    let align_len = query_seq.len().min(subject_seq.len());
    if align_len == 0 {
        return false;
    }

    let (mut num_identical, query_right_len, subject_right_len, _) =
        extend_right(query_seq, subject_seq, NEAR_IDENTICAL_MAX_SHIFT);
    if query_right_len >= query_seq.len() || subject_right_len >= subject_seq.len() {
        return (num_identical as f64) / (align_len as f64) > MIN_FRACTION_NEAR_IDENTICAL;
    }

    let (left_identical, query_left_len, subject_left_len, _) = extend_left(
        &query_seq[query_right_len..],
        &subject_seq[subject_right_len..],
        NEAR_IDENTICAL_MAX_SHIFT,
    );
    num_identical += left_identical;
    if query_left_len + query_right_len >= query_seq.len()
        || subject_left_len + subject_right_len >= subject_seq.len()
    {
        return (num_identical as f64) / (align_len as f64) > MIN_FRACTION_NEAR_IDENTICAL;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1322-1333
    // ```c
    // num_identical += s_FindNumIdentical(queryData->data + qStart + query_right_len,
    //                         query_words + qStart + query_right_len,
    //                         query_len - query_left_len - query_right_len,
    //                         seqData->data + sStart + subject_right_len,
    //                         subject_len - subject_left_len - subject_right_len,
    //                         max_shift);
    // fraction_identical = (double)num_identical / (double)align_len;
    // ```
    //
    // If the remaining middle slice is shorter than the 8-mer word size,
    // `s_FindNumIdentical` returns 0 and NCBI still evaluates the final
    // identity fraction. Do not convert a missing or exhausted hash slice
    // into an immediate `false`.
    let query_words = query_words.unwrap_or(&[]);
    let query_word_start = q_start + query_right_len;
    let query_word_slice = if query_word_start <= query_words.len() {
        &query_words[query_word_start..]
    } else {
        &[]
    };

    num_identical += find_num_identical(
        &query_seq[query_right_len..query_seq.len() - query_left_len],
        query_word_slice,
        &subject_seq[subject_right_len..subject_seq.len() - subject_left_len],
        NEAR_IDENTICAL_MAX_SHIFT,
    );

    (num_identical as f64) / (align_len as f64) > MIN_FRACTION_NEAR_IDENTICAL
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1053-1097
// ```c
// static Boolean s_preliminaryTestNearIdentical(..., double cutoff)
// {
//   ...
// }
// ```
pub fn preliminary_test_near_identical(
    query_info: &BlastCompoQueryInfo,
    window: &WindowInfo,
    align: &BlastCompoAlignment,
    cutoff: f64,
) -> bool {
    let query_length = query_info.seq.length;
    if cutoff > 0.0 {
        if (align.match_end - align.match_start + 1)
            < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
        let align_len =
            (align.query_end - align.query_start).min(align.match_end - align.match_start);
        if align_len <= 0 {
            return false;
        }
        return (align.score as f64) / (align_len as f64) >= cutoff;
    }

    if window.hspcnt != 1 {
        return false;
    }
    if (align.query_end - align.query_start) != (align.match_end - align.match_start) {
        return false;
    }
    (align.match_end - align.match_start + 1) >= query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1292
// ```c
// int
// Blast_RedoOneMatch(BlastCompo_Alignment ** alignments,
//                    Blast_RedoAlignParams * params,
//                    BlastCompo_Alignment * incoming_aligns, int hspcnt,
//                    double Lambda,
//                    BlastCompo_MatchingSequence * matchingSeq,
//                    ...
// {
//     ...
// }
// ```
pub fn blast_redo_one_match<'seq>(
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    params: &BlastRedoAlignParams,
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    query_info: &BlastCompoQueryInfo,
    lambda: f64,
    callbacks: &BlastRedoAlignCallbacks,
) -> Result<BlastRedoOneMatchResult> {
    if params.smith_waterman {
        bail!("blastp redo_alignment Smith-Waterman path is not yet ported");
    }
    if params.query_is_translated || params.subject_is_translated {
        bail!("blastp redo_alignment translated paths are not yet ported");
    }
    if params.position_based {
        bail!("blastp redo_alignment position-based path is not yet ported");
    }
    if params.alphsize as usize != BLASTAA_SIZE {
        bail!(
            "blastp redo_alignment alphsize={} is not yet ported beyond BLASTAA_SIZE={}",
            params.alphsize,
            BLASTAA_SIZE
        );
    }

    let windows = windows_from_protein_aligns(incoming_aligns, query_info, matching_seq.length)?;
    let mut alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>> = Vec::new();
    let mut pvalue_for_this_pair = None;
    let mut lambda_ratio = None;
    let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
    let mut adjusted_matrix = None;

    for window in windows {
        let query_index = usize::try_from(window.query_range.context).unwrap_or_default();
        if alignments_by_query.len() <= query_index {
            alignments_by_query.resize_with(query_index + 1, || None);
        }
        let Some(window_align_head) = window.align.as_deref() else {
            continue;
        };
        let mut range_result: Option<BlastRedoRangeResult> = None;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut hsp_index = 0usize;
        let mut num_adjustments = 0usize;
        let mut in_align = window.align.as_deref();

        while let Some(current) = in_align {
            let near_identical = if hsp_index == 0 || subject_maybe_biased {
                preliminary_test_near_identical(
                    query_info,
                    &window,
                    current,
                    params.near_identical_cutoff,
                )
            } else {
                old_near_identical
            };

            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1188-1208
                // ```c
                // status = callbacks->get_range(
                //         matchingSeq,
                //         &window->subject_range,
                //         &subject,
                //         &query_info[query_index].seq,
                //         &window->query_range,
                //         &query,
                //         query_info[query_index].words,
                //         window->align,
                //         nearIdenticalStatus,
                //         compo_adjust_mode,
                //         FALSE,
                //         &subject_maybe_biased
                // );
                // ```
                let next_range = (callbacks.get_range)(
                    matching_seq,
                    &window.subject_range,
                    &query_info.seq,
                    &window.query_range,
                    query_info.words.as_deref(),
                    window_align_head,
                    near_identical,
                    subject_maybe_biased,
                    params,
                )?;
                subject_maybe_biased = next_range.subject_maybe_biased;
                range_result = Some(next_range);
            }

            if !is_contained(current, &alignments_by_query[query_index], lambda) {
                let mut adjust_search_failed = false;
                let range = range_result
                    .as_ref()
                    .expect("range result initialized before redo_one_alignment");

                // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1219-1254
                // ```c
                // if (compo_adjust_mode != eNoCompositionBasedStats &&
                //     (subject_is_translated || hsp_index == 0
                //      || (nearIdenticalStatus != oldNearIdenticalStatus))) {
                //     Blast_AminoAcidComposition subject_composition;
                //     s_GetComposition(&subject_composition, ...);
                //     adjust_search_failed = Blast_AdjustScores(...);
                //     if (adjust_search_failed < 0) { ... }
                //     num_adjustments++;
                // }
                // ```
                if params.compo_adjust_mode != BlastCompoAdjustMode::NoCompositionBasedStats
                    && (params.subject_is_translated
                        || hsp_index == 0
                        || (near_identical != old_near_identical))
                {
                    let query_composition = if params.query_is_translated {
                        read_aa_composition(range.query.data())
                    } else {
                        query_info.composition.clone()
                    };
                    let subject_composition = read_aa_composition(range.subject.data());
                    if query_composition.num_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1228-1249
                        // ```c
                        // adjust_search_failed =
                        //     Blast_AdjustScores(..., &matrix_adjust_rule,
                        //                        ..., pvalueForThisPair,
                        //                        compositionTestIndex,
                        //                        LambdaRatio);
                        // if (adjust_search_failed < 0) { ... }
                        // ```
                        let adjusted = blast_adjust_scores(
                            &params.matrix_info,
                            &query_composition,
                            range.query.length,
                            &subject_composition,
                            range.subject.length,
                            params.compo_adjust_mode,
                            params.composition_test_index,
                        )?;
                        if let Some(adjusted) = adjusted {
                            matrix_adjust_rule = adjusted.matrix_adjust_rule;
                            adjusted_matrix = Some(adjusted.adjusted_matrix);
                            pvalue_for_this_pair = adjusted.pvalue_for_this_pair;
                            lambda_ratio = Some(adjusted.lambda_ratio);
                        } else {
                            adjust_search_failed = true;
                        }
                    }
                    num_adjustments += 1;
                }

                if !adjust_search_failed {
                    let mut new_align = (callbacks.redo_one_alignment)(
                        current,
                        matrix_adjust_rule,
                        adjusted_matrix.as_ref(),
                        &range.query,
                        &window.query_range,
                        params.ccat_query_length,
                        &range.subject,
                        &window.subject_range,
                        matching_seq.length,
                        params,
                    )?;
                    if new_align
                        .as_ref()
                        .is_some_and(|align| align.score >= params.cutoff_score)
                    {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1268-1269
                        // ```c
                        // s_WithDistinctEnds(&newAlign, &alignments[query_index],
                        //                    ..., num_adjustments == 1);
                        // ```
                        with_distinct_ends(
                            &mut new_align,
                            &mut alignments_by_query[query_index],
                            num_adjustments == 1,
                        );
                    }
                }
            }

            old_near_identical = near_identical;
            in_align = current.next.as_deref();
            hsp_index += 1;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1272-1274
    // ```c
    // s_WithDistinctEnds(&newAlign, &alignments[query_index],
    //                    free_align_context, num_adjustments == 1);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1277-1292
    // ```c
    // function_level_cleanup:
    // ...
    // return status;
    // ```
    Ok(BlastRedoOneMatchResult {
        alignments_by_query,
        pvalue_for_this_pair,
        lambda_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::BlastCompoAlignmentContext::PreliminaryHspIndex;
    use super::*;

    fn make_align(
        score: i32,
        q_start: i32,
        q_end: i32,
        s_start: i32,
        s_end: i32,
        context: usize,
    ) -> Box<BlastCompoAlignment> {
        blast_compo_alignment_new(
            score,
            EMatrixAdjustRule::DontAdjustMatrix,
            q_start,
            q_end,
            0,
            s_start,
            s_end,
            0,
            Some(PreliminaryHspIndex(context)),
        )
    }

    #[test]
    fn test_distinct_alignments_sort_orders_ncbi_style() {
        let mut head = Some(make_align(50, 20, 40, 10, 30, 0));
        if let Some(node) = head.as_mut() {
            node.next = Some(make_align(60, 30, 50, 10, 31, 1));
            if let Some(node2) = node.next.as_mut() {
                node2.next = Some(make_align(60, 25, 45, 9, 31, 2));
            }
        }

        distinct_alignments_sort(&mut head);
        let values = list_to_vec(head);
        assert_eq!(values[0].context, Some(PreliminaryHspIndex(2)));
        assert_eq!(values[1].context, Some(PreliminaryHspIndex(1)));
        assert_eq!(values[2].context, Some(PreliminaryHspIndex(0)));
    }

    #[test]
    fn test_with_distinct_ends_keeps_best_shared_endpoint() {
        let mut old = Some(make_align(80, 10, 30, 100, 120, 0));
        let mut new_align = Some(make_align(70, 10, 32, 100, 122, 1));
        with_distinct_ends(&mut new_align, &mut old, true);
        let values = list_to_vec(old);
        assert_eq!(values.len(), 1);
        assert_eq!(values[0].context, Some(PreliminaryHspIndex(0)));
    }

    #[test]
    fn test_with_distinct_ends_similar_endpoint_requires_same_frame() {
        let mut old = Some(make_align(80, 10, 30, 100, 120, 0));
        old.as_mut().unwrap().frame = 1;

        let mut new_align = Some(make_align(70, 12, 32, 102, 122, 1));
        new_align.as_mut().unwrap().frame = -1;

        with_distinct_ends(&mut new_align, &mut old, false);

        let values = list_to_vec(old);
        assert_eq!(values.len(), 2);
        assert_eq!(values[0].context, Some(PreliminaryHspIndex(1)));
        assert_eq!(values[1].context, Some(PreliminaryHspIndex(0)));
    }

    #[test]
    fn test_is_contained_requires_score_tolerance() {
        let mut existing = Some(make_align(100, 10, 50, 100, 140, 0));
        let mut inner = make_align(90, 20, 40, 110, 130, 1);
        inner.next = None;
        assert!(is_contained(&inner, &existing, 0.267));
        existing.as_mut().unwrap().score = 80;
        assert!(!is_contained(&inner, &existing, 0.267));
    }

    #[test]
    fn test_windows_from_protein_aligns_preserves_incoming_order_after_reverse_only() {
        let mut aligns = Some(make_align(10, 0, 10, 0, 10, 0));
        if let Some(node) = aligns.as_mut() {
            node.next = Some(make_align(20, 5, 15, 20, 30, 1));
        }
        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&[1, 2, 3, 4]),
            composition: BlastAminoAcidComposition::empty(),
            eff_search_space: 1.0,
            words: None,
        };

        let windows = windows_from_protein_aligns(&aligns, &query_info, 100).unwrap();
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].subject_range.begin, 0);
        assert_eq!(windows[0].subject_range.end, 100);
        let mut contexts = Vec::new();
        let mut current = windows[0].align.as_deref();
        while let Some(node) = current {
            contexts.push(node.context.clone());
            current = node.next.as_deref();
        }
        assert_eq!(
            contexts,
            vec![Some(PreliminaryHspIndex(0)), Some(PreliminaryHspIndex(1))]
        );
    }

    #[test]
    fn test_windows_from_protein_aligns_separates_queries_by_index() {
        let mut aligns = Some(make_align(10, 0, 10, 0, 10, 0));
        if let Some(node) = aligns.as_mut() {
            node.next = Some(make_align(20, 5, 15, 20, 30, 1));
            node.next.as_mut().unwrap().query_index = 1;
            node.next.as_mut().unwrap().next = Some(make_align(24, 4, 16, 40, 50, 2));
            node.next
                .as_mut()
                .unwrap()
                .next
                .as_mut()
                .unwrap()
                .query_index = 1;
        }

        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&[1, 2, 3, 4]),
            composition: BlastAminoAcidComposition::empty(),
            eff_search_space: 1.0,
            words: None,
        };

        let windows = windows_from_protein_aligns(&aligns, &query_info, 100).unwrap();
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].query_range.context, 0);
        assert_eq!(windows[1].query_range.context, 1);
        assert_eq!(windows[0].hspcnt, 1);
        assert_eq!(windows[1].hspcnt, 2);

        let mut query1_contexts = Vec::new();
        let mut current = windows[1].align.as_deref();
        while let Some(node) = current {
            query1_contexts.push(node.context.clone());
            current = node.next.as_deref();
        }
        assert_eq!(
            query1_contexts,
            vec![Some(PreliminaryHspIndex(1)), Some(PreliminaryHspIndex(2))]
        );
    }

    #[test]
    fn test_blast_redo_one_match_keeps_distinct_alignments_per_query_index() {
        let mut aligns = Some(make_align(80, 0, 40, 0, 40, 0));
        if let Some(node) = aligns.as_mut() {
            node.query_index = 2;
            node.next = Some(make_align(60, 10, 30, 50, 70, 1));
            node.next.as_mut().unwrap().query_index = 2;
        }

        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 100]),
            composition: BlastAminoAcidComposition::empty(),
            eff_search_space: 1.0,
            words: None,
        };
        let matching_seq = BlastCompoMatchingSequence::new(0, &vec![1u8; 100]);
        let params = BlastRedoAlignParams {
            matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
                crate::config::ScoringMatrix::Blosum62,
                0.3176,
            )
            .unwrap(),
            gapping_params: BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: 0,
                x_dropoff: 10,
                context: Cell::new(None),
            },
            compo_adjust_mode: BlastCompoAdjustMode::NoCompositionBasedStats,
            alphsize: BLASTAA_SIZE as i32,
            composition_test_index: 0,
            unified_p: false,
            log_k: 0.0,
            score_divisor: 1.0,
            restricted_alignment: false,
            smith_waterman: false,
            is_same_adjustment: false,
            near_identical_cutoff: 0.0,
            position_based: false,
            re_matrix_adjustment_pseudocounts: 20,
            ccat_query_length: 100,
            query_is_translated: false,
            subject_is_translated: false,
            cutoff_score: 1,
            cutoff_evalue: 10.0,
            do_link_hsps: false,
        };
        fn test_get_range<'seq>(
            _matching_seq: &BlastCompoMatchingSequence<'seq>,
            subject_range: &BlastCompoSequenceRange,
            orig_query: &BlastCompoSequenceData,
            query_range: &BlastCompoSequenceRange,
            _query_words: Option<&[u64]>,
            _align: &BlastCompoAlignment,
            _should_test_identical: bool,
            _subject_maybe_biased: bool,
            _params: &BlastRedoAlignParams,
        ) -> Result<BlastRedoRangeResult> {
            Ok(BlastRedoRangeResult {
                query: orig_query.slice_range(query_range),
                subject: BlastCompoSequenceData::from_ncbistdaa(&vec![
                    1u8;
                    (subject_range.end - subject_range.begin)
                        as usize
                ]),
                subject_maybe_biased: false,
            })
        }

        fn test_redo_one_alignment(
            incoming_align: &BlastCompoAlignment,
            matrix_adjust_rule: EMatrixAdjustRule,
            _adjusted_matrix: Option<&AdjustedProteinMatrix>,
            _query_data: &BlastCompoSequenceData,
            _query_range: &BlastCompoSequenceRange,
            _ccat_query_length: i32,
            _subject_data: &BlastCompoSequenceData,
            _subject_range: &BlastCompoSequenceRange,
            _full_subject_length: i32,
            _params: &BlastRedoAlignParams,
        ) -> Result<Option<Box<BlastCompoAlignment>>> {
            Ok(Some(blast_compo_alignment_new(
                incoming_align.score,
                matrix_adjust_rule,
                incoming_align.query_start,
                incoming_align.query_end,
                incoming_align.query_index,
                incoming_align.match_start,
                incoming_align.match_end,
                incoming_align.frame,
                incoming_align.context.clone(),
            )))
        }

        let callbacks = BlastRedoAlignCallbacks {
            calc_lambda: None,
            get_range: test_get_range,
            redo_one_alignment: test_redo_one_alignment,
            new_xdrop_align: None,
            free_align_traceback: None,
        };

        let result = blast_redo_one_match(
            &aligns,
            &params,
            &matching_seq,
            &query_info,
            0.267,
            &callbacks,
        )
        .unwrap();

        assert_eq!(result.alignments_by_query.len(), 3);
        assert!(result.alignments_by_query[0].is_none());
        assert!(result.alignments_by_query[1].is_none());
        let mut contexts = Vec::new();
        let mut current = result.alignments_by_query[2].as_deref();
        while let Some(node) = current {
            contexts.push(node.context.clone());
            current = node.next.as_deref();
        }
        assert_eq!(
            contexts,
            vec![Some(PreliminaryHspIndex(1)), Some(PreliminaryHspIndex(0))]
        );
    }

    #[test]
    fn test_preliminary_test_near_identical_uses_cutoff_branch() {
        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 60]),
            composition: BlastAminoAcidComposition::empty(),
            eff_search_space: 1.0,
            words: None,
        };
        let window = window_info_new(0, 60, 0, 0, 60, 0, Some(make_align(240, 0, 59, 0, 59, 0)));
        let align = window.align.as_deref().unwrap();
        assert!(preliminary_test_near_identical(
            &query_info,
            &window,
            align,
            3.5
        ));
        assert!(!preliminary_test_near_identical(
            &query_info,
            &window,
            align,
            5.0
        ));
    }

    #[test]
    fn test_location_compare_windows_sorts_by_query_then_subject() {
        let left = window_info_new(0, 100, 0, 0, 50, 0, None);
        let right = window_info_new(10, 100, 0, 0, 50, 0, None);
        assert_eq!(location_compare_windows(&left, &right), Ordering::Less);
    }

    #[test]
    fn test_build_query_word_hashes_matches_ncbi_shift_hash() {
        let query = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let hashes = build_query_word_hashes(&query);
        assert_eq!(hashes.len(), query.len() - NEAR_IDENTICAL_WORD_SIZE + 1);
        assert_eq!(hashes[0], get_hash(&query, NEAR_IDENTICAL_WORD_SIZE));
        assert_eq!(hashes[1], get_hash(&query[1..], NEAR_IDENTICAL_WORD_SIZE));
        assert_eq!(hashes[2], 0);
    }

    #[test]
    fn test_find_num_identical_accepts_short_hash_slice_for_9aa_window() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1160-1188
        // ```c
        // if (!query_seq || !query_hashes || !subject_seq
        //     || query_len < word_size || subject_len < word_size) {
        //     return 0;
        // }
        // for (q_pos = query_from; q_pos < query_len - word_size; q_pos++) {
        //     if (query_hashes[q_pos] == hash) {
        //         break;
        //     }
        // }
        // ```
        //
        // A 9-aa slice has only two stored 8-mer hashes; NCBI still scans it.
        let query = [1u8, 3, 4, 5, 6, 7, 8, 9, 10];
        let hashes = build_query_word_hashes(&query);
        assert_eq!(hashes.len(), 2);
        assert_eq!(
            find_num_identical(&query, &hashes, &query, NEAR_IDENTICAL_MAX_SHIFT),
            query.len()
        );
    }

    #[test]
    fn test_test_near_identical_accepts_high_identity_alignment() {
        let residues = vec![1u8; 80];
        let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
        let subject = BlastCompoSequenceData::from_ncbistdaa(&residues);
        let align = make_align(80, 0, 80, 0, 80, 0);
        let words = build_query_word_hashes(query.data());
        assert!(test_near_identical(
            &subject,
            0,
            &query,
            0,
            Some(&words),
            &align,
        ));
    }

    #[test]
    fn test_test_near_identical_rejects_low_identity_alignment() {
        let query = BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 80]);
        let subject = BlastCompoSequenceData::from_ncbistdaa(&vec![2u8; 80]);
        let align = make_align(20, 0, 80, 0, 80, 0);
        let words = build_query_word_hashes(query.data());
        assert!(!test_near_identical(
            &subject,
            0,
            &query,
            0,
            Some(&words),
            &align,
        ));
    }

    #[test]
    fn test_test_near_identical_keeps_fraction_check_when_middle_slice_has_no_words() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1314-1333
        // ```c
        // if (query_left_len + query_right_len >= query_len
        //     || subject_left_len + subject_right_len >= subject_len) {
        //     ...
        // }
        // num_identical += s_FindNumIdentical(...);
        // fraction_identical = (double)num_identical / (double)align_len;
        // ```
        //
        // A short unresolved middle slice still flows to the final fraction
        // check; NCBI does not return FALSE just because there are no 8-mer
        // hashes left to inspect.
        let mut query_residues = vec![1u8; 48];
        query_residues.extend_from_slice(&[3, 4, 5, 6]);
        query_residues.extend_from_slice(&vec![1u8; 48]);

        let mut subject_residues = vec![1u8; 48];
        subject_residues.extend_from_slice(&[7, 8, 9, 10]);
        subject_residues.extend_from_slice(&vec![1u8; 48]);

        let query = BlastCompoSequenceData::from_ncbistdaa(&query_residues);
        let subject = BlastCompoSequenceData::from_ncbistdaa(&subject_residues);
        let align = make_align(96, 0, 100, 0, 100, 0);

        assert!(test_near_identical(&subject, 0, &query, 0, None, &align));
    }
}
