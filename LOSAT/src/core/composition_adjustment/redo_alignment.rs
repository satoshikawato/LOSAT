//! Port of the non-Smith-Waterman protein and translated-subject control flow from
//! `composition_adjustment/redo_alignment.{h,c}` used by `blastp` postprocess.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c

use anyhow::{bail, Result};
use std::cell::Cell;
use std::cmp::Ordering;
use std::ffi::c_void;
use std::ptr::NonNull;

use super::x_align_index::{x_add_distinct, x_contains, x_finish, XAlignIndexes};
use crate::common::GapEditOp;
use crate::sequence::packed_subject::SubjectNt;
use crate::stats::karlin_calc::calc_lambda_of_score_probs;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::{aa_char_to_ncbistdaa, BLASTAA_SIZE};
use crate::utils::seg::{SegMasker, SegParams};

use super::adjust_scores::{
    blast_adjust_scores, blast_get_composition_range, read_aa_composition, AdjustedProteinMatrix,
    BlastAminoAcidComposition, BlastCompositionWorkspace, BlastMatrixInfo,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:550-572
/// The redo's `calc_lambda` callback of BLASTP, TBLASTN and BLASTX: NCBI's `s_CalcLambda` as
/// written (`calc_lambda_of_score_probs`: the probabilities of the scores `min_score` to
/// `max_score` as given, neither normalized nor narrowed to the nonzero scores), -1.0 where
/// `Blast_KarlinLambdaNR` gives it, which `Blast_CompositionBasedStats` turns into the lower bound
/// of the ratio (`adjust_scores::compute_lambda_ratio`). Mode 1 calls it for every pair (item SC5
/// of line NP), modes 2 and 3 where the matrix adjustment falls back to scaling.
#[inline]
pub(crate) fn redo_calc_lambda(
    probs: &[f64],
    min_score: i32,
    max_score: i32,
    lambda0: f64,
) -> Result<f64> {
    Ok(calc_lambda_of_score_probs(
        probs.to_vec(),
        min_score,
        max_score,
        lambda0,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:41-46
pub(crate) const LOCAL_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:955
pub(crate) const KAPPA_BIT_TOL: f64 = 2.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:79-84
pub const WINDOW_BORDER: i32 = 200;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1050
pub const MINIMUM_LENGTH_NEAR_IDENTICAL: i32 = 50;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1269-1272
const MIN_FRACTION_NEAR_IDENTICAL: f64 = 0.95;
const NEAR_IDENTICAL_MAX_SHIFT: usize = 8;
const NEAR_IDENTICAL_WORD_SIZE: usize = 8;
const NEAR_IDENTICAL_HASH_MASK: u64 = 0xFFFFFFFFFF;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1124-1128
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastCompoAdjustMode {
    NoCompositionBasedStats,
    CompositionBasedStats,
    CompositionMatrixAdjust,
    ForceFullMatrixAdjust,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_constants.h:59-70
impl BlastCompoAdjustMode {
    /// NCBI's value of the mode (`ECompoAdjustModes`, 0 to 3).
    pub fn ncbi_value(self) -> i32 {
        match self {
            BlastCompoAdjustMode::NoCompositionBasedStats => 0,
            BlastCompoAdjustMode::CompositionBasedStats => 1,
            BlastCompoAdjustMode::CompositionMatrixAdjust => 2,
            BlastCompoAdjustMode::ForceFullMatrixAdjust => 3,
        }
    }

    /// Whether the search runs the composition-based redo (a mode above 0).
    pub fn is_enabled(self) -> bool {
        self != BlastCompoAdjustMode::NoCompositionBasedStats
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_constants.h:76-83
// `Blast_ChooseMatrixAdjustRule` (compo_mode_condition.c) returns only `eDontAdjustMatrix`,
// `eCompoScaleOldMatrix` and `eUserSpecifiedRelEntropy`, so the three other values are not
// declared.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EMatrixAdjustRule {
    DontAdjustMatrix = -1,
    CompoScaleOldMatrix = 0,
    UserSpecifiedRelEntropy = 4,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-71
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:316-323
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1916-1926
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BlastCompoAlignmentContext {
    PreliminaryHspIndex(usize),
    EditScript(Vec<GapEditOp>),
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-71
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

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:97-112
// `decline_align` is not stored: nothing reads it (it is 0 in every caller).
#[derive(Debug)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub x_dropoff: i32,
    pub context: Cell<Option<NonNull<c_void>>>,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:111-125
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastCompoSequenceRange {
    pub begin: i32,
    pub end: i32,
    pub context: i32,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:123-130
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastCompoSequenceData {
    pub buffer: Vec<u8>,
    pub data_offset: usize,
    pub length: i32,
}

impl BlastCompoSequenceData {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
    #[inline]
    pub fn data(&self) -> &[u8] {
        let start = self.data_offset;
        let end = start + self.length as usize;
        &self.buffer[start..end]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1686-1696
    #[cfg(test)]
    pub fn slice_range(&self, range: &BlastCompoSequenceRange) -> Self {
        let begin = range.begin.max(0) as usize;
        let end = range.end.max(range.begin) as usize;
        Self::from_ncbistdaa(&self.data()[begin..end])
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1677-1704
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1638-1646
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

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:143-146
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastCompoMatchingSequence<'a> {
    pub length: i32,
    pub index: i32,
    pub data: &'a [u8],
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1494-1497,1522-1525
    pub genetic_code_id: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1518-1520
    /// The nucleotides a translated subject's windows are read from (`na_sequence`): its
    /// letters (then also `data`), or (EXPERIMENT, LOSAT_X_TBNSTORE) its packed store, whose
    /// windows are made when read. `data` for a protein.
    pub nucleotides: SubjectNt<'a>,
}

impl<'a> BlastCompoMatchingSequence<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1067-1106
    pub fn new(index: i32, data: &'a [u8]) -> Self {
        Self {
            length: data.len() as i32,
            index,
            data,
            genetic_code_id: None,
            nucleotides: SubjectNt::Letters(data),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1501
    pub fn new_translated(
        index: i32,
        nucleotides: SubjectNt<'a>,
        genetic_code_id: u8,
    ) -> Result<Self> {
        GeneticCode::try_from_id(genetic_code_id).map_err(anyhow::Error::msg)?;
        Ok(Self {
            length: i32::try_from(nucleotides.len())?,
            index,
            // the letters, when the subject has them (a store's windows are read through
            // `nucleotides` only)
            data: match nucleotides {
                SubjectNt::Letters(letters) => letters,
                SubjectNt::Packed(_) => &[],
            },
            genetic_code_id: Some(genetic_code_id),
            nucleotides,
        })
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:166-178
#[derive(Debug, Clone, PartialEq)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: BlastAminoAcidComposition,
    pub eff_search_space: f64,
    pub words: Option<Vec<u64>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:85-98
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowInfo {
    pub query_range: BlastCompoSequenceRange,
    pub subject_range: BlastCompoSequenceRange,
    pub align: Option<Box<BlastCompoAlignment>>,
    pub hspcnt: i32,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:328-357
// `unified_p`, `restricted_alignment`, `is_same_adjustment`, `re_matrix_adjustment_pseudocounts`
// and `cutoff_evalue` are not stored: nothing reads them (`composition_test_index` carries what
// the redo takes from `unified_p`).
#[derive(Debug)]
pub struct BlastRedoAlignParams {
    pub matrix_info: BlastMatrixInfo,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: BlastCompoAdjustMode,
    pub alphsize: i32,
    pub composition_test_index: i32,
    pub log_k: f64,
    pub score_divisor: f64,
    pub smith_waterman: bool,
    pub near_identical_cutoff: f64,
    pub position_based: bool,
    pub ccat_query_length: i32,
    pub query_is_translated: bool,
    pub subject_is_translated: bool,
    pub cutoff_score: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1040
    /// `do_link_hsps`: the search links HSPs with sum statistics (`hitParams->do_sum_stats`), so
    /// the redo keeps alignments by `cutoff_score` and BLASTP's evaluation links them
    /// (`blastp::kappa`, psiblast's `-sum_stats true`). BLASTX and TBLASTN link in their own
    /// evaluations.
    pub do_link_hsps: bool,
}

impl BlastRedoAlignParams {
    #[inline]
    pub fn uses_composition_based_stats(&self) -> bool {
        self.compo_adjust_mode != BlastCompoAdjustMode::NoCompositionBasedStats
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:214-236
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastRedoRangeResult {
    pub query: BlastCompoSequenceData,
    pub subject: BlastCompoSequenceData,
    pub subject_maybe_biased: bool,
    /// `read_aa_composition(subject.data())`, computed once
    /// when blastp's SEG pass masked these subject bytes (and its cache keeps them); `None` where
    /// nothing computed it. Shared, so that a cache hit copies a pointer.
    pub x_subject_composition: Option<std::sync::Arc<XSubjectComposition>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:934-936
// No NCBI counterpart: the composition NCBI reads from a subject range, kept with the range's
// bytes so that the same bytes are not read again. Two values are equal when their bits are.
#[derive(Debug, Clone)]
pub struct XSubjectComposition(pub BlastAminoAcidComposition);

impl PartialEq for XSubjectComposition {
    fn eq(&self, other: &Self) -> bool {
        self.0.num_true_amino_acids == other.0.num_true_amino_acids
            && self
                .0
                .prob
                .iter()
                .zip(&other.0.prob)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }
}

impl Eq for XSubjectComposition {}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:934-936
// The composition of a whole, untranslated subject range: the stored one when the range carries
// it (it was read from the same bytes by the same function),
// otherwise read here.
pub(crate) fn x_whole_subject_composition(
    range: &BlastRedoRangeResult,
) -> BlastAminoAcidComposition {
    if let Some(stored) = &range.x_subject_composition {
        return stored.0.clone();
    }
    read_aa_composition(range.subject.data())
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:312-323
// LOSAT has no counterpart of `new_xdrop_align` and `free_align_traceback`: LOSAT's three callers
// (blastp, blastx, tblastn) never set them, so the two fields and their function types are not
// ported (P-V3b).
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
pub struct BlastRedoAlignCallbacks {
    pub calc_lambda: Option<BlastCalcLambdaFn>,
    pub get_range: BlastGetRangeFn,
    pub redo_one_alignment: BlastRedoOneAlignmentFn,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1292
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRedoOneMatchResult {
    pub alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>>,
    pub pvalue_for_this_pair: Option<f64>,
    pub lambda_ratio: Option<f64>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1110-1125
fn get_hash(data: &[u8], word_size: usize) -> u64 {
    let mut hash = 0u64;
    for residue in data.iter().take(word_size) {
        hash <<= 5;
        hash += u64::from(*residue);
    }
    hash
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2235-2272
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:103-118
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:160-180
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:199-208
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:214-280
pub fn distinct_alignments_sort(plist: &mut Option<Box<BlastCompoAlignment>>) {
    let hspcnt = distinct_alignments_length(plist);
    distinct_alignments_sort_internal(plist, hspcnt);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:214-280
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:138-152
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:282-294
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:303-316
fn is_same_endpoint(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    new_align.frame == align.frame
        && ((align.query_start == new_align.query_start
            && align.match_start == new_align.match_start)
            || (align.query_end == new_align.query_end && align.match_end == new_align.match_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:321-351
pub(crate) fn is_similar_endpoint(
    new_align: &BlastCompoAlignment,
    align: &BlastCompoAlignment,
) -> bool {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:385-451
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:420-446
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:468-489
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:537-564
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:497-505
fn window_swap_range(window: &mut WindowInfo) {
    std::mem::swap(&mut window.query_range, &mut window.subject_range);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:569-593
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:598-621
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

/// `s_GetTranslatedLength`: the letters of `frame`'s translation of a sequence of `length`
/// (nucleotides; with a position-based matrix, RPSTBLASTN's translation buffer of the six frames,
/// whose nucleotide length `GET_NUCL_LENGTH` estimates). C's `/` (toward zero) as NCBI's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:624-635
pub(crate) fn translated_length(length: i32, frame: i32, is_pos_based: bool) -> i32 {
    if is_pos_based {
        get_translated_length(get_nucl_length(length), get_seq_frame(frame))
    } else {
        (length - frame.abs() + 1) / 3
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:491-499
/// `GET_SEQ_FRAME`: the index (0 to 5) of a frame.
pub(crate) fn get_seq_frame(frame: i32) -> i32 {
    if frame < 0 {
        frame.abs() + 2
    } else {
        frame - 1
    }
}

/// `GET_NUCL_LENGTH`: the nucleotide length of a translation buffer of `length` bytes.
pub(crate) fn get_nucl_length(length: i32) -> i32 {
    (length - 5) / 2 + 2
}

/// `GET_TRANSLATED_LENGTH`: the letters of the frame of index `f` of `length` nucleotides.
pub(crate) fn get_translated_length(length: i32, f: i32) -> i32 {
    (length - f % 3) / 3
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:80,644-765:
// if (same subject/query context && window->subject_range.end >=
// The sorts are glibc's `qsort`, a stable merge sort (see `blastp::hsp::sort_hsplist_by_score`).
// `is_pos_based`: NCBI's `params->positionBased` (RPSTBLASTN's redo, series PR, item PS12b).
pub(crate) fn windows_from_translated_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_lengths: &[i32],
    subject_nt_length: i32,
    border: i32,
    is_pos_based: bool,
) -> Result<Vec<WindowInfo>> {
    let mut windows = Vec::with_capacity(distinct_alignments_length(alignments));
    let mut current = alignments.as_deref();
    while let Some(align) = current {
        let query_index = usize::try_from(align.query_index)
            .map_err(|_| anyhow::anyhow!("NCBI translated window query index is negative"))?;
        let query_length = *query_lengths
            .get(query_index)
            .ok_or_else(|| anyhow::anyhow!("NCBI translated window query index is missing"))?;
        let translated_length = translated_length(subject_nt_length, align.frame, is_pos_based);
        windows.push(window_info_new(
            (align.match_start - border).max(0),
            (align.match_end + border).min(translated_length),
            align.frame,
            0,
            query_length,
            align.query_index,
            Some(alignment_copy(align)),
        ));
        current = align.next.as_deref();
    }
    windows.sort_by(location_compare_windows);
    let mut input = windows.into_iter().peekable();
    let mut joined = Vec::new();
    while let Some(window) = input.next() {
        if let Some(next) = input.peek_mut() {
            if window.subject_range.context == next.subject_range.context
                && window.query_range.context == next.query_range.context
                && window.subject_range.end >= next.subject_range.begin
            {
                window_info_join(next, window);
                continue;
            }
        }
        joined.push(window);
    }
    for window in &mut joined {
        distinct_alignments_sort(&mut window.align);
    }
    joined.sort_by(subject_compare_windows);
    Ok(joined)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:681-693
pub(crate) fn windows_from_translated_query_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_infos: &[BlastCompoQueryInfo],
    subject_length: i32,
    border: i32,
) -> Result<Vec<WindowInfo>> {
    let mut windows = Vec::with_capacity(distinct_alignments_length(alignments));
    let mut current = alignments.as_deref();
    while let Some(align) = current {
        let query_length = query_infos
            .get(usize::try_from(align.query_index)?)
            .ok_or_else(|| anyhow::anyhow!("BLASTX redo context is missing"))?
            .seq
            .length;
        windows.push(window_info_new(
            (align.query_start - border).max(0),
            (align.query_end + border).min(query_length),
            align.query_index,
            0,
            subject_length,
            0,
            Some(alignment_copy(align)),
        ));
        current = align.next.as_deref();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:700-751
    windows.sort_unstable_by(location_compare_windows);
    let mut input = windows.into_iter().peekable();
    let mut joined = Vec::new();
    while let Some(window) = input.next() {
        if let Some(next) = input.peek_mut() {
            if window.subject_range.context == next.subject_range.context
                && window.query_range.context == next.query_range.context
                && window.subject_range.end >= next.subject_range.begin
            {
                window_info_join(next, window);
                continue;
            }
        }
        joined.push(window);
    }
    for window in &mut joined {
        window_swap_range(window);
        distinct_alignments_sort(&mut window.align);
    }
    joined.sort_unstable_by(subject_compare_windows);
    Ok(joined)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:906-938
pub(crate) fn translated_query_composition(
    query: &BlastCompoSequenceData,
    range: &BlastCompoSequenceRange,
    align: &BlastCompoAlignment,
) -> Result<BlastAminoAcidComposition> {
    let start = usize::try_from(align.query_start - range.begin)?;
    let finish = usize::try_from(align.query_end - range.begin)?;
    let (left, right) = blast_get_composition_range(query.data(), start, finish)?;
    Ok(read_aa_composition(&query.data()[left..right]))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1504-1525;
// c++/src/algo/blast/core/blast_util.c:428-454,1141-1205:
// translation_start = frame > 0 ? 3 * range->begin
fn translated_subject_window_sequence(
    subject_nt: SubjectNt<'_>,
    range: &BlastCompoSequenceRange,
    genetic_code: &GeneticCode,
) -> Result<Vec<u8>> {
    let frame = range.context;
    if frame == 0 || frame.abs() > 3 || range.begin < 0 || range.end < range.begin {
        bail!("invalid NCBI translated subject window");
    }
    let full_length = i32::try_from(subject_nt.len())?;
    let start = if frame > 0 {
        3 * range.begin
    } else {
        full_length - 3 * range.end + frame + 1
    };
    let nucleotide_length = 3 * (range.end - range.begin) + frame.abs() - 1;
    if start < 0 || nucleotide_length < 0 || start + nucleotide_length > full_length {
        bail!("NCBI translated subject window exceeds nucleotide sequence");
    }
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject gives the window's letters in upper case;
    // the strand and codon reads below do not depend on the case.
    let mut window = Vec::new();
    let segment = subject_nt.fetch(
        start as usize,
        (start + nucleotide_length) as usize,
        &mut window,
    );
    let strand = if frame > 0 {
        segment.to_vec()
    } else {
        crate::utils::genetic_code::revcomp(segment)
    };
    let mut translated = Vec::with_capacity((range.end - range.begin) as usize + 2);
    translated.push(0);
    for codon in strand[(frame.abs() - 1) as usize..].chunks_exact(3) {
        translated.push(aa_char_to_ncbistdaa(genetic_code.get(codon)));
    }
    translated.push(0);
    if translated.len() != (range.end - range.begin) as usize + 2 {
        bail!("NCBI translated subject window length differs");
    }
    Ok(translated)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:906-938:
pub(crate) fn translated_subject_composition(
    translation: &[u8],
    range: &BlastCompoSequenceRange,
    align: &BlastCompoAlignment,
) -> Result<BlastAminoAcidComposition> {
    if translation.len() < 2 || translation[0] != 0 || translation[translation.len() - 1] != 0 {
        bail!("NCBI translated subject data lacks sentinels");
    }
    let data = &translation[1..translation.len() - 1];
    let start = usize::try_from(align.match_start - range.begin)?;
    let finish = usize::try_from(align.match_end - range.begin)?;
    let (left, right) = blast_get_composition_range(data, start, finish)?;
    Ok(read_aa_composition(&data[left..right]))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1414-1455,1530-1548:
fn mask_translated_subject_seg(translation: &mut [u8]) -> Result<bool> {
    if translation.len() < 2 || translation[0] != 0 || translation[translation.len() - 1] != 0 {
        bail!("NCBI translated subject SEG data lacks sentinels");
    }
    let data_len = translation.len() - 1;
    let data = &mut translation[1..data_len];
    // LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): the windows of one subject overlap from query to query, so
    // this masker takes the trim of a segment trimmed before (on any thread) from the memo.
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1)).x_with_trim_memo();
    let intervals = masker.mask_sequence(data);
    for interval in &intervals {
        for residue in &mut data[interval.start..interval.end] {
            *residue = 21;
        }
    }
    Ok(!intervals.is_empty())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1475-1552,1677-1704
pub fn translated_subject_get_range<'seq>(
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    subject_maybe_biased: bool,
    params: &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult> {
    if !params.subject_is_translated || params.query_is_translated {
        bail!("TBLASTN translated range requires protein query and translated subject");
    }
    let code_id = matching_seq
        .genetic_code_id
        .ok_or_else(|| anyhow::anyhow!("TBLASTN translated range requires subject genetic code"))?;
    let genetic_code = GeneticCode::try_from_id(code_id).map_err(anyhow::Error::msg)?;
    let query =
        BlastCompoSequenceData::copy_query_range_with_selenocysteine_fix(orig_query, query_range);
    let translated =
        translated_subject_window_sequence(matching_seq.nucleotides, subject_range, &genetic_code)?;
    let mut subject = BlastCompoSequenceData {
        length: i32::try_from(translated.len() - 2)?,
        buffer: translated,
        data_offset: 1,
    };
    let mut maybe_biased = subject_maybe_biased;
    if params.uses_composition_based_stats()
        && maybe_biased
        && (!should_test_identical
            || !test_near_identical(
                &subject,
                subject_range.begin,
                &query,
                query_range.begin,
                query_words,
                align,
            ))
    {
        maybe_biased = mask_translated_subject_seg(&mut subject.buffer)?;
    }
    Ok(BlastRedoRangeResult {
        query,
        subject,
        subject_maybe_biased: maybe_biased,
        x_subject_composition: None,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:738-787
pub(crate) fn windows_from_protein_aligns(
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:800-806
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:807-812
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:815-820
        alignments_rev(&mut window.align);

        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:821-836
        //
        // Protein-path windows stop at `s_AlignmentsRev`; NCBI does not call
        // `s_DistinctAlignmentsSort` here.
        packed.push(window);
    }

    packed.sort_unstable_by(subject_compare_windows);
    Ok(packed)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:965-978
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:936-1013
fn extend_right<Q: WalkLetters>(
    query_seq: Q,
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
            && query_seq.at(q_pos) == subject_seq[s_pos]
        {
            num_identical += 1;
            q_pos += 1;
            s_pos += 1;
        }

        for n in 1..max_shift {
            if q_pos + n + 1 >= query_seq.len() || s_pos + n + 1 >= subject_seq.len() || matched {
                break;
            }

            if query_seq.at(q_pos + n) == subject_seq[s_pos + n]
                && query_seq.at(q_pos + n + 1) == subject_seq[s_pos + n + 1]
            {
                q_pos += n + 2;
                s_pos += n + 2;
                num_identical += 2;
                matched = true;
            }

            if !matched
                && query_seq.at(q_pos + n) == subject_seq[s_pos]
                && query_seq.at(q_pos + n + 1) == subject_seq[s_pos + 1]
            {
                q_pos += n + 2;
                s_pos += 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            if !matched
                && query_seq.at(q_pos) == subject_seq[s_pos + n]
                && query_seq.at(q_pos + 1) == subject_seq[s_pos + n + 1]
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1025-1103
fn extend_left<Q: WalkLetters>(
    query_seq: Q,
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

        while q_pos > 0 && s_pos > 0 && query_seq.at(q_pos) == subject_seq[s_pos] {
            num_identical += 1;
            q_pos -= 1;
            s_pos -= 1;
        }

        for n in 1..max_shift {
            if q_pos <= n + 1 || s_pos <= n + 1 || matched {
                break;
            }

            if query_seq.at(q_pos - n) == subject_seq[s_pos - n]
                && query_seq.at(q_pos - n - 1) == subject_seq[s_pos - n - 1]
            {
                q_pos -= n + 2;
                s_pos -= n + 2;
                num_identical += 2;
                matched = true;
            }

            if !matched
                && query_seq.at(q_pos - n) == subject_seq[s_pos]
                && query_seq.at(q_pos - n - 1) == subject_seq[s_pos - 1]
            {
                q_pos -= n + 2;
                s_pos -= 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            if !matched
                && query_seq.at(q_pos) == subject_seq[s_pos - n]
                && query_seq.at(q_pos - 1) == subject_seq[s_pos - n - 1]
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1124-1204
fn find_num_identical<Q: WalkLetters, W: WalkWords>(
    query_seq: Q,
    query_hashes: W,
    subject_seq: &[u8],
    max_shift: usize,
) -> usize {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1160-1165
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1166-1188
    //
    // `s_FindNumIdentical` indexes `query_hashes` only up to
    // `query_len - word_size - 1`, so a query slice of length `L` needs
    // exactly `L - word_size` usable hashes here (`WalkWords::missing`).
    if query_seq.len() < NEAR_IDENTICAL_WORD_SIZE
        || query_hashes.missing(query_seq.len())
        || subject_seq.len() < NEAR_IDENTICAL_WORD_SIZE
    {
        return 0;
    }

    let mut hash = 0u64;
    num_identical_by(
        query_seq,
        subject_seq,
        max_shift,
        |s_pos, fresh_hash, query_positions| {
            if fresh_hash {
                hash = get_hash(&subject_seq[s_pos..], NEAR_IDENTICAL_WORD_SIZE);
            } else {
                hash <<= 5;
                hash &= NEAR_IDENTICAL_HASH_MASK;
                hash += u64::from(subject_seq[s_pos + NEAR_IDENTICAL_WORD_SIZE - 1]);
            }
            query_positions
                .into_iter()
                .find(|&q_pos| query_hashes.word(q_pos) == hash)
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1166-1188
/// The loop of `s_FindNumIdentical`, on `query_seq` and `subject_seq` of at least the word size.
/// `find` gets the subject position, whether the hash loop must hash the word afresh, and the query
/// positions to test in order; it returns the first of them whose word equals the subject's.
fn num_identical_by<Q: WalkLetters>(
    query_seq: Q,
    subject_seq: &[u8],
    max_shift: usize,
    mut find: impl FnMut(usize, bool, std::ops::Range<usize>) -> Option<usize>,
) -> usize {
    let mut query_from = 0usize;
    let mut subject_from = 0usize;
    let mut num_identical = 0usize;
    let mut matched_last = false;

    let last_subject_start = subject_seq.len() - NEAR_IDENTICAL_WORD_SIZE;
    let last_query_start = query_seq.len() - NEAR_IDENTICAL_WORD_SIZE;
    let mut s_pos = 0usize;
    while s_pos < last_subject_start {
        let matched_q = find(
            s_pos,
            s_pos == 0 || matched_last,
            query_from..last_query_start,
        );

        if let Some(query_start) = matched_q {
            let subject_start = s_pos;
            let (left_identical, _query_left_len, _subject_left_len, _align_len_left) = extend_left(
                query_seq.range(query_from, query_start),
                &subject_seq[subject_from..subject_start],
                max_shift,
            );
            let (right_identical, query_right_len, subject_right_len, _align_len_right) =
                extend_right(
                    query_seq.range(query_start + NEAR_IDENTICAL_WORD_SIZE, query_seq.len()),
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1153-1188 (see
// `find_num_identical`); c++/src/algo/blast/core/blast_kappa.c:2235-2272 (`build_query_word_hashes`)
/// EXPERIMENT (LOSAT_X_BXENDMEM): `find_num_identical` without the table of the query's word
/// hashes: `word_residues` are the residues of the whole query context from the first position the
/// table slice would start at. A hash is its word's letters in 5 bits each and residues are below
/// 32, so two hashes are equal exactly where the 8 letters are; the table's last entry is the one
/// its loop leaves 0 (`build_query_word_hashes`), equal only to a subject word of 8 letters 0.
fn x_find_num_identical_on_residues(
    query_seq: &[u8],
    word_residues: &[u8],
    table_entries: usize,
    last_is_zero: bool,
    subject_seq: &[u8],
    max_shift: usize,
) -> usize {
    let required_query_hashes = query_seq.len().saturating_sub(NEAR_IDENTICAL_WORD_SIZE);
    if query_seq.len() < NEAR_IDENTICAL_WORD_SIZE
        || table_entries < required_query_hashes
        || subject_seq.len() < NEAR_IDENTICAL_WORD_SIZE
    {
        return 0;
    }
    let word = |residues: &[u8], at: usize| {
        let bytes: [u8; NEAR_IDENTICAL_WORD_SIZE] = residues[at..at + NEAR_IDENTICAL_WORD_SIZE]
            .try_into()
            .expect("a slice of the word size");
        u64::from_ne_bytes(bytes)
    };
    let zero_entry = table_entries.wrapping_sub(1);
    num_identical_by(
        query_seq,
        subject_seq,
        max_shift,
        |s_pos, _, query_positions| {
            let subject_word = word(subject_seq, s_pos);
            query_positions.into_iter().find(|&q_pos| {
                if last_is_zero && q_pos == zero_entry {
                    subject_word == 0
                } else {
                    word(word_residues, q_pos) == subject_word
                }
            })
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1258-1312
pub fn test_near_identical(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    query_data: &BlastCompoSequenceData,
    query_offset: i32,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
) -> bool {
    near_identical_by(
        seq_data,
        seq_offset,
        query_data.data(),
        query_offset,
        align,
        |query_word_start, query_seq, subject_seq| {
            let query_words = query_words.unwrap_or(&[]);
            let query_word_slice = if query_word_start <= query_words.len() {
                &query_words[query_word_start..]
            } else {
                &[]
            };
            find_num_identical(
                query_seq,
                query_word_slice,
                subject_seq,
                NEAR_IDENTICAL_MAX_SHIFT,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1322-1333 (see
// `test_near_identical`)
/// EXPERIMENT (LOSAT_X_BXENDMEM): `test_near_identical` with `query_words` the table
/// `build_query_word_hashes(query_context)` would make, read from `query_context` instead
/// (`x_find_num_identical_on_residues`).
pub(crate) fn x_test_near_identical_on_residues(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    query_data: &BlastCompoSequenceData,
    query_offset: i32,
    query_context: &[u8],
    align: &BlastCompoAlignment,
) -> bool {
    near_identical_by(
        seq_data,
        seq_offset,
        query_data.data(),
        query_offset,
        align,
        |query_word_start, query_seq, subject_seq| {
            // The table has `length - word size + 1` entries; a start past them leaves no slice.
            let entries = (query_context.len() + 1).saturating_sub(NEAR_IDENTICAL_WORD_SIZE);
            let table_entries = entries.saturating_sub(query_word_start);
            x_find_num_identical_on_residues(
                query_seq,
                query_context.get(query_word_start..).unwrap_or(&[]),
                table_entries,
                // Entry 0 is hashed even in a table of one entry.
                entries > 1,
                subject_seq,
                NEAR_IDENTICAL_MAX_SHIFT,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1258-1312
/// `test_near_identical` with the query's copy and words of an RPS-BLAST redo (`NearIdenticalCopy`:
/// the profile's letters, which can be fewer than the columns of its preliminary alignments): the
/// same walk, its reads past the copy or the words reported by `copy` instead of panicking.
pub(crate) fn test_near_identical_on_copy(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    copy: NearIdenticalCopy<'_>,
    query_offset: i32,
    align: &BlastCompoAlignment,
) -> bool {
    near_identical_by(
        seq_data,
        seq_offset,
        copy.letters(),
        query_offset,
        align,
        |query_word_start, query_seq, subject_seq| {
            find_num_identical(
                query_seq,
                copy.words_from(query_word_start),
                subject_seq,
                NEAR_IDENTICAL_MAX_SHIFT,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1258-1312
/// The statements of `s_TestNearIdentical` on the query letters `query` (NCBI's `queryData->data`,
/// from its first letter); `find_identical` is its call of `s_FindNumIdentical` (the query words'
/// start, the query's middle part, the subject's middle slice).
fn near_identical_by<Q: WalkLetters>(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    query: Q,
    query_offset: i32,
    align: &BlastCompoAlignment,
    find_identical: impl FnOnce(usize, Q, &[u8]) -> usize,
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

    if s_end >= seq_data.data().len() || q_start > q_end || s_start > s_end {
        return false;
    }
    let Some(query_seq) = query.window(q_start, q_end) else {
        return false;
    };

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
        query_seq.range(query_right_len, query_seq.len()),
        &subject_seq[subject_right_len..],
        NEAR_IDENTICAL_MAX_SHIFT,
    );
    num_identical += left_identical;
    if query_left_len + query_right_len >= query_seq.len()
        || subject_left_len + subject_right_len >= subject_seq.len()
    {
        return (num_identical as f64) / (align_len as f64) > MIN_FRACTION_NEAR_IDENTICAL;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1322-1333
    //
    // If the remaining middle slice is shorter than the 8-mer word size,
    // `s_FindNumIdentical` returns 0 and NCBI still evaluates the final
    // identity fraction. Do not convert a missing or exhausted hash slice
    // into an immediate `false`.
    num_identical += find_identical(
        q_start + query_right_len,
        query_seq.range(query_right_len, query_seq.len() - query_left_len),
        &subject_seq[subject_right_len..subject_seq.len() - subject_left_len],
    );

    (num_identical as f64) / (align_len as f64) > MIN_FRACTION_NEAR_IDENTICAL
}

/// The query letters that the near-identical walk reads (`s_TestNearIdentical`,
/// `s_ExtendRight`, `s_ExtendLeft`, `s_FindNumIdentical`): a slice (BLASTP, PSIBLAST, TBLASTN,
/// BLASTX and RPS-BLAST's consistent profiles alike), or an RPS-BLAST profile's copy of fewer
/// letters than its alignments' columns (`NearIdenticalCopy`). Each function of the walk is one
/// generic function; the slice's instance is the walk on slices as before.
pub(crate) trait WalkLetters: Copy {
    /// The number of letters the walk may read (NCBI's `query_len` of the part).
    fn len(self) -> usize;
    /// Whether `len` is 0.
    #[inline(always)]
    fn is_empty(self) -> bool {
        self.len() == 0
    }
    /// The letter at `i` (below `len`).
    fn at(self, i: usize) -> u8;
    /// The letters `from..to` of these (NCBI's pointer plus `from`, with the length `to - from`).
    fn range(self, from: usize, to: usize) -> Self;
    /// The letters `start..=end` of the whole query that `s_TestNearIdentical` walks, `None` where
    /// LOSAT's slice ends before them (the walk then returns `false`).
    fn window(self, start: usize, end: usize) -> Option<Self>;
}

impl WalkLetters for &[u8] {
    #[inline(always)]
    fn len(self) -> usize {
        <[u8]>::len(self)
    }

    #[inline(always)]
    fn at(self, i: usize) -> u8 {
        self[i]
    }

    #[inline(always)]
    fn range(self, from: usize, to: usize) -> Self {
        &self[from..to]
    }

    #[inline(always)]
    fn window(self, start: usize, end: usize) -> Option<Self> {
        self.get(start..=end)
    }
}

/// The query words that `s_FindNumIdentical` reads (`s_CreateWordArray`'s table, from the walk's
/// start): a slice, or an RPS-BLAST profile's (`NearIdenticalCopy::words_from`).
pub(crate) trait WalkWords: Copy {
    /// Whether `s_FindNumIdentical` returns 0 before it reads a word of a part of `query_len`
    /// letters: NCBI's test of a null table; a slice of fewer words than its loop reads.
    fn missing(self, query_len: usize) -> bool;
    /// The word at `i`.
    fn word(self, i: usize) -> u64;
}

impl WalkWords for &[u64] {
    #[inline(always)]
    fn missing(self, query_len: usize) -> bool {
        self.len() < query_len.saturating_sub(NEAR_IDENTICAL_WORD_SIZE)
    }

    #[inline(always)]
    fn word(self, i: usize) -> u64 {
        self[i]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1671-1708
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2242-2262
/// The copy of the query that NCBI's redo hands the near-identical walk (`s_SequenceGetRange`:
/// the letters, U as C, between two zero bytes in a block of the length plus 2) and its word table
/// (`s_CreateWordArray`: length minus 7 words, none below 8 letters, where NCBI's table is a null
/// pointer), read where the walk reads them: RPS-BLAST's redo of a profile whose `.psq` letters
/// are fewer than the columns of its preliminary alignments, where NCBI's walk can read past the
/// block or the table (`rps_search::profile`, read C). A letter at or past the block's last byte
/// and a word past the table, or of the null table at a start other than 0, read as 0 and are
/// reported in `read_past` (NCBI's read is out of bounds; LOSAT then refuses the query batch).
#[derive(Clone, Copy)]
pub(crate) struct NearIdenticalCopy<'a> {
    /// The block from its second byte: the letters and the zero byte after them.
    bytes: &'a [u8],
    /// The length of the query (the copy's letters).
    length: usize,
    words: Option<&'a [u64]>,
    read_past: &'a Cell<bool>,
}

impl<'a> NearIdenticalCopy<'a> {
    /// The copy `query` (`BlastCompoSequenceData::from_ncbistdaa` of the letters, U as C) and the
    /// words of its letters (`build_query_word_hashes`, empty below 8 letters).
    pub(crate) fn new(
        query: &'a BlastCompoSequenceData,
        words: &'a [u64],
        read_past: &'a Cell<bool>,
    ) -> Self {
        Self {
            bytes: &query.buffer[query.data_offset..],
            length: query.length as usize,
            words: (!words.is_empty()).then_some(words),
            read_past,
        }
    }

    fn letters(self) -> CopyLetters<'a> {
        CopyLetters {
            bytes: self.bytes,
            len: self.length,
            read_past: self.read_past,
        }
    }

    /// The table from word `start` (NCBI's `query_words + start`).
    fn words_from(self, start: usize) -> CopyWords<'a> {
        CopyWords {
            words: self
                .words
                .and_then(|words| words.get(start..))
                .unwrap_or(&[]),
            null: self.words.is_none() && start == 0,
            read_past: self.read_past,
        }
    }
}

/// A read of the walk past a `NearIdenticalCopy` (its letters or its words): recorded in
/// `read_past`, out of the walk's loops.
#[cold]
#[inline(never)]
fn copy_read_past(read_past: &Cell<bool>) {
    read_past.set(true);
}

/// Letters of a `NearIdenticalCopy` from a start: `bytes` the block's bytes from there (none when
/// the start lies past the block), `len` the letters the walk takes.
#[derive(Clone, Copy)]
struct CopyLetters<'a> {
    bytes: &'a [u8],
    len: usize,
    read_past: &'a Cell<bool>,
}

impl WalkLetters for CopyLetters<'_> {
    #[inline]
    fn len(self) -> usize {
        self.len
    }

    #[inline]
    fn at(self, i: usize) -> u8 {
        match self.bytes.get(i) {
            Some(&letter) => letter,
            None => {
                copy_read_past(self.read_past);
                0
            }
        }
    }

    #[inline]
    fn range(self, from: usize, to: usize) -> Self {
        Self {
            bytes: self.bytes.get(from..).unwrap_or(&[]),
            len: to - from,
            ..self
        }
    }

    // NCBI's walk takes the alignment's range whatever the copy's length (`s_TestNearIdentical`).
    #[inline]
    fn window(self, start: usize, end: usize) -> Option<Self> {
        Some(self.range(start, end + 1))
    }
}

/// Words of a `NearIdenticalCopy` from a word start: `words` the table from there (none when the
/// table is null or ends before it), `null` whether NCBI's pointer is the null table itself.
#[derive(Clone, Copy)]
struct CopyWords<'a> {
    words: &'a [u64],
    null: bool,
    read_past: &'a Cell<bool>,
}

impl WalkWords for CopyWords<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1161-1166
    // Only a null table at offset 0 is a null pointer; a null table plus a start passes the test.
    #[inline]
    fn missing(self, _query_len: usize) -> bool {
        self.null
    }

    #[inline]
    fn word(self, i: usize) -> u64 {
        match self.words.get(i) {
            Some(&word) => word,
            None => {
                copy_read_past(self.read_past);
                0
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1053-1097
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1292
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1146
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3631-3641
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234,1240-1241,1256-1256,1259-1260
// No NCBI counterpart: probe of the BLASTX parallel redo. It records, per match, whether
// `Blast_AdjustScores` ran (bit 0) and whether an alignment was redone first (bit 1), so a caller
// that redoes matches out of order can tell which matrix a redo used. It is read by the BLASTX
// parallel redo and does not change any value NCBI computes.
thread_local! {
    /// What one match did with the matrix it was handed.  Bit 0: an adjustment succeeded (the
    /// matrix was replaced); bit 1: an alignment was redone before that, i.e. with the matrix left
    /// behind by an earlier match (redo_alignment.c:1232-1259 skips `Blast_AdjustScores` for a
    /// protein subject unless `hsp_index == 0`).  A caller that redoes matches out of order resets
    /// this before a match and reads it afterwards.
    pub(crate) static X_REDO_PROBE: Cell<u8> = const { Cell::new(0) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1191-1211
// No NCBI counterpart: identifies one `get_range` call of the window loop (window, HSP index, and the
// `nearIdenticalStatus` and `subject_maybe_biased` values it is called with); the arguments of the
// call are fixed by these. It does not change any value NCBI computes.
/// Key of a `get_range` result made ahead.
#[derive(Clone, Copy, PartialEq, Eq)]
struct XLaneRangeKey {
    window: usize,
    hsp_index: usize,
    should_test_identical: bool,
    subject_maybe_biased: bool,
}

// No NCBI counterpart: removes and returns the `get_range` result made ahead for this call, if any;
// it does not change any value NCBI computes.
fn x_take_lane_range(
    ranges: &mut Vec<(XLaneRangeKey, BlastRedoRangeResult)>,
    key: XLaneRangeKey,
) -> Option<BlastRedoRangeResult> {
    let position = ranges.iter().position(|(stored, _)| *stored == key)?;
    Some(ranges.swap_remove(position).1)
}

// No NCBI counterpart: number of Newton inputs gathered per TBLASTN prefetch (at least); it does not
// change any value NCBI computes. Each prefetch adds its results to the thread's glibc arena at
// once, and the space they leave seldom fits the next prefetch's, so fewer prefetches push the
// arena's top less often (line EB s26: with 16 / 32 / 64 / 128 inputs, the peak of TBLASTN
// `-subject` at 64 threads was +17 % / base / -3 % / +18 %; 128 holds too many results at once).
const X_LANES_PREFETCH_INPUTS: usize = 64;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1176-1198
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1239-1244
// Newton lanes, TBLASTN (translated subject, protein query): a replica of the
// window loop of `blast_redo_one_match_with_workspace_queries_and_matrix` from window
// `first_window` on: the same `query_info` choice, the same near-identical test and `get_range` calls
// with the same arguments (results kept in `ranges` for the loop), and the subject composition of
// every HSP. Containment is not tested (every HSP is assumed not contained), so the inputs are a
// superset of those the loop adjusts. Stops after the window in which `X_LANES_PREFETCH_INPUTS`
// inputs are reached and solves them (`x_newton_lanes::prefetch`). Errors end the replica; the loop
// meets them itself.
// No NCBI counterpart: computes ahead what the loop computes; it does not change any value NCBI
// computes.
/// Returns the first window not covered.
#[allow(clippy::too_many_arguments)]
fn x_newton_lanes_prefetch_windows<'seq>(
    windows: &[WindowInfo],
    first_window: usize,
    params: &BlastRedoAlignParams,
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    query_infos: &[BlastCompoQueryInfo],
    callbacks: &BlastRedoAlignCallbacks,
    composition_workspace: &BlastCompositionWorkspace,
    ranges: &mut Vec<(XLaneRangeKey, BlastRedoRangeResult)>,
) -> usize {
    let mut inputs = Vec::new();
    let mut shared_q = None;
    let mut window_index = first_window;
    while window_index < windows.len() && inputs.len() < X_LANES_PREFETCH_INPUTS {
        let window = &windows[window_index];
        window_index += 1;
        let Some(query_info) = usize::try_from(window.query_range.context)
            .ok()
            .and_then(|query_index| query_infos.get(query_index))
        else {
            window_index = windows.len();
            break;
        };
        let Some(window_align_head) = window.align.as_deref() else {
            continue;
        };
        let mut range_position = 0usize;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut hsp_index = 0usize;
        let mut in_align = window.align.as_deref();
        while let Some(current) = in_align {
            let near_identical = if hsp_index == 0 || subject_maybe_biased {
                preliminary_test_near_identical(
                    query_info,
                    window,
                    current,
                    params.near_identical_cutoff,
                )
            } else {
                old_near_identical
            };
            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                let key = XLaneRangeKey {
                    window: window_index - 1,
                    hsp_index,
                    should_test_identical: near_identical,
                    subject_maybe_biased,
                };
                let Ok(next_range) = (callbacks.get_range)(
                    matching_seq,
                    &window.subject_range,
                    &query_info.seq,
                    &window.query_range,
                    query_info.words.as_deref(),
                    window_align_head,
                    near_identical,
                    subject_maybe_biased,
                    params,
                ) else {
                    super::x_newton_lanes::prefetch(inputs);
                    return windows.len();
                };
                subject_maybe_biased = next_range.subject_maybe_biased;
                ranges.push((key, next_range));
                range_position = ranges.len() - 1;
            }
            let range = &ranges[range_position].1;
            if let Ok(subject_composition) = translated_subject_composition(
                &range.subject.buffer,
                &window.subject_range,
                current,
            ) {
                if query_info.composition.num_true_amino_acids != 0
                    && subject_composition.num_true_amino_acids != 0
                {
                    if let Some(input) = super::adjust_scores::x_newton_input(
                        &params.matrix_info,
                        &query_info.composition,
                        range.query.length,
                        &subject_composition,
                        range.subject.length,
                        params.compo_adjust_mode,
                        composition_workspace,
                        &mut shared_q,
                    ) {
                        inputs.push(input);
                    }
                }
            }
            old_near_identical = near_identical;
            in_align = current.next.as_deref();
            hsp_index += 1;
        }
    }
    super::x_newton_lanes::prefetch(inputs);
    window_index
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1274
pub(crate) fn blast_redo_one_match_with_workspace_queries_and_matrix<'seq>(
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    params: &BlastRedoAlignParams,
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    query_infos: &[BlastCompoQueryInfo],
    lambda: f64,
    callbacks: &BlastRedoAlignCallbacks,
    composition_workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
) -> Result<BlastRedoOneMatchResult> {
    if params.smith_waterman {
        bail!("blastp redo_alignment Smith-Waterman path is not yet ported");
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:865-889
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:879-888
    let windows = if params.query_is_translated {
        windows_from_translated_query_aligns(
            incoming_aligns,
            query_infos,
            matching_seq.length,
            WINDOW_BORDER,
        )?
    } else if params.subject_is_translated {
        if matching_seq.genetic_code_id.is_none() {
            bail!("TBLASTN redo requires a translated subject genetic code");
        }
        let query_lengths: Vec<_> = query_infos.iter().map(|info| info.seq.length).collect();
        windows_from_translated_aligns(
            incoming_aligns,
            &query_lengths,
            matching_seq.length,
            WINDOW_BORDER,
            params.position_based,
        )?
    } else {
        let query_info = query_infos
            .first()
            .ok_or_else(|| anyhow::anyhow!("BLASTP redo requires query information"))?;
        windows_from_protein_aligns(incoming_aligns, query_info, matching_seq.length)?
    };
    let mut alignments_by_query: Vec<Option<Box<BlastCompoAlignment>>> = Vec::new();
    // No NCBI counterpart: LOSAT_X_REDOCONTAINIDX moves a query's list into an index by subject
    // start once it is long (`x_align_index.rs`); the answers and the final lists are NCBI's.
    let x_index_from =
        super::x_align_index::x_redo_contain_idx().then_some(super::x_align_index::X_INDEX_FROM);
    let mut x_indexes: XAlignIndexes = Vec::new();
    let mut pvalue_for_this_pair = None;
    let mut lambda_ratio = None;
    let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1439-1444
    let mut adjusted_matrix = matrix_state.take();
    // TBLASTN Newton inputs of the coming windows solved ahead, four per vector (eight on aarch64,
    // x_newton_lanes.rs); their get_range results are kept for this loop.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1143-1143
    // No NCBI counterpart: the windows of this match are walked ahead to collect the inputs of their
    // `Blast_AdjustScores` calls; a result is used only by a call with the same input bits, and a kept
    // range only by the `get_range` call with the same window, HSP index and flags. It does not change
    // any value NCBI computes.
    let x_lanes_tblastn = params.subject_is_translated && !params.query_is_translated;
    let x_lanes_on = x_lanes_tblastn
        && params.compo_adjust_mode != BlastCompoAdjustMode::NoCompositionBasedStats;
    if x_lanes_tblastn {
        crate::utils::xreach::x_reach!("LOSAT_X_NEWTONLANES", x_lanes_on);
    }
    let mut x_next_prefetch = 0usize;
    let mut x_lane_ranges: Vec<(XLaneRangeKey, BlastRedoRangeResult)> = Vec::new();

    for (window_index, window) in windows.iter().enumerate() {
        // No NCBI counterpart: Newton lanes prefetch (see above); it does not change any value NCBI
        // computes.
        if x_lanes_on && window_index >= x_next_prefetch {
            x_next_prefetch = x_newton_lanes_prefetch_windows(
                &windows,
                window_index,
                params,
                matching_seq,
                query_infos,
                callbacks,
                composition_workspace,
                &mut x_lane_ranges,
            );
        }
        let query_index = usize::try_from(window.query_range.context)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1162-1188
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1162-1168
        let query_info = if params.subject_is_translated || params.query_is_translated {
            query_infos.get(query_index).ok_or_else(|| {
                anyhow::anyhow!("TBLASTN redo query index {query_index} is missing")
            })?
        } else {
            &query_infos[0]
        };
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
                    window,
                    current,
                    params.near_identical_cutoff,
                )
            } else {
                old_near_identical
            };

            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1188-1210
                // No NCBI counterpart: the result of this call made ahead (same window, HSP index,
                // near-identical status and bias flag) is used. It does not change any value NCBI
                // computes.
                let x_stored = if x_lanes_on {
                    x_take_lane_range(
                        &mut x_lane_ranges,
                        XLaneRangeKey {
                            window: window_index,
                            hsp_index,
                            should_test_identical: near_identical,
                            subject_maybe_biased,
                        },
                    )
                } else {
                    None
                };
                let next_range = if let Some(stored) = x_stored {
                    stored
                } else {
                    (callbacks.get_range)(
                        matching_seq,
                        &window.subject_range,
                        &query_info.seq,
                        &window.query_range,
                        query_info.words.as_deref(),
                        window_align_head,
                        near_identical,
                        subject_maybe_biased,
                        params,
                    )?
                };
                subject_maybe_biased = next_range.subject_maybe_biased;
                range_result = Some(next_range);
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1218-1224
            let query_composition = if params.query_is_translated {
                translated_query_composition(
                    &range_result.as_ref().expect("range initialized").query,
                    &window.query_range,
                    current,
                )?
            } else {
                query_info.composition.clone()
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1225-1227
            let contained = x_contains(
                &x_indexes,
                query_index,
                &alignments_by_query[query_index],
                current,
                lambda,
            );
            if !contained {
                let mut adjust_search_failed = false;
                let range = range_result
                    .as_ref()
                    .expect("range result initialized before redo_one_alignment");

                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1219-1256
                if params.compo_adjust_mode != BlastCompoAdjustMode::NoCompositionBasedStats
                    && (params.subject_is_translated
                        || hsp_index == 0
                        || (near_identical != old_near_identical))
                {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:906-938,1234-1238
                    let subject_composition = if params.subject_is_translated {
                        translated_subject_composition(
                            &range.subject.buffer,
                            &window.subject_range,
                            current,
                        )?
                    } else {
                        x_whole_subject_composition(range)
                    };
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1239-1256
                    if query_composition.num_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1228-1249
                        let adjusted = blast_adjust_scores(
                            &params.matrix_info,
                            &query_composition,
                            range.query.length,
                            &subject_composition,
                            range.subject.length,
                            params.compo_adjust_mode,
                            params.composition_test_index,
                            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1240-1247
                            composition_workspace,
                            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1239-1248
                            callbacks.calc_lambda.unwrap_or(redo_calc_lambda),
                        )?;
                        if let Some(adjusted) = adjusted {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1240-1241,1256-1256
                            // No NCBI counterpart: sets bit 0 of the probe when `Blast_AdjustScores` returned a matrix; the
                            // probe is not read by this function and does not change any value NCBI computes.
                            X_REDO_PROBE.with(|probe| probe.set(probe.get() | 1));
                            matrix_adjust_rule = adjusted.matrix_adjust_rule;
                            adjusted_matrix = Some(adjusted.adjusted_matrix);
                            pvalue_for_this_pair = adjusted.pvalue_for_this_pair;
                            lambda_ratio = Some(adjusted.lambda_ratio);
                        } else {
                            adjust_search_failed = true;
                        }
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1251-1256
                    num_adjustments += 1;
                }

                if !adjust_search_failed {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1259-1274
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1234,1259-1260
                    // No NCBI counterpart: sets bit 1 of the probe when `redo_one_alignment` runs without a
                    // `Blast_AdjustScores` call in this match (the protein-subject case `hsp_index != 0` skips it),
                    // so the redo uses the matrix left by an earlier match; it does not change any value NCBI computes.
                    X_REDO_PROBE.with(|probe| {
                        if probe.get() & 1 == 0 {
                            probe.set(probe.get() | 2);
                        }
                    });
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1268-1269,1272
                        x_add_distinct(
                            &mut x_indexes,
                            query_index,
                            &mut alignments_by_query[query_index],
                            &mut new_align,
                            num_adjustments == 1,
                            x_index_from,
                        );
                    }
                }
            }

            old_near_identical = near_identical;
            in_align = current.next.as_deref();
            hsp_index += 1;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1272-1274
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1277-1304
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1239-1259
    if x_index_from.is_some() {
        let x_indexed = x_finish(x_indexes, &mut alignments_by_query);
        crate::utils::xreach::x_reach!("LOSAT_X_REDOCONTAINIDX", x_indexed);
    }
    *matrix_state = adjusted_matrix;
    Ok(BlastRedoOneMatchResult {
        alignments_by_query,
        pvalue_for_this_pair,
        lambda_ratio,
    })
}

#[cfg(test)]
#[path = "redo_alignment_tests.rs"]
mod tests;
