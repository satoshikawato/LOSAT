//! Internal TBLASTN protein lookup and translated-subject seed stage: the scan plan and the
//! per-unit scan of the preliminary search.

use crate::algorithm::blastp::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_query_frame_with_seg_timed,
};
use crate::algorithm::tblastx::blast_aascan::{s_blast_aa_scan_subject_one_range, BlastOffsetPair};
#[cfg(test)]
use crate::algorithm::tblastx::lookup::{build_ncbi_lookup, build_ncbi_lookup_for_profile};
#[cfg(test)]
use crate::algorithm::tblastx::translation::generate_frames;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::config::ScoringMatrix;
#[cfg(test)]
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::ncbistdaa;
use crate::utils::seg::SegParams;
use anyhow::{bail, Result};

// NCBI c++/src/algo/blast/blastinput/blast_fasta_input.cpp:486-502;
// c++/src/algo/blast/api/blast_setup_cxx.cpp:655-657;
// c++/src/algo/blast/core/blast_filter.c:1241-1255;
// c++/src/algo/blast/core/blast_setup.c:614-638:
// BlastSeqLocAppend(filter_out, lcase_mask_slp);
// BlastSeqLocCombine(filter_out, 0);
// if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
// BLAST_ComplementMaskLocations(..., filter_maskloc, lookup_segments);
// query lowercase spans join the filtering locations; their complement
// defines lookup segments, and hard masking replaces those residues with X.
// Keep the unmasked copy for identity calculation after traceback.
pub(super) fn encode_tblastn_lookup_query(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> QueryFrame {
    encode_tblastn_lookup_query_with(query, seg, mask_lowercase, false)
}

/// `encode_tblastn_lookup_query`; with LOSAT_TIMING the encoding and the SEG go to the common
/// line (utils/xtiming.rs). For the driver's query preparation on the calling thread only.
pub(super) fn encode_tblastn_lookup_query_timed(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> QueryFrame {
    encode_tblastn_lookup_query_with(query, seg, mask_lowercase, true)
}

fn encode_tblastn_lookup_query_with(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    timed: bool,
) -> QueryFrame {
    let mut frame = if timed {
        encode_protein_query_frame_with_seg_timed(query, seg)
    } else {
        encode_protein_query_frame_with_seg(query, seg)
    };
    if !mask_lowercase {
        return frame;
    }
    let mut start = None;
    for (index, &residue) in query.iter().enumerate() {
        if residue.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(left) = start.take() {
            frame.seg_masks.push((left, index));
        }
    }
    if let Some(left) = start {
        frame.seg_masks.push((left, query.len()));
    }
    if frame.seg_masks.is_empty() {
        return frame;
    }
    frame.seg_masks.sort_unstable_by_key(|&(left, _)| left);
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (left, right) in frame.seg_masks.drain(..) {
        if let Some(last) = merged.last_mut() {
            if left <= last.1 {
                last.1 = last.1.max(right);
                continue;
            }
        }
        merged.push((left, right));
    }
    frame.seg_masks = merged;
    if frame.aa_seq_nomask.is_none() {
        frame.aa_seq_nomask = Some(frame.aa_seq.clone());
    }
    for &(left, right) in &frame.seg_masks {
        frame.aa_seq[1 + left..1 + right].fill(ncbistdaa::X);
    }
    frame
}

/// The runs of lower-case letters of `query`, as half-open intervals, appended to `out`.
fn push_lowercase_runs(query: &[u8], out: &mut Vec<(usize, usize)>) {
    let mut start = None;
    for (index, &residue) in query.iter().enumerate() {
        if residue.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(left) = start.take() {
            out.push((left, index));
        }
    }
    if let Some(left) = start {
        out.push((left, query.len()));
    }
}

/// The SEG intervals of `query` and, with `mask_lowercase`, its lower-case runs (half-open,
/// unsorted, not joined), for the masks of the BLAST XML report, which joins them as NCBI's
/// `BlastSeqLocCombine(filter_out, 0)` does: overlapping intervals only. The search's list
/// (`encode_tblastn_lookup_query`) also joins adjacent ones; that changes no masked letter,
/// but NCBI's XML mask loop (`s_MaskQuerySeq`) leaves a letter unmasked where one interval
/// ends on the second-last column of an alignment segment and the next begins on its last.
pub(super) fn tblastn_report_mask_intervals(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> Vec<(usize, usize)> {
    let mut intervals = encode_protein_query_frame_with_seg(query, seg).seg_masks;
    if mask_lowercase {
        push_lowercase_runs(query, &mut intervals);
    }
    intervals
}

// NCBI reference: c++/include/algo/blast/core/blast_def.h:141-150
// typedef union BlastOffsetPair {
//     struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
// } BlastOffsetPair;
// NCBI c++/src/algo/blast/core/blast_engine.c:808-812:
// subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Seed {
    pub frame: i8,
    // NCBI c++/src/algo/blast/core/blast_engine.c:241-250:
    // backup.offset identifies this chunk before WordFinder returns local offsets.
    pub chunk_offset: u32,
    pub query_offset: u32,
    pub subject_offset: u32,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
// if (backup->offset + MAX_DBSEQ_LEN < hard_ranges[hm_index].right) {
//     subject->length = MAX_DBSEQ_LEN;
//     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
// } else { subject->length = hard_ranges[hm_index].right - backup->offset; }
// NCBI c++/include/algo/blast/core/blast_gapalign.h:54:
// #define MAX_DBSEQ_LEN 5000000
// NCBI c++/include/algo/blast/core/blast_hits.h:192:
// #define DBSEQ_CHUNK_OVERLAP 100
// This schedule covers the unmasked translated frame's single hard range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct TranslatedChunk {
    pub offset: usize,
    pub length: usize,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
// backup->offset = backup->next;
// while (backup->next < backup->full_range.right) { s_GetNextSubjectChunk(...); }
pub(super) fn unmasked_translated_chunks(frame_length: usize) -> Vec<TranslatedChunk> {
    const MAX_DBSEQ_LEN: usize = 5_000_000;
    const DBSEQ_CHUNK_OVERLAP: usize = 100;
    let mut chunks = Vec::new();
    let mut offset = 0;
    while offset < frame_length {
        if offset + MAX_DBSEQ_LEN < frame_length {
            chunks.push(TranslatedChunk {
                offset,
                length: MAX_DBSEQ_LEN,
            });
            offset += MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP;
        } else {
            chunks.push(TranslatedChunk {
                offset,
                length: frame_length - offset,
            });
            break;
        }
    }
    chunks
}

// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:813-826
// BlastSeqBlkSetSeqRanges(subj, (SSeqRange*) masked_ranges.get_data(),
//                          masked_ranges.size() + 1, true, eSoftSubjMasking);
// NCBI c++/src/algo/blast/blastinput/blast_fasta_input.cpp:489-502:
// apply_mask_to_both_strands = true; PackedSeqLocToMaskedQueryRegions(...);
// NCBI c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
// if (assume_both_strands) do_pos = do_neg = true;
// if (do_pos) mqr.push_back(...); if (do_neg) mqr.push_back(...);
// NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:689-707,825-826:
// both mask copies enter masked_ranges before BlastSeqBlkSetSeqRanges.
// NCBI c++/src/algo/blast/core/blast_util.c:211-215:
// tmp[0].left = 0;
// tmp[num_seq_ranges - 1].right = seq_blk->length;
// NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
// if (context == 0) { seq_ranges[i].left = backup.seq_ranges[i].left / 3;
//                     seq_ranges[i].right = backup.seq_ranges[i].right / 3; }
// else if (context == 3) { /* reverse ranges using subject->length */ }
// NCBI reference: c++/src/algo/blast/core/masksubj.inl:51-71
// range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
// range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
fn translated_scan_ranges(
    subject: &[u8],
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    mask_lowercase: bool,
) -> Vec<(i32, i32)> {
    if !mask_lowercase || !subject.iter().any(u8::is_ascii_lowercase) {
        return vec![(0, frame_length as i32)];
    }
    let mut nucleotide_ranges = Vec::new();
    let mut left = 0;
    let mut index = 0;
    while index < subject.len() {
        if subject[index].is_ascii_lowercase() {
            let start = index;
            while index < subject.len() && subject[index].is_ascii_lowercase() {
                index += 1;
            }
            // NCBI c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
            // if (assume_both_strands) do_pos = do_neg = true;
            // if (do_pos) mqr.push_back(...); if (do_neg) mqr.push_back(...);
            for _ in 0..2 {
                nucleotide_ranges.push((left, start));
                left = index - 1;
            }
        } else {
            index += 1;
        }
    }
    nucleotide_ranges.push((left, subject.len()));
    if frame > 0 {
        nucleotide_ranges
            .into_iter()
            .map(|(left, right)| ((left / 3) as i32, (right / 3) as i32))
            .collect()
    } else {
        let length = negative_first_length as i32;
        nucleotide_ranges
            .into_iter()
            .rev()
            .map(|(left, right)| (length - (right / 3) as i32, length - (left / 3) as i32))
            .collect()
    }
}

// NCBI c++/src/algo/blast/core/blast_engine.c:259-310:
// unsplit subjects keep backup.soft_ranges; split soft-masked subjects
// advance sm_index while right < chunk start, include ranges with left < end,
// clamp only the first left and last right, and skip SUBJECT_SPLIT_NO_RANGE.
// NCBI c++/src/algo/blast/core/blast_engine.c:261-274:
// without soft masking, each chunk has one range [0, subject->length].
pub(super) fn translated_chunk_scan_ranges(
    subject: &[u8],
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    mask_lowercase: bool,
) -> Result<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    let frame_ranges = translated_scan_ranges(
        subject,
        frame,
        frame_length,
        negative_first_length,
        mask_lowercase,
    );
    let soft_masked = mask_lowercase && subject.iter().any(u8::is_ascii_lowercase);
    let mut sm_index = 0usize;
    let mut calls = Vec::new();
    for chunk in unmasked_translated_chunks(frame_length) {
        if !soft_masked {
            calls.push((chunk, vec![(0, i32::try_from(chunk.length)?)]));
            continue;
        }
        if chunk.offset == 0 && chunk.length == frame_length {
            calls.push((chunk, frame_ranges.clone()));
            continue;
        }
        let start = i32::try_from(chunk.offset)?;
        let end = i32::try_from(chunk.offset + chunk.length)?;
        let mut i = sm_index;
        while i < frame_ranges.len() && frame_ranges[i].1 < start {
            i += 1;
        }
        if i == frame_ranges.len() {
            bail!("NCBI translated soft-range sentinel is missing");
        }
        let first = i;
        while i < frame_ranges.len() && frame_ranges[i].0 < end {
            i += 1;
        }
        if i == 0 {
            bail!("NCBI translated soft-range index underflow");
        }
        sm_index = i - 1;
        if i == first {
            continue;
        }
        let mut ranges: Vec<_> = frame_ranges[first..i]
            .iter()
            .map(|&(left, right)| (left - start, right - start))
            .collect();
        ranges[0].0 = ranges[0].0.max(0);
        let last = ranges.len() - 1;
        ranges[last].1 = ranges[last].1.min(i32::try_from(chunk.length)?);
        calls.push((chunk, ranges));
    }
    Ok(calls)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ...
//     for (context=first_context; context<=last_context; context++) {
//         BlastHSPList* hsp_list_for_chunks = NULL;
//
//         if (kTranslatedSubject) {
//             Uint4 i;
//             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//             subject->sequence = translation_buffer + frame_offsets[context] + 1;
//             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
// ```
// NCBI translates the subject into its six frames once, then calls WordFinder for each frame and
// chunk with a read-only lookup table. XScanPlan holds what those calls need (the lookup table and
// the translated frames with their chunks and scan ranges) so that every unit can be scanned on its
// own.
// The scan state of the preliminary search, built once.  NCBI c++/src/algo/blast/core/
// blast_engine.c:478-500, 804-841 calls WordFinder once per unit; the lookup table is read-only
// there.
pub(super) struct XScanPlan {
    pub lookup: std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>,
    pub subject: std::sync::Arc<XSubjectPlan>,
    pub pair_capacity: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ...
//     for (context=first_context; context<=last_context; context++) {
//         BlastHSPList* hsp_list_for_chunks = NULL;
//
//         if (kTranslatedSubject) {
//             Uint4 i;
//             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//             subject->sequence = translation_buffer + frame_offsets[context] + 1;
//             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
// ```
// The six translated frames of one subject and the (frame, chunk, scan range) units in NCBI call
// order.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNSSIDE`): the part of a scan plan that depends on the
/// subject only.
pub(super) struct XSubjectPlan {
    pub frames: Vec<QueryFrame>,
    /// (index into `frames`, chunk, scan ranges) in NCBI call order.
    pub units: Vec<(usize, TranslatedChunk, Vec<(i32, i32)>)>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
// ```c
//     if (program_number != eBlastTypeRpsTblastn) {
//         if (translation_buffer) {
//             sfree(translation_buffer);
//         }
//     }
// ```
// NCBI frees a subject's translation after the subject; a dropped plan gives its frames to the
// next plan of the thread instead (LOSAT_X_TBNFRAMEREUSE).
impl Drop for XSubjectPlan {
    fn drop(&mut self) {
        for frame in self.frames.drain(..) {
            super::x_subject_translate::x_reuse_frame(frame.aa_seq);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ...
//     for (context=first_context; context<=last_context; context++) {
//         BlastHSPList* hsp_list_for_chunks = NULL;
//
//         if (kTranslatedSubject) {
//             Uint4 i;
//             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//             subject->sequence = translation_buffer + frame_offsets[context] + 1;
//             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
// ```
// NCBI translates a subject once per query batch (BLAST_GetAllTranslations in
// s_BlastSearchEngineCore). The translation depends only on the subject, the genetic code and the
// lowercase-masking option, so XSubjectSide keeps it for later batches. A reuse of the same value;
// no value is computed differently.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNSSIDE`; its rule stays: plans are kept from
/// `KEEP_FROM_BATCHES` batches on, up to `KEEP_LIMIT`): the subject plans of one search, by
/// subject index.
///
/// A plan (the six translated frames and their chunks) is a function of the
/// subject, the genetic code and the lowercase-masking option. The search
/// consists of one preliminary stage per query batch, and each of them made
/// the plan of every subject again. The owner of the search creates one value
/// for its fixed list of subjects and options; the first batch that reaches a
/// subject computes its plan, as before, and the later ones read it. Plans are
/// kept only up to a total size, beyond which they are made per use as before.
/// A search of fewer than `KEEP_FROM_BATCHES` batches gets no value
/// (`for_batches`): one batch has no later reader, and two or three give the
/// time their later readers would save for the memory of the plans. With
/// `LOSAT_X_TBNKEEPFEW` (`keep_few`) a search of more batches that searches them
/// one after another gets none either unless its plans hold no more than its
/// threads could have in flight (`kept_within_threads`): many subjects are
/// translated again for each batch, as NCBI translates them, and few long ones
/// are kept. A search that searches batches side by side keeps them
/// (`stage_d_pipeline::x_searched_side_by_side`).
///
/// Batches can reach a subject at the same time. The cell of a subject is
/// initialised by one of them while the others wait, so a plan is made once
/// and its size is added to `kept_bytes` once: `kept_bytes` is always the
/// size of the plans the cells hold.
pub(super) struct XSubjectSide {
    gencode: u8,
    mask_lowercase: bool,
    /// Per subject: unset until a batch reaches it; then the kept plan, or
    /// `None` when its plan is not kept (over the limit, or the subject could
    /// not be prepared) and every use makes its own.
    plans: Vec<std::sync::OnceLock<Option<std::sync::Arc<XSubjectPlan>>>>,
    kept_bytes: std::sync::atomic::AtomicUsize,
}

// No NCBI counterpart: the methods below keep, count and hand out the cached subject plans (a memo
// of x_build_subject); it does not change any value NCBI computes.
impl XSubjectSide {
    /// Translated residues kept at most (two bytes per subject base).
    pub(super) const KEEP_LIMIT: usize = 256 << 20;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
    // ```c
    //     if (program_number != eBlastTypeRpsTblastn) {
    //         if (translation_buffer) {
    //             sfree(translation_buffer);
    //         }
    //     }
    // ```
    // NCBI frees a subject's translation when the subject has been searched, in every query
    // batch. Below this number of batches the search does the same.
    /// Query batches (20,000 query residues each) from which a search keeps its
    /// subject plans. A kept plan saves the later batches the translation of
    /// its subject and holds the translated residues until the search ends.
    /// Timed on 28,036 subjects (40.6 Mb) at 4 and 8 threads: the plans hold
    /// 115 MiB, and a search of two, three or four batches without them took a
    /// median of 1.03 to 1.11 of the time with them. The saving does not start
    /// at this number; below it the search gives that time for the memory.
    pub(super) const KEEP_FROM_BATCHES: usize = 4;

    /// EXPERIMENT (LOSAT_X_TBNKEEPFEW): a switch of the strict set (`utils/xswitch.rs`), on
    /// unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it a search of
    /// `KEEP_FROM_BATCHES` or more batches keeps its plans only when `kept_within_threads`
    /// holds or it searches batches side by side (`stage_d_pipeline::x_searched_side_by_side`).
    /// `LOSAT_X_TBNKEEPFEWSHADOW` (`shadow`) keeps them as with `=0`.
    pub(super) fn keep_few() -> bool {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| !Self::shadow() && crate::utils::xswitch::on("LOSAT_X_TBNKEEPFEW"))
    }

    /// EXPERIMENT (LOSAT_X_TBNKEEPFEWSHADOW): set (opt-in). The plans are kept as with
    /// `LOSAT_X_TBNKEEPFEW=0`, and every read of a kept plan also makes the plan of that use, as
    /// a search that keeps none makes it, and compares the two (`plan`); a difference stops
    /// the search. `x_keep_few_print_summary` prints the count at exit.
    pub(super) fn shadow() -> bool {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var_os("LOSAT_X_TBNKEEPFEWSHADOW").is_some())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
    // ```c
    //     if (program_number != eBlastTypeRpsTblastn) {
    //         if (translation_buffer) {
    //             sfree(translation_buffer);
    //         }
    //     }
    // ```
    // NCBI frees each subject's translation after the subject, so a search holds the
    // translations of the subjects its threads are on. The kept plans are held where they are no
    // more than that.
    /// EXPERIMENT (LOSAT_X_TBNKEEPFEW): true when the plans of subjects of these lengths hold
    /// no more than `threads` plans of the longest one: no more than the threads of a search
    /// that keeps none can have in flight, each on one subject. A few long subjects (one
    /// genome, a few chromosomes) are kept, and batches searched at the same time share one
    /// plan of each instead of making one per batch. Many subjects searched one batch after
    /// another are not: without the plans the search held 115 MiB less on 28,036 16S subjects
    /// (40.6 Mb) and 309 MiB less on 57,775 zebrafish cDNAs (122 Mb), with five query batches,
    /// and took 0.980-0.994 of the time at 4, 8 and 16 threads (workstation, 2026-10-19:
    /// translating a packed subject again, `LOSAT_X_TBNXLATESIMD`, costs less than holding its
    /// plan; with `LOSAT_X_TBNFRAMEREUSE` 0.998 at one thread).
    pub(super) fn kept_within_threads(
        subject_lengths: impl IntoIterator<Item = usize>,
        threads: usize,
    ) -> bool {
        let (mut total, mut longest) = (0usize, 0usize);
        for length in subject_lengths {
            total = total.saturating_add(length);
            longest = longest.max(length);
        }
        total <= longest.saturating_mul(threads.max(1))
    }

    pub(super) fn new(subjects: usize, gencode: u8, mask_lowercase: bool) -> Self {
        Self {
            gencode,
            mask_lowercase,
            plans: (0..subjects).map(|_| std::sync::OnceLock::new()).collect(),
            kept_bytes: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// The subject plans of a search of `batches` query batches, or `None`
    /// when the search keeps none: `on` is false (`LOSAT_X_TBNKEEPFEW` keeps
    /// none, `kept_within_threads`), or there are
    /// fewer than `KEEP_FROM_BATCHES` batches. A kept plan is read by
    /// the batches after the first one that reaches its subject. A search of
    /// fewer batches makes its plans per use and drops them after the subject,
    /// like those above `KEEP_LIMIT`.
    pub(super) fn for_batches(
        on: bool,
        batches: usize,
        subjects: usize,
        gencode: u8,
        mask_lowercase: bool,
    ) -> Option<Self> {
        // The own side is recorded where a plan is kept or refused (KEEP_LIMIT).
        if on && batches < Self::KEEP_FROM_BATCHES {
            crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", false);
        }
        (on && batches >= Self::KEEP_FROM_BATCHES)
            .then(|| Self::new(subjects, gencode, mask_lowercase))
    }

    fn slot(
        &self,
        oid: usize,
        gencode: u8,
        mask_lowercase: bool,
    ) -> Option<&std::sync::OnceLock<Option<std::sync::Arc<XSubjectPlan>>>> {
        if self.gencode != gencode || self.mask_lowercase != mask_lowercase {
            return None;
        }
        self.plans.get(oid)
    }

    /// Translated residues of a plan: what `KEEP_LIMIT` counts.
    fn plan_bytes(plan: &XSubjectPlan) -> usize {
        plan.frames.iter().map(|frame| frame.aa_seq.len()).sum()
    }

    /// True when the plans of subjects of these lengths are all kept: a plan
    /// has six frames of at most `length / 3` residues and two sentinels, so
    /// `2 * length + 12` bounds `plan_bytes`.
    pub(super) fn all_fit(subject_lengths: impl IntoIterator<Item = usize>) -> bool {
        let mut total = 0usize;
        for length in subject_lengths {
            total = total.saturating_add(length.saturating_mul(2).saturating_add(12));
        }
        total <= Self::KEEP_LIMIT
    }

    /// Size of the plans held, for tests.
    #[cfg(test)]
    fn kept_bytes(&self) -> usize {
        self.kept_bytes.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Reserves room for `plan`. Called once per subject, from the
    /// initialiser of its cell. The total only changes when the plan fits, so
    /// a plan that does not fit never makes another one look too large.
    fn keep(&self, plan: &XSubjectPlan) -> bool {
        Self::reserve(&self.kept_bytes, Self::plan_bytes(plan), Self::KEEP_LIMIT)
    }

    /// Adds `bytes` to `kept` if the total stays within `limit`.
    fn reserve(kept: &std::sync::atomic::AtomicUsize, bytes: usize, limit: usize) -> bool {
        use std::sync::atomic::Ordering;
        kept.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |held| {
            held.checked_add(bytes).filter(|&total| total <= limit)
        })
        .is_ok()
    }

    /// The plan of the subject in `slot`: the kept one, made by the first
    /// caller, or a plan of the caller's own when it is not kept.
    fn plan(
        &self,
        slot: &std::sync::OnceLock<Option<std::sync::Arc<XSubjectPlan>>>,
        build: impl Fn() -> Result<std::sync::Arc<XSubjectPlan>>,
    ) -> Result<std::sync::Arc<XSubjectPlan>> {
        let mut failed = None;
        let mut own = None;
        let mut made = false;
        let kept = slot.get_or_init(|| {
            made = true;
            match build() {
                Ok(plan) if self.keep(&plan) => {
                    crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", true);
                    Some(plan)
                }
                Ok(plan) => {
                    // Over `KEEP_LIMIT`: the plan is made per use.
                    crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", false);
                    own = Some(plan);
                    None
                }
                Err(error) => {
                    failed = Some(error);
                    None
                }
            }
        });
        if let Some(error) = failed {
            return Err(error);
        }
        match (kept, own) {
            (Some(plan), _) => {
                if !made && Self::shadow() {
                    Self::x_compare_with_own(plan, build()?);
                }
                Ok(std::sync::Arc::clone(plan))
            }
            (None, Some(plan)) => Ok(plan),
            (None, None) => build(),
        }
    }

    // No NCBI counterpart: the comparison of the shadow mode; it does not change any value NCBI
    // computes.
    /// EXPERIMENT (LOSAT_X_TBNKEEPFEWSHADOW): stops the search when the plan made for a use
    /// differs from the kept plan that use reads; counts the comparison otherwise.
    fn x_compare_with_own(kept: &XSubjectPlan, own: std::sync::Arc<XSubjectPlan>) {
        if let Some(difference) = super::x_subject_translate::x_plan_fields_difference(&own, kept) {
            panic!(
                "LOSAT_X_TBNKEEPFEWSHADOW: the plan made for a use differs from the kept plan of the subject: {difference}"
            );
        }
        let residues: usize = kept.frames.iter().map(|frame| frame.aa_seq.len()).sum();
        X_KEEP_SHADOW_READS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        X_KEEP_SHADOW_RESIDUES.fetch_add(residues as u64, std::sync::atomic::Ordering::Relaxed);
    }
}

// No NCBI counterpart: counters of the shadow summary; they do not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNKEEPFEWSHADOW): the reads of kept plans compared with a plan made for
/// the use, and the translated residues in them (sentinels included).
static X_KEEP_SHADOW_READS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_KEEP_SHADOW_RESIDUES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNKEEPFEWSHADOW): the summary line, printed at exit when the shadow is
/// set (`x_subject_translate::print_summary`).
pub(super) fn x_keep_few_print_summary() {
    if !XSubjectSide::shadow() {
        return;
    }
    let reads = X_KEEP_SHADOW_READS.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!(
        "[X_SHADOW] TBNKEEPFEW kept_plan_reads_compared={reads} residues_compared={} ({})",
        X_KEEP_SHADOW_RESIDUES.load(std::sync::atomic::Ordering::Relaxed),
        if reads == 0 {
            "nothing compared: no kept plan was read again"
        } else {
            "frames and units identical"
        }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// ```c
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ...
//     for (context=first_context; context<=last_context; context++) {
//         BlastHSPList* hsp_list_for_chunks = NULL;
//
//         if (kTranslatedSubject) {
//             Uint4 i;
//             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//             subject->sequence = translation_buffer + frame_offsets[context] + 1;
//             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
// ```
// x_scan_plan does the setup that precedes NCBI's WordFinder calls for one subject (translation,
// scan ranges) and attaches the lookup table of the query set.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
// ```c
//        BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//                              &lookup_wrap->lut);
//        ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
//        BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                                  query, lookup_segments, 0);
// ```
/// `None` when the subject is too short to be searched.
pub(super) fn x_scan_plan(
    queries: &[&[u8]],
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subject's letters, or its packed store; (EXPERIMENT,
    // LOSAT_X_TBNDBSTREAM) or a `-db` sequence, read here when its plan is made.
    subject: super::x_subject_store::XPrelimSubject<'_>,
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_length: usize,
    // Where the table of this query set is kept (`XQuerySide`).
    x_lookup: &std::sync::OnceLock<
        std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>,
    >,
    // The search's subject plans and this subject's index among them.
    x_subject: Option<(&XSubjectSide, usize)>,
) -> Result<Option<XScanPlan>> {
    // An unknown genetic code is refused before the subject length is looked at.
    GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1334,1429-1432
    // return word_length * 3 + 2;
    // if (seq_arg.seq->length < min_subj_seq_length) { ... continue; }
    if subject.len() < word_length * 3 + 2 {
        return Ok(None);
    }
    // Everything below that reads the subject.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
    // ```c
    //             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
    //                                      backup.full_range.right,
    //                                      subject->gen_code_string, &translation_buffer,
    //                                      &frame_offsets, NULL);
    // ...
    //     for (context=first_context; context<=last_context; context++) {
    //         BlastHSPList* hsp_list_for_chunks = NULL;
    //
    //         if (kTranslatedSubject) {
    //             Uint4 i;
    //             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //             subject->sequence = translation_buffer + frame_offsets[context] + 1;
    //             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
    // ```
    // Translation and chunking of the subject (the shared plan memoises them).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1071-1074,1080-1090
    // ```c
    //    } else {
    //       translation_table = s_BlastGetTranslationTable(genetic_code, FALSE);
    //       translation_table_rc = s_BlastGetTranslationTable(genetic_code, TRUE);
    //    }
    // ...
    //    for (context = 0; context < NUM_FRAMES; ++context) {
    //       frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //       if (encoding == eBlastEncodingNcbi2na) {
    //          if (frame > 0) {
    //             length =
    //                BLAST_TranslateCompressedSequence(translation_table,
    //                   nucl_length, nucl_seq, frame, translation_buffer+offset);
    //          } else {
    //             length =
    //                BLAST_TranslateCompressedSequence(translation_table_rc,
    //                   nucl_length, nucl_seq, frame, translation_buffer+offset);
    // ```
    // What the shared plan memoises: the subject becomes ncbi2na in one pass and each frame is read
    // from the codon table of its strand, as BLAST_GetAllTranslations does for the compressed
    // subject (tblastn/x_subject_translate.rs). Graduated 2026-10-22 (was LOSAT_X_TBNXLATE).
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject is translated from its packed bases
    // (x_packed_subject_plan); a stage that reads the case of a subject with lower-case letters
    // gets the letters instead. The letters of a subject searched from its letters (the records
    // `run_local` is given, or a set that keeps its letters for a byte without an ncbi4na code)
    // are packed for the plan (`x_plan_store`), and the plan is made from them as for a packed
    // subject.
    let x_build_subject = || -> Result<std::sync::Arc<XSubjectPlan>> {
        use super::x_subject_store::XPrelimSubject;
        use crate::sequence::packed_subject::SubjectNt;
        match subject {
            // EXPERIMENT (LOSAT_X_TBNDBSTREAM): the plan is made from the sequence's bytes, read
            // for it (a database sequence has no lower-case letters).
            XPrelimSubject::Stream(stream, oid) => stream.plan_with(oid, |packed, length| {
                super::x_subject_translate::x_packed_bytes_plan(packed, length, db_gencode)
            }),
            XPrelimSubject::Nt(SubjectNt::Letters(letters)) => {
                match super::x_subject_store::x_plan_store(letters, mask_lowercase) {
                    Some(store) => {
                        super::x_subject_translate::x_packed_subject_plan(&store, db_gencode)
                    }
                    None => super::x_subject_translate::x_subject_plan(
                        letters,
                        db_gencode,
                        mask_lowercase,
                    ),
                }
            }
            XPrelimSubject::Nt(SubjectNt::Packed(store)) => {
                anyhow::ensure!(
                    !mask_lowercase || store.lowercase().is_empty(),
                    "TBLASTN: a packed subject with lower-case letters searched with lowercase masking"
                );
                super::x_subject_translate::x_packed_subject_plan(store, db_gencode)
            }
        }
    };
    let x_slot = x_subject.and_then(|(side, oid)| {
        side.slot(oid, db_gencode, mask_lowercase)
            .map(|slot| (side, slot))
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
    // ```c
    //             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
    //                                      backup.full_range.right,
    //                                      subject->gen_code_string, &translation_buffer,
    //                                      &frame_offsets, NULL);
    // ...
    //     for (context=first_context; context<=last_context; context++) {
    //         BlastHSPList* hsp_list_for_chunks = NULL;
    //
    //         if (kTranslatedSubject) {
    //             Uint4 i;
    //             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //             subject->sequence = translation_buffer + frame_offsets[context] + 1;
    //             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
    // ```
    // Dispatch point: with the search's subject plans the plan is the shared one (made by the
    // first batch that reaches the subject); a search of fewer than four batches makes it here.
    let subject_plan = match x_slot {
        Some((side, slot)) => side.plan(slot, x_build_subject)?,
        None => x_build_subject()?,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
    // ```c
    //        BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
    //                              &lookup_wrap->lut);
    //        ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
    //        BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
    //                                  query, lookup_segments, 0);
    // ```
    // The query-side lookup table: BlastAaLookupTableNew and BlastAaLookupIndexQuery on the masked
    // query. The table-only variant (`x_lookup_table_for_profile`, graduated 2026-10-22, was
    // LOSAT_X_LUTSPLIT) skips the per-context Karlin blocks, which the table does not read.
    // NCBI c++/src/algo/blast/core/blast_aalookup.c:227-245:
    // lookup->word_length = opt->word_size;
    // lookup->threshold = (Int4)opt->threshold;
    let x_build_lookup = || {
        let query_frames: Vec<_> = queries
            .iter()
            // NCBI blast_filter.c:1241-1255 merges lowercase query and SEG masks.
            .map(|query| vec![encode_tblastn_lookup_query(query, seg, mask_lowercase)])
            .collect();
        std::sync::Arc::new(
            crate::algorithm::tblastx::lookup::x_lookup_table_for_profile(
                &query_frames,
                threshold,
                matrix,
                word_length,
            ),
        )
    };
    // The table depends on the query set and the word-finder options only; the first subject that
    // gets here builds it (`XQuerySide`, was LOSAT_X_TBNQSIDE).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
    // ```c
    //        BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
    //                              &lookup_wrap->lut);
    //        ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
    //        BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
    //                                  query, lookup_segments, 0);
    // ```
    // NCBI builds the lookup table once per query batch: the first subject builds it and the
    // others read it.
    let lookup = std::sync::Arc::clone(x_lookup.get_or_init(x_build_lookup));
    let pair_capacity = (lookup.longest_chain.max(1) as usize) * 1024;
    i32::try_from(pair_capacity).expect("NCBI offset array fits Int4");
    Ok(Some(XScanPlan {
        lookup,
        subject: subject_plan,
        pair_capacity,
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
// ```c
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
// ...
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
// ```
// x_scan_unit runs the scansub() calls of one unit and passes every filled offset array to the
// caller in scan order, as the loop at line 506 consumes it.
/// Scan one unit and hand every filled offset
/// array to `on_pairs` in scan order, as NCBI c++/src/algo/blast/core/
/// aa_ungapped.c:492-516 consumes it (`for (i = 0; i < hits; ++i)`).
pub(super) fn x_scan_unit(
    plan: &XScanPlan,
    unit: usize,
    pairs: &mut Vec<BlastOffsetPair>,
    mut on_pairs: impl FnMut(&[BlastOffsetPair]),
) {
    let (frame_index, chunk, ranges) = &plan.subject.units[unit];
    let frame = &plan.subject.frames[*frame_index];
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1018-1019
    // BlastChooseProteinScanSubject(lookup_wrap);
    // aux_struct->WordFinder = BlastAaWordFinder;
    // The scan variant of the unit is chosen once per table, as NCBI picks one scansub.
    //
    // EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan on presence bytes, for a table of tblastn's shape
    // (without AVX2 a sparse one); it writes up to X_SCAN_BYTE_SLACK entries after the capacity.
    let scan_byte = super::x_scan_byte::x_scan_byte();
    let presence = super::x_scan_byte::x_presence(scan_byte, &plan.lookup);
    // EXPERIMENT (LOSAT_X_TBNSCANBF): the branch-free scan writes up to X_SCAN_BF_SLACK entries
    // after the capacity.
    let scan_bf = super::x_scan_bf::x_tbn_scan_bf();
    let pairs_len = if presence.is_some() {
        plan.pair_capacity + super::x_scan_byte::X_SCAN_BYTE_SLACK
    } else if scan_bf {
        plan.pair_capacity + super::x_scan_bf::X_SCAN_BF_SLACK
    } else {
        plan.pair_capacity
    };
    if pairs.len() < pairs_len {
        pairs.resize(pairs_len, BlastOffsetPair::default());
    }
    let pair_capacity = plan.pair_capacity as i32;
    // EXPERIMENT (LOSAT_X_TBNSCANBF): the scan calls of the unit, for LOSAT_X_STATS.
    let mut x_calls = 0u64;
    for (range_index, &(left, right)) in ranges.iter().enumerate() {
        // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:496-501
        // scan_range[1] = subject->seq_ranges[0].left;
        // scan_range[2] = subject->seq_ranges[0].right - wordsize;
        // if (scan_range[2] < scan_range[1])
        //     scan_range[2] = scan_range[1];
        let mut end = right - plan.lookup.word_length as i32;
        if range_index == 0 && end < left {
            end = left;
        }
        let mut scan_range = [0, left, end];
        while scan_range[1] <= scan_range[2] {
            let hits = if let Some(presence) = &presence {
                super::x_scan_byte::x_scan_one_range_byte(
                    scan_byte,
                    &plan.lookup,
                    presence,
                    &frame.aa_seq[1 + chunk.offset..],
                    pairs,
                    pair_capacity,
                    &mut scan_range,
                )
            } else if scan_bf {
                super::x_scan_bf::x_scan_one_range_bf(
                    &plan.lookup,
                    &frame.aa_seq[1 + chunk.offset..],
                    pairs,
                    pair_capacity,
                    &mut scan_range,
                )
            } else {
                s_blast_aa_scan_subject_one_range(
                    &plan.lookup,
                    &frame.aa_seq[1 + chunk.offset..],
                    pairs,
                    pair_capacity,
                    &mut scan_range,
                )
            };
            on_pairs(&pairs[..hits as usize]);
            x_calls += 1;
        }
    }
    if presence.is_some() {
        super::x_scan_byte::x_count_unit(x_calls);
    } else if super::x_scan_bf::x_stats() {
        super::x_scan_bf::x_count_unit(scan_bf, x_calls);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::encoding::ncbistdaa_to_ascii;
    use crate::algorithm::tblastn::search_gapped::{preliminary_stages, scan_profile};
    use crate::sequence::packed_subject::SubjectNt;
    use std::collections::HashMap;
    use std::fs;

    fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
        let text = fs::read_to_string(path).unwrap();
        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        for line in text.lines() {
            if let Some(id) = line.strip_prefix('>') {
                records.push((
                    id.split_whitespace().next().unwrap().to_string(),
                    Vec::new(),
                ));
            } else {
                records.last_mut().unwrap().1.extend(line.as_bytes());
            }
        }
        records
    }

    // NCBI c++/src/algo/blast/core/blast_query_info.c:68-96:
    // one context per protein query; query offsets advance by length + 1.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505:
    // scansub emits the ordered global offset pairs for each subject frame.
    // This comparison uses the unmodified candidate stream from the pinned
    // local -subject oracle; no final-output filter is applied to candidates.
    #[test]
    fn multi_query_candidate_contexts_and_order_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let frames: Vec<_> = refs
            .iter()
            .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
            .collect();
        let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (_, contexts) = build_ncbi_lookup(&frames, 13, &karlin, false);
        let actual_contexts: Vec<_> = contexts
            .iter()
            .map(|context| {
                (
                    context.q_idx as usize,
                    context.frame_base,
                    context.aa_len,
                    context.is_valid,
                )
            })
            .collect();
        let trace = fs::read_to_string(format!("{root}query_context.tsv")).unwrap();
        let expected_contexts: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[7] == "1",
                )
            })
            .collect();
        assert_eq!(actual_contexts, expected_contexts);

        let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let frame_order = [1, 2, 3, -1, -2, -3];
        let expected: Vec<Seed> = trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                Seed {
                    frame: frame_order[f[1].parse::<usize>().unwrap()],
                    chunk_offset: 0,
                    query_offset: f[3].parse().unwrap(),
                    subject_offset: f[4].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &refs,
            subject,
            1,
            scan_profile(refs.len(), None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "first multi-query candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first),
                expected.get(first),
                actual.len(),
                expected.len()
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_aalookup.c:227-245,438-465:
    // word_length and neighbor scores come from the selected options/matrix.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505:
    // the subject scanner emits global query offsets in its observed order.
    #[test]
    fn alternate_matrix_word2_candidates_and_contexts_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let frames: Vec<_> = refs
            .iter()
            .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
            .collect();
        let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum45);
        let (lookup, contexts) =
            build_ncbi_lookup_for_profile(&frames, 16, &karlin, false, ScoringMatrix::Blosum45, 2);
        assert_eq!(lookup.word_length, 2);
        let trace = fs::read_to_string(format!("{root}query_context.tsv")).unwrap();
        let expected_contexts: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse::<usize>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<usize>().unwrap(),
                    f[7] == "1",
                )
            })
            .collect();
        let actual_contexts: Vec<_> = contexts
            .iter()
            .map(|ctx| (ctx.q_idx as usize, ctx.frame_base, ctx.aa_len, ctx.is_valid))
            .collect();
        assert_eq!(actual_contexts, expected_contexts);
        let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let frame_order = [1, 2, 3, -1, -2, -3];
        let expected: Vec<Seed> = trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                Seed {
                    frame: frame_order[f[1].parse::<usize>().unwrap()],
                    chunk_offset: 0,
                    query_offset: f[3].parse().unwrap(),
                    subject_offset: f[4].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &refs,
            subject,
            1,
            scan_profile(refs.len(), None, 16, false, ScoringMatrix::Blosum45, 2),
        )
        .unwrap()
        .seeds;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!("first BLOSUM45/word-2 candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first), expected.get(first), actual.len(), expected.len());
        }
    }

    // NCBI reference: c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
    // CRandom random(base_length); ambiguous ncbi4na bases are picked once
    // for the compressed ncbi2na preliminary subject.
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:767-775,804-812
    // BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na, ...);
    // subject->sequence = translation_buffer + frame_offsets[context] + 1;
    #[test]
    fn ncbi_preliminary_subject_frames_match_all_fixture_bytes() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let records = read_fasta(&format!("{root}subjects.fna"));
        let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let mut plans = HashMap::new();
        let trace = fs::read_to_string(format!("{root}ncbi_frame_trace.tsv")).unwrap();
        let mut compared = 0;
        for line in trace.lines().skip(1) {
            let columns: Vec<_> = line.split('\t').collect();
            let subject_id = columns[1];
            let frame: i8 = columns[2].parse().unwrap();
            // The frames the search translates: the unit path's plan.
            let plan = plans.entry(subject_id).or_insert_with(|| {
                x_scan_plan(
                    &[query],
                    SubjectNt::Letters(by_id[subject_id]).into(),
                    1,
                    None,
                    13,
                    false,
                    ScoringMatrix::Blosum62,
                    3,
                    &std::sync::OnceLock::new(),
                    None,
                )
                .unwrap()
                .unwrap()
            });
            let actual = plan
                .subject
                .frames
                .iter()
                .find(|f| f.frame == frame)
                .unwrap();
            let actual_hex: String = actual.aa_seq[1..actual.aa_seq.len() - 1]
                .iter()
                .map(|base| format!("{base:02x}"))
                .collect();
            assert_eq!(
                actual.aa_len,
                columns[3].parse::<usize>().unwrap(),
                "{subject_id}/{frame}"
            );
            assert_eq!(actual_hex, columns[4], "{subject_id}/{frame}");
            compared += 1;
        }
        assert_eq!(compared, 78);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:478-505
    // scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    // for (i = 0; i < hits; ++i) { query_offset = offset_pairs[i].qs_offsets.q_off;
    //                            subject_offset = offset_pairs[i].qs_offsets.s_off; }
    #[test]
    fn ncbi_all_candidates_match_in_set_and_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let records = read_fasta(&format!("{root}subjects.fna"));
        let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                (
                    fields[1].to_string(),
                    fields[2].parse::<i8>().unwrap(),
                    fields[4].parse::<u32>().unwrap(),
                    fields[5].parse::<u32>().unwrap(),
                )
            })
            .collect();
        let mut actual = Vec::new();
        for (id, subject) in &records {
            let seeds = preliminary_stages(
                &[query],
                subject,
                1,
                scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds;
            for seed in seeds {
                actual.push((
                    id.clone(),
                    seed.frame,
                    seed.query_offset,
                    seed.subject_offset,
                ));
            }
        }
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
    // if (context == 0) { /* convert subject scan ranges to AA */ }
    // else if (context == 3) { /* reverse ranges for minus frames */ }
    // NCBI reference: c++/src/algo/blast/core/masksubj.inl:51-71
    // range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
    // range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
    #[test]
    fn ncbi_lowercase_mask_candidates_match_in_set_and_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/lowercase_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let records = read_fasta(&format!("{root}subjects.fna"));
        let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].to_string(),
                    f[2].parse::<i8>().unwrap(),
                    f[4].parse::<u32>().unwrap(),
                    f[5].parse::<u32>().unwrap(),
                )
            })
            .collect();
        let mut actual = Vec::new();
        for (id, subject) in &records {
            for seed in preliminary_stages(
                &[query],
                subject,
                1,
                scan_profile(1, None, 13, true, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds
            {
                actual.push((
                    id.clone(),
                    seed.frame,
                    seed.query_offset,
                    seed.subject_offset,
                ));
            }
        }
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!("first lowercase candidate mismatch at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                   actual.get(first), expected.get(first), actual.len(), expected.len());
        }
    }

    // NCBI reference: c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
    // CRandom random(base_length); each ambiguous 4na value consumes one GetRand().
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    #[test]
    fn ncbi_ambiguity_spectrum_frames_and_candidates_match() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        // The frames the search translates: the unit path's plan.
        let plan = x_scan_plan(
            &[query],
            SubjectNt::Letters(subject).into(),
            1,
            None,
            13,
            false,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            None,
        )
        .unwrap()
        .unwrap();
        let actual_frames = &plan.subject.frames;
        let frame_trace = fs::read_to_string(format!("{root}ncbi_frame_trace.tsv")).unwrap();
        for line in frame_trace.lines().skip(1) {
            let f: Vec<_> = line.split('\t').collect();
            let frame: i8 = f[2].parse().unwrap();
            let actual = actual_frames.iter().find(|x| x.frame == frame).unwrap();
            let hex: String = actual.aa_seq[1..actual.aa_seq.len() - 1]
                .iter()
                .map(|base| format!("{base:02x}"))
                .collect();
            assert_eq!(hex, f[4], "frame {frame}");
        }
        let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[2].parse::<i8>().unwrap(),
                    f[4].parse::<u32>().unwrap(),
                    f[5].parse::<u32>().unwrap(),
                )
            })
            .collect();
        let actual: Vec<_> = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds
        .into_iter()
        .map(|x| (x.frame, x.query_offset, x.subject_offset))
        .collect();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 181);
    }

    // NCBI reference: c++/src/algo/blast/core/blast_util.c:1080-1101
    // frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    // offset += length + 1;
    // frame_offsets[context+1] = offset;
    // New NCBI -outfmt 6 output supplies the frame and aligned subject AA.
    #[test]
    fn six_frame_translation_and_seed_coverage() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let records = read_fasta(&format!("{root}subjects.fna"));
        let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
        let output = fs::read_to_string(format!("{root}raw_isolation_fields.out")).unwrap();
        let mut frames_checked = 0;
        for line in output.lines() {
            let fields: Vec<_> = line.split('\t').collect();
            let subject_id = fields[1];
            if !(subject_id.starts_with("plus") || subject_id.starts_with("minus")) {
                continue;
            }
            let expected_frame: i8 = fields[7].parse().unwrap();
            let subject = by_id[subject_id];
            // The frames the search translates: the unit path's plan.
            let plan = x_scan_plan(
                &[query],
                SubjectNt::Letters(subject).into(),
                1,
                None,
                13,
                false,
                ScoringMatrix::Blosum62,
                3,
                &std::sync::OnceLock::new(),
                None,
            )
            .unwrap()
            .unwrap();
            let frame = plan
                .subject
                .frames
                .iter()
                .find(|f| f.frame == expected_frame)
                .unwrap();
            let translated: String = frame.aa_seq[1..frame.aa_seq.len() - 1]
                .iter()
                .map(|&aa| ncbistdaa_to_ascii(aa))
                .collect();
            assert!(translated.contains(fields[9]), "{subject_id}: {translated}");
            let seeds = preliminary_stages(
                &[query],
                subject,
                1,
                scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds;
            assert!(
                seeds.iter().any(|seed| seed.frame == expected_frame
                    && seed.query_offset == 0
                    && seed.subject_offset == 0),
                "{subject_id}"
            );
            frames_checked += 1;
        }
        assert_eq!(frames_checked, 6);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234
    // status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
    // Blast_InitHitListSortByScore(init_hitlist);
    // The comparison-only C probe captures NCBI's saved seed pair. Each
    // unambiguous saved pair must exist in the Rust scan for the same frame.
    #[test]
    fn ncbi_saved_init_seeds_exist_in_rust_word_scan() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let records = read_fasta(&format!("{root}subjects.fna"));
        let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
        let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
        let mut compared = 0;
        for line in trace.lines().skip(1) {
            let fields: Vec<_> = line.split('\t').collect();
            let subject_id = fields[1];
            if subject_id == "ambiguous" {
                continue;
            }
            let frame: i8 = fields[2].parse().unwrap();
            let q_seed: u32 = fields[4].parse().unwrap();
            let s_seed: u32 = fields[5].parse().unwrap();
            let seeds = preliminary_stages(
                &[query],
                by_id[subject_id],
                1,
                scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds;
            assert!(
                seeds.iter().any(|seed| seed.frame == frame
                    && seed.query_offset == q_seed
                    && seed.subject_offset == s_seed),
                "{subject_id}: ({q_seed},{s_seed})"
            );
            compared += 1;
        }
        assert_eq!(compared, 10);
    }

    // NCBI reference: c++/src/algo/blast/core/blast_filter.c:337-370
    // BlastSetUp_Filter(..., query_blk, ...);
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    #[test]
    fn seg_suppresses_low_complexity_query_seeds() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query_low.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}low_subject.fna"))[0].1;
        let unmasked = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        let masked = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(
                1,
                Some(&SegParams::default()),
                13,
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .seeds;
        println!(
            "fixture=low_complexity seg_off_seeds={} seg_on_seeds={}",
            unmasked.len(),
            masked.len()
        );
        assert!(!unmasked.is_empty());
        assert!(masked.is_empty());
        assert!(fs::read(format!("{root}low_query_seg_on.out"))
            .unwrap()
            .is_empty());
        assert!(!fs::read(format!("{root}low_query_seg_off.out"))
            .unwrap()
            .is_empty());
    }

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:747-775,804-841
    // BLAST_GetAllTranslations(..., subject->gen_code_string, ...);
    // subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:478-505
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    #[test]
    fn code32_api_oracle_candidates_match_in_set_and_order() {
        let fixture = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let trace_root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/code32_20260923/"
        );
        let query = &read_fasta(&format!("{fixture}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{fixture}subject_code32.fna"))[0].1;
        let trace = fs::read_to_string(format!("{trace_root}ncbi_candidates.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[0].parse::<i8>().unwrap(),
                    f[2].parse::<u32>().unwrap(),
                    f[3].parse::<u32>().unwrap(),
                )
            })
            .collect();
        let actual: Vec<_> = preliminary_stages(
            &[query],
            subject,
            32,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds
        .into_iter()
        .map(|x| (x.frame, x.query_offset, x.subject_offset))
        .collect();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 199);
    }

    // NCBI reference: c++/src/objects/seqfeat/gc.prt:340-347
    // id 32 , ncbieaa "FFLLSSSSYYWWCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    #[test]
    fn code_32_changes_search_words_and_ambiguity_resolves() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subject_code32.fna"))[0].1;
        let code1 = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        let code32 = preliminary_stages(
            &[query],
            subject,
            32,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        println!(
            "fixture=code32 code1_seeds={} code32_seeds={}",
            code1.len(),
            code32.len()
        );
        assert_ne!(code1, code32);
        assert!(code32
            .iter()
            .any(|seed| seed.frame == 1 && seed.query_offset == 0 && seed.subject_offset == 0));
        assert!(preliminary_stages(
            &[query],
            b"ATN",
            32,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3)
        )
        .is_ok());
        assert!(preliminary_stages(
            &[query],
            b"atg",
            32,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3)
        )
        .is_ok());
    }
    // NCBI c++/src/algo/blast/core/blast_query_info.c:68-96:
    // each protein query receives one context; offsets advance length + 1,
    // and the X-only third context is invalid for the lookup.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505 and
    // c++/src/algo/blast/core/blast_engine.c:804-841:
    // scansub emits every pair in frame/chunk order over both long chunks.
    #[test]
    fn long_multi_query_every_candidate_and_context_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let fixture = format!("{root}long_multi_query_20260924/");
        let queries = read_fasta(&format!("{fixture}query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let plus1 = read_fasta(&format!("{root}run_20260923/subjects.fna"))
            .into_iter()
            .find(|(id, _)| id == "plus1")
            .unwrap()
            .1;
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(&plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        let frames: Vec<_> = refs
            .iter()
            .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
            .collect();
        let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (_, contexts) = build_ncbi_lookup(&frames, 13, &karlin, false);
        let actual_contexts: Vec<_> = contexts
            .iter()
            .map(|c| (c.q_idx as usize, c.frame_base, c.aa_len, c.is_valid))
            .collect();
        let gapped_trace = fs::read_to_string(format!("{fixture}gapped_events.tsv")).unwrap();
        let expected_contexts: Vec<_> = gapped_trace
            .lines()
            .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[7] == "1",
                )
            })
            .collect();
        assert_eq!(actual_contexts, expected_contexts);
        assert_eq!(
            actual_contexts,
            [(0, 0, 120, true), (1, 121, 70, true), (2, 192, 120, false)]
        );
        let trace = fs::read_to_string(format!("{fixture}candidate.stderr")).unwrap();
        assert_eq!(
            trace.lines().filter(|l| l.starts_with("PARAM\t")).count(),
            12
        );
        assert_eq!(trace.lines().filter(|l| l.starts_with("END\t")).count(), 12);
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505:
        // scansub reads the same per-call scan
        // parameters and chunk state before returning each ordered pair.
        let params: Vec<_> = trace
            .lines()
            .filter(|l| l.starts_with("PARAM\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(params.len(), 12);
        for (call, row) in params.iter().enumerate() {
            assert_eq!(*row, (call, 16, 16, 28));
        }
        let frame_order = [1, 2, 3, -1, -2, -3];
        let chunk_rows = fs::read_to_string(format!("{fixture}frame_chunks.tsv")).unwrap();
        let observed_chunks: Vec<_> = chunk_rows
            .lines()
            .filter(|line| line.starts_with("FRAME_CHUNK\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<i8>().unwrap(),
                    f[3].parse::<usize>().unwrap(),
                )
            })
            .collect();
        assert_eq!(observed_chunks.len(), 12);
        for (call, row) in observed_chunks.iter().enumerate() {
            assert_eq!(
                *row,
                (
                    call,
                    frame_order[call / 2],
                    if call % 2 == 0 { 5_000_000 } else { 420 }
                )
            );
        }
        let mut candidate_counts = [0usize; 12];
        let expected: Vec<_> = trace
            .lines()
            .filter(|l| l.starts_with("CAND\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                candidate_counts[call] += 1;
                Seed {
                    frame: frame_order[call / 2],
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    query_offset: f[3].parse().unwrap(),
                    subject_offset: f[4].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(expected.len(), 394);
        let ends: Vec<_> = trace
            .lines()
            .filter(|l| l.starts_with("END\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<usize>().unwrap(),
                )
            })
            .collect();
        assert_eq!(ends.len(), 12);
        for (call, row) in ends.iter().enumerate() {
            assert_eq!(*row, (call, candidate_counts[call]));
        }
        assert!(expected.iter().all(|seed| seed.query_offset < 192));
        let actual = preliminary_stages(
            &refs,
            &subject,
            1,
            scan_profile(refs.len(), None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!("first long multi-query candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                   actual.get(first), expected.get(first), actual.len(), expected.len());
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:804-841:
    // for (context=first_context; context<=last_context; context++) {
    //     status = s_BlastSearchEngineOneContext(...);
    // }
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:495-505:
    // scan_range[1] = subject->seq_ranges[0].left;
    // scan_range[2] = subject->seq_ranges[0].right - wordsize;
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    #[test]
    fn long_subject_every_candidate_matches_ncbi_in_chunk_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let frame_order = [1, 2, 3, -1, -2, -3];
        for (fixture, mask_lowercase) in [
            ("long_chunk_20260924", false),
            ("masked_chunk_boundary_20260924", true),
        ] {
            let query = &read_fasta(&format!("{root}{fixture}/query.faa"))[0].1;
            let mut subject = Vec::with_capacity(15_000_962);
            for _ in 0..4_999_950 {
                subject.extend_from_slice(b"ATG");
            }
            subject.extend_from_slice(plus1);
            for _ in 0..250 {
                subject.extend_from_slice(b"ATG");
            }
            assert_eq!(subject.len(), 15_000_962);
            if mask_lowercase {
                subject[14_999_700..15_000_000].make_ascii_lowercase();
            }
            let trace = fs::read_to_string(format!(
                "{root}ordered_candidates_20260924/{fixture}.stderr"
            ))
            .unwrap();
            assert_eq!(
                trace
                    .lines()
                    .filter(|line| line.starts_with("PARAM\t"))
                    .count(),
                12
            );
            assert_eq!(
                trace
                    .lines()
                    .filter(|line| line.starts_with("END\t"))
                    .count(),
                12
            );
            let expected: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("CAND\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    let call: usize = f[1].parse().unwrap();
                    Seed {
                        frame: frame_order[call / 2],
                        chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                        query_offset: f[3].parse().unwrap(),
                        subject_offset: f[4].parse().unwrap(),
                    }
                })
                .collect();
            let actual = preliminary_stages(
                &[query.as_slice()],
                &subject,
                1,
                scan_profile(1, None, 13, mask_lowercase, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds;
            if actual != expected {
                let first = actual
                    .iter()
                    .zip(&expected)
                    .position(|(left, right)| left != right)
                    .unwrap_or(actual.len().min(expected.len()));
                panic!(
                    "{fixture}: first candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                    actual.get(first), expected.get(first), actual.len(), expected.len()
                );
            }
        }
        // NCBI blast_engine.c:246-264,478-500,804-841:
        // s_GetNextSubjectChunk skips a middle chunk with no soft range;
        // each surviving chunk calls WordFinder in frame order.
        for (fixture, early_hit) in [
            ("no_range_middle_20260924", false),
            ("no_range_two_hits_20260924", true),
        ] {
            let query = &read_fasta(&format!("{root}{fixture}/query.faa"))[0].1;
            let mut subject = if early_hit {
                let mut nt = b"ATG".repeat(4_998_450);
                nt.extend_from_slice(plus1);
                nt.push(b'A');
                nt.extend_from_slice(&b"ATG".repeat(10_000_050 - 4_998_450 - 121));
                nt.extend_from_slice(plus1);
                nt
            } else {
                let mut nt = b"ATG".repeat(10_000_050);
                nt.extend_from_slice(plus1);
                nt
            };
            subject.extend_from_slice(&b"ATG".repeat(250));
            assert_eq!(subject.len(), 30_001_262);
            subject[14_997_000..30_000_000].make_ascii_lowercase();
            let trace = fs::read_to_string(format!(
                "{root}ordered_candidates_20260924/{fixture}.stderr"
            ))
            .unwrap();
            let frames = [1, 1, 2, 2, 3, 3, -1, -1, -1, -2, -2, -2, -3, -3, -3];
            let offsets = [
                0, 9_999_800, 0, 9_999_800, 0, 9_999_800, 0, 4_999_900, 9_999_800, 0, 4_999_900,
                9_999_800, 0, 4_999_900, 9_999_800,
            ];
            let range_trace =
                fs::read_to_string(format!("{root}{fixture}/chunk_ranges.tsv")).unwrap();
            let calls: Vec<_> = range_trace
                .lines()
                .filter(|line| line.starts_with("RANGE_CALL\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (f[2].parse::<i8>().unwrap(), f[3].parse::<usize>().unwrap())
                })
                .collect();
            assert_eq!(calls.len(), 15);
            for (call, (frame, length)) in calls.iter().enumerate() {
                assert_eq!(*frame, frames[call]);
                assert_eq!(
                    *length,
                    if offsets[call] == 9_999_800 {
                        620
                    } else {
                        5_000_000
                    }
                );
            }
            assert_eq!(
                trace
                    .lines()
                    .filter(|line| line.starts_with("PARAM\t"))
                    .count(),
                15
            );
            assert_eq!(
                trace
                    .lines()
                    .filter(|line| line.starts_with("END\t"))
                    .count(),
                15
            );
            let expected: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("CAND\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    let call: usize = f[1].parse().unwrap();
                    Seed {
                        frame: frames[call],
                        chunk_offset: offsets[call],
                        query_offset: f[3].parse().unwrap(),
                        subject_offset: f[4].parse().unwrap(),
                    }
                })
                .collect();
            let actual = preliminary_stages(
                &[query.as_slice()],
                &subject,
                1,
                scan_profile(1, None, 13, true, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .seeds;
            if actual != expected {
                let first = actual
                    .iter()
                    .zip(&expected)
                    .position(|(left, right)| left != right)
                    .unwrap_or(actual.len().min(expected.len()));
                panic!(
                    "{fixture}: first candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                    actual.get(first), expected.get(first), actual.len(), expected.len()
                );
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
    // split only when offset + MAX_DBSEQ_LEN < hard-range right;
    // the next offset overlaps the first chunk by DBSEQ_CHUNK_OVERLAP.
    #[test]
    fn unmasked_long_translated_chunks_match_ncbi_wordfinder_calls() {
        let rows = unmasked_translated_chunks(5_000_320);
        assert_eq!(
            rows,
            vec![
                TranslatedChunk {
                    offset: 0,
                    length: 5_000_000
                },
                TranslatedChunk {
                    offset: 4_999_900,
                    length: 420
                },
            ]
        );
        assert_eq!(
            unmasked_translated_chunks(5_000_000),
            vec![TranslatedChunk {
                offset: 0,
                length: 5_000_000
            },]
        );
        assert_eq!(
            unmasked_translated_chunks(5_000_001),
            vec![
                TranslatedChunk {
                    offset: 0,
                    length: 5_000_000
                },
                TranslatedChunk {
                    offset: 4_999_900,
                    length: 101
                },
            ]
        );
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/long_chunk_20260924/frame_chunks.tsv"
        );
        let trace = fs::read_to_string(path).unwrap();
        let observed: Vec<(i8, usize)> = trace
            .lines()
            .filter(|line| line.starts_with("FRAME_CHUNK\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (f[2].parse().unwrap(), f[3].parse().unwrap())
            })
            .collect();
        let expected: Vec<_> = [1, 2, 3, -1, -2, -3]
            .into_iter()
            .flat_map(|frame| rows.iter().map(move |chunk| (frame, chunk.length)))
            .collect();
        assert_eq!(observed, expected);
    }

    // No NCBI counterpart: test helper (pseudo-random nucleotide subject); it does not change any
    // value NCBI computes.
    fn x_test_subject(length: usize) -> Vec<u8> {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect()
    }

    // No NCBI counterpart: test of the plan cache (one plan, charged once, when batches race); it
    // does not change any value NCBI computes.
    // Subject plans: batches that reach a subject at the same time get one
    // plan. It is made once, and the shared limit is charged for it once.
    #[test]
    fn x_subject_plan_is_made_and_charged_once_when_batches_race() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let subject = x_test_subject(9_000);
        let code = GeneticCode::try_from_id(1).unwrap();
        let side = XSubjectSide::new(3, 1, false);
        let slot = side.slot(1, 1, false).unwrap();
        let builds = AtomicUsize::new(0);
        let build = || -> Result<std::sync::Arc<XSubjectPlan>> {
            builds.fetch_add(1, Ordering::Relaxed);
            // Long enough for every thread to arrive while the first builds.
            std::thread::sleep(std::time::Duration::from_millis(50));
            Ok(std::sync::Arc::new(XSubjectPlan {
                frames: generate_frames(&subject, &code),
                units: Vec::new(),
            }))
        };
        let threads = 8;
        let barrier = std::sync::Barrier::new(threads);
        let plans: Vec<_> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    scope.spawn(|| {
                        barrier.wait();
                        side.plan(slot, &build).unwrap()
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        assert_eq!(builds.load(Ordering::Relaxed), 1);
        for plan in &plans[1..] {
            assert!(std::sync::Arc::ptr_eq(&plans[0], plan));
        }
        let bytes = XSubjectSide::plan_bytes(&plans[0]);
        // Three frames per strand of 3000, 2999 and 2999 residues, two sentinels each.
        assert_eq!(bytes, 2 * (3_000 + 2_999 + 2_999) + 12);
        assert_eq!(side.kept_bytes(), bytes);
        // A later batch reads the same plan and adds nothing.
        let again = side.plan(slot, &build).unwrap();
        assert!(std::sync::Arc::ptr_eq(&plans[0], &again));
        assert_eq!(
            (builds.load(Ordering::Relaxed), side.kept_bytes()),
            (1, bytes)
        );
        // Another gencode or masking option does not use these plans.
        assert!(side.slot(1, 4, false).is_none() && side.slot(1, 1, true).is_none());
    }

    // No NCBI counterpart: test of the plan cache through x_scan_plan; it does not change any value
    // NCBI computes.
    // The whole stage: `x_scan_plan` from eight threads on one subject.
    #[test]
    fn x_scan_plans_of_racing_batches_share_one_subject_plan() {
        let subject = x_test_subject(30_000);
        let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
        let side = XSubjectSide::new(1, 1, false);
        let threads = 8;
        let barrier = std::sync::Barrier::new(threads);
        let plans: Vec<_> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    scope.spawn(|| {
                        barrier.wait();
                        x_scan_plan(
                            &[query],
                            crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
                            1,
                            None,
                            13,
                            false,
                            ScoringMatrix::Blosum62,
                            3,
                            &std::sync::OnceLock::new(),
                            Some((&side, 0)),
                        )
                        .unwrap()
                        .unwrap()
                        .subject
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        for plan in &plans[1..] {
            assert!(std::sync::Arc::ptr_eq(&plans[0], plan));
        }
        assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes(&plans[0]));
        assert!(!plans[0].units.is_empty());
    }

    // No NCBI counterpart: test of when the plan cache exists (XSubjectSide::for_batches); it does
    // not change any value NCBI computes.
    // Subject plans: a search of fewer than four query batches keeps no
    // plan and every use makes its own (one batch has no later batch to read a
    // kept plan; two and three give the time it would save for its memory);
    // four or more batches keep the plans. The plan is the same either way.
    #[test]
    fn x_subject_plans_are_kept_only_from_four_batches() {
        assert_eq!(XSubjectSide::KEEP_FROM_BATCHES, 4);
        for batches in [0, 1, 2, 3] {
            assert!(XSubjectSide::for_batches(true, batches, 1, 1, false).is_none());
        }
        for batches in [4, 5, 9] {
            assert!(XSubjectSide::for_batches(true, batches, 1, 1, false).is_some());
        }
        for batches in [0, 1, 2, 3, 4, 9] {
            assert!(XSubjectSide::for_batches(false, batches, 1, 1, false).is_none());
        }
        let subject = x_test_subject(9_000);
        let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
        let scan = |side: Option<&XSubjectSide>| {
            x_scan_plan(
                &[query],
                crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
                1,
                None,
                13,
                false,
                ScoringMatrix::Blosum62,
                3,
                &std::sync::OnceLock::new(),
                side.map(|side| (side, 0)),
            )
            .unwrap()
            .unwrap()
            .subject
        };
        // Three batches: two uses make two plans, and nothing else holds them.
        let three = XSubjectSide::for_batches(true, 3, 1, 1, false);
        let (first, second) = (scan(three.as_ref()), scan(three.as_ref()));
        assert!(!std::sync::Arc::ptr_eq(&first, &second));
        assert_eq!(std::sync::Arc::strong_count(&first), 1);
        // Four batches: the second use reads the plan the first one made.
        let four = XSubjectSide::for_batches(true, 4, 1, 1, false).unwrap();
        let (kept, read) = (scan(Some(&four)), scan(Some(&four)));
        assert!(std::sync::Arc::ptr_eq(&kept, &read));
        assert_eq!(four.kept_bytes(), XSubjectSide::plan_bytes(&kept));
        assert!(first
            .frames
            .iter()
            .map(|frame| (frame.frame, &frame.aa_seq))
            .eq(kept.frames.iter().map(|frame| (frame.frame, &frame.aa_seq))));
        assert!(first
            .units
            .iter()
            .map(|unit| (unit.0, &unit.2))
            .eq(kept.units.iter().map(|unit| (unit.0, &unit.2))));
        assert!(!first.units.is_empty());
    }

    // No NCBI counterpart: test of when LOSAT_X_TBNKEEPFEW keeps the plans; it does not change any
    // value NCBI computes.
    // The plans are kept when they hold no more than `threads` plans of the longest subject: one
    // subject at every thread count, a few long ones at enough threads, many short ones never.
    #[test]
    fn x_subject_plans_are_kept_within_the_threads() {
        fn within(lengths: impl IntoIterator<Item = usize>, threads: usize) -> bool {
            XSubjectSide::kept_within_threads(lengths, threads)
        }
        assert!(within([4_700_000usize], 1));
        assert!(within([4_700_000usize], 32));
        assert!(!within([3_000usize, 3_000], 1));
        assert!(within([3_000usize, 3_000], 2));
        assert!(within([3_000usize, 2_999, 1], 2));
        assert!(!within([3_000usize, 3_000, 1], 2));
        // Three fly chromosomes (76.9 Mb): kept at 3 or more threads (three times the longest,
        // 28.1 Mb, is 84.3 Mb), not at 1 or 2.
        let fly = [23_513_712usize, 25_286_936, 28_110_227];
        assert!(!within(fly, 1));
        assert!(!within(fly, 2));
        assert!(within(fly, 3));
        assert!(within(fly, 16));
        // 28,036 16S records of about 1,450 bases: never at the thread counts of a search.
        assert!(!within(std::iter::repeat_n(1_450usize, 28_036), 32));
        assert!(!within(std::iter::repeat_n(1_450usize, 28_036), 1_024));
        // No subjects keep nothing; no threads count as one.
        assert!(within(std::iter::empty::<usize>(), 1));
        assert!(!within([3_000usize, 3_000], 0));
        assert!(within([usize::MAX, 1], 2));
    }

    // No NCBI counterpart: test of the size limit of the plan cache; it does not change any value
    // NCBI computes.
    // A plan that does not fit leaves the total as it was; a failure is
    // returned and nothing is kept for the subject.
    #[test]
    fn x_subject_plan_limit_counts_kept_plans_only() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let kept = AtomicUsize::new(0);
        assert!(XSubjectSide::reserve(&kept, 60, 100));
        assert!(!XSubjectSide::reserve(&kept, 60, 100));
        assert_eq!(kept.load(Ordering::Relaxed), 60);
        assert!(XSubjectSide::reserve(&kept, 40, 100));
        assert!(!XSubjectSide::reserve(&kept, 1, 100));
        assert!(XSubjectSide::reserve(&kept, 0, 100));
        assert_eq!(kept.load(Ordering::Relaxed), 100);
        assert!(!XSubjectSide::reserve(
            &AtomicUsize::new(usize::MAX),
            1,
            usize::MAX
        ));

        assert!(XSubjectSide::all_fit([1_000_000usize; 100]));
        assert!(XSubjectSide::all_fit([(XSubjectSide::KEEP_LIMIT - 12) / 2]));
        assert!(!XSubjectSide::all_fit([(XSubjectSide::KEEP_LIMIT - 12)
            / 2
            + 1]));
        assert!(!XSubjectSide::all_fit([usize::MAX, 1]));

        let side = XSubjectSide::new(1, 1, false);
        let slot = side.slot(0, 1, false).unwrap();
        let failing = || -> Result<std::sync::Arc<XSubjectPlan>> { bail!("no plan") };
        assert_eq!(
            side.plan(slot, failing).err().unwrap().to_string(),
            "no plan"
        );
        assert_eq!(side.kept_bytes(), 0);
        // The subject is not kept; the next use makes (and here fails) its own.
        assert_eq!(slot.get().map(Option::is_none), Some(true));
        assert!(side.plan(slot, failing).is_err());
    }
}
