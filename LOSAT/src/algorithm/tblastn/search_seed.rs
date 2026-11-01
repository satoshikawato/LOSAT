//! Internal TBLASTN protein lookup and translated-subject seed stage: the scan plan and the
//! per-unit scan of the preliminary search.

use crate::algorithm::blastp::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_query_frame_with_seg_timed,
};
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{XCompressedLookup, X_TBLASTN};
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
pub(crate) struct TranslatedChunk {
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
pub(crate) fn translated_chunk_scan_ranges(
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
    pub lookup: XPlanLookup,
    pub subject: std::sync::Arc<XSubjectPlan>,
    pub pair_capacity: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-108
// ```c
//        BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//                              &lookup_wrap->lut);
// ...
//    case eCompressedAaLookupTable:
//       BlastCompressedAaLookupTableNew(query, lookup_segments,
//                          (BlastCompressedAaLookupTable* *) &(lookup_wrap->lut),
//                          lookup_options, sbp);
// ```
/// The lookup table of a scan plan: the standard table (`eAaLookupTable`), or the compressed
/// table of word sizes 5 to 7 (`eCompressedAaLookupTable`, tblastn-fast).
pub(super) enum XPlanLookup {
    Aa(std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>),
    Compressed(std::sync::Arc<XCompressedLookup>),
}

/// The compressed lookup table of a query set (word sizes 5 to 7), kept for all its subjects as
/// `XQuerySide` keeps the standard table: the first subject that needs it builds it while the
/// others wait. A refusal is returned to that subject and leaves nothing behind.
#[derive(Default)]
pub(super) struct XCompressedSide {
    table: std::sync::OnceLock<std::sync::Arc<XCompressedLookup>>,
    /// Held while `table` is being built (on the builder's thread alone: no pool task waits on
    /// it while it runs others).
    in_progress: std::sync::Mutex<()>,
}

// No NCBI counterpart: builds the table once for the subjects of a query batch, as NCBI builds it
// once before them (blast_aux_priv.cpp:206-213); it does not change any value NCBI computes.
impl XCompressedSide {
    fn table(
        &self,
        build: impl FnOnce() -> Result<XCompressedLookup>,
    ) -> Result<std::sync::Arc<XCompressedLookup>> {
        if let Some(table) = self.table.get() {
            return Ok(std::sync::Arc::clone(table));
        }
        let _in_progress = self
            .in_progress
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(table) = self.table.get() {
            return Ok(std::sync::Arc::clone(table));
        }
        let built = std::sync::Arc::new(build()?);
        Ok(std::sync::Arc::clone(self.table.get_or_init(|| built)))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:104-108
// ```c
//    case eCompressedAaLookupTable:
//       BlastCompressedAaLookupTableNew(query, lookup_segments,
//                          (BlastCompressedAaLookupTable* *) &(lookup_wrap->lut),
//                          lookup_options, sbp);
// ```
/// The compressed lookup table of the query set (word sizes 5 to 7, `threshold` in its units,
/// `stage_d_pipeline::lookup_threshold`): the masked query frames of the standard table
/// (`x_scan_plan`), their concatenation and lookup segments as the standard table takes them (the
/// coordinates of the word finder's query, `x_lookup_query_for_profile`), and
/// blastp-fast's builders of NCBI's one table of the query batch (`x_tbn_build_compressed_lookup`:
/// the flat table, or the reference one with `LOSAT_X_BPCLUT=0`) on this thread.
fn x_build_compressed(
    queries: &[&[u8]],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    threshold: i32,
    word_length: usize,
) -> Result<XCompressedLookup> {
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_tbn_build_compressed_lookup;
    let query_frames: Vec<_> = queries
        .iter()
        .map(|query| vec![encode_tblastn_lookup_query(query, seg, mask_lowercase)])
        .collect();
    let (query, locations) =
        crate::algorithm::tblastx::lookup::x_lookup_query_for_profile(&query_frames, word_length);
    x_tbn_build_compressed_lookup(
        word_length,
        threshold,
        &query,
        &locations,
        &[0..locations.len()],
    )
    .ok_or_else(|| {
        anyhow::anyhow!(
            "-word_size {word_length} (the compressed lookup table) is not supported by LOSAT's TBLASTN"
        )
    })?
    .map_err(|()| {
        anyhow::anyhow!(
            "-threshold {} with the compressed lookup table of -word_size {} is not supported by LOSAT's TBLASTN: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
            f64::from(threshold) / 100.0,
            word_length
        )
    })
}

/// The compressed table of a scan plan (`x_scan_plan`, word sizes 5 to 7) and its longest chain,
/// built once for the query set (`XCompressedSide`). Cold and not inlined, so that the standard
/// table's plan keeps its code.
#[cold]
#[inline(never)]
fn x_compressed_plan_lookup(
    x_compressed: &XCompressedSide,
    queries: &[&[u8]],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    threshold: i32,
    word_length: usize,
) -> Result<(XPlanLookup, i32)> {
    let table = x_compressed
        .table(|| x_build_compressed(queries, seg, mask_lowercase, threshold, word_length))?;
    let longest_chain = match &*table {
        XCompressedLookup::Reference(table) => table.longest_chain,
        XCompressedLookup::Flat(table) => table.longest_chain,
    };
    Ok((XPlanLookup::Compressed(table), longest_chain))
}

/// The word length of a compressed table (`BlastCompressedAaLookupTable::word_length`).
fn x_compressed_word_length(table: &XCompressedLookup) -> i32 {
    match table {
        XCompressedLookup::Reference(table) => table.word_length,
        XCompressedLookup::Flat(table) => table.word_length,
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
// next plan of the thread instead.
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
/// kept only up to a total size (`keeps`: decided from the subject's length by
/// the first batch that reaches it), beyond which they are made per use as the
/// plans of a search without a side are: without residues when the subject is
/// long (LOSAT_X_TBNUNITXLATE, `x_scan_plan`), so that four batches over a
/// genome hold the frames of the chromosomes they are on no more than one does.
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
/// Batches can reach a subject at the same time. Whether its plan is kept is
/// decided by one of them while the others wait, and its bytes are added to
/// `kept_bytes` once; the kept plan is then made by one of them while the
/// others wait, so a plan is made once: `kept_bytes` is always the size of the
/// plans the cells hold.
pub(super) struct XSubjectSide {
    gencode: u8,
    mask_lowercase: bool,
    /// One cell per subject.
    slots: Vec<XSubjectSlot>,
    kept_bytes: std::sync::atomic::AtomicUsize,
    /// `KEEP_LIMIT` (a smaller limit in tests).
    limit: usize,
}

/// The cell of one subject of `XSubjectSide`: unset until a batch reaches the
/// subject.
#[derive(Default)]
struct XSubjectSlot {
    /// Whether the plan is kept (its bytes reserved within the limit), decided
    /// by the first batch that reaches the subject (`XSubjectSide::keeps`).
    kept: std::sync::OnceLock<bool>,
    /// The kept plan, made by the first batch that reaches the subject, or
    /// `None` when it could not be prepared (every later use then makes its
    /// own).
    plan: std::sync::OnceLock<Option<std::sync::Arc<XSubjectPlan>>>,
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
    /// unless it is `0`. With it a search of
    /// `KEEP_FROM_BATCHES` or more batches keeps its plans only when `kept_within_threads`
    /// holds or it searches batches side by side (`stage_d_pipeline::x_searched_side_by_side`).
    pub(super) fn keep_few() -> bool {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNKEEPFEW"))
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
    /// plan; with the frame reuse 0.998 at one thread).
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
        Self::with_limit(subjects, gencode, mask_lowercase, Self::KEEP_LIMIT)
    }

    /// `new`, with `limit` bytes of kept plans instead of `KEEP_LIMIT`.
    fn with_limit(subjects: usize, gencode: u8, mask_lowercase: bool, limit: usize) -> Self {
        Self {
            gencode,
            mask_lowercase,
            slots: (0..subjects).map(|_| XSubjectSlot::default()).collect(),
            kept_bytes: std::sync::atomic::AtomicUsize::new(0),
            limit,
        }
    }

    /// The subject plans of a search of `batches` query batches, or `None`
    /// when the search keeps none: `on` is false (`LOSAT_X_TBNKEEPFEW` keeps
    /// none, `kept_within_threads`), or there are
    /// fewer than `KEEP_FROM_BATCHES` batches. A kept plan is read by
    /// the batches after the first one that reaches its subject. A search of
    /// fewer batches makes its plans per use and drops them after the subject,
    /// like those `keeps` refuses.
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

    /// The cell of subject `oid`, when the side is one of this genetic code
    /// and masking option.
    fn slot(&self, oid: usize, gencode: u8, mask_lowercase: bool) -> Option<&XSubjectSlot> {
        if self.gencode != gencode || self.mask_lowercase != mask_lowercase {
            return None;
        }
        self.slots.get(oid)
    }

    /// Translated residues of a plan: what `KEEP_LIMIT` counts
    /// (`plan_bytes_of`, from the length; this one checks it in tests).
    #[cfg(test)]
    fn plan_bytes(plan: &XSubjectPlan) -> usize {
        plan.frames.iter().map(|frame| frame.aa_seq.len()).sum()
    }

    /// `plan_bytes` of the plan of a subject of `length` bases, from the
    /// length alone: the frames `i + 1` and `-(i + 1)` (`i + 3 <= length`)
    /// hold `(length - i) / 3` residues and two sentinels each
    /// (`x_subject_translate::x_packed_frames`, `x_frame`).
    fn plan_bytes_of(length: usize) -> usize {
        (0..3)
            .filter(|&i| i + 3 <= length)
            .map(|i| ((length - i) / 3 + 2).saturating_mul(2))
            .fold(0usize, usize::saturating_add)
    }

    /// Whether the plan of the subject in `slot`, of `length` bases, is kept:
    /// decided by the first batch that reaches the subject, which reserves
    /// the plan's bytes within the limit (`plan_bytes_of`). A plan that does
    /// not fit is made per use, as the plans of a search without a side are
    /// (`x_scan_plan`; `x_lean_subject` reads the decision before a `-db`
    /// sequence is read for its units).
    fn keeps(&self, slot: &XSubjectSlot, length: usize) -> bool {
        *slot.kept.get_or_init(|| {
            let kept = Self::reserve(&self.kept_bytes, Self::plan_bytes_of(length), self.limit);
            // The own side is recorded where a plan is kept or refused.
            crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", kept);
            kept
        })
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

    /// Adds `bytes` to `kept` if the total stays within `limit`. The total
    /// only changes when the plan fits, so a plan that does not fit never
    /// makes another one look too large.
    fn reserve(kept: &std::sync::atomic::AtomicUsize, bytes: usize, limit: usize) -> bool {
        use std::sync::atomic::Ordering;
        kept.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |held| {
            held.checked_add(bytes).filter(|&total| total <= limit)
        })
        .is_ok()
    }

    /// The kept plan of the subject in `slot` (`keeps`, of `length` bases):
    /// made by the first caller and read by the others. When it could not be
    /// made, the first caller gets the error, its bytes are given back, and
    /// every later caller makes its own.
    fn plan(
        &self,
        slot: &XSubjectSlot,
        length: usize,
        build: impl Fn() -> Result<std::sync::Arc<XSubjectPlan>>,
    ) -> Result<std::sync::Arc<XSubjectPlan>> {
        let mut failed = None;
        let kept = slot.plan.get_or_init(|| match build() {
            Ok(plan) => Some(plan),
            Err(error) => {
                self.kept_bytes.fetch_sub(
                    Self::plan_bytes_of(length),
                    std::sync::atomic::Ordering::Relaxed,
                );
                failed = Some(error);
                None
            }
        });
        if let Some(error) = failed {
            return Err(error);
        }
        match kept {
            Some(plan) => Ok(std::sync::Arc::clone(plan)),
            None => build(),
        }
    }
}

// No NCBI counterpart: the rule of LOSAT_X_TBNUNITXLATE in `x_scan_plan` (its packed arm) for a
// subject that is not read yet; it computes no value NCBI computes.
/// Whether `x_scan_plan` makes the plan of a packed subject of
/// `length` bases without residues (LOSAT_X_TBNUNITXLATE): `x_units_translate` and the plan is
/// made for its use alone (`x_own_plan`). The length is read first, so only a subject long enough
/// for a plan without residues decides here whether its plan is kept; the others decide in
/// `x_scan_plan`, after its check of the subject's length (a subject too short to be searched
/// charges nothing).
pub(super) fn x_lean_subject(
    x_subject: Option<(&XSubjectSide, usize)>,
    length: usize,
    db_gencode: u8,
    mask_lowercase: bool,
) -> bool {
    super::x_subject_translate::x_units_translate(length)
        && x_own_plan(x_subject, length, db_gencode, mask_lowercase)
}

// No NCBI counterpart: where a plan is made for one use (NCBI makes every plan so); it computes
// no value NCBI computes.
/// Whether the plan of subject `oid` of `x_subject`, of `length` bases, is made for its use
/// alone: the search keeps no plans (no side, or one of another genetic code or masking
/// option), or none of this subject (`XSubjectSide::keeps`, decided once per subject).
fn x_own_plan(
    x_subject: Option<(&XSubjectSide, usize)>,
    length: usize,
    db_gencode: u8,
    mask_lowercase: bool,
) -> bool {
    x_subject.is_none_or(|(side, oid)| {
        side.slot(oid, db_gencode, mask_lowercase)
            .is_none_or(|slot| !side.keeps(slot, length))
    })
}

/// No NCBI counterpart: the offset array of `LOSAT_X_TBNPAIRCAP` holds this many entries more
/// than the longest chain, eight times NCBI's OFFSET_ARRAY_SIZE (lookup_wrap.h:119). A scan call
/// that stops on a full array starts again in the next call at the word that did not fit
/// (`x_scan_byte.rs` builds that word's block of candidates again), which costs little per call
/// but adds up at NCBI's size: dmel_b4.faa x zebrafish cDNAs (`-db`) at one thread took 1.0041 of
/// the instructions with 4096 + the longest chain, 1.0016 with 8192 +, 1.0004 with 16384 + and
/// 1.0001 with 32768 + (against 1024 x the longest chain; D14 s1 c1, 2026-10-26). With
/// `LOSAT_X_TBNSCANRESUME` a stopped call keeps its block, and the units that take its scan have
/// NCBI's size (`OFFSET_ARRAY_SIZE`).
const X_PAIR_ROOM: usize = 8 * 4096;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.h:119
// #define OFFSET_ARRAY_SIZE 4096
/// NCBI's OFFSET_ARRAY_SIZE: the offset array of a unit scanned with `LOSAT_X_TBNSCANRESUME` holds
/// this many entries more than the longest chain (GetOffsetArraySize, lookup_wrap.c:264-266).
const OFFSET_ARRAY_SIZE: usize = 4096;

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// Whether a plan's offset array holds `X_PAIR_ROOM` + the longest chain entries (at most
/// 1024 x the longest chain, the size without the switch) (`LOSAT_X_TBNPAIRCAP`, read once).
fn x_pair_cap_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNPAIRCAP"))
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
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subject's letters, or its packed store; or a `-db`
    // sequence, read here when its plan is made.
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
    // Where the compressed table of this query set is kept (word sizes 5 to 7, `XQuerySide`).
    x_compressed: &XCompressedSide,
    // The search's subject plans and this subject's index among them.
    x_subject: Option<(&XSubjectSide, usize)>,
    // The query set's frames as the table's query (`encode_tblastn_lookup_query` of each query
    // with `seg` and `mask_lowercase`), made once for the batch (`XQuerySide`), when the caller
    // has them.
    x_frames: Option<&[QueryFrame]>,
) -> Result<Option<XScanPlan>> {
    // An unknown genetic code is refused before the subject length is looked at.
    GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1334,1429-1432
    // return word_length * 3 + 2;
    // if (seq_arg.seq->length < min_subj_seq_length) { ... continue; }
    let length = subject.len();
    if length < word_length * 3 + 2 {
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
    // The slot of the subject's kept plan when the search keeps them (the dispatch point below);
    // a plan made for this use only may have no residues (LOSAT_X_TBNUNITXLATE).
    let x_slot = x_subject.and_then(|(side, oid)| {
        side.slot(oid, db_gencode, mask_lowercase)
            .map(|slot| (side, slot))
    });
    // Whether the plan is made for this use alone: the search keeps no plans, or none of this
    // subject (`XSubjectSide::keeps`, decided once per subject; `x_lean_subject` read the same
    // decision before a `-db` sequence was read for its units).
    let x_own = x_slot.is_none_or(|(side, slot)| !side.keeps(slot, length));
    let x_build_subject = || -> Result<std::sync::Arc<XSubjectPlan>> {
        use super::x_subject_store::XPrelimSubject;
        use crate::sequence::packed_subject::SubjectNt;
        match subject {
            // The plan is made from the sequence's bytes, read
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
                // EXPERIMENT (LOSAT_X_TBNUNITXLATE): the plan of a long subject made for this
                // use only (no side, or a plan the side does not keep) has no residues; each
                // unit translates those it reads.
                if x_own {
                    let lean = super::x_subject_translate::x_units_translate(store.len());
                    crate::utils::xreach::x_reach!("LOSAT_X_TBNUNITXLATE", lean);
                    if lean {
                        return super::x_subject_translate::x_packed_lean_plan(store, db_gencode);
                    }
                }
                super::x_subject_translate::x_packed_subject_plan(store, db_gencode)
            }
        }
    };
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
    // first batch that reaches the subject); a search of fewer than four batches makes it here,
    // and so does every use of a plan the side does not keep (`XSubjectSide::keeps`).
    let subject_plan = match x_slot {
        Some((side, slot)) if !x_own => side.plan(slot, length, x_build_subject)?,
        _ => x_build_subject()?,
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
        let query_frames: Vec<_> = match x_frames {
            Some(frames) => frames.iter().map(|frame| vec![frame.clone()]).collect(),
            None => queries
                .iter()
                // NCBI blast_filter.c:1241-1255 merges lowercase query and SEG masks.
                .map(|query| vec![encode_tblastn_lookup_query(query, seg, mask_lowercase)])
                .collect(),
        };
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
    // ```c
    // CBlastOptionsLocal::SetWordSize(int ws)
    // {
    //     m_LutOpts->word_size = ws;
    //     if (m_LutOpts->lut_type == eCompressedAaLookupTable && ws <= 4)
    // 	m_LutOpts->lut_type = eAaLookupTable;
    //     else if (m_LutOpts->lut_type == eAaLookupTable && ws > 4)
    // 	m_LutOpts->lut_type = eCompressedAaLookupTable;
    // }
    // ```
    // Dispatch point: word sizes 5 to 7 take the compressed table (`x_build_compressed`).
    let (lookup, longest_chain) = if word_length > 4 {
        x_compressed_plan_lookup(
            x_compressed,
            queries,
            seg,
            mask_lowercase,
            threshold,
            word_length,
        )?
    } else {
        let lookup = std::sync::Arc::clone(x_lookup.get_or_init(x_build_lookup));
        let longest_chain = lookup.longest_chain;
        (XPlanLookup::Aa(lookup), longest_chain)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-270
    // ```c
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eCompressedAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastCompressedAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    // blast_engine.c:992,1040-1041 allocates each thread's offset array of that size and
    // aa_ungapped.c:504-505 passes it to scansub (OFFSET_ARRAY_SIZE 4096, lookup_wrap.h:119).
    // The size only decides where a scan call stops; the next call resumes there, so the pairs
    // reach the caller in the same order at any size. With `LOSAT_X_TBNPAIRCAP` the array has
    // NCBI's form with `X_PAIR_ROOM` for OFFSET_ARRAY_SIZE, never more entries than without it
    // (1024 x the longest chain, which a table whose longest chain is below 33 keeps).
    let longest_chain = longest_chain.max(1) as usize;
    let pair_capacity = if x_pair_cap_on() {
        (X_PAIR_ROOM + longest_chain).min(longest_chain.saturating_mul(1024))
    } else {
        longest_chain * 1024
    };
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
/// `subject` holds the unit's frame from the chunk's first residue on (the frame's `aa_seq`
/// from `1 + chunk.offset`, or the residues the unit translated, `x_unit_residues`).
pub(super) fn x_scan_unit(
    plan: &XScanPlan,
    unit: usize,
    subject: &[u8],
    pairs: &mut Vec<BlastOffsetPair>,
    mut on_pairs: impl FnMut(&[BlastOffsetPair]),
) {
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1018-1019
    // BlastChooseProteinScanSubject(lookup_wrap);
    // aux_struct->WordFinder = BlastAaWordFinder;
    // The scan variant of the unit is chosen once per table, as NCBI picks one scansub: the
    // compressed table's (`x_scan_compressed_call`, range by range as below), or one of the
    // standard table's below.
    let lookup = match &plan.lookup {
        XPlanLookup::Aa(lookup) => lookup,
        XPlanLookup::Compressed(table) => {
            let (_, _, ranges) = &plan.subject.units[unit];
            if pairs.len() < plan.pair_capacity {
                pairs.resize(plan.pair_capacity, BlastOffsetPair::default());
            }
            let pair_capacity = plan.pair_capacity as i32;
            let word_length = x_compressed_word_length(table);
            for (range_index, &(left, right)) in ranges.iter().enumerate() {
                let mut end = right - word_length;
                if range_index == 0 && end < left {
                    end = left;
                }
                let mut scan_range = [0, left, end];
                while scan_range[1] <= scan_range[2] {
                    let hits = x_scan_compressed_call(
                        table,
                        subject,
                        (left, right),
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                    );
                    on_pairs(&pairs[..hits as usize]);
                }
            }
            return;
        }
    };
    let (_, _, ranges) = &plan.subject.units[unit];
    //
    // EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan on presence bytes, for a table of tblastn's shape
    // (without AVX2 a sparse one); it writes up to X_SCAN_BYTE_SLACK entries after the capacity.
    let scan_byte = super::x_scan_byte::x_scan_byte();
    let presence = super::x_scan_byte::x_presence(scan_byte, lookup);
    // EXPERIMENT (LOSAT_X_TBNSCANRESUME): a call of the scan on presence bytes that stops on a full
    // offset array leaves its block in `held`, and the next call copies on from there instead of
    // making the block's candidates again; the array then has NCBI's size (OFFSET_ARRAY_SIZE + the
    // longest chain, lookup_wrap.c:264-266), never more entries than the plan's.
    let resume = presence.is_some() && super::x_scan_byte::x_scan_resume();
    let capacity = if resume {
        let longest_chain = lookup.longest_chain.max(1) as usize;
        (OFFSET_ARRAY_SIZE + longest_chain).min(plan.pair_capacity)
    } else {
        plan.pair_capacity
    };
    // The branch-free scan writes up to X_SCAN_BF_SLACK entries after the capacity.
    let pairs_len = if presence.is_some() {
        capacity + super::x_scan_byte::X_SCAN_BYTE_SLACK
    } else {
        capacity + super::x_scan_bf::X_SCAN_BF_SLACK
    };
    if pairs.len() < pairs_len {
        pairs.resize(pairs_len, BlastOffsetPair::default());
    }
    let pair_capacity = capacity as i32;
    let mut held = super::x_scan_byte::XScanHeld::new();
    // The scan calls of the unit, for LOSAT_X_STATS.
    let mut x_calls = 0u64;
    for (range_index, &(left, right)) in ranges.iter().enumerate() {
        // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:496-501
        // scan_range[1] = subject->seq_ranges[0].left;
        // scan_range[2] = subject->seq_ranges[0].right - wordsize;
        // if (scan_range[2] < scan_range[1])
        //     scan_range[2] = scan_range[1];
        let mut end = right - lookup.word_length as i32;
        if range_index == 0 && end < left {
            end = left;
        }
        let mut scan_range = [0, left, end];
        while scan_range[1] <= scan_range[2] {
            let hits = if let Some(presence) = &presence {
                if resume {
                    super::x_scan_byte::x_scan_one_range_held(
                        scan_byte,
                        lookup,
                        presence,
                        subject,
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                        &mut held,
                    )
                } else {
                    super::x_scan_byte::x_scan_one_range_byte(
                        scan_byte,
                        lookup,
                        presence,
                        subject,
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                    )
                }
            } else {
                super::x_scan_bf::x_scan_one_range_bf(
                    lookup,
                    subject,
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
        super::x_scan_bf::x_count_unit(x_calls);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:487-500 and
// c++/src/algo/blast/core/blast_aascan.c:564-568
// ```c
//     else {
//         BlastCompressedAaLookupTable *lookup =
//                         (BlastCompressedAaLookupTable *)(lookup_wrap->lut);
//         scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
//         wordsize = lookup->word_length;
//     }
//
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
// ...
//         lut->scansub_callback = (void *)s_BlastCompressedAaScanSubject;
// ```
/// One scansub call of `s_BlastCompressedAaScanSubject` (the table's `scan_subject`, which writes
/// no entry past `array_size`) on the range `(left, right)` from `scan_range`, for the compressed
/// arm of `x_scan_unit` (`s_DetermineScanningOffsets` of a list of one range ends after it, and a
/// call that stops at a full array resumes in it). Not inlined, and it takes no closure, so that
/// the standard table's scan keeps its code.
#[inline(never)]
fn x_scan_compressed_call(
    table: &XCompressedLookup,
    subject: &[u8],
    (left, right): (i32, i32),
    pairs: &mut [BlastOffsetPair],
    pair_capacity: i32,
    scan_range: &mut [i32; 3],
) -> i32 {
    let range = [(left, right)];
    match table {
        XCompressedLookup::Reference(table) => {
            table.scan_subject(subject, &range, pairs, pair_capacity, scan_range)
        }
        XCompressedLookup::Flat(table) => {
            table.scan_subject::<X_TBLASTN>(subject, &range, pairs, pair_capacity, scan_range)
        }
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
                    &Default::default(),
                    None,
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
            &Default::default(),
            None,
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
                &Default::default(),
                None,
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
    // plan. Its bytes are charged to the shared limit once (`keeps`, by the
    // first of them) and it is made once.
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
                        assert!(side.keeps(slot, subject.len()));
                        side.plan(slot, subject.len(), &build).unwrap()
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
        assert_eq!(XSubjectSide::plan_bytes_of(subject.len()), bytes);
        assert_eq!(side.kept_bytes(), bytes);
        // A later batch reads the same plan and adds nothing.
        assert!(side.keeps(slot, subject.len()));
        let again = side.plan(slot, subject.len(), &build).unwrap();
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
                            &Default::default(),
                            Some((&side, 0)),
                            None,
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
                &Default::default(),
                side.map(|side| (side, 0)),
                None,
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
        // The kept plan is the plan a use makes for itself, in every field. For subject plans
        // `aa_seq_nomask` is None and `seg_masks` is empty (and `mask_lowercase` is off here), so
        // the comparison is trivial for those two; it covers the residues, lengths and units.
        assert_eq!(
            crate::algorithm::tblastn::x_subject_translate::x_plan_fields_difference(&first, &kept),
            None
        );
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
        assert!(side.keeps(slot, 30_000));
        assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes_of(30_000));
        let failing = || -> Result<std::sync::Arc<XSubjectPlan>> { bail!("no plan") };
        assert_eq!(
            side.plan(slot, 30_000, failing).err().unwrap().to_string(),
            "no plan"
        );
        // The bytes of the plan that failed are given back; the next use makes (and here
        // fails) its own.
        assert_eq!(side.kept_bytes(), 0);
        assert_eq!(slot.plan.get().map(Option::is_none), Some(true));
        assert!(side.plan(slot, 30_000, failing).is_err());
    }

    // No NCBI counterpart: test of the size a subject's plan is charged before it is made; it does
    // not change any value NCBI computes.
    // The bytes of a plan from the subject's length alone are the bytes of the plan made from its
    // letters or its packed bases, at every length (no frame below three bases; two sentinels per
    // frame), and never above the bound `all_fit` uses.
    #[test]
    fn x_plan_bytes_of_a_length_are_the_bytes_of_its_plan() {
        use crate::algorithm::tblastn::x_subject_translate::{
            x_packed_subject_plan, x_subject_plan,
        };
        use crate::sequence::packed_subject::PackedNucleotide;
        for length in (1..=40).chain([1_000, 4_999, 100_003]) {
            let subject = x_test_subject(length);
            let bytes = XSubjectSide::plan_bytes_of(length);
            let letters = x_subject_plan(&subject, 1, false).unwrap();
            assert_eq!(
                XSubjectSide::plan_bytes(&letters),
                bytes,
                "{length} bases, letters"
            );
            let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
            let packed = x_packed_subject_plan(&store, 1).unwrap();
            assert_eq!(
                XSubjectSide::plan_bytes(&packed),
                bytes,
                "{length} bases, packed"
            );
            assert!(bytes <= 2 * length + 12, "{length} bases");
        }
        assert_eq!(XSubjectSide::plan_bytes_of(0), 0);
        assert_eq!(XSubjectSide::plan_bytes_of(2), 0);
        assert_eq!(XSubjectSide::plan_bytes_of(3), 6);
        assert_eq!(XSubjectSide::plan_bytes_of(5), 6 + 6 + 6);
        assert_eq!(XSubjectSide::plan_bytes_of(usize::MAX), usize::MAX);
        // A subject too short to be searched (x_scan_plan returns None before the plan) is not
        // charged: the pre-read decides only for a subject long enough for a lean plan, and
        // x_scan_plan decides after its length check.
        let side = XSubjectSide::new(1, 1, false);
        assert!(!x_lean_subject(Some((&side, 0)), 10, 1, false));
        assert_eq!(side.kept_bytes(), 0);
        assert!(x_scan_plan(
            &[b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ".as_slice()],
            crate::sequence::packed_subject::SubjectNt::Letters(b"ACGTACGTAC").into(),
            1,
            None,
            13,
            false,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            &Default::default(),
            Some((&side, 0)),
            None,
        )
        .unwrap()
        .is_none());
        assert_eq!(side.kept_bytes(), 0);
    }

    // No NCBI counterpart: test of the plans a side does not keep; it does not change any value
    // NCBI computes.
    // A plan the side's limit refuses is made per use, as without a side: without residues when
    // the subject is long (LOSAT_X_TBNUNITXLATE), with its frames when it is short; the decision
    // is made once per subject, before the plan, and `x_lean_subject` reads it. A plan that fits
    // is kept, with its residues, and read by the next use.
    #[test]
    fn x_plans_the_side_does_not_keep_are_made_per_use_as_without_a_side() {
        use crate::algorithm::tblastn::x_subject_translate::{
            x_plan_fields_difference, x_units_translate,
        };
        use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};
        let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
        let scan = |side: Option<&XSubjectSide>, store: &PackedNucleotide| {
            x_scan_plan(
                &[query],
                SubjectNt::Packed(store).into(),
                1,
                None,
                13,
                false,
                ScoringMatrix::Blosum62,
                3,
                &std::sync::OnceLock::new(),
                &Default::default(),
                side.map(|side| (side, 0)),
                None,
            )
            .unwrap()
            .unwrap()
            .subject
        };
        let long =
            PackedNucleotide::from_subject_letters(&x_test_subject(3_200_000)).expect("letters");
        let short =
            PackedNucleotide::from_subject_letters(&x_test_subject(9_000)).expect("letters");
        for (store, lean) in [(&long, x_units_translate(3_200_000)), (&short, false)] {
            let length = store.len();
            // A limit of no bytes refuses every plan: each use makes its own, the plan a use
            // makes without a side in every field.
            let side = XSubjectSide::with_limit(1, 1, false, 0);
            assert_eq!(x_lean_subject(Some((&side, 0)), length, 1, false), lean);
            let (first, second) = (scan(Some(&side), store), scan(Some(&side), store));
            assert!(!std::sync::Arc::ptr_eq(&first, &second));
            assert_eq!(side.kept_bytes(), 0);
            assert!(!side.keeps(side.slot(0, 1, false).unwrap(), length));
            assert_eq!(first.frames[0].aa_seq.is_empty(), lean, "{length} bases");
            assert_eq!(x_plan_fields_difference(&first, &scan(None, store)), None);
            // The default limit keeps both: the plan has its residues and the next use reads it.
            let side = XSubjectSide::new(1, 1, false);
            assert!(!x_lean_subject(Some((&side, 0)), length, 1, false));
            let (kept, read) = (scan(Some(&side), store), scan(Some(&side), store));
            assert!(std::sync::Arc::ptr_eq(&kept, &read));
            assert!(!kept.frames[0].aa_seq.is_empty());
            assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes(&kept));
            assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes_of(length));
        }
    }

    // A small deterministic generator (xorshift64) for the compressed scan tests.
    struct XRng(u64);

    impl XRng {
        fn below(&mut self, n: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 % n as u64) as usize
        }
    }

    /// The hits of cell `index` of a compressed table, in NCBI's scan order.
    fn x_cell_hits(table: &XCompressedLookup, index: usize) -> Vec<i32> {
        match table {
            XCompressedLookup::Reference(table) => {
                if table.backbone[index].num_used == 0 {
                    Vec::new()
                } else {
                    table.query_offsets_in_scan_order(index)
                }
            }
            XCompressedLookup::Flat(table) => table.cell(index).to_vec(),
        }
    }

    /// The pairs of every unit of `plan`: what `x_scan_unit` hands to its caller, and what NCBI's
    /// compressed scan finds (blast_aascan.c:264-409): at each offset of each range whose word has
    /// a compressed letter at every position, the hits of the word's cell (index sum c_k A^k, the
    /// index the scan's preshift forms), in order.
    fn x_compressed_unit_pairs(plan: &XScanPlan) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
        let XPlanLookup::Compressed(table) = &plan.lookup else {
            panic!("a compressed plan")
        };
        let (word, alphabet, scaled) = match &**table {
            XCompressedLookup::Reference(t) => (
                t.word_length,
                t.compressed_alphabet_size,
                t.scaled_compress_table_u8,
            ),
            XCompressedLookup::Flat(t) => (
                t.word_length,
                t.compressed_alphabet_size,
                t.scaled_compress_table_u8,
            ),
        };
        let top = alphabet.pow(word as u32 - 1);
        let mut scanned = Vec::new();
        let mut expected = Vec::new();
        let mut pairs = Vec::new();
        for (unit, (frame_index, chunk, ranges)) in plan.subject.units.iter().enumerate() {
            let residues = &plan.subject.frames[*frame_index].aa_seq[1 + chunk.offset..];
            x_scan_unit(plan, unit, residues, &mut pairs, |filled| {
                scanned.extend(filled.iter().map(|pair| (pair.q_off, pair.s_off)));
            });
            for (range_index, &(left, right)) in ranges.iter().enumerate() {
                let mut end = right - word;
                if range_index == 0 && end < left {
                    end = left;
                }
                for s in left..=end.max(left - 1) {
                    let letters: Vec<i32> = (0..word)
                        .map(|k| {
                            let letter = residues.get((s + k) as usize).copied().unwrap_or(0);
                            scaled[usize::from(letter)]
                        })
                        .collect();
                    if letters.iter().any(|&c| c < 0) {
                        continue;
                    }
                    let index: i32 = letters
                        .iter()
                        .enumerate()
                        .map(|(k, &c)| c / top * alphabet.pow(k as u32))
                        .sum();
                    for hit in x_cell_hits(table, index as usize) {
                        expected.push((hit as u32, s as u32));
                    }
                }
            }
        }
        (scanned, expected)
    }

    fn x_random_protein(rng: &mut XRng, length: usize) -> Vec<u8> {
        const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYXBZU*";
        (0..length)
            .map(|_| LETTERS[rng.below(LETTERS.len())])
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:264-409 and
    // c++/src/algo/blast/core/aa_ungapped.c:494-516
    // ```c
    //     while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
    // ...
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                                   offset_pairs, array_size, scan_range);
    // ```
    // The unit scan of the compressed table (word sizes 5 to 7) finds the hits of NCBI's scan, in
    // its order, with any offset array (one that holds only the longest cell resumes at every
    // cell), on the flat table and, for word size 5, on the reference table, which is the flat
    // one; a threshold of 100 puts only the exact words in the table.
    #[test]
    fn x_compressed_unit_scan_finds_ncbis_hits_in_order() {
        use crate::algorithm::tblastx::lookup::compressed::try_build_blosum62_compressed_lookup;
        use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
            x_assert_same_compressed_table, x_try_build_flat_compressed_lookup, XFlatBuild,
        };
        let mut rng = XRng(0x9e37_79b9_7f4a_7c15);
        let queries: Vec<Vec<u8>> = (0..3)
            .map(|i| x_random_protein(&mut rng, 40 + 23 * i))
            .collect();
        let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        const BASES: &[u8] = b"ACGTACGTACGTACGTACGTN";
        let mut subject: Vec<u8> = (0..2400).map(|_| BASES[rng.below(BASES.len())]).collect();
        // Fragments of the queries, back-translated (one codon for each letter, `NNN` for the
        // others) into the subject, forward and reverse-complemented, so that the exact words of
        // a threshold of 100 have hits too.
        let codon = |letter: u8| -> &'static [u8] {
            match letter {
                b'A' => b"GCT",
                b'C' => b"TGT",
                b'D' => b"GAT",
                b'E' => b"GAA",
                b'F' => b"TTT",
                b'G' => b"GGT",
                b'H' => b"CAT",
                b'I' => b"ATT",
                b'K' => b"AAA",
                b'L' => b"CTG",
                b'M' => b"ATG",
                b'N' => b"AAT",
                b'P' => b"CCG",
                b'Q' => b"CAG",
                b'R' => b"CGT",
                b'S' => b"TCT",
                b'T' => b"ACT",
                b'V' => b"GTT",
                b'W' => b"TGG",
                b'Y' => b"TAT",
                _ => b"NNN",
            }
        };
        for (index, query) in queries.iter().enumerate() {
            let mut fragment: Vec<u8> = query[5..35]
                .iter()
                .flat_map(|&l| codon(l).to_vec())
                .collect();
            if index == 1 {
                fragment.reverse();
                for base in &mut fragment {
                    *base = match *base {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        other => other,
                    };
                }
            }
            let at = 300 + 700 * index;
            subject.splice(at..at, fragment);
        }
        for (word, threshold) in [
            (5usize, 1100),
            (5, 2000),
            (5, 10000),
            (6, 2100),
            (6, 10000),
            (7, 2025),
            (7, 10000),
        ] {
            let frames: Vec<_> = query_refs
                .iter()
                .map(|query| vec![encode_tblastn_lookup_query(query, None, false)])
                .collect();
            let (query, locations) =
                crate::algorithm::tblastx::lookup::x_lookup_query_for_profile(&frames, word);
            let tables = [0..locations.len()];
            let XFlatBuild::Table(flat) =
                x_try_build_flat_compressed_lookup(word, threshold, &query, &locations, &tables)
                    .unwrap()
            else {
                panic!("word {word} threshold {threshold}: no flat table")
            };
            if threshold == 10000 {
                // Only the exact words: one hit per query offset whose word has compressed letters.
                let hits: usize = (0..flat.starts.len() - 1).map(|i| flat.cell(i).len()).sum();
                let words: usize = locations
                    .iter()
                    .map(|&(from, to)| {
                        (from..=to - word as i32 + 1)
                            .filter(|&o| {
                                (0..word).all(|k| {
                                    flat.scaled_compress_table_u8
                                        [usize::from(query[(o as usize) + k])]
                                        >= 0
                                })
                            })
                            .count()
                    })
                    .sum();
                assert_eq!(hits, words, "word {word}: exact words only");
            }
            // The side holds the flat table, so that the plan takes it whatever LOSAT_X_BPCLUT says.
            let flat = std::sync::Arc::new(XCompressedLookup::Flat(flat));
            let compressed = XCompressedSide::default();
            assert!(compressed.table.set(std::sync::Arc::clone(&flat)).is_ok());
            let mut plan = x_scan_plan(
                &query_refs,
                SubjectNt::Letters(&subject).into(),
                1,
                None,
                threshold,
                false,
                ScoringMatrix::Blosum62,
                word,
                &std::sync::OnceLock::new(),
                &compressed,
                None,
                None,
            )
            .unwrap()
            .unwrap();
            let XPlanLookup::Compressed(planned) = &plan.lookup else {
                panic!("word {word}: the plan of the compressed table")
            };
            assert!(std::sync::Arc::ptr_eq(planned, &flat));
            let mut tables = vec![flat];
            if word == 5 {
                let reference = try_build_blosum62_compressed_lookup(
                    word,
                    threshold,
                    &query,
                    &locations,
                    &[0..locations.len()],
                )
                .unwrap()
                .unwrap();
                let XCompressedLookup::Flat(flat) = &*tables[0] else {
                    unreachable!()
                };
                x_assert_same_compressed_table(flat, &reference);
                tables.push(std::sync::Arc::new(XCompressedLookup::Reference(reference)));
            }
            let planned_capacity = plan.pair_capacity;
            let mut first = None;
            for table in tables {
                let longest = match &*table {
                    XCompressedLookup::Reference(t) => t.longest_chain,
                    XCompressedLookup::Flat(t) => t.longest_chain,
                };
                plan.lookup = XPlanLookup::Compressed(table);
                for capacity in [planned_capacity, longest.max(1) as usize] {
                    plan.pair_capacity = capacity;
                    let (scanned, expected) = x_compressed_unit_pairs(&plan);
                    assert_eq!(scanned, expected, "word {word} threshold {threshold}");
                    assert!(
                        !scanned.is_empty(),
                        "word {word} threshold {threshold}: no hits"
                    );
                    match &first {
                        None => first = Some(scanned),
                        Some(first) => assert_eq!(&scanned, first, "word {word}"),
                    }
                }
            }
        }
    }
}
