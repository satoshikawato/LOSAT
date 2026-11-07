//! The chunks of a translated subject frame and the scan ranges of each (`TranslatedChunk`,
//! `translated_chunk_scan_ranges`): the schedule of NCBI's `MAX_DBSEQ_LEN` subject chunks.
//!
//! Moved from `tblastn/search_seed.rs` (SW6 item 3) so that `tblastx` and `tblastn` both use it
//! without `tblastx` depending on `tblastn`.

use anyhow::{bail, Result};

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
pub(crate) fn unmasked_translated_chunks(frame_length: usize) -> Vec<TranslatedChunk> {
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
