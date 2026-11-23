//! The chunks of a translated subject frame and the scan ranges of each (`TranslatedChunk`,
//! `translated_chunk_scan_ranges`): the schedule of NCBI's `MAX_DBSEQ_LEN` subject chunks.
//!
//! Moved from `tblastn/search_seed.rs` (SW6 item 3) so that `tblastx` and `tblastn` both use it
//! without `tblastx` depending on `tblastn`.

use anyhow::{bail, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-264:
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54:
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:192:
// This schedule covers the unmasked translated frame's single hard range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct TranslatedChunk {
    pub offset: usize,
    pub length: usize,
    /// The overlap that the merge after the chunk reads (`Blast_HSPListsMerge`'s
    /// `chunk_overlap_size`): 0 where the chunk starts the hard range at `hm_index`, else
    /// `DBSEQ_CHUNK_OVERLAP` ([`TranslatedChunk::new`] for a frame without hard ranges;
    /// `core::subject_split::HardChunk::overlap` for one with `-db_hard_mask`).
    pub overlap: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:578-584
// After each chunk the merge takes no overlap when the chunk's offset is the left of the hard
// range at `hm_index`, else the chunk overlap. A frame without hard ranges has one, the frame
// itself: its first chunk starts it, and the later chunks of a long frame do not (the last one
// reads the range after the frame's only one, which LOSAT takes as not that chunk's offset, as
// `core::subject_split::HardChunk::overlap_past_list` does).
impl TranslatedChunk {
    /// A chunk `[offset, offset + length)` of a frame without hard ranges (no `-db_hard_mask`).
    pub(crate) const fn new(offset: usize, length: usize) -> Self {
        use crate::core::subject_split::DBSEQ_CHUNK_OVERLAP;
        TranslatedChunk {
            offset,
            length,
            overlap: if offset == 0 { 0 } else { DBSEQ_CHUNK_OVERLAP },
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-264:
// Advance through overlapping subject chunks via `s_GetNextSubjectChunk` until the full range ends.
pub(crate) fn unmasked_translated_chunks(frame_length: usize) -> Vec<TranslatedChunk> {
    use crate::core::subject_split::{DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN};
    let mut chunks = Vec::new();
    let mut offset = 0;
    while offset < frame_length {
        if offset + MAX_DBSEQ_LEN < frame_length {
            chunks.push(TranslatedChunk::new(offset, MAX_DBSEQ_LEN));
            offset += MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP;
        } else {
            chunks.push(TranslatedChunk::new(offset, frame_length - offset));
            break;
        }
    }
    chunks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:813-826
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:489-502:
// Apply lowercase masking to both strands via `PackedSeqLocToMaskedQueryRegions`.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:689-707,825-826:
// both mask copies enter masked_ranges before BlastSeqBlkSetSeqRanges.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-215
// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:51-59
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
            for _ in 0..2 {
                nucleotide_ranges.push((left as i32, start as i32));
                left = index - 1;
            }
        } else {
            index += 1;
        }
    }
    nucleotide_ranges.push((left as i32, subject.len() as i32));
    translated_frame_ranges(&nucleotide_ranges, frame, negative_first_length)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:814-830 and :79
// The nucleotide ranges of a translated subject (`seq_ranges`, the complement of its masks) are
// mapped to protein ranges once per strand, by integer division by the codon length and with no
// frame offset: at the first positive context for the three positive frames, and at the first
// negative context for the three negative ones, in reverse order and reflected by the length of
// frame -1, which frames -2 and -3 keep even where they are one residue shorter.
/// The protein ranges of frame `frame` (+1..+3, -1..-3) of a subject whose nucleotide ranges are
/// `nucleotide_ranges` (in the order NCBI holds them), where `negative_first_length` is the number
/// of residues of frame -1: `[left / 3, right / 3)` each for a positive frame, and for a negative
/// one `[n - right / 3, n - left / 3)` with `n = negative_first_length`, last range first.
pub(crate) fn translated_frame_ranges(
    nucleotide_ranges: &[(i32, i32)],
    frame: i8,
    negative_first_length: usize,
) -> Vec<(i32, i32)> {
    if frame > 0 {
        nucleotide_ranges
            .iter()
            .map(|&(left, right)| (left / 3, right / 3))
            .collect()
    } else {
        let length = negative_first_length as i32;
        nucleotide_ranges
            .iter()
            .rev()
            .map(|&(left, right)| (length - right / 3, length - left / 3))
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:259-310:
// unsplit subjects keep backup.soft_ranges; split soft-masked subjects
// advance sm_index while right < chunk start, include ranges with left < end,
// clamp only the first left and last right, and skip SUBJECT_SPLIT_NO_RANGE.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:261-274:
// without soft masking, each chunk has one range [0, subject->length].
pub(crate) fn translated_chunk_scan_ranges(
    subject: &[u8],
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    mask_lowercase: bool,
) -> Result<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    let soft_masked = mask_lowercase && subject.iter().any(u8::is_ascii_lowercase);
    if !soft_masked {
        return unmasked_translated_chunks(frame_length)
            .into_iter()
            .map(|chunk| Ok((chunk, vec![(0, i32::try_from(chunk.length)?)])))
            .collect();
    }
    let frame_ranges = translated_scan_ranges(
        subject,
        frame,
        frame_length,
        negative_first_length,
        mask_lowercase,
    );
    soft_chunk_scan_ranges(&frame_ranges, frame_length)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:259-310 (cited at
// `translated_chunk_scan_ranges`)
/// The chunks of a frame of `frame_length` residues whose soft ranges are `frame_ranges`, each with
/// the soft ranges it scans (a chunk without one, NCBI's `SUBJECT_SPLIT_NO_RANGE`, is left out).
fn soft_chunk_scan_ranges(
    frame_ranges: &[(i32, i32)],
    frame_length: usize,
) -> Result<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    let mut sm_index = 0usize;
    let mut calls = Vec::new();
    for chunk in unmasked_translated_chunks(frame_length) {
        if chunk.offset == 0 && chunk.length == frame_length {
            calls.push((chunk, frame_ranges.to_vec()));
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321, 398-405 and
// c++/src/algo/blast/core/blast_engine.c:757-762, 814-830, 476-486, 221-311
// A database sequence of a search with `-db_soft_mask` or `-db_hard_mask` comes with the complement
// of its stored masks as its ranges (one range, the whole sequence, where it has none) and the mask
// kind. The frames of a translated subject map those ranges (`translated_frame_ranges`) and each
// frame is cut by them as a protein subject is: with soft masks, the ranges are the scan ranges of
// the frame's chunks of `MAX_DBSEQ_LEN`; with hard masks, each range is a chunk of its own (cut
// again past `MAX_DBSEQ_LEN`) that scans its letters from its start. NCBI's "no chunking" branch
// scans the whole frame, `[0, frame_length)`, from the chunk at the frame's start
// (blast_engine.c:150-152, 162, 264-268), so a whole frame one residue shorter than its strand's
// first frame, whose one chunk ends on its end sentinel, scans no word there. Where that chunk ends
// before the frame's end at a mask, LOSAT scans the chunk alone, as NCBI's normal chunks do for the
// near input, the mask shortened by the fewest bases after which a letter follows the chunk in
// every frame; the masked bases stay masked in the traceback and the statistics (an approved NCBI
// defect, PD-LOSAT-NCBI-DEFECTS exception 21: the Owner's answers Q1 = B and Q2 = B of M4b,
// 2026-11-21; NCBI's scan past the chunk gives rows that differ from the near input's, in TBLASTX
// the same row once for each word it reads in the mask, or a crash).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:145-215, 264-276
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-500, 525-545
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// Where the masks leave a frame an empty first part (a mask within the sequence's first two bases
// for the plus frames, or within the one or two bases after its last whole codon for the minus
// frames, with no other mask in the frame but at its other end; or an empty mask at 0), the empty
// part's chunk takes the hard branch, which allocates the subject's list of ranges, and the chunk
// from 0 after it takes the "no chunking" branch, which points that list at the backup's full
// range: NCBI searches the frame whole (or, shorter than the frame, as above), then frees that
// pointer when it restores the subject and crashes, in every run. The empty chunk adds nothing (its
// one word, if any, is a window or more behind every later hit on its diagonal), so LOSAT searches
// the frame's chunks as they come, the masked bases masked as for every other mask (an approved
// NCBI defect, PD-LOSAT-NCBI-DEFECTS exception 27: the Owner's answers Q4 = B of M4b, 2026-11-22,
// and Q6 = B, 2026-11-23).
/// The chunks of frame `frame` (+1..+3, -1..-3, of `frame_length` residues) of a database sequence
/// with the nucleotide ranges `nucleotide_ranges` (the complement of its stored masks,
/// `core::subject_split::complement_masks`), each with its scan ranges, in a search with
/// `-db_hard_mask` (`hard`) or `-db_soft_mask`. `negative_first_length` is the number of residues
/// of frame -1; `split` is the chunker, whose buffer of ranges is reused. A chunk of a hard range
/// may end one residue past a frame that is one residue shorter than its strand's first frame (the
/// frame's end sentinel), as NCBI's does. A hard chunk scans its own letters, `[0, length)`, and
/// the chunk of a whole frame the frame's, `[0, frame_length)`.
pub(crate) fn translated_db_mask_chunks(
    nucleotide_ranges: &[(i32, i32)],
    hard: bool,
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    split: &mut crate::core::subject_split::HardSplitState,
) -> Result<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    use crate::core::subject_split::DBSEQ_CHUNK_OVERLAP;
    let frame_ranges = translated_frame_ranges(nucleotide_ranges, frame, negative_first_length);
    if !hard {
        return soft_chunk_scan_ranges(&frame_ranges, frame_length);
    }
    let mut calls = Vec::new();
    if frame_length == 0 {
        return Ok(calls);
    }
    split.start(frame_length, frame_ranges, false);
    while let Some(chunk) = split.next_chunk(DBSEQ_CHUNK_OVERLAP) {
        // The protein programs' chunks have no residual: the scan starts at the chunk's start. The
        // "no chunking" chunk of a whole frame stops at the frame's end, not on its sentinel.
        let right = if chunk.full_range {
            chunk.length.min(frame_length)
        } else {
            chunk.length
        };
        calls.push((
            TranslatedChunk {
                offset: chunk.offset,
                length: chunk.length,
                overlap: chunk.overlap,
            },
            vec![(0, i32::try_from(right)?)],
        ));
    }
    split.end();
    Ok(calls)
}

#[cfg(test)]
#[path = "translated_chunks_tests.rs"]
mod tests;
