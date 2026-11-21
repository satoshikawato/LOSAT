//! The chunks of a database sequence searched with hard ranges (`-db_hard_mask`): NCBI's
//! `SubjectSplitStruct` and `s_GetNextSubjectChunk` (core/blast_engine.c:122-310) with
//! `mask_type == eHardSubjMasking`, and the overlap of the merge after each chunk
//! (blast_engine.c:579-584), as arithmetic on the hard ranges.
//!
//! The one owner of the hard chunks (series NP, M4a-2b review S1): BLASTN's preliminary stage
//! searches these chunks (`algorithm::blastn::blast_engine::run::next_hard_chunk` turns each into
//! the engine's chunk), and LOSAT's refusals of `-db_hard_mask` read them
//! (`blastinput::db_mask::HardMaskShape::of_scan`), so the two cannot disagree. The chunks of a subject
//! without hard ranges (unmasked or soft-masked, `run::SubjectSplitState::next_chunk`) keep their
//! own copy of the same arithmetic on `&full_range`, for the default path's code: the one
//! function for both cost the chunk loop's closure 280,200 more instructions (14 per subject) in
//! a 20,000-subject megablast search (callgrind, release builds; review 2 S1 of M4a-2b), and
//! `run::tests::hard_masked_chunks_follow_ncbi_s_get_next_subject_chunk` holds both copies to one
//! C model of `s_GetNextSubjectChunk`. [`MAX_DBSEQ_LEN`] and [`DBSEQ_CHUNK_OVERLAP`] are the one
//! copy of NCBI's two macros for every program.

use crate::core::blast_encoding::COMPRESSION_RATIO;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
/// The longest chunk of a subject (NCBI's `MAX_DBSEQ_LEN`).
pub const MAX_DBSEQ_LEN: usize = 5_000_000;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:190-192
/// The overlap of a subject's consecutive chunks within one range (NCBI's `DBSEQ_CHUNK_OVERLAP`).
pub const DBSEQ_CHUNK_OVERLAP: usize = 100;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:197-218
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:280-306
// `get_data()` starts one `TSeqPos` before the first mask, so its `SSeqRange` pairs are
// (the spare value, start 0), (end 0, start 1), ..., (end n-1, the spare value after the
// last): with `tmp[0].left = 0` and the last `right` the length, the complement of the
// masks in the order they are stored, the `TSeqPos` values read as `Int4`. Unsorted or
// overlapping masks give inverted members, which the scan skips (`masksubj.inl:43-60`);
// nothing is sorted or merged (decision D3 of M4a-2).
/// The `seq_ranges` of a database sequence of `length` letters with the stored masks `masks`
/// (`[start, end)`, in the order stored): their complement in that order, one range more than
/// the masks, into `ranges`, which it replaces. A search with `-db_soft_mask` takes them as the
/// sequence's soft ranges, one with `-db_hard_mask` as its hard ranges (blast_engine.c:165-176).
pub fn complement_masks(masks: &[(u32, u32)], length: usize, ranges: &mut Vec<(i32, i32)>) {
    ranges.clear();
    let mut left = 0i32;
    for &(start, end) in masks {
        ranges.push((left, start as i32));
        left = end as i32;
    }
    ranges.push((left, length as i32));
}

/// One chunk of a subject's hard ranges ([`HardSplitState::next_chunk`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HardChunk {
    /// `backup->offset`: the chunk's first letter in the subject, a multiple of
    /// `COMPRESSION_RATIO` in BLASTN (`next` itself in the protein programs).
    pub offset: usize,
    /// `subject->length`: the chunk's letters.
    pub length: usize,
    /// `residual`: the letters before the chunk's hard range (`next % COMPRESSION_RATIO`),
    /// where its scan starts (`seq_ranges[0].left`).
    pub residual: usize,
    /// NCBI's "no chunking" branch: NCBI's `seq_ranges` of the chunk are the soft ranges, the
    /// full range `[0, subject length)`, also where the chunk is shorter (a mask that runs to the
    /// end after a chunk from the start, SURVEY §7.1 of M4a-2). LOSAT's search scans every chunk
    /// from its residual (`[0, length)` here): where the chunk is shorter than the subject that
    /// is the approved NCBI defect of `blastinput::db_mask::TAIL_MASK_SEARCHED_WORD`, below whose
    /// word `HardMaskShape::TailMask` refuses the sequence (series NP, M4a-2c).
    pub full_range: bool,
    /// The hard range the chunk is cut from (`hm_index` when it is cut).
    pub range: usize,
    /// The overlap that the merge after the chunk takes: 0 where the chunk starts at the left of
    /// the hard range at `hm_index`, else `chunk_overlap`.
    pub overlap: usize,
    /// The chunk ends the last hard range: NCBI reads the overlap's range one past the list
    /// (blast_engine.c:579-580; another sequence's mask end, or unset memory, SURVEY §7.3 of
    /// M4a-2). The overlap is then LOSAT's: the full range's end, which no chunk starts at, gives
    /// `chunk_overlap`, as the unmasked search's `(&full_range)[1]`.
    pub overlap_past_list: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:122-143
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:171-179
/// `SubjectSplitStruct` of a subject with hard ranges: its full range's end, its hard ranges
/// (the complement of its masks in the order stored, [`complement_masks`]), `hm_index` and
/// `next`. Its soft ranges are the full range (blast_engine.c:162-164). The first hard range
/// starts at 0 (`BlastSeqBlkSetSeqRanges` sets `tmp[0].left = 0`), so `next` starts at 0. A
/// search's worker keeps one, its ranges empty outside a hard-masked sequence.
#[derive(Debug, Default)]
pub struct HardSplitState {
    full_right: i32,
    /// The hard ranges of the sequence being cut (empty outside one; the buffer is kept).
    pub ranges: Vec<(i32, i32)>,
    hm_index: usize,
    next: i32,
    /// `is_nucleotide` (`Blast_ProgramIsNucleotide`): BLASTN's chunks start at a multiple of
    /// `COMPRESSION_RATIO`, those of the protein programs (`residual` 0) at the range's left.
    nucleotide: bool,
}

impl HardSplitState {
    /// The hard ranges `ranges` (never empty: the complement has one range more than the masks)
    /// of a subject of `subject_len` letters, from its first chunk; `nucleotide` is NCBI's
    /// `is_nucleotide` of the program (BLASTN; not the protein programs, whose letters are not
    /// packed: blast_engine.c:448, `Blast_ProgramIsNucleotide`).
    pub fn start(&mut self, subject_len: usize, ranges: Vec<(i32, i32)>, nucleotide: bool) {
        debug_assert!(!ranges.is_empty());
        self.full_right = subject_len as i32;
        self.ranges = ranges;
        self.hm_index = 0;
        self.next = 0;
        self.nucleotide = nucleotide;
    }

    /// The end of the subject's chunks: no hard ranges (the buffer is kept).
    pub fn end(&mut self) {
        self.ranges.clear();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, with
    // `mask_type == eHardSubjMasking`
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:578-584
    /// The next chunk of the subject (`None`: `SUBJECT_SPLIT_DONE`): one chunk per hard range (a
    /// range above `MAX_DBSEQ_LEN` in overlapping chunks of `MAX_DBSEQ_LEN` letters), from `next -
    /// next % 4` (from `next` in a protein program), of `right - offset` letters (0 for a hard
    /// range of no letters: a mask at 0, masks that touch), scanned from the residual (NCBI scans
    /// the full range in the "no chunking" branch, [`HardChunk::full_range`]); with the overlap of
    /// the merge after it ([`HardChunk::overlap`], [`HardChunk::overlap_past_list`]). NCBI's
    /// `Int4` arithmetic: the
    /// hard ranges of sorted masks never go back, so no length is negative (masks out of order
    /// are refused before the search, `HardMaskShape::Unsorted`); a negative one would be clamped
    /// to 0, where NCBI moves the diagonal offset back by it.
    pub fn next_chunk(&mut self, chunk_overlap: usize) -> Option<HardChunk> {
        if self.next >= self.full_right {
            return None;
        }
        let residual = if self.nucleotide {
            self.next % COMPRESSION_RATIO as i32
        } else {
            0
        };
        let offset = self.next - residual;
        let range = self.hm_index;
        let hard_right = self.ranges[range].1;
        let length = if offset + (MAX_DBSEQ_LEN as i32) < hard_right {
            self.next = offset + MAX_DBSEQ_LEN as i32 - chunk_overlap as i32;
            MAX_DBSEQ_LEN as i32
        } else {
            debug_assert!(
                hard_right >= offset,
                "a hard range that ends at {hard_right}, before its chunk's offset {offset}"
            );
            self.hm_index += 1;
            self.next = self
                .ranges
                .get(self.hm_index)
                .map_or(self.full_right, |&(left, _)| left);
            (hard_right - offset).max(0)
        };
        let full_range = offset == 0 && residual == 0 && self.next == self.full_right;
        let next_left = self.ranges.get(self.hm_index).map(|&(left, _)| left);
        Some(HardChunk {
            offset: offset as usize,
            length: length as usize,
            residual: residual as usize,
            full_range,
            range,
            overlap: if offset == next_left.unwrap_or(self.full_right) {
                0
            } else {
                chunk_overlap
            },
            overlap_past_list: next_left.is_none(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunks(masks: &[(u32, u32)], length: usize) -> Vec<HardChunk> {
        chunks_of(masks, length, true)
    }

    fn chunks_of(masks: &[(u32, u32)], length: usize, nucleotide: bool) -> Vec<HardChunk> {
        let mut ranges = Vec::new();
        complement_masks(masks, length, &mut ranges);
        let mut split = HardSplitState::default();
        split.start(length, ranges, nucleotide);
        std::iter::from_fn(|| split.next_chunk(DBSEQ_CHUNK_OVERLAP)).collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 579-584 (quoted
    // at `HardSplitState::next_chunk`)
    // The fields that the refusals read: the "no chunking" branch of a chunk from 0 before a mask
    // that runs to the end (after empty masks at 0 as well, review F1), the range of each chunk,
    // and the overlap read past the list after the last range's chunk (a short last range, and the
    // last full chunk of a long one, review F3).
    #[test]
    fn hard_chunks_name_their_range_and_the_overlap_read_past_the_list() {
        let found = chunks(&[(0, 0), (0, 0), (1000, 1200)], 1200);
        assert_eq!(found.len(), 3);
        assert_eq!((found[2].offset, found[2].length), (0, 1000));
        assert!(found[2].full_range && found[2].range == 2 && !found[2].overlap_past_list);
        let found = chunks(&[(1001, 1002)], 1200);
        assert_eq!(
            found
                .iter()
                .map(|c| (c.offset, c.length, c.residual, c.range, c.overlap_past_list))
                .collect::<Vec<_>>(),
            [(0, 1001, 0, 0, false), (1000, 200, 2, 1, true)]
        );
        let found = chunks(&[(1001, 1002)], 5_300_000);
        assert_eq!(
            found
                .iter()
                .map(|c| (c.offset, c.length, c.range, c.overlap, c.overlap_past_list))
                .collect::<Vec<_>>(),
            [
                (0, 1001, 0, 100, false),
                (1000, 5_000_000, 1, 100, false),
                (5_000_900, 299_100, 1, 100, true),
            ]
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 578-584 (quoted
    // at `HardSplitState::next_chunk`)
    // A protein program's chunks start at their hard range's left (`residual` 0), so the overlap
    // that the merge after a chunk reads is the next range's left against that offset; the first
    // mask running to the end takes the "no chunking" branch for the chunk before it.
    #[test]
    fn protein_chunks_start_at_their_range() {
        let found = chunks_of(&[(150, 190)], 400, false);
        assert_eq!(
            found
                .iter()
                .map(|c| (c.offset, c.length, c.residual, c.range, c.overlap))
                .collect::<Vec<_>>(),
            [(0, 150, 0, 0, 100), (190, 210, 0, 1, 100)]
        );
        assert!(!found[0].full_range && !found[1].full_range);
        let found = chunks_of(&[(300, 400)], 400, false);
        assert_eq!(found.len(), 1);
        assert!(found[0].full_range && found[0].length == 300);
        let found = chunks_of(&[(0, 400)], 400, false);
        assert_eq!(found.len(), 1);
        assert!(found[0].full_range && found[0].length == 0);
    }
}
