//! Alignment algorithms for BLASTN
//!
//! This module contains:
//! - Greedy alignment algorithms for high-identity sequences
//! - Gapped extension algorithms with affine gap penalties
//! - Alignment statistics calculation
//! - Utility functions (GCD, coordinate conversion, etc.)

pub mod gapped;
pub mod greedy;
pub mod statistics;
pub mod utilities;
// No NCBI counterpart: the fences that a thread's alignments met (LOSAT_X_BNTBRANGE traceback,
// preliminary greedy windows).
pub(crate) mod x_fence;
// No NCBI counterpart: the preliminary greedy alignment on a window of the packed subject, one
// base per byte, or on a one-byte copy of the whole subject (see its doc for the NCBI references).
pub mod x_greedy_prelim;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:385-390
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:466-469
/// NCBI's pointer into a larger buffer with the `Int4` length NCBI gives with it: the pointer is
/// `buf` from `start`, the length `len`. NCBI's `Int4` offsets relative to the pointer (`at`) read
/// the buffer outside `0..len` too, before the pointer for a negative offset and past `len` for one
/// at `len` or more, as NCBI's pointer arithmetic does; a greedy extension reads `len - off`
/// letters to the right of an offset and `off` to its left, none for a count below 0
/// (`right_of`, `left_of`). Two of them:
/// - `ContextQuery`: the query of one context, `query + query_offset` with `query_length`, in the
///   concatenated query (every context, with the sentinel between two of them; NCBI's `sequence`
///   or `sequence_nomask`), read outside the context by an HSP from a hit that passes an end of
///   its context (remaining-gaps X1b-o);
/// - `AdjustedSubject`: the subject range of a traceback, `subject + start_shift` with
///   `adjusted_s_length` from `AdjustSubjectRange`, in the whole subject. On a subject of 90000
///   letters or more a traceback start more than 3000 letters past an end of its context makes
///   the start in the range negative (left of the context) or the range shorter than the start,
///   so that `len - off` is negative (right of it), and the shift can pass the subject's end
///   (remaining-gaps X1b-o review 2).
#[derive(Clone, Copy, Debug)]
pub struct SeqPointer<'a> {
    pub buf: &'a [u8],
    pub start: usize,
    pub len: i32,
}

/// NCBI's `query` and `query_length` of one context (`SeqPointer`).
pub type ContextQuery<'a> = SeqPointer<'a>;

/// NCBI's `adjusted_subject` and `adjusted_s_length` of a traceback (`SeqPointer`).
pub type AdjustedSubject<'a> = SeqPointer<'a>;

impl<'a> SeqPointer<'a> {
    /// The place in `buf` of NCBI's offset `off` from the pointer (`ptr + off`).
    #[inline(always)]
    pub fn at(&self, off: i32) -> usize {
        (self.start as isize + off as isize) as usize
    }

    /// What an extension to the right of offset `off` reads: the letters from `at(off)` and their
    /// count, NCBI's `q_avail` or `s_avail` (`len - off`), 0 for a count below 0.
    #[inline(always)]
    pub fn right_of(&self, off: i32) -> (&'a [u8], i32) {
        let count = (self.len - off).max(0);
        let at = self.at(off);
        (&self.buf[at..at + count as usize], count)
    }

    /// What an extension to the left of offset `off` reads, backward from its end: the `off`
    /// letters before `at(off)` and their count, NCBI's `q_off` or `s_off`, 0 for an offset below 0.
    #[inline(always)]
    pub fn left_of(&self, off: i32) -> (&'a [u8], i32) {
        let count = off.max(0);
        let at = self.at(off);
        (&self.buf[at - count as usize..at], count)
    }
}

// Re-export public APIs
pub use greedy::{
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_with_traceback,
    GreedyAlignScratch,
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3305 (BlastGetOffsetsForGappedAlignment)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3323-3389 (BlastGetStartForGappedAlignmentNucl)
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
// The traceback extension reads its matrix with `get_unchecked` for letters below 16 (BLASTNA,
// which `run_local` encodes first), so it is the crate's only.
pub use gapped::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    build_blastna_matrix,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
};
pub(crate) use gapped::{extend_gapped_heuristic_with_traceback_with_scratch, stats_from_edit_ops};
