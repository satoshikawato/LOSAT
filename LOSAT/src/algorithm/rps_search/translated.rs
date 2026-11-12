//! The frames of RPSTBLASTN's translated queries as the RPS engine reads them (series PR, item
//! PS12a, DECISIONS 38 D1): the six frames of each nucleotide query are the "subjects" of the
//! reversed search, one context each, with the frame's statistics. RPSTBLASTN's batch loop
//! (`algorithm::rpstblastn::search`) builds them from BLASTX's set-up of the query batch; nothing
//! here reads a BLASTX type. NCBI's `eBlastTypeRpsTblastn` is `PSSM_SUBJECT_MASK | eBlastTypeBlastx`.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:85
//! ```c
//!     eBlastTypeRpsTblastn    = (PSSM_SUBJECT_MASK | eBlastTypeBlastx),
//! ```

use anyhow::{bail, Result};

use super::engine::RpsSearch;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::KarlinParams;

/// The six contexts of a translated query (`NUM_FRAMES`).
pub(crate) const NUM_FRAMES: usize = 6;

/// The RPS search of one volume for the frames of a batch's queries: each query's hit list, in
/// the batch's order. The frames of a query are the subjects of the reversed search (NCBI's
/// frame loop of `s_RPSPreliminarySearchEngine`, blast_engine.c:757-846); the frames' results
/// combine by score before the traceback, which aligns each HSP on its frame.
pub(crate) fn search_volume(
    _search: &RpsSearch,
    _files: &RpsVolumeFiles,
    _frames: &[RpsFrameQuery<'_>],
) -> Result<Vec<Option<BlastpHitList>>> {
    bail!("the RPS search of translated queries is not supported by LOSAT's RPSTBLASTN yet")
}

/// One frame (context) of a user query as the RPS search reads it: NCBI's `subject->frame`,
/// `subject->sequence` and `subject->length` of the frame loop, and the frame's statistics.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-813
/// ```c
///     for (context=first_context; context<=last_context; context++) {
///         BlastHSPList* hsp_list_for_chunks = NULL;
///
///         if (kTranslatedSubject) {
///             Uint4 i;
///             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
///             subject->sequence = translation_buffer + frame_offsets[context] + 1;
///             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
///             if (subject->length > 0) stat_length = subject->length;
/// ```
#[derive(Debug, Clone, Copy)]
#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
pub(crate) struct RpsFrame<'a> {
    /// +1, +2, +3, -1, -2 or -3.
    pub(crate) frame: i8,
    /// The frame's letters between two sentinels as searched (SEG and lower-case masks); empty
    /// for a frame without letters (the other strand of `-strand`, a query shorter than the
    /// frame), which the search does not scan.
    pub(crate) aa_seq: &'a [u8],
    /// The frame's letters as read (between two sentinels) where the searched ones differ: the
    /// reports' identities count on them, as rpsblast's (`traceback::formatter_identities`).
    pub(crate) aa_seq_nomask: Option<&'a [u8]>,
    /// Whether the frame's Karlin-Altschul block exists (an invalid frame gets no HSP).
    pub(crate) valid: bool,
    /// `kbp_std[context]`: the frame's composition with BLOSUM62, the ideal block where its
    /// lambda is the larger (`check_ideal`).
    pub(crate) kbp_std: KarlinParams,
    /// `eff_searchsp` of the context.
    pub(crate) search_space: i64,
}

#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
impl RpsFrame<'_> {
    /// The frame's letters without the sentinels (`subject->length`).
    pub(crate) fn length(&self) -> usize {
        self.aa_seq.len().saturating_sub(2)
    }
}

/// The six frames of one user query, in NCBI's context order (+1, +2, +3, -1, -2, -3).
#[derive(Debug, Clone)]
#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
pub(crate) struct RpsFrameQuery<'a> {
    pub(crate) frames: [RpsFrame<'a>; NUM_FRAMES],
}

#[allow(
    dead_code,
    reason = "the frames' search (`search_volume`) reads them once the engine runs frames (item PS12a's second part)"
)]
impl RpsFrameQuery<'_> {
    /// `stat_length` of the preliminary e-values: the length of the last frame that has letters
    /// (frame -3 with both strands, +3 with `-strand plus`), not the query's length; 0 when no
    /// frame has letters.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:813
    /// ```c
    ///             if (subject->length > 0) stat_length = subject->length;
    /// ```
    pub(crate) fn prelim_stat_length(&self) -> usize {
        self.frames
            .iter()
            .rev()
            .map(RpsFrame::length)
            .find(|&length| length > 0)
            .unwrap_or(0)
    }
}

/// BLASTX's context of a frame within its query (`BLAST_FrameToContext` of blastx): +1, +2, +3
/// are 0, 1, 2; -1, -2, -3 are 3, 4, 5.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1211-1220
/// ```c
/// Int4 BLAST_FrameToContext(Int2 frame, EBlastProgramType program)
/// {
///     if (Blast_QueryIsTranslated(program) ||
///         Blast_SubjectIsTranslated(program)) {
///         ASSERT(frame >= -3 && frame <= 3 && frame != 0);
///         if (frame > 0) {
///             return frame - 1;
///         } else {
///             return 2 - frame;
///         }
/// ```
pub(crate) fn frame_context(frame: i32) -> usize {
    debug_assert!((-3..=3).contains(&frame) && frame != 0);
    if frame > 0 {
        (frame - 1) as usize
    } else {
        (2 - frame) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(frame: i8, letters: &[u8]) -> RpsFrame<'_> {
        RpsFrame {
            frame,
            aa_seq: letters,
            aa_seq_nomask: None,
            valid: !letters.is_empty(),
            kbp_std: KarlinParams::default(),
            search_space: 0,
        }
    }

    // blast_engine.c:813 (quoted at `prelim_stat_length`): the length of the last frame with
    // letters, whatever its validity; 0 when no frame has letters.
    #[test]
    fn the_preliminary_stat_length_is_the_last_frame_with_letters() {
        let (five, four, none): (&[u8], &[u8], &[u8]) =
            (&[0, 1, 2, 3, 4, 5, 0], &[0, 1, 2, 3, 4, 0], &[]);
        let query = |lengths: [&'static [u8]; 6]| RpsFrameQuery {
            frames: std::array::from_fn(|i| frame([1, 2, 3, -1, -2, -3][i], lengths[i])),
        };
        assert_eq!(
            query([five, five, four, five, five, four]).prelim_stat_length(),
            4
        );
        // -strand plus: the minus frames have no letters.
        assert_eq!(
            query([five, four, four, none, none, none]).prelim_stat_length(),
            4
        );
        assert_eq!(
            query([five, five, five, five, four, none]).prelim_stat_length(),
            4
        );
        let mut invalid = query([five, five, five, five, five, five]);
        invalid.frames[5].valid = false;
        assert_eq!(invalid.prelim_stat_length(), 5);
        assert_eq!(query([none; 6]).prelim_stat_length(), 0);
    }

    // blast_util.c:1211-1220 (quoted at `frame_context`).
    #[test]
    fn frames_have_blastxs_contexts() {
        let contexts: Vec<usize> = [1, 2, 3, -1, -2, -3].map(frame_context).to_vec();
        assert_eq!(contexts, [0, 1, 2, 3, 4, 5]);
    }
}
