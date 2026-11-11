//! Where the scores of a TBLASTN search come from: the matrix (every FASTA query) or the query's
//! PSSM (`-in_pssm`, PSI-TBLASTN; line PS, item PS7, DECISIONS 35 Q7). The stage functions that
//! score letters (`search_gapped::score_gapped_chunk`, `search_gapped::x_traceback_in_windows`, and
//! the two per-subject functions that carry the type) are generic over `TbnScore`; the word
//! finder's two-hit extension takes the PSSM's walks (`XPssmWalks`) through the walks that
//! `search_init::XTwoHit::pairs_on` is already generic over. The search picks the source once per
//! subject (`stage_d_pipeline::run_preliminary_subjects_into`) and once per subject list (the
//! traceback, `stage_d_pipeline::finish_local_mode0_no_sum_stats`), on a value that is the same for
//! the whole query batch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3410-3420
//! ```c
//!     for (index1=q_start; index1<hsp_end; index1++) {
//!         if (!(positionBased))
//!             score += sbp->matrix->data[*query_var][*subject_var];
//!         else
//!             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
//! ```
//! NCBI's code is the same for both, with a branch on `positionBased` (or `use_pssm`) where it
//! reads a score: `MatrixScore` holds tblastn's calls as they were, each `#[inline(always)]`, so
//! that the default instance is tblastn's code; `PssmScore` calls the PSSM kernels of
//! `algorithm::pssm_search`, which port the other branch.

use super::search_init::{XTwoHit, XUngapWalks};
use crate::algorithm::blastp::extension::BlastpTwoHitUngappedResult;
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_in_window, blastp_get_start_for_gapped_alignment,
    blastp_score_only_gapped_alignment_with_scratch, BlastpGapAlignResult,
    BlastpGappedAlignmentMode, BlastpScoreOnlyGappedAlignment, GapAlignScratch,
};
use crate::algorithm::pssm_search::gapalign::{
    pssm_gapped_alignment_with_traceback, pssm_get_offsets_for_gapped_alignment,
    pssm_get_start_for_gapped_alignment, pssm_score_only_gapped_alignment, PssmGapScratch,
};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::pssm_search::ungapped::pssm_extend_two_hit;
use crate::config::ScoringMatrix;

/// The score source of the stage functions of a TBLASTN search: the matrix (`MatrixScore`) or the
/// query's PSSM (`PssmScore`). Its methods are the calls of those functions that read a score;
/// the query offsets are those of the one context of a PSSM query (its row indices).
pub(super) trait TbnScore: Copy + Send + Sync {
    /// The gapped-alignment memory that a unit of the preliminary stage keeps
    /// (`GapAlignScratch` for the matrix's kernels).
    type Scratch: Send;
    /// The two-hit extension of the PSSM (`XPssmWalks`); none for the matrix, whose word finder
    /// picks its own walks.
    type Walks: XUngapWalks;

    /// The walks of the two-hit extension, or `None` for the matrix's own.
    fn x_walks(self) -> Option<Self::Walks>;

    /// The word finder of a unit (`XTwoHit::with_table`'s), made for this source.
    fn x_two_hit(two_hit: XTwoHit<'_>) -> XTwoHit<'_>;

    /// A unit's gapped-alignment memory from the one a thread kept (`XUnitBuffers`).
    fn unit_scratch(kept: GapAlignScratch) -> Self::Scratch;

    /// What the thread keeps of a unit's gapped-alignment memory.
    fn kept_scratch(scratch: Self::Scratch) -> GapAlignScratch;

    /// `BlastGetStartForGappedAlignment` of an ungapped HSP of `length` letters
    /// (blast_gapalign.c:3393-3439).
    fn gapped_start(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        length: usize,
        s_start: usize,
        matrix: ScoringMatrix,
    ) -> usize;

    /// The score-only gapped alignment of the preliminary stage (`s_BlastProtGappedAlignment`
    /// without restricted alignment, blast_gapalign.c:4209-4347).
    #[allow(clippy::too_many_arguments)]
    fn score_only(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        s_start: usize,
        matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        scratch: &mut Self::Scratch,
    ) -> BlastpScoreOnlyGappedAlignment;

    /// `BlastGetOffsetsForGappedAlignment` of the traceback (blast_gapalign.c:3248-3320).
    #[allow(clippy::too_many_arguments)]
    fn offsets(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        matrix: ScoringMatrix,
    ) -> Option<(usize, usize)>;

    /// `BLAST_GappedAlignmentWithTraceback` on a translated window (blast_gapalign.c:4599-4718):
    /// `subject` holds `subject_length` letters and then the window's sentinel; `fence_hit` is set
    /// when the alignment reached a `FENCE_SENTRY`.
    #[allow(clippy::too_many_arguments)]
    fn traceback(
        self,
        query: &[u8],
        subject: &[u8],
        subject_length: usize,
        q_start: usize,
        s_start: usize,
        matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        fence_hit: &mut bool,
    ) -> Option<BlastpGapAlignResult>;
}

/// The scores of the matrix (`sbp->matrix`): tblastn's calls.
#[derive(Clone, Copy)]
pub(super) struct MatrixScore;

impl TbnScore for MatrixScore {
    type Scratch = GapAlignScratch;
    type Walks = std::convert::Infallible;

    #[inline(always)]
    fn x_walks(self) -> Option<Self::Walks> {
        None
    }

    #[inline(always)]
    fn x_two_hit(two_hit: XTwoHit<'_>) -> XTwoHit<'_> {
        two_hit
    }

    #[inline(always)]
    fn unit_scratch(kept: GapAlignScratch) -> Self::Scratch {
        kept
    }

    #[inline(always)]
    fn kept_scratch(scratch: Self::Scratch) -> GapAlignScratch {
        scratch
    }

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
    // max_offset = BlastGetStartForGappedAlignment(query_tmp.sequence,
    //     subject->sequence, gap_align->sbp, q_start, length, s_start, length);
    #[inline(always)]
    fn gapped_start(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        length: usize,
        s_start: usize,
        matrix: ScoringMatrix,
    ) -> usize {
        blastp_get_start_for_gapped_alignment(
            query, subject, q_start, length, s_start, length, matrix,
        )
    }

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3953:
    // status = s_BlastProtGappedAlignment(...);
    #[inline(always)]
    fn score_only(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        s_start: usize,
        matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        scratch: &mut Self::Scratch,
    ) -> BlastpScoreOnlyGappedAlignment {
        blastp_score_only_gapped_alignment_with_scratch(
            query,
            subject,
            q_start,
            s_start,
            matrix,
            gap_open,
            gap_extend,
            x_drop,
            BlastpGappedAlignmentMode::Exact,
            scratch,
        )
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:436-449:
    // retval = BlastGetOffsetsForGappedAlignment(...);
    #[inline(always)]
    fn offsets(
        self,
        query: &[u8],
        subject: &[u8],
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        matrix: ScoringMatrix,
    ) -> Option<(usize, usize)> {
        super::search_gapped::blast_get_offsets_for_gapped_alignment_protein(
            query, subject, q_start, q_end, s_start, s_end, matrix,
        )
    }

    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:509-512
    // ```c
    //            BLAST_GappedAlignmentWithTraceback(program_number, query,
    //                  adjusted_subject, gap_align, score_params, q_start, s_start,
    //                  query_length, adjusted_s_length,
    //                  fence_hit);
    // ```
    #[inline(always)]
    fn traceback(
        self,
        query: &[u8],
        subject: &[u8],
        subject_length: usize,
        q_start: usize,
        s_start: usize,
        matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        fence_hit: &mut bool,
    ) -> Option<BlastpGapAlignResult> {
        let mut scratch = GapAlignScratch::new();
        blast_gapped_alignment_with_traceback_in_window(
            query,
            subject,
            subject_length,
            q_start,
            s_start,
            matrix,
            gap_open,
            gap_extend,
            x_drop,
            &mut scratch,
            Some(fence_hit),
        )
    }
}

/// The scores of the query's PSSM (`sbp->psi_matrix->pssm->data[query position][residue]`,
/// `positionBased`): the PSSM kernels of `algorithm::pssm_search`. The query letters are not
/// read for scores (the identities read them, `protein_identities_from_edit_ops`).
#[derive(Clone, Copy)]
pub(super) struct PssmScore<'a>(pub(super) &'a QueryPssm);

impl<'a> TbnScore for PssmScore<'a> {
    type Scratch = (GapAlignScratch, PssmGapScratch);
    type Walks = XPssmWalks<'a>;

    fn x_walks(self) -> Option<Self::Walks> {
        Some(XPssmWalks(self.0))
    }

    fn x_two_hit(two_hit: XTwoHit<'_>) -> XTwoHit<'_> {
        two_hit.x_with_walks()
    }

    fn unit_scratch(kept: GapAlignScratch) -> Self::Scratch {
        (kept, PssmGapScratch::new())
    }

    fn kept_scratch(scratch: Self::Scratch) -> GapAlignScratch {
        scratch.0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3410-3420 (quoted at
    // `pssm_search::gapalign::pssm_get_start_for_gapped_alignment`)
    fn gapped_start(
        self,
        _query: &[u8],
        subject: &[u8],
        q_start: usize,
        length: usize,
        s_start: usize,
        _matrix: ScoringMatrix,
    ) -> usize {
        pssm_get_start_for_gapped_alignment(
            self.0,
            subject,
            q_start as i32,
            length as i32,
            s_start as i32,
        ) as usize
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4245-4347 (quoted at
    // `pssm_search::gapalign::pssm_score_only_gapped_alignment`; `subject->length` is the unit's
    // letters, the slice)
    fn score_only(
        self,
        _query: &[u8],
        subject: &[u8],
        q_start: usize,
        s_start: usize,
        _matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        scratch: &mut Self::Scratch,
    ) -> BlastpScoreOnlyGappedAlignment {
        let gapped = pssm_score_only_gapped_alignment(
            self.0,
            subject,
            subject.len(),
            q_start,
            s_start,
            gap_open,
            gap_extend,
            x_drop,
            &mut scratch.1,
        );
        BlastpScoreOnlyGappedAlignment {
            query_start: gapped.query_start,
            query_stop: gapped.query_stop,
            subject_start: gapped.subject_start,
            subject_stop: gapped.subject_stop,
            score: gapped.score,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3320 (quoted at
    // `pssm_search::gapalign::pssm_get_offsets_for_gapped_alignment`)
    fn offsets(
        self,
        _query: &[u8],
        subject: &[u8],
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        _matrix: ScoringMatrix,
    ) -> Option<(usize, usize)> {
        let (q, s) = pssm_get_offsets_for_gapped_alignment(
            self.0,
            subject,
            i32::try_from(q_start).ok()?,
            i32::try_from(q_end).ok()?,
            i32::try_from(s_start).ok()?,
            i32::try_from(s_end).ok()?,
        )?;
        Some((usize::try_from(q).ok()?, usize::try_from(s).ok()?))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4599-4718 (quoted at
    // `pssm_search::gapalign::pssm_gapped_alignment_with_traceback`, which reports the fence)
    fn traceback(
        self,
        _query: &[u8],
        subject: &[u8],
        subject_length: usize,
        q_start: usize,
        s_start: usize,
        _matrix: ScoringMatrix,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        fence_hit: &mut bool,
    ) -> Option<BlastpGapAlignResult> {
        let traced = pssm_gapped_alignment_with_traceback(
            self.0,
            subject,
            subject_length,
            q_start,
            s_start,
            gap_open,
            gap_extend,
            x_drop,
            &mut PssmGapScratch::new(),
        );
        *fence_hit = traced.fence_hit;
        Some(BlastpGapAlignResult {
            query_start: traced.query_start,
            query_stop: traced.query_stop,
            subject_start: traced.subject_start,
            subject_stop: traced.subject_stop,
            score: traced.score,
            edit_script: traced.edit_script,
        })
    }
}

/// The two-hit extension of a PSSM query (`use_pssm` of `s_BlastAaExtendTwoHit`), as the walks of
/// `XTwoHit::pairs_on`: every extension is the PSSM's (`Some`), none falls back to the matrix's
/// reference (the finder is made with `XTwoHit::x_with_walks`). The query offsets of the word
/// finder are the PSSM's row indices (one context, at offset 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1104-1157 (quoted at
/// `pssm_search::ungapped::pssm_extend_two_hit`)
#[derive(Clone, Copy)]
pub(super) struct XPssmWalks<'a>(pub(super) &'a QueryPssm);

impl XUngapWalks for XPssmWalks<'_> {
    fn two_hit(
        self,
        _query: &[u8],
        subject: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
        word_size: usize,
    ) -> Option<Option<BlastpTwoHitUngappedResult>> {
        Some(Some(pssm_extend_two_hit(
            self.0,
            subject,
            s_left_off as i32,
            s_right_off as i32,
            q_right_off as i32,
            dropoff,
            word_size as i32,
        )))
    }
}

/// The matrix has no walks of its own here (`MatrixScore::x_walks`).
impl XUngapWalks for std::convert::Infallible {
    fn two_hit(
        self,
        _query: &[u8],
        _subject: &[u8],
        _s_left_off: usize,
        _s_right_off: usize,
        _q_right_off: usize,
        _dropoff: i32,
        _word_size: usize,
    ) -> Option<Option<BlastpTwoHitUngappedResult>> {
        match self {}
    }
}
