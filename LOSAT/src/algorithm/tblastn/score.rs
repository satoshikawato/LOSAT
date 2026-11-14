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
use crate::algorithm::blastp::extension::{
    extend_one_hit_for, BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult,
};
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_in_window, blastp_get_start_for_gapped_alignment,
    blastp_score_only_gapped_alignment_with_scratch, gapped_start_int4, BlastpGapAlignResult,
    BlastpGappedAlignmentMode, BlastpScoreOnlyGappedAlignment, GapAlignScratch,
};
use crate::algorithm::pssm_search::gapalign::{
    pssm_gapped_alignment_with_traceback, pssm_get_offsets_for_gapped_alignment,
    pssm_get_start_for_gapped_alignment, pssm_score_only_gapped_alignment, PssmGapScratch,
};
use crate::algorithm::pssm_search::rows::{PssmRows, QueryPssm};
use crate::algorithm::pssm_search::ungapped::{pssm_extend_one_hit, pssm_extend_two_hit};
use crate::config::ScoringMatrix;
use crate::utils::matrix::BLASTAA_SIZE;
use anyhow::Result;

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
    /// The memory of one traceback (`GapAlignScratch` for the matrix's kernel).
    type TracebackScratch;
    /// Whether the window of `gapped_start_int4` reads the query's letters: the matrix's does
    /// (`matrix->data[*query_var][*subject_var]`), a PSSM's reads its rows instead
    /// (`pssm->data[index1][*subject_var]`, blast_gapalign.c:3413-3416).
    const WINDOW_READS_QUERY: bool;

    /// A traceback's memory, made where the traceback of an HSP makes it.
    fn traceback_scratch() -> Self::TracebackScratch;

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

    /// `gapped_start` of a one-hit HSP whose Int4 `length` is negative (`-window_size 0`,
    /// aa_ungapped.c:1054,1083), which NCBI passes as Uint4 (blast_gapalign.c:3394): its window of
    /// HSP_MAX_WINDOW letters is scored. `query_window` and `subject_window` are the bytes NCBI
    /// reads there, from the query's `q_start` and the subject's `s_start` on, in NCBI's buffers
    /// (past a context or a frame: its NULLB, then the next one's letters;
    /// `search_gapped::XOneHitBuffers`).
    fn gapped_start_int4(
        self,
        query_window: &[u8; 11],
        subject_window: &[u8; 11],
        q_start: usize,
        length: i32,
        matrix: ScoringMatrix,
    ) -> Result<usize>;

    /// `score_only` of a gapped start at or past the end of its chunk (`s_start >=
    /// subject_length`, NCBI's `subject->length`), which only a one-hit HSP of negative length
    /// gives (`-window_size 0`): `subject` holds the bytes of NCBI's buffer from the chunk's first
    /// letter to `s_start`. NCBI's left extension reads back from `s_start`, and the right one
    /// does not run (blast_gapalign.c:4257-4258,4295).
    #[allow(clippy::too_many_arguments)]
    fn score_only_past_chunk(
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
        scratch: &mut Self::Scratch,
    ) -> Result<BlastpScoreOnlyGappedAlignment>;

    /// `s_BlastAaExtendOneHit` (aa_ungapped.c:1020-1083) of the one-hit word finder
    /// (`-window_size 0`, `XTwoHit::x_pairs_one_hit`) at the word `(q_off, s_off)`.
    #[allow(clippy::too_many_arguments)]
    fn one_hit(
        self,
        query: &[u8],
        subject: &[u8],
        q_off: usize,
        s_off: usize,
        dropoff: i32,
        word_size: usize,
        matrix: ScoringMatrix,
    ) -> Option<BlastpOneHitUngappedResult>;

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
        scratch: &mut Self::TracebackScratch,
        fence_hit: &mut bool,
    ) -> Option<BlastpGapAlignResult>;
}

/// The scores of the matrix (`sbp->matrix`): tblastn's calls.
#[derive(Clone, Copy)]
pub(super) struct MatrixScore;

impl TbnScore for MatrixScore {
    type Scratch = GapAlignScratch;
    type Walks = std::convert::Infallible;
    type TracebackScratch = GapAlignScratch;
    const WINDOW_READS_QUERY: bool = true;

    #[inline(always)]
    fn traceback_scratch() -> Self::TracebackScratch {
        GapAlignScratch::new()
    }

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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083 and
    // c++/src/algo/blast/core/blast_gapalign.c:3394,3401,3405-3413,3419-3420,3435-3437
    // ```c
    //     init_hit_width = q_right_off - q_left_off + 1;
    //     *hsp_len = left_disp + right_disp + init_hit_width;
    // ...
    //    const BlastScoreBlk* sbp, Uint4 q_start, Uint4 q_length,
    //     if (q_length <= HSP_MAX_WINDOW) {
    //     hsp_end = q_start + HSP_MAX_WINDOW;
    // ...
    //             score += sbp->matrix->data[*query_var][*subject_var];
    // ...
    //     hsp_end = q_start + MIN(q_length, s_length);
    //     for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
    // ...
    //        max_offset -= HSP_MAX_WINDOW/2;
    //     else
    //        max_offset = q_start;
    // ```
    // A negative Int4 length read as Uint4 is above HSP_MAX_WINDOW, and q_start + MIN(q_length,
    // s_length) wraps back to q_start + length, so only the first window is scored: its sum over
    // the bytes NCBI reads gives q_start + 5 when above 0, else q_start (`gapped_start_int4` of
    // `blastp::gapalign`, which BLASTX shares). BLASTP's own function for this
    // (`blastp::blast_engine`) reads its sequences with the letters past them as the sentinel
    // (PD-LOSAT-NCBI-DEFECTS 1.3 item 3) and stays where BLASTP's gapped stage inlines it.
    fn gapped_start_int4(
        self,
        query_window: &[u8; 11],
        subject_window: &[u8; 11],
        q_start: usize,
        length: i32,
        matrix: ScoringMatrix,
    ) -> Result<usize> {
        debug_assert!(length < 0);
        Ok(gapped_start_int4(
            query_window,
            subject_window,
            q_start,
            matrix,
        ))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4171-4189,4257-4258,4295
    // ```c
    //    if (subject_length < MAX_SUBJECT_OFFSET) {
    //       *start_shift = 0;
    //       return;
    //    }
    // ...
    //    if (s_offset <= max_extension_left) {
    //       *start_shift = 0;
    //    } else {
    //       *start_shift = s_offset - max_extension_left;
    //       *subject_offset_ptr = max_extension_left;
    //    }
    //
    //    *subject_length_ptr =
    //       MIN(subject_length, s_offset + max_extension_right) - *start_shift;
    // ...
    //    AdjustSubjectRange(&s_length, &subject_length, q_length, query_length,
    //                       &subject_shift);
    // ...
    //    if (q_length < query_length && s_length < subject_length) {
    // ```
    // NCBI's `subject->length` enters the alignment only through `AdjustSubjectRange` and the
    // test of the right extension. With `s_off >= subject->length`, the alignment on the letters
    // up to `s_off` (a `subject->length` of `s_off + 1`) gives NCBI's values whenever both lengths
    // are on the same side of MAX_SUBJECT_OFFSET: the same shift (it depends on `s_offset` alone),
    // a right extension that does not run in either case (`s_length` is not below the adjusted
    // length), and the same left extension. That holds for the chunks a gapped start can pass
    // (MAX_DBSEQ_LEN letters, followed by letters of their frame); for the matrix, the start of the
    // last chunk of a frame stays inside it (the NULLB after the frame scores BLAST_SCORE_MIN, so
    // a window that holds it gives `q_start`). The other case is refused as impossible.
    #[cold]
    fn score_only_past_chunk(
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
        scratch: &mut Self::Scratch,
    ) -> Result<BlastpScoreOnlyGappedAlignment> {
        const MAX_SUBJECT_OFFSET: usize = 90_000;
        anyhow::ensure!(
            s_start >= subject_length
                && subject.len() == s_start + 1
                && ((subject_length >= MAX_SUBJECT_OFFSET) == (subject.len() >= MAX_SUBJECT_OFFSET)),
            "TBLASTN -window_size 0: a gapped start past a chunk of {} letters at letter {} (internal error)",
            subject_length,
            s_start
        );
        Ok(blastp_score_only_gapped_alignment_with_scratch(
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
        ))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:785-790
    // ```c
    //                 score = s_BlastAaExtendOneHit(matrix, subject, query,
    //                                               subject_offset, query_offset,
    //                                               cutoffs->x_dropoff,
    //                                               &hsp_q, &hsp_s, &hsp_len,
    //                                               wordsize, use_pssm, &s_last_off);
    // ```
    #[inline(always)]
    fn one_hit(
        self,
        query: &[u8],
        subject: &[u8],
        q_off: usize,
        s_off: usize,
        dropoff: i32,
        word_size: usize,
        matrix: ScoringMatrix,
    ) -> Option<BlastpOneHitUngappedResult> {
        extend_one_hit_for(matrix, query, subject, q_off, s_off, dropoff, word_size)
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
        scratch: &mut Self::TracebackScratch,
        fence_hit: &mut bool,
    ) -> Option<BlastpGapAlignResult> {
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
            scratch,
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
    type TracebackScratch = PssmGapScratch;
    const WINDOW_READS_QUERY: bool = false;

    fn traceback_scratch() -> Self::TracebackScratch {
        PssmGapScratch::new()
    }

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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3394,3405-3413
    // ```c
    //    const BlastScoreBlk* sbp, Uint4 q_start, Uint4 q_length,
    // ...
    //     hsp_end = q_start + HSP_MAX_WINDOW;
    // ...
    //             score += sbp->psi_matrix->pssm->data[index1][*subject_var];
    // ```
    // A negative length is not below HSP_MAX_WINDOW as Uint4, so the first window of
    // HSP_MAX_WINDOW positions is scored (`pssm_get_start_for_gapped_alignment` takes the Int4
    // length) on the bytes NCBI reads (`subject_window`: past the frame its NULLB, whose score is
    // the PSSM's '-' column, then the next frame's letters). A window past the PSSM's last row
    // reads past NCBI's matrix (`data` has one row per query position): no result to port, so the
    // search stops (PD-LOSAT-NCBI-DEFECTS, rule 3).
    fn gapped_start_int4(
        self,
        _query_window: &[u8; 11],
        subject_window: &[u8; 11],
        q_start: usize,
        length: i32,
        _matrix: ScoringMatrix,
    ) -> Result<usize> {
        const HSP_MAX_WINDOW: usize = 11;
        if q_start + HSP_MAX_WINDOW > self.0.rows.len() {
            anyhow::bail!(
                "TBLASTN -in_pssm with -window_size 0: a one-hit alignment of negative length at query position {} would read past the PSSM (NCBI BLAST+ reads past its matrix there), which LOSAT does not support",
                q_start + 1
            );
        }
        Ok(pssm_get_start_for_gapped_alignment(
            self.0,
            subject_window,
            i32::try_from(q_start)?,
            length,
            0,
        ) as usize)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4245-4347 (quoted at
    // `pssm_search::gapalign::pssm_score_only_gapped_alignment`, which takes `subject->length`)
    #[cold]
    fn score_only_past_chunk(
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
        scratch: &mut Self::Scratch,
    ) -> Result<BlastpScoreOnlyGappedAlignment> {
        anyhow::ensure!(
            s_start >= subject_length && subject.len() == s_start + 1,
            "TBLASTN -window_size 0: a gapped start past a chunk of {} letters at letter {} (internal error)",
            subject_length,
            s_start
        );
        let gapped = pssm_score_only_gapped_alignment(
            self.0,
            subject,
            subject_length,
            q_start,
            s_start,
            gap_open,
            gap_extend,
            x_drop,
            &mut scratch.1,
        );
        Ok(BlastpScoreOnlyGappedAlignment {
            query_start: gapped.query_start,
            query_stop: gapped.query_stop,
            subject_start: gapped.subject_start,
            subject_stop: gapped.subject_stop,
            score: gapped.score,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:785-790 (quoted at
    // `pssm_search::ungapped::pssm_extend_one_hit`, the `use_pssm` branch)
    // ```c
    //                 score = s_BlastAaExtendOneHit(matrix, subject, query,
    //                                               subject_offset, query_offset,
    // ...
    //                                               wordsize, use_pssm, &s_last_off);
    // ```
    // The query offsets are the PSSM's row indices (one context at offset 0).
    fn one_hit(
        self,
        _query: &[u8],
        subject: &[u8],
        q_off: usize,
        s_off: usize,
        dropoff: i32,
        word_size: usize,
        _matrix: ScoringMatrix,
    ) -> Option<BlastpOneHitUngappedResult> {
        Some(pssm_extend_one_hit(
            &XOneHitRows(self.0),
            subject,
            i32::try_from(s_off).ok()?,
            i32::try_from(q_off).ok()?,
            dropoff,
            i32::try_from(word_size).ok()?,
        ))
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
        scratch: &mut Self::TracebackScratch,
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
            scratch,
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

/// No NCBI counterpart: a type of its own for the PSSM's rows in TBLASTN's one-hit word finder
/// (`PssmScore::one_hit`), so that `pssm_extend_one_hit` is compiled for it apart from BLASTP's PSSM
/// stage, which calls it for `QueryPssm` and keeps its code; the rows read are `QueryPssm`'s.
struct XOneHitRows<'a>(&'a QueryPssm);

impl PssmRows for XOneHitRows<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        self.0.row(pos)
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.0.len()
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
