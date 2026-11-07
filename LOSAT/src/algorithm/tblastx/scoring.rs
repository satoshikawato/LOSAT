//! The substitution scores of one TBLASTX search: BLOSUM62, or the matrix that `-matrix` names.
//!
//! NCBI reference: c++/src/algo/blast/core/blast_setup.c:383-387
//! ```c
//!         sbp->read_in_matrix = TRUE;
//!         BLAST_ScoreSetAmbigRes(sbp, 'X');
//!         sbp->name = BLAST_StrToUpper(scoring_options->matrix);
//!     }
//!     status = Blast_ScoreBlkMatrixFill(sbp, get_path);
//! ```
//!
//! The search reads its matrix through [`TbxScore`], chosen once before the search
//! (`blast_engine::run_impl::run_in_pool`). BLOSUM62's instance ([`Blosum62Score`], no data)
//! calls what the search called before `-matrix`: `blosum62_score` and the paths built on
//! BLOSUM62's scores (the batch table of TBXLUT, the table walks of TBXUNGAP with the sentinel
//! walks of TBXUNGAPSENT, the AVX2 reevaluation). The other matrices' instance ([`MatrixScore`])
//! reads a table of the matrix: PAIRLUT builds its group tables with the matrix's scores and the
//! units of PAIRPAR run wherever `decide` takes the batch, the reference frame loop builds NCBI's
//! table with them (`build_ncbi_lookup_for_profile`), and both two-hit extensions take the table
//! walks of TBXUNGAP with a byte table of the matrix ([`TbxScore::table_two_hit`]).
//!
//! The other instance also runs every search of a word length other than 3 (`-word_size` 2 or
//! 4, NCBI's 2..=4; BLOSUM62 too): it holds the word length ([`TbxScore::word`]), which BLOSUM62's
//! instance fixes at 3. Its two-hit extensions take the table walks and the sentinel walks with
//! the word loop over that length (`two_hit_in2`, `two_hit_in4`, `two_hit_sent_in::<W>` for `W`
//! in 2..=4, the matrix's byte table), and the reference two-hit extension, whose loop over the
//! word reads the length, for the calls those do not cover. The unit path takes the words whose
//! group tables the 16-bit direct builder makes (`x_pair_lut::x_direct_covers`: 2 and 3
//! letters), its units the walks of their word length ([`TbxScore::sent_two_hit_of`]); a search
//! of 4-letter words takes the reference frame loop with NCBI's table of that word length
//! (`build_ncbi_lookup_for_profile`).

use super::extension::{
    extend_hit_two_hit_sent, extend_hit_two_hit_with, matrix_walk_table, two_hit_in, two_hit_in2,
    two_hit_in4, two_hit_sent_in, MatrixWalkTable, TwoHit,
};
use crate::config::ScoringMatrix;
use crate::utils::matrix::{blosum62_score, protein_score, BLASTAA_SIZE};

/// The scores of a search's matrix, `matrix[a][b]` of NCBI's score block.
///
/// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:847-849
/// ```c
///     for (i = 0; i < n; i++) {
///         score += matrix[q[i]][s[i]];
/// ```
pub(crate) trait TbxScore: Copy + Send + Sync {
    /// `true` for BLOSUM62's instance: the search keeps the paths built on BLOSUM62's scores.
    const BLOSUM62: bool;

    /// The score of the NCBISTDAA letters `a` (query) and `b` (subject).
    fn score(self, a: u8, b: u8) -> i32;

    /// The matrix.
    fn matrix(self) -> ScoringMatrix;

    /// The word length of the search's lookup table (`-word_size`): 3 for BLOSUM62's instance.
    ///
    /// NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:237
    /// ```c
    ///     lookup->word_length = opt->word_size;
    /// ```
    #[inline(always)]
    fn word(self) -> usize {
        3
    }

    /// Whether the matrix's extensions may take the walks that the sentinels stop
    /// (LOSAT_X_TBXUNGAPSENT): its byte table exists and stops at a sentinel.
    fn x_sent_walks(self) -> bool {
        false
    }

    /// `extend_hit_two_hit_with(self, &q_full[1..q_full.len() - 1], &s_full[1..s_full.len() - 1],
    /// .., false)` with the walks the sentinels stop.
    ///
    /// SAFETY: `self.x_sent_walks()`, `q_full` is a context buffer `sent_contexts` accepted and
    /// `s_full` a frame buffer `sent_subject` accepted.
    // NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
    #[inline(always)]
    unsafe fn sent_two_hit(
        self,
        q_full: &[u8],
        s_full: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> TwoHit {
        extend_hit_two_hit_with(
            self,
            &q_full[1..q_full.len() - 1],
            &s_full[1..s_full.len() - 1],
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
            false,
        )
    }

    /// `sent_two_hit` for a search of `W`-letter words, with the word length chosen by the caller
    /// instead of at every call: the unit path's extensions (`x_hit`, the batched `resolve`), whose
    /// units take the instance of the search's word length (the dispatches of x_tbx_hitbf.rs, such
    /// as `scan_unit_or_bf`, once per unit). BLOSUM62's instance has 3-letter words only.
    ///
    /// SAFETY: those of `sent_two_hit`, and `self.word() == W`.
    // NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (`word_size` W).
    #[inline(always)]
    unsafe fn sent_two_hit_of<const W: usize>(
        self,
        q_full: &[u8],
        s_full: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> TwoHit {
        // The units check it once per unit (`scan_unit_cells`).
        debug_assert_eq!(W, self.word());
        // SAFETY: the caller's.
        unsafe { self.sent_two_hit(q_full, s_full, s_left_off, s_right_off, q_right_off, x_drop) }
    }

    /// `extend_hit_two_hit_reference(self, ..)` (debug output off) from the table walks of
    /// TBXUNGAP with the matrix's byte table, or `None` when the reference has to run. BLOSUM62's
    /// instance takes `extend_hit_two_hit` instead (`extend_hit_two_hit_with`).
    // NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
    #[inline(always)]
    fn table_two_hit(
        self,
        _q_seq: &[u8],
        _s_seq: &[u8],
        _s_left_off: usize,
        _s_right_off: usize,
        _q_right_off: usize,
        _x_drop: i32,
    ) -> Option<TwoHit> {
        None
    }
}

/// BLOSUM62, read from `blosum62_score`'s table at its static address.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Blosum62Score;

impl TbxScore for Blosum62Score {
    const BLOSUM62: bool = true;

    #[inline(always)]
    fn score(self, a: u8, b: u8) -> i32 {
        blosum62_score(a, b)
    }

    #[inline(always)]
    fn matrix(self) -> ScoringMatrix {
        ScoringMatrix::Blosum62
    }

    #[inline(always)]
    fn x_sent_walks(self) -> bool {
        true
    }

    #[inline(always)]
    unsafe fn sent_two_hit(
        self,
        q_full: &[u8],
        s_full: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> TwoHit {
        // SAFETY: the caller's.
        unsafe {
            extend_hit_two_hit_sent(q_full, s_full, s_left_off, s_right_off, q_right_off, x_drop)
        }
    }
}

/// The table of one matrix's scores for the NCBISTDAA letters, as NCBI fills its score block
/// (`protein_score`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:1559-1592
/// ```c
///     for (i = 0; i < sbp->alphabet_size; i++) {
///         for (j = 0; j < sbp->alphabet_size; j++) {
///             matrix[i][j] = BLAST_SCORE_MIN;
/// ...
///             matrix[i][j] = NCBISM_GetScore((const SNCBIPackedScoreMatrix*) psm,
///                                            i, j);
/// ...
///         matrix[u_index][i] = matrix[c_index][i];
///         matrix[i][u_index] = matrix[i][c_index];
///         matrix[o_index][i] = matrix[x_index][i];
///         matrix[i][o_index] = matrix[i][x_index];
/// ```
pub(crate) struct MatrixScores {
    matrix: ScoringMatrix,
    /// `-word_size` (`TbxScore::word`).
    word: usize,
    scores: Box<[i32; BLASTAA_SIZE * BLASTAA_SIZE]>,
    /// The byte table of the two-hit walks (`TbxScore::table_two_hit`), a static of the matrix;
    /// `None` if a score does not fit it.
    walk: Option<&'static MatrixWalkTable>,
}

impl MatrixScores {
    /// The scores of `matrix` for a search of 3-letter words.
    pub(crate) fn new(matrix: ScoringMatrix) -> Self {
        Self::with_word(matrix, 3)
    }

    /// The scores of `matrix` for a search of `word`-letter words (`-word_size`).
    pub(crate) fn with_word(matrix: ScoringMatrix, word: usize) -> Self {
        let mut scores = Box::new([0i32; BLASTAA_SIZE * BLASTAA_SIZE]);
        for a in 0..BLASTAA_SIZE {
            for b in 0..BLASTAA_SIZE {
                scores[a * BLASTAA_SIZE + b] = protein_score(matrix, a as u8, b as u8);
            }
        }
        let walk = matrix_walk_table(matrix, |a, b| {
            MatrixScore {
                matrix,
                word,
                scores: &scores,
                walk: None,
            }
            .score(a, b)
        });
        Self {
            matrix,
            word,
            scores,
            walk,
        }
    }

    /// The scores as the search reads them.
    pub(crate) fn scorer(&self) -> MatrixScore<'_> {
        MatrixScore {
            matrix: self.matrix,
            word: self.word,
            scores: &self.scores,
            walk: self.walk,
        }
    }
}

/// A matrix other than BLOSUM62, or any matrix with a word length other than 3, read from its
/// [`MatrixScores`].
#[derive(Clone, Copy)]
pub(crate) struct MatrixScore<'a> {
    matrix: ScoringMatrix,
    word: usize,
    scores: &'a [i32; BLASTAA_SIZE * BLASTAA_SIZE],
    walk: Option<&'static MatrixWalkTable>,
}

impl TbxScore for MatrixScore<'_> {
    const BLOSUM62: bool = false;

    // The letters of a search are below BLASTAA_SIZE; other bytes are read as
    // `blosum62_score` reads them (`blosum62_score_ncbistdaa_direct`).
    #[inline(always)]
    fn score(self, a: u8, b: u8) -> i32 {
        if (a as usize) < BLASTAA_SIZE && (b as usize) < BLASTAA_SIZE {
            return self.scores[a as usize * BLASTAA_SIZE + b as usize];
        }
        protein_score(self.matrix, a, b)
    }

    fn matrix(self) -> ScoringMatrix {
        self.matrix
    }

    #[inline(always)]
    fn word(self) -> usize {
        self.word
    }

    // The walks take words of 2..=4 letters (`sent_two_hit`).
    #[inline(always)]
    fn x_sent_walks(self) -> bool {
        (2..=4).contains(&self.word) && self.walk.is_some_and(MatrixWalkTable::x_sentinels_stop)
    }

    #[inline(always)]
    unsafe fn sent_two_hit(
        self,
        q_full: &[u8],
        s_full: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> TwoHit {
        // SAFETY: the caller's; `x_sent_walks` holds, so the table exists and stops at a sentinel
        // (and the word has 2..=4 letters).
        let walks = unsafe {
            let t = self.walk.unwrap_unchecked();
            // The walks with the word loop over `self.word` letters (`TbxScore::word`), as in
            // `table_two_hit`.
            match self.word {
                3 => two_hit_sent_in::<3>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                2 => two_hit_sent_in::<2>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                4 => two_hit_sent_in::<4>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                _ => None,
            }
        };
        match walks {
            Some(r) => r,
            None => x_matrix_sent_refused(
                self,
                q_full,
                s_full,
                s_left_off,
                s_right_off,
                q_right_off,
                x_drop,
            ),
        }
    }

    // The walks of `W`-letter words (`two_hit_sent_in::<W>`) without the choice in
    // `sent_two_hit`; a `W` without walks (an instance the caller never runs) takes the reference
    // (`x_matrix_sent_refused`).
    #[inline(always)]
    unsafe fn sent_two_hit_of<const W: usize>(
        self,
        q_full: &[u8],
        s_full: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> TwoHit {
        debug_assert_eq!(self.word, W);
        // SAFETY: the caller's; `x_sent_walks` holds, so the table exists and stops at a sentinel.
        let walks = unsafe {
            let t = self.walk.unwrap_unchecked();
            match W {
                2 => two_hit_sent_in::<2>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                3 => two_hit_sent_in::<3>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                4 => two_hit_sent_in::<4>(
                    t,
                    q_full,
                    s_full,
                    s_left_off,
                    s_right_off,
                    q_right_off,
                    x_drop,
                ),
                _ => None,
            }
        };
        match walks {
            Some(r) => r,
            None => x_matrix_sent_refused(
                self,
                q_full,
                s_full,
                s_left_off,
                s_right_off,
                q_right_off,
                x_drop,
            ),
        }
    }

    #[inline(always)]
    fn table_two_hit(
        self,
        q_seq: &[u8],
        s_seq: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> Option<TwoHit> {
        // The walks with the word loop over `self.word` letters (`TbxScore::word`); another word
        // length takes the reference.
        // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:1104-1108 (`word_size`, the lookup
        // table's word length)
        // ```c
        //     /* find one beyond the position (up to word_size-1 letters to the right)
        //        that gives the largest starting score */
        //     /* Use "one beyond" to make the numbering consistent with how it's done
        //        for BlastAaExtendOneHit and the "Extend" functions called here. */
        //     for (i = 0; i < word_size; i++) {
        // ```
        let t = self.walk?;
        match self.word {
            3 => two_hit_in(
                t,
                q_seq,
                s_seq,
                s_left_off,
                s_right_off,
                q_right_off,
                x_drop,
            ),
            2 => two_hit_in2(
                t,
                q_seq,
                s_seq,
                s_left_off,
                s_right_off,
                q_right_off,
                x_drop,
            ),
            4 => two_hit_in4(
                t,
                q_seq,
                s_seq,
                s_left_off,
                s_right_off,
                q_right_off,
                x_drop,
            ),
            _ => None,
        }
    }
}

/// The calls `MatrixScore::sent_two_hit` leaves to `extend_hit_two_hit_with` (an X-drop outside
/// 1..=127, a word that does not fit), out of the callers' code.
// NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit`.
#[cold]
#[inline(never)]
fn x_matrix_sent_refused(
    sc: MatrixScore<'_>,
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> TwoHit {
    extend_hit_two_hit_with(
        sc,
        &q_full[1..q_full.len() - 1],
        &s_full[1..s_full.len() - 1],
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
        false,
    )
}

/// A matrix other than BLOSUM62 in the reports' identities and positives
/// (`report::set_displayed_identities`, `report::engine_counts`, `report::PairwiseRows`), read
/// from `protein_score`. The reports pass BLOSUM62 as [`Blosum62Score`].
impl TbxScore for ScoringMatrix {
    const BLOSUM62: bool = false;

    fn score(self, a: u8, b: u8) -> i32 {
        protein_score(self, a, b)
    }

    fn matrix(self) -> ScoringMatrix {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // No NCBI counterpart: the table of BLOSUM62 holds the scores the search read before
    // `-matrix`, for every two bytes.
    #[test]
    fn the_table_of_blosum62_holds_blosum62_score() {
        let table = MatrixScores::new(ScoringMatrix::Blosum62);
        let scorer = table.scorer();
        for a in 0..=u8::MAX {
            for b in 0..=u8::MAX {
                assert_eq!(scorer.score(a, b), Blosum62Score.score(a, b), "{a} {b}");
            }
        }
    }

    // NCBI reference: c++/src/util/tables/sm_pam30.c (the PAM30 scores of A/A and W/W) and
    // blast_stat.c:1559-1592 (the gap letter scores BLAST_SCORE_MIN).
    #[test]
    fn the_table_of_another_matrix_holds_its_scores() {
        use crate::utils::matrix::ncbistdaa::{A, GAP, W};
        let table = MatrixScores::new(ScoringMatrix::Pam30);
        let scorer = table.scorer();
        assert_eq!(scorer.score(A, A), 6);
        assert_eq!(scorer.score(W, W), 13);
        assert_eq!(scorer.score(A, W), -13);
        assert_eq!(scorer.score(GAP, A), crate::utils::matrix::BLAST_SCORE_MIN);
        for a in 0..BLASTAA_SIZE as u8 {
            for b in 0..BLASTAA_SIZE as u8 {
                assert_eq!(
                    scorer.score(a, b),
                    protein_score(ScoringMatrix::Pam30, a, b)
                );
            }
        }
    }
}
