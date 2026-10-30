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

use super::extension::{matrix_walk_table, two_hit_in, MatrixWalkTable, TwoHit};
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
    scores: Box<[i32; BLASTAA_SIZE * BLASTAA_SIZE]>,
    /// The byte table of the two-hit walks (`TbxScore::table_two_hit`), a static of the matrix;
    /// `None` if a score does not fit it.
    walk: Option<&'static MatrixWalkTable>,
}

impl MatrixScores {
    pub(crate) fn new(matrix: ScoringMatrix) -> Self {
        let mut scores = Box::new([0i32; BLASTAA_SIZE * BLASTAA_SIZE]);
        for a in 0..BLASTAA_SIZE {
            for b in 0..BLASTAA_SIZE {
                scores[a * BLASTAA_SIZE + b] = protein_score(matrix, a as u8, b as u8);
            }
        }
        let walk = matrix_walk_table(matrix, |a, b| {
            MatrixScore {
                matrix,
                scores: &scores,
                walk: None,
            }
            .score(a, b)
        });
        Self {
            matrix,
            scores,
            walk,
        }
    }

    /// The scores as the search reads them.
    pub(crate) fn scorer(&self) -> MatrixScore<'_> {
        MatrixScore {
            matrix: self.matrix,
            scores: &self.scores,
            walk: self.walk,
        }
    }
}

/// A matrix other than BLOSUM62, read from its [`MatrixScores`].
#[derive(Clone, Copy)]
pub(crate) struct MatrixScore<'a> {
    matrix: ScoringMatrix,
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
    fn table_two_hit(
        self,
        q_seq: &[u8],
        s_seq: &[u8],
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        x_drop: i32,
    ) -> Option<TwoHit> {
        two_hit_in(
            self.walk?,
            q_seq,
            s_seq,
            s_left_off,
            s_right_off,
            q_right_off,
            x_drop,
        )
    }
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
