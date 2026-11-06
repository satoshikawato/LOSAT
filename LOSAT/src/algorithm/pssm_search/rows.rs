//! The rows of a position-specific score matrix, as the PSSM kernels read them: one row of
//! `BLASTAA_SIZE` scores for each query position.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:131-144
//! ```c++
//!     // Assign the matrix scores/frequency ratios
//!     const size_t kQueryLength = pssm->GetPssm().GetNumColumns();
//!     score_blk->psi_matrix = SPsiBlastScoreMatrixNew(kQueryLength);
//!
//!     // Get the scores
//!     bool missing_scores = false;
//!     try {
//!         unique_ptr< CNcbiMatrix<int> > scores
//!             (CScorematPssmConverter::GetScores(*pssm));
//!         ...
//!         for (TSeqPos c = 0; c < scores->GetCols(); c++) {
//!             for (TSeqPos r = 0; r < scores->GetRows(); r++) {
//!                 score_blk->psi_matrix->pssm->data[c][r] = (*scores)(r, c);
//!             }
//!         }
//! ```
//! NCBI's `psi_matrix->pssm->data[position][residue]` (`Int4**`): the kernels index it by the
//! absolute query position and the subject residue (NCBIstdaa).

use crate::core::blast_util::BLASTAA_SIZE;

/// The source of the score rows of a PSSM search: the query's PSSM (psiblast `-in_pssm`), and
/// later the rows of other PSSM searches (tblastn `-in_pssm`, the profiles of RPS-BLAST).
pub(crate) trait PssmRows: Sync {
    /// The scores of query position `pos` against each residue.
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE];
    /// The number of rows (the query length).
    fn len(&self) -> usize;
}

/// The PSSM of a query: `rows[position][residue]`, as NCBI's `psi_matrix->pssm->data`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct QueryPssm {
    pub(crate) rows: Vec<[i32; BLASTAA_SIZE]>,
}

impl PssmRows for QueryPssm {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        &self.rows[pos]
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.rows.len()
    }
}
