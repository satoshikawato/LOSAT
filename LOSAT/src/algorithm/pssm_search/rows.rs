//! The rows of a position-specific score matrix, as the PSSM kernels read them: one row of
//! `BLASTAA_SIZE` scores for each query position.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:129-144
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
    /// Called by the extension to the right (`pssm_extend_right`) when it scored every position
    /// up to its bound without stopping: `end` is the query position after the last one scored,
    /// `letter` the subject's letter after it (`None` at the subject's end). Where NCBI's bound
    /// lies past this one (RPS-BLAST's concatenated length above its rows, `rps_search::RpsRows`)
    /// and the subject has a letter left, NCBI scores one more row: its score for `letter`, which
    /// the extension adds as NCBI's loop does. `None` (every other PSSM, whose bound is its rows)
    /// ends the extension there.
    #[inline(always)]
    fn past_end(&self, _end: i32, _letter: Option<u8>) -> Option<i32> {
        None
    }
    /// Called after the extension scored the row of `past_end` and went on: `end` is the query
    /// position after it, `subject_left` whether the subject has letters after it (NCBI's loop
    /// then reads past `past_end`'s row too).
    #[inline(always)]
    fn past_end_continues(&self, _end: i32, _subject_left: bool) {}
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
