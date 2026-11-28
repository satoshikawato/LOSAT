//! The statistics of a PSSM search's e-values and hit-saving cutoff (psiblast's PSSM rounds:
//! `-in_pssm` and every round after the first; RPS-BLAST's redo): Spouge's finite-size correction
//! with the Gumbel block, or, without it (`-searchsp` at any value, or the environment variable
//! `OLD_FSC`), Karlin-Altschul's simple forms on the context's effective search space.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:921-923
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307

use crate::stats::karlin::evalue_from_raw_score;
use crate::stats::{blast_spouge_etos, blast_spouge_stoe, BlastGumbelBlk, KarlinParams};

/// What the e-values and the hit-saving cutoff of a PSSM search read besides the Karlin block:
/// the Gumbel block, and without it the context's effective search space.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EvalueStats<'a> {
    /// `sbp->gbp`: `None` with `-searchsp` or `OLD_FSC`.
    pub(crate) gumbel: Option<&'a BlastGumbelBlk>,
    /// `query_info->contexts[context].eff_searchsp` (NCBI's `Int8`). `e_value` and `hit_cutoff`
    /// read it only without the Gumbel block; other readers (the Smith-Waterman redo's
    /// `s_EvalueFromScore`, redo_alignment.c:944) may read it with the block too.
    pub(crate) search_space: i64,
}

impl<'a> EvalueStats<'a> {
    /// The Gumbel block, if any, and the context's effective search space.
    pub(crate) fn new(gumbel: Option<&'a BlastGumbelBlk>, search_space: i64) -> Self {
        Self {
            gumbel,
            search_space,
        }
    }

    /// With the Gumbel block, for a caller that has no search space (0: `e_value` and
    /// `hit_cutoff` do not read it).
    pub(crate) fn spouge(gumbel: &'a BlastGumbelBlk) -> Self {
        Self::new(Some(gumbel), 0)
    }

    /// One HSP's e-value from its score: `m` and `n` are Spouge's lengths in NCBI's order for the
    /// program (RPS-BLAST passes the subject's first, `isRPS`), `kbp` the block that NCBI reads
    /// (`kbp_gap`: `kbp_gap_psi` in a PSSM search, the redo copy's scaled block in the redo).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1872-1890
    pub(crate) fn e_value(&self, score: i32, kbp: &KarlinParams, m: i32, n: i32) -> f64 {
        match self.gumbel {
            Some(gumbel) => blast_spouge_stoe(score, kbp, gumbel, m, n),
            None => evalue_from_raw_score(score, kbp, self.search_space as f64),
        }
    }

    /// The hit-saving cutoff `cutoff_score_max` before any scale: `prelim_evalue` is
    /// `params->prelim_evalue` (`-evalue`, times 5 for the composition-based statistics' modes 2
    /// and 3 only with the Gumbel block), `m` the query's length and `n` the subject length of the
    /// cutoff (`avg_subject_length`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:915-945
    /// Without the Gumbel block, `BLAST_Cutoffs` from the `new_cutoff` of 1
    /// (`stats::cutoffs::cutoff_score`): a block with a negative value keeps the cutoff of 1
    /// (`BLAST_Cutoffs` returns at once for a block of -1s, and `BlastKarlinEtoS_simple` gives
    /// `BLAST_SCORE_MIN` for any negative Lambda, K or H), as does an e-value not above 0
    /// (`-evalue -nan`, which NCBI's option check lets pass).
    pub(crate) fn hit_cutoff(&self, prelim_evalue: f64, kbp: &KarlinParams, m: i32, n: i32) -> i32 {
        match self.gumbel {
            Some(gumbel) => blast_spouge_etos(prelim_evalue, kbp, gumbel, m, n),
            None => crate::stats::cutoffs::cutoff_score(
                prelim_evalue,
                kbp,
                self.search_space,
                None,
                crate::stats::cutoffs::CutoffConversion::PROTEIN,
            ),
        }
    }
}

#[cfg(test)]
#[path = "evalue_tests.rs"]
mod tests;
