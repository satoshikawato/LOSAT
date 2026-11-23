//! The statistics of a PSSM search's e-values and hit-saving cutoff (psiblast's PSSM rounds:
//! `-in_pssm` and every round after the first; RPS-BLAST's redo): Spouge's finite-size correction
//! with the Gumbel block, or, without it (`-searchsp` at any value, or the environment variable
//! `OLD_FSC`), Karlin-Altschul's simple forms on the context's effective search space.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:921-923
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307

use crate::stats::{
    blast_spouge_etos, blast_spouge_stoe, BlastGumbelBlk, KarlinParams, SpougeCutoffOverflow,
};
use crate::utils::glibc_fma_libm::GlibcLibm;

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
            None => karlin_stoe_simple(score, kbp, self.search_space),
        }
    }

    /// The hit-saving cutoff `cutoff_score_max` before any scale: `prelim_evalue` is
    /// `params->prelim_evalue` (`-evalue`, times 5 for the composition-based statistics' modes 2
    /// and 3 only with the Gumbel block), `m` the query's length and `n` the subject length of the
    /// cutoff (`avg_subject_length`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:915-945
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4097-4129
    /// `BlastKarlinEtoS_simple` is `cutoff_score_from_evalue`; its value is raised to the
    /// `new_cutoff` of 1 (`-evalue` is above 0: NCBI's options validation refuses the others).
    /// A block with a negative value keeps the cutoff of 1:
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4101-4102
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4054-4057
    /// (`es` is then `BLAST_SCORE_MIN`, not above the cutoff of 1.) With the Gumbel block, a search
    /// for the score that passes NCBI's `Int4` is `SpougeCutoffOverflow` (`blast_spouge_etos`).
    pub(crate) fn hit_cutoff(
        &self,
        prelim_evalue: f64,
        kbp: &KarlinParams,
        m: i32,
        n: i32,
    ) -> Result<i32, SpougeCutoffOverflow> {
        Ok(match self.gumbel {
            Some(gumbel) => blast_spouge_etos(prelim_evalue, kbp, gumbel, m, n)?,
            None if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 => 1,
            None => crate::algorithm::tblastx::ncbi_cutoffs::cutoff_score_from_evalue(
                prelim_evalue,
                self.search_space,
                kbp,
            )
            .max(1),
        })
    }
}

/// `BLAST_KarlinStoE_simple`: the e-value of `score` in the search space `searchsp`, -1 for a
/// block without valid values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4156-4171
/// `kbp->logK` is the natural log of `kbp->K` in every block a PSSM search reads
/// (psiblast_aux_priv.cpp:121, blast_kappa.c:2127).
fn karlin_stoe_simple(score: i32, kbp: &KarlinParams, searchsp: i64) -> f64 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    searchsp as f64 * (-kbp.lambda * score as f64 + kbp.k.glibc_ln()).glibc_exp()
}

#[cfg(test)]
#[path = "evalue_tests.rs"]
mod tests;
