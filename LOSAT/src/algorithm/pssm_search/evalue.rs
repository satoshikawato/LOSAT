//! The statistics of a PSSM search's e-values and hit-saving cutoff (psiblast's PSSM rounds:
//! `-in_pssm` and every round after the first; RPS-BLAST's redo): Spouge's finite-size correction
//! with the Gumbel block, or, without it (`-searchsp` at any value, or the environment variable
//! `OLD_FSC`), Karlin-Altschul's simple forms on the context's effective search space.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:921-923
//! ```c
//!     /* FSCOLD: to switch back to original FSC, comment out the following line */
//!     use_old_fsc = getenv("OLD_FSC");
//!     if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
//! ```c++
//!     if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace]) {
//!         CNcbiEnvironment env;
//!         env.Set("OLD_FSC", "true");
//!         opt.SetEffectiveSearchSpace(args[kArgEffSearchSpace].AsInt8());
//!     }
//! ```

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
    /// ```c
    ///       if (sbp->gbp) {
    ///           /* Only try Spouge's method if gumbel parameters are available */
    ///           if (!isRPS) {
    ///               hsp->evalue =
    ///                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
    ///                                query_info->contexts[hsp->context].query_length,
    ///                                subject_length);
    ///           } else {
    ///               /* for RPS blast, query and subject is swapped */
    ///               hsp->evalue =
    ///                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
    ///                                subject_length,
    ///                                query_info->contexts[hsp->context].query_length);
    ///           }
    ///       } else {
    ///           /* Get effective search space from the query information block */
    ///           hsp->evalue =
    ///               BLAST_KarlinStoE_simple(score, kbp[kbp_context],
    ///                                query_info->contexts[hsp->context].eff_searchsp);
    ///       }
    /// ```
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
    /// ```c
    ///          Int8 searchsp;
    ///          Int4 new_cutoff = 1;
    ///          double evalue = options->expect_value;
    ///          ...
    ///          searchsp = query_info->contexts[context].eff_searchsp;
    ///          ...
    ///          if (sbp->gbp && sbp->gbp->filled) {
    ///              ...
    ///              new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
    ///                          query_info->contexts[context].query_length,
    ///                          avg_subject_length);
    ///          } else {
    ///              BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
    ///          }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4097-4129
    /// ```c
    ///    Int4  s = *S, es;
    ///    ...
    ///    es = 1;
    ///    esave = e;
    ///    if (e > 0.)
    ///    {
    ///         ...
    ///         es = BlastKarlinEtoS_simple(e, kbp, searchsp);
    ///    }
    ///    /*
    ///    Pick the larger cutoff score between the user's choice
    ///    and that calculated from the value of E.
    ///    */
    ///    if (es > s) {
    ///       s_changed = TRUE;
    ///       *S = s = es;
    ///    }
    /// ```
    /// `BlastKarlinEtoS_simple` is `cutoff_score_from_evalue`; its value is raised to the
    /// `new_cutoff` of 1 (`-evalue` is above 0: NCBI's options validation refuses the others).
    /// A block with a negative value keeps the cutoff of 1:
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4101-4102
    /// ```c
    ///    if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.)
    ///       return 1;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4054-4057
    /// ```c
    ///    if (Lambda < 0. || K < 0. || H < 0.0)
    ///    {
    ///       return BLAST_SCORE_MIN;
    ///    }
    /// ```
    /// (`es` is then `BLAST_SCORE_MIN`, not above the cutoff of 1.)
    pub(crate) fn hit_cutoff(&self, prelim_evalue: f64, kbp: &KarlinParams, m: i32, n: i32) -> i32 {
        match self.gumbel {
            Some(gumbel) => blast_spouge_etos(prelim_evalue, kbp, gumbel, m, n),
            None if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 => 1,
            None => crate::algorithm::tblastx::ncbi_cutoffs::cutoff_score_from_evalue(
                prelim_evalue,
                self.search_space,
                kbp,
            )
            .max(1),
        }
    }
}

/// `BLAST_KarlinStoE_simple`: the e-value of `score` in the search space `searchsp`, -1 for a
/// block without valid values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4156-4171
/// ```c
/// double
/// BLAST_KarlinStoE_simple(Int4 S,
///       Blast_KarlinBlk* kbp,
///       Int8  searchsp)   /* size of search space. */
/// {
///    double   Lambda, K, H; /* parameters for Karlin statistics */
///
///    Lambda = kbp->Lambda;
///    K = kbp->K;
///    H = kbp->H;
///    if (Lambda < 0. || K < 0. || H < 0.) {
///       return -1.;
///    }
///
///    return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);
/// }
/// ```
/// `kbp->logK` is `log(kbp->K)` in every block a PSSM search reads (psiblast_aux_priv.cpp:121,
/// blast_kappa.c:2127).
fn karlin_stoe_simple(score: i32, kbp: &KarlinParams, searchsp: i64) -> f64 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    searchsp as f64 * (-kbp.lambda * score as f64 + kbp.k.ln()).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ProteinScoringSpec;
    use crate::stats::lookup_protein_gumbel_params;

    /// BLOSUM62 11/1's gapped block (blast_stat.c's table: Lambda 0.267, K 0.041, H 0.14).
    fn blosum62_gapped() -> KarlinParams {
        KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        }
    }

    #[test]
    fn the_gumbel_block_gives_spouges_values() {
        let gumbel =
            lookup_protein_gumbel_params(&ProteinScoringSpec::default(), 3_234_567).unwrap();
        let stats = EvalueStats::spouge(&gumbel);
        let kbp = blosum62_gapped();
        for (score, m, n) in [(48, 302, 250), (25, 302, 10), (120, 87, 1_400)] {
            assert_eq!(
                stats.e_value(score, &kbp, m, n).to_bits(),
                blast_spouge_stoe(score, &kbp, &gumbel, m, n).to_bits()
            );
            // RPS-BLAST's swapped roles reach Spouge's function as given.
            assert_eq!(
                stats.e_value(score, &kbp, n, m).to_bits(),
                blast_spouge_stoe(score, &kbp, &gumbel, n, m).to_bits()
            );
        }
        for (evalue, m, n) in [(50.0, 302, 10), (10.0, 87, 1_400)] {
            assert_eq!(
                stats.hit_cutoff(evalue, &kbp, m, n),
                blast_spouge_etos(evalue, &kbp, &gumbel, m, n)
            );
        }
    }

    // blast_stat.c:4157-4171: `searchsp * exp(-Lambda S + log K)`, whatever the lengths; -1 for
    // a block with a negative value.
    #[test]
    fn without_the_gumbel_block_the_e_value_reads_the_search_space() {
        let kbp = blosum62_gapped();
        let stats = EvalueStats::new(None, 10_000);
        let e = stats.e_value(48, &kbp, 302, 250);
        assert_eq!(e, stats.e_value(48, &kbp, 1, 1));
        assert!((e - 0.0011139502021664525).abs() < 1e-18, "{e:e}");
        assert_eq!(
            e.to_bits(),
            (10_000f64 * (-0.267 * 48f64 + 0.041f64.ln()).exp()).to_bits()
        );
        let max = EvalueStats::new(None, i64::MAX);
        assert_eq!(
            max.e_value(48, &kbp, 302, 250),
            9_223_372_036_854_775_807f64 * (-0.267 * 48f64 + 0.041f64.ln()).exp()
        );
        for invalid in [
            KarlinParams {
                lambda: -1.0,
                ..kbp
            },
            KarlinParams { k: -1.0, ..kbp },
            KarlinParams { h: -1.0, ..kbp },
        ] {
            assert_eq!(stats.e_value(48, &invalid, 302, 250), -1.0);
        }
    }

    // blast_parameters.c:915-945 and blast_stat.c:4097-4129: `ceil(ln(K searchsp / E) / Lambda)`,
    // at least 1 (`new_cutoff`), whatever the lengths.
    #[test]
    fn without_the_gumbel_block_the_hit_cutoff_reads_the_search_space() {
        let kbp = blosum62_gapped();
        let cutoff = |search_space, evalue| {
            EvalueStats::new(None, search_space).hit_cutoff(evalue, &kbp, 302, 10)
        };
        // ln(0.041 * 10000 / 10) / 0.267 = 13.2...
        assert_eq!(cutoff(10_000, 10.0), 14);
        // A search space of 1 (`-searchsp 1`) gives a negative score: the cutoff stays at 1.
        assert_eq!(cutoff(1, 10.0), 1);
        // The largest `-searchsp`: ln(0.041 * 2^63 / 10) / 0.267 = 142.9...
        assert_eq!(cutoff(i64::MAX, 10.0), 143);
        assert_eq!(
            EvalueStats::new(None, 10_000).hit_cutoff(10.0, &kbp, 1, 1_000_000),
            14
        );
        // blast_stat.c:4054-4057, 4101-4102: a block with a negative value keeps the cutoff 1.
        for invalid in [
            KarlinParams {
                lambda: -1.0,
                ..kbp
            },
            KarlinParams { k: -0.5, ..kbp },
            KarlinParams { h: -0.14, ..kbp },
        ] {
            assert_eq!(
                EvalueStats::new(None, i64::MAX).hit_cutoff(10.0, &invalid, 302, 10),
                1
            );
        }
        // The largest search space with an E-value below kSmallFloat: the quotient is past the
        // double, the cutoff the finite ceil((ln K + ln searchsp - ln E) / Lambda)
        // (`cutoff_score_from_evalue`, line PS DECISIONS 28).
        assert_eq!(cutoff(i64::MAX, 1e-300), 2713);
    }
}
