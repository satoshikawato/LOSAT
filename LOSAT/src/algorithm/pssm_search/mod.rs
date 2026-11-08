//! The search of a query given as a position-specific score matrix (PSSM): psiblast `-in_pssm`
//! (series PS, items PS4a, `-comp_based_stats 0`, and PS4b, the position-based composition-based
//! statistics of modes 1 to 3).
//!
//! The kernels read the PSSM's rows by absolute query position where blastp's read the matrix
//! row of the query letter (`PssmRows`): the lookup table (`lookup.rs`), the ungapped extensions
//! (`ungapped.rs`), the gapped start, score-only and traceback alignments (`gapalign.rs`), all
//! scalar as NCBI's PSSM code is. The stages (`prelim.rs`, `traceback0.rs`, `redo1.rs`) run them for the one
//! query of a PSSM search; blastp's per-run driver calls them where a PSSM search differs
//! (`blastp::blast_engine`, the `psi_search` hooks) and does the rest (the subject sets and
//! databases, the hit lists, the reports) as for blastp.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2826
//! ```c
//!    sbp->kbp = Blast_QueryIsPssm(program) ? sbp->kbp_psi : sbp->kbp_std;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:113-126
//! ```c
//!     if (Blast_QueryIsPssm(program_number)) {
//!         ...
//!                 Blast_KarlinBlkCopy(sbp->kbp_gap_psi[index],
//!                                     sbp->kbp_gap_std[index]);
//!         ...
//!         sbp->kbp_gap = sbp->kbp_gap_psi;
//! ```

pub(crate) mod gapalign;
pub(crate) mod lookup;
pub(crate) mod prelim;
pub(crate) mod redo1;
pub(crate) mod rows;
pub(crate) mod traceback0;
pub(crate) mod ungapped;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::stats::KarlinParams;

/// The statistical values of a PSSM's `finalData`, each `None` where NCBI's getter gives
/// `CPssm::kInvalidStat` (the member is absent, or 0.0).
///
/// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:101-155
/// ```c++
/// double
/// CPssm::GetLambda() const
/// {
///     if ( !CanGetFinalData() || !GetFinalData().CanGetLambda() ) {
///         return kInvalidStat;
///     }
///     return GetFinalData().GetLambda();
/// }
/// ...
/// const double CPssm::kInvalidStat = 0.0;
/// ```
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct PssmStats {
    pub(crate) lambda_ungapped: Option<f64>,
    pub(crate) kappa_ungapped: Option<f64>,
    pub(crate) h_ungapped: Option<f64>,
    pub(crate) lambda: Option<f64>,
    pub(crate) kappa: Option<f64>,
    pub(crate) h: Option<f64>,
}

/// A query given as a PSSM (psiblast `-in_pssm`): its rows, the statistical values of its file,
/// the query record made from its Bioseq, and NCBI's validation error of the PSSM, which NCBI
/// raises when the search starts (after the outfmt 0 prolog).
#[derive(Debug, Clone)]
pub struct PsiSearch {
    pub(crate) pssm: rows::QueryPssm,
    /// The PSSM's frequency ratios, 28 per position (`psi_matrix->freq_ratios`), which the
    /// position-based composition-based statistics read; empty when the PSSM has none.
    pub(crate) freq_ratios: Vec<[f64; crate::core::blast_util::BLASTAA_SIZE]>,
    pub(crate) stats: PssmStats,
    pub(crate) query: FastaRecord,
    pub(crate) validation_error: Option<String>,
    /// The warnings of `PsiBlastSetupScoreBlock` for the search's options (the frequency ratios
    /// all zero, the composition-based statistics reset to mode 1), which the query's warning
    /// line carries (`algorithm::psiblast::in_pssm::setup_messages`).
    pub(crate) setup_messages: Vec<String>,
    /// The title of the query's warning label when it is not the record's: a PSSM query's
    /// (`CBlastQuerySourceBioseqSet::GetTitle`, `algorithm::psiblast::in_pssm`); `None` for a FASTA
    /// query's later rounds, whose label keeps the record's title.
    pub(crate) warning_title: Option<Vec<u8>>,
}

impl PsiSearch {
    /// `kbp_psi[0]`: the file's ungapped values, or those of `kbp_std` (the query's composition
    /// with BLOSUM62, which NCBI also computes for `kbp_psi` before it overwrites them).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:89-107
    /// ```c++
    ///     // Assign the ungapped Karlin-Altschul block
    ///     if (pssm->GetPssm().GetLambdaUngapped() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_psi[0]->Lambda = pssm->GetPssm().GetLambdaUngapped();
    ///     } else if (score_blk->kbp_std[0]->Lambda > 0.0) {
    ///         score_blk->kbp_psi[0]->Lambda = score_blk->kbp_std[0]->Lambda;
    ///     }
    ///
    ///     if (pssm->GetPssm().GetKappaUngapped() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_psi[0]->K = pssm->GetPssm().GetKappaUngapped();
    ///     } else if (score_blk->kbp_std[0]->K > 0.0) {
    ///         score_blk->kbp_psi[0]->K = score_blk->kbp_std[0]->K;
    ///     }
    ///     score_blk->kbp_psi[0]->logK = log(score_blk->kbp_psi[0]->K);
    ///
    ///     if (pssm->GetPssm().GetHUngapped() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_psi[0]->H = pssm->GetPssm().GetHUngapped();
    ///     } else if (score_blk->kbp_std[0]->K > 0.0) {
    ///         score_blk->kbp_psi[0]->H = score_blk->kbp_std[0]->H;
    ///     }
    /// ```
    /// `kbp_psi[0]` holds `kbp_std[0]`'s values before (`Blast_ScoreBlkKbpUngappedCalc`,
    /// blast_stat.c:2799-2808: the same calculation), so either branch leaves those values.
    pub(crate) fn ungapped_karlin(&self, kbp_std: &KarlinParams) -> KarlinParams {
        KarlinParams {
            lambda: self.stats.lambda_ungapped.unwrap_or(kbp_std.lambda),
            k: self.stats.kappa_ungapped.unwrap_or(kbp_std.k),
            h: self.stats.h_ungapped.unwrap_or(kbp_std.h),
            ..*kbp_std
        }
    }

    /// `kbp_gap_psi[0]`: the file's gapped values, or those of `kbp_gap_std` (the table's values
    /// of BLOSUM62 with the search's gap costs, which `kbp_gap_psi` is a copy of).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:109-127
    /// ```c++
    ///     // Assign the gapped Karlin-Altschul block
    ///     if (pssm->GetPssm().GetLambda() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_gap_psi[0]->Lambda = pssm->GetPssm().GetLambda();
    ///     } else if (score_blk->kbp_gap_std[0]->Lambda > 0.0) {
    ///         score_blk->kbp_gap_psi[0]->Lambda = score_blk->kbp_gap_std[0]->Lambda;
    ///     }
    ///
    ///     if (pssm->GetPssm().GetKappa() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_gap_psi[0]->K = pssm->GetPssm().GetKappa();
    ///     } else if (score_blk->kbp_gap_std[0]->K > 0.0) {
    ///         score_blk->kbp_gap_psi[0]->K = score_blk->kbp_gap_std[0]->K;
    ///     }
    ///     score_blk->kbp_gap_psi[0]->logK = log(score_blk->kbp_gap_psi[0]->K);
    ///
    ///     if (pssm->GetPssm().GetH() != CPssm::kInvalidStat) {
    ///         score_blk->kbp_gap_psi[0]->H = pssm->GetPssm().GetH();
    ///     } else if (score_blk->kbp_gap_std[0]->H > 0.0) {
    ///         score_blk->kbp_gap_psi[0]->H = score_blk->kbp_gap_std[0]->H;
    ///     }
    /// ```
    pub(crate) fn gapped_karlin(&self, kbp_gap_std: &KarlinParams) -> KarlinParams {
        KarlinParams {
            lambda: self.stats.lambda.unwrap_or(kbp_gap_std.lambda),
            k: self.stats.kappa.unwrap_or(kbp_gap_std.k),
            h: self.stats.h.unwrap_or(kbp_gap_std.h),
            ..*kbp_gap_std
        }
    }

    /// The query record made from the PSSM's Bioseq.
    pub fn query(&self) -> &FastaRecord {
        &self.query
    }

    /// The PSSM search of a later round of the FASTA query `query` (psiblast's rounds,
    /// `algorithm::psiblast::rounds`): the engine's PSSM of the round before, in full precision
    /// (NCBI's `CPssmWithParameters` in memory, never read back from text): its scores, its
    /// frequency ratios and the six statistical values of `x_PSIMatrix2Asn1`, each `None` where
    /// it is 0.0 (`CPssm::kInvalidStat`), as `psi_search_of` reads a file's; and the warnings of
    /// `PsiBlastSetupScoreBlock` for the round's options (`setup_messages`). The query stays the
    /// FASTA query: `CPsiBlastImpl::SetPssm` replaces only the PSSM.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:173-181
    /// ```c++
    /// CPsiBlastImpl::SetPssm(CConstRef<objects::CPssmWithParameters> pssm)
    /// {
    ///     if (pssm.Empty()) {
    ///         NCBI_THROW(CBlastException, eInvalidArgument,
    ///                    "Setting empty reference for pssm");
    ///     }
    ///     CPsiBlastValidate::Pssm(*pssm, true);
    ///     m_Pssm.Reset(const_cast<CPssmWithParameters*>(&*pssm));
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:648-660 (quoted at
    /// `algorithm::psiblast::pssm::pssm_with_parameters`: the scores and values that the next
    /// round's `PsiBlastSetupScoreBlock` reads)
    pub(crate) fn from_engine(
        pssm: &crate::core::psi::PsiPssm,
        query: FastaRecord,
        setup_messages: Vec<String>,
    ) -> Self {
        let stat = |value: f64| (value != 0.0).then_some(value);
        Self {
            pssm: rows::QueryPssm {
                rows: pssm.scores.clone(),
            },
            freq_ratios: pssm
                .diagnostics
                .frequency_ratios
                .clone()
                .unwrap_or_default(),
            stats: PssmStats {
                lambda_ungapped: stat(pssm.ungapped.lambda),
                kappa_ungapped: stat(pssm.ungapped.k),
                h_ungapped: stat(pssm.ungapped.h),
                lambda: stat(pssm.gapped.lambda),
                kappa: stat(pssm.gapped.k),
                h: stat(pssm.gapped.h),
            },
            query,
            validation_error: None,
            setup_messages,
            warning_title: None,
        }
    }
}
