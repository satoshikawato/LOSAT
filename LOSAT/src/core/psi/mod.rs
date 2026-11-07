//! PSI-BLAST's PSSM engine: the port of NCBI's `algo/blast/core/blast_psi.c`
//! (`PSICreatePssmWithDiagnostics`) and of the stages in `blast_psi_priv.c` (`stages.rs`).
//!
//! The engine builds a position-specific scoring matrix from a multiple sequence alignment whose
//! first row is the query. It is reached only where a PSSM is requested (psiblast's
//! `-save_pssm_after_last_round` with `-out_pssm` or `-out_ascii_pssm`, `algorithm::psiblast::pssm`);
//! no search path calls it. It returns what `s_PSISavePssm` and `_PSISaveDiagnostics` keep: the
//! scores, the Karlin-Altschul blocks of the scaling, and the diagnostics of psiblast's request
//! (`PsiDiagnosticsRequest`): the frequency ratios always, and with `-out_ascii_pssm` the
//! information content, the weighted residue frequencies, the relative weight of gapless real
//! matches to pseudocounts, sigma, the interval sizes and the number of matching sequences of each
//! position.

mod stages;

use crate::utils::matrix::BLASTAA_SIZE;

/// One cell of the multiple sequence alignment (`PSIMsaCell`).
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_psi.h:49-53
// ```c
// typedef struct PSIMsaCell {
//     Uint1   letter;             /**< Preferred letter at this position, in
//                                   ncbistdaa encoding */
//     Boolean is_aligned;         /**< Is this letter part of the alignment? */
// } PSIMsaCell;
// ```
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct PsiMsaCell {
    pub(crate) letter: u8,
    pub(crate) is_aligned: bool,
}

/// The engine's input (`PSIMsa`): `rows[0]` is the query, `rows[1..]` the aligned sequences;
/// every row has `query_length` cells.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:462-485
// ```c
// PSIMsaNew(const PSIMsaDimensions* dimensions)
// {
//     ...
//     retval->data = (PSIMsaCell**) _PSIAllocateMatrix(dimensions->num_seqs + 1,
//                                                     dimensions->query_length,
//                                                     sizeof(PSIMsaCell));
// ```
#[derive(Debug, Clone)]
pub(crate) struct PsiMsa {
    pub(crate) query_length: usize,
    pub(crate) rows: Vec<Vec<PsiMsaCell>>,
}

impl PsiMsa {
    /// `PSIMsaNew`: `num_seqs + 1` rows of empty cells (letter 0, not aligned).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:493-498
    // ```c
    //         for (s = 0; s < dimensions->num_seqs + 1; s++) {
    //             for (p = 0; p < dimensions->query_length; p++) {
    //                 retval->data[s][p].letter = 0;
    //                 retval->data[s][p].is_aligned = FALSE;
    //             }
    //         }
    // ```
    pub(crate) fn new(query_length: usize, num_seqs: usize) -> Self {
        Self {
            query_length,
            rows: vec![vec![PsiMsaCell::default(); query_length]; num_seqs + 1],
        }
    }
}

/// A Karlin-Altschul block (`Blast_KarlinBlk`): the Lambda, K and H that the PSSM keeps.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:66-72
// ```c
// typedef struct Blast_KarlinBlk {
//       double   Lambda; /**< Lambda value used in statistics */
//       double   K; /**< K value used in statistics */
//       double   logK; /**< natural log of K value used in statistics */
//       double   H; /**< H value used in statistics */
//       double   paramC;  /**< for use in seed. */
//    } Blast_KarlinBlk;
// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct KarlinBlk {
    pub(crate) lambda: f64,
    pub(crate) k: f64,
    pub(crate) h: f64,
}

/// The score block the engine reads (`BlastScoreBlk` as `CPssmEngine::x_InitializeScoreBlock`
/// sets it up for psiblast and `SetUngappedStatisticalParams` updates it): the matrix
/// (`sbp->matrix->data`), the ideal block (`sbp->kbp_ideal`), the standard gapped K
/// (`sbp->kbp_gap_std[0]->K`), the gapped PSI block before the scaling (`sbp->kbp_gap_psi[0]`),
/// and the matrix's frequency ratios and background frequencies (`_PSIMatrixFrequencyRatiosNew`,
/// `Blast_GetMatrixBackgroundFreq`).
pub(crate) struct PsiScoreBlk {
    /// `sbp->matrix->data[i][j]`, NCBIstdaa by NCBIstdaa.
    pub(crate) matrix: [[i32; BLASTAA_SIZE]; BLASTAA_SIZE],
    pub(crate) ideal_lambda: f64,
    /// `sbp->kbp_ideal->K`.
    pub(crate) ideal_k: f64,
    /// `sbp->kbp_gap_std[0]->K`: the gapped K of the matrix and gap costs
    /// (`Blast_KarlinBlkGappedCalc`).
    pub(crate) gap_std_k: f64,
    /// `sbp->kbp_gap_psi[0]` before the scaling: a copy of `kbp_gap_std[0]`
    /// (`Blast_ScoreBlkKbpGappedCalc`, blast_setup.c:114-120) that
    /// `SetUngappedStatisticalParams` replaces by the round's PSI gapped block
    /// (`set_ungapped_statistical_params`). The scaling changes only its K.
    pub(crate) gap_psi: KarlinBlk,
    /// `SFreqRatios::data`, NCBIstdaa by NCBIstdaa.
    pub(crate) freq_ratios: &'static [[f64; BLASTAA_SIZE]; BLASTAA_SIZE],
    /// `SFreqRatios::bit_scale_factor`.
    pub(crate) bit_scale_factor: i32,
    /// The background frequencies in the order A R N D C Q E G H I L K M F P S T W Y V.
    pub(crate) background: &'static [f64; 20],
}

impl PsiScoreBlk {
    /// `CPssmEngine::SetUngappedStatisticalParams` with the round's ancillary data, whose PSI
    /// gapped block is `psi_gapped`. Its PSI ungapped block replaces `kbp_psi[0]`, which the
    /// scaling computes again before anything reads it (`_PSIUpdateLambdaK` at every factor of
    /// `_PSIScaleMatrix`), so it is not kept.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:483-509
    // ```c++
    //     if (ancillary_data->GetPsiUngappedKarlinBlk()) {
    //         _ASSERT(m_ScoreBlk->kbp_psi && m_ScoreBlk->kbp_psi[0]);
    //         m_ScoreBlk->kbp_psi[0]->Lambda =
    //             ancillary_data->GetPsiUngappedKarlinBlk()->Lambda;
    //         ...
    //     if (ancillary_data->GetPsiGappedKarlinBlk()) {
    //         _ASSERT(m_ScoreBlk->kbp_gap_psi && m_ScoreBlk->kbp_gap_psi[0]);
    //         m_ScoreBlk->kbp_gap_psi[0]->Lambda =
    //             ancillary_data->GetPsiGappedKarlinBlk()->Lambda;
    //         m_ScoreBlk->kbp_gap_psi[0]->K =
    //             ancillary_data->GetPsiGappedKarlinBlk()->K;
    //         m_ScoreBlk->kbp_gap_psi[0]->logK = log(m_ScoreBlk->kbp_gap_psi[0]->K);
    //         m_ScoreBlk->kbp_gap_psi[0]->H =
    //             ancillary_data->GetPsiGappedKarlinBlk()->H;
    //     }
    // ```
    pub(crate) fn set_ungapped_statistical_params(&mut self, psi_gapped: KarlinBlk) {
        self.gap_psi = psi_gapped;
    }
}

/// The diagnostics the engine is asked for (`PSIDiagnosticsRequest`), as psiblast asks for them.
/// psiblast never asks for the residue frequencies or the independent observations, which are
/// not modelled.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:591-608
// ```c
// PSIDiagnosticsRequestNewEx(Boolean save_ascii_pssm)
// {
//     PSIDiagnosticsRequest* retval = PSIDiagnosticsRequestNew();
//     if ( !retval ) {
//         return NULL;
//     }
//
//     retval->frequency_ratios = TRUE;
//     if (save_ascii_pssm) {
//         retval->information_content = TRUE;
//         retval->weighted_residue_frequencies = TRUE;
//         retval->gapless_column_weights = TRUE;
//         retval->sigma = TRUE;
//         retval->interval_sizes = TRUE;
//         retval->num_matching_seqs = TRUE;
//     }
//     return retval;
// }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:175-176
// ```c++
//     CPSIDiagnosticsRequest
//         diags(PSIDiagnosticsRequestNewEx(m_CmdLineArgs->SaveAsciiPssm()));
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PsiDiagnosticsRequest {
    pub(crate) information_content: bool,
    pub(crate) weighted_residue_frequencies: bool,
    pub(crate) frequency_ratios: bool,
    pub(crate) gapless_column_weights: bool,
    pub(crate) sigma: bool,
    pub(crate) interval_sizes: bool,
    pub(crate) num_matching_seqs: bool,
}

impl PsiDiagnosticsRequest {
    /// `PSIDiagnosticsRequestNewEx(save_ascii_pssm)`.
    pub(crate) fn new_ex(save_ascii_pssm: bool) -> Self {
        Self {
            information_content: save_ascii_pssm,
            weighted_residue_frequencies: save_ascii_pssm,
            frequency_ratios: true,
            gapless_column_weights: save_ascii_pssm,
            sigma: save_ascii_pssm,
            interval_sizes: save_ascii_pssm,
            num_matching_seqs: save_ascii_pssm,
        }
    }
}

/// The diagnostics that were asked for (`PSIDiagnosticsResponse`), by query position; `None`
/// where the request did not ask.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct PsiDiagnostics {
    /// `information_content`.
    pub(crate) information_content: Option<Vec<f64>>,
    /// `weighted_residue_freqs[position][residue]`.
    pub(crate) weighted_residue_freqs: Option<Vec<[f64; BLASTAA_SIZE]>>,
    /// `frequency_ratios[position][residue]`.
    pub(crate) frequency_ratios: Option<Vec<[f64; BLASTAA_SIZE]>>,
    /// `gapless_column_weights`.
    pub(crate) gapless_column_weights: Option<Vec<f64>>,
    /// `sigma`.
    pub(crate) sigma: Option<Vec<f64>>,
    /// `interval_sizes`.
    pub(crate) interval_sizes: Option<Vec<u32>>,
    /// `num_matching_seqs`.
    pub(crate) num_matching_seqs: Option<Vec<u32>>,
}

/// The PSSM (`PSIMatrix`) and its diagnostics.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_psi.h:150-160
// ```c
// typedef struct PSIMatrix {
//     Uint4   ncols;      /**< Number of columns in PSSM (query_length) */
//     Uint4   nrows;      /**< Number of rows in PSSM (alphabet_size) */
//     int**   pssm;       /**< Position-specific score matrix */
//     double  lambda;     /**< Lambda Karlin-Altschul parameter */
//     double  kappa;      /**< Kappa Karlin-Altschul parameter */
//     double  h;          /**< H Karlin-Altschul parameter */
//     double  ung_lambda; /**< Ungapped Lambda Karlin-Altschul parameter */
//     double  ung_kappa;  /**< Ungapped Kappa Karlin-Altschul parameter */
//     double  ung_h;      /**< Ungapped H Karlin-Altschul parameter */
// } PSIMatrix;
// ```
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PsiPssm {
    /// `pssm[position][residue]`.
    pub(crate) scores: Vec<[i32; BLASTAA_SIZE]>,
    /// `lambda`, `kappa`, `h`: `kbp_gap_psi[0]` after the scaling.
    pub(crate) gapped: KarlinBlk,
    /// `ung_lambda`, `ung_kappa`, `ung_h`: `kbp_psi[0]` after the scaling.
    pub(crate) ungapped: KarlinBlk,
    pub(crate) diagnostics: PsiDiagnostics,
}

/// The engine's failures that its stages can return for psiblast's input (`PSIERR_*`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.h:366-392
// ```c
// #define PSI_SUCCESS             (0)
// #define PSIERR_BADPARAM         (-1)
// #define PSIERR_OUTOFMEM         (-2)
// #define PSIERR_BADSEQWEIGHTS    (-3)
// #define PSIERR_NOFREQRATIOS     (-4)
// #define PSIERR_POSITIVEAVGSCORE (-5)
// #define PSIERR_NOALIGNEDSEQS    (-6)
// #define PSIERR_GAPINQUERY       (-7)
// #define PSIERR_UNALIGNEDCOLUMN  (-8)
// #define PSIERR_COLUMNOFGAPS     (-9)
// #define PSIERR_STARTINGGAP      (-10)
// #define PSIERR_ENDINGGAP        (-11)
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PsiError {
    BadSeqWeights,
    PositiveAvgScore,
    NoAlignedSeqs,
    GapInQuery,
    UnalignedColumn,
    ColumnOfGaps,
    StartingGap,
    EndingGap,
}

/// `PSICreatePssmWithDiagnostics` with psiblast's options (`PSIBlastOptionsNew`: no structure
/// group compatibility mode, no IMPALA scaling, unaligned positions not ignored) and the
/// pseudocount of `-pseudocount` (already validated, `s_Validate`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:104-221
// ```c
//     packed_msa = _PSIPackedMsaNew(msap);
//     status = _PSIPurgeBiasedSegments(packed_msa);
//     ...
//     msa = _PSIMsaNew(packed_msa, (Uint4) sbp->alphabet_size);
//     aligned_block = _PSIAlignedBlockNew(msa->dimensions->query_length);
//     seq_weights = _PSISequenceWeightsNew(msa->dimensions, sbp);
//     internal_pssm = _PSIInternalPssmDataNew(msa->dimensions->query_length,
//                                             (Uint4) sbp->alphabet_size);
//     ...
//         status = _PSIValidateMSA(msa, options->ignore_unaligned_positions);
//     ...
//     status = _PSIComputeAlignmentBlocks(msa, aligned_block);
//     ...
//     status = _PSIComputeSequenceWeights(msa, aligned_block,
//                                         options->nsg_compatibility_mode,
//                                         seq_weights);
//     ...
//     status = _PSIComputeFreqRatios(msa, seq_weights, sbp, aligned_block,
//                                    options->pseudo_count,
//                                    options->nsg_compatibility_mode,
//                                    internal_pssm);
//     ...
//     status = _PSICreateAndScalePssmFromFrequencyRatios
//         (internal_pssm, msa->query, msa->dimensions->query_length,
//          seq_weights->std_prob, sbp, options->impala_scaling_factor);
//     ...
//     s_PSISavePssm(internal_pssm, sbp, *pssm);
//     ...
//         status = _PSISaveDiagnostics(msa, aligned_block, seq_weights,
//                                      internal_pssm, *diagnostics);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1566-1572
// ```c
//    options->inclusion_ethresh = PSI_INCLUSION_ETHRESH;
//    options->pseudo_count = PSI_PSEUDO_COUNT_CONST;
//    options->use_best_alignment = TRUE;
//
//    options->nsg_compatibility_mode = FALSE;
//    options->impala_scaling_factor = kPSSM_NoImpalaScaling;
//    options->ignore_unaligned_positions = FALSE;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:440-457
// ```c
// s_PSISavePssm(const _PSIInternalPssmData* internal_pssm,
//              const BlastScoreBlk* sbp,
//              PSIMatrix* pssm)
// {
//     ...
//     _PSICopyMatrix_int(pssm->pssm, internal_pssm->pssm,
//                        pssm->ncols, pssm->nrows);
//
//     pssm->lambda = sbp->kbp_gap_psi[0]->Lambda;
//     pssm->kappa = sbp->kbp_gap_psi[0]->K;
//     pssm->h = sbp->kbp_gap_psi[0]->H;
//     pssm->ung_lambda = sbp->kbp_psi[0]->Lambda;
//     pssm->ung_kappa = sbp->kbp_psi[0]->K;
//     pssm->ung_h = sbp->kbp_psi[0]->H;
// }
// ```
// `kbp_psi[0]` is the block of the scaling's last `_PSIUpdateLambdaK` (`stages::scale_matrix`),
// and `kbp_gap_psi[0]` the score block's (`PsiScoreBlk::gap_psi`) with the K that the same call
// sets (`stages::gapped_psi_k`).
pub(crate) fn create_pssm_with_diagnostics(
    msa: &PsiMsa,
    pseudo_count: i32,
    sbp: &PsiScoreBlk,
    request: PsiDiagnosticsRequest,
) -> Result<PsiPssm, PsiError> {
    let mut packed = stages::PackedMsa::new(msa);
    stages::purge_biased_segments(&mut packed);
    let mut msa = stages::Msa::new(&packed);
    let query_length = msa.query_length;
    let mut aligned_block = stages::AlignedBlock::new(query_length);
    let mut seq_weights = stages::SequenceWeights::new(msa.num_seqs, query_length);
    let mut internal = stages::InternalPssm::new(query_length);
    stages::validate_msa(&msa)?;
    stages::compute_alignment_blocks(&mut msa, &mut aligned_block);
    stages::compute_sequence_weights(&msa, &aligned_block, &mut seq_weights)?;
    stages::compute_freq_ratios(
        &msa,
        &mut seq_weights,
        sbp,
        &aligned_block,
        pseudo_count,
        &mut internal,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:385-414
    // ```c
    //     status = _PSIConvertFreqRatiosToPSSM(internal_pssm, query, sbp, std_prob);
    //     ...
    //     if (impala_scaling_factor == kPSSM_NoImpalaScaling) {
    //         status = _PSIScaleMatrix(query, std_prob, internal_pssm, sbp);
    // ```
    stages::convert_freq_ratios_to_pssm(&mut internal, &msa.query, sbp, &seq_weights.std_prob);
    let ungapped = stages::scale_matrix(&msa.query, &seq_weights.std_prob, &mut internal, sbp)?;
    let gapped = KarlinBlk {
        k: stages::gapped_psi_k(ungapped, sbp),
        ..sbp.gap_psi
    };
    let diagnostics =
        stages::save_diagnostics(&msa, &aligned_block, &seq_weights, &internal, request);
    Ok(PsiPssm {
        scores: internal.into_scores(),
        gapped,
        ungapped,
        diagnostics,
    })
}

/// The scaled PSSM of the position-based composition-based statistics (`s_ScalePosMatrix`,
/// `stages::scale_pos_matrix`): the redo's `startMatrix`, and `kbp_psi[0]` and the gapped K of
/// the factor-1 PSSM, which `_PSIUpdateLambdaK` leaves in the redo's copy of the score block
/// (`kbp_gap_psi[0]->K`, `stages::gapped_psi_k`).
pub(crate) struct PosScaledMatrix {
    pub(crate) start_matrix: Vec<[i32; BLASTAA_SIZE]>,
    pub(crate) kbp_psi: KarlinBlk,
    pub(crate) gapped_k: f64,
}

/// `s_ScalePosMatrix` for `query` (NCBIstdaa), the PSSM's frequency ratios `freq_ratios` (28
/// rows each) and the score block `sbp`, scaled by `scaling_factor` (`SCALING_FACTOR`, 32).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2227
// ```c
//     if (self->positionBased) {
//         status = s_GetPosBasedStartFreqRatios(self->startFreqRatios,
//                                               queryBlk->length,
//                                               queryBlk->sequence,
//                                               matrixName,
//                                               sbp->psi_matrix->freq_ratios);
//         if (status == 0) {
//             status = s_ScalePosMatrix(self->startMatrix, matrixName,
//                                       sbp->psi_matrix->freq_ratios,
//                                       queryBlk->sequence,
//                                       queryBlk->length, sbp, scale_factor);
//             self->ungappedLambda = sbp->kbp_psi[0]->Lambda / scale_factor;
//         }
// ```
pub(crate) fn scale_pos_matrix(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
    scaling_factor: f64,
) -> PosScaledMatrix {
    let (start_matrix, kbp_psi) =
        stages::scale_pos_matrix(query, freq_ratios, sbp, scaling_factor);
    PosScaledMatrix {
        start_matrix,
        kbp_psi,
        gapped_k: stages::gapped_psi_k(kbp_psi, sbp),
    }
}

/// The PSSM of a query from its frequency ratios (`PSICreatePssmFromFrequencyRatios`,
/// `stages::pssm_from_freq_ratios`): the scores and `kbp_gap_psi[0]` after the scaling (the
/// lambda, kappa and h that `s_PSISavePssm` keeps). `PsiBlastComputePssmScores` gives a PSSM of
/// frequency ratios only these scores and values.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:486-498
// ```c++
//     pssm->SetPssm().SetFinalData().SetScores() =
//         pssm_with_scores->GetPssm().GetFinalData().GetScores();
//     pssm->SetPssm().SetFinalData().SetLambda() =
//         pssm_with_scores->GetPssm().GetFinalData().GetLambda();
//     pssm->SetPssm().SetFinalData().SetKappa() =
//         pssm_with_scores->GetPssm().GetFinalData().GetKappa();
//     pssm->SetPssm().SetFinalData().SetH() =
//         pssm_with_scores->GetPssm().GetFinalData().GetH();
// ```
pub(crate) fn pssm_from_freq_ratios(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
) -> Result<(Vec<[i32; BLASTAA_SIZE]>, KarlinBlk), PsiError> {
    let (scores, ungapped) = stages::pssm_from_freq_ratios(query, freq_ratios, sbp)?;
    Ok((
        scores,
        KarlinBlk {
            k: stages::gapped_psi_k(ungapped, sbp),
            ..sbp.gap_psi
        },
    ))
}
