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
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct PsiMsaCell {
    pub(crate) letter: u8,
    pub(crate) is_aligned: bool,
}

/// The engine's input (`PSIMsa`): `rows[0]` is the query, `rows[1..]` the aligned sequences;
/// every row has `query_length` cells.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:462-485
#[derive(Debug, Clone)]
pub(crate) struct PsiMsa {
    pub(crate) query_length: usize,
    pub(crate) rows: Vec<Vec<PsiMsaCell>>,
}

impl PsiMsa {
    /// `PSIMsaNew`: `num_seqs + 1` rows of empty cells (letter 0, not aligned).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:493-498
    pub(crate) fn new(query_length: usize, num_seqs: usize) -> Self {
        Self {
            query_length,
            rows: vec![vec![PsiMsaCell::default(); query_length]; num_seqs + 1],
        }
    }
}

/// A Karlin-Altschul block (`Blast_KarlinBlk`): the Lambda, K and H that the PSSM keeps.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:66-72
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
    pub(crate) fn set_ungapped_statistical_params(&mut self, psi_gapped: KarlinBlk) {
        self.gap_psi = psi_gapped;
    }
}

/// The diagnostics the engine is asked for (`PSIDiagnosticsRequest`), as psiblast asks for them.
/// psiblast never asks for the residue frequencies or the independent observations, which are
/// not modelled.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:591-608
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:175-176
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
    /// `PSIDiagnosticsRequestNewEx` with `save_ascii_pssm`.
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
    /// `PSIERR_UNKNOWN` (-255): a column without weight in the structure group's compatibility
    /// mode (`_PSIComputeFreqRatios`).
    Unknown,
    /// The stage of the frequency ratios of a matrix without frequency data (IDENTITY), which
    /// NCBI's `_PSIComputeFreqRatios` reaches with a null pointer and crashes (SIGSEGV); the
    /// engine stops here and the caller refuses the run (`algorithm::psiblast::in_msa`).
    NoMatrixFrequencyData,
    /// The scaling's first `_PSIUpdateLambdaK` with a null `kbp_psi[0]` (the query's ungapped
    /// Karlin-Altschul block could not be computed, [`PsiOptions::ungapped_block_missing`]),
    /// where NCBI crashes (SIGSEGV) before any other error of the scaling.
    NoUngappedKarlinBlock,
}

/// The `PSIBlastOptions` that the engine's stages read besides the pseudocount: the options of the
/// rounds' PSSMs have neither (`PSIBlastOptionsNew`), those of `-in_msa` (`CPsiBlastInputClustalW`,
/// `x_CreatePssmFromMsa`) set `ignore_unaligned_positions` and, with `-ignore_msa_master`,
/// `nsg_compatibility_mode`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct PsiOptions {
    pub(crate) ignore_unaligned_positions: bool,
    pub(crate) nsg_compatibility_mode: bool,
    /// The score block's `kbp_psi[0]` is null (the query's ungapped Karlin-Altschul parameters
    /// could not be computed: `Blast_ScoreBlkKbpUngappedCalc`'s failure goes unreported,
    /// blast_setup.c:353-370); the scaling then ends with [`PsiError::NoUngappedKarlinBlock`].
    pub(crate) ungapped_block_missing: bool,
}

/// `PSICreatePssmWithDiagnostics` with psiblast's options (`PSIBlastOptionsNew`: no structure
/// group compatibility mode, no IMPALA scaling, unaligned positions not ignored) and the
/// pseudocount of `-pseudocount` (already validated, `s_Validate`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:104-221
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1566-1572
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:440-457
// `kbp_psi[0]` is the block of the scaling's last `_PSIUpdateLambdaK` (`stages::scale_matrix`),
// and `kbp_gap_psi[0]` the score block's (`PsiScoreBlk::gap_psi`) with the K that the same call
// sets (`stages::gapped_psi_k`).
pub(crate) fn create_pssm_with_diagnostics(
    msa: &PsiMsa,
    pseudo_count: i32,
    sbp: &PsiScoreBlk,
    request: PsiDiagnosticsRequest,
) -> Result<PsiPssm, PsiError> {
    create_pssm_with_options(msa, pseudo_count, PsiOptions::default(), Some(sbp), request)
}

/// `PSICreatePssmWithDiagnostics` with the options `options`; `sbp` is `None` for a matrix
/// without frequency data (the engine then stops at the stage of the frequency ratios with
/// [`PsiError::NoMatrixFrequencyData`], where NCBI crashes).
pub(crate) fn create_pssm_with_options(
    msa: &PsiMsa,
    pseudo_count: i32,
    options: PsiOptions,
    sbp: Option<&PsiScoreBlk>,
    request: PsiDiagnosticsRequest,
) -> Result<PsiPssm, PsiError> {
    let mut packed = stages::PackedMsa::new(msa);
    stages::purge_biased_segments(&mut packed);
    let mut msa = stages::Msa::new(&packed);
    let query_length = msa.query_length;
    let mut aligned_block = stages::AlignedBlock::new(query_length);
    let mut seq_weights = stages::SequenceWeights::new(msa.num_seqs, query_length);
    let mut internal = stages::InternalPssm::new(query_length);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:151-156
    if options.nsg_compatibility_mode {
        msa.structure_group_customization();
        stages::validate_msa_structure_group(&msa)?;
    } else {
        stages::validate_msa(&msa, options.ignore_unaligned_positions)?;
    }
    stages::compute_alignment_blocks(&mut msa, &mut aligned_block);
    stages::compute_sequence_weights(
        &msa,
        &aligned_block,
        options.nsg_compatibility_mode,
        &mut seq_weights,
    )?;
    let Some(sbp) = sbp else {
        return Err(PsiError::NoMatrixFrequencyData);
    };
    stages::compute_freq_ratios(
        &msa,
        &mut seq_weights,
        sbp,
        &aligned_block,
        pseudo_count,
        options.nsg_compatibility_mode,
        &mut internal,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:385-414
    stages::convert_freq_ratios_to_pssm(&mut internal, &msa.query, sbp, &seq_weights.std_prob);
    if options.ungapped_block_missing {
        return Err(PsiError::NoUngappedKarlinBlock);
    }
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
/// `stages::scale_pos_matrix`): the redo's `startMatrix`, and `kbp_psi[0]` and the gapped K that
/// `_PSIUpdateLambdaK` leaves in the redo's copy of the score block (`kbp_gap_psi[0]->K`,
/// `stages::gapped_psi_k`): those of the factor-1 PSSM, or, where NCBI's `Blast_ScoreFreqNew`
/// gives no score frequencies for it, the copy's `kbp_psi[0]` as it was and its K.
pub(crate) struct PosScaledMatrix {
    pub(crate) start_matrix: Vec<[i32; BLASTAA_SIZE]>,
    pub(crate) kbp_psi: KarlinBlk,
    pub(crate) gapped_k: f64,
}

/// `s_ScalePosMatrix` for `query` (NCBIstdaa), the PSSM's frequency ratios `freq_ratios` (28
/// rows each) and the score block `sbp`, scaled by `scaling_factor` (`SCALING_FACTOR`, 32);
/// `kbp_psi` is the search's `kbp_psi[0]` (`PsiBlastSetupScoreBlock`'s), which the redo's copy
/// holds before the call.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2215-2227
pub(crate) fn scale_pos_matrix(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
    scaling_factor: f64,
    kbp_psi: KarlinBlk,
) -> PosScaledMatrix {
    let (start_matrix, factor_one) =
        stages::scale_pos_matrix(query, freq_ratios, sbp, scaling_factor);
    let kbp_psi = factor_one.unwrap_or(kbp_psi);
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
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:486-499
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
