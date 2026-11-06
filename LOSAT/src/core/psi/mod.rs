//! PSI-BLAST's PSSM engine: the port of NCBI's `algo/blast/core/blast_psi.c`
//! (`PSICreatePssmWithDiagnostics`) and of the stages in `blast_psi_priv.c` (`stages.rs`).
//!
//! The engine builds a position-specific scoring matrix from a multiple sequence alignment whose
//! first row is the query. It is reached only where a PSSM is requested (psiblast's
//! `-save_pssm_after_last_round` with `-out_ascii_pssm`, `algorithm::psiblast::pssm`); no search
//! path calls it. It keeps only the diagnostics that the ASCII PSSM reads
//! (`PrintAsciiPssm`, blastfmtutil.cpp:172-327): the information content, the weighted residue
//! frequencies, the relative weight of gapless real matches to pseudocounts and the number of
//! matching sequences of each position.

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

/// The score block the engine reads (`BlastScoreBlk` as `CPssmEngine::x_InitializeScoreBlock`
/// sets it up for psiblast): the matrix (`sbp->matrix->data`), the ideal Lambda
/// (`sbp->kbp_ideal->Lambda`), and the matrix's frequency ratios and background frequencies
/// (`_PSIMatrixFrequencyRatiosNew`, `Blast_GetMatrixBackgroundFreq`).
pub(crate) struct PsiScoreBlk {
    /// `sbp->matrix->data[i][j]`, NCBIstdaa by NCBIstdaa.
    pub(crate) matrix: [[i32; BLASTAA_SIZE]; BLASTAA_SIZE],
    pub(crate) ideal_lambda: f64,
    /// `SFreqRatios::data`, NCBIstdaa by NCBIstdaa.
    pub(crate) freq_ratios: &'static [[f64; BLASTAA_SIZE]; BLASTAA_SIZE],
    /// `SFreqRatios::bit_scale_factor`.
    pub(crate) bit_scale_factor: i32,
    /// The background frequencies in the order A R N D C Q E G H I L K M F P S T W Y V.
    pub(crate) background: &'static [f64; 20],
}

/// The PSSM and the diagnostics the ASCII PSSM reads, by query position.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PsiPssm {
    /// `PSIMatrix::pssm[position][residue]`.
    pub(crate) scores: Vec<[i32; BLASTAA_SIZE]>,
    /// `PSIDiagnosticsResponse::information_content`.
    pub(crate) information_content: Vec<f64>,
    /// `PSIDiagnosticsResponse::weighted_residue_freqs[position][residue]`.
    pub(crate) weighted_residue_freqs: Vec<[f64; BLASTAA_SIZE]>,
    /// `PSIDiagnosticsResponse::gapless_column_weights`.
    pub(crate) gapless_column_weights: Vec<f64>,
    /// `PSIDiagnosticsResponse::num_matching_seqs`.
    pub(crate) num_matching_seqs: Vec<u32>,
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
// The Karlin-Altschul blocks that the scaling updates (`kbp_psi`, `kbp_gap_psi`) are not kept:
// the ASCII PSSM prints those of the round's ancillary data (`PrintAsciiPssm`).
pub(crate) fn create_pssm_with_diagnostics(
    msa: &PsiMsa,
    pseudo_count: i32,
    sbp: &PsiScoreBlk,
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
    stages::scale_matrix(&msa.query, &seq_weights.std_prob, &mut internal, sbp)?;
    Ok(stages::save_diagnostics(
        &msa,
        &aligned_block,
        &seq_weights,
        internal,
    ))
}
