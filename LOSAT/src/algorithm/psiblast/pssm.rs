//! psiblast's PSSM after the last round (`-num_iterations 1 -save_pssm_after_last_round
//! -out_ascii_pssm F`): the input of NCBI's PSSM engine from a query's round-1 alignments
//! (`CPsiBlastInputData`), the engine's score block (`CPssmEngine`), the ASCII PSSM
//! (`CBlastFormatUtil::PrintAsciiPssm`) and its file (`CAutoOutputFileReset`). The engine itself
//! is `core::psi`. The blastp engine reaches this module only when a PSSM is requested
//! (`blast_engine.rs` `x_psi_ascii_pssms`).

use std::io::Write as _;

use anyhow::Result;

use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::core::psi::{self, PsiError, PsiMsa, PsiMsaCell, PsiPssm, PsiScoreBlk};
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

/// psiblast's PSSM options as given on the command line.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PsiblastPssmOptions {
    pub save_pssm_after_last_round: bool,
    pub save_each_pssm: bool,
    pub out_ascii_pssm: Option<String>,
    pub pseudocount: i32,
}

/// A run that writes the ASCII PSSM of each query whose round-1 results have alignments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PssmRequest {
    /// `-out_ascii_pssm`, as given.
    pub ascii_path: String,
    /// `-save_each_pssm`: the files are `F.1`, `F.2`, ... (`CAutoOutputFileReset` with versions).
    pub save_each: bool,
    /// `-pseudocount` (`PSIBlastOptions::pseudo_count`).
    pub pseudocount: i32,
}

impl PsiblastPssmOptions {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1283-1303
    // ```c++
    //         if (args.Exist(kArgSaveLastPssm) && args[kArgSaveLastPssm] &&
    //             (!args.Exist(kArgPSIOutputChkPntFile) ||
    //              !args[kArgPSIOutputChkPntFile]) &&
    //             (!args.Exist(kArgAsciiPssmOutputFile) ||
    //              !args[kArgAsciiPssmOutputFile])) {
    //
    //             NCBI_THROW(CInputException, eInvalidInput, kArgSaveLastPssm +
    //                        " option requires " + kArgPSIOutputChkPntFile + " or " +
    //                        kArgAsciiPssmOutputFile);
    //         }
    //
    //         if (args.Exist(kArgSaveAllPssms) && args[kArgSaveAllPssms] &&
    //             (!args.Exist(kArgPSIOutputChkPntFile) ||
    //              !args[kArgPSIOutputChkPntFile]) &&
    //             (!args.Exist(kArgAsciiPssmOutputFile) ||
    //              !args[kArgAsciiPssmOutputFile])) {
    //
    //             NCBI_THROW(CInputException, eInvalidInput, kArgSaveAllPssms +
    //                        " option requires " + kArgPSIOutputChkPntFile + " or " +
    //                        kArgAsciiPssmOutputFile);
    //         }
    // ```
    // `CPsiBlastArgs::ExtractAlgorithmOptions` runs after the composition-based statistics
    // handler's (psiblast_args.cpp:102-110, the order of `m_Args`) and before `Validate`
    // (blast_args.cpp:3826-3832); `-out_pssm` is refused by the parser (`cli.rs`
    // `is_unported_psiblast_arg`), so only `-out_ascii_pssm` can satisfy it.
    pub(crate) fn check(&self) -> Result<()> {
        if self.save_pssm_after_last_round && self.out_ascii_pssm.is_none() {
            return Err(crate::blastinput::app::options_error(
                "save_pssm_after_last_round option requires out_pssm or out_ascii_pssm",
            ));
        }
        if self.save_each_pssm && self.out_ascii_pssm.is_none() {
            return Err(crate::blastinput::app::options_error(
                "save_each_pssm option requires out_pssm or out_ascii_pssm",
            ));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-321
    // ```c++
    //                         if (itr || m_CmdLineArgs->GetSaveLastPssm()) {
    //                         CConstRef<CBioseq> seq =
    //                          s_GetQueryBioseq(query, scope, pssm);
    //                         pssm =
    //                          ComputePssmForNextIteration(*seq, aln, psi_opts, scope,
    //                                      results_1st_query.GetAncillaryData());
    //     ...
    //        if (m_CmdLineArgs->GetSaveLastPssm()) {
    //            SavePssmToFile(pssm, &itr);
    //        }
    // ```
    // With one round (`itr` is false after it), the PSSM is computed only with
    // `-save_pssm_after_last_round`, and `SavePssmToFile` writes only `-out_ascii_pssm`'s file
    // here (`SaveAsciiPssm`).
    pub(crate) fn request(&self) -> Option<PssmRequest> {
        match (&self.out_ascii_pssm, self.save_pssm_after_last_round) {
            (Some(path), true) => Some(PssmRequest {
                ascii_path: path.clone(),
                save_each: self.save_each_pssm,
                pseudocount: self.pseudocount,
            }),
            _ => None,
        }
    }
}

/// The failure of the PSSM of the searched query `query` (an error that NCBI's psiblast raises
/// after that query's report, when it computes its PSSM): the run reports the queries up to it
/// without the epilog, then fails with `error` (`blast_engine.rs` `search`).
#[derive(Debug)]
pub(crate) struct PssmStop {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
    /// The PSSMs of the queries before it, which NCBI has written.
    pub(crate) pssms: Vec<Vec<u8>>,
}

impl std::fmt::Display for PssmStop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for PssmStop {}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-220
// ```c++
//     catch (const blast::CBlastException& e) {                               \
//         const string& msg = e.GetMsg();                                     \
//         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
//             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
//             exit_code = BLAST_INPUT_ERROR;                                  \
// ```
fn blast_options_error(message: &str) -> anyhow::Error {
    crate::cli::NativeError {
        exit: 1,
        message: format!("BLAST options error: {message}\n"),
    }
    .into()
}

/// One Seq-align of a query's round-1 results, in the order of the Seq-align set: the subject
/// (`hsp.s_idx`), the e-value, the 0-based starts in the query and subject records and the
/// alignment's operations.
pub(crate) struct PssmAlignment<'a> {
    pub(crate) subject: u32,
    pub(crate) e_value: f64,
    pub(crate) q_start: usize,
    pub(crate) s_start: usize,
    pub(crate) q_length: usize,
    pub(crate) ops: Option<&'a [GapEditOp]>,
}

/// The Karlin-Altschul blocks of the round's ancillary data (`CBlastAncillaryData`): `(K,
/// Lambda)` of the ungapped and gapped blocks. In the first round the PSI blocks are the
/// standard ones.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PssmAncillary {
    pub(crate) ungapped: (f64, f64),
    pub(crate) gapped: (f64, f64),
}

/// The ASCII PSSM of one query (`ComputePssmForNextIteration` and `PrintAsciiPssm`), from its
/// record's letters (the whole record, without masks), its Seq-aligns, the letters of their
/// subjects (`subject(s_idx)`: the record's letters in NCBIstdaa) and the round's ancillary data.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast.cpp:108-145
// ```c++
//     CPSIBlastOptions opts;
//     PSIBlastOptionsNew(&opts);
//     opts->pseudo_count = opts_handle.GetPseudoCount();
//     opts->inclusion_ethresh = opts_handle.GetInclusionThreshold();
//     ...
//     const SBlastSequence query_seq =
//         query_source.GetBlastSequence(0, eBlastEncodingProtein,
//                                       eNa_strand_unknown,
//                                       eSentinels, &warnings);
//     ...
//     CPsiBlastInputData input(query_seq.data.get()+1,    // skip sentinel
//                              query_seq.length-2,        // don't count sentinels
//                              alignment, database_scope,
//                              *opts.Get(),
//                              opts_handle.GetMatrixName(),
//                              ...
//     CPssmEngine engine(&input);
//     engine.SetUngappedStatisticalParams(ancillary_data);
//     CRef<CPssmWithParameters> retval(engine.Run());
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:376-399
// ```c++
//     m_PssmInput->Process();
//     s_Validate(m_PssmInput);
//
//     CPSIMatrix pssm;
//     CPSIDiagnosticsResponse diagnostics;
//     int status =
//         PSICreatePssmWithDiagnostics(m_PssmInput->GetData(),
//                                      m_PssmInput->GetOptions(),
//                                      m_ScoreBlk,
//                                      m_PssmInput->GetDiagnosticsRequest(),
//                                      &pssm,
//                                      &diagnostics);
//     if (status != PSI_SUCCESS) {
//         // FIXME: need to use core level perror-like facility
//         string msg = x_ErrorCodeToString(status);
//         NCBI_THROW(CBlastException, eCoreBlastError, msg);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:131-148
// ```c++
//     Blast_Message* errors = NULL;
//     if (PSIBlastOptionsValidate(pssm_input_msa->GetOptions(), &errors) != 0) {
//         string msg("IPssmInputData returns invalid PSIBlastOptions: ");
//         msg += string(errors->message);
//         errors = Blast_MessageFree(errors);
//         NCBI_THROW(CBlastException, eInvalidOptions, msg);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1586-1590
// ```c
//     if (psi_options->pseudo_count < 0) {
//         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//                            "Pseudo count must be greater than or equal to 0");
//         return retval;
//     }
// ```
// The inclusion threshold is checked before the search (blastp's `Validate`). The ungapped
// statistical parameters that `SetUngappedStatisticalParams` gives the score block are not read
// by the ASCII PSSM, which prints the ancillary data's.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn ascii_pssm(
    request: &PssmRequest,
    matrix: ScoringMatrix,
    inclusion_ethresh: f64,
    query_letters: &[u8],
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Vec<u8>,
    ancillary: PssmAncillary,
) -> Result<Vec<u8>> {
    let sbp = score_block(matrix)?;
    let query = crate::algorithm::blastp::encoding::encode_protein_sequence(query_letters);
    let msa = input_msa(
        &query.aa_seq[1..=query.aa_len],
        alignments,
        subject,
        inclusion_ethresh,
    );
    if request.pseudocount < 0 {
        return Err(blast_options_error(
            "IPssmInputData returns invalid PSIBlastOptions: Pseudo count must be greater than or equal to 0",
        ));
    }
    let pssm = psi::create_pssm_with_diagnostics(&msa, request.pseudocount, &sbp)
        .map_err(|error| crate::blastinput::app::engine_error(error_code_to_string(error)))?;
    Ok(print_ascii_pssm(
        &pssm,
        &ncbistdaa_letters(query_letters),
        ancillary,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:218-282
// ```c++
//     case PSIERR_BADSEQWEIGHTS:
//         retval = "Error computing sequence weights";
//         break;
//     ...
//     case PSIERR_POSITIVEAVGSCORE:
//         retval = "PSSM has positive average score";
//         break;
//
//     case PSIERR_NOALIGNEDSEQS:
//         retval = "No sequences left after purging biased sequences in ";
//         retval += "multiple sequence alignment";
//         break;
//
//     case PSIERR_GAPINQUERY:
//         retval = "Gap found in query sequence";
//         break;
//
//     case PSIERR_UNALIGNEDCOLUMN:
//         retval = "Found column with no sequences aligned in it";
//         break;
//
//     case PSIERR_COLUMNOFGAPS:
//         retval = "Found column with only GAP residues";
//         break;
//
//     case PSIERR_STARTINGGAP:
//         retval = "Found flanking gap at start of alignment";
//         break;
//
//     case PSIERR_ENDINGGAP:
//         retval = "Found flanking gap at end of alignment";
//         break;
// ```
fn error_code_to_string(error: PsiError) -> &'static str {
    match error {
        PsiError::BadSeqWeights => "Error computing sequence weights",
        PsiError::PositiveAvgScore => "PSSM has positive average score",
        PsiError::NoAlignedSeqs => {
            "No sequences left after purging biased sequences in multiple sequence alignment"
        }
        PsiError::GapInQuery => "Gap found in query sequence",
        PsiError::UnalignedColumn => "Found column with no sequences aligned in it",
        PsiError::ColumnOfGaps => "Found column with only GAP residues",
        PsiError::StartingGap => "Found flanking gap at start of alignment",
        PsiError::EndingGap => "Found flanking gap at end of alignment",
    }
}

/// The score block of `CPssmEngine::x_InitializeScoreBlock` that the engine reads.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:525-577
// ```c++
//     const EBlastProgramType kProgramType = eBlastTypePsiBlast;
//     ...
//     BlastScoringOptionsSetMatrix(opts, matrix_name);
//     ...
//     status = BlastSetup_ScoreBlkInit(query_blk,
//                                      query_info,
//                                      opts,
//                                      kProgramType,
//                                      &retval,
//                                      kScaleFactor,
//                                      &errors,
//                                      &BlastFindMatrixPath);
//     ...
//     _ASSERT(retval->kbp_ideal);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2843-2848
// ```c
//    stdrfp = Blast_ResFreqNew(sbp);
//    Blast_ResFreqStdComp(sbp, stdrfp);
//    sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
//    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
//    sbp->kbp_ideal = Blast_KarlinBlkNew();
//    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1506-1528 (the matrix's
// `loscore` and `hiscore`: its smallest and largest scores between BLAST_SCORE_MIN and
// BLAST_SCORE_MAX)
// ```c
//          if (matrix[index1][index2] <= BLAST_SCORE_MIN ||
//              matrix[index1][index2] >= BLAST_SCORE_MAX)
//             continue;
//          if (sbp->loscore > matrix[index1][index2])
//             sbp->loscore = matrix[index1][index2];
//          if (sbp->hiscore < matrix[index1][index2])
//             sbp->hiscore = matrix[index1][index2];
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/matrix_freq_ratios.c:1661-1736
// ```c
//     if ( !strcasecmp(matrix_name, "BLOSUM62") ||
//          !strcasecmp(matrix_name, "BLOSUM62_20")) {
//         ...
//         retval->bit_scale_factor = 2;
//     ...
//     } else if ( !strcasecmp(matrix_name, "BLOSUM45") ) {
//         ...
//         retval->bit_scale_factor = 3;
//     } else if ( !strcasecmp(matrix_name, "BLOSUM80") ) {
//         ...
//         retval->bit_scale_factor = 2;
//     ... (BLOSUM50, BLOSUM90, PAM30, PAM70, PAM250: 2)
//     } else {
//         retval = _PSIMatrixFrequencyRatiosFree(retval);
//     }
// ```
pub(crate) fn score_block(matrix: ScoringMatrix) -> Result<PsiScoreBlk> {
    // psiblast's `Validate` refuses IDENTITY, the one matrix without frequency ratios
    // (`stats::protein_options`, the standard matrices only).
    let data = crate::core::composition_adjustment::adjust_scores::compo_matrix_data(matrix)
        .ok_or_else(|| anyhow::anyhow!("no frequency ratios for {}", matrix.ncbi_name()))?;
    let mut table = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];
    let (mut loscore, mut hiscore) = (i32::MAX, i32::MIN);
    for (i, row) in table.iter_mut().enumerate() {
        for (j, score) in row.iter_mut().enumerate() {
            *score = crate::utils::matrix::protein_score(matrix, i as u8, j as u8);
            if *score <= BLAST_SCORE_MIN || *score >= i16::MAX as i32 {
                continue;
            }
            loscore = loscore.min(*score);
            hiscore = hiscore.max(*score);
        }
    }
    let std_comp = crate::stats::karlin_calc::compute_std_aa_composition();
    let sfp = crate::stats::karlin_calc::compute_score_freq_profile_for_matrix(
        &std_comp, &std_comp, loscore, hiscore, matrix,
    );
    let ideal = crate::stats::karlin_calc::compute_karlin_params_ungapped(&sfp)
        .map_err(|error| crate::blastinput::app::engine_error(&error))?;
    Ok(PsiScoreBlk {
        matrix: table,
        ideal_lambda: ideal.lambda,
        freq_ratios: data.freq_ratios,
        bit_scale_factor: if matrix == ScoringMatrix::Blosum45 {
            3
        } else {
            2
        },
        background: data.background,
    })
}

/// The engine's input from the query and its Seq-aligns (`CPsiBlastInputData::Process`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:122-141
// ```c++
//     // Update the number of aligned sequences
//     m_MsaDimensions.num_seqs = x_CountAndSelectQualifyingAlignments();
//
//     // Create multiple alignment data structure and populate with query
//     // sequence
//     m_Msa = PSIMsaNew(&m_MsaDimensions);
//     ...
//     x_CopyQueryToMsa();
//     x_ExtractAlignmentData();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:529-565
// ```c++
// double GetLowestEvalue(const objects::CDense_seg::TScores& scores,
//                        double* bit_score /* = NULL */)
// {
//     double retval = BLAST_EXPECT_VALUE;
//     ...
//     ITERATE(CDense_seg::TScores, i, scores) {
//         if ( (tmp = s_GetEvalue(**i)) < retval) {
//             retval = tmp;
//         }
//     ...
//     ITERATE(CSeq_align_set::Tdata, hsp, alignments.Get()) {
//         // Look for HSP with score less than inclusion_ethresh
//         double e = GetLowestEvalue((*hsp)->GetScore());
//         if (e < evalue_inclusion_threshold) {
//             CSeq_id_Handle sid =
//                 CSeq_id_Handle::GetHandle((*hsp)->GetSeq_id(1));
//             output.insert(sid);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:283-323
// ```c++
//     for (unsigned int i = 0; i < GetQueryLength(); i++) {
//         m_Msa->data[kQueryIndex][i].letter = m_Query[i];
//         m_Msa->data[kQueryIndex][i].is_aligned = true;
//     }
//     ...
//     unsigned int msa_index = kQueryIndex + 1;
//
//     CSeq_id* last_sid=NULL;
//
//     // For each HSP...
//     ITERATE(CSeq_align_set::Tdata, itr, m_SeqAlignSet->Get()) {
//
//         double bit_score;
//         double evalue = GetLowestEvalue((*itr)->GetScore(), &bit_score);
//         CSeq_id* current_sid = const_cast<CSeq_id*> (&(*itr)->GetSeq_id(1));
//
//         // ... below the e-value inclusion threshold
//         if (evalue < m_Opts.inclusion_ethresh) {
//             // Increment msa_index (if appropriate) after all CDense_seg for a
//             // given target sequence have been processed.
//             if (last_sid && !current_sid->Match(*last_sid)) {
//                 msa_index++;
//             }
//             _ASSERT(msa_index < GetNumAlignedSequences() + 1);
//             const CDense_seg& seg = (*itr)->GetSegs().GetDenseg();
//             x_ProcessDenseg(seg, msa_index, evalue, bit_score);
//         }
//         last_sid = current_sid;
//     }
// ```
// A Seq-align names its subject by its Seq-id; LOSAT's names it by the subject's index, which
// the Seq-ids of different subjects do not share in the runs LOSAT accepts. An HSP after a
// subject that is not counted (an HSP above the threshold before the first one below it) moves
// `msa_index` past the rows NCBI allocates; NCBI writes there out of bounds, LOSAT drops it.
fn input_msa(
    query: &[u8],
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Vec<u8>,
    inclusion_ethresh: f64,
) -> PsiMsa {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:158
    // ```c
    // #define BLAST_EXPECT_VALUE 10.0 /**< by default, alignments whose expect
    // ```
    const BLAST_EXPECT_VALUE: f64 = 10.0;
    let lowest_evalue = |e_value: f64| {
        if e_value < BLAST_EXPECT_VALUE {
            e_value
        } else {
            BLAST_EXPECT_VALUE
        }
    };
    let mut qualifying: Vec<u32> = alignments
        .iter()
        .filter(|alignment| lowest_evalue(alignment.e_value) < inclusion_ethresh)
        .map(|alignment| alignment.subject)
        .collect();
    qualifying.sort_unstable();
    qualifying.dedup();
    let mut msa = PsiMsa::new(query.len(), qualifying.len());
    for (cell, &letter) in msa.rows[0].iter_mut().zip(query) {
        *cell = PsiMsaCell {
            letter,
            is_aligned: true,
        };
    }
    let mut msa_index = 1;
    let mut last_subject: Option<u32> = None;
    for alignment in alignments {
        if lowest_evalue(alignment.e_value) < inclusion_ethresh {
            if last_subject.is_some_and(|last| last != alignment.subject) {
                msa_index += 1;
            }
            if let Some(row) = msa.rows.get_mut(msa_index) {
                process_denseg(row, alignment, &subject(alignment.subject));
            }
        }
        last_subject = Some(alignment.subject);
    }
    msa
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:363-408
// ```c++
//     for (int segmt_idx = 0; segmt_idx < kNumSegments; segmt_idx++) {
//
//         TSeqPos query_offset = starts[query_index];
//         TSeqPos subject_offset = starts[subj_index];
//         ...
//         if (query_offset == GAP_IN_ALIGNMENT) {
//
//             // gap in query, just skip residues on subject sequence
//             subj_seq_idx += lengths[segmt_idx];
//             continue;
//
//         } else if (subject_offset == GAP_IN_ALIGNMENT) {
//
//             // gap in subject, initialize appropriately
//             for (TSeqPos i = 0; i < lengths[segmt_idx]; i++) {
//                 PSIMsaCell& msa_cell = m_Msa->data[msa_index][query_offset++];
//                 if ( !msa_cell.is_aligned ) {
//                     msa_cell.letter = GAP;
//                     msa_cell.is_aligned = true;
//                 }
//             }
//
//         } else {
//
//             // Aligned segments without any gaps
//             for (TSeqPos i = 0; i < lengths[segmt_idx]; i++, subj_seq_idx++) {
//                 PSIMsaCell& msa_cell =
//                     m_Msa->data[msa_index][query_offset++];
//                 if ( !msa_cell.is_aligned ) {
//                     msa_cell.letter = static_cast<Uint1>(seq[subj_seq_idx]);
//                     msa_cell.is_aligned = true;
//                 }
//             }
//         }
//
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:425-447
// ```c++
//     for (int i = 0; i < kNumSegments; i++) {
//
//         if (starts[subj_index] != (TSignedSeqPos)GAP_IN_ALIGNMENT) {
//             if ( !subj_start_found ) {
//                 subj_start = starts[subj_index];
//                 subj_start_found = true;
//             }
//             subjlen += lengths[i];
//         }
//
//         subj_index += kDimensions;
//     }
//     ...
//     CSeq_loc seqloc(const_cast<CSeq_id&>(*ds.GetIds().back()), subj_start,
//                     subj_start+subjlen-1);
//     ...
//         CSeqVector sv(seqloc, scope);
//         sv.SetCoding(CSeq_data::e_Ncbistdaa);
//         sv.GetSeqData(0, kInvalidSeqPos, sequence_data);
// ```
// The subject's letters from the first aligned one on: an HSP's operations walk them from
// `s_start` (an HSP starts with aligned letters). An HSP without operations is ungapped.
fn process_denseg(row: &mut [PsiMsaCell], alignment: &PssmAlignment<'_>, subject: &[u8]) {
    let whole = [GapEditOp::Sub(alignment.q_length as u32)];
    let ops = match alignment.ops {
        Some(ops) if !ops.is_empty() => ops,
        _ => &whole[..],
    };
    let mut q = alignment.q_start;
    let mut s = alignment.s_start;
    for op in ops {
        match *op {
            GapEditOp::Del(n) => s += n as usize,
            GapEditOp::Ins(n) => {
                for cell in &mut row[q..q + n as usize] {
                    if !cell.is_aligned {
                        *cell = PsiMsaCell {
                            letter: ncbistdaa::GAP,
                            is_aligned: true,
                        };
                    }
                }
                q += n as usize;
            }
            GapEditOp::Sub(n) => {
                for (cell, &letter) in row[q..q + n as usize]
                    .iter_mut()
                    .zip(&subject[s..s + n as usize])
                {
                    if !cell.is_aligned {
                        *cell = PsiMsaCell {
                            letter,
                            is_aligned: true,
                        };
                    }
                }
                q += n as usize;
                s += n as usize;
            }
        }
    }
}

/// A record's letters in NCBIstdaa as its Bioseq holds them (`CSeqVector`,
/// `CPssm::GetQuerySequenceData`): upper case, and O stays O (only the search reads it as X).
pub(crate) fn ncbistdaa_letters(letters: &[u8]) -> Vec<u8> {
    letters
        .iter()
        .map(|&letter| aa_char_to_ncbistdaa(letter.to_ascii_uppercase()))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:145-166
// ```c++
// static int RESIDUE_ORDER[] = {
//      1,     /* A */
//     16,     /* R */
//     13,     /* N */
//      4,     /* D */
//      3,     /* C */
//     15,     /* Q */
//      5,     /* E */
//      7,     /* G */
//      8,     /* H */
//      9,     /* I */
//     11,     /* L */
//     10,     /* K */
//     12,     /* M */
//      6,     /* F */
//     14,     /* P */
//     17,     /* S */
//     18,     /* T */
//     20,     /* W */
//     22,     /* Y */
//     19      /* V */
// };
// ```
const RESIDUE_ORDER: [usize; 20] = [
    1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-118
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'};
// ```
const NCBISTDAA_TO_AMINOACID: &[u8; BLASTAA_SIZE] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

/// A double as the C++ stream writes it with `fixed` and `setprecision(precision)` (the C
/// library's `%.*f`: correctly rounded, as Rust's formatting; `inf` and `nan` as glibc spells
/// them).
fn fixed(value: f64, precision: usize) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .to_string();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    format!("{value:.precision$}")
}

/// `CBlastFormatUtil::PrintAsciiPssm` of a computed PSSM.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:191-327
// ```c++
//     if (pssm_calculation_done) {
//         out << "\nLast position-specific scoring matrix computed, weighted ";
//         out << "observed percentages rounded down, information per position, ";
//         out << "and relative weight of gapless real matches to pseudocounts\n";
//     ...
//     int max_score = 0;
//     ITERATE (TNcbiMatrixInt::TData, it, pssm->GetData()) {
//         if (*it <= BLAST_SCORE_MIN) {
//             continue;
//         }
//         if (*it > max_score) {
//             max_score = *it;
//         }
//         if (-*it > max_score) {
//             max_score = -*it;
//         }
//     }
//     int num_digits = 0;
//     while (max_score > 0) {
//         max_score /= 10;
//         num_digits++;
//     }
//     int width = num_digits + 2;
//
//     out << "         ";
//     for (size_t c = 0; c < DIM(RESIDUE_ORDER); c++) {
//         out << setw(width) << NCBISTDAA_TO_AMINOACID[RESIDUE_ORDER[c]];
//     }
//     if (pssm_calculation_done) {
//         for (size_t c = 0; c < DIM(RESIDUE_ORDER); c++) {
//             out << "   " << NCBISTDAA_TO_AMINOACID[RESIDUE_ORDER[c]];
//         }
//     }
//     ...
//     out << fixed;
//     for (SIZE_TYPE i = 0; i < kQueryLength; i++) {
//         out << "\n" << setw(5) << (i+1) << " " <<
//             NCBISTDAA_TO_AMINOACID[(int)query_seq[i]] << "  ";
//         for (SIZE_TYPE c = 0; c < DIM(RESIDUE_ORDER); c++) {
//             if ((*pssm)(RESIDUE_ORDER[c], i) == BLAST_SCORE_MIN) {
//                 out << "-I ";
//             } else {
//                 out << setw(width) << (*pssm)(RESIDUE_ORDER[c], i);
//             }
//         }
//         out << " ";
//         if (pssm_calculation_done) {
//             for (SIZE_TYPE c = 0; c < DIM(RESIDUE_ORDER); c++) {
//                 if ((*pssm)(RESIDUE_ORDER[c], i) != BLAST_SCORE_MIN) {
//                     double value = 100;
//                     value *= (*weighted_res_freqs)(RESIDUE_ORDER[c], i);
//                     // round to the nearest integer
//                     value = (int)(value + (value > 0. ? 0.5 : -0.5));
//                     out << setw(4) << (int)value;
//                 }
//             }
//             out << "  " << setprecision(2) << info_content[i] << " ";
//             if ((num_matching_seqs[i] > 1) && (query_seq[i] != kXResidue)) {
//                 out << setprecision(2) << gapless_col_weights[i];
//             } else {
//                 out << "    0.00";
//             }
//         }
//     }
//     ...
//     out << "\n\n" << setprecision(4);
//     out << "                      K         Lambda\n";
//     if (ungapped_kbp) {
//         out << "Standard Ungapped    "
//             << ungapped_kbp->K << "     "
//             << ungapped_kbp->Lambda << "\n";
//     }
//     if (gapped_kbp) {
//         out << "Standard Gapped      "
//             << gapped_kbp->K << "     "
//             << gapped_kbp->Lambda << "\n";
//     }
//     if (psi_ungapped_kbp) {
//         out << "PSI Ungapped         "
//             << psi_ungapped_kbp->K << "     "
//             << psi_ungapped_kbp->Lambda << "\n";
//     }
//     if (psi_gapped_kbp) {
//         out << "PSI Gapped           "
//             << psi_gapped_kbp->K << "     "
//             << psi_gapped_kbp->Lambda << "\n";
//     }
// ```
// The scores matrix of `GetScores` has the 28 rows of NCBIstdaa (pssm_engine.cpp
// `x_PSIMatrix2Asn1` stores every row). The query's letters are those of its Bioseq
// (`GetQuerySequenceData`, blastfmtutil.cpp:253-255), in which an O stays O.
pub(crate) fn print_ascii_pssm(pssm: &PsiPssm, query: &[u8], ancillary: PssmAncillary) -> Vec<u8> {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push_str("\nLast position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n");
    let mut max_score: i32 = 0;
    for &score in pssm.scores.iter().flatten() {
        if score <= BLAST_SCORE_MIN {
            continue;
        }
        if score > max_score {
            max_score = score;
        }
        if -score > max_score {
            max_score = -score;
        }
    }
    let mut num_digits = 0;
    while max_score > 0 {
        max_score /= 10;
        num_digits += 1;
    }
    let width = num_digits + 2;
    out.push_str("         ");
    for &c in &RESIDUE_ORDER {
        let _ = write!(out, "{:>width$}", NCBISTDAA_TO_AMINOACID[c] as char);
    }
    for &c in &RESIDUE_ORDER {
        let _ = write!(out, "   {}", NCBISTDAA_TO_AMINOACID[c] as char);
    }
    for (i, &residue) in query.iter().enumerate() {
        let scores = &pssm.scores[i];
        let _ = write!(
            out,
            "\n{:>5} {}  ",
            i + 1,
            NCBISTDAA_TO_AMINOACID[residue as usize] as char
        );
        for &c in &RESIDUE_ORDER {
            if scores[c] == BLAST_SCORE_MIN {
                out.push_str("-I ");
            } else {
                let _ = write!(out, "{:>width$}", scores[c]);
            }
        }
        out.push(' ');
        for &c in &RESIDUE_ORDER {
            if scores[c] != BLAST_SCORE_MIN {
                let mut value = 100.0;
                value *= pssm.weighted_residue_freqs[i][c];
                value = (value + if value > 0.0 { 0.5 } else { -0.5 }) as i32 as f64;
                let _ = write!(out, "{:>4}", value as i32);
            }
        }
        let _ = write!(out, "  {} ", fixed(pssm.information_content[i], 2));
        if pssm.num_matching_seqs[i] > 1 && residue != ncbistdaa::X {
            out.push_str(&fixed(pssm.gapless_column_weights[i], 2));
        } else {
            out.push_str("    0.00");
        }
    }
    out.push_str("\n\n");
    out.push_str("                      K         Lambda\n");
    let (ungapped_k, ungapped_lambda) = ancillary.ungapped;
    let (gapped_k, gapped_lambda) = ancillary.gapped;
    for (label, k, lambda) in [
        ("Standard Ungapped    ", ungapped_k, ungapped_lambda),
        ("Standard Gapped      ", gapped_k, gapped_lambda),
        ("PSI Ungapped         ", ungapped_k, ungapped_lambda),
        ("PSI Gapped           ", gapped_k, gapped_lambda),
    ] {
        let _ = writeln!(out, "{label}{}     {}", fixed(k, 4), fixed(lambda, 4));
    }
    out.into_bytes()
}

/// Writes the ASCII PSSMs of the queries whose round-1 results have alignments, in the order of
/// the queries: each write opens the file again (`CAutoOutputFileReset::GetStream`), so without
/// `-save_each_pssm` the last PSSM stays, and with it the writes are `F.1`, `F.2`, ...
/// A file that cannot be opened or written is not reported, as NCBI's `ofstream` is not checked.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:398-407
// ```c++
//     if (m_CmdLineArgs->SaveAsciiPssm() &&
//         (itr == NULL || // this is true in the case of remote PSI-BLAST
//          itr->GetIterationNumber() >= 1)) {
//         ...
//         CBlastFormatUtil::PrintAsciiPssm(*pssm,
//                                          m_AncillaryData,
//                                          *m_CmdLineArgs->GetAsciiPssmStream());
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:50-66
// ```c++
// CAutoOutputFileReset::GetStream()
// {
//     string filename = m_FileName;
//     if (!m_Version) {
//         CFile file_deleter(m_FileName);
//         if (file_deleter.Exists()) {
//             _TRACE("Deleting " << m_FileName);
//             file_deleter.Remove();
//         }
//     }
//     else {
//         filename = m_FileName + "." + NStr::IntToString(m_Version);
//         m_Version++;
//     }
//
//     m_FileStream.reset(new ofstream(filename.c_str()));
//     return m_FileStream.get();
// }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:54-55
// ```c++
//     CAutoOutputFileReset(const string& file_name, bool use_versions = false)
//         : m_FileName(file_name), m_Version(use_versions ? 1 : 0) {}
// ```
#[cold]
#[inline(never)]
pub(crate) fn write_ascii_pssm_files<'a>(
    request: &PssmRequest,
    pssms: impl Iterator<Item = &'a [u8]>,
) {
    let mut version: u32 = if request.save_each { 1 } else { 0 };
    for pssm in pssms {
        let filename = if version == 0 {
            let path = std::path::Path::new(&request.ascii_path);
            if path.is_file() {
                let _ = std::fs::remove_file(path);
            }
            request.ascii_path.clone()
        } else {
            let filename = format!("{}.{}", request.ascii_path, version);
            version += 1;
            filename
        };
        if let Ok(mut file) = std::fs::File::create(&filename) {
            let _ = file.write_all(pssm);
        }
    }
}

#[cfg(test)]
#[path = "pssm_tests.rs"]
mod tests;
