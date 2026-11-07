//! psiblast's PSSM after the last round (`-num_iterations 1 -save_pssm_after_last_round` with
//! `-out_pssm F` and/or `-out_ascii_pssm F`): the input of NCBI's PSSM engine from a query's
//! round-1 alignments (`CPsiBlastInputData`), the engine's score block (`CPssmEngine`), the
//! ASN.1 PSSM (`x_PSIMatrix2Asn1`, written by `report::serial::pssm`), the ASCII PSSM
//! (`CBlastFormatUtil::PrintAsciiPssm`) and their files (`CAutoOutputFileReset`). The engine
//! itself is `core::psi`. The blastp engine reaches this module only when a PSSM is requested
//! (`blast_engine.rs` `psiblast_pssms`).

use std::io::Write as _;

use anyhow::Result;

use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::core::psi::{
    self, KarlinBlk, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiMsaCell, PsiPssm, PsiScoreBlk,
};
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmSeqId, PssmWithParameters, SeqData, SeqInst,
};
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

/// psiblast's PSSM options as given on the command line.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PsiblastPssmOptions {
    pub save_pssm_after_last_round: bool,
    pub save_each_pssm: bool,
    pub out_pssm: Option<String>,
    pub out_ascii_pssm: Option<String>,
    pub pseudocount: i32,
    /// `-in_pssm`: the query as a PSSM file (`algorithm::psiblast::in_pssm`).
    pub in_pssm: Option<std::path::PathBuf>,
    /// `-num_iterations` when it is not 1 (0: until convergence): the rounds of a FASTA query
    /// (`algorithm::psiblast::rounds`); `None` for one round.
    pub rounds: Option<u32>,
}

/// A run that writes the PSSM of each query whose round-1 results have alignments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PssmRequest {
    /// `-out_pssm`, as given: the ASN.1 text PSSM (`SaveCheckpoint`).
    pub asn_path: Option<String>,
    /// `-out_ascii_pssm`, as given (`SaveAsciiPssm`).
    pub ascii_path: Option<String>,
    /// `-save_each_pssm`: the files are `F.1`, `F.2`, ... (`CAutoOutputFileReset` with versions).
    pub save_each: bool,
    /// `-pseudocount` (`PSIBlastOptions::pseudo_count`).
    pub pseudocount: i32,
}

impl PssmRequest {
    /// The diagnostics the engine computes (`PSIDiagnosticsRequestNewEx(SaveAsciiPssm())`).
    pub(crate) fn diagnostics(&self) -> PsiDiagnosticsRequest {
        PsiDiagnosticsRequest::new_ex(self.ascii_path.is_some())
    }
}

impl PsiblastPssmOptions {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1283-1303
    // ```c
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
    // (blast_args.cpp:3826-3832).
    pub(crate) fn check(&self) -> Result<()> {
        let no_file = self.out_pssm.is_none() && self.out_ascii_pssm.is_none();
        if self.save_pssm_after_last_round && no_file {
            return Err(crate::blastinput::app::options_error(
                "save_pssm_after_last_round option requires out_pssm or out_ascii_pssm",
            ));
        }
        if self.save_each_pssm && no_file {
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
    // `-save_pssm_after_last_round`, and `SavePssmToFile` writes the files of `-out_pssm` and
    // `-out_ascii_pssm` (`SaveCheckpoint`, `SaveAsciiPssm`).
    // With more rounds the driver writes the files (`rounds::RoundFiles`), not the search.
    pub(crate) fn request(&self) -> Option<PssmRequest> {
        if self.rounds.is_some()
            || !self.save_pssm_after_last_round
            || (self.out_pssm.is_none() && self.out_ascii_pssm.is_none())
        {
            return None;
        }
        Some(PssmRequest {
            asn_path: self.out_pssm.clone(),
            ascii_path: self.out_ascii_pssm.clone(),
            save_each: self.save_each_pssm,
            pseudocount: self.pseudocount,
        })
    }
}

/// The failure of the PSSM of the searched query `query` (an error that NCBI's psiblast raises
/// after that query's report, when it computes or saves its PSSM): the run reports the queries up
/// to it without the epilog, then fails with `error` (`blast_engine.rs` `search`).
#[derive(Debug)]
pub(crate) struct PssmStop {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
    /// The PSSM files of the queries before it, which NCBI has written: LOSAT writes them after
    /// the reports of the queries up to it (`PssmFiles::commit`).
    pub(crate) files: PssmFiles,
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

/// The Karlin-Altschul blocks of the round's ancillary data (`CBlastAncillaryData`, the blocks of
/// the query's context): `kbp_std`, `kbp_gap`, `kbp_psi` and `kbp_gap_psi`. In the first round
/// the PSI blocks are copies of the standard ones (`kbp_psi` and `kbp_gap_psi` of the search's
/// score block are its standard blocks, blast_setup.c:114-120); in a round that searched a PSSM
/// they are the PSSM's, and `kbp_gap` is `kbp_gap_psi` (blast_setup.c:113-126, quoted at
/// `algorithm::pssm_search`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:99-111
// ```c++
//     if (sbp->kbp_std) {
//         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
//     }
//     if (sbp->kbp_gap) {
//         s_InitializeKarlinBlk(sbp->kbp_gap[ctx_index], &m_GappedKarlinBlk);
//     }
//     if (sbp->kbp_psi) {
//         s_InitializeKarlinBlk(sbp->kbp_psi[ctx_index], &m_PsiUngappedKarlinBlk);
//     }
//     if (sbp->kbp_gap_psi) {
//         s_InitializeKarlinBlk(sbp->kbp_gap_psi[ctx_index],
//                               &m_PsiGappedKarlinBlk);
//     }
// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct PssmAncillary {
    pub(crate) ungapped: KarlinBlk,
    pub(crate) gapped: KarlinBlk,
    pub(crate) psi_ungapped: KarlinBlk,
    pub(crate) psi_gapped: KarlinBlk,
}

/// The query of a PSSM: its letters (the whole record, without masks), its Seq-align id
/// (`Query_<n>`) and its title.
pub(crate) struct PssmQueryRecord<'a> {
    pub(crate) letters: &'a [u8],
    pub(crate) id: &'a str,
    pub(crate) title: &'a [u8],
}

/// The scoring of the search that the PSSM records: the matrix as named (`-matrix`, any case)
/// and the gap costs.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PssmScoring<'a> {
    pub(crate) matrix: ScoringMatrix,
    pub(crate) matrix_name: &'a str,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
}

/// The files of one query's PSSM: the ASN.1 text with `-out_pssm`, the ASCII PSSM with
/// `-out_ascii_pssm`.
#[derive(Debug)]
pub(crate) struct QueryPssm {
    pub(crate) asn: Option<Vec<u8>>,
    pub(crate) ascii: Option<Vec<u8>>,
    /// The engine's PSSM, which the next round searches (`algorithm::psiblast::rounds`).
    pub(crate) engine: PsiPssm,
}

/// The PSSM of one query (`ComputePssmForNextIteration`) and its files: the ASN.1 text
/// (`x_PSIMatrix2Asn1`, `PsiBlastAddAncillaryPssmData` and `MSerial_AsnText`) and the ASCII
/// PSSM (`PrintAsciiPssm`), as requested. The inputs are the query, its Seq-aligns, the letters
/// of their subjects (`subject(s_idx)`: the record's letters in NCBIstdaa), the round's
/// ancillary data, and the gapped block of the matrix and gap costs (`gapped`, the score block's
/// `kbp_gap_std[0]`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast.cpp:110-150
// ```c++
//     CPSIBlastOptions opts;
//     PSIBlastOptionsNew(&opts);
//     opts->pseudo_count = opts_handle.GetPseudoCount();
//     opts->inclusion_ethresh = opts_handle.GetInclusionThreshold();
//
//     string query_descr = NcbiEmptyString;
//
//     if (query.IsSetDescr()) {
//          const CBioseq::TDescr::Tdata& data = query.GetDescr().Get();
//          ITERATE(CBioseq::TDescr::Tdata, iter, data) {
//              if((*iter)->IsTitle()) {
//                  query_descr += (*iter)->GetTitle();
//              }
//          }
//     }
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
//                              opts_handle.GetGapOpeningCost(),
//                              opts_handle.GetGapExtensionCost(),
//                              diagnostics_request,
//                              query_descr);
//
//     CPssmEngine engine(&input);
//     engine.SetUngappedStatisticalParams(ancillary_data);
//     CRef<CPssmWithParameters> retval(engine.Run());
//
//     PsiBlastAddAncillaryPssmData(*retval,
//                                   opts_handle.GetGapOpeningCost(),
//                                   opts_handle.GetGapExtensionCost());
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:380-404
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
//
//     // Convert core BLAST matrix structure into ASN.1 score matrix object
//     CRef<CPssmWithParameters> retval;
//     retval = x_PSIMatrix2Asn1(pssm, m_PssmInput->GetMatrixName(),
//                               m_PssmInput->GetOptions(), diagnostics);
//     CRef<CBioseq> query = m_PssmInput->GetQueryForPssm();
//     if (query.NotEmpty()) {
//         retval->SetQuery().SetSeq(*query);
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1577-1600
// ```c
// Int2 PSIBlastOptionsValidate(const PSIBlastOptions* psi_options,
//                              Blast_Message** blast_msg)
// {
//     Int2 retval = 1;    /* assume failure */
//
//     if ( !psi_options ) {
//         return retval;
//     }
//
//     if (psi_options->pseudo_count < 0) {
//         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//                            "Pseudo count must be greater than or equal to 0");
//         return retval;
//     }
//
//     if (psi_options->inclusion_ethresh <= 0.0) {
//         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//                            "Inclusion threshold must be greater than 0");
//         return retval;
//     }
//
//     retval = 0;
//     return retval;
// }
// ```
// Every check of `PSIBlastOptionsValidate`, in its order (the options are never NULL here). The
// search before it does not check the inclusion threshold: nothing else calls this function, so
// a run with `-inclusion_ethresh 0` reports the query first. The ASCII PSSM prints the ancillary
// data's blocks; the ASN.1 PSSM keeps those of the scaling (`core::psi`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:894-903
// ```c++
//     for (i = 0; i < sv.size(); i++) {
//         // Change unsupported residues to X
//         if (sv[i] == AMINOACID_TO_NCBISTDAA[(int)'O']) {
//             replaced_residues.push_back(i);
//             *buf_var++ = AMINOACID_TO_NCBISTDAA[(int)'X'];
//         } else if (!s_IsValidResidue(sv[i])) {
//             invalid_residues.push_back(i);
//         } else {
//             *buf_var++ = sv[i];
//         }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:98-99 and 166-171
// ```c++
//     m_Query = new Uint1[query_length];
//     memcpy((void*) m_Query, (void*) query, query_length);
//     ...
//     // set the sequence data in ncbistdaa format
//     CNCBIstdaa& seq = m_QueryBioseq->SetInst().SetSeq_data().SetNcbistdaa();
//     seq.Set().reserve(GetQueryLength());
//     for (TSeqPos i = 0; i < GetQueryLength(); i++) {
//         seq.Set().push_back(m_Query[i]);
//     }
// ```
// The MSA's query row and the PSSM's query Bioseq, which `PrintAsciiPssm` prints, are both the
// search's buffer (`GetBlastSequence` with `eBlastEncodingProtein`): an O is X in both. That is
// the only letter the buffer changes from the record's NCBIstdaa letters, and
// `encode_protein_sequence` makes the same change, so the letters are the buffer's (NCBI
// 2.17.0 and LOSAT give the same PSSMs for a query with lower-case letters, a lower-case o, and
// a `-` that the FASTA reader drops). The query has one title (the FASTA reader's), empty when
// the definition line has none, and `CPsiBlastInputData` writes it as the Bioseq's title also
// when it is empty.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn query_pssm(
    request: &PssmRequest,
    scoring: PssmScoring<'_>,
    gapped: KarlinBlk,
    inclusion_ethresh: f64,
    query: PssmQueryRecord<'_>,
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Vec<u8>,
    ancillary: PssmAncillary,
) -> Result<QueryPssm> {
    let mut sbp = score_block(scoring.matrix, gapped)?;
    sbp.set_ungapped_statistical_params(ancillary.psi_gapped);
    let encoded = crate::algorithm::blastp::encoding::encode_protein_sequence(query.letters);
    let letters = &encoded.aa_seq[1..=encoded.aa_len];
    let msa = input_msa(letters, alignments, subject, inclusion_ethresh)?;
    if request.pseudocount < 0 {
        return Err(blast_options_error(
            "IPssmInputData returns invalid PSIBlastOptions: Pseudo count must be greater than or equal to 0",
        ));
    }
    if inclusion_ethresh <= 0.0 {
        return Err(blast_options_error(
            "IPssmInputData returns invalid PSIBlastOptions: Inclusion threshold must be greater than 0",
        ));
    }
    let pssm =
        psi::create_pssm_with_diagnostics(&msa, request.pseudocount, &sbp, request.diagnostics())
            .map_err(|error| crate::blastinput::app::engine_error(error_code_to_string(error)))?;
    let asn = match request.asn_path {
        Some(_) => {
            let value = pssm_with_parameters(&pssm, letters, &query, scoring, request.pseudocount);
            Some(
                crate::report::serial::pssm::asn_text(&value)
                    .map_err(|error| anyhow::anyhow!("{error}"))?,
            )
        }
        None => None,
    };
    let ascii = request
        .ascii_path
        .as_ref()
        .map(|_| print_ascii_pssm(&pssm, letters, ancillary));
    Ok(QueryPssm {
        asn,
        ascii,
        engine: pssm,
    })
}

/// The `PssmWithParameters` of a PSSM (`x_PSIMatrix2Asn1` with its query, then
/// `PsiBlastAddAncillaryPssmData`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:632-661
// ```c++
//     CRef<CPssmWithParameters> retval(new CPssmWithParameters);
//
//     // Record the parameters
//     string mtx(matrix_name);
//     mtx = NStr::ToUpper(mtx); // save the matrix name in all capital letters
//     retval->SetParams().SetRpsdbparams().SetMatrixName(mtx);
//     if (opts) {
//         retval->SetParams().SetPseudocount(opts->pseudo_count);
//     }
//
//     CPssm& asn1_pssm = retval->SetPssm();
//     asn1_pssm.SetIsProtein(true);
//     // number of rows is alphabet size
//     asn1_pssm.SetNumRows(pssm->nrows);
//     // number of columns is query length
//     asn1_pssm.SetNumColumns(pssm->ncols);
//     asn1_pssm.SetByRow(false);  // this is the default
//
//     asn1_pssm.SetLambda(pssm->lambda);
//     asn1_pssm.SetKappa(pssm->kappa);
//     asn1_pssm.SetH(pssm->h);
//     asn1_pssm.SetLambdaUngapped(pssm->ung_lambda);
//     asn1_pssm.SetKappaUngapped(pssm->ung_kappa);
//     asn1_pssm.SetHUngapped(pssm->ung_h);
//     if (asn1_pssm.GetByRow() == false) {
//         for (unsigned int i = 0; i < pssm->ncols; i++) {
//             for (unsigned int j = 0; j < pssm->nrows; j++) {
//                 asn1_pssm.SetFinalData().SetScores().
//                     push_back(pssm->pssm[i][j]);
//             }
//         }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:684-778
// ```c++
//     if (diagnostics->information_content) {
//         CPssmIntermediateData::TInformationContent& info_content =
//             asn1_pssm.SetIntermediateData().SetInformationContent();
//         for (Uint4 i = 0; i < diagnostics->query_length; i++) {
//             info_content.push_back(diagnostics->information_content[i]);
//         }
//     }
//     ...
//     if (diagnostics->weighted_residue_freqs) {
//         CPssmIntermediateData::TWeightedResFreqsPerPos& wres_freqs =
//             asn1_pssm.SetIntermediateData().SetWeightedResFreqsPerPos();
//         if (asn1_pssm.GetByRow() == false) {
//             for (unsigned int i = 0; i < pssm->ncols; i++) {
//                 for (unsigned int j = 0; j < pssm->nrows; j++) {
//                     wres_freqs.
//                         push_back(diagnostics->weighted_residue_freqs[i][j]);
//     ...
//     if (diagnostics->frequency_ratios) {
//         CPssmIntermediateData::TFreqRatios& freq_ratios =
//             asn1_pssm.SetIntermediateData().SetFreqRatios();
//         if (asn1_pssm.GetByRow() == false) {
//             for (unsigned int i = 0; i < pssm->ncols; i++) {
//                 for (unsigned int j = 0; j < pssm->nrows; j++) {
//                     freq_ratios.push_back(diagnostics->frequency_ratios[i][j]);
//     ...
//     if (diagnostics->gapless_column_weights) {
//         ...
//             gcw.push_back(diagnostics->gapless_column_weights[i]);
//     ...
//     if (diagnostics->sigma) {
//         ...
//             sigma.push_back(diagnostics->sigma[i]);
//     ...
//     if (diagnostics->interval_sizes) {
//         ...
//             interval_sizes.push_back(diagnostics->interval_sizes[i]);
//     ...
//     if (diagnostics->num_matching_seqs) {
//         ...
//             num_matching_seqs.push_back(diagnostics->num_matching_seqs[i]);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:426-433
// ```c++
// PsiBlastAddAncillaryPssmData(objects::CPssmWithParameters& pssm,
//                          int gap_open,
//                          int gap_extend)
// {
//     _ASSERT(pssm.GetParams().GetRpsdbparams().IsSetMatrixName());
//     pssm.SetParams().SetRpsdbparams().SetGapOpen(gap_open);
//     pssm.SetParams().SetRpsdbparams().SetGapExtend(gap_extend);
// }
// ```
// psiblast's options have no IMPALA scaling (`scalingFactor` is not set), and its diagnostics
// request never asks for the residue frequencies or the independent observations. The interval
// sizes and numbers of matching sequences are `Uint4` pushed to a `vector<int>`.
fn pssm_with_parameters(
    pssm: &PsiPssm,
    letters: &[u8],
    query: &PssmQueryRecord<'_>,
    scoring: PssmScoring<'_>,
    pseudocount: i32,
) -> PssmWithParameters {
    let by_column =
        |rows: &[[f64; BLASTAA_SIZE]]| -> Vec<f64> { rows.iter().flatten().copied().collect() };
    let diagnostics = &pssm.diagnostics;
    let intermediate_data = PssmIntermediateData {
        res_freqs_per_pos: None,
        weighted_res_freqs_per_pos: diagnostics.weighted_residue_freqs.as_deref().map(by_column),
        freq_ratios: diagnostics
            .frequency_ratios
            .as_deref()
            .map(by_column)
            .unwrap_or_default(),
        information_content: diagnostics.information_content.clone(),
        gapless_column_weights: diagnostics.gapless_column_weights.clone(),
        sigma: diagnostics.sigma.clone(),
        interval_sizes: diagnostics
            .interval_sizes
            .as_ref()
            .map(|sizes| sizes.iter().map(|&size| size as i32).collect()),
        num_matching_seqs: diagnostics
            .num_matching_seqs
            .as_ref()
            .map(|counts| counts.iter().map(|&count| count as i32).collect()),
        num_indept_obsr: None,
    };
    PssmWithParameters {
        pssm: Pssm {
            is_protein: Some(true),
            identifier: None,
            num_rows: BLASTAA_SIZE as i32,
            num_columns: letters.len() as i32,
            row_labels: None,
            by_row: Some(false),
            query: Some(PssmQuery {
                ids: vec![PssmSeqId::Local(ObjectId::Str(
                    query.id.as_bytes().to_vec(),
                ))],
                titles: Some(vec![query.title.to_vec()]),
                inst: SeqInst {
                    repr: ("raw", 2),
                    mol: ("aa", 3),
                    length: Some(letters.len() as i32),
                    seq_data: Some(SeqData::Ncbistdaa(letters.to_vec())),
                },
            }),
            intermediate_data: Some(intermediate_data),
            final_data: Some(PssmFinalData {
                scores: pssm.scores.iter().flatten().copied().collect(),
                lambda: pssm.gapped.lambda,
                kappa: pssm.gapped.k,
                h: pssm.gapped.h,
                scaling_factor: None,
                lambda_ungapped: Some(pssm.ungapped.lambda),
                kappa_ungapped: Some(pssm.ungapped.k),
                h_ungapped: Some(pssm.ungapped.h),
                word_score_threshold: None,
            }),
        },
        params: Some(PssmParameters {
            pseudocount: Some(pseudocount),
            rpsdbparams: Some(FormatRpsDbParameters {
                matrix_name: scoring.matrix_name.to_ascii_uppercase().into_bytes(),
                gap_open: Some(scoring.gap_open),
                gap_extend: Some(scoring.gap_extend),
            }),
            bit_score_thresh: None,
            bit_score_reporting_thresh: None,
        }),
    }
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:104-120
// ```c
//         } else {
//             retval =
//                 Blast_KarlinBlkGappedCalc(sbp->kbp_gap_std[index],
//                     scoring_options->gap_open, scoring_options->gap_extend,
//                     sbp->name, error_return);
//         }
//         ...
//         /* For right now, copy the contents from kbp_gap_std to
//          * kbp_gap_psi (as in old code - BLASTSetUpSearchInternalByLoc) */
//         if (program != eBlastTypeBlastn && program != eBlastTypeMapping) {
//             sbp->kbp_gap_psi[index] = Blast_KarlinBlkNew();
//             Blast_KarlinBlkCopy(sbp->kbp_gap_psi[index],
//                                 sbp->kbp_gap_std[index]);
//         }
// ```
// `gapped` is the gapped block of the matrix and gap costs (`Blast_KarlinBlkGappedCalc`), the
// search's: `kbp_gap_std[0]`, copied to `kbp_gap_psi[0]`.
pub(crate) fn score_block(matrix: ScoringMatrix, gapped: KarlinBlk) -> Result<PsiScoreBlk> {
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
        ideal_k: ideal.k,
        gap_std_k: gapped.k,
        gap_psi: gapped,
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
// the Seq-ids of different subjects do not share in the runs LOSAT accepts.
// `msa_index` also moves on after a subject that is not counted, so an HSP below the threshold
// whose subject comes after one without such an HSP, before the first counted subject, can take
// `msa_index` past the rows NCBI allocates (one per counted subject). NCBI then writes out of
// bounds (`_ASSERT` only), a failure without a checkable result: LOSAT refuses the PSSM there.
// The hit lists are sorted by their best e-value (`s_EvalueCompareHSPLists`), which puts the
// counted subjects first, but two e-values below 1.0e-180 compare as equal and those lists are
// then ordered by score, so with an inclusion threshold below 1.0e-180 the counted subjects are
// not shown to come first:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1390-1395 and 3095-3102
// ```c
// s_EvalueComp(double evalue1, double evalue2)
// {
//     const double epsilon = 1.0e-180;
//     if (evalue1 < epsilon && evalue2 < epsilon) {
//         return 0;
//     }
//     ...
//    if ((retval = s_EvalueComp(h1->best_evalue,
//                                    h2->best_evalue)) != 0)
//       return retval;
//
//    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score)
//       return -1;
//    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score)
//       return 1;
// ```
fn input_msa(
    query: &[u8],
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Vec<u8>,
    inclusion_ethresh: f64,
) -> Result<PsiMsa> {
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
            let Some(row) = msa.rows.get_mut(msa_index) else {
                anyhow::bail!(
                    "an alignment past the last row of the PSSM's multiple alignment (a subject included by -inclusion_ethresh after one that is not, where NCBI BLAST+'s psiblast writes out of bounds) is not supported by LOSAT's PSIBLAST"
                );
            };
            process_denseg(row, alignment, &subject(alignment.subject));
        }
        last_subject = Some(alignment.subject);
    }
    Ok(msa)
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

/// A subject record's letters in NCBIstdaa as its Bioseq holds them (the `CSeqVector` of
/// `CPsiBlastInputData::x_GetSubjectSequence`, quoted at `process_denseg`): upper case, and O
/// stays O. The query's letters are not these: the PSSM's query is the search's buffer, in which
/// O is X (`query_pssm`).
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
// `x_PSIMatrix2Asn1` stores every row). The query's letters are those of the PSSM's query Bioseq
// (`GetQuerySequenceData`, blastfmtutil.cpp:253-255), which `CPsiBlastInputData` builds from the
// search's buffer, in which an O is X (quoted at `query_pssm`): `query` is that buffer, so an O
// of the record prints as X and takes the X branch of the gapless weight. The diagnostics it
// reads are those that `-out_ascii_pssm` asks for (`PssmRequest::diagnostics`).
pub(crate) fn print_ascii_pssm(pssm: &PsiPssm, query: &[u8], ancillary: PssmAncillary) -> Vec<u8> {
    use std::fmt::Write as _;
    const REQUESTED: &str = "a diagnostic that -out_ascii_pssm requests";
    let diagnostics = &pssm.diagnostics;
    let weighted_residue_freqs = diagnostics
        .weighted_residue_freqs
        .as_ref()
        .expect(REQUESTED);
    let information_content = diagnostics.information_content.as_ref().expect(REQUESTED);
    let num_matching_seqs = diagnostics.num_matching_seqs.as_ref().expect(REQUESTED);
    let gapless_column_weights = diagnostics
        .gapless_column_weights
        .as_ref()
        .expect(REQUESTED);
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
                value *= weighted_residue_freqs[i][c];
                value = (value + if value > 0.0 { 0.5 } else { -0.5 }) as i32 as f64;
                let _ = write!(out, "{:>4}", value as i32);
            }
        }
        let _ = write!(out, "  {} ", fixed(information_content[i], 2));
        if num_matching_seqs[i] > 1 && residue != ncbistdaa::X {
            out.push_str(&fixed(gapless_column_weights[i], 2));
        } else {
            out.push_str("    0.00");
        }
    }
    out.push_str("\n\n");
    out.push_str("                      K         Lambda\n");
    for (label, block) in [
        ("Standard Ungapped    ", ancillary.ungapped),
        ("Standard Gapped      ", ancillary.gapped),
        ("PSI Ungapped         ", ancillary.psi_ungapped),
        ("PSI Gapped           ", ancillary.psi_gapped),
    ] {
        let _ = writeln!(
            out,
            "{label}{}     {}",
            fixed(block.k, 4),
            fixed(block.lambda, 4)
        );
    }
    out.into_bytes()
}

/// With `-out_pssm`, the refusal of a query whose title has a byte outside the visible ASCII
/// range (0x20 to 0x7E), before the search: the PSSM's Bioseq holds the title as a
/// VisibleString, and NCBI's ASN.1 text stream writes `#` for each such byte with an error message
/// of three lines for each (`Error: (803.7) [psiblast] PssmWithParameters.pssm.query.seq.descr.E.title`,
/// `Bad char [0xC3] in string at line 12`, the title; NCBI 2.17.0, a UTF-8 title, status 0),
/// which LOSAT does not reproduce (as outfmt 11's titles, `report::archive::check_titles`). The
/// FASTA reader cuts a title at its first byte below 0x20, so the bytes refused are those from
/// 0x7F on.
// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
// ```c++
// void CObjectOStreamAsn::WriteString(const char* ptr, size_t length)
// {
//     m_Output.PutChar('"');
//     CTempString original(ptr, length);
//     size_t valid=0;
//     while ( length > 0 ) {
//         char c = *ptr;
//         if ( x_FixCharsMethod() != eFNP_Allow ) {
//             if ( !GoodVisibleChar(c) ) {
//                 if (valid == 0) {
// #if SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_WRITING
//                     valid = CUtf8::EvaluateSymbolLength(CTempString(ptr,length));
// #endif
//                     if (valid == 0) {
//                         c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, original, x_FixCharsSubst());
//                     }
// ```
// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:456-459
// ```c++
// bool GoodVisibleChar(char c)
// {
//     return c >= ' ' && c <= '~';
// }
// ```
pub(crate) fn check_pssm_titles(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
) -> Result<()> {
    if let Some(record) = records.iter().find(|record| {
        record
            .title
            .iter()
            .any(|byte| !(0x20..=0x7E).contains(byte))
    }) {
        anyhow::bail!(
            "the title of the query {} has a character the PSSM file cannot write as NCBI does (outside the visible ASCII range), which is not supported by LOSAT's PSIBLAST with -out_pssm",
            record.local_id
        );
    }
    Ok(())
}

/// `CAutoOutputFileReset`: the file of one PSSM option, opened again at each write
/// (`GetStream`). Without `-save_each_pssm` an existing file of the name is deleted and the file
/// made again, so the last PSSM stays; with it the writes are `F.1`, `F.2`, ...
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
#[derive(Debug)]
pub(crate) struct OutputFileReset {
    file_name: String,
    version: u32,
}

impl OutputFileReset {
    pub(crate) fn new(file_name: &str, use_versions: bool) -> Self {
        Self {
            file_name: file_name.to_string(),
            version: u32::from(use_versions),
        }
    }

    /// `GetStream`'s file: `F`, or `F.<version>` (the version counts the writes).
    pub(crate) fn next_file(&mut self) -> PssmFile {
        if self.version == 0 {
            return PssmFile {
                path: std::path::PathBuf::from(&self.file_name),
                replaced: true,
            };
        }
        let path = std::path::PathBuf::from(format!("{}.{}", self.file_name, self.version));
        self.version += 1;
        PssmFile {
            path,
            replaced: false,
        }
    }
}

/// One file of `CAutoOutputFileReset::GetStream`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PssmFile {
    pub(crate) path: std::path::PathBuf,
    /// Version 0 (no `-save_each_pssm`): an existing regular file of the name is removed before
    /// the file is opened, so the file is made again; a numbered file is only truncated, which
    /// keeps an existing file's inode and mode.
    pub(crate) replaced: bool,
}

impl PssmFile {
    /// `GetStream`: with version 0, `CFile::Exists` (a regular file, through symbolic links) and
    /// `Remove` (the name itself, `remove(3)`; a failure is neither checked nor, by default,
    /// logged), then `ofstream` (made or truncated, through symbolic links).
    // NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:4039-4042
    // ```c++
    // bool CFile::Exists(void) const
    // {
    //     return IsFile();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2595-2610
    // ```c++
    // bool CDirEntry::Remove(TRemoveFlags flags) const
    // {
    //     // Is this a directory ? (and processing not entry only)
    //     if ( (flags & (fDir_All | fDir_Recursive)) != eEntryOnly  &&  IsDir(eIgnoreLinks) ) {
    //         return CDir(GetPath()).Remove(flags);
    //     }
    //     // Other entries
    //     return RemoveEntry(flags);
    // }
    //
    //
    // bool CDirEntry::RemoveEntry(TRemoveFlags flags) const
    // {
    //
    //     if ( NcbiSys_remove(_T_XCSTRING(GetPath())) != 0 ) {
    // ```
    pub(crate) fn open(&self) -> std::io::Result<std::fs::File> {
        if self.replaced && self.path.is_file() {
            let _ = std::fs::remove_file(&self.path);
        }
        std::fs::File::create(&self.path)
    }
}

/// What a PSSM file's name is before the run, for NCBI's `GetStream` on it (`PssmFile::open`).
#[derive(Debug)]
enum TargetKind {
    /// No file: `ofstream` makes it in its directory.
    Absent,
    /// A regular file, also through a symbolic link.
    Regular,
    /// A character device (`/dev/null`, `/dev/full`): opened for writing as it is; `null`: it is
    /// `/dev/null`.
    #[cfg_attr(
        not(unix),
        expect(dead_code, reason = "constructed only on unix (special_target_kind)")
    )]
    Device { null: bool },
    /// A directory or a socket: `ofstream` cannot open it.
    Unopenable,
    /// A named pipe, a block device or a symbolic link to no file: not supported.
    Unsupported(&'static str),
    /// A name whose directories cannot be read: `ofstream` fails the same way.
    Unreachable(std::io::Error),
}

fn target_kind(path: &std::path::Path) -> TargetKind {
    match std::fs::metadata(path) {
        Ok(metadata) if metadata.is_file() => TargetKind::Regular,
        Ok(metadata) => special_target_kind(&metadata),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            if std::fs::symlink_metadata(path).is_ok() {
                TargetKind::Unsupported("a symbolic link to a missing file")
            } else {
                TargetKind::Absent
            }
        }
        Err(error) => TargetKind::Unreachable(error),
    }
}

#[cfg(unix)]
fn special_target_kind(metadata: &std::fs::Metadata) -> TargetKind {
    use std::os::unix::fs::FileTypeExt;
    let file_type = metadata.file_type();
    if file_type.is_fifo() {
        TargetKind::Unsupported("a named pipe")
    } else if file_type.is_block_device() {
        TargetKind::Unsupported("a block device")
    } else if file_type.is_char_device() {
        TargetKind::Device {
            null: device_number(metadata).is_some_and(|rdev| Some(rdev) == null_device()),
        }
    } else {
        TargetKind::Unopenable
    }
}

#[cfg(not(unix))]
fn special_target_kind(_metadata: &std::fs::Metadata) -> TargetKind {
    TargetKind::Unopenable
}

/// The device number of a character device.
#[cfg(unix)]
fn device_number(metadata: &std::fs::Metadata) -> Option<u64> {
    use std::os::unix::fs::{FileTypeExt, MetadataExt};
    metadata
        .file_type()
        .is_char_device()
        .then(|| metadata.rdev())
}

/// `/dev/null`'s device number.
#[cfg(unix)]
fn null_device() -> Option<u64> {
    std::fs::metadata("/dev/null")
        .ok()
        .as_ref()
        .and_then(device_number)
}

/// Whether the test open of an existing file failed with `EPERM`: an immutable or append-only file
/// (`open(2)`), whose name cannot be removed either (`unlink(2)`, `EPERM`).
#[cfg(unix)]
fn immutable(error: &std::io::Error) -> bool {
    rustix::io::Errno::from_io_error(error) == Some(rustix::io::Errno::PERM)
}

#[cfg(not(unix))]
fn immutable(_error: &std::io::Error) -> bool {
    false
}

/// Whether removing `path` succeeds, given that its directory takes a new file (the staging file
/// made there by the user `uid`): in a directory with the sticky bit only the owner of the file or
/// of the directory, or root, removes a name (`unlink(2)`, `EPERM`).
#[cfg(unix)]
fn removable(path: &std::path::Path, directory: &std::path::Path, uid: u32) -> bool {
    use std::os::unix::fs::MetadataExt;
    let (Ok(entry), Ok(directory)) = (
        std::fs::symlink_metadata(path),
        std::fs::metadata(directory),
    ) else {
        return false;
    };
    directory.mode() & 0o1000 == 0 || uid == 0 || entry.uid() == uid || directory.uid() == uid
}

/// Elsewhere a file that cannot be opened for writing is taken as one that cannot be removed
/// either (a refusal where NCBI might remove it).
#[cfg(not(unix))]
fn removable(_path: &std::path::Path, _directory: &std::path::Path, _uid: u32) -> bool {
    false
}

#[cfg(unix)]
fn owner(metadata: &std::fs::Metadata) -> u32 {
    std::os::unix::fs::MetadataExt::uid(metadata)
}

#[cfg(not(unix))]
fn owner(_metadata: &std::fs::Metadata) -> u32 {
    0
}

/// The PSSM files of a psiblast run (`-out_pssm`, `-out_ascii_pssm`): made with the PSSMs, before
/// the reports (`add`), and written after them (`commit`), each query's ASN.1 file then its ASCII
/// one, as psiblast's `SavePssmToFile` after the query's report. The ASCII PSSMs, whose writes
/// NCBI does not check, are kept in memory.
///
/// NCBI ends a run whose ASN.1 file cannot be opened or written after that query's report, so LOSAT
/// must know before its first report byte whether each file can be written; and a run that fails
/// before its reports are out (a report write error, a LOSAT refusal) must leave every file as it
/// was, as NCBI's, which stops before `SavePssmToFile`. So `add` writes the ASN.1 bytes to a
/// staging file beside the target (`.losat-<pid>-<n>`, a new file in the target's directory, whose
/// name is short so that any name the target can have is no reason for its failure), which shows
/// that the directory takes a file and the file system the bytes, and foresees NCBI's open of the
/// target without changing it; only `commit` writes the targets.
/// Dropped without `commit`, `PssmFiles` removes the staging files: no target has changed. A
/// failure in `add` is the refusal of a file that NCBI cannot write either (its stop after that
/// query's report, `PssmStop`) or that LOSAT cannot write as NCBI does; no partial file is left.
///
/// The rule for a target, from what `GetStream` does with it (`PssmFile::open`):
/// - no file of the name: made in its directory, which takes the staging file;
/// - a regular file (also through a symbolic link): removed and made again (version 0) or
///   truncated in place. LOSAT opens it for writing (without truncation); if that fails with
///   `EPERM` (an immutable or append-only file, whose name cannot be removed either), or fails and
///   the name cannot be removed (`removable`) or the file is numbered, NCBI's open fails there.
///   An existing file that opens for writing but beside which no staging file can be made (a
///   directory that takes no new file) is written in place, as NCBI's version 0, which cannot remove
///   its name, and a numbered file do: its bytes are kept in memory for `commit` (only the last
///   query's with version 0), which cannot show before the reports that they fit;
/// - a character device: opened for writing and written zero bytes, which a device that takes no
///   bytes refuses (`/dev/full`: ENOSPC), so the run stops at that query as NCBI's. `/dev/null`
///   gets nothing more; the PSSM of every query for another device is kept in memory (a few MB per
///   query of some thousand residues) and written by `commit`. A device that takes zero bytes but
///   not the PSSM's would fail only after the reports;
/// - a directory or a socket: NCBI's open fails, and so does LOSAT's;
/// - a named pipe, a block device, a symbolic link to no file: refused (also before the search,
///   `check_pssm_files`).
///
/// Not foreseen: a file system that fills between `add` and `commit`, or that has no room for the
/// bytes of a file written in place (LOSAT's peak use is the old file, at most two staged PSSMs of
/// version 0 and every numbered one); such a failure comes after the reports (`commit`).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-407
// ```c++
// CPsiBlastApp::SavePssmToFile(CRef<CPssmWithParameters> pssm,
//                              const CPsiBlastIterationState* itr)
// {
//     if (pssm.Empty()) {
//         return;
//     }
//
//     if (m_CmdLineArgs->SaveCheckpoint() &&
//         (itr == NULL || // this is true in the case of remote PSI-BLAST
//          itr->GetIterationNumber() >= 1)) {
//         *m_CmdLineArgs->GetCheckpointStream() << MSerial_AsnText << *pssm;
//     }
//
//     if (m_CmdLineArgs->SaveAsciiPssm() &&
//         (itr == NULL || // this is true in the case of remote PSI-BLAST
//          itr->GetIterationNumber() >= 1)) {
//         ...
//         CBlastFormatUtil::PrintAsciiPssm(*pssm,
//                                          m_AncillaryData,
//                                          *m_CmdLineArgs->GetAsciiPssmStream());
//     }
// ```
#[derive(Debug)]
pub(crate) struct PssmFiles {
    asn: Option<OutputFileReset>,
    ascii: Option<OutputFileReset>,
    pending: Vec<PendingPssm>,
    /// The entry whose staging file holds the last version-0 ASN.1 PSSM, which the next one
    /// replaces: of NCBI's writes of one name only the last one's bytes stay, made or truncated
    /// as each of them would make or truncate it.
    last_replaced: Option<usize>,
    staging_count: u64,
}

/// One query's PSSM files, to be written after the reports.
#[derive(Debug)]
struct PendingPssm {
    asn: Option<(PssmFile, Staged)>,
    ascii: Option<Vec<u8>>,
}

#[derive(Debug)]
enum Staged {
    /// The bytes in a staging file beside the target.
    File(std::path::PathBuf),
    /// The bytes of an existing file in a directory that takes no new file, written in place.
    InPlace(Vec<u8>),
    /// The bytes for a character device.
    Device(Vec<u8>),
    /// Nothing: the target is `/dev/null`.
    Null,
    /// Written, or replaced by a later query's file of the same name.
    Done,
}

impl PssmFiles {
    pub(crate) fn new(request: &PssmRequest) -> Self {
        let stream =
            |path: Option<&str>| path.map(|path| OutputFileReset::new(path, request.save_each));
        Self {
            asn: stream(request.asn_path.as_deref()),
            ascii: stream(request.ascii_path.as_deref()),
            pending: Vec::new(),
            last_replaced: None,
            staging_count: 0,
        }
    }

    /// One query's PSSM: its ASN.1 file staged and its ASCII PSSM kept. `Err` is LOSAT's refusal
    /// of an ASN.1 file that cannot be written as NCBI does (NCBI stops after that query's report).
    // NCBI's stream writes through a 4 KiB buffer (`COStreamBuffer`); the first write that fails
    // throws a `CIOException` (`eWrite`, "write fault"), the object stream reports it once as an
    // error of its own (`SetFailFlags`: `Error: (802.5) [psiblast] CObjectOStream: error at byte
    // N: stack is empty: write fault`, N being where the buffer first overflowed), and psiblast's
    // `CATCH_ALL` ends the run with status 0, without the epilog and the later queries. LOSAT
    // stops at the same place, with its own refusal instead of NCBI's message and status, and
    // without NCBI's partial file.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-904 and 921-924
    // ```c++
    // void COStreamBuffer::FlushBuffer(bool fullBuffer)
    //     THROWS1((CIOException))
    // {
    //     ...
    //     if ( count != 0 ) {
    //         if ( !m_Output.write(m_Buffer, count) ) {
    //             m_Error = "write fault";
    //             NCBI_THROW(CIOException,eWrite,m_Error);
    //         }
    //     ...
    // char* COStreamBuffer::DoReserve(size_t count)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     FlushBuffer(false);
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384
    // ```c++
    // CObjectOStream::TFailFlags CObjectOStream::SetFailFlags(TFailFlags flags,
    //                                                         const char* message)
    // {
    //     TFailFlags old = m_Fail;
    //     m_Fail |= flags;
    //     if ( !old && flags ) {
    //         // first fail
    //         ERR_POST_X(5, "CObjectOStream: error at "<<
    //                    GetPosition()<<": "<<GetStackTrace() << ": " << message);
    //     }
    //     return old;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    // ```c++
    //     catch (const CIOException& e) {                                         \
    //         if (e.GetErrCode() == CIOException::eFlush) {                       \
    //             LOG_POST(Error << "BLAST failed to write output: " << e.GetMsg());\
    //             exit_code = BLAST_OUTPUT_ERROR;                                 \
    //         }                                                                   \
    //     }                                                                       \
    // ```
    pub(crate) fn add(&mut self, asn: Option<&[u8]>, ascii: Option<Vec<u8>>) -> Result<()> {
        let file = asn.and(self.asn.as_mut()).map(OutputFileReset::next_file);
        let asn = match (file, asn) {
            (Some(file), Some(bytes)) => {
                let staged = self.stage(&file, bytes)?;
                Some((file, staged))
            }
            _ => None,
        };
        if let Some((file, Staged::File(_) | Staged::InPlace(_))) = &asn {
            if file.replaced {
                if let Some(previous) = self.last_replaced.replace(self.pending.len()) {
                    if let Some((_, staged)) = self.pending[previous].asn.as_mut() {
                        if let Staged::File(staging) = std::mem::replace(staged, Staged::Done) {
                            let _ = std::fs::remove_file(staging);
                        }
                    }
                }
            }
        }
        self.pending.push(PendingPssm { asn, ascii });
        Ok(())
    }

    /// The ASN.1 bytes of `file`, staged (`PssmFiles`).
    fn stage(&mut self, file: &PssmFile, asn: &[u8]) -> Result<Staged> {
        let path = &file.path;
        let cannot_write = |error: std::io::Error| {
            anyhow::anyhow!(
                "a -out_pssm file that cannot be written ({}: {error}), where NCBI BLAST+'s psiblast stops after the query's report with a serial-stream write error, is not supported by LOSAT's PSIBLAST",
                path.display()
            )
        };
        let existing = match target_kind(path) {
            TargetKind::Absent => false,
            TargetKind::Regular => true,
            TargetKind::Device { null } => {
                return test_device(path, null, asn).map_err(cannot_write)
            }
            TargetKind::Unopenable => return test_device(path, false, asn).map_err(cannot_write),
            TargetKind::Unsupported(kind) => return Err(unsupported_target(path, kind)),
            TargetKind::Unreachable(error) => return Err(cannot_write(error)),
        };
        if path.file_name().is_none() {
            return Err(cannot_write(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "not a file name",
            )));
        }
        let directory = match path.parent() {
            Some(parent) if !parent.as_os_str().is_empty() => parent,
            _ => std::path::Path::new("."),
        };
        let (staging, uid) = match self.write_staging(directory, asn) {
            Ok(staged) => staged,
            Err(_) if existing => {
                // Where the directory takes no new file, NCBI's version 0 cannot remove the file
                // either and opens it in place, as a numbered file.
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(path)
                    .map_err(cannot_write)?;
                return Ok(Staged::InPlace(asn.to_vec()));
            }
            Err(error) => return Err(cannot_write(error)),
        };
        if existing {
            if let Err(error) = std::fs::OpenOptions::new().write(true).open(path) {
                if immutable(&error) || !(file.replaced && removable(path, directory, uid)) {
                    let _ = std::fs::remove_file(&staging);
                    return Err(cannot_write(error));
                }
            }
        }
        Ok(Staged::File(staging))
    }

    /// Writes `bytes` to a new staging file `.losat-<pid>-<n>` in the target's `directory`; its
    /// path and owner.
    fn write_staging(
        &mut self,
        directory: &std::path::Path,
        bytes: &[u8],
    ) -> std::io::Result<(std::path::PathBuf, u32)> {
        loop {
            self.staging_count += 1;
            let staging = directory.join(format!(
                ".losat-{}-{}",
                std::process::id(),
                self.staging_count
            ));
            let mut file = match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&staging)
            {
                Ok(file) => file,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            };
            let written = file
                .write_all(bytes)
                .and_then(|()| file.metadata())
                .map(|metadata| owner(&metadata));
            return match written {
                Ok(uid) => Ok((staging, uid)),
                Err(error) => {
                    drop(file);
                    let _ = std::fs::remove_file(&staging);
                    Err(error)
                }
            };
        }
    }

    /// After the reports: each query's ASN.1 file (`GetStream` on the target, then the staged
    /// bytes), then its ASCII file (not checked, as NCBI's `ofstream`). `Err`: an ASN.1 file that
    /// `add` found writable failed, after the reports (see `PssmFiles`).
    pub(crate) fn commit(mut self) -> Result<()> {
        for index in 0..self.pending.len() {
            if let Some((file, staged)) = self.pending[index].asn.as_mut() {
                let bytes = match std::mem::replace(staged, Staged::Done) {
                    Staged::File(staging) => {
                        let bytes = std::fs::read(&staging);
                        let _ = std::fs::remove_file(&staging);
                        Some(bytes)
                    }
                    Staged::InPlace(bytes) | Staged::Device(bytes) => Some(Ok(bytes)),
                    Staged::Null | Staged::Done => None,
                };
                if let Some(bytes) = bytes {
                    bytes
                        .and_then(|bytes| file.open()?.write_all(&bytes))
                        .map_err(|error| {
                            anyhow::anyhow!(
                                "a -out_pssm file that could not be written after the reports ({}: {error}), although LOSAT found it writable before them, is not supported by LOSAT's PSIBLAST",
                                file.path.display()
                            )
                        })?;
                }
            }
            if let (Some(stream), Some(ascii)) =
                (self.ascii.as_mut(), self.pending[index].ascii.take())
            {
                if let Ok(mut out) = stream.next_file().open() {
                    let _ = out.write_all(&ascii);
                }
            }
        }
        Ok(())
    }
}

impl Drop for PssmFiles {
    // Out of line: it runs where the blastp engine's function that holds the files
    // (`run_resolved_in_pool`) returns.
    #[cold]
    #[inline(never)]
    fn drop(&mut self) {
        for pending in &self.pending {
            if let Some((_, Staged::File(staging))) = &pending.asn {
                let _ = std::fs::remove_file(staging);
            }
        }
    }
}

/// A character device's test (`PssmFiles`): opened for writing and written zero bytes; what
/// `commit` writes to it.
fn test_device(path: &std::path::Path, null: bool, asn: &[u8]) -> std::io::Result<Staged> {
    let mut device = std::fs::OpenOptions::new().write(true).open(path)?;
    // A write of zero bytes is passed to the device (`write(2)`).
    let written = device.write(&[])?;
    debug_assert_eq!(written, 0);
    Ok(if null {
        Staged::Null
    } else {
        Staged::Device(asn.to_vec())
    })
}

fn unsupported_target(path: &std::path::Path, kind: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "an -out_pssm file that is {kind} ({}) is not supported by LOSAT's PSIBLAST",
        path.display()
    )
}

/// Before the search: refuses PSSM files that LOSAT cannot write as NCBI does. NCBI writes each
/// query's files after its report through streams of their own, so a PSSM file that is also the
/// standard output, the standard error or the `-out` file gets its bytes among the reports' in an
/// order that their buffers make, which LOSAT does not reproduce (it writes the PSSM files after
/// all reports). Every name the run can write is checked: `F`, or with `-save_each_pssm` `F.1` to
/// `F.<number of queries>`; a name that does not exist cannot be one of those streams. The same
/// file is the same device and inode (`metadata` follows `/dev/stdout` to the stream's file),
/// or for a character device the same device number (`/dev/null` excepted: no order of its bytes
/// shows); without inode numbers, the same canonical path as `-out`, or one of the names of the
/// standard streams (`ReportStreams`, taken once). Also refused, for `-out_pssm`: a named pipe, a
/// block device and a symbolic link to no file (`PssmFiles`).
pub(crate) fn check_pssm_files(
    request: &PssmRequest,
    queries: usize,
    out: Option<&std::path::Path>,
) -> Result<()> {
    let check = PssmFileCheck::new(out);
    for (option, path) in [
        ("-out_pssm", request.asn_path.as_deref()),
        ("-out_ascii_pssm", request.ascii_path.as_deref()),
    ] {
        let Some(path) = path else {
            continue;
        };
        let names: Vec<std::path::PathBuf> = if request.save_each {
            (1..=queries)
                .map(|n| std::path::PathBuf::from(format!("{path}.{n}")))
                .collect()
        } else {
            vec![std::path::PathBuf::from(path)]
        };
        for name in &names {
            check.check(option, name)?;
        }
    }
    Ok(())
}

/// `check_pssm_files`' check of one name, also made by the round driver before each write
/// (`rounds::RoundFiles`), whose numbered names are not known before the search.
pub(crate) struct PssmFileCheck {
    streams: ReportStreams,
}

impl PssmFileCheck {
    pub(crate) fn new(out: Option<&std::path::Path>) -> Self {
        Self {
            streams: ReportStreams::new(out),
        }
    }

    /// The refusal of the file `name` of `option` (`-out_pssm` or `-out_ascii_pssm`) where it
    /// is a report stream, or for `-out_pssm` a named pipe, a block device or a symbolic link to
    /// no file.
    pub(crate) fn check(&self, option: &str, name: &std::path::Path) -> Result<()> {
        if let Some(stream) = self.streams.stream_of(name) {
            anyhow::bail!(
                "an {option} file that is also {stream} ({}) is not supported by LOSAT's PSIBLAST",
                name.display()
            );
        }
        if option == "-out_pssm" {
            if let TargetKind::Unsupported(kind) = target_kind(name) {
                return Err(unsupported_target(name, kind));
            }
        }
        Ok(())
    }
}

/// Whether `path` names a regular file (also through a symbolic link), which a failed write
/// leaves partly written (`rounds::RoundFiles`).
pub(crate) fn is_regular_file(path: &std::path::Path) -> bool {
    matches!(target_kind(path), TargetKind::Regular)
}

/// The report streams (`check_pssm_files`): the standard output, the standard error (`fstat` of a
/// copy of each descriptor) and the `-out` file, as files.
#[cfg(unix)]
struct ReportStreams {
    streams: Vec<(&'static str, std::fs::Metadata)>,
    null: Option<u64>,
}

#[cfg(unix)]
impl ReportStreams {
    fn new(out: Option<&std::path::Path>) -> Self {
        use std::os::fd::AsFd;
        let stream = |fd: std::os::fd::BorrowedFd<'_>| {
            fd.try_clone_to_owned()
                .ok()
                .and_then(|fd| std::fs::File::from(fd).metadata().ok())
        };
        let streams = [
            ("the standard output", stream(std::io::stdout().as_fd())),
            ("the standard error", stream(std::io::stderr().as_fd())),
            (
                "the -out file",
                out.and_then(|out| std::fs::metadata(out).ok()),
            ),
        ]
        .into_iter()
        .filter_map(|(name, metadata)| Some((name, metadata?)))
        .collect();
        Self {
            streams,
            null: null_device(),
        }
    }

    /// The report stream that the file `name` is; `/dev/null`, where no order of bytes shows, is
    /// none.
    fn stream_of(&self, name: &std::path::Path) -> Option<&'static str> {
        use std::os::unix::fs::MetadataExt;
        let target = std::fs::metadata(name).ok()?;
        let device = device_number(&target);
        if device.is_some() && device == self.null {
            return None;
        }
        self.streams.iter().find_map(|(stream, metadata)| {
            ((metadata.dev(), metadata.ino()) == (target.dev(), target.ino())
                || device.is_some_and(|rdev| device_number(metadata) == Some(rdev)))
            .then_some(*stream)
        })
    }
}

#[cfg(not(unix))]
struct ReportStreams {
    out: Option<std::path::PathBuf>,
}

#[cfg(not(unix))]
impl ReportStreams {
    fn new(out: Option<&std::path::Path>) -> Self {
        Self {
            out: out.and_then(|out| std::fs::canonicalize(out).ok()),
        }
    }

    fn stream_of(&self, name: &std::path::Path) -> Option<&'static str> {
        match name.to_str() {
            Some("/dev/stdout" | "/dev/fd/1" | "/proc/self/fd/1") => {
                return Some("the standard output")
            }
            Some("/dev/stderr" | "/dev/fd/2" | "/proc/self/fd/2") => {
                return Some("the standard error")
            }
            _ => {}
        }
        let target = std::fs::canonicalize(name).ok()?;
        (self.out.as_ref() == Some(&target)).then_some("the -out file")
    }
}

#[cfg(test)]
#[path = "pssm_tests.rs"]
mod tests;
