//! psiblast `-in_pssm` with rounds or with the PSSM file options (series PS, item PS5b): the input
//! PSSM as NCBI's `CPssmWithParameters` holds it when round 1 starts (after `CPsiBlastImpl`'s
//! constructor, which computes the scores of a PSSM of frequency ratios only), which
//! `SavePssmToFile` writes then; the same object after round 1's K/A write-back, which it writes
//! again after a round 1 without alignments with `-save_pssm_after_last_round`; the ASCII PSSM of a
//! file PSSM (`ExtractPssmAncillaryData`: the PSI blocks only, the short form or, with
//! `informationContent`, the long one); the Seq-id that the PSSMs of the later rounds keep; and
//! LOSAT's refusals of what NCBI would write back that LOSAT's reader does not keep, or that NCBI
//! cannot write.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-408
//! ```c++
//! CPsiBlastApp::SavePssmToFile(CRef<CPssmWithParameters> pssm,
//!                              const CPsiBlastIterationState* itr)
//! {
//!     if (pssm.Empty()) {
//!         return;
//!     }
//!
//!     if (m_CmdLineArgs->SaveCheckpoint() &&
//!         (itr == NULL || // this is true in the case of remote PSI-BLAST
//!          itr->GetIterationNumber() >= 1)) {
//!         *m_CmdLineArgs->GetCheckpointStream() << MSerial_AsnText << *pssm;
//!     }
//!
//!     if (m_CmdLineArgs->SaveAsciiPssm() &&
//!         (itr == NULL || // this is true in the case of remote PSI-BLAST
//!          itr->GetIterationNumber() >= 1)) {
//!         if (m_AncillaryData.Empty() && pssm.NotEmpty()) {
//!             _ASSERT(itr->GetIterationNumber() == 1);
//!             m_AncillaryData = ExtractPssmAncillaryData(*pssm);
//!         }
//!         CBlastFormatUtil::PrintAsciiPssm(*pssm,
//!                                          m_AncillaryData,
//!                                          *m_CmdLineArgs->GetAsciiPssmStream());
//!     }
//! }
//! ```
//! With `-in_pssm` the loop's first `SavePssmToFile(pssm, &itr)` writes the input PSSM, the object
//! that `CPsiBlast` searches (`CPsiBlast(pssm, db_adapter, psi_opts)` holds the same `CRef`); its
//! ASCII PSSM takes the ancillary data of the PSSM itself, which stays for the files written after
//! the loop until the next PSSM is computed (`ComputePssmForNextIteration` replaces
//! `m_AncillaryData`).

use anyhow::{bail, Result};

use crate::core::blast_util::{BLASTAA_SIZE, BLAST_SCORE_MIN};
use crate::core::psi::KarlinBlk;
use crate::report::serial::pssm::{PssmWithParameters, SeqData};
use crate::seqid::ParsedSeqId;

use super::pssm::{AsciiPssmDiagnostics, PssmRequest};
#[cfg(test)]
use crate::algorithm::pssm_search::in_pssm::with_computed_scores;
use crate::algorithm::pssm_search::in_pssm::InPssmFile;

/// The query PSSM of psiblast's rounds with `-in_pssm` (`rounds::PsiRoundPlan::input`): the
/// first Seq-id of its query, which the PSSM of each later round keeps (the first Seq-align's query
/// id: psiblast_app.cpp:539 makes the query's Seq-loc of `GetFirstId()`),
/// and, when a PSSM file is written, the object that `SavePssmToFile` writes.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct InPssmRounds {
    pub(crate) id: ParsedSeqId,
    object: Option<PssmWithParameters>,
}

/// The input PSSM of a run of psiblast's rounds or PSSM files with `-in_pssm` (`request`: the
/// files), or LOSAT's refusal before the search of a file NCBI would write back in a way LOSAT
/// does not:
/// - with `-out_pssm`, a member that the reader skipped (NCBI's object keeps it and writes it:
///   the descriptors of the query other than titles, its `annot`, the members of its `Seq-inst`
///   other than `repr`, `mol`, `length` and `seq-data`, the `constraints` and `annotatedSites` of
///   the parameters; NCBI 2.17.0 probes of PS5b), and a PSSM of frequency ratios only without
///   `rpsdbparams` (NCBI writes part of the file, then stops with a C++ exception: "cannot write
///   unassigned member matrixName", rc 255);
/// - with `-out_ascii_pssm`, a PSSM with `informationContent` (the long form) whose
///   `weightedResFreqsPerPos`, `informationContent`, `numMatchingSeqs` or the
///   `gaplessColumnWeights` it reads do not cover every column (`PrintAsciiPssm` dereferences
///   the missing matrix or reads past the vectors: undefined behaviour).
///
/// A PSSM that `CPsiBlastValidate::Pssm` refuses gets NCBI's error when the search starts, before
/// any file (the batch stop, `blast_engine::psiblast_stop_at_invalid_pssm`).
pub(crate) fn rounds_input(file: InPssmFile, request: &PssmRequest) -> Result<InPssmRounds> {
    let id = match file
        .object
        .pssm
        .query
        .as_ref()
        .map(|query| query.ids.as_slice())
    {
        Some([id, ..]) => id.clone(),
        // `in_pssm::psi_search_and_object` has refused a query without a Seq-id.
        _ => bail!("a PSSM query without a Seq-id is not supported by LOSAT's PSIBLAST"),
    };
    if !file.valid {
        return Ok(InPssmRounds { id, object: None });
    }
    if request.asn_path.is_some() {
        if !file.skipped.is_empty() {
            bail!(
                "a PSSM file with members that LOSAT's reader does not keep ({}), which NCBI BLAST+ writes back with -out_pssm, is not supported by LOSAT's PSIBLAST",
                file.skipped.join(", ")
            );
        }
        let rpsdbparams = file
            .object
            .params
            .as_ref()
            .is_some_and(|params| params.rpsdbparams.is_some());
        if file.computed && !rpsdbparams {
            bail!(
                "a PSSM of frequency ratios only without params.rpsdbparams with -out_pssm (NCBI BLAST+ cannot write the matrix name it lacks and stops with a C++ exception) is not supported by LOSAT's PSIBLAST"
            );
        }
    }
    if request.ascii_path.is_some() {
        check_long_ascii(&file.object)?;
    }
    let files = request.asn_path.is_some() || request.ascii_path.is_some();
    Ok(InPssmRounds {
        id,
        object: files.then_some(file.object),
    })
}

/// The refusal of the long ASCII form of a PSSM file that NCBI cannot print (`rounds_input`).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:186-188, 278-301 (quoted
// at `pssm::print_ascii_pssm`: the long form reads `(*weighted_res_freqs)(r, i)`,
// `info_content[i]`, `num_matching_seqs[i]` and, where that is above 1 and the letter is not X,
// `gapless_col_weights[i]`; `GetWeightedResidueFrequencies` gives NULL without the member)
fn check_long_ascii(object: &PssmWithParameters) -> Result<()> {
    let pssm = &object.pssm;
    let Some(data) = pssm.intermediate_data.as_ref() else {
        return Ok(());
    };
    let info = data.information_content.as_deref().unwrap_or(&[]);
    if info.is_empty() {
        return Ok(());
    }
    let columns = usize::try_from(pssm.num_columns).unwrap_or(0);
    let rows_n = usize::try_from(pssm.num_rows).unwrap_or(0);
    let cells = columns * rows_n;
    let weighted = data.weighted_res_freqs_per_pos.as_deref().unwrap_or(&[]);
    let matching = data.num_matching_seqs.as_deref().unwrap_or(&[]);
    let gapless = data.gapless_column_weights.as_deref().unwrap_or(&[]);
    let letters = query_codes(object);
    let gapless_short = (gapless.len()..columns.min(matching.len())).any(|i| {
        matching[i] > 1 && letters.get(i).copied() != Some(crate::utils::matrix::ncbistdaa::X)
    });
    // `Convert2Matrix` reads the first numRows × numColumns values (a longer list's others are
    // never read; `_ASSERT` only): fewer is a read past the list.
    if weighted.len() < cells || info.len() < columns || matching.len() < columns || gapless_short {
        bail!(
            "an ASCII PSSM (-out_ascii_pssm) of a PSSM file whose informationContent comes without weightedResFreqsPerPos, informationContent, numMatchingSeqs or gaplessColumnWeights for every column (NCBI BLAST+ reads past them) is not supported by LOSAT's PSIBLAST"
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:274-282
    // ```c++
    //             for (SIZE_TYPE c = 0; c < DIM(RESIDUE_ORDER); c++) {
    //                 if ((*pssm)(RESIDUE_ORDER[c], i) != BLAST_SCORE_MIN) {
    //                     double value = 100;
    //                     value *= (*weighted_res_freqs)(RESIDUE_ORDER[c], i);
    //                     // round to the nearest integer
    //                     value = (int)(value + (value > 0. ? 0.5 : -0.5));
    //                     out << setw(4) << (int)value;
    // ```
    // A weighted frequency whose rounded percentage is not an `int` (or is NaN) is undefined
    // behaviour in that conversion (x86-64 prints -2147483648): LOSAT's refusal, for the cells
    // NCBI prints (the engine's frequencies lie in [0, 1]; only a file reaches this).
    let by_row = pssm.by_row == Some(true);
    let scores = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        pssm.final_data
            .as_ref()
            .map_or(&[][..], |data| data.scores.as_slice()),
        rows_n,
        columns,
        by_row,
        BLAST_SCORE_MIN,
    );
    let weighted = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        weighted, rows_n, columns, by_row, 0.0,
    );
    let outside_int = scores.iter().zip(&weighted).any(|(scores, weighted)| {
        super::pssm::RESIDUE_ORDER.iter().any(|&r| {
            let value = 100.0 * weighted[r];
            let rounded = value + if value > 0.0 { 0.5 } else { -0.5 };
            scores[r] != BLAST_SCORE_MIN
                && !(rounded > -2_147_483_649.0 && rounded < 2_147_483_648.0)
        })
    });
    if outside_int {
        bail!(
            "an ASCII PSSM (-out_ascii_pssm) of a PSSM file with a weighted residue frequency whose percentage is not an int or is NaN (NCBI BLAST+'s conversion is undefined behaviour) is not supported by LOSAT's PSIBLAST"
        );
    }
    Ok(())
}

/// The query's letters as NCBIstdaa codes (`CPssm::GetQuerySequenceData`, which converts an
/// IUPACaa or NCBIeaa `Seq-data`; `in_pssm` admits only letters of NCBIstdaa).
fn query_codes(object: &PssmWithParameters) -> Vec<u8> {
    const LETTERS: &[u8; BLASTAA_SIZE] =
        crate::algorithm::pssm_search::in_pssm::NCBISTDAA_TO_AMINOACID;
    match object
        .pssm
        .query
        .as_ref()
        .and_then(|query| query.inst.seq_data.as_ref())
    {
        Some(SeqData::Ncbistdaa(codes)) => codes.clone(),
        Some(SeqData::Iupacaa(text)) | Some(SeqData::Ncbieaa(text)) => text
            .iter()
            .map(|letter| LETTERS.iter().position(|c| c == letter).unwrap_or(0) as u8)
            .collect(),
        None => Vec::new(),
    }
}

impl InPssmRounds {
    /// The files of the input PSSM at the start of round 1 (`SavePssmToFile`): its ASN.1 text and
    /// its ASCII PSSM, as requested.
    pub(crate) fn files(
        &self,
        request: &PssmRequest,
    ) -> Result<(Option<Vec<u8>>, Option<Vec<u8>>)> {
        let Some(object) = self.object.as_ref() else {
            return Ok((None, None));
        };
        let asn = match request.asn_path {
            Some(_) => Some(
                crate::report::serial::pssm::asn_text(object)
                    .map_err(|error| anyhow::anyhow!("{error}"))?,
            ),
            None => None,
        };
        let ascii = request.ascii_path.as_ref().map(|_| ascii_pssm(object));
        Ok((asn, ascii))
    }

    /// The ASN.1 text of the input PSSM after round 1's K/A write-back (`ungapped`: `kbp_psi[0]`,
    /// `gapped`: `kbp_gap_psi[0]` of the round's search), which `-save_pssm_after_last_round`
    /// writes after a round 1 without alignments.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:150-168
    /// ```c++
    ///     // Save the K&A values be as they might have been modified in the
    ///     // composition adjustment library
    ///     if (m_Pssm.NotEmpty()) {
    ///         CPssm& pssm = m_Pssm->SetPssm();
    ///         pssm.SetLambda
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_gap_psi[0]->Lambda);
    ///         pssm.SetKappa
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_gap_psi[0]->K);
    ///         pssm.SetH
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_gap_psi[0]->H);
    ///         pssm.SetLambdaUngapped
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_psi[0]->Lambda);
    ///         pssm.SetKappaUngapped
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_psi[0]->K);
    ///         pssm.SetHUngapped
    ///             (core_data->m_ScoreBlk->GetPointer()->kbp_psi[0]->H);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:157-191 (`CPssm::SetLambda`
    /// and the others set the members of `finalData`)
    /// The redo of the composition-based statistics works on copies of the score block
    /// (`s_BlastScoreBlk_Copy` in `Blast_RedoAlignmentCore_MT`, blast_kappa.c:3254-3259), so the
    /// blocks are those that `PsiBlastSetupScoreBlock` set for the search.
    pub(crate) fn written_back(
        &self,
        request: &PssmRequest,
        ungapped: KarlinBlk,
        gapped: KarlinBlk,
    ) -> Result<Option<Vec<u8>>> {
        let (Some(object), Some(_)) = (self.object.as_ref(), request.asn_path.as_ref()) else {
            return Ok(None);
        };
        // The six values are written in place of the object's (no copy of the object).
        let karlin = crate::report::serial::pssm::FinalKarlin {
            lambda: gapped.lambda,
            kappa: gapped.k,
            h: gapped.h,
            lambda_ungapped: ungapped.lambda,
            kappa_ungapped: ungapped.k,
            h_ungapped: ungapped.h,
        };
        crate::report::serial::pssm::asn_text_with_karlin(object, Some(karlin))
            .map(Some)
            .map_err(|error| anyhow::anyhow!("{error}"))
    }
}

/// `PrintAsciiPssm(*pssm, ExtractPssmAncillaryData(*pssm))` of a file PSSM: the scores of the file
/// (`GetScores`), the query's letters (`GetQuerySequenceData`), the long form with the file's
/// diagnostics when it has `informationContent`, and the PSI blocks of its `finalData` (a member
/// it lacks is `CPssm::kInvalidStat`, 0).
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:792-804
// ```c++
// ExtractPssmAncillaryData(const CPssmWithParameters& pssm)
// {
//     _ASSERT(pssm.CanGetPssm());
//     pair<double, double> lambda, k, h;
//     lambda.first = pssm.GetPssm().GetLambdaUngapped();
//     lambda.second = pssm.GetPssm().GetLambda();
//     k.first = pssm.GetPssm().GetKappaUngapped();
//     k.second = pssm.GetPssm().GetKappa();
//     h.first = pssm.GetPssm().GetHUngapped();
//     h.second = pssm.GetPssm().GetH();
//     return CRef<CBlastAncillaryData>(new CBlastAncillaryData(lambda, k, h, 0,
//                                                              true));
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:117-134
// ```c++
// CBlastAncillaryData::CBlastAncillaryData(pair<double, double> lambda,
//                                          pair<double, double> k,
//                                          pair<double, double> h,
//                                          Int8 effective_search_space,
//                                          bool is_psiblast /* = false */)
// : m_GumbelBlk(0), m_UngappedKarlinBlk(0), m_GappedKarlinBlk(0), m_PsiUngappedKarlinBlk(0),
//   m_PsiGappedKarlinBlk(0), m_SearchSpace(0), m_LengthAdjustment(0)
// {
//     if (is_psiblast) {
//         m_PsiGappedKarlinBlk = Blast_KarlinBlkNew();
//         m_PsiGappedKarlinBlk->Lambda = lambda.second;
//         m_PsiGappedKarlinBlk->K = k.second;
//         m_PsiGappedKarlinBlk->H = h.second;
//
//         m_PsiUngappedKarlinBlk = Blast_KarlinBlkNew();
//         m_PsiUngappedKarlinBlk->Lambda = lambda.first;
//         m_PsiUngappedKarlinBlk->K = k.first;
//         m_PsiUngappedKarlinBlk->H = h.first;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:186-196
// ```c++
//     // We use whether the information content is available to assume whether
//     // the PSSM computation was done or not
//     bool pssm_calculation_done = info_content.empty() ? false : true;
//
//     if (pssm_calculation_done) {
//         out << "\nLast position-specific scoring matrix computed, weighted ";
//         out << "observed percentages rounded down, information per position, ";
//         out << "and relative weight of gapless real matches to pseudocounts\n";
//     } else {
//         out << "\nLast position-specific scoring matrix computed\n";
//     }
// ```
fn ascii_pssm(object: &PssmWithParameters) -> Vec<u8> {
    let pssm = &object.pssm;
    let rows_n = usize::try_from(pssm.num_rows).unwrap_or(0);
    let columns = usize::try_from(pssm.num_columns).unwrap_or(0);
    let by_row = pssm.by_row == Some(true);
    let final_data = pssm.final_data.as_ref();
    let scores = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        final_data.map_or(&[][..], |data| data.scores.as_slice()),
        rows_n,
        columns,
        by_row,
        BLAST_SCORE_MIN,
    );
    let letters = query_codes(object);
    let data = pssm.intermediate_data.as_ref();
    let info = data.and_then(|data| data.information_content.as_deref());
    let weighted;
    let long = match (data, info) {
        (Some(data), Some(info)) if !info.is_empty() => {
            weighted = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
                data.weighted_res_freqs_per_pos.as_deref().unwrap_or(&[]),
                rows_n,
                columns,
                by_row,
                0.0,
            );
            Some(AsciiPssmDiagnostics {
                weighted_residue_freqs: &weighted,
                information_content: info,
                num_matching_seqs: data.num_matching_seqs.as_deref().unwrap_or(&[]),
                gapless_column_weights: data.gapless_column_weights.as_deref().unwrap_or(&[]),
            })
        }
        _ => None,
    };
    let stat = |value: Option<f64>| value.unwrap_or(0.0);
    let psi_ungapped = KarlinBlk {
        lambda: stat(final_data.and_then(|data| data.lambda_ungapped)),
        k: stat(final_data.and_then(|data| data.kappa_ungapped)),
        h: stat(final_data.and_then(|data| data.h_ungapped)),
    };
    let psi_gapped = KarlinBlk {
        lambda: stat(final_data.map(|data| data.lambda)),
        k: stat(final_data.map(|data| data.kappa)),
        h: stat(final_data.map(|data| data.h)),
    };
    super::pssm::write_ascii_pssm(
        &scores,
        &letters,
        long,
        &[
            ("PSI Ungapped         ", psi_ungapped),
            ("PSI Gapped           ", psi_gapped),
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::serial::pssm::{
        FormatRpsDbParameters, ObjectId, Pssm, PssmIntermediateData, PssmParameters, PssmQuery,
        SeqInst,
    };

    /// A PSSM of `columns` positions of the letter M whose query has the Seq-id `id`.
    fn pssm(id: ObjectId, columns: usize, rows: usize) -> PssmWithParameters {
        PssmWithParameters {
            pssm: Pssm {
                is_protein: Some(true),
                identifier: None,
                num_rows: rows as i32,
                num_columns: columns as i32,
                row_labels: None,
                by_row: Some(false),
                query: Some(PssmQuery {
                    ids: vec![ParsedSeqId::Local(id)],
                    titles: None,
                    inst: SeqInst {
                        repr: ("raw", 2),
                        mol: ("aa", 3),
                        length: Some(columns as i32),
                        seq_data: Some(SeqData::Ncbistdaa(vec![12; columns])),
                    },
                }),
                intermediate_data: Some(PssmIntermediateData {
                    freq_ratios: (0..columns * rows).map(|i| i as f64).collect(),
                    ..PssmIntermediateData::default()
                }),
                final_data: None,
            },
            params: Some(PssmParameters {
                pseudocount: Some(0),
                rpsdbparams: Some(FormatRpsDbParameters {
                    matrix_name: b"BLOSUM62".to_vec(),
                    gap_open: None,
                    gap_extend: None,
                }),
                bit_score_thresh: None,
                bit_score_reporting_thresh: None,
            }),
        }
    }

    fn request(asn: bool, ascii: bool) -> PssmRequest {
        PssmRequest {
            asn_path: asn.then(|| "f.asn".to_string()),
            ascii_path: ascii.then(|| "f.txt".to_string()),
            save_each: false,
            pseudocount: 0,
        }
    }

    const GAPPED: KarlinBlk = KarlinBlk {
        lambda: 0.267,
        k: 0.0425,
        h: 0.14,
    };

    // psiblast_aux_priv.cpp:440-502 (quoted at `with_computed_scores`): the ratios of a PSSM of
    // fewer rows get zeros up to 28 rows in each column, the scores are the engine's by column, the
    // gapped values the engine's, and the parameters get the gap costs.
    #[test]
    fn a_ratio_only_pssm_gets_the_constructors_scores() {
        let mut value = pssm(ObjectId::Id(7), 2, 26);
        let rows = vec![[1; BLASTAA_SIZE], [2; BLASTAA_SIZE]];
        with_computed_scores(&mut value, &rows, GAPPED);
        let matrix = &value.pssm;
        assert_eq!(matrix.num_rows, 28);
        let ratios = &matrix.intermediate_data.as_ref().unwrap().freq_ratios;
        assert_eq!(ratios.len(), 56);
        assert_eq!(&ratios[24..30], &[24.0, 25.0, 0.0, 0.0, 26.0, 27.0]);
        assert_eq!(&ratios[54..], &[0.0, 0.0]);
        let data = matrix.final_data.as_ref().unwrap();
        assert_eq!(data.scores.len(), 56);
        assert_eq!((data.scores[27], data.scores[28]), (1, 2));
        assert_eq!((data.lambda, data.kappa, data.h), (0.267, 0.0425, 0.14));
        assert_eq!(data.lambda_ungapped, None);
        let rps = value.params.unwrap().rpsdbparams.unwrap();
        assert_eq!((rps.gap_open, rps.gap_extend), (Some(11), Some(1)));
    }

    // The refusals of `rounds_input`, and what it keeps: the query's Seq-id (an integer stays
    // one), the object only when a file is written.
    #[test]
    fn rounds_input_keeps_the_id_and_refuses_what_ncbi_writes_otherwise() {
        let file = |skipped: Vec<&'static str>, computed: bool, params: bool| {
            let mut object = pssm(ObjectId::Id(7), 2, 28);
            if !params {
                object.params = None;
            }
            InPssmFile {
                object,
                skipped,
                computed,
                valid: true,
            }
        };
        let input = rounds_input(file(vec![], false, true), &request(false, false)).unwrap();
        assert_eq!(
            (input.id.clone(), input.object.is_none()),
            (ParsedSeqId::Local(ObjectId::Id(7)), true)
        );
        let input = rounds_input(file(vec![], false, true), &request(true, false)).unwrap();
        assert!(input.object.is_some());
        let error = rounds_input(
            file(vec!["Seqdesc", "topology"], false, true),
            &request(true, false),
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains("(Seqdesc, topology)"), "{error}");
        // Without -out_pssm the skipped members are never written.
        rounds_input(file(vec!["annot"], false, true), &request(false, true)).unwrap();
        let error = rounds_input(file(vec![], true, false), &request(true, false))
            .unwrap_err()
            .to_string();
        assert!(error.contains("without params.rpsdbparams"), "{error}");
        rounds_input(file(vec![], false, false), &request(true, false)).unwrap();
        // A PSSM that NCBI's checks refuse writes nothing.
        let mut invalid = file(vec!["annot"], true, false);
        invalid.valid = false;
        assert!(rounds_input(invalid, &request(true, true))
            .unwrap()
            .object
            .is_none());
        // The long ASCII form needs every column it reads.
        let mut long = file(vec![], false, true);
        let data = long.object.pssm.intermediate_data.as_mut().unwrap();
        data.information_content = Some(vec![0.5, 0.5]);
        data.num_matching_seqs = Some(vec![2, 2]);
        data.gapless_column_weights = Some(vec![1.0]);
        data.weighted_res_freqs_per_pos = Some(vec![0.0; 56]);
        let error = rounds_input(long, &request(false, true))
            .unwrap_err()
            .to_string();
        assert!(error.contains("informationContent"), "{error}");
    }

    // The review of PS5b, N1 and S1: a longer `weightedResFreqsPerPos` is read as NCBI's
    // `Convert2Matrix` reads it (its first cells); a printed weighted frequency whose rounded
    // percentage is not an `int`, or is NaN, is refused (NCBI's conversion is undefined behaviour).
    #[test]
    fn the_long_ascii_form_reads_the_first_cells_and_refuses_weights_outside_int() {
        let long = |weighted: Vec<f64>| {
            let mut object = pssm(ObjectId::Id(1), 1, 28);
            let mut row = [-1; BLASTAA_SIZE];
            row[0] = BLAST_SCORE_MIN;
            with_computed_scores(&mut object, &[row], GAPPED);
            let data = object.pssm.intermediate_data.as_mut().unwrap();
            data.information_content = Some(vec![0.5]);
            data.num_matching_seqs = Some(vec![1]);
            data.gapless_column_weights = Some(vec![0.0]);
            data.weighted_res_freqs_per_pos = Some(weighted);
            check_long_ascii(&object)
        };
        let mut weighted = vec![0.25; 28];
        long(weighted.clone()).unwrap();
        weighted.extend([7.0, 8.0]);
        long(weighted.clone()).unwrap();
        long(vec![0.25; 27]).unwrap_err();
        // Row 0 (`-`) is not printed: its value is not read.
        let mut unprinted = vec![0.25; 28];
        unprinted[0] = f64::NAN;
        long(unprinted).unwrap();
        for value in [1e8, -1e8, f64::NAN, f64::INFINITY] {
            let mut weighted = vec![0.25; 28];
            weighted[1] = value;
            let error = long(weighted).unwrap_err().to_string();
            assert!(error.contains("not an int or is NaN"), "{value}: {error}");
        }
        let mut edge = vec![0.25; 28];
        edge[1] = 21_474_836.47;
        long(edge).unwrap();
    }

    // psiblast_impl.cpp:150-168 (quoted at `written_back`): the six values of `finalData`.
    #[test]
    fn the_write_back_sets_the_searched_blocks() {
        let mut object = pssm(ObjectId::Str(b"myq".to_vec()), 1, 28);
        with_computed_scores(&mut object, &[[0; BLASTAA_SIZE]], GAPPED);
        let input = InPssmRounds {
            id: ParsedSeqId::Local(ObjectId::Str(b"myq".to_vec())),
            object: Some(object),
        };
        let ungapped = KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            h: 0.401,
        };
        let text = input
            .written_back(&request(true, false), ungapped, GAPPED)
            .unwrap()
            .unwrap();
        let read = crate::blastinput::pssm_text::read_pssm_with_parameters(&text).unwrap();
        let data = read.pssm.final_data.unwrap();
        assert_eq!(
            (data.lambda_ungapped, data.kappa_ungapped, data.h_ungapped),
            (Some(0.3176), Some(0.134), Some(0.401))
        );
        assert_eq!((data.lambda, data.kappa, data.h), (0.267, 0.0425, 0.14));
        assert_eq!(
            input
                .written_back(&request(false, true), ungapped, GAPPED)
                .unwrap(),
            None
        );
    }

    // blastfmtutil.cpp:186-196 (quoted at `ascii_pssm`): the short form of a PSSM without its
    // information content, with the PSI blocks of its file (0 for those it lacks).
    #[test]
    fn a_file_pssm_without_information_content_prints_the_short_form() {
        let mut object = pssm(ObjectId::Id(1), 1, 28);
        let mut row = [-1; BLASTAA_SIZE];
        row[1] = 4;
        row[0] = BLAST_SCORE_MIN;
        with_computed_scores(&mut object, &[row], GAPPED);
        let text = String::from_utf8(ascii_pssm(&object)).unwrap();
        let expected = "\nLast position-specific scoring matrix computed\n           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n    1 M    4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 \n\n                      K         Lambda\nPSI Ungapped         0.0000     0.0000\nPSI Gapped           0.0425     0.2670\n";
        assert_eq!(text, expected);
    }
}
