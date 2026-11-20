//! psiblast `-in_msa F`, `-msa_master_idx N` and `-ignore_msa_master` (series PS, item PS6): the
//! PSSM that NCBI's option processing makes from an alignment file, and the search of it as
//! `-in_pssm` searches a PSSM file.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1253-1275
//! ```c++
//! CRef<CPssmWithParameters>
//! CPsiBlastArgs::x_CreatePssmFromMsa(CNcbiIstream& input_stream,
//!                                    CBlastOptions& opt, bool save_ascii_pssm,
//!                                    unsigned int msa_master_idx,
//!                                    bool ignore_pssm_tmplt_seq)
//! {
//!     // FIXME get these from CBlastOptions
//!     CPSIBlastOptions psiblast_opts;
//!     PSIBlastOptionsNew(&psiblast_opts);
//!     psiblast_opts->nsg_compatibility_mode = ignore_pssm_tmplt_seq;
//!
//!     CPSIDiagnosticsRequest diags(PSIDiagnosticsRequestNewEx(save_ascii_pssm));
//!     CPsiBlastInputClustalW pssm_input(input_stream, *psiblast_opts,
//!                                       opt.GetMatrixName(), diags, NULL, 0,
//!                                       opt.GetGapOpeningCost(),
//!                                       opt.GetGapExtensionCost(),
//!                                       msa_master_idx);
//!     CPssmEngine pssm_engine(&pssm_input);
//!     return pssm_engine.Run();
//! }
//! ```
//! The reader is `blastinput::aln_reader`; this module is `CPsiBlastInputClustalW`
//! (`x_ExtractQueryFromMsa`, `x_ValidateQueryInMsa`, `Process`) and the use of `CPssmEngine` on
//! it (the score block of its constructor, then `Run`), with the PSSM engine of `core::psi`. The
//! object that comes out is searched by the one PSSM search (`pssm_search::in_pssm`), with the
//! rounds, files, strategy and refusals of `-in_pssm`.

use std::path::{Path, PathBuf};

use anyhow::Result;

use crate::blastinput::aln_reader::{read_alignment, Alignment, AlnBioseq, AlnError, Mol};
use crate::blastinput::app::engine_error;
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::core::psi::{
    create_pssm_with_options, KarlinBlk, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiOptions,
    PsiScoreBlk,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::pssm::{PssmQuery, SeqData, SeqEntry, SeqInst};
use crate::utils::matrix::aa_char_to_ncbistdaa;

use super::pssm::{blast_options_error, error_code_to_string, msa_pssm_with_parameters};

/// `-in_msa` with its two companions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MsaInput {
    /// `-in_msa`: the alignment file (`-` is the standard input).
    pub path: PathBuf,
    /// `-msa_master_idx`: the 1-based number of the row that is the query.
    pub master_index: Option<i32>,
    /// `-ignore_msa_master`: the PSSM engine's structure-group mode, which leaves the master out
    /// of the sequence counts.
    pub ignore_master: bool,
}

/// What the PSSM's score block takes from the options at that point of option processing: the
/// matrix and the gap costs (the generic search handler's, which comes before the PSI-BLAST one).
#[derive(Debug, Clone, Copy)]
pub(crate) struct MsaScoring<'a> {
    pub(crate) matrix_name: &'a str,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// `-out_ascii_pssm` is given (`PSIDiagnosticsRequestNewEx(SaveAsciiPssm())`).
    pub(crate) save_ascii_pssm: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
// ```c++
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
fn query_error(message: &str) -> anyhow::Error {
    crate::cli::NativeError {
        exit: 1,
        message: format!("BLAST query error: {message}\n"),
    }
    .into()
}

/// The bytes of `-in_msa`: a file that cannot be opened is NCBI's "File is not accessible" (the
/// open of `CArg_Ios::x_Open`, quoted at `pssm_search::in_pssm::read_in_pssm`); one that opens and
/// cannot be read (a directory) gives the reader no bytes, which no format guess accepts.
fn read_input(path: &Path) -> Result<Vec<u8>> {
    if path.as_os_str() == "-" {
        let mut bytes = Vec::new();
        return std::io::Read::read_to_end(&mut std::io::stdin().lock(), &mut bytes)
            .map(|_| bytes)
            .map_err(|_| crate::cli::inaccessible("in_msa", path));
    }
    std::fs::File::open(path)
        .map(|mut file| {
            let mut bytes = Vec::new();
            let _ = std::io::Read::read_to_end(&mut file, &mut bytes);
            bytes
        })
        .map_err(|_| crate::cli::inaccessible("in_msa", path))
}

/// LOSAT's refusal of what the alignment reader does not read yet.
fn refused(what: impl std::fmt::Display) -> anyhow::Error {
    anyhow::anyhow!("-in_msa: {what} is not supported by LOSAT's PSIBLAST yet")
}

/// The letters of a Bioseq as NCBIstdaa after `CSeqportUtil::Convert` from IUPACaa
/// (`CIupacaaToStdaa`): the letters A to Z; any other byte (`*` among them) is 0xff.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1690 (the table's
/// upper-case half: A=1, B=2, ... I=9, J=0x1b, K=0xa, ... O=0x1a, ... U=0x18, ... Z=0x17; all
/// other entries 0xff)
fn iupacaa_to_ncbistdaa(letter: u8) -> u8 {
    if letter.is_ascii_uppercase() {
        aa_char_to_ncbistdaa(letter)
    } else {
        0xff
    }
}

/// `CPsiBlastInputClustalW::x_ExtractQueryForPssm`: the first Bioseq, in the file's order, of
/// the query's length whose letters as NCBIstdaa are the query's; `None` where no Bioseq is
/// (a `*` in the query: IUPACaa has none). A Bioseq of that length that is a nucleotide one stops
/// the run (`CSeqportUtil::Convert` cannot convert it).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:178-199
/// ```c++
///     for (CTypeIterator<CBioseq> itr(Begin(*m_SeqEntry)); itr; ++itr) {
///         _ASSERT(itr->IsAa());
///         if (itr->GetLength() != GetQueryLength()) {
///             continue;
///         }
///         // let's check the sequence data
///         CNCBIstdaa sequence;
///         s_GetQuerySequenceData(*itr, GetQueryLength(), sequence);
///         if (s_AreSequencesEqual(sequence, m_Query.get())) {
///             m_QueryBioseq.Reset(&*itr);
///             break;
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:221-226 (the throw of a
/// conversion that no function covers; NCBI BLAST+ 2.17.0, PS6 corpus `e.len1`, `n.q`)
/// ```c++
///     NCBI_THROW(CSeqUtilException, eInvalidCoding,
///                "Unknown conversion.");
/// ```
fn query_bioseq<'a>(
    bioseqs: &'a [AlnBioseq],
    query: &[u8],
) -> Result<Option<&'a AlnBioseq>, NcbiException> {
    for bioseq in bioseqs {
        if bioseq.letters.len() != query.len() {
            continue;
        }
        if bioseq.mol == Mol::Nucleotide {
            return Err(NcbiException::new(
                ThrowSite {
                    file: "src/util/sequtil/sequtil_convert_imp.cpp",
                    line: 225,
                    kind: "CSeqUtilException::Invalid coding",
                    function: Some("ncbi::CSeqConvert_imp::Convert"),
                },
                "Unknown conversion.",
            ));
        }
        if bioseq
            .letters
            .iter()
            .zip(query)
            .all(|(&letter, &code)| iupacaa_to_ncbistdaa(letter) == code)
        {
            return Ok(Some(bioseq));
        }
    }
    Ok(None)
}

/// `x_ExtractQueryFromMsa`: the master row's letters without gaps, as NCBIstdaa.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:264-306
/// ```c++
///     if (msa_master_idx >= m_AsciiMsa.size()) {
///         CNcbiOstrstream oss;
///         oss << "Invalid master sequence index, please use a value between 1 "
///             << "and " << m_AsciiMsa.size();
///         NCBI_THROW(CBlastException, eInvalidOptions,
///                    CNcbiOstrstreamToString(oss));
///     }
///     const string& kQuery = m_AsciiMsa.at(msa_master_idx);
///     ...
///         m_Query.get()[query_idx] = AMINOACID_TO_NCBISTDAA[toupper(*residue)];
/// ```
fn extract_query(rows: &[Vec<u8>], master_index: Option<i32>) -> Result<Vec<u8>> {
    // `AsInteger() - 1` of a value of 1 or more.
    let index = master_index.map_or(0, |index| (index - 1) as usize);
    if index >= rows.len() {
        return Err(blast_options_error(&format!(
            "Invalid master sequence index, please use a value between 1 and {}",
            rows.len()
        )));
    }
    Ok(rows[index]
        .iter()
        .filter(|&&residue| residue != b'-')
        .map(|&residue| aa_char_to_ncbistdaa(residue.to_ascii_uppercase()))
        .collect())
}

/// `x_ValidateQueryInMsa`: the first row whose letters start with the query's (a row that goes on
/// after the query is one too) is swapped with the first row; none is NCBI's options error.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:228-262
/// ```c++
///     const size_t kAligmentLength = m_AsciiMsa.front().size();
///     const char kMaskingRes = NCBISTDAA_TO_AMINOACID[kProtMask];
///     ...
///             /* Selenocysteines are replaced by X's in query; test for this
///              * possibility */
///             if (query_res == kMaskingRes && kCurrentRes == 'U') {
///                 query_res = kCurrentRes;
///             }
///     ...
///         if (query_idx == GetQueryLength()) {
///             break;
///         }
///     ...
///         string msg("No sequence in the multiple sequence alignment provided ");
///         msg += "matches the query sequence";
///         NCBI_THROW(CBlastException, eInvalidOptions, msg);
/// ```
fn validate_query(rows: &mut [Vec<u8>], query: &[u8]) -> Result<()> {
    let masking = crate::algorithm::pssm_search::in_pssm::NCBISTDAA_TO_AMINOACID[21];
    let found = rows.iter().position(|row| {
        let mut query_index = 0;
        for &residue in row.iter() {
            if query_index >= query.len() {
                break;
            }
            if residue == b'-' {
                continue;
            }
            let mut query_residue = crate::algorithm::pssm_search::in_pssm::NCBISTDAA_TO_AMINOACID
                [query[query_index] as usize];
            let current = residue.to_ascii_uppercase();
            if query_residue == masking && current == b'U' {
                query_residue = current;
            }
            if query_residue != current {
                break;
            }
            query_index += 1;
        }
        query_index == query.len()
    });
    match found {
        Some(index) => {
            rows.swap(0, index);
            Ok(())
        }
        None => Err(blast_options_error(
            "No sequence in the multiple sequence alignment provided matches the query sequence",
        )),
    }
}

/// `x_CopyQueryToMsa` and `x_ExtractAlignmentData`: the engine's alignment of the rows (the query
/// first, in the file's case: a lower-case residue of the query is not aligned), with the flanking
/// gaps and the internal gaps of 10 or more of each other row unaligned.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:308-405
/// ```c++
///         m_Msa->data[kQueryIndex][query_idx].letter = m_Query.get()[query_idx];
///         m_Msa->data[kQueryIndex][query_idx].is_aligned =
///             (isupper(*residue) ? true : false);
///     ...
///             m_Msa->data[seq_index][query_idx].letter =
///                 AMINOACID_TO_NCBISTDAA[(int) toupper(kCurrentRes)];
///             m_Msa->data[seq_index][query_idx].is_aligned = true;
///     ...
///     // set flanking gaps and long internal gaps as unaligned
///     const int kGapResidue = 0;
///     const int kLongGapLen = 10;
///     for (; seq_index < m_Msa->dimensions->num_seqs + 1; seq_index++) {
///         size_t i = 0;
///
///         // find left flanking gaps
///         while (i < m_Msa->dimensions->query_length &&
///                m_Msa->data[seq_index][i].letter == kGapResidue) {
///             m_Msa->data[seq_index][i].is_aligned = false;
///             i++;
///         }
///
///         // find long internal gaps
///         while (i < m_Msa->dimensions->query_length) {
///             while (i < m_Msa->dimensions->query_length &&
///                    m_Msa->data[seq_index][i].letter != kGapResidue) {
///                 i++;
///             }
///
///             int k = static_cast<int>(i) + 1;
///             while (k < m_Msa->dimensions->query_length &&
///                    m_Msa->data[seq_index][k].letter == kGapResidue) {
///                 k++;
///             }
///
///             if (k - i >= kLongGapLen) {
///                 for (int j=static_cast<int>(i);j < k;j++) {
///                     m_Msa->data[seq_index][j].is_aligned = false;
///                 }
///             }
///
///             i = k;
///         }
///
///         // find right flanking gaps
///         int k = m_Msa->dimensions->query_length - 1;
///         while (k >= 0 && m_Msa->data[seq_index][k].letter == kGapResidue) {
///             m_Msa->data[seq_index][k].is_aligned = false;
///             k--;
///         }
///     }
/// ```
fn engine_msa(rows: &[Vec<u8>], query: &[u8]) -> PsiMsa {
    const LONG_GAP: usize = 10;
    let length = query.len();
    let mut msa = PsiMsa::new(length, rows.len() - 1);
    let mut query_index = 0;
    for &residue in &rows[0] {
        if residue == b'-' {
            continue;
        }
        let cell = &mut msa.rows[0][query_index];
        cell.letter = query[query_index];
        cell.is_aligned = residue.is_ascii_uppercase();
        query_index += 1;
    }
    for (seq, row) in rows.iter().enumerate().skip(1) {
        let mut query_index = 0;
        for (column, &master) in rows[0].iter().enumerate() {
            if master == b'-' {
                continue;
            }
            let cell = &mut msa.rows[seq][query_index];
            cell.letter = aa_char_to_ncbistdaa(row[column].to_ascii_uppercase());
            cell.is_aligned = true;
            query_index += 1;
        }
        let cells = &mut msa.rows[seq];
        let mut i = 0;
        while i < length && cells[i].letter == 0 {
            cells[i].is_aligned = false;
            i += 1;
        }
        while i < length {
            while i < length && cells[i].letter != 0 {
                i += 1;
            }
            let mut k = i + 1;
            while k < length && cells[k].letter == 0 {
                k += 1;
            }
            if k - i >= LONG_GAP {
                for cell in &mut cells[i..k] {
                    cell.is_aligned = false;
                }
            }
            i = k;
        }
        let mut k = length;
        while k > 0 && cells[k - 1].letter == 0 {
            cells[k - 1].is_aligned = false;
            k -= 1;
        }
    }
    msa
}

/// The score block of `CPssmEngine::x_InitializeScoreBlock` for the options' matrix and gap
/// costs: NCBI's errors of a matrix it does not know and of gap costs that the matrix's table
/// does not hold, with the gapped block `kbp_gap_psi` copies. `None` for IDENTITY, whose table has
/// no frequency ratios (the engine stops there).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:327-342 and 395-410
/// ```c
///     status = Blast_ScoreBlkMatrixInit(program_number, scoring_options, sbp, get_path);
///     if (status) {
///         Blast_PerrorWithLocation(blast_message, status, kBlastMessageNoContext);
///         return status;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3527-3560
/// ```c
///       else if (status == 2)
///       {
///          snprintf(buffer, sizeof(buffer), "Gap existence and extension values of %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, matrix_name);
///          Blast_MessageWrite(error_return, eBlastSevError, kBlastMessageNoContext, buffer);
///          BlastKarlinReportAllowedValues(matrix_name, error_return);
///       }
/// ```
/// (The message that `CPssmEngine` throws is the first one written; the allowed values follow it
/// as further messages. NCBI BLAST+ 2.17.0, PS6 corpus `b.gap.bad`, `o.ok.matrix`.)
fn score_block(scoring: &MsaScoring<'_>) -> Result<(ScoringMatrix, Option<PsiScoreBlk>)> {
    let Ok(matrix) = scoring.matrix_name.parse::<ScoringMatrix>() else {
        return Err(engine_error("Unknown error code -1"));
    };
    let spec = ProteinScoringSpec {
        matrix,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
    };
    if !crate::stats::tables::protein_scoring_supported(&spec) {
        return Err(engine_error(&format!(
            "Gap existence and extension values of {} and {} not supported for {}",
            scoring.gap_open,
            scoring.gap_extend,
            scoring.matrix_name.to_ascii_uppercase()
        )));
    }
    let params = crate::stats::tables::lookup_protein_params(&spec);
    let gapped = KarlinBlk {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
    };
    let sbp = if matrix == ScoringMatrix::Identity {
        None
    } else {
        Some(crate::algorithm::pssm_search::in_pssm::score_block(
            matrix, gapped,
        )?)
    };
    Ok((matrix, sbp))
}

/// Whether the ungapped Karlin-Altschul parameters of the query's composition can be computed:
/// NCBI's score block leaves the query's context invalid where they cannot (a query of `X`
/// alone, or of letters that no positive score pairs with the standard composition), and the
/// engine then reads a null block.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2776-2800
/// ```c
///       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
///       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
///       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
///       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
///       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
///           sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
///           sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:353-370 (the ungapped
/// calculation's status is replaced by the gapped one's, so the failure goes unreported and the
/// run goes on with `kbp_psi[0]` null; `s_PSISavePssm` and the scaling read it: NCBI BLAST+ 2.17.0
/// ends with SIGSEGV, PS6 corpus `r.ch.q.allX`, `x.in.allx.60`, `x.in.allO`)
/// ```c
///        status = Blast_ScoreBlkKbpUngappedCalc(program_number, sbp, query_blk->sequence,
///                query_info, blast_message);
///
///        if (scoring_options->gapped_calculation) {
///           status =
///               Blast_ScoreBlkKbpGappedCalc(sbp, scoring_options, program_number,
///                                           query_info, blast_message);
/// ```
fn ungapped_karlin_valid(matrix: ScoringMatrix, query: &[u8]) -> bool {
    let (_, loscore, hiscore) = crate::algorithm::pssm_search::in_pssm::score_table(matrix);
    let mut sentinelled = Vec::with_capacity(query.len() + 2);
    sentinelled.push(0);
    sentinelled.extend_from_slice(query);
    sentinelled.push(0);
    let composition = crate::stats::karlin_calc::compute_aa_composition(&sentinelled, query.len());
    let standard = crate::stats::karlin_calc::compute_std_aa_composition();
    let profile = crate::stats::karlin_calc::compute_score_freq_profile_for_matrix(
        &composition,
        &standard,
        loscore,
        hiscore,
        matrix,
    );
    crate::stats::karlin_calc::compute_karlin_params_ungapped(&profile).is_ok()
}

/// The PSSM that NCBI makes of the alignment `input` (`x_CreatePssmFromMsa`), or its error at the
/// point of the run where NCBI raises it, or LOSAT's refusal of what NCBI crashes on or what the
/// alignment reader does not read.
///
/// The PSSM has no query Seq-entry where no Bioseq of the alignment is the query (a `*` in the
/// master): `None` in the second place.
pub(crate) fn pssm_from_msa(
    input: &MsaInput,
    scoring: &MsaScoring<'_>,
) -> Result<crate::report::serial::pssm::PssmWithParameters> {
    let bytes = read_input(&input.path)?;
    let alignment: Alignment = match read_alignment(&bytes) {
        Ok(alignment) => alignment,
        Err(AlnError::NotAvailable) => {
            return Err(query_error(
                crate::blastinput::aln_reader::SEQ_ENTRY_NOT_AVAILABLE,
            ))
        }
        Err(AlnError::InvalidMolType) => {
            return Err(query_error(crate::blastinput::aln_reader::INVALID_MOL_TYPE))
        }
        Err(AlnError::Unsupported(format)) => {
            return Err(refused(format!("an alignment file in {format} format")))
        }
        Err(AlnError::Refused(what)) => return Err(refused(what)),
    };
    let mut rows = alignment.rows.clone();
    let query = extract_query(&rows, input.master_index)?;
    validate_query(&mut rows, &query)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:64-84, 197-200
    // ```c++
    //     if (pssm_input_msa->GetQueryLength() == 0) {
    //         NCBI_THROW(CPssmEngineException, eNullInputData,
    //            "Query length provided by IPssmInputData is 0");
    //     }
    // ...
    // CPssmEngine::CPssmEngine(IPssmInputData* input)
    //     : m_PssmInput(input), m_PssmInputFreqRatios(NULL)
    // {
    //     s_CheckAgainstNullData(input);
    //     x_InitializeScoreBlock(x_GetQuery(), x_GetQueryLength(), x_GetMatrixName(),
    //            x_GetGapExistence(), x_GetGapExtension());
    // }
    // ```
    if query.is_empty() {
        return Err(engine_error("Query length provided by IPssmInputData is 0"));
    }
    let (matrix, sbp) = score_block(scoring)?;
    // `Process` before the engine: the query's Bioseq is looked for last, and may stop the run.
    let msa = engine_msa(&rows, &query);
    let bioseq = query_bioseq(&alignment.bioseqs, &query)
        .map_err(|exception| anyhow::Error::from(crate::cli::NativeError::exception(&exception)))?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:366-394
    // ```c++
    //     m_PssmInput->Process();
    //     s_Validate(m_PssmInput);
    //     ...
    //     int status =
    //         PSICreatePssmWithDiagnostics(m_PssmInput->GetData(),
    //                                      m_PssmInput->GetOptions(),
    //                                      m_ScoreBlk,
    //                                      m_PssmInput->GetDiagnosticsRequest(),
    //                                      &pssm,
    //                                      &diagnostics);
    //     if (status != PSI_SUCCESS) {
    //         string msg = x_ErrorCodeToString(status);
    //         NCBI_THROW(CBlastException, eCoreBlastError, msg);
    //     }
    // ```
    // The options are `PSIBlastOptionsNew`'s (pseudocount 0 whatever `-pseudocount` says, which
    // acts on the later rounds) with the unaligned positions ignored and the structure-group mode
    // of `-ignore_msa_master`.
    let pssm = create_pssm_with_options(
        &msa,
        0,
        PsiOptions {
            ignore_unaligned_positions: true,
            nsg_compatibility_mode: input.ignore_master,
            ungapped_block_missing: !ungapped_karlin_valid(matrix, &query),
        },
        sbp.as_ref(),
        PsiDiagnosticsRequest::new_ex(scoring.save_ascii_pssm),
    )
    .map_err(|error| match error {
        PsiError::NoMatrixFrequencyData => anyhow::anyhow!(
            "-in_msa with -matrix {} is not supported by LOSAT's PSIBLAST (NCBI BLAST+ 2.17.0 stops with a crash: the matrix has no frequency ratios)",
            scoring.matrix_name.to_ascii_uppercase()
        ),
        PsiError::NoUngappedKarlinBlock => anyhow::anyhow!(
            "-in_msa with a master sequence whose ungapped Karlin-Altschul parameters cannot be computed (a sequence of X alone, for one) is not supported by LOSAT's PSIBLAST (NCBI BLAST+ 2.17.0 stops with a crash)"
        ),
        other => engine_error(error_code_to_string(other)),
    })?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:395-403
    // ```c++
    //     retval = x_PSIMatrix2Asn1(pssm, m_PssmInput->GetMatrixName(),
    //                               m_PssmInput->GetOptions(), diagnostics);
    //     CRef<CBioseq> query = m_PssmInput->GetQueryForPssm();
    //     if (query.NotEmpty()) {
    //         retval->SetQuery().SetSeq(*query);
    //     }
    // ```
    // The Bioseq is the alignment's own: its ids, its letters as IUPACaa, no title.
    let entry = bioseq.map(|bioseq| {
        SeqEntry::Seq(Box::new(PssmQuery {
            ids: bioseq.ids.clone(),
            descr: None,
            inst: SeqInst::new(
                ("raw", 2),
                ("aa", 3),
                Some(bioseq.letters.len() as u32),
                Some(SeqData::Iupacaa(bioseq.letters.clone())),
            ),
            annot: None,
        }))
    });
    Ok(msa_pssm_with_parameters(
        &pssm,
        query.len(),
        entry,
        scoring.matrix_name,
    ))
}

/// The PSSM search of `-in_msa`: the PSSM of [`pssm_from_msa`] through the one PSSM search
/// (`pssm_search::in_pssm::psi_search_and_object`), as a PSSM file's.
pub(crate) fn read_in_msa(
    input: &MsaInput,
    scoring: &MsaScoring<'_>,
    reader: crate::algorithm::pssm_search::in_pssm::PssmReader<'_>,
) -> Result<(
    crate::algorithm::pssm_search::PsiSearch,
    crate::algorithm::pssm_search::in_pssm::InPssmFile,
)> {
    let pssm = pssm_from_msa(input, scoring)?;
    let (mut search, file) = if pssm.pssm.query.is_none() {
        crate::algorithm::pssm_search::in_pssm::psi_search_unset_query(pssm)
    } else {
        crate::algorithm::pssm_search::in_pssm::psi_search_and_object(pssm, reader)?
    };
    search.from_msa = true;
    Ok((search, file))
}

#[cfg(test)]
#[path = "in_msa_tests.rs"]
mod tests;
