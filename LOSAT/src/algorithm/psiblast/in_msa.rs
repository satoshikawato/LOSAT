//! psiblast `-in_msa F`, `-msa_master_idx N` and `-ignore_msa_master` (series PS, item PS6): the
//! PSSM that NCBI's option processing makes from an alignment file, and the search of it as
//! `-in_pssm` searches a PSSM file.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1252-1271
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

use std::path::PathBuf;

use anyhow::Result;

use crate::blastinput::aln_reader::{read_alignment, Alignment, AlnBioseq, AlnError, Mol};
use crate::blastinput::app::{engine_error, query_error};
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

// What the PSSM engine reads from `CPsiBlastInputClustalW`: an empty matrix name and a gap cost
// of 0 are "not given" and stand for the defaults of `IPssmInputData` (BLOSUM62, 11, 1), whatever
// the options hold. So `-gapopen 0` builds the PSSM with 11, and the search's own option check
// then reports the options' 0 (NCBI BLAST+ 2.17.0, review probes `z.go0*`, `z.ge0*` of line PS
// item PS6; `-matrix ''`: `PS/logs/fix6/probe_mx`).
//
// NCBI reference (598d8ae6): c++/include/algo/blast/api/msa_pssm_input.hpp:104-122
// ```c++
//     const char* GetMatrixName() {
//         return m_MatrixName.empty()
//             ? IPssmInputData::GetMatrixName()
//             : m_MatrixName.c_str();
//     }
//
//    /// Obtain the gap existence value to use when building the PSSM
//     int GetGapExistence() {
//          return m_GapExistence
//          ? m_GapExistence
//          : IPssmInputData::GetGapExistence();
//     }
//
//     /// Obtain the gap extension value to use when building the PSSM
//     int GetGapExtension() {
//          return m_GapExtension
//          ? m_GapExtension
//          : IPssmInputData::GetGapExtension();
//     }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/pssm_input.hpp:68-80
// ```c++
//     virtual const char* GetMatrixName() {
//         return BLAST_DEFAULT_MATRIX;
//     }
//
//     /// Obtain the gap existence value for the underlying matrix used to build the PSSM.
//     virtual int GetGapExistence() {
//         return BLAST_GAP_OPEN_PROT;
//     }
//
//     /// Obtain the gap extension value for the underlying matrix used to build the PSSM.
//     virtual int GetGapExtension() {
//         return BLAST_GAP_EXTN_PROT;
//     }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:77
// ```c
// #define BLAST_DEFAULT_MATRIX "BLOSUM62"
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:84-85
// ```c
// #define BLAST_GAP_OPEN_PROT 11 /**< default gap open penalty (all
//                                     protein searches) */
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:92-93
// ```c
// #define BLAST_GAP_EXTN_PROT 1 /**< default gap open penalty (all
//                                    protein searches) */
// ```
const BLAST_DEFAULT_MATRIX: &str = "BLOSUM62";
const BLAST_GAP_OPEN_PROT: i32 = 11;
const BLAST_GAP_EXTN_PROT: i32 = 1;

impl MsaScoring<'_> {
    /// `GetMatrixName()`: the options' matrix, or `BLAST_DEFAULT_MATRIX` for an empty name.
    fn engine_matrix_name(&self) -> &str {
        if self.matrix_name.is_empty() {
            BLAST_DEFAULT_MATRIX
        } else {
            self.matrix_name
        }
    }

    /// `GetGapExistence()`: the options' gap opening cost, or `BLAST_GAP_OPEN_PROT` for 0.
    fn engine_gap_existence(&self) -> i32 {
        if self.gap_open != 0 {
            self.gap_open
        } else {
            BLAST_GAP_OPEN_PROT
        }
    }

    /// `GetGapExtension()`: the options' gap extension cost, or `BLAST_GAP_EXTN_PROT` for 0.
    fn engine_gap_extension(&self) -> i32 {
        if self.gap_extend != 0 {
            self.gap_extend
        } else {
            BLAST_GAP_EXTN_PROT
        }
    }
}

/// LOSAT's refusal of what the alignment reader does not read yet.
fn refused(what: impl std::fmt::Display) -> anyhow::Error {
    anyhow::anyhow!("-in_msa: {what} is not supported by LOSAT's PSIBLAST yet")
}

/// The letters of a Bioseq as NCBIstdaa after `CSeqportUtil::Convert` from IUPACaa
/// (`CIupacaaToStdaa`): the letters A to Z; any other byte (`*` among them) is 0xff.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1686
/// ```c++
/// const Uint1 CIupacaaToStdaa::scm_Table[256] = {
/// ...
///     0xff, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7,
///     0x8, 0x9, 0x1b, 0xa, 0xb, 0xc, 0xd, 0x1a,
///     0xe, 0xf, 0x10, 0x11, 0x12, 0x18, 0x13, 0x14,
///     0x15, 0x16, 0x17, 0xff, 0xff, 0xff, 0xff, 0xff,
/// ```
/// (The table's upper-case half: A=1, B=2, ... I=9, J=0x1b, K=0xa, ... O=0x1a, ... U=0x18, ...
/// Z=0x17; all other entries 0xff.)
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:189-201
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
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:225
/// ```c++
///     NCBI_THROW(CSeqUtilException, eInvalidCoding, "Unknown conversion.");
/// ```
/// (The throw of a conversion that no function covers; NCBI BLAST+ 2.17.0, PS6 corpus `e.len1`, `n.q`.)
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:276-300
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
/// Returns the place that row had in the file (0-based).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:228-270
/// ```c++
///     const size_t kAligmentLength = m_AsciiMsa.front().size();
///     const char kMaskingRes = NCBISTDAA_TO_AMINOACID[kProtMask];
///     _ASSERT( !m_AsciiMsa.empty() );
///
///     size_t seq_idx = 0;
///     for (; seq_idx < m_AsciiMsa.size(); seq_idx++) {
///         size_t query_idx = 0;
///         for (size_t align_idx = 0;
///              align_idx < kAligmentLength && query_idx < GetQueryLength();
///              align_idx++) {
///             if (m_AsciiMsa[seq_idx][align_idx] == kGapChar) {
///                 continue;
///             }
///             char query_res = NCBISTDAA_TO_AMINOACID[m_Query.get()[query_idx]];
///             const char kCurrentRes = toupper(m_AsciiMsa[seq_idx][align_idx]);
///             /* Selenocysteines are replaced by X's in query; test for this
///              * possibility */
///             if (query_res == kMaskingRes && kCurrentRes == 'U') {
///                 query_res = kCurrentRes;
///             }
///             if (query_res != kCurrentRes) {
///                 break;  // character mismatch
///             } else {
///                 query_idx++;
///             }
///         }
///
///         if (query_idx == GetQueryLength()) {
///             break;
///         }
///     }
///
///     if (seq_idx < m_AsciiMsa.size()) {
///         // If the query was found at position seq_idx, swap it with the first
///         // element in the m_AsciiMsa vector
///         for (size_t align_idx = 0; align_idx < kAligmentLength; align_idx++) {
///             swap(m_AsciiMsa[seq_idx][align_idx], m_AsciiMsa.front()[align_idx]);
///         }
///     } else {
///         string msg("No sequence in the multiple sequence alignment provided ");
///         msg += "matches the query sequence";
///         NCBI_THROW(CBlastException, eInvalidOptions, msg);
///     }
/// ```
fn validate_query(rows: &mut [Vec<u8>], query: &[u8]) -> Result<usize> {
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
            Ok(index)
        }
        None => Err(blast_options_error(
            "No sequence in the multiple sequence alignment provided matches the query sequence",
        )),
    }
}

/// The rows after `x_ValidateQueryInMsa`'s swap whose first row has exactly the query's letters:
/// the only input of [`engine_msa`], and made only by [`QueryRows::new`], so the engine's
/// alignment is never filled from a first row with more letters than the query.
///
/// The first row is the first one whose letters start with the query's, so it can have more:
/// with `-msa_master_idx`, a master that is the start of a sequence above it in the file. NCBI
/// then fills its alignment (`query_length` cells a row) by that longer row's letters and reads
/// past the query and writes past every row (the two functions quoted at [`engine_msa`]; NCBI
/// BLAST+ 2.17.0, review probes `p.pfx.*`, `u.pfx.*` of line PS item PS6: exit 0 with a PSSM of
/// the longer row's columns when it is longer by one or two letters, `malloc(): invalid size`
/// and SIGABRT when by many). That result depends on the allocator, so LOSAT refuses the run
/// (`pssm_from_msa`; DECISIONS 75.1 of line PS).
struct QueryRows<'a> {
    /// The first row after the swap (the query's letters with the row's gaps, in the file's case).
    front: &'a [u8],
    /// The rows after it.
    others: &'a [Vec<u8>],
    /// The query's letters as NCBIstdaa: as many as `front` has letters.
    query: &'a [u8],
}

impl<'a> QueryRows<'a> {
    /// `None` where there is no row or the first row has not exactly `query.len()` letters.
    fn new(rows: &'a [Vec<u8>], query: &'a [u8]) -> Option<Self> {
        let (front, others) = rows.split_first()?;
        let letters = front.iter().filter(|&&residue| residue != b'-').count();
        (letters == query.len()).then_some(Self {
            front,
            others,
            query,
        })
    }
}

/// `x_CopyQueryToMsa` and `x_ExtractAlignmentData`: the engine's alignment of the rows (the query
/// first, in the file's case: a lower-case residue of the query is not aligned), with the flanking
/// gaps and the internal gaps of 10 or more of each other row unaligned.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:314-325
/// ```c++
///     const string& ascii_query = m_AsciiMsa.front();
///
///     unsigned int query_idx = 0;
///     ITERATE(string, residue, ascii_query) {
///         if (*residue == kGapChar) {
///             continue;
///         }
///         m_Msa->data[kQueryIndex][query_idx].letter = m_Query.get()[query_idx];
///         m_Msa->data[kQueryIndex][query_idx].is_aligned =
///             (isupper(*residue) ? true : false);
///         query_idx++;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:335-396
/// ```c++
///     size_t seq_index = kQueryIndex + 1;
///     for (; seq_index < m_AsciiMsa.size(); seq_index++) {
///         size_t query_idx = 0;
///         for (size_t align_idx = 0; align_idx < kAlignmentLength; align_idx++) {
///             if (m_AsciiMsa.front()[align_idx] == kGapChar) {
///                 continue;
///             }
///             _ASSERT(toupper(m_AsciiMsa.front()[align_idx]) ==
///                     NCBISTDAA_TO_AMINOACID[m_Query.get()[query_idx]]);
///             const char kCurrentRes = m_AsciiMsa[seq_index][align_idx];
///             _ASSERT(isalpha(kCurrentRes) || kCurrentRes == kGapChar);
///             m_Msa->data[seq_index][query_idx].letter =
///                 AMINOACID_TO_NCBISTDAA[(int) toupper(kCurrentRes)];
///             m_Msa->data[seq_index][query_idx].is_aligned = true;
///             query_idx++;
///         }
///     }
///
///     // set flanking gaps and long internal gaps as unaligned
///     seq_index = kQueryIndex + 1;
///     const int kGapResidue = 0;
///     const int kLongGapLen = 10;
///     for (; seq_index < m_Msa->dimensions->num_seqs + 1; seq_index++) {
///         size_t i = 0;
///
///         // find left flanking gaps
///         while (i < m_Msa->dimensions->query_length &&
///                m_Msa->data[seq_index][i].letter == kGapResidue) {
///
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
/// The cells of a row are taken by the letters of the first row in step (`zip`), never by an
/// index: [`QueryRows`] holds a first row of exactly `query.len()` letters, the number of cells
/// of each row, so every cell is written once and none is written past the row.
fn engine_msa(aligned: &QueryRows<'_>) -> PsiMsa {
    const LONG_GAP: usize = 10;
    let QueryRows {
        front,
        others,
        query,
    } = *aligned;
    let mut msa = PsiMsa::new(query.len(), others.len());
    let mut msa_rows = msa.rows.iter_mut();
    let residues = front.iter().filter(|&&residue| residue != b'-');
    let master_cells = msa_rows.next().into_iter().flatten();
    for ((cell, &residue), &letter) in master_cells.zip(residues).zip(query) {
        cell.letter = letter;
        cell.is_aligned = residue.is_ascii_uppercase();
    }
    for (cells, row) in msa_rows.zip(others) {
        let letters = front
            .iter()
            .zip(row)
            .filter(|(&master, _)| master != b'-')
            .map(|(_, &residue)| residue);
        for (cell, residue) in cells.iter_mut().zip(letters) {
            cell.letter = aa_char_to_ncbistdaa(residue.to_ascii_uppercase());
            cell.is_aligned = true;
        }
        let length = cells.len();
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

/// The score block of `CPssmEngine::x_InitializeScoreBlock` for the matrix and gap costs that the
/// engine reads ([`MsaScoring`]'s accessors) and the query: NCBI's errors of a matrix it does not
/// know and of gap costs that the matrix's table does not hold, with the gapped block
/// `kbp_gap_psi` copies. `None` for IDENTITY, whose table has no frequency ratios (the engine
/// stops there).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:498-520
/// ```c
///     status = Blast_ScoreBlkMatrixInit(program_number, scoring_options, sbp, get_path);
///     if (status) {
///         Blast_PerrorWithLocation(blast_message, status, kBlastMessageNoContext);
///         return status;
///     }
///
///     /* Fills in block for gapped blast. */
///     if (Blast_ProgramIsPhiBlast(program_number)) {
///        status = s_PHIScoreBlkFill(sbp, scoring_options, blast_message, get_path);
///     } else if (Blast_ProgramIsMapping(program_number)) {
///         /* Create fake score blocks for each query without computing base
///            frequencies. We do not compute e-values for mapping, so the
///            KA statistics are not needed, but the data structures are checked
///            and used in BLAST engine. */
///         status = s_JumperScoreBlkFill(sbp, query_info, blast_message);
///     } else {
///        status = Blast_ScoreBlkKbpUngappedCalc(program_number, sbp, query_blk->sequence,
///                query_info, blast_message);
///
///        if (scoring_options->gapped_calculation) {
///           status =
///               Blast_ScoreBlkKbpGappedCalc(sbp, scoring_options, program_number,
///                                           query_info, blast_message);
/// ```
/// The ungapped calculation writes its message for a query whose parameters it cannot compute
/// and its status is then replaced by the gapped calculation's:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2783-2790
/// ```c
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
///           sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
///           sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
///           if (!Blast_QueryIsTranslated(program) ) {
///              Blast_MessageWrite(blast_message, eBlastSevWarning, context,
///              kBlastErrMsg_CantCalculateUngappedKAParams);
///           }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3557-3562
/// ```c
///       else if (status == 2)
///       {
///          snprintf(buffer, sizeof(buffer), "Gap existence and extension values of %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, matrix_name);
///          Blast_MessageWrite(error_return, eBlastSevError, kBlastMessageNoContext, buffer);
///          BlastKarlinReportAllowedValues(matrix_name, error_return);
///       }
/// ```
/// The engine throws the first message that was written, so with gap costs that the matrix does
/// not have the error is the ungapped calculation's for such a query and the gap costs' for any
/// other (whose allowed values follow as further messages); with gap costs that the matrix has,
/// the status is 0 and the engine goes on without the query's block (`ungapped_karlin_valid`).
/// NCBI BLAST+ 2.17.0, PS6 corpus `b.gap.bad`, `o.ok.matrix`; review probes `k.X20.badgap`,
/// `z.x20.go0.ge5`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:566-570
/// ```c++
///     if (status != 0) {
///         retval = BlastScoreBlkFree(retval);
///         if (errors) {
///             string msg(errors->message);
///             errors = Blast_MessageFree(errors);
/// ```
fn score_block(
    scoring: &MsaScoring<'_>,
    query: &[u8],
) -> Result<(ScoringMatrix, Option<PsiScoreBlk>)> {
    let matrix_name = scoring.engine_matrix_name();
    let Ok(matrix) = matrix_name.parse::<ScoringMatrix>() else {
        return Err(engine_error("Unknown error code -1"));
    };
    let spec = ProteinScoringSpec {
        matrix,
        gap_open: scoring.engine_gap_existence(),
        gap_extend: scoring.engine_gap_extension(),
    };
    if !crate::stats::tables::protein_scoring_supported(&spec) {
        if !ungapped_karlin_valid(matrix, query) {
            return Err(engine_error(
                crate::report::query_warnings::INVALID_QUERY_MESSAGE,
            ));
        }
        return Err(engine_error(&format!(
            "Gap existence and extension values of {} and {} not supported for {}",
            spec.gap_open,
            spec.gap_extend,
            matrix_name.to_ascii_uppercase()
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2786
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:514-520 (the ungapped
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

/// LOSAT's refusal of a master sequence that is the start of a longer sequence above it in the
/// file ([`QueryRows`]; DECISIONS 75.1 of line PS): `row` is that longer sequence's row, counted
/// from 1 as `-msa_master_idx` counts. The first row is its own first match, so the master of such
/// a run is one that `-msa_master_idx` names. With the master above that row, the master is the
/// first row that starts with its letters, and NCBI's alignment has its length.
fn longer_row_above_master(master_index: Option<i32>, row: usize) -> anyhow::Error {
    let master = match master_index {
        Some(index) => format!("-msa_master_idx {index}"),
        None => "its first sequence as the master".to_string(),
    };
    anyhow::anyhow!(
        "-in_msa with {master}, a master sequence that is the start of the longer sequence of row {row} above it, is not supported by LOSAT's PSIBLAST (NCBI BLAST+ 2.17.0 reads and writes past its alignment buffers there); put the master sequence above row {row} in the alignment file"
    )
}

/// LOSAT's refusal of a run whose alignment reader reached a state that no input gives (S3 of
/// the PS6 review: the reader's two preconditions on the Seq-ids, `aln_reader::verify_alignment_info`).
fn internal_error(what: &str) -> anyhow::Error {
    anyhow::anyhow!("LOSAT internal error: -in_msa: {what}")
}

/// The PSSM that NCBI makes of the alignment `input` (`x_CreatePssmFromMsa`), or its error at the
/// point of the run where NCBI raises it, or LOSAT's refusal of what NCBI crashes on or what the
/// alignment reader does not read.
///
/// The PSSM has no query Seq-entry where no Bioseq of the alignment is the query (a `*` in the
/// master).
pub(crate) fn pssm_from_msa(
    input: &MsaInput,
    scoring: &MsaScoring<'_>,
) -> Result<crate::report::serial::pssm::PssmWithParameters> {
    let bytes = crate::algorithm::pssm_search::in_pssm::read_input_argument("in_msa", &input.path)?;
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
        Err(AlnError::Internal(what)) => return Err(internal_error(what)),
    };
    let mut rows = alignment.rows.clone();
    let query = extract_query(&rows, input.master_index)?;
    let front = validate_query(&mut rows, &query)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:86-89
    // ```c++
    //     if (pssm_input_msa->GetQueryLength() == 0) {
    //         NCBI_THROW(CPssmEngineException, eNullInputData,
    //            "Query length provided by IPssmInputData is 0");
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:188-194
    // ```c++
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
    let (matrix, sbp) = score_block(scoring, &query)?;
    // `Process` comes after the score block (the matrix and gap errors above are NCBI's first):
    // it fills the engine's alignment by the first row's letters, which NCBI does past its
    // buffers where that row has more letters than the query ([`QueryRows`]). LOSAT refuses there.
    let Some(aligned) = QueryRows::new(&rows, &query) else {
        return Err(longer_row_above_master(input.master_index, front + 1));
    };
    // `Process` before the engine: the query's Bioseq is looked for last, and may stop the run.
    let msa = engine_msa(&aligned);
    let bioseq = query_bioseq(&alignment.bioseqs, &query)
        .map_err(|exception| anyhow::Error::from(crate::cli::NativeError::exception(&exception)))?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:380-396
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
            scoring.engine_matrix_name().to_ascii_uppercase()
        ),
        PsiError::NoUngappedKarlinBlock => anyhow::anyhow!(
            "-in_msa with a master sequence whose ungapped Karlin-Altschul parameters cannot be computed (a sequence of X alone, for one) is not supported by LOSAT's PSIBLAST (NCBI BLAST+ 2.17.0 stops with a crash)"
        ),
        other => engine_error(error_code_to_string(other)),
    })?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:400-405
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
        scoring.engine_matrix_name(),
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
        crate::algorithm::pssm_search::in_pssm::psi_search_and_object_of(
            pssm,
            reader,
            crate::algorithm::pssm_search::in_pssm::IDS_OF_THE_MASTER_SEQUENCE,
        )?
    };
    search.from_msa = true;
    Ok((search, file))
}

#[cfg(test)]
#[path = "in_msa_tests.rs"]
mod tests;
