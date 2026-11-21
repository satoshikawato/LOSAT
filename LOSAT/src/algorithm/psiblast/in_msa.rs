//! psiblast `-in_msa F`, `-msa_master_idx N` and `-ignore_msa_master` (series PS, item PS6): the
//! PSSM that NCBI's option processing makes from an alignment file, and the search of it as
//! `-in_pssm` searches a PSSM file.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1252-1271
//! The reader is `blastinput::aln_reader`; this module is `CPsiBlastInputClustalW`
//! (`x_ExtractQueryFromMsa`, `x_ValidateQueryInMsa`, `Process`) and the use of `CPssmEngine` on
//! it (the score block of its constructor, then `Run`), with the PSSM engine of `core::psi`. The
//! object that comes out is searched by the one PSSM search (`pssm_search::in_pssm`), with the
//! rounds, files, strategy and refusals of `-in_pssm`.

use std::path::PathBuf;

use anyhow::Result;

use crate::blastinput::aln_reader::{read_alignment, Alignment, AlnBioseq, AlnError, Mol, SetMol};
use crate::blastinput::app::{engine_error, query_error};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::core::psi::{
    create_pssm_with_options, KarlinBlk, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiOptions,
    PsiScoreBlk,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::pssm::{PssmQuery, SeqData, SeqEntry, SeqInst, Seqdesc};
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
// NCBI reference (598d8ae6): c++/include/algo/blast/api/pssm_input.hpp:68-80
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:77
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:84-85
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:92-93
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

/// LOSAT's refusal of an alignment that NCBI's reader or its PSSM input takes where the result
/// is undefined, or that carries what LOSAT does not write: `cause` says what.
fn refused_cause(cause: &str) -> anyhow::Error {
    anyhow::anyhow!("-in_msa with this alignment is not supported by LOSAT's PSIBLAST: {cause}")
}

/// Posts the warnings of the reader (`ERR_POST(Warning << ..)`, prefix `psiblast`) to the
/// diagnostics stream.
fn post_warnings(diagnostics: &mut dyn std::io::Write, warnings: &[Vec<u8>]) -> Result<()> {
    for message in warnings {
        diagnostics.write_all(&crate::ncbi_exception::post_line(
            crate::ncbi_exception::Severity::Warning,
            None,
            Some("psiblast"),
            message,
        ))?;
    }
    Ok(())
}

/// The letters of a Bioseq as NCBIstdaa after `CSeqportUtil::Convert` from IUPACaa
/// (`CIupacaaToStdaa`): the letters A to Z; any other byte (`*` among them) is 0xff.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1686
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
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:225
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
fn validate_query(rows: &mut [Vec<u8>], query: &[u8]) -> Result<usize> {
    let masking = crate::algorithm::pssm_search::in_pssm::NCBISTDAA_TO_AMINOACID[21];
    // The columns that NCBI looks at and swaps are those of the first row of the file
    // (`kAligmentLength`); a row of another length is read as far as it goes (where NCBI reads
    // the end of the string, a NUL, which is no match).
    let window = rows[0].len();
    let found = rows.iter().position(|row| {
        let mut query_index = 0;
        for &residue in row.iter().take(window) {
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
        Some(0) => Ok(0),
        Some(index) => {
            // NCBI swaps the first `window` letters of the two rows by index: a row with fewer
            // letters is read and written past its end.
            if rows[index].len() < window {
                return Err(refused_cause(SHORT_ROW_SWAPPED));
            }
            let (head, tail) = rows.split_at_mut(index);
            head[0][..window].swap_with_slice(&mut tail[0][..window]);
            Ok(index)
        }
        None => Err(blast_options_error(
            "No sequence in the multiple sequence alignment provided matches the query sequence",
        )),
    }
}

/// The cause of the refusal of a row that holds the query and is shorter than the first row of the
/// file: NCBI swaps it with the first row column by column through the first row's length, past
/// the end of the shorter row's string.
const SHORT_ROW_SWAPPED: &str = "a row that holds the query is shorter than the first row of the file, and NCBI BLAST+ 2.17.0 swaps the two rows column by column through the first row's length, reading and writing past the shorter row's end";

/// The cause of the refusal of a row shorter than the first row where the first row has a
/// letter beyond the end of the shorter one: NCBI reads that row's string past its end.
pub(crate) const SHORT_ROW_READ: &str = "a row of the alignment is shorter than the first row and the first row has a letter past its end, and NCBI BLAST+ 2.17.0 reads the shorter row's string past its end for it";

/// The first row after the swap that has a letter of the first row past the end of its own
/// string plus one (the string's terminator is a readable NUL, the bytes after it are not).
fn short_row_read(front: &[u8], others: &[Vec<u8>]) -> bool {
    others.iter().any(|row| {
        front
            .iter()
            .skip(row.len() + 1)
            .any(|&residue| residue != b'-')
    })
}

/// The rows after `x_ValidateQueryInMsa`'s swap whose first row has exactly the query's letters:
/// the only input of [`engine_msa`], and made only by [`QueryRows::new`], so the engine's
/// alignment is never filled from a first row with more letters than the query.
///
/// The first row is the first one whose letters start with the query's, so it can have more:
/// with `-msa_master_idx`, a master that is the start of a sequence above it in the file. NCBI
/// then fills its alignment (`query_length` cells a row) by that longer row's letters and reads
/// past the query and writes past every row (the two functions cited at [`engine_msa`]; NCBI
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:335-396
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
/// The ungapped calculation writes its message for a query whose parameters it cannot compute
/// and its status is then replaced by the gapped calculation's:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2783-2790
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3557-3562
/// The engine throws the first message that was written, so with gap costs that the matrix does
/// not have the error is the ungapped calculation's for such a query and the gap costs' for any
/// other (whose allowed values follow as further messages); with gap costs that the matrix has,
/// the status is 0 and the engine goes on without the query's block (`ungapped_karlin_valid`).
/// NCBI BLAST+ 2.17.0, PS6 corpus `b.gap.bad`, `o.ok.matrix`; review probes `k.X20.badgap`,
/// `z.x20.go0.ge5`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:566-570
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:514-520 (the ungapped
/// calculation's status is replaced by the gapped one's, so the failure goes unreported and the
/// run goes on with `kbp_psi[0]` null; `s_PSISavePssm` and the scaling read it: NCBI BLAST+ 2.17.0
/// ends with SIGSEGV, PS6 corpus `r.ch.q.allX`, `x.in.allx.60`, `x.in.allO`)
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
    diagnostics: &mut dyn std::io::Write,
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
        Err(AlnError::Modifier(failure)) => {
            post_warnings(diagnostics, &failure.warnings)?;
            return Err(query_error(&String::from_utf8_lossy(&failure.message)));
        }
        Err(AlnError::ModifierException(warnings, exception)) => {
            post_warnings(diagnostics, &warnings)?;
            return Err(anyhow::Error::from(crate::cli::NativeError::exception(
                &exception,
            )));
        }
        Err(AlnError::StdException(what)) => {
            return Err(crate::cli::NativeError {
                exit: 255,
                message: format!("Error: {what}\n").into_bytes(),
            }
            .into())
        }
        Err(AlnError::Undefined(cause)) => return Err(refused_cause(cause)),
        Err(AlnError::Exception(exception)) => {
            return Err(anyhow::Error::from(crate::cli::NativeError::exception(
                &exception,
            )))
        }
        Err(AlnError::Internal(what)) => return Err(internal_error(what)),
    };
    post_warnings(diagnostics, &alignment.warnings)?;
    let mut rows = alignment.rows.clone();
    let query = extract_query(&rows, input.master_index)?;
    let front = validate_query(&mut rows, &query)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:86-89
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:188-194
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
    if short_row_read(aligned.front, aligned.others) {
        return Err(refused_cause(SHORT_ROW_READ));
    }
    // `Process` before the engine: the query's Bioseq is looked for last, and may stop the run.
    let msa = engine_msa(&aligned);
    let bioseq = query_bioseq(&alignment.bioseqs, &query)
        .map_err(|exception| anyhow::Error::from(crate::cli::NativeError::exception(&exception)))?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:380-396
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
    // The Bioseq is the alignment's own: its ids, its letters as IUPACaa, its title (the text of the
    // defline outside its modifiers).
    let entry = bioseq.map(|bioseq| {
        let mut descr: Vec<Seqdesc> = bioseq
            .content
            .descriptors()
            .into_iter()
            .map(Seqdesc::Other)
            .collect();
        if let Some(title) = &bioseq.title {
            descr.push(Seqdesc::Title(title.clone()));
        }
        let mol = match bioseq.set_mol {
            None | Some(SetMol::Aa) => ("aa", 3),
            Some(SetMol::Dna) => ("dna", 1),
            Some(SetMol::Rna) => ("rna", 2),
            Some(SetMol::Na) => ("na", 4),
            Some(SetMol::Other) => ("other", 255),
        };
        let mut inst = SeqInst::new(
            ("raw", 2),
            mol,
            Some(bioseq.letters.len() as u32),
            Some(SeqData::Iupacaa(bioseq.letters.clone())),
        );
        inst.hist = bioseq.content.hist();
        inst.topology = bioseq.content.topology.map(|name| match name {
            "linear" => 1,
            "circular" => 2,
            "tandem" => 3,
            _ => 255,
        });
        inst.strand = bioseq.content.strand.map(|name| match name {
            "ss" => 1,
            "ds" => 2,
            "mixed" => 3,
            _ => 255,
        });
        SeqEntry::Seq(Box::new(PssmQuery {
            ids: bioseq.ids.clone(),
            descr: (!descr.is_empty()).then_some(descr),
            inst,
            annot: bioseq.content.annot(&bioseq.ids),
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
    diagnostics: &mut dyn std::io::Write,
) -> Result<(
    crate::algorithm::pssm_search::PsiSearch,
    crate::algorithm::pssm_search::in_pssm::InPssmFile,
)> {
    let pssm = pssm_from_msa(input, scoring, diagnostics)?;
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
