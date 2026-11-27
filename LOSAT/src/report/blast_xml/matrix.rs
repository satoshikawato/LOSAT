//! The score matrix that outfmt 5 counts positives with (`CBlastFormattingMatrix`).
//!
//! This is **not** outfmt 0's matrix (`utils::matrix::protein_display_score`): outfmt 0 shows
//! `+` from `k_PSymbol` ("ARNDCQEGHILKMFPSTWYVBZX", align_format_util.cpp:88) and a J, O or U
//! pair as no positive; the XML counts every pair whose score in the packed table is above 0,
//! J included.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:185-232
//! (`kMatrixCols` = 28, the NCBIstdaa codes; `NCBISM_Unpack` gives a letter outside the
//! packed symbols the table's `defscore`, util/tables/raw_scoremat.c:81-118.)
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:392-425
//! So the matrix is indexed by IUPACaa letters; `-` and `*` have no IUPACaa letter (INT_MIN),
//! U and O are letters outside the packed symbols (BLOSUM62's `defscore` -4), and J has the
//! packed table's J row.

use crate::config::ScoringMatrix;
use crate::utils::matrix::standard_protein_matrix;
use crate::wording::ENGINE_ERROR;

/// The packed table's symbols, in the order of its rows (NCBI `s_Blosum62PSM`:
/// A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, J, Z, X, *).
const SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// The letters of the NCBIstdaa codes 1-27 that `CSeqportUtil::Convert` gives an IUPACaa
/// letter: all but `*` (code 25), and `-` (code 0).
const IUPACAA_LETTERS: &[u8; 26] = b"ABCDEFGHIKLMNPQRSTVWXYZUOJ";

/// The scores of `CBlastFormattingMatrix` for a pair of letters, `i32::MIN` for a pair with a
/// letter that has no entry.
pub struct FormattingMatrix {
    scores: Box<[i32]>,
}

/// The packed matrix of `data4xmlformat.cpp`'s `x_FillScoreMatrix` by name, compared with
/// `strcmp` (so `blosum62` has none and a protein program's XML stops with "unsupported score
/// matrix").
fn packed_matrix(name: impl AsRef<[u8]>) -> Option<ScoringMatrix> {
    let name = name.as_ref();
    ScoringMatrix::ALL
        .into_iter()
        .find(|matrix| matrix.ncbi_name().as_bytes() == name)
}

/// Whether `outfmt` is a format of the BLAST XML family, whose report data takes the formatting
/// matrix: BLAST XML (5) and the XML2 and JSON documents (16, 15, 14, 13), read as NCBI reads
/// `-outfmt` (`parse_formatting_string`).
pub fn is_structured_xml(outfmt: impl AsRef<[u8]>) -> bool {
    use crate::blastinput::app::ReportFormat;
    crate::blastinput::app::parse_formatting_string(outfmt)
        .ok()
        .and_then(|choice| ReportFormat::from_number(choice.number))
        .is_some_and(|format| {
            format == ReportFormat::Xml || crate::report::blast_xml2::Xml2Kind::of(format).is_some()
        })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2392
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2204-2212
/// Where NCBI's run stops when its BLAST XML family report data does not know the matrix name
/// (`xml_formatting_failure`): the report data, and with it the error, is made where NCBI writes
/// a query's report (NCBI BLAST+ 2.17.0, probe K/k4f of line RG, 221 runs).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XmlMatrixFailure {
    /// BLAST XML (outfmt 5), and the XML2 and JSON documents of a `-db` search of a program that
    /// does not iterate: each result's report is written as it is formatted, so the batch of the
    /// run's first searched query fails.
    FirstQuery,
    /// The XML2 and JSON documents of a `-subject` search (`m_IsBl2Seq`) or of an iterative
    /// program (`m_IsIterative`: psiblast, tblastn `-in_pssm`; `blast_xml2::reports_held`): a
    /// query's report is held until a result arrives whose query Seq-id does not match the held
    /// one (blast_format.cpp:2377); a result of the held Seq-id is accumulated with it. The report
    /// data is made, and fails, when the result of the run's first searched query whose Seq-id
    /// does not match the first's arrives (the run's query, `query_batch::QueryScopes::
    /// xml_failing_query`, the one owner of it: without `-parse_deflines` every query has a local
    /// id of its own, so the second), in its batch, after its reader messages and search strategy
    /// request and before anything of its reports. Without such a query (or a run that searched no
    /// further), the run fails at its end (`PrintEpilog`) instead of writing the epilog, unless
    /// another error ends it first. (psiblast's PHI-BLAST search, `-phi_pattern`, is not
    /// iterative: its `-db` reports are not held. tblastn's `-in_pssm` search has one query.)
    SecondQueryOrEnd,
}

/// Whether and where a protein program's run of the formats `outfmts` with the score matrix
/// `matrix_name` (as typed) stops in its formats: a BLAST XML family format
/// (`is_structured_xml`) makes report data, whose matrix is found by `strcmp` of the name
/// (`packed_matrix`). The programs give their facts: a search of subject records
/// (`subject_search`: `-subject`) and an iterative program (`iterative`: psiblast, tblastn
/// `-in_pssm`), whose XML2 and JSON reports are held (`blast_xml2::reports_held`). With several
/// formats (the library), the first to fail stops the run. The one decision behind the error
/// (`unsupported_score_matrix`) and the search strategy's last request
/// (`StrategyExport::stop_after`) of BLASTP (psiblast too), TBLASTN, TBLASTX and BLASTX; the query
/// whose batch fails is the run's (`query_batch::QueryScopes::xml_failing_query`).
pub fn xml_formatting_failure<O: AsRef<[u8]>>(
    outfmts: impl IntoIterator<Item = O>,
    matrix_name: impl AsRef<[u8]>,
    subject_search: bool,
    iterative: bool,
) -> Option<XmlMatrixFailure> {
    if packed_matrix(matrix_name).is_some() {
        return None;
    }
    let holds_back = crate::report::blast_xml2::reports_held(subject_search, iterative);
    outfmts
        .into_iter()
        .filter(|outfmt| is_structured_xml(outfmt.as_ref()))
        .map(|outfmt| {
            if holds_back && !is_blast_xml(outfmt.as_ref()) {
                XmlMatrixFailure::SecondQueryOrEnd
            } else {
                XmlMatrixFailure::FirstQuery
            }
        })
        // The batch of the run's first searched query fails before any held report is made.
        .min_by_key(|&failure| failure == XmlMatrixFailure::SecondQueryOrEnd)
}

/// Whether `outfmt` is BLAST XML (5), whose report is never held (blast_format.cpp:707-743).
fn is_blast_xml(outfmt: impl AsRef<[u8]>) -> bool {
    crate::blastinput::app::parse_formatting_string(outfmt)
        .ok()
        .and_then(|choice| crate::blastinput::app::ReportFormat::from_number(choice.number))
        == Some(crate::blastinput::app::ReportFormat::Xml)
}

/// NCBI's error when it formats a result in BLAST XML with a score matrix that its XML data does
/// not know (`xml_formatting_failure`): `CATCH_ALL`'s engine error, exit status 3.
pub fn unsupported_score_matrix() -> anyhow::Error {
    crate::cli::NativeError {
        exit: 3,
        message: format!("{ENGINE_ERROR} unsupported score matrix\n").into(),
    }
    .into()
}

impl FormattingMatrix {
    /// The matrix of the score matrix `name` (`CCmdLineBlastXMLReportData::x_FillScoreMatrix`).
    pub fn new(name: impl AsRef<[u8]>) -> Option<Self> {
        let matrix = packed_matrix(name)?;
        let mut scores = vec![i32::MIN; 256 * 256].into_boxed_slice();
        for &row in IUPACAA_LETTERS {
            for &col in IUPACAA_LETTERS {
                scores[usize::from(row) * 256 + usize::from(col)] = packed_score(matrix, row, col);
            }
        }
        Some(Self { scores })
    }

    /// The element of NCBI's `matrix` at row `query` and column `subject` (its `operator()`).
    pub fn score(&self, query: u8, subject: u8) -> i32 {
        self.scores[usize::from(query) * 256 + usize::from(subject)]
    }

    /// Whether the pair counts as a positive (`score > 0`).
    #[inline]
    pub fn positive(&self, query: u8, subject: u8) -> bool {
        self.score(query, subject) > 0
    }
}

/// `NCBISM_Unpack` of a packed matrix for two letters: the packed table, or its `defscore`
/// (BLOSUM62 -4, BLOSUM45 -5, ...) when either letter is not one of the packed symbols.
fn packed_score(matrix: ScoringMatrix, row: u8, col: u8) -> i32 {
    let table = standard_protein_matrix(matrix);
    let index = |letter: u8| SYMBOLS.iter().position(|&symbol| symbol == letter);
    match (index(row), index(col)) {
        (Some(r), Some(c)) => i32::from(table.scores[r * SYMBOLS.len() + c]),
        _ => table.defscore,
    }
}

#[cfg(test)]
#[path = "matrix_tests.rs"]
mod tests;
