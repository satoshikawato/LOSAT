//! The score matrix that outfmt 5 counts positives with (`CBlastFormattingMatrix`).
//!
//! This is **not** outfmt 0's matrix (`utils::matrix::protein_display_score`): outfmt 0 shows
//! `+` from `k_PSymbol` ("ARNDCQEGHILKMFPSTWYVBZX", align_format_util.cpp:88) and a J, O or U
//! pair as no positive; the XML counts every pair whose score in the packed table is above 0,
//! J included.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:185-232
//! ```c++
//!     if (strcmp(matrix_name, "BLOSUM45") == 0)
//!         packed_matrix = &NCBISM_Blosum45;
//!     ...
//!     else if (strcmp(matrix_name, "IDENTITY") == 0)
//!         	packed_matrix = &NCBISM_Identity;
//!     else {
//!         string prog_name = Blast_ProgramNameFromType(
//!                                            m_Options.GetProgramType());
//!         if (prog_name != "blastn" && prog_name != "megablast") {
//!             NCBI_THROW(blast::CBlastException, eInvalidArgument,
//!                         "unsupported score matrix");
//!         }
//!     }
//!
//!     if (packed_matrix) {
//!         SNCBIFullScoreMatrix m;
//!
//!         NCBISM_Unpack(packed_matrix, &m);
//!
//!         for (unsigned int i = 0; i < kMatrixCols; i++) {
//!             for (unsigned int j = 0; j < kMatrixCols; j++) {
//!                 m_Matrix[i][j] = m.s[i][j];
//! ```
//! (`kMatrixCols` = 28, the NCBIstdaa codes; `NCBISM_Unpack` gives a letter outside the
//! packed symbols the table's `defscore`, util/tables/raw_scoremat.c:81-118.)
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:392-425
//! ```c++
//! CBlastFormattingMatrix::CBlastFormattingMatrix(int** data, unsigned int nrows,
//!                                                unsigned int ncols)
//! {
//!     const int kAsciiSize = 256;
//!     Resize(kAsciiSize, kAsciiSize, INT_MIN);
//!     ...
//!     // Convert to IUPACaa using the CSeqportUtil::Convert method.
//!     ...
//!             if (iupacaa_values[row] >= 0 && iupacaa_values[col] >= 0) {
//!                 (*this)((int)iupacaa_values[row], (int)iupacaa_values[col]) =
//!                     data[row][col];
//! ```
//! So the matrix is indexed by IUPACaa letters; `-` and `*` have no IUPACaa letter (INT_MIN),
//! U and O are letters outside the packed symbols (BLOSUM62's `defscore` -4), and J has the
//! packed table's J row.

use crate::config::ScoringMatrix;
use crate::utils::matrix::standard_protein_matrix;

/// The packed table's symbols, in the order of its rows (NCBI `s_Blosum62PSM`:
/// A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, J, Z, X, *).
const SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// The letters of the NCBIstdaa codes 1-27 that `CSeqportUtil::Convert` gives an IUPACaa
/// letter: all but `*` (code 25), and `-` (code 0).
const IUPACAA_LETTERS: &[u8; 26] = b"ABCDEFGHIKLMNPQRSTVWXYZUOJ";

/// `CBlastFormattingMatrix(letter, letter)` as scores, `i32::MIN` for a pair with a letter that
/// has no entry.
pub struct FormattingMatrix {
    scores: Box<[i32]>,
}

/// The packed matrix of `data4xmlformat.cpp`'s `x_FillScoreMatrix` by name, compared with
/// `strcmp` (so `blosum62` has none and a protein program's XML stops with "unsupported score
/// matrix").
fn packed_matrix(name: &str) -> Option<ScoringMatrix> {
    ScoringMatrix::ALL
        .into_iter()
        .find(|matrix| matrix.ncbi_name() == name)
}

/// Whether `outfmt` is a format of the BLAST XML family, whose report data takes the formatting
/// matrix: BLAST XML (5) and the XML2 and JSON documents (16, 15, 14, 13), read as NCBI reads
/// `-outfmt` (`parse_formatting_string`).
pub fn is_structured_xml(outfmt: &str) -> bool {
    use crate::blastinput::app::ReportFormat;
    crate::blastinput::app::parse_formatting_string(outfmt)
        .ok()
        .and_then(|choice| ReportFormat::from_number(choice.number))
        .is_some_and(|format| {
            format == ReportFormat::Xml || crate::report::blast_xml2::Xml2Kind::of(format).is_some()
        })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
// ```c++
//     else {
//         string prog_name = Blast_ProgramNameFromType(
//                                            m_Options.GetProgramType());
//         if (prog_name != "blastn" && prog_name != "megablast") {
//             NCBI_THROW(blast::CBlastException, eInvalidArgument,
//                         "unsupported score matrix");
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2392
// ```c++
// 	if((m_IsIterative) || (m_IsBl2Seq)) {
// 		if(m_AccumulatedResults.empty()) {
// 			...
// 			m_AccumulatedResults.push_back(res);
// 			...
// 		else {
// 			CConstRef<CSeq_id> query_id = results.GetSeqId();
// 			if(m_AccumulatedResults[0].GetSeqId()->Match(*query_id)) {
// 				m_AccumulatedResults.push_back(res);
// 			}
// 			else {
// 				CRef <CBlastSearchQuery> q = m_AccumulatedQueries->GetBlastSearchQuery(0);
// 				if(m_IsBl2Seq) {
// 				    CCmdLineBlastXML2ReportData report_data(q, m_AccumulatedResults, m_Options,
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2204-2212
// ```c++
// CBlastFormat::PrintEpilog(const blast::CBlastOptions& options)
// {
//     if ((m_FormatType == CFormattingArgs::eXml2) || (m_FormatType == CFormattingArgs::eJson) ||
//         (m_FormatType == CFormattingArgs::eXml2_S) || (m_FormatType == CFormattingArgs::eJson_S)) {
//     	if(!m_AccumulatedResults.empty()) {
//     		CRef <CBlastSearchQuery> q = m_AccumulatedQueries->GetBlastSearchQuery(0);
//     		if(m_IsBl2Seq) {
//     			CCmdLineBlastXML2ReportData report_data(q, m_AccumulatedResults, m_Options,
// ```
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
    /// program (`m_IsIterative`: psiblast, tblastn `-in_pssm`): a query's report is held until a
    /// result of the next query arrives, so the batch of the run's second searched query fails,
    /// after its reader messages and search strategy request and before anything of its reports;
    /// a run of one searched query fails at its end (`PrintEpilog`) instead of writing the
    /// epilog, unless another error ends it first.
    SecondQueryOrEnd,
}

impl XmlMatrixFailure {
    /// The run's searched query (0-based, in the run's order) whose batch fails as NCBI formats
    /// it.
    pub fn query(self) -> usize {
        match self {
            Self::FirstQuery => 0,
            Self::SecondQueryOrEnd => 1,
        }
    }

    /// Whether NCBI fails while it formats the batch of the run's searched queries `batch`.
    pub fn fails_in(self, batch: std::ops::Range<usize>) -> bool {
        batch.contains(&self.query())
    }

    /// Whether NCBI fails at the end of a run whose batches searched `searched` queries and
    /// ended without an error (its epilog makes the held report of the one query).
    pub fn fails_at_end(self, searched: usize) -> bool {
        self == Self::SecondQueryOrEnd && searched == 1
    }
}

/// Whether and where a protein program's run of the formats `outfmts` with the score matrix
/// `matrix_name` (as typed) stops in its formats: a BLAST XML family format
/// (`is_structured_xml`) makes report data, whose matrix is found by `strcmp` of the name
/// (`packed_matrix`). `holds_back` is NCBI's `m_IsBl2Seq || m_IsIterative`: a `-subject`
/// search, psiblast, tblastn `-in_pssm`. With several formats (the library), the first to fail
/// stops the run. The one decision behind the error (`unsupported_score_matrix`) and the search
/// strategy's last request (`StrategyExport::stop_after`) of BLASTP (psiblast too), TBLASTN,
/// TBLASTX and BLASTX.
pub fn xml_formatting_failure<'a>(
    outfmts: impl IntoIterator<Item = &'a str>,
    matrix_name: &str,
    holds_back: bool,
) -> Option<XmlMatrixFailure> {
    if packed_matrix(matrix_name).is_some() {
        return None;
    }
    outfmts
        .into_iter()
        .filter(|outfmt| is_structured_xml(outfmt))
        .map(|outfmt| {
            if holds_back && !is_blast_xml(outfmt) {
                XmlMatrixFailure::SecondQueryOrEnd
            } else {
                XmlMatrixFailure::FirstQuery
            }
        })
        .min_by_key(|failure| failure.query())
}

/// Whether `outfmt` is BLAST XML (5), whose report is never held (blast_format.cpp:707-743).
fn is_blast_xml(outfmt: &str) -> bool {
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
        message: "BLAST engine error: unsupported score matrix\n".to_string(),
    }
    .into()
}

impl FormattingMatrix {
    /// The matrix of the score matrix `name` (`CCmdLineBlastXMLReportData::x_FillScoreMatrix`).
    pub fn new(name: &str) -> Option<Self> {
        let matrix = packed_matrix(name)?;
        let mut scores = vec![i32::MIN; 256 * 256].into_boxed_slice();
        for &row in IUPACAA_LETTERS {
            for &col in IUPACAA_LETTERS {
                scores[usize::from(row) * 256 + usize::from(col)] = packed_score(matrix, row, col);
            }
        }
        Some(Self { scores })
    }

    /// `(*matrix)(query, subject)`.
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
mod tests {
    use super::*;

    /// The stop of a run's XML formats (BLAST XML, the XML2 and JSON documents) with a matrix name
    /// that NCBI's formatting table does not spell: every other format, and every name it spells,
    /// formats.
    #[test]
    fn the_xml_formats_fail_for_a_name_ncbi_does_not_spell() {
        for outfmt in ["5", "13", "14", "15", "16", " 5"] {
            assert!(is_structured_xml(outfmt), "{outfmt}");
            for holds_back in [false, true] {
                let failure = xml_formatting_failure([outfmt], "blosum62", holds_back);
                let held = holds_back && outfmt.trim() != "5";
                assert_eq!(
                    failure,
                    Some(if held {
                        XmlMatrixFailure::SecondQueryOrEnd
                    } else {
                        XmlMatrixFailure::FirstQuery
                    }),
                    "{outfmt} {holds_back}"
                );
                assert_eq!(
                    xml_formatting_failure([outfmt], "BLOSUM62", holds_back),
                    None,
                    "{outfmt}"
                );
            }
        }
        for outfmt in ["0", "6", "7 qseqid", "8", "11", "12", "17", "abc"] {
            assert!(!is_structured_xml(outfmt), "{outfmt}");
            assert_eq!(
                xml_formatting_failure([outfmt], "blosum62", true),
                None,
                "{outfmt}"
            );
        }
        assert_eq!(
            xml_formatting_failure(["6", "16"], "blosum62", false),
            Some(XmlMatrixFailure::FirstQuery)
        );
        // Several formats: the first to fail stops the run.
        assert_eq!(
            xml_formatting_failure(["16", "15"], "blosum62", true),
            Some(XmlMatrixFailure::SecondQueryOrEnd)
        );
        assert_eq!(
            xml_formatting_failure(["16", "5"], "blosum62", true),
            Some(XmlMatrixFailure::FirstQuery)
        );
        let error = unsupported_score_matrix();
        let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
        assert_eq!(
            (native.exit, native.message.as_str()),
            (3, "BLAST engine error: unsupported score matrix\n")
        );
    }

    /// Where the run stops (probe K/k4f, NCBI 2.17.0): the batch of the first searched query, or
    /// of the second when the report is held (`-subject` 13-16: 3 batches of one query wrote
    /// the requests of Query_1 and Query_2), or the end of a run of one query.
    #[test]
    fn a_held_report_fails_with_the_second_query_or_at_the_end() {
        let first = XmlMatrixFailure::FirstQuery;
        let held = XmlMatrixFailure::SecondQueryOrEnd;
        assert!(first.fails_in(0..1) && first.fails_in(0..3) && !first.fails_in(1..2));
        assert!(!held.fails_in(0..1) && held.fails_in(1..2) && held.fails_in(0..3));
        assert!(!held.fails_in(2..3) && !held.fails_in(1..1));
        assert!(held.fails_at_end(1) && !held.fails_at_end(0) && !held.fails_at_end(2));
        assert!(!first.fails_at_end(1));
    }

    #[test]
    fn engine_table_is_the_packed_ncbi_table_with_its_symbol_order() {
        // sm_blosum62.c: the A row, the J row and the V column.
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'A', b'A'), 4);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'A', b'R'), -1);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'J', b'I'), 3);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'J', b'V'), 2);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'J', b'L'), 3);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'*', b'*'), 1);
        assert_eq!(packed_score(ScoringMatrix::Blosum62, b'*', b'A'), -4);
    }

    #[test]
    fn positives_of_rare_letters_follow_the_packed_table() {
        let matrix = FormattingMatrix::new("BLOSUM62").expect("BLOSUM62");
        // J is a packed symbol: J/V and I/J are positive in the XML (and not in outfmt 0).
        assert!(matrix.positive(b'J', b'V'));
        assert!(matrix.positive(b'I', b'J'));
        assert!(matrix.positive(b'L', b'J'));
        assert!(!matrix.positive(b'J', b'A'));
        // U and O are not packed symbols: `defscore` -4, so never positive unless equal
        // (equal letters are counted before the matrix is read).
        assert_eq!(matrix.score(b'U', b'C'), -4);
        assert!(!matrix.positive(b'U', b'C'));
        assert!(!matrix.positive(b'O', b'X'));
        assert_eq!(matrix.score(b'U', b'U'), -4);
        // '-' and '*' have no IUPACaa letter.
        assert_eq!(matrix.score(b'*', b'*'), i32::MIN);
        assert_eq!(matrix.score(b'-', b'A'), i32::MIN);
        assert!(!matrix.positive(b'*', b'A'));
        // ordinary pairs: I/V +3, B/D +4, Z/E +4, X/X -1, E/D +2.
        assert!(matrix.positive(b'I', b'V'));
        assert!(matrix.positive(b'B', b'D'));
        assert!(matrix.positive(b'Z', b'E'));
        assert!(!matrix.positive(b'X', b'A'));
        assert!(matrix.positive(b'E', b'D'));
        assert!(!matrix.positive(b'A', b'W'));
    }

    #[test]
    fn the_core_agrees_with_the_engine_matrix() {
        let matrix = FormattingMatrix::new("BLOSUM62").expect("BLOSUM62");
        for (r, &row) in SYMBOLS.iter().enumerate() {
            for (c, &col) in SYMBOLS.iter().enumerate() {
                if row == b'*' || col == b'*' {
                    continue;
                }
                assert_eq!(
                    matrix.score(row, col),
                    i32::from(crate::utils::matrix::BLOSUM62[r * 25 + c]),
                    "{} {}",
                    row as char,
                    col as char
                );
            }
        }
        assert!(FormattingMatrix::new("blosum62").is_none());
        assert!(FormattingMatrix::new("BLOSUM62_20").is_none());
        assert!(FormattingMatrix::new("IDENTITYX").is_none());
    }

    #[test]
    fn blosum45_is_the_packed_table_with_defscore_minus_five() {
        // sm_blosum45.c: A/A 5, A/R -2, J/I 4, J/V 3, the `*` diagonal 1, defscore -5.
        let blosum45 = ScoringMatrix::Blosum45;
        assert_eq!(packed_score(blosum45, b'A', b'A'), 5);
        assert_eq!(packed_score(blosum45, b'A', b'R'), -2);
        assert_eq!(packed_score(blosum45, b'J', b'I'), 4);
        assert_eq!(packed_score(blosum45, b'*', b'*'), 1);
        assert_eq!(packed_score(blosum45, b'*', b'A'), -5);
        let matrix = FormattingMatrix::new("BLOSUM45").expect("BLOSUM45");
        assert_eq!(matrix.score(b'U', b'C'), -5);
        assert_eq!(matrix.score(b'A', b'A'), 5);
        assert!(matrix.positive(b'I', b'V'));
        assert!(!matrix.positive(b'A', b'R'));
    }

    #[test]
    fn all_nine_ncbi_matrices_have_a_formatting_matrix() {
        // data4xmlformat.cpp:195-212: the nine built-in matrices by their exact names.
        for matrix in ScoringMatrix::ALL {
            let formatting = FormattingMatrix::new(matrix.ncbi_name()).expect("built in");
            let table = standard_protein_matrix(matrix);
            assert_eq!(formatting.score(b'A', b'A'), i32::from(table.scores[0]));
            assert_eq!(formatting.score(b'U', b'A'), table.defscore, "{matrix:?}");
        }
        // sm_identity.c: 9 on the diagonal but X/X, -5 elsewhere and defscore -5.
        let identity = FormattingMatrix::new("IDENTITY").expect("IDENTITY");
        assert_eq!(identity.score(b'W', b'W'), 9);
        assert_eq!(identity.score(b'X', b'X'), -5);
        assert_eq!(identity.score(b'I', b'V'), -5);
        assert!(!identity.positive(b'I', b'V'));
        // sm_pam30.c: the A row starts 6, -7, -4; defscore -17.
        let pam30 = FormattingMatrix::new("PAM30").expect("PAM30");
        assert_eq!(
            (
                pam30.score(b'A', b'A'),
                pam30.score(b'A', b'R'),
                pam30.score(b'A', b'N')
            ),
            (6, -7, -4)
        );
        assert_eq!(pam30.score(b'O', b'A'), -17);
    }
}
