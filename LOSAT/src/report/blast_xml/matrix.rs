//! The score matrix that outfmt 5 counts positives with (`CBlastFormattingMatrix`).
//!
//! This is **not** outfmt 0's matrix (`utils::matrix::protein_display_score`): outfmt 0 shows
//! `+` from `k_PSymbol` ("ARNDCQEGHILKMFPSTWYVBZX", align_format_util.cpp:88) and a J, O or U
//! pair as no positive; the XML counts every pair whose score in the packed table is above 0,
//! J included.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:185-232
//! ```c++
//!     if (strcmp(matrix_name, "BLOSUM62") == 0)
//!         packed_matrix = &NCBISM_Blosum62;
//!     ...
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

/// The packed matrix of `data4xmlformat.cpp`'s `x_FillScoreMatrix` by name (exactly as typed)
/// that LOSAT runs; none for a name it does not know (BLOSUM50, 80, 90, PAM and IDENTITY
/// are NCBI's, but no LOSAT program with outfmt 5 runs them).
fn packed_matrix(name: &str) -> Option<ScoringMatrix> {
    match name {
        "BLOSUM45" => Some(ScoringMatrix::Blosum45),
        "BLOSUM62" => Some(ScoringMatrix::Blosum62),
        _ => None,
    }
}

impl FormattingMatrix {
    /// The matrix of the score matrix `name` (`CCmdLineBlastXMLReportData::x_FillScoreMatrix`);
    /// BLOSUM62 and BLOSUM45 run in LOSAT.
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
/// (BLOSUM62 -4, BLOSUM45 -5) when either letter is not one of the packed symbols.
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
        assert!(FormattingMatrix::new("BLOSUM80").is_none());
        assert!(FormattingMatrix::new("blosum62").is_none());
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
}
