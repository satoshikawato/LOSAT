//! NCBI's option functions over the protein matrix tables (`protein_tables.rs`): the
//! table lookup of `Blast_KarlinBlkGappedLoadFromTables`, its two error messages, the
//! preferred gap costs of a matrix and the suggested word threshold and window.

use super::protein_tables::{ProteinMatrixTable, ProteinTableRow, PROTEIN_MATRIX_TABLES};
use crate::config::ScoringMatrix;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
pub fn blast_nint(x: f64) -> i64 {
    let x = x + if x >= 0.0 { 0.5 } else { -0.5 };
    x as i64
}

/// The table of the matrix `name` (case-insensitive, `strcasecmp`), IDENTITY only when
/// `standard_only` is false.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3593-3605
pub fn matrix_table(name: &str, standard_only: bool) -> Option<&'static ProteinMatrixTable> {
    let tables = if standard_only {
        &PROTEIN_MATRIX_TABLES[..PROTEIN_MATRIX_TABLES.len() - 1]
    } else {
        &PROTEIN_MATRIX_TABLES[..]
    };
    tables
        .iter()
        .find(|table| table.name.eq_ignore_ascii_case(name))
}

/// The table of `matrix`: the one whose name is `matrix.ncbi_name()` (`strcasecmp`), found by
/// its index in `PROTEIN_MATRIX_TABLES` instead of the name.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2958-2992 (the order of the tables)
pub fn scoring_matrix_table(matrix: ScoringMatrix) -> &'static ProteinMatrixTable {
    let index = match matrix {
        ScoringMatrix::Blosum80 => 0,
        ScoringMatrix::Blosum62 => 1,
        ScoringMatrix::Blosum50 => 2,
        ScoringMatrix::Blosum45 => 3,
        ScoringMatrix::Pam250 => 4,
        ScoringMatrix::Blosum90 => 5,
        ScoringMatrix::Pam30 => 6,
        ScoringMatrix::Pam70 => 7,
        ScoringMatrix::Identity => 8,
    };
    &PROTEIN_MATRIX_TABLES[index]
}

/// The first row of `table` with the gap costs, the search of the table loaders.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3612-3616
pub fn gap_costs_row(
    table: &'static ProteinMatrixTable,
    gap_open: i32,
    gap_extend: i32,
) -> Option<&'static ProteinTableRow> {
    table.values.iter().find(|row| {
        blast_nint(row.0[0]) == i64::from(gap_open) && blast_nint(row.0[1]) == i64::from(gap_extend)
    })
}

/// Why `Blast_KarlinBlkGappedLoadFromTables` fails: 1, the matrix is not in the tables;
/// 2, the gap costs are not in the matrix's table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableStatus {
    UnknownMatrix = 1,
    UnknownGapCosts = 2,
}

/// The table row of the gap costs of the matrix `name`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3608-3641
pub fn karlin_blk_gapped_load_from_tables(
    gap_open: i32,
    gap_extend: i32,
    name: &str,
    standard_only: bool,
) -> Result<&'static ProteinTableRow, TableStatus> {
    let table = matrix_table(name, standard_only).ok_or(TableStatus::UnknownMatrix)?;
    gap_costs_row(table, gap_open, gap_extend).ok_or(TableStatus::UnknownGapCosts)
}

/// `snprintf` into what is left of a fixed buffer: at most `size - 1` bytes are written
/// (none when `size` is 0), and the would-be length is returned.
fn snprintf(buffer: &mut Vec<u8>, size: i32, text: &str) -> i32 {
    if size > 0 {
        let room = (size - 1) as usize;
        buffer.extend_from_slice(&text.as_bytes()[..text.len().min(room)]);
    }
    text.len() as i32
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3760-3788
pub fn print_matrix_message(matrix_name: &str, standard_only: bool) -> String {
    let mut buffer = Vec::new();
    let mut buffer_sz: i32 = 1024;
    let out_sz = snprintf(
        &mut buffer,
        buffer_sz,
        &format!("{matrix_name} is not a supported matrix, supported matrices are:\n"),
    );
    buffer_sz = (buffer_sz - (1 + out_sz)).max(0);
    let count = if standard_only {
        PROTEIN_MATRIX_TABLES.len() - 1
    } else {
        PROTEIN_MATRIX_TABLES.len()
    };
    for table in &PROTEIN_MATRIX_TABLES[..count] {
        let out_sz = snprintf(&mut buffer, buffer_sz, &format!("{} \n", table.name));
        buffer_sz = (buffer_sz - out_sz).max(0);
    }
    String::from_utf8_lossy(&buffer).into_owned()
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3792-3839
pub fn print_allowed_values(matrix_name: &str, gap_open: i32, gap_extend: i32) -> String {
    let mut buffer = Vec::new();
    let mut buffer_sz: i32 = 2048;
    let out_sz = snprintf(
        &mut buffer,
        buffer_sz,
        &format!(
            "Gap existence and extension values of {gap_open} and {gap_extend} not supported for {matrix_name}\nsupported values are:\n"
        ),
    );
    buffer_sz = (buffer_sz - (1 + out_sz)).max(0);
    if let Some(table) = matrix_table(matrix_name, false) {
        for row in table.values {
            let line = if blast_nint(row.0[2]) == 32767 {
                format!("{}, {}\n", blast_nint(row.0[0]), blast_nint(row.0[1]))
            } else {
                format!(
                    "{}, {}, {}\n",
                    blast_nint(row.0[0]),
                    blast_nint(row.0[1]),
                    blast_nint(row.0[2])
                )
            };
            let out_sz = snprintf(&mut buffer, buffer_sz, &line);
            buffer_sz = (buffer_sz - out_sz).max(0);
        }
    }
    String::from_utf8_lossy(&buffer).into_owned()
}

/// The best gap costs of the matrix `name`, the row whose preference is
/// `BLAST_MATRIX_BEST` after the first (ungapped) row; `None` for a matrix that is not in
/// the tables (NCBI then leaves its caller's values, 0 and 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3374-3396
/// `Blast_GetMatrixValues` converts the costs with `(Int4)` (blast_stat.c:3068-3070).
pub fn protein_gap_existence_extend_params(name: &str) -> Option<(i32, i32)> {
    let table = matrix_table(name, false)?;
    let mut costs = None;
    for (row, &pref) in table.values.iter().zip(table.prefs).skip(1) {
        if pref == 2 {
            costs = Some((row.0[0] as i32, row.0[1] as i32));
            break;
        }
    }
    costs
}

/// The program family that the suggested threshold depends on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SuggestionProgram {
    /// blastp: neither sequence is translated.
    Protein,
    /// tblastn, tblastx: the subject is translated.
    TranslatedSubject,
    /// blastx: only the query is translated.
    TranslatedQuery,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1174-1209
pub fn suggested_threshold(program: SuggestionProgram, matrix_name: &str) -> f64 {
    let is = |name: &str| matrix_name.eq_ignore_ascii_case(name);
    let mut threshold = if is("BLOSUM62") {
        11.0
    } else if is("BLOSUM45") {
        14.0
    } else if is("BLOSUM62_20") {
        100.0
    } else if is("BLOSUM80") {
        12.0
    } else if is("PAM30") {
        16.0
    } else if is("PAM70") {
        14.0
    } else if is("IDENTITY") {
        27.0
    } else {
        11.0
    };
    match program {
        SuggestionProgram::TranslatedSubject => threshold += 2.0,
        SuggestionProgram::TranslatedQuery => threshold += 1.0,
        SuggestionProgram::Protein => {}
    }
    threshold
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1211-1236
pub fn suggested_window_size(matrix_name: &str) -> i32 {
    let is = |name: &str| matrix_name.eq_ignore_ascii_case(name);
    if is("BLOSUM45") {
        60
    } else if is("BLOSUM80") {
        25
    } else if is("PAM30") {
        15
    } else if is("PAM70") {
        20
    } else {
        40
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI BLAST+ 2.17.0 `blastp -matrix FOO` (BlastScoringOptionsValidate status 1).
    #[test]
    fn matrix_message_lists_ncbis_matrices() {
        assert_eq!(
            print_matrix_message("FOO", false),
            "FOO is not a supported matrix, supported matrices are:\nBLOSUM80 \nBLOSUM62 \nBLOSUM50 \nBLOSUM45 \nPAM250 \nBLOSUM90 \nPAM30 \nPAM70 \nIDENTITY \n"
        );
        assert!(!print_matrix_message("FOO", true).contains("IDENTITY"));
        // The 1024-byte buffer: a long name truncates the message.
        let long = "M".repeat(2000);
        assert_eq!(print_matrix_message(&long, false).len(), 1023);
    }

    // NCBI BLAST+ 2.17.0 `blastp -matrix PAM30 -gapopen 1 -gapextend 1` (status 2).
    #[test]
    fn allowed_values_list_the_table_in_order() {
        assert_eq!(
            print_allowed_values("PAM30", 1, 1),
            "Gap existence and extension values of 1 and 1 not supported for PAM30\nsupported values are:\n32767, 32767\n7, 2\n6, 2\n5, 2\n10, 1\n9, 1\n8, 1\n15, 3\n14, 2\n14, 1\n13, 3\n"
        );
    }

    #[test]
    fn tables_are_found_without_case_and_with_identity_only_for_nonstandard() {
        assert!(karlin_blk_gapped_load_from_tables(11, 1, "blosum62", false).is_ok());
        assert_eq!(
            karlin_blk_gapped_load_from_tables(11, 1, "IDENTITY", true),
            Err(TableStatus::UnknownMatrix)
        );
        assert_eq!(
            karlin_blk_gapped_load_from_tables(1, 1, "BLOSUM62", false),
            Err(TableStatus::UnknownGapCosts)
        );
        assert!(karlin_blk_gapped_load_from_tables(32767, 32767, "PAM70", false).is_ok());
    }

    #[test]
    fn preferred_gap_costs_are_the_best_rows() {
        let costs: Vec<_> = [
            "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250",
            "IDENTITY",
        ]
        .iter()
        .map(|name| protein_gap_existence_extend_params(name).unwrap())
        .collect();
        assert_eq!(
            costs,
            vec![
                (14, 2),
                (13, 2),
                (11, 1),
                (10, 1),
                (10, 1),
                (9, 1),
                (10, 1),
                (15, 2),
                (15, 2)
            ]
        );
        assert_eq!(protein_gap_existence_extend_params("FOO"), None);
    }

    #[test]
    fn suggestions_follow_blast_options_c() {
        assert_eq!(
            suggested_threshold(SuggestionProgram::Protein, "pam30"),
            16.0
        );
        assert_eq!(
            suggested_threshold(SuggestionProgram::TranslatedSubject, "BLOSUM62"),
            13.0
        );
        assert_eq!(suggested_threshold(SuggestionProgram::Protein, "FOO"), 11.0);
        assert_eq!(suggested_window_size("blosum80"), 25);
        assert_eq!(suggested_window_size("IDENTITY"), 40);
    }
}
