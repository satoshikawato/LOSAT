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
pub fn matrix_table(
    name: impl AsRef<[u8]>,
    standard_only: bool,
) -> Option<&'static ProteinMatrixTable> {
    let name = name.as_ref();
    let tables = if standard_only {
        &PROTEIN_MATRIX_TABLES[..PROTEIN_MATRIX_TABLES.len() - 1]
    } else {
        &PROTEIN_MATRIX_TABLES[..]
    };
    tables
        .iter()
        .find(|table| table.name.as_bytes().eq_ignore_ascii_case(name))
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
    name: impl AsRef<[u8]>,
    standard_only: bool,
) -> Result<&'static ProteinTableRow, TableStatus> {
    let table = matrix_table(name, standard_only).ok_or(TableStatus::UnknownMatrix)?;
    gap_costs_row(table, gap_open, gap_extend).ok_or(TableStatus::UnknownGapCosts)
}

/// `snprintf` into what is left of a fixed buffer: at most `size - 1` bytes are written
/// (none when `size` is 0), and the would-be length is returned.
fn snprintf(buffer: &mut Vec<u8>, size: i32, text: &[u8]) -> i32 {
    if size > 0 {
        let room = (size - 1) as usize;
        buffer.extend_from_slice(&text[..text.len().min(room)]);
    }
    text.len() as i32
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3760-3788
/// The message's bytes: the buffer can end inside a character of the name (NCBI BLAST+ 2.17.0,
/// probe R9 of line LIB-B: 1,000 `a` and 30 `é` end with the byte 0xC3).
pub fn print_matrix_message(matrix_name: impl AsRef<[u8]>, standard_only: bool) -> Vec<u8> {
    let mut buffer = Vec::new();
    let mut buffer_sz: i32 = 1024;
    let first = [
        matrix_name.as_ref(),
        b" is not a supported matrix, supported matrices are:\n",
    ]
    .concat();
    let out_sz = snprintf(&mut buffer, buffer_sz, &first);
    buffer_sz = (buffer_sz - (1 + out_sz)).max(0);
    let count = if standard_only {
        PROTEIN_MATRIX_TABLES.len() - 1
    } else {
        PROTEIN_MATRIX_TABLES.len()
    };
    for table in &PROTEIN_MATRIX_TABLES[..count] {
        let line = format!("{} \n", table.name);
        let out_sz = snprintf(&mut buffer, buffer_sz, line.as_bytes());
        buffer_sz = (buffer_sz - out_sz).max(0);
    }
    buffer
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3792-3839
/// The message's bytes ([`print_matrix_message`]).
pub fn print_allowed_values(
    matrix_name: impl AsRef<[u8]>,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<u8> {
    let matrix_name = matrix_name.as_ref();
    let mut buffer = Vec::new();
    let mut buffer_sz: i32 = 2048;
    let first = [
        format!(
            "Gap existence and extension values of {gap_open} and {gap_extend} not supported for "
        )
        .as_bytes(),
        matrix_name,
        b"\nsupported values are:\n",
    ]
    .concat();
    let out_sz = snprintf(&mut buffer, buffer_sz, &first);
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
            let out_sz = snprintf(&mut buffer, buffer_sz, line.as_bytes());
            buffer_sz = (buffer_sz - out_sz).max(0);
        }
    }
    buffer
}

/// The best gap costs of the matrix `name`, the row whose preference is
/// `BLAST_MATRIX_BEST` after the first (ungapped) row; `None` for a matrix that is not in
/// the tables (NCBI then leaves its caller's values, 0 and 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3374-3396
/// `Blast_GetMatrixValues` converts the costs with `(Int4)` (blast_stat.c:3068-3070).
pub fn protein_gap_existence_extend_params(name: impl AsRef<[u8]>) -> Option<(i32, i32)> {
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
pub fn suggested_threshold(program: SuggestionProgram, matrix_name: impl AsRef<[u8]>) -> f64 {
    let matrix_name = matrix_name.as_ref();
    let is = |name: &str| matrix_name.eq_ignore_ascii_case(name.as_bytes());
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
pub fn suggested_window_size(matrix_name: impl AsRef<[u8]>) -> i32 {
    let matrix_name = matrix_name.as_ref();
    let is = |name: &str| matrix_name.eq_ignore_ascii_case(name.as_bytes());
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
#[path = "protein_options_tests.rs"]
mod tests;
