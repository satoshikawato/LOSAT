//! NCBI's option functions over the protein matrix tables (`protein_tables.rs`): the
//! table lookup of `Blast_KarlinBlkGappedLoadFromTables`, its two error messages, the
//! preferred gap costs of a matrix and the suggested word threshold and window.

use super::protein_tables::{ProteinMatrixTable, ProteinTableRow, PROTEIN_MATRIX_TABLES};

/// NCBI reference: c++/src/algo/blast/core/ncbi_math.c:437-441
/// ```c
/// long BLAST_Nint(double x)
/// {
///    x += (x >= 0. ? 0.5 : -0.5);
///    return (long)x;
/// }
/// ```
pub fn blast_nint(x: f64) -> i64 {
    let x = x + if x >= 0.0 { 0.5 } else { -0.5 };
    x as i64
}

/// The table of the matrix `name` (case-insensitive, `strcasecmp`), IDENTITY only when
/// `standard_only` is false.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3593-3605
/// ```c
///    vnp = head = BlastLoadMatrixValues(standard_only);
///    while (vnp)
///    {
///       matrix_info = vnp->ptr;
///       if (strcasecmp(matrix_info->name, matrix_name) == 0)
///       {
///          values = matrix_info->values;
///          max_number_values = matrix_info->max_number_values;
///          found_matrix = TRUE;
///          break;
///       }
///       vnp = vnp->next;
///    }
/// ```
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

/// Why `Blast_KarlinBlkGappedLoadFromTables` fails: 1, the matrix is not in the tables;
/// 2, the gap costs are not in the matrix's table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableStatus {
    UnknownMatrix = 1,
    UnknownGapCosts = 2,
}

/// The table row of the gap costs of the matrix `name`.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3608-3641
/// ```c
///    if (found_matrix)
///    {
///                 Boolean found_values=FALSE;
///            Int4 index;
///       for (index=0; index<max_number_values; index++)
///       {
///          if (BLAST_Nint(values[index][0]) == gap_open &&
///             BLAST_Nint(values[index][1]) == gap_extend)
///          {
/// ...
///             found_values = TRUE;
///             break;
///          }
///       }
///
///       if (found_values == TRUE)
///       {
///          status = 0;
///       }
///       else
///       {
///          status = 2;
///       }
///    }
///    else
///    {
///       status = 1;
///    }
/// ```
pub fn karlin_blk_gapped_load_from_tables(
    gap_open: i32,
    gap_extend: i32,
    name: &str,
    standard_only: bool,
) -> Result<&'static ProteinTableRow, TableStatus> {
    let table = matrix_table(name, standard_only).ok_or(TableStatus::UnknownMatrix)?;
    table
        .values
        .iter()
        .find(|row| {
            blast_nint(row.0[0]) == i64::from(gap_open)
                && blast_nint(row.0[1]) == i64::from(gap_extend)
        })
        .ok_or(TableStatus::UnknownGapCosts)
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

/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3760-3788
/// ```c
/// #define BUF_SZ_1024 (1024)
///    int buffer_sz= BUF_SZ_1024 ;
/// ...
///    out_sz = snprintf(ptr, (size_t)buffer_sz, "%s is not a supported matrix, supported matrices are:\n", matrix_name);
///    buffer_sz -= ( 1 + out_sz ); if ( buffer_sz < 0 ) buffer_sz = 0; // decrease out buffer size first, "+1" for null termination , once.
///
///    ptr += strlen(ptr);
///
///         vnp = head = BlastLoadMatrixValues(standard_only);
///
///         while (vnp)
///         {
///          matrix_info = vnp->ptr;
///          out_sz = snprintf(ptr,(size_t)buffer_sz, "%s \n", matrix_info->name);
/// 	 buffer_sz -= out_sz ; if ( buffer_sz < 0 ) buffer_sz = 0;
///       ptr += strlen(ptr);
///       vnp = vnp->next;
///         }
/// ```
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

/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3792-3839
/// ```c
/// #define BUF_SZ_2048 (2048)
/// ...
///    out_sz = snprintf(ptr, (size_t)buffer_sz,"Gap existence and extension values of %ld and %ld not supported for %s\nsupported values are:\n",
///       (long) gap_open, (long) gap_extend, matrix_name);
///
///    buffer_sz -= ( 1 + out_sz ); if ( buffer_sz < 0 ) buffer_sz = 0; // decrease out buffer size first, "+1" for null termination , once.
///    ptr += strlen(ptr);  // advance writing buffer position
///
///    vnp = head = BlastLoadMatrixValues(FALSE);
/// ...
///       for (index=0; index<max_number_values; index++)
///       {
///          if (BLAST_Nint(values[index][2]) == INT2_MAX)
///             out_sz = snprintf(ptr, (size_t)buffer_sz, "%ld, %ld\n", (long) BLAST_Nint(values[index][0]), (long) BLAST_Nint(values[index][1]));
///          else
///             out_sz = snprintf(ptr, (size_t)buffer_sz, "%ld, %ld, %ld\n", (long) BLAST_Nint(values[index][0]), (long) BLAST_Nint(values[index][1]), (long) BLAST_Nint(values[index][2]));
/// 	 buffer_sz -= out_sz ; if ( buffer_sz < 0 ) buffer_sz = 0; // decrease out buffer size first, "+1" for null termination , once.
///          ptr += strlen(ptr);
///       }
/// ```
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
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3374-3396
/// ```c
/// Int2 BLAST_GetProteinGapExistenceExtendParams(const char* matrixName,
///                                        Int4* gap_existence,
///                                        Int4* gap_extension)
/// {
///      Int4* gapOpen_arr,* gapExtend_arr,* pref_flags;
///      Int4 i; /*loop index*/
///      Int2 num_values = Blast_GetMatrixValues(matrixName, &gapOpen_arr,
///        &gapExtend_arr, NULL, NULL, NULL,  NULL, NULL, &pref_flags);
///
///      if (num_values <= 0)
///          return -1;
///
///      for(i = 1; i < num_values; i++) {
///          if(pref_flags[i]==BLAST_MATRIX_BEST) {
///              (*gap_existence) = gapOpen_arr[i];
///              (*gap_extension) = gapExtend_arr[i];
///              break;
///          }
///      }
/// ```
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

/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1174-1209
/// ```c
/// Int2 BLAST_GetSuggestedThreshold(EBlastProgramType program_number, const char* matrixName, double* threshold)
/// {
///
///     const double kB62_threshold = 11;
/// ...
///     if(strcasecmp(matrixName, "BLOSUM62") == 0)
///         *threshold = kB62_threshold;
///     else if(strcasecmp(matrixName, "BLOSUM45") == 0)
///         *threshold = 14;
///     else if(strcasecmp(matrixName, "BLOSUM62_20") == 0)
///         *threshold = 100;
///     else if(strcasecmp(matrixName, "BLOSUM80") == 0)
///         *threshold = 12;
///     else if(strcasecmp(matrixName, "PAM30") == 0)
///         *threshold = 16;
///     else if(strcasecmp(matrixName, "PAM70") == 0)
///         *threshold = 14;
///     else if(strcasecmp(matrixName, "IDENTITY") == 0)
///         *threshold = 27;
///     else
///         *threshold = kB62_threshold;
///
///     if (Blast_SubjectIsTranslated(program_number) == TRUE)
///         *threshold += 2;  /* Covers tblastn, tblastx, psi-tblastn rpstblastn. */
///     else if (Blast_QueryIsTranslated(program_number) == TRUE)
///         *threshold += 1;
/// ```
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

/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1211-1236
/// ```c
///     if(strcasecmp(matrixName, "BLOSUM62") == 0)
///         *window_size = kB62_windowsize;
///     else if(strcasecmp(matrixName, "BLOSUM45") == 0)
///         *window_size = 60;
///     else if(strcasecmp(matrixName, "BLOSUM80") == 0)
///         *window_size = 25;
///     else if(strcasecmp(matrixName, "PAM30") == 0)
///         *window_size = 15;
///     else if(strcasecmp(matrixName, "PAM70") == 0)
///         *window_size = 20;
///     else
///         *window_size = kB62_windowsize;
/// ```
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

/// The options that NCBI's `BLAST_ValidateOptions` checks for blastp and tblastn, as their
/// option handlers set them.
pub struct ProteinOptionsCheck<'a> {
    /// Not `-ungapped`.
    pub gapped: bool,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// The matrix name as typed.
    pub matrix_name: &'a str,
    pub threshold: f64,
    pub word_size: i32,
    /// Whether the compressed-alphabet lookup table is used (`eCompressedAaLookupTable`).
    pub compressed_lookup: bool,
    pub evalue: f64,
}

/// NCBI's `BLAST_ValidateOptions` for blastp and tblastn, in its order, with NCBI's
/// messages (`CInputException`): the scoring options, the lookup table options, the hit
/// saving options and the word size of the identity matrix. LOSAT has no greedy or
/// Smith-Waterman extension option, `-xdrop_ungap` or cutoff score option for these
/// programs, so NCBI's other checks cannot fail.
///
/// NCBI reference: c++/src/algo/blast/core/blast_options.c:1759-1776
/// ```c
///    if ((status = BlastExtensionOptionsValidate(program_number, ext_options,
///                                                blast_msg)) != 0)
///        return status;
///    if ((status = BlastScoringOptionsValidate(program_number, score_options,
///                                                blast_msg)) != 0)
///        return status;
///    if ((status = LookupTableOptionsValidate(program_number,
///                     lookup_options, blast_msg)) != 0)
///        return status;
///    if ((status = BlastInitialWordOptionsValidate(program_number,
///                     word_options, blast_msg)) != 0)
///        return status;
///    if ((status = BlastHitSavingOptionsValidate(program_number, hit_options,
///                                                blast_msg)) != 0)
///        return status;
/// ```
pub fn validate_protein_options(check: &ProteinOptionsCheck<'_>) -> anyhow::Result<()> {
    use crate::blastinput::app::options_error;
    // NCBI reference: c++/src/algo/blast/core/blast_options.c:910-936 (the scoring
    // check of BLAST_ValidateOptions; blastp and tblastn may use IDENTITY)
    // ```c
    //                 if (options->gapped_calculation && !Blast_ProgramIsRpsBlast(program_number))
    //                 {
    //                     Int2 status=0;
    //                     Boolean std_matrix_only =
    //                         (program_number != eBlastTypeBlastp &&
    //                          program_number != eBlastTypeTblastn);
    //                     if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
    //                           options->gap_extend, options->matrix, std_matrix_only)) != 0)
    //                      {
    // 			if (status == 1)
    // 			{
    // ...
    // 				buffer = BLAST_PrintMatrixMessage(options->matrix,
    //                                                   std_matrix_only);
    // ...
    // 			else if (status == 2)
    // 			{
    // ...
    // 				buffer = BLAST_PrintAllowedValues(options->matrix,
    //                         options->gap_open, options->gap_extend);
    // ```
    if check.gapped {
        match karlin_blk_gapped_load_from_tables(
            check.gap_open,
            check.gap_extend,
            check.matrix_name,
            false,
        ) {
            Ok(_) => {}
            Err(TableStatus::UnknownMatrix) => {
                return Err(options_error(&print_matrix_message(
                    check.matrix_name,
                    false,
                )))
            }
            Err(TableStatus::UnknownGapCosts) => {
                return Err(options_error(&print_allowed_values(
                    check.matrix_name,
                    check.gap_open,
                    check.gap_extend,
                )))
            }
        }
    }
    // NCBI reference: c++/src/algo/blast/core/blast_options.c:1303-1395 (the lookup
    // check of BLAST_ValidateOptions for blastp and tblastn)
    // ```c
    //     if (program_number != eBlastTypeBlastn &&
    //         program_number != eBlastTypeMapping &&
    //         (!Blast_ProgramIsRpsBlast(program_number)) &&
    //         options->threshold <= 0)
    //     {
    //         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "Non-zero threshold required");
    // ...
    //             if (options->word_size > 7) {
    //                 Blast_MessageWrite(blast_msg, eBlastSevError,
    //                                    kBlastMessageNoContext,
    //                                    "Word-size must be less than "
    //                                    "8 for a tblastn, blastp or blastx search");
    // ...
    //         if (options->word_size > 5 &&
    //             options->lut_type != eCompressedAaLookupTable) {
    //            Blast_MessageWrite(blast_msg, eBlastSevError,
    //                               kBlastMessageNoContext,
    //                               "Blastp, Blastx or Tblastn with word size"
    //                               " > 5 requires a "
    //                               "compressed alphabet lookup table");
    //            return BLASTERR_OPTION_VALUE_INVALID;
    //         }
    //         else if (options->lut_type == eCompressedAaLookupTable &&
    //                  options->word_size != 5 && options->word_size != 6 &&
    //                  options->word_size != 7) {
    //            Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "Compressed alphabet lookup table requires "
    //                          "word size 5, 6 or 7");
    // ```
    // (`threshold <= 0` is false for NaN, as in C.)
    if check.threshold <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    if check.word_size <= 0 {
        return Err(options_error("Word-size must be greater than zero"));
    }
    if check.word_size > 7 {
        return Err(options_error(
            "Word-size must be less than 8 for a tblastn, blastp or blastx search",
        ));
    }
    if check.word_size > 5 && !check.compressed_lookup {
        return Err(options_error(
            "Blastp, Blastx or Tblastn with word size > 5 requires a compressed alphabet lookup table",
        ));
    }
    if check.compressed_lookup && !(5..=7).contains(&check.word_size) {
        return Err(options_error(
            "Compressed alphabet lookup table requires word size 5, 6 or 7",
        ));
    }
    // NCBI reference: c++/src/algo/blast/core/blast_options.c:1518-1521 (the hit
    // saving check; LOSAT's blastp and tblastn have no cutoff score option)
    // ```c
    // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //          "expect value or cutoff score must be greater than zero");
    // ```
    if check.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    // NCBI reference: c++/src/algo/blast/core/blast_options.c:1783-1795
    // ```c
    //        char* matrix = BLAST_StrToUpper(score_options->matrix);
    //        Boolean is_identity = strcmp(matrix, "IDENTITY") == 0;
    // ...
    //        if (lookup_options->word_size > 5 && is_identity) {
    //
    //            Blast_MessageWrite(blast_msg, eBlastSevError,
    //                               kBlastMessageNoContext,
    //                               "Word size larger than 5 is not supported for "
    //                               "the identity scoring matrix");
    // ```
    if check.word_size > 5 && check.matrix_name.eq_ignore_ascii_case("IDENTITY") {
        return Err(options_error(
            "Word size larger than 5 is not supported for the identity scoring matrix",
        ));
    }
    Ok(())
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
