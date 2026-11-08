//! NCBI's `BLAST_ValidateOptions` for blastp, tblastn and psiblast (`validate_protein_options`):
//! the option checks over the protein matrix tables of `stats::protein_options`.

use crate::stats::protein_options::{
    karlin_blk_gapped_load_from_tables, print_allowed_values, print_matrix_message, TableStatus,
};

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
    /// The hit list size (`hit_options->hitlist_size`).
    pub hitlist_size: usize,
    pub evalue: f64,
    /// The program is psiblast (`eBlastTypePsiBlast`) rather than blastp or tblastn.
    pub psiblast: bool,
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
    // psiblast is neither blastp nor tblastn: its matrices are the standard ones
    // (`std_matrix_only`).
    let std_matrix_only = check.psiblast;
    if check.gapped {
        match karlin_blk_gapped_load_from_tables(
            check.gap_open,
            check.gap_extend,
            check.matrix_name,
            std_matrix_only,
        ) {
            Ok(_) => {}
            Err(TableStatus::UnknownMatrix) => {
                return Err(options_error(&print_matrix_message(
                    check.matrix_name,
                    std_matrix_only,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1350-1357
    // ```c
    //         else if (program_number == eBlastTypePsiBlast &&
    //                  options->word_size > 4) {
    //                 Blast_MessageWrite(blast_msg, eBlastSevError,
    //                                kBlastMessageNoContext,
    //                                "Word-size must be less "
    //                                "than 5 for psiblast");
    //                 return BLASTERR_OPTION_VALUE_INVALID;
    //         }
    // ```
    // The checks of the compressed lookup table below are blastp's, tblastn's and blastx's
    // only (blast_options.c:1376-1378).
    if check.psiblast {
        if check.word_size > 4 {
            return Err(options_error("Word-size must be less than 5 for psiblast"));
        }
    } else if check.word_size > 7 {
        return Err(options_error(
            "Word-size must be less than 8 for a tblastn, blastp or blastx search",
        ));
    }
    if !check.psiblast && check.word_size > 5 && !check.compressed_lookup {
        return Err(options_error(
            "Blastp, Blastx or Tblastn with word size > 5 requires a compressed alphabet lookup table",
        ));
    }
    if !check.psiblast && check.compressed_lookup && !(5..=7).contains(&check.word_size) {
        return Err(options_error(
            "Compressed alphabet lookup table requires word size 5, 6 or 7",
        ));
    }
    // NCBI reference: c++/src/algo/blast/core/blast_options.c:1511-1517 (the hit saving
    // check of the hit list size, before the e-value's)
    // ```c
    // 	if (options->hitlist_size < 1)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "No hits are being saved");
    // 		return BLASTERR_OPTION_VALUE_INVALID;
    // 	}
    // ```
    if check.hitlist_size < 1 {
        return Err(options_error("No hits are being saved"));
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

    // blast_options.c:1336-1357 and 913-915: psiblast's word size is below 5 (NCBI 2.17.0
    // `psiblast -word_size 5`), its matrices are the standard ones, and the compressed lookup
    // table's checks are blastp's, tblastn's and blastx's.
    #[test]
    fn psiblast_takes_word_sizes_below_5_and_the_standard_matrices() {
        let check = |word_size: i32, compressed_lookup: bool, matrix_name: &str, psiblast: bool| {
            validate_protein_options(&ProteinOptionsCheck {
                gapped: true,
                gap_open: 11,
                gap_extend: 1,
                matrix_name,
                threshold: 11.0,
                word_size,
                compressed_lookup,
                hitlist_size: 500,
                evalue: 10.0,
                psiblast,
            })
            .map_err(|error| error.to_string())
        };
        assert!(check(4, false, "BLOSUM62", true).is_ok());
        for (word_size, compressed) in [(5, true), (6, true), (8, true), (5, false)] {
            let error = check(word_size, compressed, "BLOSUM62", true).unwrap_err();
            assert!(
                error.contains("Word-size must be less than 5 for psiblast"),
                "{word_size}: {error}"
            );
        }
        assert!(check(5, true, "BLOSUM62", false).is_ok());
        assert!(check(8, true, "BLOSUM62", false)
            .unwrap_err()
            .contains("Word-size must be less than 8 for a tblastn, blastp or blastx search"));
        let identity = check(3, false, "IDENTITY", true).unwrap_err();
        assert!(
            identity.contains("IDENTITY is not a supported matrix"),
            "{identity}"
        );
        assert!(!identity.contains("\nIDENTITY \n"), "{identity}");
    }
}
