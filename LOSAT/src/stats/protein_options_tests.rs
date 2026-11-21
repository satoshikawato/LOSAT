use super::*;

// NCBI BLAST+ 2.17.0 `blastp -matrix FOO` (BlastScoringOptionsValidate status 1).
#[test]
fn matrix_message_lists_ncbis_matrices() {
    assert_eq!(
        print_matrix_message("FOO", false),
        b"FOO is not a supported matrix, supported matrices are:\nBLOSUM80 \nBLOSUM62 \nBLOSUM50 \nBLOSUM45 \nPAM250 \nBLOSUM90 \nPAM30 \nPAM70 \nIDENTITY \n"
    );
    assert!(!print_matrix_message("FOO", true)
        .windows(8)
        .any(|w| w == b"IDENTITY"));
    // The 1024-byte buffer: a long name truncates the message.
    let long = "M".repeat(2000);
    assert_eq!(print_matrix_message(&long, false).len(), 1023);
    // ... inside a character, whose first byte stays (probe R9 of line LIB-B).
    let cut = format!("{}{}", "a".repeat(1000), "\u{e9}".repeat(30));
    let message = print_matrix_message(&cut, false);
    assert_eq!(message.len(), 1023);
    assert_eq!(&message[1020..], b"\xc3\xa9\xc3");
}

// NCBI BLAST+ 2.17.0 `blastp -matrix PAM30 -gapopen 1 -gapextend 1` (status 2).
#[test]
fn allowed_values_list_the_table_in_order() {
    assert_eq!(
        print_allowed_values("PAM30", 1, 1),
        b"Gap existence and extension values of 1 and 1 not supported for PAM30\nsupported values are:\n32767, 32767\n7, 2\n6, 2\n5, 2\n10, 1\n9, 1\n8, 1\n15, 3\n14, 2\n14, 1\n13, 3\n"
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
