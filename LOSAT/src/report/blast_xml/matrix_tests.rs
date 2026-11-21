use super::*;

/// The stop of a run's XML formats (BLAST XML, the XML2 and JSON documents) with a matrix name
/// that NCBI's formatting table does not spell: every other format, and every name it spells,
/// formats.
#[test]
fn the_xml_formats_fail_for_a_name_ncbi_does_not_spell() {
    // The facts of the search: a search of subject records, an iterative program.
    let searches = [(false, false), (true, false), (false, true), (true, true)];
    for outfmt in ["5", "13", "14", "15", "16", " 5"] {
        assert!(is_structured_xml(outfmt), "{outfmt}");
        for (subject_search, iterative) in searches {
            let failure = xml_formatting_failure([outfmt], "blosum62", subject_search, iterative);
            let held = (subject_search || iterative) && outfmt.trim() != "5";
            assert_eq!(
                failure,
                Some(if held {
                    XmlMatrixFailure::SecondQueryOrEnd
                } else {
                    XmlMatrixFailure::FirstQuery
                }),
                "{outfmt} {subject_search} {iterative}"
            );
            assert_eq!(
                xml_formatting_failure([outfmt], "BLOSUM62", subject_search, iterative),
                None,
                "{outfmt}"
            );
        }
    }
    for outfmt in ["0", "6", "7 qseqid", "8", "11", "12", "17", "abc"] {
        assert!(!is_structured_xml(outfmt), "{outfmt}");
        assert_eq!(
            xml_formatting_failure([outfmt], "blosum62", true, true),
            None,
            "{outfmt}"
        );
    }
    assert_eq!(
        xml_formatting_failure(["6", "16"], "blosum62", false, false),
        Some(XmlMatrixFailure::FirstQuery)
    );
    // Several formats: the first to fail stops the run.
    assert_eq!(
        xml_formatting_failure(["16", "15"], "blosum62", true, false),
        Some(XmlMatrixFailure::SecondQueryOrEnd)
    );
    assert_eq!(
        xml_formatting_failure(["16", "5"], "blosum62", true, false),
        Some(XmlMatrixFailure::FirstQuery)
    );
    let error = unsupported_score_matrix();
    let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
    assert_eq!(
        (native.exit, native.message.as_slice()),
        (
            3,
            b"LOSAT engine error: unsupported score matrix\n".as_slice()
        )
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
