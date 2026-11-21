use super::*;

#[test]
fn test_ncbistdaa_encoding() {
    assert_eq!(aa_char_to_ncbistdaa(b'A'), 1);
    assert_eq!(aa_char_to_ncbistdaa(b'*'), 25);
    assert_eq!(aa_char_to_ncbistdaa(b'X'), 21);
    assert_eq!(aa_char_to_ncbistdaa(b'J'), 27);
}

#[test]
fn test_standard_matrix_conversion() {
    assert_eq!(ncbistdaa_to_blosum62(1), 0);
    assert_eq!(ncbistdaa_to_blosum62(25), 24);
    assert_eq!(ncbistdaa_to_blosum62(21), 23);
    assert_eq!(ncbistdaa_to_blosum62(ncbistdaa::U), blosum62_order::C);
    assert_eq!(ncbistdaa_to_blosum62(ncbistdaa::O), blosum62_order::X);
}

#[test]
fn test_blosum62_scores() {
    assert_eq!(blosum62_score(ncbistdaa::A, ncbistdaa::A), 4);
    assert_eq!(blosum62_score(ncbistdaa::STOP, ncbistdaa::STOP), 1);
    assert_eq!(blosum62_score(ncbistdaa::X, ncbistdaa::X), -1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1592
#[test]
fn test_blosum62_ncbistdaa_direct_matches_generic_for_blastaa_alphabet() {
    for q in 0..BLASTAA_SIZE as u8 {
        for s in 0..BLASTAA_SIZE as u8 {
            assert_eq!(
                blosum62_score_ncbistdaa_direct(q, s),
                protein_score(ScoringMatrix::Blosum62, q, s),
                "direct BLOSUM62 score mismatch for NCBISTDAA {q}/{s}"
            );
        }
    }
}

// No NCBI counterpart: the rows give the scores of the functions they are built from for
// every two bytes (the bytes past the alphabet read one row).
#[test]
fn the_score_rows_hold_the_scores_of_every_two_bytes() {
    for a in 0..=u8::MAX {
        for b in 0..=u8::MAX {
            assert_eq!(
                BLOSUM62_NCBISTDAA_ROWS.score(a, b),
                blosum62_score_ncbistdaa_direct(a, b),
                "BLOSUM62 {a} {b}"
            );
            for matrix in ScoringMatrix::ALL {
                assert_eq!(
                    protein_score_rows(matrix).score(a, b),
                    protein_score(matrix, a, b),
                    "{matrix:?} {a} {b}"
                );
            }
        }
    }
    assert!(std::ptr::eq(
        protein_score_rows(ScoringMatrix::Pam30),
        protein_score_rows(ScoringMatrix::Pam30)
    ));
}

#[test]
fn test_blosum62_ncbistdaa_score_row_matches_direct_scores() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:562-563
    for q in 0..BLASTAA_SIZE as u8 {
        let row = blosum62_ncbistdaa_score_row(q as usize);
        for s in 0..BLASTAA_SIZE as u8 {
            assert_eq!(row[s as usize], blosum62_score_ncbistdaa_direct(q, s));
        }
    }
}

#[test]
fn test_non_blosum62_standard_matrices() {
    assert_eq!(
        protein_score(ScoringMatrix::Blosum45, ncbistdaa::A, ncbistdaa::A),
        5
    );
    assert_eq!(
        protein_score(ScoringMatrix::Blosum80, ncbistdaa::W, ncbistdaa::W),
        11
    );
    assert_eq!(
        protein_score(ScoringMatrix::Pam30, ncbistdaa::A, ncbistdaa::A),
        6
    );
    assert_eq!(
        protein_score(ScoringMatrix::Pam250, ncbistdaa::A, ncbistdaa::A),
        2
    );
}

#[test]
fn test_matrix_defscores_follow_ncbi() {
    assert_eq!(protein_defscore(ScoringMatrix::Blosum45), -5);
    assert_eq!(protein_defscore(ScoringMatrix::Blosum62), -4);
    assert_eq!(protein_defscore(ScoringMatrix::Blosum80), -6);
    assert_eq!(protein_defscore(ScoringMatrix::Pam30), -17);
    assert_eq!(protein_defscore(ScoringMatrix::Pam250), -8);
    assert_eq!(protein_defscore(ScoringMatrix::Identity), -5);
}

// NCBI reference (598d8ae6): c++/src/util/tables/sm_identity.c:38-96 (9 on the diagonal, -5
// elsewhere, the X row -5 throughout; defscore -5)
#[test]
fn identity_matrix_is_ncbis_packed_table() {
    let identity = standard_protein_matrix(ScoringMatrix::Identity);
    let symbols = b"ARNDCQEGHILKMFPSTWYVBJZX*";
    for (r, &row) in symbols.iter().enumerate() {
        for c in 0..symbols.len() {
            let expected = if r == c && row != b'X' { 9 } else { -5 };
            assert_eq!(
                identity.scores[r * PROTEIN_MATRIX_SIZE + c],
                expected,
                "{} {}",
                row as char,
                symbols[c] as char
            );
        }
    }
    assert_eq!(
        protein_score(ScoringMatrix::Identity, ncbistdaa::W, ncbistdaa::W),
        9
    );
    assert_eq!(
        protein_score(ScoringMatrix::Identity, ncbistdaa::X, ncbistdaa::X),
        -5
    );
    assert_eq!(
        protein_score(ScoringMatrix::Identity, ncbistdaa::A, ncbistdaa::GAP),
        BLAST_SCORE_MIN
    );
}

#[test]
fn test_sentinel_uses_blast_score_min() {
    assert_eq!(
        protein_score(ScoringMatrix::Blosum45, ncbistdaa::GAP, ncbistdaa::A),
        BLAST_SCORE_MIN
    );
    assert_eq!(
        protein_score(ScoringMatrix::Pam30, ncbistdaa::A, ncbistdaa::GAP),
        BLAST_SCORE_MIN
    );
    assert_eq!(
        protein_score(ScoringMatrix::Blosum62, ncbistdaa::GAP, ncbistdaa::GAP),
        BLAST_SCORE_MIN
    );
}

#[test]
fn test_u_and_o_follow_ncbi_matrix_aliases() {
    assert_eq!(
        protein_score(ScoringMatrix::Blosum62, ncbistdaa::U, ncbistdaa::A),
        protein_score(ScoringMatrix::Blosum62, ncbistdaa::C, ncbistdaa::A)
    );
    assert_eq!(
        protein_score(ScoringMatrix::Blosum62, ncbistdaa::O, ncbistdaa::W),
        protein_score(ScoringMatrix::Blosum62, ncbistdaa::X, ncbistdaa::W)
    );
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
#[test]
fn display_matrix_matches_pinned_cpp_all_ascii_matrices() {
    let expected = include_str!("../../tests/unit/blastx_stage_e_ascii_matrix_expected.tsv");
    let mut rows = 0;
    for line in expected.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        let matrix = match fields[0] {
            "BLOSUM45" => ScoringMatrix::Blosum45,
            "BLOSUM50" => ScoringMatrix::Blosum50,
            "BLOSUM62" => ScoringMatrix::Blosum62,
            "BLOSUM80" => ScoringMatrix::Blosum80,
            "BLOSUM90" => ScoringMatrix::Blosum90,
            "PAM30" => ScoringMatrix::Pam30,
            "PAM70" => ScoringMatrix::Pam70,
            "PAM250" => ScoringMatrix::Pam250,
            _ => panic!("unregistered independent matrix"),
        };
        let query: u8 = fields[1].parse().unwrap();
        let subject: u8 = fields[2].parse().unwrap();
        let score: i32 = fields[3].parse().unwrap();
        assert_eq!(
            protein_display_score(matrix, query, subject),
            score,
            "{line}"
        );
        rows += 1;
    }
    assert_eq!(rows, 131072);
}

// The slice search that the byte table in protein_display_score replaced.
fn protein_display_score_slice_search(matrix: ScoringMatrix, query: u8, subject: u8) -> i32 {
    const SYMBOLS: &[u8] = b"ARNDCQEGHILKMFPSTWYVBZX";
    if (query == b'U' && matches!(subject, b'U' | b'C')) || (query == b'C' && subject == b'U') {
        return protein_score(
            matrix,
            aa_char_to_ncbistdaa(b'C'),
            aa_char_to_ncbistdaa(b'C'),
        );
    }
    if query == b'*' && subject == b'*' {
        return 1;
    }
    if (query == b'*' && SYMBOLS.contains(&subject))
        || (subject == b'*' && SYMBOLS.contains(&query))
    {
        return -4;
    }
    if !SYMBOLS.contains(&query) || !SYMBOLS.contains(&subject) {
        return -1000;
    }
    protein_score(
        matrix,
        aa_char_to_ncbistdaa(query),
        aa_char_to_ncbistdaa(subject),
    )
}

#[test]
fn display_score_byte_table_equals_slice_search_for_all_byte_pairs() {
    for matrix in [
        ScoringMatrix::Blosum45,
        ScoringMatrix::Blosum50,
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum80,
        ScoringMatrix::Blosum90,
        ScoringMatrix::Pam30,
        ScoringMatrix::Pam70,
        ScoringMatrix::Pam250,
        ScoringMatrix::Identity,
    ] {
        for query in 0..=u8::MAX {
            for subject in 0..=u8::MAX {
                assert_eq!(
                    protein_display_score(matrix, query, subject),
                    protein_display_score_slice_search(matrix, query, subject),
                    "{matrix:?} {query} {subject}"
                );
            }
        }
    }
}
