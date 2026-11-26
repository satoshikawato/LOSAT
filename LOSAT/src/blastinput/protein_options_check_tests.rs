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
            matrix_name: matrix_name.as_bytes(),
            threshold: 11.0,
            word_size,
            compressed_lookup,
            hitlist_size: 500,
            evalue: 10.0,
            psiblast,
            psi_tblastn: false,
            x_dropoff: None,
            culling_limit: None,
            best_hit_overhang: None,
            best_hit_score_edge: None,
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
