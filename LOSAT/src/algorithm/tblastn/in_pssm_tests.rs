use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1011-1015: the code 32 is not NCBI's.
#[test]
fn the_genetic_code_32_is_refused_with_a_pssm() {
    assert!(check_db_gencode(1).is_ok());
    assert!(check_db_gencode(33).is_ok());
    assert!(check_db_gencode(32)
        .unwrap_err()
        .to_string()
        .contains("-db_gencode 32 with -in_pssm"));
}
