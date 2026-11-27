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

// Review of PS15b-1, N3: PSI-TBLASTN has one name (`ProgramType::name`, `Blast_ProgramNameFromType`),
// which its reports upper-case: the program line of outfmt 0 and 7 and BLAST XML's version.
#[test]
fn psi_tblastn_has_one_name() {
    let psi = program_type(true).expect("PSI-TBLASTN's");
    assert_eq!(psi.name(), "psitblastn");
    assert_eq!(PROGRAM, psi.name().to_ascii_uppercase());
    assert!(super::super::tabular::program_version(true).starts_with(&format!("{PROGRAM} ")));
    assert_eq!(program_type(false), None);
    assert!(!crate::report::record::carries_query_masks(Some(psi)));
    assert!(crate::report::record::carries_query_masks(None));
}
