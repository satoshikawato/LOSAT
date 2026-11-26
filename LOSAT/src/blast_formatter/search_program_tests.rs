//! The Blast Type that blast_formatter's formatter writes for each program of an archive.

use super::*;

// NCBI BLAST+ 2.17.0 blast_formatter on archives of each program (line RG item F3, fixture rows
// of `LOSAT/tests/blast_formatter`; DELTA-BLAST and PHI-BLAST: rows `of8_delta`, `of8_phi`): the
// handle's type names the program, DELTA-BLAST is named on its own.
#[test]
fn each_program_has_ncbi_s_blast_type() {
    let cases = [
        (EProgram::Blastn, ("blastn", 1)),
        (EProgram::Megablast, ("blastn", 1)),
        (EProgram::DiscMegablast, ("blastn", 1)),
        (EProgram::VecScreen, ("blastn", 1)),
        (EProgram::Blastp, ("blastp", 2)),
        (EProgram::Blastx, ("blastx", 3)),
        (EProgram::Tblastn, ("tblastn", 4)),
        (EProgram::Tblastx, ("tblastx", 5)),
        (EProgram::RpsBlast, ("rpsblast", 6)),
        (EProgram::RpsTblastn, ("rpstblastn", 7)),
        (EProgram::PsiBlast, ("psiblast", 10)),
        (EProgram::PsiTblastn, ("psitblastn", 11)),
        (EProgram::DeltaBlast, ("deltablast", 14)),
    ];
    for (program, expected) in cases {
        assert_eq!(blast_type(program).unwrap(), expected, "{program:?}");
    }
    // DELTA-BLAST's handle is PSI-BLAST's type: protein subjects, no SAM.
    let delta = program_type(EProgram::DeltaBlast);
    assert!(!delta.subject_nucleotide && !delta.blastn);
    assert!(program_type(EProgram::Megablast).blastn);
}

// The PHI-BLAST names are no program of `ProgramNameToEnum`: NCBI's engine error, exit status 3;
// the mapper is refused.
#[test]
fn phi_blast_is_ncbi_s_engine_error_and_the_mapper_is_refused() {
    for (program, name) in [
        (EProgram::PhiBlastp, "phiblastp"),
        (EProgram::PhiBlastn, "phiblastn"),
    ] {
        let error = blast_type(program).unwrap_err();
        let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
        assert_eq!(native.exit, 3, "{name}");
        assert_eq!(
            String::from_utf8_lossy(&native.message),
            format!("LOSAT engine error: Program type '{name}' not supported\n")
        );
    }
    let error = blast_type(EProgram::Mapper).unwrap_err();
    assert_eq!(
        error.downcast_ref::<Refusal>(),
        Some(&Refusal::Program("mapr2g"))
    );
}
