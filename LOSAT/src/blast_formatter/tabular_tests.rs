//! blast_formatter's tabular driver: the formatter's matrix, the options' matrix name and the
//! ids that the rows show for a local id.

use super::*;

#[test]
fn the_formatter_matrix_is_ncbis_standard_matrix_of_the_name() {
    assert_eq!(formatter_matrix(None), None);
    assert_eq!(formatter_matrix(Some(b"")), None);
    assert_eq!(
        formatter_matrix(Some(b"BLOSUM62")),
        Some(ScoringMatrix::Blosum62)
    );
    assert_eq!(
        formatter_matrix(Some(b"blosum45")),
        Some(ScoringMatrix::Blosum45)
    );
    assert_eq!(formatter_matrix(Some(b"PAM30")), Some(ScoringMatrix::Pam30));
    assert_eq!(
        formatter_matrix(Some(b"pam250")),
        Some(ScoringMatrix::Pam250)
    );
    assert_eq!(formatter_matrix(Some(b"BLOSUM62x")), None);
    assert_eq!(formatter_matrix(Some(b"PAM40")), None);
    assert_eq!(
        formatter_matrix(Some(b"IDENTITYxyz")),
        Some(ScoringMatrix::Identity)
    );
    assert_eq!(formatter_matrix(Some(b"other")), None);
}

#[test]
fn a_local_id_shows_the_first_word_of_the_title() {
    let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
    let word = |ids: Vec<ParsedSeqId>| match &ids[0] {
        ParsedSeqId::Local(ObjectId::Str(text)) => String::from_utf8(text.clone()).unwrap(),
        other => panic!("{other:?}"),
    };
    assert_eq!(word(replaced_ids([&local], b"mq1 a query", false)), "mq1");
    assert_eq!(
        word(replaced_ids([&local], b"mq1 a query", true)),
        "Query_1"
    );
    assert_eq!(word(replaced_ids([&local], b"", false)), "Query_1");
    assert_eq!(word(replaced_ids([&local], b" leading", false)), "Query_1");
    let numbered = ParsedSeqId::Local(ObjectId::Id(7));
    assert_eq!(word(replaced_ids([&numbered], b"", false)), "7");
}

#[test]
fn only_the_covered_programs_are_written() {
    for program in [
        EProgram::Blastn,
        EProgram::Megablast,
        EProgram::Blastp,
        EProgram::PsiBlast,
        EProgram::DeltaBlast,
        EProgram::PhiBlastp,
        EProgram::RpsBlast,
        EProgram::Blastx,
        EProgram::Tblastn,
        EProgram::Tblastx,
        EProgram::RpsTblastn,
        EProgram::PsiTblastn,
    ] {
        assert!(covered(program), "{program:?}");
    }
    for program in [EProgram::PhiBlastn] {
        assert!(!covered(program), "{program:?}");
    }
}
