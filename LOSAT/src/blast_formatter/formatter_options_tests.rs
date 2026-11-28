use super::*;
use crate::report::archive::read::options::SeqRange;

fn options(program: EProgram, setters: Vec<Setter>) -> BuiltOptions {
    BuiltOptions {
        program,
        setters,
        entrez_query: None,
        first_db_seq: None,
        final_db_seq: None,
        gi_list: None,
        negative_gi_list: None,
        taxid_list: None,
        negative_taxid_list: None,
        db_filtering_algorithm_id: None,
        db_filtering_algorithm_key: None,
        subject_masking_type: None,
        query_masks: None,
        query_range: SeqRange::EMPTY,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:136-219
// The factory defaults of each program's handle, which the epilog prints where the archive
// sets nothing (set A: blastn 2 -3 5 2, megablast 1 -2 0 0, dc-megablast window 40, rpsblast
// without a word threshold, the protein handles' 11 1 11 40).
#[test]
fn the_values_are_the_handles_defaults() {
    let values = |program| {
        let v = FormatterOptions::of(&options(program, Vec::new()));
        (
            v.reward,
            v.penalty,
            v.gap_open,
            v.gap_extend,
            v.word_threshold,
            v.window_size,
        )
    };
    assert_eq!(values(EProgram::Blastn), (2, -3, 5, 2, 0, 0));
    assert_eq!(values(EProgram::Megablast), (1, -2, 0, 0, 0, 0));
    assert_eq!(values(EProgram::DiscMegablast), (2, -3, 5, 2, 0, 40));
    let (_, _, open, extend, threshold, window) = values(EProgram::RpsBlast);
    assert_eq!((open, extend, threshold, window), (11, 1, 0, 40));
    for program in [
        EProgram::Blastp,
        EProgram::PsiBlast,
        EProgram::DeltaBlast,
        EProgram::PhiBlastp,
    ] {
        let (_, _, open, extend, threshold, window) = values(program);
        assert_eq!(
            (open, extend, threshold, window),
            (11, 1, 11, 40),
            "{program:?}"
        );
    }
}

// The archive's calls on the factory options, in order, the last one of each kind kept
// (blastp-short's PAM30 9 1 16 15 of set A; an ungapped search; sum statistics).
#[test]
fn the_values_take_the_archives_calls() {
    let v = FormatterOptions::of(&options(
        EProgram::Blastp,
        vec![
            Setter::GapOpeningCost(10),
            Setter::GapOpeningCost(9),
            Setter::WordThreshold(16),
            Setter::WindowSize(15),
            Setter::GappedMode(false),
            Setter::SumStatisticsMode(true),
        ],
    ));
    assert_eq!(
        (
            v.gap_open,
            v.gap_extend,
            v.word_threshold,
            v.window_size,
            v.gapped,
            v.sum_statistics
        ),
        (9, 1, 16, 15, false, true)
    );
    let v = FormatterOptions::of(&options(
        EProgram::Blastn,
        vec![Setter::MatchReward(1), Setter::MismatchPenalty(-3)],
    ));
    assert_eq!((v.reward, v.penalty, v.gapped), (1, -3, true));
}

// The matrix name: the last `MatrixName` call, else BLOSUM62 for the protein handles and none
// (NULL) for the nucleotide handles and RPS-BLAST's.
#[test]
fn the_matrix_name_is_the_last_call_or_the_handles() {
    let name = |program, setters| FormatterOptions::of(&options(program, setters)).matrix_name;
    assert_eq!(
        name(EProgram::Blastp, Vec::new()),
        Some(b"BLOSUM62".to_vec())
    );
    assert_eq!(
        name(EProgram::DeltaBlast, Vec::new()),
        Some(b"BLOSUM62".to_vec())
    );
    assert_eq!(name(EProgram::Blastn, Vec::new()), None);
    assert_eq!(name(EProgram::RpsBlast, Vec::new()), None);
    assert_eq!(
        name(
            EProgram::Blastp,
            vec![
                Setter::MatrixName("PAM30".into()),
                Setter::MatrixName("BLOSUM80".into())
            ]
        ),
        Some(b"BLOSUM80".to_vec())
    );
}
