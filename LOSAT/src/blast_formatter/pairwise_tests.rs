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
fn epilog_values_are_the_handles_defaults() {
    let values = |program| {
        let v = EpilogValues::of(&options(program, Vec::new()));
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
fn epilog_values_take_the_archives_calls() {
    let v = EpilogValues::of(&options(
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
    let v = EpilogValues::of(&options(
        EProgram::Blastn,
        vec![Setter::MatchReward(1), Setter::MismatchPenalty(-3)],
    ));
    assert_eq!((v.reward, v.penalty, v.gapped), (1, -3, true));
}

// A stored gap residue is a `-` at a column where the row has letters; a gap column's `-` is
// not one.
#[test]
fn stored_gap_columns_are_letters_shown_as_gaps() {
    let ds = PairDenseg {
        starts: vec![0, 0, -1, 3, 3, 5],
        lens: vec![3, 2, 2],
        strands: [Strand::Unknown, Strand::Unknown],
        first_minus: false,
    };
    // Query row: letters, then a gap of 2, then letters; its `-` at column 1 is a residue.
    assert_eq!(stored_gap_columns(&ds, b"A-C--DE", 0), vec![1]);
    // Subject row: letters everywhere; its `-` at columns 3 and 6 are residues.
    assert_eq!(stored_gap_columns(&ds, b"ACD-EF-", 1), vec![3, 6]);
    assert!(stored_gap_columns(&ds, b"ACDEFGH", 1).is_empty());
}

// The subject letters that the rows read: from the first to the last letter of row 1.
#[test]
fn subject_range_spans_row_one() {
    let ds = PairDenseg {
        starts: vec![10, 20, 13, -1, 15, 23],
        lens: vec![3, 2, 4],
        strands: [Strand::Plus, Strand::Plus],
        first_minus: false,
    };
    assert_eq!(subject_range(&ds), Some((20, 26)));
}
