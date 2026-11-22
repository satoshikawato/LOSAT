use super::*;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2131-2155 and c++/src/objtools/align_format/align_format_util.cpp:3563-3581 (cited above and at
// `protein_display_score`): a masked capital compares as its capital; a PSSM query's
// lower-case letter is identical to no subject letter and scores -1000.
#[test]
fn a_masked_letter_compares_as_its_capital_and_a_bioseq_letter_as_it_is() {
    let matrix = ScoringMatrix::Blosum62;
    assert_eq!(
        display_positives(b"aKv", b"AKI", matrix, QueryCase::Masked),
        3
    );
    assert_eq!(
        display_positives(b"aKv", b"AKI", matrix, QueryCase::Bioseq),
        1
    );
    assert_eq!(
        display_positives(b"-A", b"-a", matrix, QueryCase::Bioseq),
        2
    );
    let pssm = QueryLetters::new(b"AKV", Some(b"aKv"));
    assert_eq!(
        (pssm.get(0), pssm.get(1), pssm.get(3)),
        (Some(b'a'), Some(b'K'), None)
    );
    assert_eq!(pssm.case(), QueryCase::Bioseq);
    let record = QueryLetters::new(b"aKv", None);
    assert_eq!((record.get(0), record.get(2)), (Some(b'A'), Some(b'V')));
    assert_eq!(record.case(), QueryCase::Masked);
    assert_eq!(
        QueryCase::Masked.compared_letters()[usize::from(b'k')],
        b'K'
    );
    assert_eq!(
        QueryCase::Bioseq.compared_letters()[usize::from(b'k')],
        b'k'
    );
    assert_eq!(
        QueryCase::Masked.compared_letters()[usize::from(b'-')],
        b'-'
    );
    assert!(pssm.shows_other_letters());
    assert!(QueryLetters::new(b"A-K", None).shows_other_letters());
    assert!(!record.shows_other_letters());
}

// A gap character is a residue of its own side's sequence: the query reads one at a
// substitution and an insertion, the subject at a substitution and a deletion.
#[test]
fn the_gap_residues_lie_at_their_sides_columns() {
    let ops = [
        GapEditOp::Sub(2),
        GapEditOp::Ins(1),
        GapEditOp::Del(1),
        GapEditOp::Sub(1),
    ];
    // Query `XA--K` read from 1: its `-` at positions 2 and 3 lie at columns 1 and 2 (the
    // insertion).
    assert_eq!(
        stored_gap_columns(b"XA--K", b'-', 1, &ops, RowSide::Query),
        [1, 2]
    );
    // Subject residues 0, 9, 0, 5 at columns 0, 1, 3, 4: code 0 at 0 and 3.
    assert_eq!(
        stored_gap_columns(&[0, 9, 0, 5], 0, 0, &ops, RowSide::Subject),
        [0, 3]
    );
    assert!(stored_gap_columns(b"MKVLA", b'-', 0, &ops, RowSide::Query).is_empty());
}
