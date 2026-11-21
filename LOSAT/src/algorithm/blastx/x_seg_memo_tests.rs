use super::*;

#[test]
fn memo_returns_reference_result_for_equal_bytes_and_recomputes_otherwise() {
    use std::cell::Cell;
    let calls = Cell::new(0usize);
    let reference = |s: &[u8]| {
        calls.set(calls.get() + 1);
        vec![MaskedInterval::new(0, s.len())]
    };
    let a = [1u8, 2, 3, 4];
    let b = [1u8, 2, 3, 5];
    let first = x_subject_seg(&a, reference);
    let again = x_subject_seg(&a, reference);
    assert_eq!(first, again);
    assert_eq!(calls.get(), 1);
    let other = x_subject_seg(&b, reference);
    assert_eq!(other, vec![MaskedInterval::new(0, 4)]);
    assert_eq!(calls.get(), 2);
    // A different length with the same prefix is a different subject.
    let shorter = x_subject_seg(&a[..3], reference);
    assert_eq!(shorter, vec![MaskedInterval::new(0, 3)]);
    assert_eq!(calls.get(), 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
#[test]
fn memo_of_the_redo_masker_equals_the_masker() {
    use crate::utils::matrix::aa_char_to_ncbistdaa;
    use crate::utils::seg::{SegMasker, SegParams};
    let masker =
        |data: &[u8]| SegMasker::with_params(&SegParams::new(10, 1.8, 2.1)).mask_sequence(data);
    // Test protein: made from the first residues of E. coli AsnA (RefSeq NP_418200.1) and other letters.
    let subjects: Vec<Vec<u8>> = [
        "MKTAYIAKQRQISFVKSHFSRQAAAAAAAAAAAAAAAAAPPPPPPPGGGGGSLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ",
        "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYG",
    ]
    .iter()
    .map(|s| s.bytes().map(aa_char_to_ncbistdaa).collect())
    .collect();
    assert!(
        !masker(&subjects[0]).is_empty(),
        "the fixture must be masked"
    );
    for subject in subjects
        .iter()
        .chain(subjects.iter())
        .chain(subjects.iter().rev())
    {
        assert_eq!(x_subject_seg(subject, masker), masker(subject));
        assert_eq!(x_subject_seg(subject, masker), masker(subject));
    }
}
