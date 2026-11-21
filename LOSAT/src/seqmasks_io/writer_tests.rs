use super::*;

#[test]
fn fasta_letters_follow_ncbis_one_step_mask_iterator() {
    // A mask list where one position passes two intervals: NCBI advances the iterator
    // once per position, so the second interval [3, 3] after [0, 1] is reached at 2 and
    // its letter 3 is masked; a third interval right after is reached only at its own
    // positions.
    let mut out = Vec::new();
    fasta_letters(b"ACGTACGT", &[(0, 1), (3, 3), (4, 5)], false, &mut out);
    assert_eq!(out, b"acGtacGT\n");
    let mut out = Vec::new();
    fasta_letters(&[b'A'; 61], &[(59, 60)], true, &mut out);
    assert_eq!(&out[..], [&[b'A'; 59][..], b"N\nN\n"].concat());
}
