use super::*;

#[test]
fn sequence_lines_are_the_width_with_a_shorter_last() {
    let mut out = Vec::new();
    write_fasta_lines(b"ACGTACGT", 3, &mut out);
    assert_eq!(out, b"ACG\nTAC\nGT\n");
    out.clear();
    write_fasta_lines(b"ACGTAC", 3, &mut out);
    assert_eq!(out, b"ACG\nTAC\n");
}
