use super::*;

#[test]
fn protein_tables_follow_ncbi() {
    assert_eq!(
        protein_fasta_letters(&[0, 1, 21, 24, 25, 26, 27]),
        b"-AXU*OJ"
    );
    assert_eq!(
        protein_string_letters(&[0, 1, 21, 24, 25, 26, 27]),
        [255, b'A', b'X', b'U', 255, b'O', b'J']
    );
    for code in [28, 31, 32, 100, 201, 255] {
        assert_eq!(protein_fasta_letters(&[code]), [255], "{code}");
        assert_eq!(protein_string_letters(&[code]), [255], "{code}");
    }
}
