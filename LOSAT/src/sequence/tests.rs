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

/// The lines written a part at a time are those of the letters written whole, whatever the
/// parts.
#[test]
fn lines_written_in_parts_are_the_lines_of_the_whole() {
    let letters: Vec<u8> = (0..300).map(|i| b"ACGT"[i % 4]).collect();
    for count in [0, 1, 69, 70, 71, 139, 140, 141, 300] {
        let mut whole = Vec::new();
        write_fasta_lines(&letters[..count], 70, &mut whole);
        for part in [1, 3, 69, 70, 71, 200] {
            let mut out = Vec::new();
            let mut lines = FastaLines::new(70);
            for chunk in letters[..count].chunks(part) {
                lines.write(chunk, &mut out);
            }
            lines.end(&mut out);
            assert_eq!(out, whole, "{count} letters in parts of {part}");
        }
    }
}
