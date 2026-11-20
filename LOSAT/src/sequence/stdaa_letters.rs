//! The letters of a protein's NCBIstdaa codes, as NCBI's `CSeqVector` converts them: ncbieaa
//! (blastdbcmd's FASTA path) and IUPACaa (`eCoding_Iupac`: blastdbcmd's `-outfmt %s` and `%h`,
//! segmasker's writers of a database's sequences).

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1822-1850
// ```c++
// const Uint1 CStdaaToEaa::scm_Table[32] = {
//     '-', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
//     'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
//     'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z',
//     'U', '*', 'O', 'J', 0, 0, 0, 0
// };
// ```
/// NCBIstdaa to ncbieaa, the letters of the FASTA path (`CSeqVector` with `Ncbieaa`): the
/// table's 28 codes.
const STDAA_TO_EAA: [u8; 28] = *b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1785-1813
// ```c++
// const Uint1 CStdaaToIupacaa::scm_Table[32] = {
//     255, 'A', 'B', 'C', 'D', 'E', 'F', 'G',
//     'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
//     'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z',
//     'U', 255, 'O', 'J', 255, 255, 255, 255
// };
// ```
/// NCBIstdaa to iupacaa, the letters of `GetSequenceAsString` (`-outfmt %s`, `%h`): the gap
/// (0) and the stop (25) are the byte 255, as NCBI's table has them.
const STDAA_TO_IUPACAA: [u8; 28] = [
    255, b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q',
    b'R', b'S', b'T', b'V', b'W', b'X', b'Y', b'Z', b'U', 255, b'O', b'J',
];

/// The letter of `code` in `table`; a code that is no NCBIstdaa residue (28 and above, which
/// makeblastdb does not write) is the byte 255 in both conversions, whatever its low bits (NCBI
/// BLAST+ 2.17.0 blastdbcmd `%f` and `%s`, segmasker `-outfmt fasta`, with `.psq` bytes 28, 31,
/// 32, 100, 201 and 255; review of line RG item K2, M1).
fn letter(table: &[u8; 28], code: u8) -> u8 {
    table.get(usize::from(code)).copied().unwrap_or(255)
}

/// The FASTA letters of a protein's NCBIstdaa residues.
pub(crate) fn protein_fasta_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| letter(&STDAA_TO_EAA, code))
        .collect()
}

/// The `GetSequenceAsString` letters of a protein's NCBIstdaa residues.
pub(crate) fn protein_string_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| letter(&STDAA_TO_IUPACAA, code))
        .collect()
}

#[cfg(test)]
mod tests {
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
}
