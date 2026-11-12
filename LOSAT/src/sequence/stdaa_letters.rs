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
/// table's 28 codes; a code above 27 (not written by makeblastdb) reads as NUL as NCBI's
/// table has it.
const STDAA_TO_EAA: [u8; 32] = *b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ\0\0\0\0";

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
const STDAA_TO_IUPACAA: [u8; 32] = [
    255, b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q',
    b'R', b'S', b'T', b'V', b'W', b'X', b'Y', b'Z', b'U', 255, b'O', b'J', 255, 255, 255, 255,
];

/// The FASTA letters of a protein's NCBIstdaa residues.
pub(crate) fn protein_fasta_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| STDAA_TO_EAA[usize::from(code & 31)])
        .collect()
}

/// The `GetSequenceAsString` letters of a protein's NCBIstdaa residues.
pub(crate) fn protein_string_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| STDAA_TO_IUPACAA[usize::from(code & 31)])
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
    }
}
