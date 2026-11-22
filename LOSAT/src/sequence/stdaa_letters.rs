//! The letters of a protein's NCBIstdaa codes, as NCBI's `CSeqVector` converts them: ncbieaa
//! (blastdbcmd's FASTA path) and IUPACaa (`eCoding_Iupac`: blastdbcmd's `-outfmt %s` and `%h`,
//! segmasker's writers of a database's sequences).

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1822-1850
/// NCBIstdaa to ncbieaa, the letters of the FASTA path (`CSeqVector` with `Ncbieaa`): the
/// table's 28 codes ([`all_codes`]).
const STDAA_TO_EAA: [u8; 256] = all_codes(b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ");

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1785-1813
/// NCBIstdaa to iupacaa, the letters of `GetSequenceAsString` (`-outfmt %s`, `%h`): the gap
/// (0) and the stop (25) are the byte 255, as NCBI's table has them ([`all_codes`]).
const STDAA_TO_IUPACAA: [u8; 256] = all_codes(&[
    255, b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q',
    b'R', b'S', b'T', b'V', b'W', b'X', b'Y', b'Z', b'U', 255, b'O', b'J',
]);

/// The letters of the 28 NCBIstdaa codes `residues` for every byte: a code that is no NCBIstdaa
/// residue (28 and above, which makeblastdb does not write) is the byte 255 in both conversions,
/// whatever its low bits (NCBI BLAST+ 2.17.0 blastdbcmd `%f` and `%s`, segmasker `-outfmt fasta`,
/// with `.psq` bytes 28, 31, 32, 100, 201 and 255; review of line RG item K2, M1). A table of
/// every byte makes the letter one load, without a branch per residue.
const fn all_codes(residues: &[u8; 28]) -> [u8; 256] {
    let mut table = [255u8; 256];
    let mut code = 0;
    while code < residues.len() {
        table[code] = residues[code];
        code += 1;
    }
    table
}

/// The letter of `code` in `table`.
fn letter(table: &[u8; 256], code: u8) -> u8 {
    table[usize::from(code)]
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
#[path = "stdaa_letters_tests.rs"]
mod tests;
