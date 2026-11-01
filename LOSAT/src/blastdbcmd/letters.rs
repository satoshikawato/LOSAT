//! The letters of the sequences blastdbcmd writes: NCBI's conversions of the stored
//! codes (NCBIstdaa, the 2-bit and 4-bit nucleotide codes) to the letters of its FASTA
//! and `-outfmt %s` paths, the reverse complement of `-strand minus`, and the CRC-32 of
//! `-outfmt %h`.

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
pub(super) fn protein_fasta_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| STDAA_TO_EAA[usize::from(code & 31)])
        .collect()
}

/// The `GetSequenceAsString` letters of a protein's NCBIstdaa residues.
pub(super) fn protein_string_letters(stdaa: &[u8]) -> Vec<u8> {
    stdaa
        .iter()
        .map(|&code| STDAA_TO_IUPACAA[usize::from(code & 31)])
        .collect()
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:2149-2178
// ```c++
// const Uint1 CIupacnaCmp::scm_Table[256] = {
//     ...
//     'T', 'V', 'G', 'H', 'D', 'E', 'F', 'C',    // A B C D E F G H
//     'I', 'J', 'M', 'L', 'K', 'N', 'O', 'P',    // I J K L M N O P
//     'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X',    // Q R S T U V W X
//     'R', 'Z', ...
//     't', 'v', 'g', 'h', 'd', 'e', 'f', 'c',    // a b c d e f g h
//     ...
// ```
/// The complement of one IUPAC letter (`CIupacnaCmp`): `A`/`T`, `C`/`G`, `M`/`K`, `R`/`Y`,
/// `V`/`B`, `H`/`D`, `U` to `A`, the case kept, every other byte itself. The FASTA path
/// complements the stored 4-bit codes (the Seq-code table of ncbi4na), which gives the same
/// letters for the letters a database holds (`-ACMGRSVTWYHKDBN`: the gap and `N`, `S`, `W`
/// are their own complement).
fn complement(letter: u8) -> u8 {
    match letter {
        b'A' => b'T',
        b'B' => b'V',
        b'C' => b'G',
        b'D' => b'H',
        b'G' => b'C',
        b'H' => b'D',
        b'K' => b'M',
        b'M' => b'K',
        b'R' => b'Y',
        b'T' => b'A',
        b'U' => b'A',
        b'V' => b'B',
        b'Y' => b'R',
        b'a' => b't',
        b'b' => b'v',
        b'c' => b'g',
        b'd' => b'h',
        b'g' => b'c',
        b'h' => b'd',
        b'k' => b'm',
        b'm' => b'k',
        b'r' => b'y',
        b't' => b'a',
        b'u' => b'a',
        b'v' => b'b',
        b'y' => b'r',
        other => other,
    }
}

/// Reverses `letters` and complements each (`CSeqManip::ReverseComplement` of IUPAC letters;
/// the FASTA path's `fReverseStrand`).
pub(super) fn reverse_complement(letters: &mut [u8]) {
    letters.reverse();
    for letter in letters.iter_mut() {
        *letter = complement(*letter);
    }
}

/// The CRC-32 of zlib (`CChecksum::eCRC32ZIP`): the IEEE polynomial reflected
/// (`0xEDB88320`), the register starting at and finally inverted by `0xFFFFFFFF`.
fn crc32_zip(bytes: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:1037-1046
// ```c++
// Uint4 CBlastSeqUtil::GetSeqHash(const char* buffer, int length)
// {
//     CChecksum crc(CChecksum::eCRC32ZIP);
//
//     for(int ii = 0; ii < length; ii++) {
//         if (buffer[ii] != '\n')
//             crc.AddChars(buffer+ii,1);
//     }
//     return (crc.GetChecksum() ^ (0xFFFFFFFFL));
// }
// ```
/// `CBlastSeqUtil::GetSeqHash` of the letters of a whole sequence (the `GetSequenceAsString`
/// letters), the line feeds skipped.
pub(super) fn sequence_hash(letters: &[u8]) -> u32 {
    if letters.contains(&b'\n') {
        let kept: Vec<u8> = letters.iter().copied().filter(|&b| b != b'\n').collect();
        return !crc32_zip(&kept);
    }
    !crc32_zip(letters)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:265-272
// ```c++
// string CBlastDB_SeqFormatter::x_GetSeqHash(CSeqDB::TOID oid)
// {
// 	string seq;
// 	m_BlastDb.GetSequenceAsString(oid, seq);
// 	CNcbiOstrstream os;
// 	os << std::showbase << std::hex << std::uppercase << CBlastSeqUtil::GetSeqHash(seq.c_str(), seq.size());
// 	return CNcbiOstrstreamToString(os);
// }
// ```
/// The text of `-outfmt %h`: the hash in upper-case hexadecimal with the prefix `0X`
/// (`showbase`), which the stream leaves out for zero.
pub(super) fn hash_text(hash: u32) -> String {
    if hash == 0 {
        "0".to_string()
    } else {
        format!("0X{hash:X}")
    }
}

/// The ranges `masks` (the end exclusive) of a sequence in lower case, in `letters`, which are
/// the sequence's `range` (`from`, `to` inclusive) or all of it: a mask is clipped to them, and
/// one that ends before it starts is empty.
pub(super) fn apply_masks(letters: &mut [u8], masks: &[(u32, u32)], range: Option<(u32, u32)>) {
    let from = range.map_or(0, |(from, _)| from) as usize;
    for &(first, second) in masks {
        let start = (first as usize).max(from);
        let end = (second as usize).min(from + letters.len());
        if start < end {
            letters[start - from..end - from].make_ascii_lowercase();
        }
    }
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

    #[test]
    fn reverse_complement_keeps_case_and_ambiguity() {
        let mut letters = b"ACGTMRWSYKVHDBN-acgtu".to_vec();
        reverse_complement(&mut letters);
        assert_eq!(letters, b"aacgt-NVHDBMRSWYKACGT");
    }

    #[test]
    fn crc32_is_zlib_s() {
        // The CRC-32 of "123456789" is the standard check value 0xCBF43926.
        assert_eq!(crc32_zip(b"123456789"), 0xCBF4_3926);
        assert_eq!(sequence_hash(b"1234\n56789"), !0xCBF4_3926);
        assert_eq!(hash_text(0), "0");
        assert_eq!(hash_text(0x9E5E_1FA), "0X9E5E1FA");
    }

    #[test]
    fn masks_are_lower_case_inside_the_letters() {
        let masked = |length: usize, masks: &[(u32, u32)], range| {
            let mut letters = vec![b'A'; length];
            apply_masks(&mut letters, masks, range);
            String::from_utf8(letters).unwrap()
        };
        assert_eq!(masked(10, &[(2, 4), (8, 10)], None), "AAaaAAAAaa");
        assert_eq!(masked(6, &[(2, 4), (8, 10)], Some((3, 8))), "aAAAAa");
        assert_eq!(masked(5, &[(2, 4)], Some((4, 8))), "AAAAA");
        assert_eq!(masked(4, &[(2, 100), (3, 2), (1, 1)], None), "AAaa");
    }
}
