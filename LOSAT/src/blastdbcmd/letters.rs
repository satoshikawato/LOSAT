//! The letters of the sequences blastdbcmd writes: the reverse complement of `-strand minus`,
//! and the CRC-32 of `-outfmt %h`. The letters of a protein's NCBIstdaa codes are
//! `sequence::stdaa_letters`.

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:2149-2178
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
#[path = "letters_tests.rs"]
mod tests;
