//! The Bioseq of a query or subject record in a BLAST archive.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:597-601
//! ```c++
//!     // (CFastaReader::ParseDefLine ... the raw defline is kept in a user object)
//!     CRef<CUser_object> uo(new CUser_object);
//!     uo->SetType().SetStr("CFastaReader");
//!     uo->AddField("DefLine", NStr::PrintableString(defline));
//! ```
//! A record is `id { local str "Query_N" }` (`Subject_N` for `-subject`), `descr` with the
//! `CFastaReader` user object (the raw defline line) and the `title` (absent when it is
//! empty), and `inst` (`repr raw`, `mol na|aa`, the length and the data: `ncbi2na` when every
//! nucleotide is A, C, G or T, else `ncbi4na` for the whole sequence; protein `ncbieaa`; an
//! empty record is `repr virtual` with no data).

use crate::report::seqid::SeqId;
use crate::report::serial::encoder::Encoder;

/// A record to write as a Bioseq.
#[derive(Clone, Copy, Debug)]
pub struct BioseqInput<'a> {
    /// `Query_N` or `Subject_N`.
    pub id: &'a str,
    /// The defline line as read, with its `>` and without the line end.
    pub defline: &'a [u8],
    /// The title (the defline after the `>` up to the first control byte, trimmed).
    pub title: &'a [u8],
    /// The letters as the reader keeps them (lower case where the lower-case mask covers
    /// them; a nucleotide `U` is already `T`).
    pub letters: &'a [u8],
    pub protein: bool,
}

/// `NStr::PrintableString(s, 0)`: a byte outside 0x20..0x7e, a backslash or a quote is
/// written as an escape.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3732-3845
/// (`\t \r \n \a \b \f \v` as letters, `\\` `\'` `\"` backslashed, a `?` next to a `?` as
/// `\?`, other bytes as reduced octal: `\303\251`, `\1`).
pub fn printable_string(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len() + 8);
    for (index, &byte) in text.iter().enumerate() {
        let escape = match byte {
            7 => Some(b'a'),
            8 => Some(b'b'),
            12 => Some(b'f'),
            13 => Some(b'r'),
            9 => Some(b't'),
            11 => Some(b'v'),
            10 => {
                out.extend_from_slice(b"\\n\\\n");
                continue;
            }
            b'\\' | b'\'' | b'"' => Some(byte),
            b'?' if (index > 0 && text[index - 1] == b'?')
                || (index + 1 < text.len() && text[index + 1] == b'?') =>
            {
                Some(b'?')
            }
            _ => None,
        };
        if let Some(escape) = escape {
            out.push(b'\\');
            out.push(escape);
        } else if !(0x20..0x7f).contains(&byte) {
            let mut reduce = index + 1 == text.len() || !(b'0'..=b'7').contains(&text[index + 1]);
            out.push(b'\\');
            let (v1, v2, v3) = (byte >> 6, (byte >> 3) & 7, byte & 7);
            if v1 != 0 || !reduce {
                out.push(b'0' + v1);
                reduce = false;
            }
            if v2 != 0 || !reduce {
                out.push(b'0' + v2);
            }
            out.push(b'0' + v3);
        } else {
            out.push(byte);
        }
    }
    out
}

/// `NCBI4na` code of an upper-case nucleotide letter (`-ACMGRSVTWYHKDBN`); a letter outside
/// the table is `N`.
fn ncbi4na_code(letter: u8) -> u8 {
    match letter {
        b'-' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// The packed data of a nucleotide sequence: `ncbi2na` (`true`) when every letter is A, C,
/// G or T, else `ncbi4na`.
pub fn pack_nucleotides(letters: &[u8]) -> (bool, Vec<u8>) {
    let two_bit = |letter: u8| match letter.to_ascii_uppercase() {
        b'A' => Some(0u8),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    if letters.iter().all(|&letter| two_bit(letter).is_some()) {
        let mut packed = vec![0u8; letters.len().div_ceil(4)];
        for (index, &letter) in letters.iter().enumerate() {
            let code = two_bit(letter).unwrap_or(0);
            packed[index / 4] |= code << (6 - 2 * (index % 4));
        }
        (true, packed)
    } else {
        let mut packed = vec![0u8; letters.len().div_ceil(2)];
        for (index, &letter) in letters.iter().enumerate() {
            let code = ncbi4na_code(letter.to_ascii_uppercase());
            packed[index / 2] |= code << (4 - 4 * (index % 2));
        }
        (false, packed)
    }
}

/// Writes the Bioseq as the element it is in a `SEQUENCE OF` / `SET OF` (the caller has
/// called `element()`; a Seq-entry wrapper is the caller's too).
pub fn write_bioseq<E: Encoder>(enc: &mut E, bioseq: &BioseqInput<'_>) {
    enc.begin_sequence("Bioseq");
    enc.member("id", 0);
    enc.begin_list(true);
    enc.element();
    SeqId::Local(bioseq.id).write(enc);
    enc.end_list();
    enc.end_member();

    enc.member("descr", 1);
    enc.begin_named("Seq-descr");
    enc.begin_list(true);
    enc.element();
    enc.variant("Seqdesc", "user", 13);
    enc.begin_sequence("User-object");
    enc.member("type", 1);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(b"CFastaReader");
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    enc.begin_list(false);
    enc.element();
    enc.begin_sequence("User-field");
    enc.member("label", 0);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(b"DefLine");
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    enc.variant("User-field_data", "str", 0);
    enc.visible_str(&printable_string(bioseq.defline));
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
    if !bioseq.title.is_empty() {
        enc.element();
        enc.variant("Seqdesc", "title", 4);
        enc.visible_str(bioseq.title);
        enc.end_variant();
    }
    enc.end_list();
    enc.end_named();
    enc.end_member();

    enc.member("inst", 2);
    enc.begin_sequence("Seq-inst");
    enc.member("repr", 0);
    let empty = bioseq.letters.is_empty();
    if empty {
        enc.enumerated("virtual", 1);
    } else {
        enc.enumerated("raw", 2);
    }
    enc.end_member();
    enc.member("mol", 1);
    if bioseq.protein {
        enc.enumerated("aa", 3);
    } else {
        enc.enumerated("na", 4);
    }
    enc.end_member();
    enc.member("length", 2);
    enc.int(bioseq.letters.len() as i32);
    enc.end_member();
    if !empty {
        enc.member("seq-data", 6);
        if bioseq.protein {
            enc.variant("Seq-data", "ncbieaa", 7);
            let upper: Vec<u8> = bioseq
                .letters
                .iter()
                .map(|letter| letter.to_ascii_uppercase())
                .collect();
            enc.begin_named("NCBIeaa");
            enc.string_store(&upper);
            enc.end_named();
        } else {
            let (two_bit, packed) = pack_nucleotides(bioseq.letters);
            let name = if two_bit { "NCBI2na" } else { "NCBI4na" };
            if two_bit {
                enc.variant("Seq-data", "ncbi2na", 2);
            } else {
                enc.variant("Seq-data", "ncbi4na", 3);
            }
            enc.begin_named(name);
            enc.octets(&packed);
            enc.end_named();
        }
        enc.end_variant();
        enc.end_member();
    }
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn printable_string_follows_ncbi() {
        assert_eq!(printable_string(b">plain text"), b">plain text");
        assert_eq!(printable_string(b"a\tb\\c'd\"e"), b"a\\tb\\\\c\\'d\\\"e");
        assert_eq!(printable_string(b"what??"), b"what\\?\\?");
        assert_eq!(printable_string(b"one?two"), b"one?two");
        assert_eq!(printable_string(b"caf\xc3\xa9"), b"caf\\303\\251");
        assert_eq!(printable_string(b"\x01x"), b"\\1x");
        assert_eq!(printable_string(b"\x011"), b"\\0011");
        assert_eq!(printable_string(b"\x07\x08\x0c\r\x0b"), b"\\a\\b\\f\\r\\v");
    }

    #[test]
    fn nucleotides_pack_into_two_or_four_bits() {
        assert_eq!(pack_nucleotides(b"ACGT"), (true, vec![0b0001_1011]));
        assert_eq!(pack_nucleotides(b"acgta"), (true, vec![0b0001_1011, 0]));
        // N makes the whole sequence ncbi4na: A=1 C=2 N=15 T=8.
        assert_eq!(pack_nucleotides(b"ACNT"), (false, vec![0x12, 0xf8]));
        assert_eq!(pack_nucleotides(b"RYK"), (false, vec![0x5a, 0xc0]));
    }
}
