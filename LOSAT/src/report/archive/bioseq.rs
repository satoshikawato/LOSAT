//! The Bioseq of a query or subject record in a BLAST archive.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:597-603
//! ```c++
//!     if ( !TestFlag(fNoUserObjs) ) {
//!         // store the raw defline in a User-object for reference
//!         CRef<CSeqdesc> desc(new CSeqdesc);
//!         desc->SetUser().SetType().SetStr("CFastaReader");
//!         desc->SetUser().AddField("DefLine", NStr::PrintableString(s));
//!         m_CurrentSeq->SetDescr().Set().push_back(desc);
//!     }
//! ```
//! A record is `id { local str "Query_N" }` (`Subject_N` for `-subject`), `descr` with the
//! `CFastaReader` user object (the raw defline line) and the `title` (absent when it is
//! empty), and `inst` (`repr raw`, `mol na|aa`, the length and the data: `ncbi2na` when every
//! nucleotide is A, C, G or T, else `ncbi4na` for the whole sequence; protein `ncbieaa`; an
//! empty record is `repr virtual` with no data).

use crate::report::seqid::SeqId;
use crate::report::serial::encoder::Encoder;

/// A record to write as a Bioseq.
#[derive(Clone, Debug)]
pub struct BioseqInput<'a> {
    /// `Query_N` or `Subject_N`.
    pub id: &'a str,
    /// With `-parse_deflines`: the record's Seq-ids, written instead of `id`.
    pub ids: Option<&'a [crate::seqid::ParsedSeqId]>,
    /// With `-parse_deflines` and an id that is not local: the query as a Seq-loc.
    pub seq_loc: Option<QuerySeqLoc<'a>>,
    /// The defline line as read, with its `>` and without the line end (`>` and the title for a
    /// record that kept none).
    pub defline: std::borrow::Cow<'a, [u8]>,
    /// The title (the defline after the `>` up to the first control byte, trimmed).
    pub title: &'a [u8],
    /// The letters as the reader keeps them (lower case where the lower-case mask covers
    /// them; a nucleotide `U` is already `T`).
    pub letters: &'a [u8],
    pub protein: bool,
}

/// The strand of a query's Seq-loc (`CBlastFastaInputSource::x_FastaToSeqLoc`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqLocStrand {
    Unknown,
    Plus,
    Minus,
    Both,
}

impl SeqLocStrand {
    /// The strand a nucleotide program reads with `-strand value`.
    pub fn of_option(value: &str) -> Self {
        match value {
            "plus" => Self::Plus,
            "minus" => Self::Minus,
            _ => Self::Both,
        }
    }
}

/// A query as the Seq-loc of the archive's `seq-loc-list`: the searched interval and the
/// record's Seq-align id (`CBlastInputReader` with `-parse_deflines`, an id that is not local).
#[derive(Clone, Copy, Debug)]
pub struct QuerySeqLoc<'a> {
    pub from: i32,
    pub to: i32,
    pub strand: SeqLocStrand,
    pub id: &'a crate::seqid::ParsedSeqId,
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

/// Writes the query as the element it is in the `seq-loc-list` (the caller has called
/// `element()`).
pub fn write_seq_loc<E: Encoder>(enc: &mut E, loc: &QuerySeqLoc<'_>) {
    enc.variant("Seq-loc", "int", 3);
    enc.begin_sequence("Seq-interval");
    enc.member("from", 0);
    enc.int(loc.from);
    enc.end_member();
    enc.member("to", 1);
    // `TSeqPos`: `seqlen - 1` of an empty record is written unsigned (x_FastaToSeqLoc).
    if loc.to < 0 {
        enc.big_int(i64::from(u32::MAX));
    } else {
        enc.int(loc.to);
    }
    enc.end_member();
    enc.member("strand", 2);
    match loc.strand {
        SeqLocStrand::Unknown => enc.enumerated("unknown", 0),
        SeqLocStrand::Plus => enc.enumerated("plus", 1),
        SeqLocStrand::Minus => enc.enumerated("minus", 2),
        SeqLocStrand::Both => enc.enumerated("both", 3),
    }
    enc.end_member();
    enc.member("id", 3);
    SeqId::Db(loc.id).write(enc);
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
}

/// Writes the Bioseq as the element it is in a `SEQUENCE OF` / `SET OF` (the caller has
/// called `element()`; a Seq-entry wrapper is the caller's too).
pub fn write_bioseq<E: Encoder>(enc: &mut E, bioseq: &BioseqInput<'_>) {
    enc.begin_sequence("Bioseq");
    enc.member("id", 0);
    enc.begin_list(true);
    match bioseq.ids {
        Some(ids) => {
            for id in ids {
                enc.element();
                SeqId::Db(id).write(enc);
            }
        }
        None => {
            enc.element();
            SeqId::Local(bioseq.id).write(enc);
        }
    }
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
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:601
    // ```c++
    //         desc->SetUser().AddField("DefLine", NStr::PrintableString(s));
    // ```
    // The default mode `fNewLine_Quote` writes a line feed as `\n` (ncbistr.hpp:2730).
    enc.visible_str(crate::ncbi_str::printable_string(&bioseq.defline).as_bytes());
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

    /// The `DefLine` field is `NStr::PrintableString(s)` of the raw line (fasta.cpp:601) in the
    /// default `fNewLine_Quote` mode: a line feed, which only a record the library builds can
    /// hold, is written as `\n`.
    #[test]
    fn the_defline_field_quotes_a_line_feed() {
        let mut buf = Vec::new();
        let mut enc = crate::report::serial::asn_text::AsnTextEncoder::new(&mut buf);
        let bioseq = BioseqInput {
            id: "Query_1",
            ids: None,
            seq_loc: None,
            defline: std::borrow::Cow::Borrowed(b">a\nb"),
            title: b"a",
            letters: b"MKV",
            protein: true,
        };
        write_bioseq(&mut enc, &bioseq);
        assert_eq!(enc.error(), None);
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains(r#"str ">a\nb""#), "{text}");
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
