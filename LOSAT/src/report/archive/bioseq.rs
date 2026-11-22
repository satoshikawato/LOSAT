//! The Bioseq of a query or subject record in a BLAST archive.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:597-603
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

impl BioseqInput<'_> {
    /// Whether `encoding` writes the strings of the Bioseq (`write_bioseq`) and of its Seq-loc
    /// (`write_seq_loc`) as they are (`ArchiveEncoding::writes_text`): the title and the
    /// Seq-ids. (The `DefLine` field is `NStr::PrintableString` of the line: visible.) The masks
    /// and the Seq-aligns of the record use the same Seq-ids.
    pub fn written_as_is(&self, encoding: super::ArchiveEncoding) -> bool {
        use super::seq_id_written_as_is;
        encoding.writes_text(self.title)
            && match self.ids {
                Some(ids) => ids
                    .iter()
                    .all(|id| seq_id_written_as_is(encoding, SeqId::Db(id))),
                None => seq_id_written_as_is(encoding, SeqId::Local(self.id)),
            }
            && self
                .seq_loc
                .is_none_or(|loc| seq_id_written_as_is(encoding, SeqId::Db(loc.id)))
    }
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
    // `TSeqPos`: for an empty record, one less than `seqlen` is written unsigned (x_FastaToSeqLoc).
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
    // The default mode `fNewLine_Quote` writes a line feed as `\n` (ncbistr.hpp:2730).
    enc.visible_str(&crate::ncbi_str::printable_string(&bioseq.defline));
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
#[path = "bioseq_tests.rs"]
mod tests;
