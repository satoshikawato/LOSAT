//! BER decoding of the header record of a BLAST database (`.phr` / `.nhr`): the reader of
//! the encoding (`Reader`, shared with the ordinal-id decoder of `seqdb::header`) and the
//! decoder of a whole `Blast-def-line-set` with every member and every Seq-id alternative
//! (`decode_full_set`).
//!
//! NCBI writes the record with `MSerial_AsnBinary`: EXPLICIT context tags, indefinite length
//! (`80 ... 00 00`) on every constructed value, definite length on INTEGER and VisibleString,
//! minimal two's-complement integers. NCBI's reader (`CObjectIStreamAsnBinary`) reads both
//! length forms; so does this one.
//!
//! ASN.1 definitions (NCBI 598d8ae6): `c++/src/objects/blastdb/blastdb.asn:31-47`
//! (`Blast-def-line-set`, `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-74`
//! (`Seq-id` and its member types), `c++/src/objects/general/general.asn:54-80` (`Date`,
//! `Dbtag`, `Object-id`), `c++/src/objects/biblio/biblio.asn:119-124` (`Id-pat`).
//!
//! NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:32-47
//! ```text
//! Blast-def-line-set ::= SEQUENCE OF Blast-def-line  -- all deflines for an entry
//!
//! Blast-def-line ::= SEQUENCE {
//!     title VisibleString OPTIONAL,             -- simple title
//!     seqid SEQUENCE OF Seq-id,                 -- Regular NCBI Seq-Id
//!     -- taxonomy id, should not be set if multiple taxids are available
//!     -- (see links below)
//!     taxid  INTEGER OPTIONAL,
//!     memberships SEQUENCE OF INTEGER OPTIONAL, -- bit arrays
//!     -- Repurposed to store the (multiple) taxIDs associated with WP proteins.
//!     links SEQUENCE OF INTEGER OPTIONAL,
//!     -- In proteins this stores the PIG, in nucleotides this stores the "origin
//!     -- GIs" (if one sequence is described in the ASN.1 as an identical copy of
//!     -- another)
//!     other-info SEQUENCE OF INTEGER OPTIONAL
//! }
//! ```

use super::defline_set::{BlastDefline, DeflineSet};
use super::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};

// NCBI reference (598d8ae6): c++/include/serial/impl/objstrasnb.hpp:71-130
// ```c++
//         eContextSpecific    = 2 << 6,
//         ...
//         eConstructed        = 1 << 5,
//         ...
//         eInteger            = 2,
//         ...
//         eSequence           = 16,
//         ...
//         eVisibleString      = 26,
//         ...
//         eLongTag            = 31,
//         ...
//         eContainterTagByte      = TByte(eConstructed) | TByte(eSequence),
//         eIndefiniteLengthByte   = 0x80,
//         eEndOfContentsByte      = 0,
// ```
const TAG_INTEGER: u8 = 2;
const TAG_VISIBLE_STRING: u8 = 26;
pub(crate) const TAG_SEQUENCE: u8 = (1 << 5) | 16;
pub(crate) const TAG_CONTEXT_CONSTRUCTED: u8 = (2 << 6) | (1 << 5);
pub(crate) const TAG_NUMBER_MASK: u8 = 31;
const INDEFINITE_LENGTH_BYTE: u8 = 0x80;
/// Deepest nesting accepted inside a skipped member.
const MAX_SKIP_DEPTH: u32 = 8;

/// Member tags of `Blast-def-line` (blastdb.asn:34-46, in order).
pub(crate) const DEFLINE_TITLE: u8 = 0;
pub(crate) const DEFLINE_SEQID: u8 = 1;
pub(crate) const DEFLINE_TAXID: u8 = 2;
const DEFLINE_MEMBERSHIPS: u8 = 3;
const DEFLINE_LINKS: u8 = 4;
const DEFLINE_OTHER_INFO: u8 = 5;
pub(crate) const DEFLINE_LAST_MEMBER: u8 = DEFLINE_OTHER_INFO;
/// The names of the members of `Blast-def-line`, by tag.
const DEFLINE_MEMBERS: [&str; 6] = [
    "title",
    "seqid",
    "taxid",
    "memberships",
    "links",
    "other-info",
];

/// An error of the BER decoding: the record is not a valid encoding of the type. `pos` is the
/// byte of the record that the reader had reached.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BerError {
    pub pos: usize,
    pub message: String,
}

impl std::fmt::Display for BerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "header record, byte {}: {}", self.pos, self.message)
    }
}

impl std::error::Error for BerError {}

/// Where a constructed value ends.
#[derive(Clone, Copy)]
pub(crate) enum End {
    /// Indefinite length: two zero bytes close the value.
    Marker,
    /// Definite length: the value ends at this position.
    At(usize),
}

/// A reader of one BER-encoded record.
pub(crate) struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Reader { data, pos: 0 }
    }

    /// Whether every byte of the record was read.
    pub(crate) fn at_record_end(&self) -> bool {
        self.pos == self.data.len()
    }

    pub(crate) fn error(&self, message: impl Into<String>) -> BerError {
        BerError {
            pos: self.pos,
            message: message.into(),
        }
    }

    fn byte(&mut self) -> Result<u8, BerError> {
        let Some(&byte) = self.data.get(self.pos) else {
            return Err(self.error("the record ends inside a value"));
        };
        self.pos += 1;
        Ok(byte)
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], BerError> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.data.len());
        let Some(end) = end else {
            return Err(self.error("the record ends inside a value"));
        };
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    /// Reads one identifier byte. Tag numbers above 30 (the long form) do not
    /// occur in these types.
    pub(crate) fn tag(&mut self) -> Result<u8, BerError> {
        let tag = self.byte()?;
        if tag & TAG_NUMBER_MASK == TAG_NUMBER_MASK {
            return Err(self.error("long-form tag"));
        }
        Ok(tag)
    }

    pub(crate) fn expect_tag(&mut self, expected: u8, what: &str) -> Result<(), BerError> {
        let tag = self.tag()?;
        if tag != expected {
            return Err(self.error(format!(
                "tag 0x{tag:02x} where {what} (0x{expected:02x}) is expected"
            )));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:394-411
    // ```c++
    // size_t CObjectIStreamAsnBinary::ReadLengthLong(TByte byte)
    // {
    //     size_t lengthLength = byte - 0x80;
    //     if ( lengthLength == 0 ) {
    //         ThrowError(fFormatError, "unexpected indefinite length");
    //     }
    //     if ( lengthLength > sizeof(size_t) ) {
    //         ThrowError(fOverflow, "length overflow");
    //     }
    //     size_t length = m_Input.GetChar() & 0xff;
    //     if ( length == 0 ) {
    //         ThrowError(fFormatError, "illegal length start");
    //     }
    //     while ( --lengthLength > 0 ) {
    //         length = (length << 8) | (m_Input.GetChar() & 0xff);
    //     }
    //     return StartTagData(length);
    // }
    // ```
    /// Reads a definite length (short or long form).
    fn definite_length(&mut self, first: u8) -> Result<usize, BerError> {
        if first < INDEFINITE_LENGTH_BYTE {
            return Ok(usize::from(first));
        }
        let mut length_length = usize::from(first - INDEFINITE_LENGTH_BYTE);
        if length_length == 0 {
            return Err(self.error("unexpected indefinite length"));
        }
        if length_length > std::mem::size_of::<usize>() {
            return Err(self.error("length overflow"));
        }
        let mut length = usize::from(self.byte()?);
        if length == 0 {
            return Err(self.error("illegal length start"));
        }
        length_length -= 1;
        while length_length > 0 {
            length = (length << 8) | usize::from(self.byte()?);
            length_length -= 1;
        }
        Ok(length)
    }

    /// Reads the length of a constructed value. NCBI writes these with
    /// indefinite length and reads both forms (objistrasnb.inl:277-306).
    pub(crate) fn open(&mut self) -> Result<End, BerError> {
        let first = self.byte()?;
        if first == INDEFINITE_LENGTH_BYTE {
            return Ok(End::Marker);
        }
        let length = self.definite_length(first)?;
        let end = self
            .pos
            .checked_add(length)
            .filter(|&end| end <= self.data.len());
        match end {
            Some(end) => Ok(End::At(end)),
            None => Err(self.error("a value is longer than the record")),
        }
    }

    /// Whether the constructed value opened with `end` has no further
    /// element; consumes the end marker of an indefinite length.
    pub(crate) fn at_end(&mut self, end: End) -> Result<bool, BerError> {
        match end {
            End::At(limit) => {
                if self.pos > limit {
                    return Err(self.error("a value passes the end of its container"));
                }
                Ok(self.pos == limit)
            }
            End::Marker => {
                if self.pos >= self.data.len() {
                    return Err(self.error("the record ends inside a value"));
                }
                if self.data.get(self.pos) == Some(&0) {
                    if self.data.get(self.pos + 1) != Some(&0) {
                        return Err(self.error("end-of-contents marker cut short"));
                    }
                    self.pos += 2;
                    return Ok(true);
                }
                Ok(false)
            }
        }
    }

    pub(crate) fn close(&mut self, end: End, what: &str) -> Result<(), BerError> {
        if self.at_end(end)? {
            Ok(())
        } else {
            Err(self.error(format!("unexpected data after {what}")))
        }
    }

    /// Skips one value of any type. `depth` bounds the nesting: the skipped
    /// members of `Blast-def-line` are an INTEGER or a SEQUENCE OF INTEGER.
    pub(crate) fn skip_value(&mut self, depth: u32) -> Result<(), BerError> {
        let tag = self.tag()?;
        if tag & (1 << 5) != 0 {
            if depth >= MAX_SKIP_DEPTH {
                return Err(self.error("values nested too deeply"));
            }
            let end = self.open()?;
            while !self.at_end(end)? {
                self.skip_value(depth + 1)?;
            }
        } else {
            let first = self.byte()?;
            let length = self.definite_length(first)?;
            self.take(length)?;
        }
        Ok(())
    }

    // NCBI's reader also takes a UTF8String (tag 0x0c) where a VisibleString is expected,
    // with a warning (objistrasnb.cpp:276-312, `SERIAL_READ_ANY_VISIBLESTRING_TAG` 1).
    // makeblastdb writes only VisibleStrings, so both decoders refuse another tag.
    pub(crate) fn visible_string(&mut self) -> Result<&'a [u8], BerError> {
        self.expect_tag(TAG_VISIBLE_STRING, "a VisibleString")?;
        let first = self.byte()?;
        let length = self.definite_length(first)?;
        self.take(length)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:476-512
    // ```c++
    // void ReadStdSigned(CObjectIStreamAsnBinary& in, T& data)
    // {
    //     size_t length = in.ReadShortLength();
    //     if ( length == 0 ) {
    //         in.ThrowError(in.fFormatError, "zero length of number");
    //     }
    //     T n;
    //     if ( length > sizeof(data) ) {
    //         // skip
    //         --length;
    //         Int1 c = in.ReadSByte();
    //         if ( c != 0 && c != -1 ) {
    //             in.ThrowError(in.fOverflow, "overflow error");
    //         }
    //         ...
    //     }
    //     else {
    //         --length;
    //         n = in.ReadSByte();
    //     }
    //     while ( length > 0 ) {
    //         --length;
    //         n = (n << 8) | in.ReadByte();
    //     }
    //     data = n;
    //     in.EndOfTag();
    // }
    // ```
    // An INTEGER wider than eight content bytes is rejected instead of being
    // reduced as above; each caller checks the range of its member.
    pub(crate) fn integer(&mut self) -> Result<i64, BerError> {
        self.expect_tag(TAG_INTEGER, "an INTEGER")?;
        let first = self.byte()?;
        if first >= INDEFINITE_LENGTH_BYTE {
            return Err(self.error("long length of number"));
        }
        let length = usize::from(first);
        if length == 0 {
            return Err(self.error("zero length of number"));
        }
        if length > std::mem::size_of::<i64>() {
            return Err(self.error("overflow error"));
        }
        let bytes = self.take(length)?;
        let mut n = i64::from(bytes[0] as i8);
        for &byte in &bytes[1..] {
            n = (n << 8) | i64::from(byte);
        }
        Ok(n)
    }

    /// An INTEGER member that NCBI reads into an `int`.
    fn int(&mut self, what: &str) -> Result<i32, BerError> {
        let value = self.integer()?;
        i32::try_from(value)
            .map_err(|_| self.error(format!("{what} {value} does not fit a 32-bit integer")))
    }

    pub(crate) fn open_context(&mut self, number: u8, what: &str) -> Result<End, BerError> {
        self.expect_tag(TAG_CONTEXT_CONSTRUCTED | number, what)?;
        self.open()
    }

    /// A VisibleString as NCBI's reader gives it (`fix_visible_chars`), and the number of
    /// bytes it replaced.
    fn fixed_string(&mut self) -> Result<(Vec<u8>, usize), BerError> {
        let stored = self.visible_string()?;
        Ok((fix_visible_chars(stored), replaced_bytes(stored)))
    }

    /// The members of a SEQUENCE whose identifier was read and whose length `end` gives:
    /// `each(reader, member)` reads the value of each member present (inside its
    /// EXPLICIT tag), in order. Member numbers must increase and name one of `names`, the
    /// members of `type_name` in order.
    fn members(
        &mut self,
        end: End,
        type_name: &str,
        names: &[&str],
        mut each: impl FnMut(&mut Self, u8) -> Result<(), BerError>,
    ) -> Result<(), BerError> {
        let mut next_member = 0;
        while !self.at_end(end)? {
            let tag = self.tag()?;
            let member = tag & TAG_NUMBER_MASK;
            if tag & !TAG_NUMBER_MASK != TAG_CONTEXT_CONSTRUCTED
                || member < next_member
                || usize::from(member) >= names.len()
            {
                return Err(self.error(format!("tag 0x{tag:02x} is not a member of {type_name}")));
            }
            next_member = member + 1;
            let member_end = self.open()?;
            each(self, member)?;
            if !self.at_end(member_end)? {
                let name = names[usize::from(member)];
                return Err(self.error(format!("unexpected data after {type_name}.{name}")));
            }
        }
        Ok(())
    }

    /// Opens a SEQUENCE value of `type_name`.
    fn open_sequence(&mut self, type_name: &str) -> Result<End, BerError> {
        self.expect_tag(TAG_SEQUENCE, type_name)?;
        self.open()
    }

    /// Reads the identifier and length of a CHOICE alternative (`[n]` EXPLICIT) of
    /// `type_name` with at most `last + 1` alternatives.
    fn open_alternative(&mut self, type_name: &str, last: u8) -> Result<(u8, End), BerError> {
        let tag = self.tag()?;
        let alternative = tag & TAG_NUMBER_MASK;
        if tag & !TAG_NUMBER_MASK != TAG_CONTEXT_CONSTRUCTED || alternative > last {
            return Err(self.error(format!(
                "tag 0x{tag:02x} is not an alternative of {type_name}"
            )));
        }
        Ok((alternative, self.open()?))
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:724-738
// ```c++
//         for(string::iterator i = str.begin(); i < to; ++i) {
//             if ( BadVisibleChar(*i) ) {
//                 ...
//                 {
//                     char c = ReplaceVisibleChar(*i, fix_method, this, str, subst);
//                     if (c != '\0') {
//                         *i = subst;
// ```
// A VisibleString is read with the method `eFNP_ReplaceAndWarn` unless the
// parameter SERIAL / WRONG_CHARS_READ says otherwise (objistr.cpp:302-303):
// every byte outside ' ' ..= '~' (`GoodVisibleChar`, objistr.inl:456-459)
// becomes the substitute '#', UTF-8 included
// (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0, objistr.hpp:1163).
// NCBI also writes a warning for each such string; this reader reports none.
// makeblastdb stores titles already repaired this way.
pub(crate) fn fix_visible_chars(value: &[u8]) -> Vec<u8> {
    value
        .iter()
        .map(|&c| if (b' '..=b'~').contains(&c) { c } else { b'#' })
        .collect()
}

/// The number of bytes of a stored VisibleString that `fix_visible_chars` replaces.
pub(crate) fn replaced_bytes(value: &[u8]) -> usize {
    value.iter().filter(|c| !(b' '..=b'~').contains(c)).count()
}

/// The `taxid` of a definition line as LOSAT keeps it: NCBI's `TTaxId` is an `int`; a value
/// outside it, which makeblastdb does not write, reads as `i32::MAX`.
pub(crate) fn taxid_value(value: i64) -> i32 {
    i32::try_from(value).unwrap_or(i32::MAX)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2208-2218
// ```c++
//     CTempString raw_data = x_GetHdrAsn1Binary(oid);
//
//     if (! raw_data.size()) {
//         return bdls;
//     }
//
//     // Now create an ASN.1 object from the memory chunk provided here.
//
//     CObjectIStreamAsnBinary inpstr(raw_data.data(), raw_data.size());
//
//     bdls.Reset(new objects::CBlast_def_line_set);
// ```
/// Decodes the `Blast-def-line-set` of one sequence with every member: the titles, the
/// Seq-ids of every alternative, the taxids, memberships, links and other-info. An empty
/// record gives no definition line, as in NCBI. The `general` ids `BL_ORD_ID` keep the
/// volume-relative tag that the record stores (`SeqDb::defline_set` adds the volume start).
pub fn decode_full_set(record: &[u8]) -> Result<DeflineSet, BerError> {
    let mut deflines = Vec::new();
    if record.is_empty() {
        return Ok(DeflineSet { deflines });
    }
    let mut ber = Reader::new(record);
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line-set")?;
    let set_end = ber.open()?;
    while !ber.at_end(set_end)? {
        deflines.push(decode_defline(&mut ber)?);
    }
    if !ber.at_record_end() {
        return Err(ber.error("unexpected data after Blast-def-line-set"));
    }
    Ok(DeflineSet { deflines })
}

fn decode_defline(ber: &mut Reader<'_>) -> Result<BlastDefline, BerError> {
    let end = ber.open_sequence("Blast-def-line")?;
    let mut defline = BlastDefline::default();
    let mut has_seqid = false;
    ber.members(end, "Blast-def-line", &DEFLINE_MEMBERS, |ber, member| {
        match member {
            DEFLINE_TITLE => {
                let (title, replaced) = ber.fixed_string()?;
                defline.title = Some(title);
                defline.replaced_title_bytes = replaced;
            }
            DEFLINE_SEQID => {
                has_seqid = true;
                let list_end = ber.open_sequence("the Seq-id list")?;
                while !ber.at_end(list_end)? {
                    let id = decode_seq_id(ber, &mut defline.replaced_id_bytes)?;
                    defline.ids.push(id);
                }
            }
            DEFLINE_TAXID => defline.taxid = Some(taxid_value(ber.integer()?)),
            DEFLINE_MEMBERSHIPS => defline.memberships = Some(int_list(ber, "memberships")?),
            DEFLINE_LINKS => defline.links = Some(int_list(ber, "links")?),
            _ => defline.other_info = Some(int_list(ber, "other-info")?),
        }
        Ok(())
    })?;
    if !has_seqid {
        return Err(ber.error("a Blast-def-line has no seqid"));
    }
    Ok(defline)
}

/// A `SEQUENCE OF INTEGER` member of `Blast-def-line` (NCBI's `list<int>`).
fn int_list(ber: &mut Reader<'_>, what: &str) -> Result<Vec<i32>, BerError> {
    let end = ber.open_sequence(what)?;
    let mut values = Vec::new();
    while !ber.at_end(end)? {
        values.push(ber.int(what)?);
    }
    Ok(values)
}

/// Decodes one `Seq-id`; `replaced` counts the bytes of its strings that NCBI's reader
/// replaces (`fix_visible_chars`).
fn decode_seq_id(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<ParsedSeqId, BerError> {
    let (alternative, end) = ber.open_alternative("Seq-id", 19)?;
    let id = match alternative {
        0 => ParsedSeqId::Local(decode_object_id(ber, replaced)?),
        1 => ParsedSeqId::Gibbsq(ber.int("gibbsq")?),
        2 => ParsedSeqId::Gibbmt(ber.int("gibbmt")?),
        3 => ParsedSeqId::Giim(decode_giimport_id(ber, replaced)?),
        4 => ParsedSeqId::Genbank(decode_textseq_id(ber, replaced)?),
        5 => ParsedSeqId::Embl(decode_textseq_id(ber, replaced)?),
        6 => ParsedSeqId::Pir(decode_textseq_id(ber, replaced)?),
        7 => ParsedSeqId::Swissprot(decode_textseq_id(ber, replaced)?),
        8 => ParsedSeqId::Patent(decode_patent_seq_id(ber, replaced)?),
        9 => ParsedSeqId::Other(decode_textseq_id(ber, replaced)?),
        10 => ParsedSeqId::General(decode_dbtag(ber, replaced)?),
        11 => ParsedSeqId::Gi(ber.integer()?),
        12 => ParsedSeqId::Ddbj(decode_textseq_id(ber, replaced)?),
        13 => ParsedSeqId::Prf(decode_textseq_id(ber, replaced)?),
        14 => ParsedSeqId::Pdb(decode_pdb_seq_id(ber, replaced)?),
        15 => ParsedSeqId::Tpg(decode_textseq_id(ber, replaced)?),
        16 => ParsedSeqId::Tpe(decode_textseq_id(ber, replaced)?),
        17 => ParsedSeqId::Tpd(decode_textseq_id(ber, replaced)?),
        18 => ParsedSeqId::Gpipe(decode_textseq_id(ber, replaced)?),
        _ => ParsedSeqId::NamedAnnotTrack(decode_textseq_id(ber, replaced)?),
    };
    ber.close(end, "Seq-id")?;
    Ok(id)
}

fn string(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<Vec<u8>, BerError> {
    let (value, count) = ber.fixed_string()?;
    *replaced += count;
    Ok(value)
}

fn decode_object_id(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<ObjectId, BerError> {
    let (alternative, end) = ber.open_alternative("Object-id", 1)?;
    let id = match alternative {
        0 => ObjectId::Id(ber.int("Object-id.id")?),
        _ => ObjectId::Str(string(ber, replaced)?),
    };
    ber.close(end, "Object-id")?;
    Ok(id)
}

fn decode_dbtag(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<Dbtag, BerError> {
    let end = ber.open_sequence("Dbtag")?;
    let (mut db, mut tag) = (None, None);
    ber.members(end, "Dbtag", &["db", "tag"], |ber, member| {
        match member {
            0 => db = Some(string(ber, replaced)?),
            _ => tag = Some(decode_object_id(ber, replaced)?),
        }
        Ok(())
    })?;
    match (db, tag) {
        (Some(db), Some(tag)) => Ok(Dbtag { db, tag }),
        _ => Err(ber.error("a Dbtag has no db or no tag")),
    }
}

fn decode_textseq_id(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<TextseqId, BerError> {
    let end = ber.open_sequence("Textseq-id")?;
    let mut id = TextseqId::default();
    let names = ["name", "accession", "release", "version"];
    ber.members(end, "Textseq-id", &names, |ber, member| {
        match member {
            0 => id.name = Some(string(ber, replaced)?),
            1 => id.accession = Some(string(ber, replaced)?),
            2 => id.release = Some(string(ber, replaced)?),
            _ => id.version = Some(ber.int("Textseq-id.version")?),
        }
        Ok(())
    })?;
    Ok(id)
}

fn decode_giimport_id(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<GiimportId, BerError> {
    let end = ber.open_sequence("Giimport-id")?;
    let (mut id, mut db, mut release) = (None, None, None);
    ber.members(
        end,
        "Giimport-id",
        &["id", "db", "release"],
        |ber, member| {
            match member {
                0 => id = Some(ber.int("Giimport-id.id")?),
                1 => db = Some(string(ber, replaced)?),
                _ => release = Some(string(ber, replaced)?),
            }
            Ok(())
        },
    )?;
    let Some(id) = id else {
        return Err(ber.error("a Giimport-id has no id"));
    };
    Ok(GiimportId { id, db, release })
}

fn decode_patent_seq_id(
    ber: &mut Reader<'_>,
    replaced: &mut usize,
) -> Result<PatentSeqId, BerError> {
    let end = ber.open_sequence("Patent-seq-id")?;
    let (mut seqid, mut cit) = (None, None);
    ber.members(end, "Patent-seq-id", &["seqid", "cit"], |ber, member| {
        match member {
            0 => seqid = Some(ber.int("Patent-seq-id.seqid")?),
            _ => cit = Some(decode_id_pat(ber, replaced)?),
        }
        Ok(())
    })?;
    match (seqid, cit) {
        (Some(seqid), Some(cit)) => Ok(PatentSeqId { seqid, cit }),
        _ => Err(ber.error("a Patent-seq-id has no seqid or no cit")),
    }
}

fn decode_id_pat(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<IdPat, BerError> {
    let end = ber.open_sequence("Id-pat")?;
    let (mut country, mut id, mut doc_type) = (None, None, None);
    ber.members(
        end,
        "Id-pat",
        &["country", "id", "doc-type"],
        |ber, member| {
            match member {
                0 => country = Some(string(ber, replaced)?),
                1 => {
                    let (alternative, choice_end) = ber.open_alternative("Id-pat.id", 1)?;
                    let number = string(ber, replaced)?;
                    id = Some(match alternative {
                        0 => PatentNumber::Number(number),
                        _ => PatentNumber::AppNumber(number),
                    });
                    ber.close(choice_end, "Id-pat.id")?;
                }
                _ => doc_type = Some(string(ber, replaced)?),
            }
            Ok(())
        },
    )?;
    match (country, id) {
        (Some(country), Some(id)) => Ok(IdPat {
            country,
            id,
            doc_type,
        }),
        _ => Err(ber.error("an Id-pat has no country or no id")),
    }
}

fn decode_pdb_seq_id(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<PdbSeqId, BerError> {
    let end = ber.open_sequence("PDB-seq-id")?;
    let (mut mol, mut chain, mut rel, mut chain_id) = (None, None, None, None);
    let names = ["mol", "chain", "rel", "chain-id"];
    ber.members(end, "PDB-seq-id", &names, |ber, member| {
        match member {
            0 => mol = Some(string(ber, replaced)?),
            1 => chain = Some(ber.int("PDB-seq-id.chain")?),
            2 => rel = Some(decode_date(ber, replaced)?),
            _ => chain_id = Some(string(ber, replaced)?),
        }
        Ok(())
    })?;
    let Some(mol) = mol else {
        return Err(ber.error("a PDB-seq-id has no mol"));
    };
    Ok(PdbSeqId {
        mol,
        chain,
        rel,
        chain_id,
    })
}

fn decode_date(ber: &mut Reader<'_>, replaced: &mut usize) -> Result<Date, BerError> {
    let (alternative, end) = ber.open_alternative("Date", 1)?;
    let date = match alternative {
        0 => Date::Str(string(ber, replaced)?),
        _ => {
            let std_end = ber.open_sequence("Date-std")?;
            let (mut year, mut month, mut day, mut season) = (None, None, None, None);
            let (mut hour, mut minute, mut second) = (None, None, None);
            let names = ["year", "month", "day", "season", "hour", "minute", "second"];
            ber.members(std_end, "Date-std", &names, |ber, member| {
                match member {
                    0 => year = Some(ber.int("Date-std.year")?),
                    1 => month = Some(ber.int("Date-std.month")?),
                    2 => day = Some(ber.int("Date-std.day")?),
                    3 => season = Some(string(ber, replaced)?),
                    4 => hour = Some(ber.int("Date-std.hour")?),
                    5 => minute = Some(ber.int("Date-std.minute")?),
                    _ => second = Some(ber.int("Date-std.second")?),
                }
                Ok(())
            })?;
            let Some(year) = year else {
                return Err(ber.error("a Date-std has no year"));
            };
            Date::Std(DateStd {
                year,
                month,
                day,
                season,
                hour,
                minute,
                second,
            })
        }
    };
    ber.close(end, "Date")?;
    Ok(date)
}
