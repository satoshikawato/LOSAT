//! Header records of a database volume (`.phr` / `.nhr`).
//!
//! Each sequence has one BER-encoded `Blast-def-line-set`. LOSAT decodes the
//! part a database built without `-parse_seqids` uses: the title and the
//! `gnl|BL_ORD_ID|n` id. Any other Seq-id means the database stores parsed
//! Seq-ids, which LOSAT does not support yet; that is an explicit error.
//!
//! ASN.1 definitions (NCBI 598d8ae6):
//! `c++/src/objects/blastdb/blastdb.asn:31-47` (`Blast-def-line-set`,
//! `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-45` (`Seq-id`),
//! `c++/src/objects/general/general.asn:71-80` (`Dbtag`, `Object-id`).

use std::path::Path;

use super::SeqDbError;

/// One definition line of a sequence in a database without parsed Seq-ids.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Defline {
    /// The title: the whole FASTA definition line after `>`.
    pub title: Vec<u8>,
    /// `n` of `gnl|BL_ORD_ID|n` as stored, that is, relative to the volume.
    pub volume_ordinal: u32,
    /// Number of stored title bytes outside ' ' ..= '~' that were read as '#'
    /// (`fix_visible_chars`). NCBI posts an error message for each of them every
    /// time it reads the definition line; makeblastdb stores none.
    pub replaced_title_bytes: usize,
    /// `taxid` of the definition line, 0 unless `decode_defline_taxids` read it (a value
    /// above `i32::MAX`, which makeblastdb does not write, reads as `i32::MAX`).
    pub taxid: i32,
}

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
const TAG_SEQUENCE: u8 = (1 << 5) | 16;
const TAG_CONTEXT_CONSTRUCTED: u8 = (2 << 6) | (1 << 5);
const TAG_NUMBER_MASK: u8 = 31;
const INDEFINITE_LENGTH_BYTE: u8 = 0x80;
/// Deepest nesting accepted inside a skipped member.
const MAX_SKIP_DEPTH: u32 = 8;

/// Member tags of `Blast-def-line` (blastdb.asn:34-46, in order).
const DEFLINE_TITLE: u8 = 0;
const DEFLINE_SEQID: u8 = 1;
const DEFLINE_TAXID: u8 = 2;
const DEFLINE_LAST_MEMBER: u8 = 5;
/// `general` is the eleventh alternative of `Seq-id` (seqloc.asn:25-35).
const SEQID_GENERAL: u8 = 10;
/// Members of `Dbtag` and alternatives of `Object-id` (general.asn:71-80).
const DBTAG_DB: u8 = 0;
const DBTAG_TAG: u8 = 1;
const OBJECT_ID_ID: u8 = 0;

/// Where a constructed value ends.
#[derive(Clone, Copy)]
enum End {
    /// Indefinite length: two zero bytes close the value.
    Marker,
    /// Definite length: the value ends at this position.
    At(usize),
}

struct Ber<'a> {
    data: &'a [u8],
    pos: usize,
    path: &'a Path,
}

impl<'a> Ber<'a> {
    fn format(&self, message: impl Into<String>) -> SeqDbError {
        SeqDbError::Format {
            path: self.path.to_path_buf(),
            message: format!("header record, byte {}: {}", self.pos, message.into()),
        }
    }

    fn unsupported(&self, message: impl Into<String>) -> SeqDbError {
        SeqDbError::Unsupported {
            path: self.path.to_path_buf(),
            message: message.into(),
        }
    }

    fn byte(&mut self) -> Result<u8, SeqDbError> {
        let Some(&byte) = self.data.get(self.pos) else {
            return Err(self.format("the record ends inside a value"));
        };
        self.pos += 1;
        Ok(byte)
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], SeqDbError> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.data.len());
        let Some(end) = end else {
            return Err(self.format("the record ends inside a value"));
        };
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    /// Reads one identifier byte. Tag numbers above 30 (the long form) do not
    /// occur in these types.
    fn tag(&mut self) -> Result<u8, SeqDbError> {
        let tag = self.byte()?;
        if tag & TAG_NUMBER_MASK == TAG_NUMBER_MASK {
            return Err(self.format("long-form tag"));
        }
        Ok(tag)
    }

    fn expect_tag(&mut self, expected: u8, what: &str) -> Result<(), SeqDbError> {
        let tag = self.tag()?;
        if tag != expected {
            return Err(self.format(format!(
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
    fn definite_length(&mut self, first: u8) -> Result<usize, SeqDbError> {
        if first < INDEFINITE_LENGTH_BYTE {
            return Ok(usize::from(first));
        }
        let mut length_length = usize::from(first - INDEFINITE_LENGTH_BYTE);
        if length_length == 0 {
            return Err(self.format("unexpected indefinite length"));
        }
        if length_length > std::mem::size_of::<usize>() {
            return Err(self.format("length overflow"));
        }
        let mut length = usize::from(self.byte()?);
        if length == 0 {
            return Err(self.format("illegal length start"));
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
    fn open(&mut self) -> Result<End, SeqDbError> {
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
            None => Err(self.format("a value is longer than the record")),
        }
    }

    /// Whether the constructed value opened with `end` has no further
    /// element; consumes the end marker of an indefinite length.
    fn at_end(&mut self, end: End) -> Result<bool, SeqDbError> {
        match end {
            End::At(limit) => {
                if self.pos > limit {
                    return Err(self.format("a value passes the end of its container"));
                }
                Ok(self.pos == limit)
            }
            End::Marker => {
                if self.pos >= self.data.len() {
                    return Err(self.format("the record ends inside a value"));
                }
                if self.data.get(self.pos) == Some(&0) {
                    if self.data.get(self.pos + 1) != Some(&0) {
                        return Err(self.format("end-of-contents marker cut short"));
                    }
                    self.pos += 2;
                    return Ok(true);
                }
                Ok(false)
            }
        }
    }

    fn close(&mut self, end: End, what: &str) -> Result<(), SeqDbError> {
        if self.at_end(end)? {
            Ok(())
        } else {
            Err(self.format(format!("unexpected data after {what}")))
        }
    }

    /// Skips one value of any type. `depth` bounds the nesting: the skipped
    /// members of `Blast-def-line` are an INTEGER or a SEQUENCE OF INTEGER.
    fn skip_value(&mut self, depth: u32) -> Result<(), SeqDbError> {
        let tag = self.tag()?;
        if tag & (1 << 5) != 0 {
            if depth >= MAX_SKIP_DEPTH {
                return Err(self.format("values nested too deeply"));
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

    fn visible_string(&mut self) -> Result<&'a [u8], SeqDbError> {
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
    // The ordinal id is a non-negative `int`; an INTEGER wider than eight
    // content bytes is rejected instead of being reduced as above.
    fn integer(&mut self) -> Result<i64, SeqDbError> {
        self.expect_tag(TAG_INTEGER, "an INTEGER")?;
        let first = self.byte()?;
        if first >= INDEFINITE_LENGTH_BYTE {
            return Err(self.format("long length of number"));
        }
        let length = usize::from(first);
        if length == 0 {
            return Err(self.format("zero length of number"));
        }
        if length > std::mem::size_of::<i64>() {
            return Err(self.format("overflow error"));
        }
        let bytes = self.take(length)?;
        let mut n = i64::from(bytes[0] as i8);
        for &byte in &bytes[1..] {
            n = (n << 8) | i64::from(byte);
        }
        Ok(n)
    }

    fn open_context(&mut self, number: u8, what: &str) -> Result<End, SeqDbError> {
        self.expect_tag(TAG_CONTEXT_CONSTRUCTED | number, what)?;
        self.open()
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
fn fix_visible_chars(value: &[u8]) -> Vec<u8> {
    value
        .iter()
        .map(|&c| if (b' '..=b'~').contains(&c) { c } else { b'#' })
        .collect()
}

const PARSED_SEQIDS: &str = "the database stores parsed Seq-ids (makeblastdb -parse_seqids, or a \
     database distributed by NCBI); LOSAT reads databases built without -parse_seqids only";

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
/// Decodes the `Blast-def-line-set` of one sequence. An empty record gives no
/// definition line, as in NCBI.
pub fn decode_defline_set(record: &[u8], path: &Path) -> Result<Vec<Defline>, SeqDbError> {
    decode_deflines(record, path, false)
}

/// The `taxid` of each definition line of one sequence's record (0 for none), which the XML2
/// and JSON formats show (`Hit` description). It is decoded apart from `decode_defline_set`,
/// which skips it, so that no other report reads it.
pub fn decode_defline_taxids(record: &[u8], path: &Path) -> Result<Vec<i32>, SeqDbError> {
    Ok(decode_deflines(record, path, true)?
        .into_iter()
        .map(|defline| defline.taxid)
        .collect())
}

fn decode_deflines(
    record: &[u8],
    path: &Path,
    want_taxid: bool,
) -> Result<Vec<Defline>, SeqDbError> {
    let mut deflines = Vec::new();
    if record.is_empty() {
        return Ok(deflines);
    }
    let mut ber = Ber {
        data: record,
        pos: 0,
        path,
    };
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line-set")?;
    let set_end = ber.open()?;
    while !ber.at_end(set_end)? {
        deflines.push(decode_defline(&mut ber, want_taxid)?);
    }
    if ber.pos != record.len() {
        return Err(ber.format("unexpected data after Blast-def-line-set"));
    }
    Ok(deflines)
}

fn decode_defline(ber: &mut Ber<'_>, want_taxid: bool) -> Result<Defline, SeqDbError> {
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line")?;
    let defline_end = ber.open()?;
    let mut title = Vec::new();
    let mut replaced_title_bytes = 0;
    let mut volume_ordinal = None;
    let mut taxid = 0;
    let mut next_member = DEFLINE_TITLE;
    while !ber.at_end(defline_end)? {
        let tag = ber.tag()?;
        let member = tag & TAG_NUMBER_MASK;
        if tag & !TAG_NUMBER_MASK != TAG_CONTEXT_CONSTRUCTED
            || member < next_member
            || member > DEFLINE_LAST_MEMBER
        {
            return Err(ber.format(format!("tag 0x{tag:02x} is not a member of Blast-def-line")));
        }
        next_member = member + 1;
        let member_end = ber.open()?;
        match member {
            DEFLINE_TITLE => {
                let stored = ber.visible_string()?;
                replaced_title_bytes = stored.iter().filter(|c| !(b' '..=b'~').contains(c)).count();
                title = fix_visible_chars(stored);
                ber.close(member_end, "the title")?;
            }
            DEFLINE_SEQID => {
                ber.expect_tag(TAG_SEQUENCE, "the Seq-id list")?;
                let list_end = ber.open()?;
                while !ber.at_end(list_end)? {
                    if volume_ordinal.is_some() {
                        return Err(ber.unsupported(PARSED_SEQIDS));
                    }
                    volume_ordinal = Some(decode_ordinal_id(ber)?);
                }
                ber.close(member_end, "the Seq-id list")?;
            }
            // The taxid is read only when asked for (`decode_defline_taxids`).
            DEFLINE_TAXID if want_taxid => {
                taxid = i32::try_from(ber.integer()?).unwrap_or(i32::MAX);
                ber.close(member_end, "the taxid")?;
            }
            // taxid, memberships, links, other-info: not used by LOSAT.
            _ => {
                while !ber.at_end(member_end)? {
                    ber.skip_value(0)?;
                }
            }
        }
    }
    let Some(volume_ordinal) = volume_ordinal else {
        return Err(ber.unsupported(PARSED_SEQIDS));
    };
    Ok(Defline {
        title,
        volume_ordinal,
        replaced_title_bytes,
        taxid,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2231-2236
// ```c++
//                 if (seqid.Which() == CSeq_id::e_General) {
//                     CDbtag & dbt = seqid.SetGeneral();
//
//                     if (dbt.GetDb() == "BL_ORD_ID") {
//                         int vol_oid = dbt.GetTag().GetId();
//                         dbt.SetTag().SetId(m_VolStart + vol_oid);
// ```
/// Decodes one `Seq-id`, which must be `general { db "BL_ORD_ID", tag id n }`.
fn decode_ordinal_id(ber: &mut Ber<'_>) -> Result<u32, SeqDbError> {
    let tag = ber.tag()?;
    if tag != TAG_CONTEXT_CONSTRUCTED | SEQID_GENERAL {
        return Err(ber.unsupported(PARSED_SEQIDS));
    }
    let choice_end = ber.open()?;
    ber.expect_tag(TAG_SEQUENCE, "Dbtag")?;
    let dbtag_end = ber.open()?;

    let db_end = ber.open_context(DBTAG_DB, "Dbtag.db")?;
    let db = ber.visible_string()?;
    ber.close(db_end, "Dbtag.db")?;
    if db != b"BL_ORD_ID" {
        return Err(ber.unsupported(PARSED_SEQIDS));
    }

    let tag_end = ber.open_context(DBTAG_TAG, "Dbtag.tag")?;
    let object_id = ber.tag()?;
    if object_id != TAG_CONTEXT_CONSTRUCTED | OBJECT_ID_ID {
        return Err(ber.unsupported(PARSED_SEQIDS));
    }
    let id_end = ber.open()?;
    let id = ber.integer()?;
    ber.close(id_end, "Object-id.id")?;
    ber.close(tag_end, "Dbtag.tag")?;
    ber.close(dbtag_end, "Dbtag")?;
    ber.close(choice_end, "Seq-id")?;

    match u32::try_from(id) {
        Ok(ordinal) if ordinal <= i32::MAX as u32 => Ok(ordinal),
        _ => Err(ber.format(format!("BL_ORD_ID {id} is not a valid ordinal id"))),
    }
}
