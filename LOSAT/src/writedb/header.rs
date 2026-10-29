//! Header record of one sequence (`.phr` / `.nhr`): the BER encoding of a
//! `Blast-def-line-set`. Without `-parse_seqids` it has one definition line
//! with the title, the id `gnl|BL_ORD_ID|n` and the taxid 0 (`defline_set`);
//! with it, the definition lines of the sequence with every Seq-id alternative
//! (`encode_defline_set`). `crate::seqdb::header` and `crate::seqid::ber`
//! decode the same bytes.
//!
//! ASN.1 definitions (NCBI 598d8ae6):
//! `c++/src/objects/blastdb/blastdb.asn:32-47` (`Blast-def-line-set`,
//! `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-45` (`Seq-id`),
//! `c++/src/objects/general/general.asn:71-80` (`Dbtag`, `Object-id`).

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
//         eContainterTagByte      = TByte(eConstructed) | TByte(eSequence),
//         eIndefiniteLengthByte   = 0x80,
//         eEndOfContentsByte      = 0,
// ```
const TAG_INTEGER: u8 = 2;
const TAG_VISIBLE_STRING: u8 = 26;
const TAG_SEQUENCE: u8 = (1 << 5) | 16;
const TAG_CONTEXT_CONSTRUCTED: u8 = (2 << 6) | (1 << 5);
const INDEFINITE_LENGTH_BYTE: u8 = 0x80;

use crate::seqid::model::{Date, ObjectId, ParsedSeqId, PatentNumber, TextseqId};
use crate::seqid::DeflineSet;

/// Member tags of `Blast-def-line` (blastdb.asn:34-46, in order).
const DEFLINE_TITLE: u8 = 0;
const DEFLINE_SEQID: u8 = 1;
const DEFLINE_TAXID: u8 = 2;
/// `general` is the eleventh alternative of `Seq-id` (seqloc.asn:25-35).
const SEQID_GENERAL: u8 = 10;
/// Members of `Dbtag` and alternatives of `Object-id` (general.asn:71-80).
const DBTAG_DB: u8 = 0;
const DBTAG_TAG: u8 = 1;
const OBJECT_ID_ID: u8 = 0;

/// A byte of a title that a VisibleString cannot hold, as NCBI's serializer
/// reports it ("Bad char [0xC3] in string at byte 21").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BadTitleChar {
    /// The byte of the title.
    pub byte: u8,
    /// Its position in the header record (`CObjectStack::GetPosition`: the
    /// stream offset at which the replacement is written).
    pub stream_offset: usize,
}

fn open(out: &mut Vec<u8>, tag: u8) {
    out.push(tag);
    out.push(INDEFINITE_LENGTH_BYTE);
}

fn close(out: &mut Vec<u8>) {
    out.extend_from_slice(&[0, 0]);
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:393-425
// ```c++
// void CObjectOStreamAsnBinary::WriteLongLength(size_t length)
// {
//     // long form
//     size_t count;
//     if ( length <= 0xffU )
//         count = 1;
//     else if ( length <= 0xffffU )
//         count = 2;
//     else if ( length <= 0xffffffU )
//         count = 3;
//     else {
//         count = sizeof(length);
//         if ( sizeof(length) > 4 ) {
//             for ( size_t shift = (count-1)*8;
//                   count > 0; --count, shift -= 8 ) {
//                 if ( Uint1(length >> shift) != 0 )
//                     break;
//             }
//         }
//     }
//     WriteByte(TByte(0x80 + count));
//     WriteBytesOf(length, count);
// }
//
// inline
// void CObjectOStreamAsnBinary::WriteLength(size_t length)
// {
//     if ( length <= 127 )
//         WriteShortLength(length);
//     else
//         WriteLongLength(length);
// }
// ```
fn write_length(out: &mut Vec<u8>, length: usize) {
    if length <= 127 {
        out.push(length as u8);
        return;
    }
    let bytes = (length as u64).to_be_bytes();
    let count = bytes.iter().skip_while(|&&byte| byte == 0).count();
    out.push(0x80 + count as u8);
    out.extend_from_slice(&bytes[bytes.len() - count..]);
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:533-554
// ```c++
// void CObjectOStreamAsnBinary::WriteNumberValue(Int4 data)
// {
//     size_t length;
//     if ( data >= Int4(-0x80) && data <= Int4(0x7f) ) {
//         // one byte
//         length = 1;
//     }
//     else if ( data >= Int4(-0x8000) && data <= Int4(0x7fff) ) {
//         // two bytes
//         length = 2;
//     }
//     else if ( data >= Int4(-0x800000) && data <= Int4(0x7fffff) ) {
//         // three bytes
//         length = 3;
//     }
//     else {
//         // full length signed
//         length = sizeof(data);
//     }
//     WriteShortLength(length);
//     WriteBytesOf(data, length);
// }
// ```
fn write_int4(out: &mut Vec<u8>, data: i32) {
    out.push(TAG_INTEGER);
    let length = if (-0x80..=0x7f).contains(&data) {
        1
    } else if (-0x8000..=0x7fff).contains(&data) {
        2
    } else if (-0x80_0000..=0x7f_ffff).contains(&data) {
        3
    } else {
        4
    };
    out.push(length as u8);
    out.extend_from_slice(&data.to_be_bytes()[4 - length..]);
}

// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:455-459
// ```c++
// inline
// bool GoodVisibleChar(char c)
// {
//     return c >= ' ' && c <= '~';
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:842-871
// ```c++
//     WriteLength(length - fixed);
//     if ( type == eStringTypeVisible && x_FixCharsMethod() != eFNP_Allow ) {
//         size_t done = 0;
//         for ( size_t i = 0; i < length; ++i ) {
//             char c = str[i];
//             if ( !GoodVisibleChar(c) ) {
//                 ...
//                 if ( i > done ) {
//                     WriteBytes(str.data() + done, i - done);
//                 }
//                 c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, str, x_FixCharsSubst());
//                 if (c != 0) {
//                     WriteByte(c);
//                 }
//                 done = i + 1;
//             }
//         }
//         if ( done < length ) {
//             WriteBytes(str.data() + done, length - done);
//         }
//     }
//     else {
//         WriteBytes(str.data(), length);
//     }
// ```
// makeblastdb runs with the default method (replace and warn) and the
// substitute '#': every byte outside ' '..'~' is written as '#' and reported
// with the stream position of the replacement.
fn write_visible_string(out: &mut Vec<u8>, text: &[u8], bad: &mut Vec<BadTitleChar>) {
    out.push(TAG_VISIBLE_STRING);
    write_length(out, text.len());
    for &byte in text {
        if (b' '..=b'~').contains(&byte) {
            out.push(byte);
        } else {
            bad.push(BadTitleChar {
                byte,
                stream_offset: out.len(),
            });
            out.push(b'#');
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:866-893
// ```c++
//     if (OID>=0) {
//         // Re-inject the BL_ORD_ID
//         CRef<CSeq_id> gnl_id(new CSeq_id);
//         gnl_id->SetGeneral().SetDb("BL_ORD_ID");
//         gnl_id->SetGeneral().SetTag().SetId(OID);
//         CRef<CBlast_def_line_set> bdls = s_EditDeflineSet(deflines);
//         bdls->Set().front()->SetSeqid().front() = gnl_id;
//
//         deflines.Reset(&* bdls);
//     }
//     ...
//     if (bin_hdr.empty() || OID>=0) {
//         // Compress the deflines to binary.
//
//         CNcbiOstrstream oss;
//         oss << MSerial_AsnBinary << *deflines;
//         bin_hdr = CNcbiOstrstreamToString(oss);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/taxid_set.cpp:130-136
// ```c++
// void
// CTaxIdSet::FixTaxId(CRef<objects::CBlast_def_line_set> deflines)
// {
//     NON_CONST_ITERATE(CBlast_def_line_set::Tdata, itr, deflines->Set()) {
//         (*itr)->SetTaxid(x_SelectBestTaxid(**itr));
//     }
// }
// ```
// Without `-taxid` the best taxid is 0 (taxid_set.cpp:84-90, 123-128: the unset global id
// is returned as it is unless the definition line has one); with it, the global taxid
// (`TaxIdSet::best_taxid`). `-taxid_map` requires `-parse_seqids`.
/// The header record of the sequence number `oid` of its volume: the
/// definition line with `title`, `gnl|BL_ORD_ID|oid` and `taxid`. Bytes of
/// the title that a VisibleString cannot hold are written as `#` and
/// appended to `bad`.
pub fn defline_set(title: &[u8], oid: u32, taxid: i32, bad: &mut Vec<BadTitleChar>) -> Vec<u8> {
    let mut out = Vec::with_capacity(64 + title.len());
    // Blast-def-line-set ::= SEQUENCE OF Blast-def-line
    open(&mut out, TAG_SEQUENCE);
    // Blast-def-line ::= SEQUENCE
    open(&mut out, TAG_SEQUENCE);
    // title VisibleString OPTIONAL
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_TITLE);
    write_visible_string(&mut out, title, bad);
    close(&mut out);
    // seqid SEQUENCE OF Seq-id
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_SEQID);
    open(&mut out, TAG_SEQUENCE);
    // Seq-id ::= CHOICE { ..., general Dbtag, ... }
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | SEQID_GENERAL);
    // Dbtag ::= SEQUENCE { db VisibleString, tag Object-id }
    open(&mut out, TAG_SEQUENCE);
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DBTAG_DB);
    write_visible_string(&mut out, b"BL_ORD_ID", &mut Vec::new());
    close(&mut out);
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DBTAG_TAG);
    // Object-id ::= CHOICE { id INTEGER, str VisibleString }
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | OBJECT_ID_ID);
    write_int4(&mut out, oid as i32);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    // taxid INTEGER OPTIONAL
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_TAXID);
    write_int4(&mut out, taxid);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    out
}

/// A byte of a VisibleString of a parsed definition line that the serializer replaces
/// (`ReplaceVisibleChar` posts "Bad char [0xE9] in string at byte 38", the member's
/// stack path before it and the whole string after it).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadChar {
    /// `CObjectStack::GetStackTrace()`: the path of the member, as
    /// `Blast-def-line-set.E.seqid.E.local.str`.
    pub path: String,
    pub byte: u8,
    /// The stream offset of the replacement (`GetPosition`).
    pub stream_offset: usize,
    /// The string that holds it.
    pub string: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:556-597
// ```c++
// void CObjectOStreamAsnBinary::WriteNumberValue(Int8 data)
// {
//     size_t length;
//     if ( data >= -Int8(0x80) && data <= Int8(0x7f) ) {
//         // one byte
//         length = 1;
//     }
//     else if ( data >= Int8(-0x8000) && data <= Int8(0x7fff) ) {
//         // two bytes
//         length = 2;
//     }
//     ...
//     else if ( data >= NCBI_CONST_INT8(-0x80000000000000) &&
//               data <= NCBI_CONST_INT8( 0x7fffffffffffff) ) {
//         // four bytes
//         length = 7;
//     }
//     else {
//         // full length signed
//         length = sizeof(data);
//     }
//     WriteShortLength(length);
//     WriteBytesOf(data, length);
// }
// ```
// `WriteInt8` writes the universal INTEGER tag first (objostrasnb.cpp:716-724; no
// `m_CStyleBigInt` in makeblastdb's stream).
fn write_int8(out: &mut Vec<u8>, data: i64) {
    out.push(TAG_INTEGER);
    let length = (1..8usize)
        .find(|&bytes| {
            let bits = 8 * bytes as u32 - 1;
            (-(1i64 << bits)..(1i64 << bits)).contains(&data)
        })
        .unwrap_or(8);
    out.push(length as u8);
    out.extend_from_slice(&data.to_be_bytes()[8 - length..]);
}

/// The BER writer of a whole `Blast-def-line-set` (`oss << MSerial_AsnBinary <<
/// *deflines`), with the member path of each string for the serializer's diagnostics.
struct SetWriter<'a> {
    out: Vec<u8>,
    path: Vec<&'static str>,
    bad: &'a mut Vec<BadChar>,
}

impl SetWriter<'_> {
    fn open(&mut self, tag: u8) {
        open(&mut self.out, tag);
    }

    fn close(&mut self) {
        close(&mut self.out);
    }

    /// A member or an alternative `[number]` around `body`, named `name` in the path.
    fn context(&mut self, number: u8, name: &'static str, body: impl FnOnce(&mut Self)) {
        self.path.push(name);
        self.open(TAG_CONTEXT_CONSTRUCTED | number);
        body(self);
        self.close();
        self.path.pop();
    }

    fn string(&mut self, text: &[u8]) {
        let mut bad = Vec::new();
        write_visible_string(&mut self.out, text, &mut bad);
        for bad_char in bad {
            self.bad.push(BadChar {
                path: self.path.join("."),
                byte: bad_char.byte,
                stream_offset: bad_char.stream_offset,
                string: text.to_vec(),
            });
        }
    }

    fn string_member(&mut self, number: u8, name: &'static str, text: &[u8]) {
        self.context(number, name, |w| w.string(text));
    }

    fn int_member(&mut self, number: u8, name: &'static str, value: i32) {
        self.context(number, name, |w| write_int4(&mut w.out, value));
    }

    fn object_id(&mut self, tag: &ObjectId) {
        match tag {
            ObjectId::Id(id) => self.int_member(0, "id", *id),
            ObjectId::Str(text) => self.string_member(1, "str", text),
        }
    }

    fn textseq(&mut self, id: &TextseqId) {
        self.open(TAG_SEQUENCE);
        if let Some(name) = &id.name {
            self.string_member(0, "name", name);
        }
        if let Some(accession) = &id.accession {
            self.string_member(1, "accession", accession);
        }
        if let Some(release) = &id.release {
            self.string_member(2, "release", release);
        }
        if let Some(version) = id.version {
            self.int_member(3, "version", version);
        }
        self.close();
    }

    fn date(&mut self, date: &Date) {
        match date {
            Date::Str(text) => self.string_member(0, "str", text),
            Date::Std(std) => self.context(1, "std", |w| {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "year", std.year);
                if let Some(month) = std.month {
                    w.int_member(1, "month", month);
                }
                if let Some(day) = std.day {
                    w.int_member(2, "day", day);
                }
                if let Some(season) = &std.season {
                    w.string_member(3, "season", season);
                }
                if let Some(hour) = std.hour {
                    w.int_member(4, "hour", hour);
                }
                if let Some(minute) = std.minute {
                    w.int_member(5, "minute", minute);
                }
                if let Some(second) = std.second {
                    w.int_member(6, "second", second);
                }
                w.close();
            }),
        }
    }

    // The Seq-id CHOICE: each alternative is `[n]` around its value (seqloc.asn:24-45,
    // `seqid::model`), the members of a SEQUENCE in ASN.1 order, only those set.
    fn seq_id(&mut self, id: &ParsedSeqId) {
        let number = id.alternative();
        let name = id.alternative_name();
        self.context(number, name, |w| match id {
            ParsedSeqId::Local(tag) => w.object_id(tag),
            ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => {
                write_int4(&mut w.out, *value)
            }
            ParsedSeqId::Giim(giim) => {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "id", giim.id);
                if let Some(db) = &giim.db {
                    w.string_member(1, "db", db);
                }
                if let Some(release) = &giim.release {
                    w.string_member(2, "release", release);
                }
                w.close();
            }
            ParsedSeqId::Patent(patent) => {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "seqid", patent.seqid);
                w.context(1, "cit", |w| {
                    w.open(TAG_SEQUENCE);
                    w.string_member(0, "country", &patent.cit.country);
                    w.context(1, "id", |w| match &patent.cit.id {
                        PatentNumber::Number(number) => w.string_member(0, "number", number),
                        PatentNumber::AppNumber(number) => w.string_member(1, "app-number", number),
                    });
                    if let Some(doc_type) = &patent.cit.doc_type {
                        w.string_member(2, "doc-type", doc_type);
                    }
                    w.close();
                });
                w.close();
            }
            ParsedSeqId::General(dbtag) => {
                w.open(TAG_SEQUENCE);
                w.string_member(DBTAG_DB, "db", &dbtag.db);
                w.context(DBTAG_TAG, "tag", |w| w.object_id(&dbtag.tag));
                w.close();
            }
            ParsedSeqId::Gi(gi) => write_int8(&mut w.out, *gi),
            ParsedSeqId::Pdb(pdb) => {
                w.open(TAG_SEQUENCE);
                w.string_member(0, "mol", &pdb.mol);
                if let Some(chain) = pdb.chain {
                    w.int_member(1, "chain", chain);
                }
                if let Some(rel) = &pdb.rel {
                    w.context(2, "rel", |w| w.date(rel));
                }
                if let Some(chain_id) = &pdb.chain_id {
                    w.string_member(3, "chain-id", chain_id);
                }
                w.close();
            }
            _ => w.textseq(
                id.textseq()
                    .expect("the other alternatives are Textseq-ids"),
            ),
        });
    }

    fn int_list(&mut self, number: u8, name: &'static str, values: &[i32]) {
        self.context(number, name, |w| {
            w.open(TAG_SEQUENCE);
            for &value in values {
                write_int4(&mut w.out, value);
            }
            w.close();
        });
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:887-893
// ```c++
//     if (bin_hdr.empty() || OID>=0) {
//         // Compress the deflines to binary.
//
//         CNcbiOstrstream oss;
//         oss << MSerial_AsnBinary << *deflines;
//         bin_hdr = CNcbiOstrstreamToString(oss);
//     }
// ```
// The writing of blastdb.asn:32-47 (quoted in `seqid::ber`): SEQUENCE OF Blast-def-line,
// each a SEQUENCE of the members that are set, in order; `title` is a VisibleString,
// `taxid` an INTEGER, the lists SEQUENCE OF INTEGER. `seqid::ber::decode_full_set` reads it
// back.
/// The header record of a sequence whose definition lines carry parsed Seq-ids
/// (`-parse_seqids`): `set` as makeblastdb serializes it. Bytes of the strings that a
/// VisibleString cannot hold are written as `#` and appended to `bad`.
pub fn encode_defline_set(set: &DeflineSet, bad: &mut Vec<BadChar>) -> Vec<u8> {
    let mut w = SetWriter {
        out: Vec::with_capacity(64),
        path: vec!["Blast-def-line-set"],
        bad,
    };
    w.open(TAG_SEQUENCE);
    for defline in &set.deflines {
        w.path.push("E");
        w.open(TAG_SEQUENCE);
        if let Some(title) = &defline.title {
            w.string_member(DEFLINE_TITLE, "title", title);
        }
        w.context(DEFLINE_SEQID, "seqid", |w| {
            w.open(TAG_SEQUENCE);
            for id in &defline.ids {
                w.path.push("E");
                w.seq_id(id);
                w.path.pop();
            }
            w.close();
        });
        if let Some(taxid) = defline.taxid {
            w.int_member(DEFLINE_TAXID, "taxid", taxid);
        }
        if let Some(values) = &defline.memberships {
            w.int_list(3, "memberships", values);
        }
        if let Some(values) = &defline.links {
            w.int_list(4, "links", values);
        }
        if let Some(values) = &defline.other_info {
            w.int_list(5, "other-info", values);
        }
        w.close();
        w.path.pop();
    }
    w.close();
    w.out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|byte| format!("{byte:02x}")).collect()
    }

    // Records that makeblastdb 2.17.0 wrote (2026-10-12).
    #[test]
    fn the_record_of_a_short_title() {
        let mut bad = Vec::new();
        assert_eq!(
            hex(&defline_set(b"a one", 0, 0, &mut bad)),
            "30803080a0801a0561206f6e650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            hex(&defline_set(b"c three", 1, 0, &mut bad)),
            "30803080a0801a07632074687265650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020101000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            hex(&defline_set(b"", 0, 0, &mut bad)),
            "30803080a0801a000000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
        );
        assert!(bad.is_empty());
    }

    #[test]
    fn bytes_outside_the_visible_range_become_number_signs() {
        let mut bad = Vec::new();
        let record = defline_set(
            b"a   tab  ctl \xc3\xa9 utf8 \x80hi \x7f del",
            0,
            0,
            &mut bad,
        );
        assert_eq!(
            hex(&record),
            "30803080a0801a1e61202020746162202063746c20232320757466382023686920232064656c0000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
        );
        // "Bad char [0xC3] in string at byte 21", [0xA9] 22, [0x80] 29, [0x7F] 33
        let reported: Vec<(u8, usize)> = bad.iter().map(|b| (b.byte, b.stream_offset)).collect();
        assert_eq!(reported, [(0xC3, 21), (0xA9, 22), (0x80, 29), (0x7F, 33)]);
    }

    #[test]
    fn lengths_and_integers_take_the_shortest_form() {
        let mut out = Vec::new();
        write_length(&mut out, 127);
        write_length(&mut out, 128);
        write_length(&mut out, 255);
        write_length(&mut out, 256);
        write_length(&mut out, 5000);
        write_length(&mut out, 0x1_0000);
        assert_eq!(hex(&out), "7f818081ff82010082138883010000");
        let mut out = Vec::new();
        for value in [0, 127, 128, 32767, 32768, 8_388_607, 8_388_608, i32::MAX] {
            write_int4(&mut out, value);
        }
        assert_eq!(
            hex(&out),
            concat!(
                "020100",
                "02017f",
                "02020080",
                "02027fff",
                "0203008000",
                "02037fffff",
                "020400800000",
                "02047fffffff"
            )
        );
    }

    #[test]
    fn the_reader_decodes_what_the_writer_encodes() {
        let path = std::path::Path::new("test.phr");
        for (title, oid) in [
            (&b"sp|P12345|ABC_HUMAN desc"[..], 0u32),
            (&[b'x'; 5000][..], 128),
            (&b""[..], 70_000),
            (&b"t"[..], i32::MAX as u32),
        ] {
            let record = defline_set(title, oid, 0, &mut Vec::new());
            let deflines = crate::seqdb::header::decode_defline_set(&record, path).unwrap();
            assert_eq!(deflines.len(), 1);
            assert_eq!(deflines[0].title, title);
            assert_eq!(deflines[0].volume_ordinal, oid);
        }
    }

    mod parsed {
        use super::super::*;
        use crate::seqid::model::{DateStd, Dbtag, GiimportId, IdPat, PatentSeqId, PdbSeqId};
        use crate::seqid::BlastDefline;

        fn hex(bytes: &[u8]) -> String {
            bytes.iter().map(|byte| format!("{byte:02x}")).collect()
        }

        fn text(s: &str) -> Option<Vec<u8>> {
            Some(s.as_bytes().to_vec())
        }

        fn textseq(
            acc: Option<&str>,
            name: Option<&str>,
            release: Option<&str>,
            version: Option<i32>,
        ) -> TextseqId {
            TextseqId {
                name: name.and_then(text),
                accession: acc.and_then(text),
                release: release.and_then(text),
                version,
            }
        }

        fn one(title: &str, ids: Vec<ParsedSeqId>) -> DeflineSet {
            DeflineSet {
                deflines: vec![BlastDefline {
                    title: text(title),
                    ids,
                    taxid: Some(0),
                    ..BlastDefline::default()
                }],
            }
        }

        fn pdb(mol: &str, chain: Option<i32>, chain_id: Option<&str>) -> ParsedSeqId {
            ParsedSeqId::Pdb(PdbSeqId {
                mol: mol.as_bytes().to_vec(),
                chain,
                rel: None,
                chain_id: chain_id.and_then(text),
            })
        }

        fn encode(set: &DeflineSet) -> String {
            let mut bad = Vec::new();
            let record = encode_defline_set(set, &mut bad);
            assert!(bad.is_empty());
            assert_eq!(&crate::seqid::ber::decode_full_set(&record).unwrap(), set);
            hex(&record)
        }

        // Records of the .phr files makeblastdb 2.17.0 wrote with -parse_seqids (survey B
        // fixtures p01_types, p24_gi_range and p02_multi, v4; the v5 files are equal).
        #[test]
        fn records_of_makeblastdb() {
            assert_eq!(
                encode(&one(
                    "Hemoglobin subunit alpha OS=Homo sapiens",
                    vec![ParsedSeqId::Swissprot(textseq(Some("P69905"), Some("HBA_HUMAN"), Some("reviewed"), None))]
                )),
                "30803080a0801a2848656d6f676c6f62696e20737562756e697420616c706861204f533d486f6d6f2073617069656e730000a1803080a7803080a0801a094842415f48554d414e0000a1801a065036393930350000a2801a08726576696577656400000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one(
                    "gi plus ref",
                    vec![ParsedSeqId::Gi(123456), ParsedSeqId::Other(textseq(Some("NP_000550"), None, None, Some(2)))]
                )),
                "30803080a0801a0b676920706c7573207265660000a1803080ab80020301e2400000a9803080a1801a094e505f3030303535300000a38002010200000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("pdb chain A", vec![pdb("1ABC", Some(65), Some("A"))])),
                "30803080a0801a0b70646220636861696e20410000a1803080ae803080a0801a04314142430000a1800201410000a3801a014100000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("pdb two-letter chain", vec![pdb("2XYZ", None, Some("AA"))])),
                "30803080a0801a147064622074776f2d6c657474657220636861696e0000a1803080ae803080a0801a043258595a0000a3801a02414100000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("pdb no chain", vec![pdb("3DEF", None, None)])),
                "30803080a0801a0c706462206e6f20636861696e0000a1803080ae803080a0801a043344454600000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("local numeric", vec![ParsedSeqId::Local(ObjectId::Id(123))])),
                "30803080a0801a0d6c6f63616c206e756d657269630000a1803080a080a08002017b0000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one(
                    "general numeric tag",
                    vec![ParsedSeqId::General(Dbtag { db: b"mydb".to_vec(), tag: ObjectId::Id(456) })]
                )),
                "30803080a0801a1367656e6572616c206e756d65726963207461670000a1803080aa803080a0801a046d7964620000a180a080020201c8000000000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one(
                    "patent",
                    vec![ParsedSeqId::Patent(PatentSeqId {
                        seqid: 1,
                        cit: IdPat {
                            country: b"US".to_vec(),
                            id: PatentNumber::Number(b"1234567".to_vec()),
                            doc_type: None,
                        },
                    })]
                )),
                "30803080a0801a06706174656e740000a1803080a8803080a0800201010000a1803080a0801a0255530000a180a0801a073132333435363700000000000000000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("gi 2^31", vec![
                    ParsedSeqId::Gi(2_147_483_648),
                    ParsedSeqId::Other(textseq(Some("NP_000011"), None, None, Some(1))),
                ])),
                "30803080a0801a07676920325e33310000a1803080ab80020500800000000000a9803080a1801a094e505f3030303031310000a38002010100000000000000000000a280020100000000000000"
            );
            assert_eq!(
                encode(&one("gi 2^32-1", vec![ParsedSeqId::Gi(4_294_967_295)])),
                "30803080a0801a09676920325e33322d310000a1803080ab80020500ffffffff000000000000a280020100000000000000"
            );
            let refseq = |acc: &str| ParsedSeqId::Other(textseq(Some(acc), None, None, Some(1)));
            let three = DeflineSet {
                deflines: [
                    ("first title", refseq("NP_100003")),
                    ("second title", refseq("NP_100004")),
                    ("third title", ParsedSeqId::Gi(1005)),
                ]
                .into_iter()
                .map(|(title, id)| BlastDefline {
                    title: text(title),
                    ids: vec![id],
                    taxid: Some(0),
                    ..BlastDefline::default()
                })
                .collect(),
            };
            assert_eq!(
                encode(&three),
                "30803080a0801a0b6669727374207469746c650000a1803080a9803080a1801a094e505f3130303030330000a38002010100000000000000000000a280020100000000003080a0801a0c7365636f6e64207469746c650000a1803080a9803080a1801a094e505f3130303030340000a38002010100000000000000000000a280020100000000003080a0801a0b7468697264207469746c650000a1803080ab80020203ed000000000000a280020100000000000000"
            );
        }

        #[test]
        fn every_alternative_reads_back() {
            let ts = textseq(Some("A1"), Some("n"), Some("r"), Some(3));
            let ids = vec![
                ParsedSeqId::Local(ObjectId::Str(b"x".to_vec())),
                ParsedSeqId::Gibbsq(7),
                ParsedSeqId::Gibbmt(-8),
                ParsedSeqId::Giim(GiimportId {
                    id: 9,
                    db: text("d"),
                    release: text("r"),
                }),
                ParsedSeqId::Genbank(ts.clone()),
                ParsedSeqId::Embl(ts.clone()),
                ParsedSeqId::Pir(ts.clone()),
                ParsedSeqId::Swissprot(ts.clone()),
                ParsedSeqId::Patent(PatentSeqId {
                    seqid: 2,
                    cit: IdPat {
                        country: b"EP".to_vec(),
                        id: PatentNumber::AppNumber(b"9".to_vec()),
                        doc_type: text("A"),
                    },
                }),
                ParsedSeqId::Other(ts.clone()),
                ParsedSeqId::General(Dbtag {
                    db: b"db".to_vec(),
                    tag: ObjectId::Str(b"t".to_vec()),
                }),
                ParsedSeqId::Gi(i64::MAX),
                ParsedSeqId::Gi(-1),
                ParsedSeqId::Ddbj(ts.clone()),
                ParsedSeqId::Prf(ts.clone()),
                ParsedSeqId::Pdb(PdbSeqId {
                    mol: b"1ABC".to_vec(),
                    chain: Some(32),
                    rel: Some(Date::Std(DateStd {
                        year: 2020,
                        month: Some(1),
                        day: None,
                        season: text("s"),
                        hour: None,
                        minute: Some(5),
                        second: None,
                    })),
                    chain_id: None,
                }),
                ParsedSeqId::Pdb(PdbSeqId {
                    mol: b"2ABC".to_vec(),
                    chain: None,
                    rel: Some(Date::Str(b"x".to_vec())),
                    chain_id: None,
                }),
                ParsedSeqId::Tpg(ts.clone()),
                ParsedSeqId::Tpe(ts.clone()),
                ParsedSeqId::Tpd(ts.clone()),
                ParsedSeqId::Gpipe(ts.clone()),
                ParsedSeqId::NamedAnnotTrack(ts.clone()),
            ];
            let set = DeflineSet {
                deflines: vec![
                    BlastDefline {
                        title: text(&"t".repeat(300)),
                        ids,
                        taxid: Some(9606),
                        memberships: Some(vec![1, 2]),
                        links: Some(vec![]),
                        other_info: Some(vec![-5]),
                        ..BlastDefline::default()
                    },
                    BlastDefline {
                        title: None,
                        ids: vec![ParsedSeqId::Gi(1)],
                        ..BlastDefline::default()
                    },
                ],
            };
            encode(&set);
        }

        // makeblastdb 2.17.0 on survey B fixture p23_nonascii (record 1, `>lcl|\xe9abc
        // latin-1 first byte`; record 3, `>ref|NP_000001.1| title with \xe9 byte`).
        #[test]
        fn bad_characters_are_reported_with_their_member() {
            let mut bad = Vec::new();
            let set = one(
                "latin-1 first byte",
                vec![ParsedSeqId::Local(ObjectId::Str(b"\xe9abc".to_vec()))],
            );
            let record = encode_defline_set(&set, &mut bad);
            assert_eq!(
                bad,
                [BadChar {
                    path: "Blast-def-line-set.E.seqid.E.local.str".to_string(),
                    byte: 0xe9,
                    stream_offset: 38,
                    string: b"\xe9abc".to_vec(),
                }]
            );
            assert_eq!(record[38], b'#');
            bad.clear();
            let mut set = one(
                "",
                vec![ParsedSeqId::Other(textseq(
                    Some("NP_000001"),
                    None,
                    None,
                    Some(1),
                ))],
            );
            set.deflines[0].title = Some(b"title with \xe9 byte".to_vec());
            encode_defline_set(&set, &mut bad);
            assert_eq!(bad.len(), 1);
            assert_eq!(
                (bad[0].path.as_str(), bad[0].stream_offset),
                ("Blast-def-line-set.E.title", 19)
            );
        }
    }
}
