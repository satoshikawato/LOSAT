//! Header record of one sequence (`.phr` / `.nhr`): the BER encoding of a
//! `Blast-def-line-set` with one definition line that carries the title, the
//! id `gnl|BL_ORD_ID|n` and the taxid 0, as makeblastdb writes it without
//! `-parse_seqids`. `crate::seqdb::header` decodes the same bytes.
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
// Without `-taxid` and `-taxid_map` the best taxid is 0
// (taxid_set.cpp:84-90, 123-128: the unset global id is returned as it is
// unless the definition line has one).
/// The header record of the sequence number `oid` of its volume: the
/// definition line with `title`, `gnl|BL_ORD_ID|oid` and taxid 0. Bytes of
/// the title that a VisibleString cannot hold are written as `#` and
/// appended to `bad`.
pub fn defline_set(title: &[u8], oid: u32, bad: &mut Vec<BadTitleChar>) -> Vec<u8> {
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
    write_int4(&mut out, 0);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    out
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
            hex(&defline_set(b"a one", 0, &mut bad)),
            "30803080a0801a0561206f6e650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            hex(&defline_set(b"c three", 1, &mut bad)),
            "30803080a0801a07632074687265650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020101000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            hex(&defline_set(b"", 0, &mut bad)),
            "30803080a0801a000000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
        );
        assert!(bad.is_empty());
    }

    #[test]
    fn bytes_outside_the_visible_range_become_number_signs() {
        let mut bad = Vec::new();
        let record = defline_set(b"a   tab  ctl \xc3\xa9 utf8 \x80hi \x7f del", 0, &mut bad);
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
            let record = defline_set(title, oid, &mut Vec::new());
            let deflines = crate::seqdb::header::decode_defline_set(&record, path).unwrap();
            assert_eq!(deflines.len(), 1);
            assert_eq!(deflines[0].title, title);
            assert_eq!(deflines[0].volume_ordinal, oid);
        }
    }
}
