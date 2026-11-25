use super::{
    decode_full_set, fix_visible_chars, read_std_signed, replaced_bytes, Reader, StdSignedError,
};
use crate::seqdb::header::decode_deflines_as;
use std::path::Path;

/// The header record of one sequence without parsed Seq-ids as makeblastdb writes it (the
/// title `t`, `gnl|BL_ORD_ID|0`), with `member` where makeblastdb writes the taxid `[2]`
/// and `after` (members 3 to 5) after it.
fn record_with(member: &[u8], after: &[u8]) -> Vec<u8> {
    let mut record = vec![0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x01, b't', 0, 0];
    record.extend_from_slice(&[0xa1, 0x80, 0x30, 0x80, 0xaa, 0x80, 0x30, 0x80, 0xa0, 0x80]);
    record.extend_from_slice(&[0x1a, 0x09]);
    record.extend_from_slice(b"BL_ORD_ID");
    record.extend_from_slice(&[0, 0, 0xa1, 0x80, 0xa0, 0x80, 0x02, 0x01, 0x00]);
    record.extend_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    record.extend_from_slice(member);
    record.extend_from_slice(after);
    record.extend_from_slice(&[0, 0, 0, 0]);
    record
}

/// `record_with` the inner bytes `inner` of an indefinite `[2]`.
fn taxid_record(inner: &[u8]) -> Vec<u8> {
    let mut member = vec![0xa2, 0x80];
    member.extend_from_slice(inner);
    member.extend_from_slice(&[0, 0]);
    record_with(&member, &[])
}

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).unwrap())
        .collect()
}

/// The taxid that both decoders read from `record` (`None`: the line has none), or NCBI's text
/// of the exception that both give (`NcbiException::ncbi_what`, the report that LOSAT's
/// `what()` maps to its words): its places, codes and messages, oldest first.
fn read_taxid(record: &[u8]) -> Result<Option<i32>, String> {
    let full = decode_full_set(record).map(|set| set.deflines[0].taxid);
    let ordinal = decode_deflines_as(record, Path::new("t.phr"), 0, true);
    match (full, ordinal) {
        (Ok(taxid), Ok(deflines)) => {
            assert_eq!(deflines[0].taxid, taxid.unwrap_or(0));
            Ok(taxid)
        }
        (Err(error), Err(ordinal)) => {
            let exception = error.exception().expect("NCBI's report");
            assert_eq!(ordinal.ncbi_exception(), Some(&exception));
            Err(String::from_utf8_lossy(&exception.ncbi_what()).into_owned())
        }
        (full, ordinal) => panic!("the decoders differ: {full:?} / {ordinal:?}"),
    }
}

// The 21 encodings of the taxid of probe RG-T5c (remaining-gaps-20261112/T/taxid_t5c.md):
// NCBI BLAST+ 2.17.0 reads the first nine (blastdbcmd `%T`, the reports of the record) and
// fails on the others at the record wherever it reads the header (`ReadInt4`:
// `ExpectIntegerTag`, `ReadShortLength`, `ReadStdSigned` into an `int`).
#[test]
fn every_taxid_encoding_reads_as_ncbis_read_int4_reads_it() {
    let read = [
        ("base", "02022586", Some(9606)),
        ("app42", "42022586", Some(9606)),
        ("sx9", "0209000000000000002586", Some(9606)),
        ("sx8", "02080000000000002586", Some(9606)),
        ("nonmin4", "020400002586", Some(9606)),
        ("neg1", "0201ff", Some(-1)),
    ];
    for (name, inner, taxid) in read {
        assert_eq!(read_taxid(&taxid_record(&hex(inner))), Ok(taxid), "{name}");
    }
    // `[2]` with a definite length; no `[2]`.
    let definite = record_with(&hex("a20402022586"), &[]);
    assert_eq!(read_taxid(&definite), Ok(Some(9606)), "defl");
    assert_eq!(read_taxid(&record_with(&[], &[])), Ok(None), "absent");

    // NCBI's reports (T/t5d/survey/errors_raw.txt): the reader's exception, then the member's
    // backlog with the frames; the taxid's INTEGER is at byte 55 of this record.
    let backlog = "line 752: Error: (CSerialException::{code}) \
                   ncbi::CMemberInfoFunctions::ReadWithSetFlagMember() - error while reading \
                   taxid ( at Blast-def-line-set.[].[].taxid)\n";
    let failed = [
        (
            "ovf9",
            "0209010000000000002586",
            "line 488",
            "eOverflow",
            "byte 58: overflow error",
        ),
        (
            "big5",
            "02050080000000",
            "line 499",
            "eOverflow",
            "byte 59: overflow error",
        ),
        (
            "nbig5",
            "0205ff7fffffff",
            "line 499",
            "eOverflow",
            "byte 59: overflow error",
        ),
        (
            "zero",
            "0200",
            "line 480",
            "eFormatError",
            "byte 57: zero length of number",
        ),
        (
            "longlen",
            "0281022586",
            "line 202",
            "eFormatError",
            "UnexpectedLongLength() - byte 57: ShortLength expected",
        ),
        (
            "indef",
            "028025860000",
            "line 202",
            "eFormatError",
            "UnexpectedLongLength() - byte 57: ShortLength expected",
        ),
        (
            "vstr",
            "1a024142",
            "line 377",
            "eFormatError",
            "byte 55: unexpected tag: VisibleString (26), should be: Integer (2)\n",
        ),
        (
            "octet",
            "04022586",
            "line 377",
            "eFormatError",
            "byte 55: unexpected tag: OctetString (4), should be: Integer (2)\n",
        ),
        (
            "enum",
            "0a022586",
            "line 377",
            "eFormatError",
            "byte 55: unexpected tag: Enumerated (10), should be: Integer (2)\n",
        ),
        (
            "cons",
            "3080020225860000",
            "line 377",
            "eFormatError",
            "byte 55: unexpected tag: constructed/Sequence (48), should be: Integer (2)\n",
        ),
        (
            "empty",
            "",
            "line 377",
            "eFormatError",
            "byte 55: unexpected tag: None (0), should be: Integer (2)\n",
        ),
    ];
    for (name, inner, line, code, message) in failed {
        let report = read_taxid(&taxid_record(&hex(inner))).expect_err(name);
        assert!(
            report.contains(line) && report.contains(message),
            "{name}: {report}"
        );
        let backlog = backlog.replace("{code}", code);
        assert!(report.ends_with(&backlog), "{name}: {report}");
    }
    // Data after the INTEGER, where the member's end marker is expected; `[2]` as a primitive
    // tag before the INTEGER: one exception, no backlog.
    let extra = read_taxid(&taxid_record(&hex("02022586020105"))).unwrap_err();
    assert!(
        extra.ends_with(
            "UnexpectedContinuation() - byte 59: EndOfContentsByte expected ( at \
             Blast-def-line-set.[].[].taxid)\n"
        ),
        "{extra}"
    );
    let primitive = read_taxid(&record_with(&hex("8202022586"), &[])).unwrap_err();
    assert!(
        primitive.ends_with(
            "UnexpectedTagClassByte() - byte 53: unexpected tag: contextspecific/2 (130), \
             should be contextspecific/constructed/0 (160) ( at \
             Blast-def-line-set.[].[].seqid)\n"
        ),
        "{primitive}"
    );
}

/// The other INTEGER members are read at the width of their NCBI type: the lists of
/// `Blast-def-line` (`list<int>`) by both decoders, the ordinal id (`Object-id.id`, an `int`),
/// the GI (`Int8`).
#[test]
fn integer_members_read_at_the_width_of_their_ncbi_type() {
    let taxid = hex("a280020225860000");
    // memberships `[3] SEQUENCE OF INTEGER`: a value within `int`, then one outside it.
    let lists = record_with(&taxid, &hex("a380308002010500000000"));
    assert_eq!(read_taxid(&lists), Ok(Some(9606)));
    let wide = record_with(&taxid, &hex("a38030800205008000000000000000"));
    let wide = read_taxid(&wide).unwrap_err();
    assert!(
        wide.contains("overflow error") && wide.contains("error while reading memberships"),
        "{wide}"
    );
    assert_eq!(
        decode_full_set(&lists).unwrap().deflines[0].memberships,
        Some(vec![5])
    );
    // The ordinal id: the application tag and sign extension read as NCBI reads them.
    let ordinal = |integer: &str| {
        let mut record = record_with(&taxid, &[]);
        let at = record
            .windows(3)
            .position(|window| window == [0x02, 0x01, 0x00])
            .unwrap();
        record.splice(at..at + 3, hex(integer));
        decode_deflines_as(&record, Path::new("t.phr"), 0, true)
            .map(|lines| lines[0].volume_ordinal)
    };
    assert_eq!(ordinal("420100").unwrap(), 0);
    assert_eq!(ordinal("02050000000000").unwrap(), 0);
    let error = ordinal("02050080000000").unwrap_err().to_string();
    assert!(error.contains("overflow error"), "{error}");
    // The GI, an `Int8`: nine bytes of sign extension; a value outside `Int8`.
    let gi = |integer: &str| Reader::new(&hex(integer)).int8();
    assert_eq!(gi("0209007fffffffffffffff"), Ok(i64::MAX));
    assert_eq!(gi("0209ff8000000000000000"), Ok(i64::MIN));
    assert_eq!(gi("420500ffffffff"), Ok(0xffff_ffff));
    let overflow = gi("020900ffffffffffffffff")
        .unwrap_err()
        .exception()
        .unwrap();
    assert_eq!(overflow.get_msg(), b"byte 4: overflow error");
}

/// `read_std_signed` stops where NCBI's `ReadStdSigned` throws: at the first byte that is
/// no sign extension, before it needs a byte that is not there.
#[test]
fn read_std_signed_stops_where_ncbis_throws() {
    let read = |bytes: &[u8], length: usize, width: usize| {
        let mut taken = 0;
        let result = read_std_signed(length, width, || {
            let byte = bytes.get(taken).copied().ok_or(taken);
            taken += 1;
            byte
        });
        let result = match result {
            Ok(value) => Ok(value),
            Err(StdSignedError::ZeroLength) => Err("zero"),
            Err(StdSignedError::Overflow { line: 488 }) => Err("sign"),
            Err(StdSignedError::Overflow { line: 493 }) => Err("extension"),
            Err(StdSignedError::Overflow { line: 499 }) => Err("change"),
            Err(StdSignedError::Overflow { .. }) => Err("overflow"),
            Err(StdSignedError::Read(_)) => Err("end"),
        };
        (result, taken)
    };
    assert_eq!(read(&[], 0, 4), (Err("zero"), 0));
    assert_eq!(read(&[0x80], 1, 4), (Ok(-128), 1));
    assert_eq!(read(&[0xff, 0xff, 0xff, 0xff, 0xfe], 5, 4), (Ok(-2), 5));
    assert_eq!(read(&[0xff, 0xff, 0xff, 0xff, 0xff], 5, 4), (Ok(-1), 5));
    assert_eq!(
        read(&[0x7f, 0xff, 0xff, 0xff], 4, 4),
        (Ok(i64::from(i32::MAX)), 4)
    );
    // The first byte is no sign; a later byte differs from it; the sign changes.
    assert_eq!(read(&[0x01], 5, 4), (Err("sign"), 1));
    assert_eq!(read(&[0x00, 0x01], 6, 4), (Err("extension"), 2));
    assert_eq!(read(&[0xff, 0x7f], 5, 4), (Err("change"), 2));
    // The data ends before the value does.
    assert_eq!(read(&[0x00, 0x00], 6, 4), (Err("end"), 3));
    assert_eq!(read(&[0x25], 2, 4), (Err("end"), 2));
}

/// `replaced_bytes` counts the bytes that `fix_visible_chars` replaces, for strings of every
/// byte value and of lengths about its pieces of 255 bytes, every byte replaced included.
#[test]
fn replaced_bytes_counts_the_replaced_bytes() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for len in [
        0, 1, 31, 32, 254, 255, 256, 509, 510, 511, 765, 766, 1_000, 4_097,
    ] {
        let random: Vec<u8> = (0..len).map(|_| next() as u8).collect();
        let mixed: Vec<u8> = (0..len)
            .map(|i| if i % 3 == 0 { 0x7f } else { b'a' })
            .collect();
        for value in [
            random,
            mixed,
            vec![0x80; len],
            vec![b'~'; len],
            vec![0x1f; len],
        ] {
            let replaced = fix_visible_chars(&value)
                .iter()
                .zip(&value)
                .filter(|(fixed, stored)| fixed != stored)
                .count();
            assert_eq!(replaced_bytes(&value), replaced, "length {len}");
        }
    }
}

#[test]
fn replaced_char_posts_and_negative_length_eof_match_ncbi() {
    // Minimal Blast-def-line-set with one Blast-def-line whose title holds 0x01 and 0x80,
    // followed by a missing mandatory `seqid` member (simulating a truncated header record).
    let mut truncated = vec![
        0x30, 0x80, // Blast-def-line-set SEQUENCE (indefinite)
        0x30, 0x80, // Blast-def-line SEQUENCE (indefinite)
        0xa0, 0x80, // [0] title
        0x1a, 0x05, b'a', 0x01, 0x80, b'h', b'a', // VisibleString "a\x01\x80ha"
        0x00, 0x00, // end title
        0x00, 0x00, // end Blast-def-line (missing mandatory seqid!)
        0x00, 0x00, // end set
    ];
    let posts = super::replaced_char_posts(&truncated);
    assert_eq!(posts.len(), 2);
    assert_eq!(posts[0].trace, "Blast-def-line-set.E.title");
    assert_eq!(posts[0].byte, 0x01);
    assert_eq!(posts[0].string.as_ref(), b"a\x01\x80ha");
    assert_eq!(posts[1].trace, "Blast-def-line-set.E.title");
    assert_eq!(posts[1].byte, 0x80);
    assert_eq!(posts[1].string.as_ref(), b"a#\x80ha");

    // Clean title has no posts.
    truncated[9] = b'l';
    truncated[10] = b'p';
    assert!(super::replaced_char_posts(&truncated).is_empty());

    // Negative header length EOF exception matches NCBI's strbuffer.cpp:358 + objistr.cpp:604.
    let eof = super::negative_length_eof();
    let text = String::from_utf8(eof.ncbi_what()).unwrap();
    assert!(text.contains("strbuffer.cpp\", line 358"), "{text}");
    assert!(text.contains("CEofException::eEof"), "{text}");
    assert!(text.contains("objistr.cpp\", line 604"), "{text}");
    assert!(
        text.contains("Frame type= eFrameNamed, Object type= Blast-def-line-set"),
        "{text}"
    );
}
