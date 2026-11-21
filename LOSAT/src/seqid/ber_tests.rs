use super::{decode_full_set, read_std_signed, Reader, StdSignedError};
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

/// The taxid that both decoders read from `record` (`None`: the line has none), or the
/// message of their error, which must be the same.
fn read_taxid(record: &[u8]) -> Result<Option<i32>, String> {
    let full = decode_full_set(record).map(|set| set.deflines[0].taxid);
    let ordinal = decode_deflines_as(record, Path::new("t.phr"), 0, true);
    match (full, ordinal) {
        (Ok(taxid), Ok(deflines)) => {
            assert_eq!(deflines[0].taxid, taxid.unwrap_or(0));
            Ok(taxid)
        }
        (Err(error), Err(ordinal)) => {
            let at = format!("header record, byte {}: ", error.pos);
            assert!(ordinal.to_string().contains(&at), "{error} / {ordinal}");
            Err(error.message)
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

    let failed = [
        ("ovf9", "0209010000000000002586", "overflow error"),
        ("big5", "02050080000000", "overflow error"),
        ("nbig5", "0205ff7fffffff", "overflow error"),
        ("zero", "0200", "zero length of number"),
        ("longlen", "0281022586", "ShortLength expected"),
        ("indef", "028025860000", "ShortLength expected"),
        (
            "vstr",
            "1a024142",
            "tag 0x1a where an INTEGER (0x02) is expected",
        ),
        (
            "octet",
            "04022586",
            "tag 0x04 where an INTEGER (0x02) is expected",
        ),
        (
            "enum",
            "0a022586",
            "tag 0x0a where an INTEGER (0x02) is expected",
        ),
        (
            "cons",
            "3080020225860000",
            "tag 0x30 where an INTEGER (0x02) is expected",
        ),
        ("extra", "02022586020105", "unexpected data after"),
        ("empty", "", "tag 0x00 where an INTEGER (0x02) is expected"),
    ];
    for (name, inner, message) in failed {
        let error = read_taxid(&taxid_record(&hex(inner))).expect_err(name);
        assert!(error.contains(message), "{name}: {error}");
    }
    // `[2]` as a primitive tag before the INTEGER.
    let primitive = record_with(&hex("8202022586"), &[]);
    assert!(read_taxid(&primitive).is_err(), "prim");
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
    assert_eq!(read_taxid(&wide), Err("overflow error".to_string()));
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
    assert_eq!(
        gi("020900ffffffffffffffff").unwrap_err().message,
        "overflow error"
    );
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
    assert_eq!(read(&[0x01], 5, 4), (Err("overflow"), 1));
    assert_eq!(read(&[0x00, 0x01], 6, 4), (Err("overflow"), 2));
    assert_eq!(read(&[0xff, 0x7f], 5, 4), (Err("overflow"), 2));
    // The data ends before the value does.
    assert_eq!(read(&[0x00, 0x00], 6, 4), (Err("end"), 3));
    assert_eq!(read(&[0x25], 2, 4), (Err("end"), 2));
}
