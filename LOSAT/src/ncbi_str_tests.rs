use super::*;
use crate::utils::byte_search::tests::{search_needles, search_texts};

// NCBI BLAST+ 2.17.0, probes of line RG D1a-3 (outfmt 0 prolog and epilog, blastdbcmd
// `-info`) and the copy of `s_PrintBase10Uint8` in D/logs/survey_d1a3/sim/commas.c.
#[test]
fn commas_lose_digits_where_ncbis_chunks_do() {
    assert_eq!(with_commas(0), "0");
    assert_eq!(with_commas(999), "999");
    assert_eq!(with_commas(1_000), "1,000");
    assert_eq!(with_commas(186_063_409), "186,063,409");
    assert_eq!(with_commas(u64::from(u32::MAX)), "4,294,967,295");
    assert_eq!(with_commas(1 << 32), "4,294,967,296");
    assert_eq!(
        with_commas(1_000_000_000_000_000_000),
        "1,000,000,000,000,000,000"
    );
    // p5_pr_L4294e15#0: the last value printed whole.
    assert_eq!(
        with_commas(4_294_967_295_999_999_999),
        "4,294,967,295,999,999,999"
    );
    // p5_pr_L4295e15#0 and #6.
    assert_eq!(
        with_commas(4_294_967_296_000_000_000),
        "494,967,296,000,000,000"
    );
    // p5_neg_m1#28 (blastdbcmd, `UInt8ToString` of 2^64 - 1).
    assert_eq!(with_commas(u64::MAX), "1,846,744,073,709,551,615");
    // p5_pr_Lmax#0.
    assert_eq!(int8_with_commas(i64::MAX), "923,372,036,854,775,807");
    // p5_neg_2p63#16 and p5_neg_2p63p1#16.
    assert_eq!(int8_with_commas(i64::MIN), "-923,372,036,854,775,808");
    assert_eq!(int8_with_commas(i64::MIN + 1), "-923,372,036,854,775,807");
    // p5_neg_m1#16 and p5_pr_sum2#0 (two names of 5e18: the Int8 sum wraps).
    assert_eq!(int8_with_commas(-1), "-1");
    assert_eq!(
        int8_with_commas(10_000_000_000_000_000_000u64 as i64),
        "-846,744,073,709,551,616"
    );
    // p5_pr_sum2#0 (one name of 5e18 in the epilog).
    assert_eq!(
        int8_with_commas(5_000_000_000_000_000_000),
        "500,000,000,000,000,000"
    );
    assert_eq!(
        int8_with_commas(13_000_000_000_000_000_000u64 as i64),
        "-546,744,073,709,551,616"
    );
    assert_eq!(
        with_commas(12_345_678_901_234_567_890),
        "1,245,678,901,234,567,890"
    );
    assert_eq!(with_commas(4_294_967_296_999), "4,294,967,296,999");
    assert_eq!(int8_with_commas(-1_294_967_296), "-1,294,967,296");
    assert_eq!(int8_with_commas(4_000_000_000), "4,000,000,000");
    // p5_pr_Nm5#0, p5_pr_N2p31#0, p5_pr_N3g#0.
    assert_eq!(int_with_commas(-5), "-5");
    assert_eq!(int_with_commas(i32::MIN), "-2,147,483,648");
    assert_eq!(int_with_commas(-1_294_967_296), "-1,294,967,296");
    assert_eq!(int_with_commas(0), "0");
    assert_eq!(int_with_commas(23), "23");
}

/// The message, the line and the position of an error (NCBI's `(m_Pos = <n>)`, which LOSAT
/// writes as `(position <n>)`).
fn failure(e: NcbiException) -> (String, u32, String) {
    let what = String::from_utf8(e.what()).unwrap();
    let pos = what
        .rsplit_once(crate::wording::POSITION)
        .map(|(_, rest)| rest.trim_end_matches(")\n").to_string())
        .unwrap_or_default();
    (
        String::from_utf8(e.get_msg().to_vec()).unwrap(),
        e.site().line,
        pos,
    )
}

// makeblastdb 2.17.0 (the official Linux x64 release), probe of 2026-11-08 (cases
// mapsz_*_2seq, txn_*_2seq; taxid_map_*): the line and position of each conversion error.
#[test]
fn errors_are_thrown_where_ncbi_throws_them() {
    let cases: [(Result<(), NcbiException>, &str, u32, &str); 9] = [
        (
            string_to_uint8(b"abc").map(drop),
            "Cannot convert string 'abc' to Uint8",
            908,
            "0",
        ),
        (
            string_to_uint8(b"1x").map(drop),
            "Cannot convert string '1x' to Uint8",
            908,
            "1",
        ),
        (
            string_to_uint8(b"+").map(drop),
            "Cannot convert string '+' to Uint8",
            895,
            "1",
        ),
        (
            string_to_uint8(b"99999999999999999999").map(drop),
            "Cannot convert string '99999999999999999999' to Uint8",
            912,
            "19",
        ),
        (
            string_to_int(b"abc", 0).map(drop),
            "Cannot convert string 'abc' to Int8",
            861,
            "0",
        ),
        (
            string_to_int(b"-1x", 0).map(drop),
            "Cannot convert string '-1x' to Int8",
            870,
            "2",
        ),
        (
            string_to_int(b"99999999999", 0).map(drop),
            "Cannot convert string '99999999999' to int, overflow",
            640,
            "0",
        ),
        (
            string_to_int(b"99999999999999999999", 0).map(drop),
            "Cannot convert string '99999999999999999999' to Int8, overflow",
            852,
            "18",
        ),
        (
            string_to_uint(b"4294967296").map(drop),
            "Cannot convert string '4294967296' to unsigned int, overflow",
            652,
            "0",
        ),
    ];
    for (result, message, line, pos) in cases {
        assert_eq!(
            failure(result.unwrap_err()),
            (message.to_string(), line, pos.to_string()),
            "{message}"
        );
    }
    let e = string_to_uint8(b"abc").unwrap_err();
    assert_eq!(
        e.site(),
        ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 908,
            kind: "CStringException::eConvert",
            function: Some("ncbi::NStr::StringToUInt8"),
        }
    );
}

// NStr::StringToInt / StringToUInt / StringToUInt8 with the default flags (makeblastdb
// 2.17.0: MAX_LMDB_TXN_ENTRY="+3" is 3, " 3" and "abc" stop it, "3000000000" overflows
// int, BLASTDB_LMDB_MAP_SIZE="-5" stops it).
#[test]
fn values_convert_as_ncbi_converts_them() {
    assert_eq!(string_to_int(b"+3", 0), Ok(3));
    assert_eq!(string_to_int(b"-1", 0), Ok(-1));
    assert_eq!(string_to_int(b"-2147483648", 0), Ok(i32::MIN));
    assert!(string_to_int(b" 3", 0).is_err());
    assert_eq!(string_to_int(b" 3", ALLOW_LEADING_SPACES), Ok(3));
    assert!(string_to_int(b"3 ", ALLOW_LEADING_SPACES).is_err());
    assert!(string_to_int(b"", 0).is_err());
    assert_eq!(string_to_int8(b"-9223372036854775808", 0), Ok(i64::MIN));
    assert_eq!(string_to_int8(b"12\x0034", 0), Ok(12));
    // Base 16 (ncbistr.cpp:684-712, 765-795): `0x` or `0X` is skipped, a letter of either case
    // is a digit, and a `0x` with no digit after it is 0.
    assert_eq!(string_to_int8_base(b"0x10", 0, 16), Ok(16));
    assert_eq!(string_to_int8_base(b"0XfF", 0, 16), Ok(255));
    assert_eq!(string_to_int8_base(b"ff", 0, 16), Ok(255));
    assert_eq!(string_to_int8_base(b"0x", 0, 16), Ok(0));
    assert_eq!(string_to_int8_base(b"-0x10", 0, 16), Ok(-16));
    assert_eq!(
        string_to_int8_base(b"0x7fffffffffffffff", 0, 16),
        Ok(i64::MAX)
    );
    assert_eq!(
        failure(string_to_int8_base(b"0x8000000000000000", 0, 16).unwrap_err()),
        (
            "Cannot convert string '0x8000000000000000' to Int8, overflow".to_string(),
            852,
            "17".to_string()
        )
    );
    assert_eq!(
        failure(string_to_int8_base(b"0x1g", 0, 16).unwrap_err()),
        (
            "Cannot convert string '0x1g' to Int8".to_string(),
            870,
            "3".to_string()
        )
    );
    // fAllowLeadingSpaces | fAllowTrailingSpaces (a taxid list's entry); C `isspace`
    // includes the vertical tab.
    let both = ALLOW_LEADING_SPACES | ALLOW_TRAILING_SPACES;
    assert_eq!(string_to_int(b" \x0b-5\t\r ", both), Ok(-5));
    assert_eq!(string_to_int(b"7 \x00x", both), Ok(7));
    for bad in [&b"1 2"[..], b"9606.0", b"+", b" ", b"0x10", b"3 -"] {
        assert!(string_to_int(bad, both).is_err(), "{bad:?}");
    }
    assert!(string_to_int(b"3 ", ALLOW_TRAILING_SPACES).is_ok());
    assert!(string_to_int(b" 3", ALLOW_TRAILING_SPACES).is_err());
    assert_eq!(string_to_uint(b"4294967295"), Ok(u32::MAX));
    assert_eq!(string_to_uint(b"+0"), Ok(0));
    assert!(string_to_uint8(b"-5").is_err());
    assert!(string_to_uint8(b"").is_err());
    assert!(string_to_uint8(b" 1").is_err());
    assert!(string_to_uint8(b"\xb11").is_err());
    assert_eq!(string_to_uint8(b"+40960"), Ok(40960));
    assert_eq!(string_to_uint8(b"18446744073709551615"), Ok(u64::MAX));
    assert_eq!(
        failure(string_to_uint8(b"18446744073709551616").unwrap_err()),
        (
            "Cannot convert string '18446744073709551616' to Uint8".to_string(),
            912,
            "19".to_string()
        )
    );
}

// `NStr::StringToInt(str, fConvErr_NoThrow | fAllowTrailingSymbols)` (the version of
// `CSeqDBVol::x_CheckVersions`, `seqdb::isam::check_versions`): the bytes after the digits
// are skipped, NUL included; `fConvErr_NoThrow` gives 0 for an error.
#[test]
fn trailing_symbols_skip_the_rest_of_the_string() {
    let no_throw = |s: &[u8]| string_to_int(s, ALLOW_TRAILING_SYMBOLS).unwrap_or(0);
    for (text, value) in [
        (&b"12"[..], 12),
        (b"12|HBA_HUMAN", 12),
        (b"-5x", -5),
        (b"+7", 7),
        (b"3 4", 3),
        (b"9\x00x", 9),
        (b"", 0),
        (b"x1", 0),
        (b" 1", 0),
        (b"2147483647", i32::MAX),
        (b"2147483648", 0),
        (b"-2147483648", i32::MIN),
        (b"99999999999999999999", 0),
        (b"-9223372036854775808", 0),
    ] {
        assert_eq!(no_throw(text), value, "{text:?}");
    }
    // `fAllowTrailingSpaces` alone still stops at a symbol.
    assert!(string_to_int(b"12|x", ALLOW_TRAILING_SPACES).is_err());
    assert_eq!(string_to_int(b"12 \t", ALLOW_TRAILING_SYMBOLS), Ok(12));
    assert_eq!(string_to_int8(b"-12abc", ALLOW_TRAILING_SYMBOLS), Ok(-12));
    // The overflow is still found before the trailing symbols.
    assert_eq!(
        failure(string_to_int8(b"99999999999999999999x", ALLOW_TRAILING_SYMBOLS).unwrap_err()),
        (
            "Cannot convert string '99999999999999999999x' to Int8, overflow".to_string(),
            852,
            "18".to_string()
        )
    );
}

// `NStr::StringToInt(str)` with the default flags, as `seqdb::mask` reads the mask
// algorithm ids and `s_IsNumericId` tests a program: a sign and digits only, up to the
// first NUL.
#[test]
fn default_flags_take_a_sign_and_digits_only() {
    for (text, value) in [
        (&b"0"[..], Some(0)),
        (b"30", Some(30)),
        (b"+30", Some(30)),
        (b"-5", Some(-5)),
        (b"007", Some(7)),
        (b"2147483647", Some(i32::MAX)),
        (b"-2147483648", Some(i32::MIN)),
        (b"2147483648", None),
        (b"-2147483649", None),
        (b"99999999999999999999999999", None),
        (b"", None),
        (b"+", None),
        (b"-", None),
        (b" 5", None),
        (b"5 ", None),
        (b"5x", None),
        (b"0x10", None),
        (b"1,000", None),
        (b"--5", None),
        // The text ends at a NUL, as the C string does.
        (b"5\x00x", Some(5)),
        (b"\x005", None),
    ] {
        assert_eq!(string_to_int(text, 0).ok(), value, "{text:?}");
    }
    // `NStr::StringToNumeric<TGi>(acc)` of `CSeq_id::Set` (`seqid::parse`).
    assert_eq!(string_to_int8(b"4294967296", 0).ok(), Some(4_294_967_296));
    assert_eq!(string_to_int8(b"-5", 0).ok(), Some(-5));
    assert_eq!(string_to_int8(b"9223372036854775808", 0).ok(), None);
    assert_eq!(
        string_to_int8(b"-9223372036854775808", 0).ok(),
        Some(i64::MIN)
    );
    assert_eq!(string_to_int8(b"12\x0034", 0).ok(), Some(12));
    // `NStr::StringToUInt(str)` of the string ISAM's data pages (`seqdb::isam`).
    assert_eq!(string_to_uint(b"+5"), Ok(5));
    assert!(string_to_uint(b"-5").is_err());
    assert!(string_to_uint(b"4294967296").is_err());
    assert!(string_to_uint(b"").is_err());
    assert!(string_to_uint(b"5\x00").is_err());
}

// `NStr::StringToNumeric<TTaxId>(id, fAllowLeadingSpaces | fAllowTrailingSpaces)` of a taxid
// list entry (`blastinput::db_filter`, `blastdbcmd::taxonomy`).
#[test]
fn taxid_list_entries_take_spaces_around_the_number() {
    let both = ALLOW_LEADING_SPACES | ALLOW_TRAILING_SPACES;
    for (text, value) in [
        (&b"9606"[..], Some(9606)),
        (b"+9606", Some(9606)),
        (b"  -5 ", Some(-5)),
        (b"2147483647", Some(i32::MAX)),
        (b"-2147483648", Some(i32::MIN)),
        (b"\x0b42\x0b", Some(42)),
        (b"abc", None),
        (b"# human", None),
        (b"2147483648", None),
        (b"9606.0", None),
        (b"+", None),
        (b"1 2", None),
        (b"0x10", None),
    ] {
        assert_eq!(
            string_to_int(text, both).ok(),
            value,
            "{}",
            String::from_utf8_lossy(text)
        );
    }
}

/// The taxid of a `-taxid_map` line (`writedb::taxid_set`): `StringToNumeric<TTaxId>` with
/// leading spaces allowed; its errors are NCBI's (t03_taxmap_badnum, x_taxmap_trail,
/// x_taxmap_minus, x_taxmap_overflow).
#[test]
fn taxid_map_values_take_leading_spaces() {
    let taxid = |s: &[u8]| string_to_int(s, ALLOW_LEADING_SPACES);
    assert_eq!(taxid(b"9606"), Ok(9606));
    assert_eq!(taxid(b"+9"), Ok(9));
    assert_eq!(taxid(b"-8"), Ok(-8));
    assert_eq!(taxid(b"\x0c19"), Ok(19));
    assert_eq!(taxid(b"12\x0034"), Ok(12));
    for (text, line, function, message) in [
        (
            &b"abc"[..],
            861,
            "StringToInt8",
            "Cannot convert string 'abc' to Int8 (m_Pos = 0)",
        ),
        (
            b"12ab",
            870,
            "StringToInt8",
            "Cannot convert string '12ab' to Int8 (m_Pos = 2)",
        ),
        (
            b"-",
            861,
            "StringToInt8",
            "Cannot convert string '-' to Int8 (m_Pos = 1)",
        ),
        (
            b"3000000000",
            640,
            "StringToInt",
            "Cannot convert string '3000000000' to int, overflow (m_Pos = 0)",
        ),
    ] {
        let error = taxid(text).unwrap_err();
        assert_eq!(error.site().line, line);
        assert_eq!(
            error.site().function.map(str::to_string),
            Some(format!("ncbi::NStr::{function}"))
        );
        // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
        assert_eq!(
            String::from_utf8(error.what()).unwrap(),
            crate::wording::losat_text(&format!(
                "NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line {line}: Error: \
                     (CStringException::eConvert) ncbi::NStr::{function}() - {message}\n"
            ))
        );
    }
}

/// `NStr::StringToNonNegativeInt` (`CObject_id::SetStrOrId`, the versions of `CSeq_id`,
/// blastdbcmd's mask algorithm of a batch line): digits after an optional `+`, at most
/// `INT_MAX`, else -1; it reads every byte, so a NUL is not a digit.
#[test]
fn non_negative_ints_are_digits_only() {
    for (text, value) in [
        (&b"+12"[..], 12),
        (b"+", -1),
        (b"", -1),
        (b"0", 0),
        (b"2147483647", i32::MAX),
        (b"2147483648", -1),
        (b"9999999999", -1),
        (b"-1", -1),
        (b" 1", -1),
        (b"1\x00", -1),
        (b"\xb11", -1),
        (b"++1", -1),
    ] {
        assert_eq!(string_to_non_negative_int(text), value, "{text:?}");
    }
}

// ncbistr.cpp:3736-3844 (`s_PrintableString`, cited at `printable_string`): C escapes for
// the quotes, the backslash, the control characters and `??` (`?` only next to another `?`),
// octal for the other bytes outside 0x20..0x7e, in as few digits as the next byte allows
// (blastdbcmd's messages and rpsbproc's data loader quote with it; makeblastdb's titles are
// kept in it).
#[test]
fn printable_string_escapes_as_ncbi() {
    let printable =
        |text: &str| String::from_utf8(printable_string(text.as_bytes()).into_owned()).unwrap();
    assert_eq!(printable("a'b\"c\\d"), "a\\'b\\\"c\\\\d");
    assert_eq!(printable("x\ty"), "x\\ty");
    assert_eq!(printable("1\n2"), "1\\n2");
    assert_eq!(printable("a??b?c"), "a\\?\\?b?c");
    assert_eq!(printable("\u{1}9"), "\\19");
    assert_eq!(printable("\u{1}1"), "\\0011");
    assert_eq!(printable("é"), "\\303\\251");
    assert_eq!(printable("a&b"), "a&b");
    assert_eq!(printable("\u{7f}"), "\\177");
    assert_eq!(&*printable_string(b">plain title 1"), b">plain title 1");
    assert!(matches!(printable_string(b">a&b?c"), Cow::Borrowed(_)));
    assert_eq!(&*printable_string(b"'\"\\"), b"\\'\\\"\\\\");
    assert_eq!(&*printable_string(b"??x?"), b"\\?\\?x?");
    assert_eq!(
        &*printable_string(b"\x07\x08\x0c\r\x0b\n"),
        b"\\a\\b\\f\\r\\v\\n"
    );
    // Octal: as few digits as the next byte allows.
    assert_eq!(&*printable_string(b"\x01x"), b"\\1x");
    assert_eq!(&*printable_string(b"\x017"), b"\\0017");
    assert_eq!(&*printable_string(b"\x00"), b"\\0");
    assert_eq!(&*printable_string(b"\x1f8"), b"\\378");
    assert_eq!(&*printable_string(b"\x7f\x80\xff"), b"\\177\\200\\377");
}

#[test]
fn parse_escapes_reads_printable_string_back() {
    let all: Vec<u8> = (0..=255u8).chain([b'1', b'?', b'?', 0, b'7']).collect();
    assert_eq!(parse_escapes(&printable_string(&all)).unwrap(), all);
}

/// `NStr::IsBlank` is C `isspace` over the whole string: the six spaces of the C locale
/// (vertical tab and form feed included), and an empty string, are blank; NUL is not.
#[test]
fn is_blank_takes_the_spaces_of_c_isspace() {
    assert!(is_blank(b""));
    assert!(is_blank(b" \t\n\x0b\x0c\r"));
    for text in [&b"\0"[..], b" 9 ", b"\x85", b"\xa0"] {
        assert!(!is_blank(text), "{text:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892 (`NStr::ParseEscapes`,
// cited at `parse_escapes`): the C escapes, `\x` with every hexadecimal digit after it (the
// value cut to a byte), up to three octal digits (cut to a byte), an escaped line feed that
// is no byte, any other byte for itself; and the three errors.
#[test]
fn parse_escapes_reads_ncbis_escapes() {
    assert_eq!(parse_escapes(b"a\\tb\\n").unwrap(), b"a\tb\n");
    assert_eq!(parse_escapes(b"a\\x41\\x4a4bz").unwrap(), b"aA\x4bz");
    assert_eq!(parse_escapes(b"\\x000000041").unwrap(), b"A");
    assert_eq!(parse_escapes(b"\\101\\1012\\8\\q").unwrap(), b"AA28q");
    assert_eq!(parse_escapes(b"\\777\\0\\08").unwrap(), b"\xff\0\08");
    assert_eq!(
        parse_escapes(b"\\a\\b\\f\\r\\v\\\\\\\"").unwrap(),
        b"\x07\x08\x0c\r\x0b\\\""
    );
    assert_eq!(parse_escapes(b"a\\\nb").unwrap(), b"ab");
    assert_eq!(parse_escapes(b"").unwrap(), b"");
    assert_eq!(parse_escapes(b"a\\"), Err("Unterminated escape sequence"));
    for bad in [&b"a\\xg"[..], b"\\xZZ", b"\\x"] {
        assert_eq!(
            parse_escapes(bad),
            Err("\\x followed by no hexadecimal digits"),
            "{bad:?}"
        );
    }
    assert_eq!(parse_escapes(b"\\xffffffff").unwrap(), b"\xff");
    assert!(parse_escapes(b"\\x100000000").is_err());
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3145-3180 (cited at
// `truncate_spaces`).
#[test]
fn truncate_spaces_cuts_the_c_spaces_at_both_ends() {
    assert_eq!(truncate_spaces(b""), b"");
    assert_eq!(truncate_spaces(b" \t\n\x0b\x0c\r"), b"");
    assert_eq!(truncate_spaces(b"\x0b a b\r\n"), b"a b");
    assert_eq!(truncate_spaces(b"a"), b"a");
    assert_eq!(truncate_spaces(b"\0 \xa0"), b"\0 \xa0");
    assert_eq!(truncate_spaces(b" \x0b\x0cid x\t\r "), b"id x");
    assert_eq!(truncate_spaces(b"\xa0id\xa0"), b"\xa0id\xa0");
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3145-3180 (cited at
// `truncate_spaces`): `eTrunc_Begin` and `eTrunc_End` cut one end each, and the `string` form
// cuts what the bytes' form cuts.
#[test]
fn truncate_spaces_cuts_one_end_for_begin_and_end() {
    assert_eq!(truncate_spaces_begin(b""), b"");
    assert_eq!(truncate_spaces_begin(b" \t\n\x0b\x0c\r"), b"");
    assert_eq!(truncate_spaces_begin(b"\x0b a b\r\n"), b"a b\r\n");
    assert_eq!(truncate_spaces_begin(b"\xa0 a"), b"\xa0 a");
    assert_eq!(truncate_spaces_end(b""), b"");
    assert_eq!(truncate_spaces_end(b" \t\n\x0b\x0c\r"), b"");
    assert_eq!(truncate_spaces_end(b"\x0b a b\r\n"), b"\x0b a b");
    assert_eq!(truncate_spaces_end(b"a \xa0"), b"a \xa0");
    assert_eq!(truncate_spaces_str(""), "");
    assert_eq!(truncate_spaces_str(" \t\n\x0b\x0c\r"), "");
    assert_eq!(truncate_spaces_str("\x0b é b\r\n"), "é b");
    assert_eq!(truncate_spaces_str("\u{a0}a\u{85}"), "\u{a0}a\u{85}");
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2817-2839 (cited at
// `string_to_bool`).
#[test]
fn ncbis_booleans() {
    for (text, value) in [
        ("1", Some(true)),
        ("TRUE", Some(true)),
        ("t", Some(true)),
        ("On", Some(true)),
        ("0", Some(false)),
        ("No", Some(false)),
        ("off", Some(false)),
        ("", None),
        ("2", None),
    ] {
        assert_eq!(string_to_bool(text), value, "{text:?}");
    }
}

// NCBI BLAST+ 2.17.0 (line NP M3d, `M/m3d/ws/run2.sh`): the parameters (base 0) and units
// (base 16) of ascii unit counts, their messages, lines and positions.
#[test]
fn base_conversions_are_ncbis() {
    assert_eq!(string_to_uint_base(b"0xad", 0).unwrap(), 173);
    assert_eq!(string_to_uint_base(b"0XAD", 0).unwrap(), 173);
    assert_eq!(string_to_uint_base(b"0255", 0).unwrap(), 173);
    assert_eq!(string_to_uint_base(b"+173", 0).unwrap(), 173);
    assert_eq!(string_to_uint_base(b"0x", 0).unwrap(), 0);
    assert_eq!(string_to_uint_base(b"0", 0).unwrap(), 0);
    assert_eq!(string_to_uint_base(b"3fd00", 16).unwrap(), 0x3fd00);
    assert_eq!(string_to_uint_base(b"0x3FD00", 16).unwrap(), 0x3fd00);
    assert_eq!(string_to_uint_base(b"173", 10).unwrap(), 173);
    let failure = |e: NcbiException| {
        let what = String::from_utf8(e.what()).unwrap();
        let pos = what.rsplit_once(crate::wording::POSITION).unwrap().1;
        (
            String::from_utf8(e.get_msg().to_vec()).unwrap(),
            e.site().line,
            pos.trim_end_matches(")\n").to_string(),
        )
    };
    let cases: [(&[u8], u32, &str, u32, &str); 10] = [
        (b"abc", 0, "Cannot convert string 'abc' to Uint8", 973, "0"),
        (b"-5", 0, "Cannot convert string '-5' to Uint8", 973, "0"),
        (b"08", 0, "Cannot convert string '08' to Uint8", 980, "1"),
        (
            b"173 ",
            0,
            "Cannot convert string '173 ' to Uint8",
            980,
            "3",
        ),
        (
            b"173\r",
            0,
            "Cannot convert string '173\\r' to Uint8",
            980,
            "3",
        ),
        (
            b"17\x003",
            0,
            "Cannot convert string '17\\0003' to Uint8",
            980,
            "2",
        ),
        (
            b"99999999999999999999",
            0,
            "Cannot convert string '99999999999999999999' to Uint8, overflow",
            964,
            "19",
        ),
        (
            b"0x10000000000000000",
            0,
            "Cannot convert string '0x10000000000000000' to Uint8, overflow",
            964,
            "18",
        ),
        (b"1g", 16, "Cannot convert string '1g' to Uint8", 980, "1"),
        (
            b"11111111111111111",
            16,
            "Cannot convert string '11111111111111111' to Uint8, overflow",
            964,
            "16",
        ),
    ];
    for (input, base, message, line, pos) in cases {
        assert_eq!(
            failure(string_to_uint8_base(input, base).unwrap_err()),
            (message.to_string(), line, pos.to_string()),
            "{input:?}"
        );
    }
    for (input, base) in [(&b"4294967296"[..], 0), (b"1ffffffff", 16)] {
        assert_eq!(
            failure(string_to_uint_base(input, base).unwrap_err()),
            (
                format!(
                    "Cannot convert string '{}' to unsigned int, overflow",
                    String::from_utf8_lossy(input)
                ),
                652,
                "0".to_string()
            )
        );
    }
}

#[test]
fn c_atoi_is_glibcs() {
    assert_eq!(c_atoi(b" 42xyz"), 42);
    assert_eq!(c_atoi(b"\t\r\x0b+7"), 7);
    assert_eq!(c_atoi(b"-3"), -3);
    assert_eq!(c_atoi(b""), 0);
    assert_eq!(c_atoi(b"x1"), 0);
    assert_eq!(c_atoi(b"- 1"), 0);
    assert_eq!(c_atoi(b"9\x001"), 9);
    // `long` truncated to `int`: 99999999999 = 0x174876E7FF.
    assert_eq!(c_atoi(b"99999999999"), 0x4876_E7FF);
    assert_eq!(c_atoi(b"4294967297"), 1);
    // Beyond `long`: LONG_MAX (-1 as `int`) and LONG_MIN (0).
    assert_eq!(c_atoi(b"99999999999999999999"), -1);
    assert_eq!(c_atoi(b"-99999999999999999999"), 0);
    assert_eq!(c_atoi(b"-9223372036854775808"), 0);
    assert_eq!(c_atoi(b"-9223372036854775807"), 1);
}

/// `NStr::Replace` as dev had it before `find` (a slice compare at every byte), the reference
/// of the tests.
fn replace_by_windows(src: &[u8], search: &[u8], replacement: &[u8]) -> Vec<u8> {
    let mut dst = Vec::new();
    let mut pos = 0;
    while let Some(found) = src[pos..]
        .windows(search.len())
        .position(|window| window == search)
    {
        dst.extend_from_slice(&src[pos..pos + found]);
        dst.extend_from_slice(replacement);
        pos += found + search.len();
    }
    dst.extend_from_slice(&src[pos..]);
    dst
}

#[test]
fn replace_is_the_window_search() {
    let texts = search_texts();
    for search in search_needles() {
        for replacement in [&b""[..], b"X", b"   ", b"a", b"aa", b"\\t"] {
            for text in &texts {
                let expected = replace_by_windows(text, &search, replacement);
                let got = replace(text, &search, replacement);
                assert_eq!(*got, expected, "{text:?} {search:?} {replacement:?}");
                let unchanged = search == replacement || find(text, &search).is_none();
                assert_eq!(
                    matches!(got, Cow::Borrowed(_)),
                    unchanged,
                    "{text:?} {search:?} {replacement:?}"
                );
                let mut in_place = text.clone();
                replace_in_place(&mut in_place, &search, replacement);
                assert_eq!(in_place, expected, "{text:?} {search:?} {replacement:?}");
            }
        }
    }
}

/// Occurrences from the left, at the ends, overlapping, and in a replacement (never read again).
#[test]
fn replace_takes_occurrences_from_the_left() {
    let r = |src: &[u8], search: &[u8], replacement: &[u8]| {
        replace(src, search, replacement).into_owned()
    };
    assert_eq!(r(b"", b"\\t", b"   "), b"");
    assert_eq!(r(b"\\t", b"\\t", b"   "), b"   ");
    assert_eq!(r(b"\\ta\\t", b"\\t", b"   "), b"   a   ");
    assert_eq!(r(b"a\\\\tb", b"\\t", b"   "), b"a\\   b");
    assert_eq!(r(b"aaa", b"aa", b"X"), b"Xa");
    assert_eq!(r(b"aaaa", b"aa", b"a"), b"aa");
    assert_eq!(r(b"abab", b"ab", b"abab"), b"abababab");
    assert_eq!(r(b"x. [. [", b". [", b" ["), b"x [ [");
    assert_eq!(r(b"a. . [b", b". [", b" ["), b"a.  [b");
    assert_eq!(r(b"b\x01\x01\x01", b"\x01", b" "), b"b   ");
    assert_eq!(r(b"a", b"ab", b"X"), b"a");
}
