use super::*;

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

// NCBI's rpsbproc printed `Cannot convert string 'abc' to double (m_Pos = 0)` for a bit score
// of `bitscore_specific.txt` (PS16a probe 4), read `0x10` as 16 and `0X1p3` as 8 and refused
// `0x` at position 1 (PS16a review oracle A4-1 cases `bs12`, `bs14`, `bs13`).
#[test]
fn values_and_errors_are_ncbis() {
    assert_eq!(string_to_double(b"55.5"), Ok(55.5));
    assert_eq!(string_to_double(b"-1e3"), Ok(-1000.0));
    assert_eq!(string_to_double(b"0x10"), Ok(16.0));
    assert_eq!(string_to_double(b"0X1p3"), Ok(8.0));
    assert_eq!(string_to_double(b"12\x00x"), Ok(12.0));
    let cases: [(&[u8], &str, u32, &str); 7] = [
        (b"abc", "Cannot convert string 'abc' to double", 1317, "0"),
        (b"", "Cannot convert string '' to double", 1317, "0"),
        (b"-x", "Cannot convert string '-x' to double", 1343, "0"),
        (
            b"12.5\r",
            "Cannot convert string '12.5\\r' to double",
            1374,
            "4",
        ),
        (b"1e", "Cannot convert string '1e' to double", 1374, "1"),
        (b"0x", "Cannot convert string '0x' to double", 1374, "1"),
        // the message quotes the C string, up to its first 0 byte
        (
            b"1x\x00y",
            "Cannot convert string '1x' to double",
            1374,
            "1",
        ),
    ];
    for (text, message, line, pos) in cases {
        assert_eq!(
            failure(string_to_double(text).unwrap_err()),
            (message.to_string(), line, pos.to_string()),
            "{text:?}"
        );
    }
    assert_eq!(
        string_to_double(b"abc").unwrap_err().site(),
        ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 1317,
            kind: "CStringException::eConvert",
            function: Some("ncbi::s_StringToDouble"),
        }
    );
}

// `LOSAT/tests/unit/string_to_double_posix_or_local_vectors.tsv`: NCBI's own
// `NStr::StringToDouble(text, fDecimalPosixOrLocal)` (its corelib linked to a probe), for
// the texts that end in an exponent mark (`StringToDoublePosix` reads them past `strtod`,
// in x87 extended precision), the edges of glibc's `strtod` that `CArg_Double` reads, and
// the errors.
#[test]
fn posix_or_local_values_are_ncbis() {
    let table = include_str!("../../tests/unit/string_to_double_posix_or_local_vectors.tsv");
    let mut rows = 0;
    for row in table.lines().filter(|row| !row.starts_with('#')) {
        let (text, bits) = row.rsplit_once('\t').unwrap();
        match (string_to_double_posix_or_local(text.as_bytes()), bits) {
            (Err(_), "ERR") => {}
            (Ok(value), bits) if bits != "ERR" => {
                assert_eq!(format!("{:016x}", value.to_bits()), bits, "{text:?}")
            }
            (other, bits) => panic!("{text:?}: {other:?}, NCBI {bits}"),
        }
        rows += 1;
    }
    assert_eq!(rows, 422);
}

/// The four conversions of a text as NCBI's probe prints them: `StringToDouble(text,
/// fDecimalPosixOrLocal)`, `StringToDouble(text)`, `StringToDoublePosix(text, &end, 0)` and
/// with `fDecimalPosixFinite`, each the bits of the double or `ERR`, the last two with
/// `:end`.
fn four_conversions(text: &[u8]) -> [String; 4] {
    let bits = |value: f64| format!("{:016x}", value.to_bits());
    let posix = |finite| {
        let read = string_to_double_posix(text, finite);
        format!("{}:{}", bits(read.value), read.end)
    };
    [
        string_to_double_posix_or_local(text).map_or("ERR".into(), bits),
        string_to_double(text).map_or("ERR".into(), bits),
        posix(false),
        posix(true),
    ]
}

// `LOSAT/tests/unit/string_to_double_posix_vectors.tsv`: NCBI's own `StringToDouble` with
// both flag sets and `StringToDoublePosix` with both (its corelib linked to a probe,
// `swebok-improve-20261105/sw8/nfix/probe/x87`): mantissas of more than 18 digits and more
// than 15 places (`third`, `mul2`), exponents at the ends of NCBI's ranges, NaN (NCBI's
// `HUGE_VAL/HUGE_VAL` has the sign bit), and random texts. The probe compared 900,184 texts
// (these, 900,000 random ones and mantissas of up to 5,000 digits): all equal but the one
// whose `int` arithmetic overflows in NCBI (below).
#[test]
fn posix_values_and_ends_are_ncbis() {
    let table = include_str!("../../tests/unit/string_to_double_posix_vectors.tsv");
    let mut rows = 0;
    for row in table.lines().filter(|row| !row.starts_with('#')) {
        let mut fields = row.split('\t');
        let text = fields.next().unwrap();
        let ncbi: Vec<&str> = fields.collect();
        assert_eq!(
            four_conversions(text.as_bytes()),
            ncbi.as_slice(),
            "{text:?}"
        );
        rows += 1;
    }
    assert_eq!(rows, 590);
    // Texts of their own: white space first (`StringToDoublePosix` skips it, `StringToDouble`
    // refuses it), the empty text, and long mantissas, whose x87 value passes the largest
    // 80-bit value from 4,933 digits on (the infinity; scaled by at least 10^-614 it stays
    // above the largest double, so NCBI's result is the infinity or, below 10^-614, zero).
    let repeat = |byte: &str, count: usize| byte.repeat(count);
    let cases: [(String, [&str; 4]); 12] = [
        (
            "\x0b1.5".into(),
            ["ERR", "ERR", "3ff8000000000000:4", "3ff8000000000000:4"],
        ),
        (
            "\x0c-2.5e".into(),
            ["ERR", "ERR", "c004000000000000:6", "c004000000000000:6"],
        ),
        (
            "".into(),
            ["ERR", "ERR", "0000000000000000:0", "0000000000000000:0"],
        ),
        (
            format!("1{}", repeat("0", 4931)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "7ff0000000000000:4932",
                "7ff0000000000000:4932",
            ],
        ),
        (
            format!("1{}", repeat("0", 4933)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "7ff0000000000000:4934",
                "7ff0000000000000:4934",
            ],
        ),
        (
            format!("1{}e", repeat("0", 5000)),
            [
                "7ff0000000000000",
                "ERR",
                "7ff0000000000000:5002",
                "7ff0000000000000:5002",
            ],
        ),
        (
            format!("{}e1", repeat("9", 5000)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "7ff0000000000000:5002",
                "7fefffffffffffff:5002",
            ],
        ),
        (
            format!("{}e-1", repeat("9", 5000)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "7ff0000000000000:5003",
                "7ff0000000000000:5003",
            ],
        ),
        (
            format!("{}.{}e", repeat("9", 5000), repeat("9", 700)),
            [
                "0000000000000000",
                "ERR",
                "0000000000000000:5702",
                "0010000000000000:5702",
            ],
        ),
        (
            format!("{}.{}e", repeat("9", 5000), repeat("9", 600)),
            [
                "7ff0000000000000",
                "ERR",
                "7ff0000000000000:5602",
                "7ff0000000000000:5602",
            ],
        ),
        (
            format!("{}e-600", repeat("5", 4950)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "7ff0000000000000:4955",
                "7ff0000000000000:4955",
            ],
        ),
        (
            format!("{}e-616", repeat("5", 4950)),
            [
                "7ff0000000000000",
                "7ff0000000000000",
                "0000000000000000:4955",
                "0010000000000000:4955",
            ],
        ),
    ];
    for (text, ncbi) in cases {
        assert_eq!(
            four_conversions(text.as_bytes()),
            ncbi,
            "{}",
            &text[..text.len().min(20)]
        );
    }
    // NCBI's `exponent - expvalue` (-2 - 2147483647) overflows `int`: undefined behaviour,
    // which NCBI's corelib wraps to 2147483647 (the infinity; DBL_MAX with
    // `fDecimalPosixFinite`). LOSAT's widened exponent is -2147483649, below -614: 0, and
    // DBL_MIN with `fDecimalPosixFinite` (the Owner's rule of 2026-10-10 for NCBI's signed
    // overflows). `strtod` reads the whole text, so `CArg_Double` takes its 0 as NCBI does.
    assert_eq!(
        four_conversions(b"1.55e-2147483647"),
        [
            "0000000000000000",
            "0000000000000000",
            "0000000000000000:16",
            "0010000000000000:16"
        ]
    );
    assert_eq!(
        four_conversions(b"-1.55e-2147483647"),
        [
            "8000000000000000",
            "8000000000000000",
            "8000000000000000:17",
            "8010000000000000:17"
        ]
    );
    // `1e2147483648` and longer exponents: NCBI's guard keeps 214748364 where `expvalue * 10`
    // wraps, LOSAT the whole (saturated) exponent; both are beyond 616 (the infinity).
    assert_eq!(
        four_conversions(b"1e99999999999999999999999999"),
        [
            "7ff0000000000000",
            "7ff0000000000000",
            "7ff0000000000000:28",
            "7fefffffffffffff:28"
        ]
    );
}

// `LOSAT/tests/unit/strtod_hex_vectors.tsv`: glibc 2.39's `strtod` of hexadecimal texts
// (edge cases of rounding, subnormals, overflow and the ends of the form, then random ones),
// made with `ctypes` by `psiblast-20261025/PR/tmp/impl16a/strtodhex.py`.
#[test]
fn hex_numbers_read_as_glibcs_strtod() {
    let table = include_str!("../../tests/unit/strtod_hex_vectors.tsv");
    let mut rows = 0;
    for row in table.lines().filter(|row| !row.starts_with('#')) {
        let mut fields = row.split('\t');
        let (text, bits, used) = (
            fields.next().unwrap(),
            fields.next().unwrap(),
            fields.next().unwrap(),
        );
        let (value, end) = strtod_prefix(text.as_bytes());
        assert_eq!(
            (format!("{:016x}", value.to_bits()), end.to_string()),
            (bits.to_string(), used.to_string()),
            "{text:?}"
        );
        rows += 1;
    }
    assert_eq!(rows, 360);
}

// glibc's `strtod` reads these to their end and stops before the end of these (the SEG
// locut and hicut of the argument reader, `blastinput::app::read_seg_option`, test
// `string_to_double` with it, which also needs a digit, a point or a sign first).
#[test]
fn strtod_reads_the_whole_of_these_only() {
    let whole = |text: &str| {
        let (_, used) = strtod_prefix(text.as_bytes());
        used != 0 && used == text.len()
    };
    for text in [
        "1",
        "+1",
        "-.5",
        "5.",
        "1e5",
        "1E-5",
        "0x10",
        "0X1p3",
        "0x.8p-1",
        "0x1.",
        "inf",
        "+INFINITY",
        "-nan",
        "nan(x_1)",
        "1e400",
    ] {
        assert!(whole(text), "{text}");
    }
    for text in [
        "1e", "1e+", "2.5x", "2,2", ".", "+.", ".e1", "0x", "0x.", "0x1p", "1_0", "--5", "2.5.5",
        "nan(", "nan(-)", "infx", "",
    ] {
        assert!(!whole(text), "{text}");
    }
}
