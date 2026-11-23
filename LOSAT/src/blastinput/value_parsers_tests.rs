use super::*;

#[test]
fn integers_are_read_as_ncbi_reads_them() {
    for (value, expected) in [
        ("7", 7),
        ("+7", 7),
        ("-07", -7),
        ("0x1F", 31),
        ("0XaB", 171),
        ("0x", 0),
        ("2147483647", i32::MAX),
        ("-2147483648", i32::MIN),
    ] {
        assert_eq!(ncbi_integer(value), Ok(expected), "{value}");
    }
    for value in [
        "",
        " 7",
        "7 ",
        "+",
        "-0x5",
        "0x+5",
        "0x5g",
        "1_000",
        "2147483648",
        "0x80000000",
        "00x5",
        "7.0",
    ] {
        assert!(ncbi_integer(value).is_err(), "{value:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129 (`s_StringToInt8`)
// An Int8 argument (`-dbsize`) reads what an integer argument reads, in 64 bits.
#[test]
fn int8_arguments_are_read_as_ncbi_reads_them() {
    for (value, expected) in [
        ("7", 7),
        ("+7", 7),
        ("-07", -7),
        ("0", 0),
        ("0x1F", 31),
        ("0XaB", 171),
        ("0x", 0),
        ("3000000000", 3_000_000_000),
        ("-3000000000", -3_000_000_000),
        ("9223372036854775807", i64::MAX),
        ("-9223372036854775808", i64::MIN),
        ("0x7FFFFFFFFFFFFFFF", i64::MAX),
    ] {
        assert_eq!(ncbi_int8(value), Ok(expected), "{value}");
    }
    for value in [
        "",
        " 7",
        "7 ",
        "+",
        "-0x5",
        "0x5g",
        "1_000",
        "9223372036854775808",
        "-9223372036854775809",
        "0x8000000000000000",
        "7.0",
        "1e3",
        "abc",
    ] {
        assert!(ncbi_int8(value).is_err(), "{value:?}");
    }
}

// LOSAT's messages for the integers NCBI cannot convert (`s_StringToInt8`, then the range of
// `int` for `CArg_Integer`): a byte after `0x` that is not a hexadecimal digit is tested before
// the base 16 overflow.
#[test]
fn integer_arguments_name_what_failed() {
    for value in ["0x10000000000000000g", "9223372036854775808", "abc", "0x5g"] {
        assert_eq!(
            ncbi_int8(value),
            Err("expected an integer".to_string()),
            "{value}"
        );
    }
    assert_eq!(
        ncbi_int8("0x8000000000000000"),
        Err("Integer value is out of range".to_string())
    );
    for value in [
        "2147483648",
        "0x80000000",
        "-2147483649",
        "0x8000000000000000",
    ] {
        assert_eq!(
            ncbi_integer(value),
            Err("Integer value is out of range".to_string()),
            "{value}"
        );
    }
    assert_eq!(
        ncbi_integer("9223372036854775808"),
        Err("expected an integer".to_string())
    );
}

// `CArgAllowGeneticCodeInteger::Verify` (`NStr::StringToInt` of `value`, blast_args.cpp:998-1007):
// the number of a code of 0 to 255 written without a minus sign.
#[test]
fn genetic_codes_are_read_as_ncbi_reads_them() {
    for (value, expected) in [("1", 1), ("+5", 5), ("005", 5), ("33", 33), ("255", 255)] {
        assert_eq!(genetic_code_number(value), Some(expected), "{value}");
    }
    for value in [
        "-0",
        "-5",
        "256",
        "2147483648",
        "",
        " 5",
        "5 ",
        "5.0",
        "0x",
        "0x5",
    ] {
        assert_eq!(genetic_code_number(value), None, "{value:?}");
    }
}

// NCBI 2.17.0's `CArg_Double`: the bits of `NStr::StringToDouble` of `value`
// with `fDecimalPosixOrLocal` from NCBI's corelib (swebok-improve-20261105/sw8/nfix/probe/x87,
// `s2dp`), and `blastn -evalue <value>` (exit 0, or 1 with the usage for the errors).
#[test]
fn reals_are_read_as_ncbis_carg_double() {
    for (value, bits) in [
        ("1e-5", 0x3ee4_f8b5_88e3_68f1),
        (".5", 0x3fe0_0000_0000_0000),
        ("+2", 0x4000_0000_0000_0000),
        ("-0", 0x8000_0000_0000_0000),
        // glibc's hexadecimal numbers
        ("0x10", 0x4030_0000_0000_0000),
        ("0x1p-3", 0x3fc0_0000_0000_0000),
        ("-0x.8", 0xbfe0_0000_0000_0000),
        // `StringToDoublePosix` reads a final exponent mark, in x87 extended precision
        ("1e", 0x3ff0_0000_0000_0000),
        ("2.5E", 0x4004_0000_0000_0000),
        ("-3e", 0xc008_0000_0000_0000),
        ("12345678901234567890e", 0x43e5_6a95_319d_63e1),
        ("0.0000000000000001e", 0x3c9c_d2b2_97d8_89bc),
        // the infinities, and NaN with glibc's sign and payload
        ("1e400", 0x7ff0_0000_0000_0000),
        ("+Infinity", 0x7ff0_0000_0000_0000),
        // `StringToDoublePosix` reads an `inf` after a point, which `strtod` does not
        (".inf", 0x7ff0_0000_0000_0000),
        ("-.Infinity", 0xfff0_0000_0000_0000),
        ("-nan", 0xfff8_0000_0000_0000),
        ("-nan(1)", 0xfff8_0000_0000_0001),
        ("+nan(5)", 0x7ff8_0000_0000_0005),
        ("-NaN()", 0xfff8_0000_0000_0000),
        ("+nan(x_9)", 0x7ff8_0000_0000_0000),
    ] {
        let read = ncbi_double(value).unwrap_or_else(|e| panic!("{value}: {e}"));
        assert_eq!(read.to_bits(), bits, "{value}");
    }
    for value in [
        "+nan(",
        "+nan(1",
        "+nan(-)",
        "+nan(1)x",
        "+-nan(1)",
        "nan(1)",
        "nan(5)",
        ",5",
        ",",
        "1x",
        "1e+",
        "1ex",
        "--5",
        "2,5",
        "inf",
        "nan",
        " 1",
        "",
        "e5",
        "0x",
        "1.2.3",
        // A character of several bytes across the fourth byte (S11 audit F-1).
        "-xx€",
        "-in\u{fffd}",
        "+\u{fffd}\u{fffd}",
        "-1éé",
    ] {
        assert_eq!(
            ncbi_double(value),
            Err("expected a number".into()),
            "{value:?}"
        );
    }
}

// NCBI 2.17.0: `-threshold`'s constraint reads the value with the default flags (`blastp
// -threshold 1e` and `12345678901234567890e`: "Illegal value, expected >=0", exit 1; `0x10`
// searches), `-perc_identity`'s with `fDecimalPosixOrLocal` (`1e` and `0x10` search, `0x1p7`
// is illegal): NCBI's argument reader checks them (`cli::reader`), the conversion reads the number.
#[test]
fn each_constraint_reads_the_value_with_its_own_conversion() {
    let read = |program: &str, name: &str, value: &str| {
        let argv = [
            "losat", program, "-query", "q", "-subject", "s", name, value,
        ];
        crate::cli::try_parse_from(argv)
    };
    let illegal = |program: &str, name: &str, value: &str| match read(program, name, value) {
        Err(crate::cli::ArgvError::Report(stderr)) => {
            String::from_utf8_lossy(&stderr).contains("Illegal value")
        }
        _ => false,
    };
    for value in ["0x10", "11", "1e400"] {
        assert!(read("blastp", "-threshold", value).is_ok(), "{value}");
    }
    for value in [
        "1e",
        "2.5E",
        "12345678901234567890e",
        "-nan",
        "+nan(5)",
        "-1",
        "-0x1",
    ] {
        assert!(illegal("blastp", "-threshold", value), "{value}");
    }
    assert!(!illegal("blastp", "-threshold", "1x"));
    assert!(read("blastp", "-threshold", "1x").is_err());
    assert_eq!(ncbi_double("0x10"), Ok(16.0));
    for value in ["1e", "0x10"] {
        assert!(read("blastn", "-perc_identity", value).is_ok(), "{value}");
    }
    assert!(read("blastp", "-qcov_hsp_perc", "99.5E").is_ok());
    for value in ["0x1p7", "1e400", "-nan", "-0.5e"] {
        assert!(illegal("blastn", "-perc_identity", value), "{value}");
    }
}
