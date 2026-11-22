//! NCBI's argument reader (`args`) on the descriptions of a small program: the conversions of
//! `eInteger` keys, the defaults, the flags and the streams of the file keys. The texts are
//! NCBI BLAST+ 2.17.0's (dustmasker probe of line RG, item K1: `LOSAT/tests/dustmasker/ncbi`,
//! and probe B of line BDC for blastdbcmd).

use super::*;
use crate::ncbi_app::no_dependencies;

/// A program with a key of each kind.
static TOOL: Tool = Tool {
    program: "t",
    refusal_name: "T",
    help: "",
    descs: &[
        key("n", Kind::Integer, Some("7"), None),
        key("m", Kind::Integer, None, Some(Constraint::AtLeast(1))),
        key("b", Kind::Integer, None, Some(Constraint::Between(0, 4))),
        key("in", Kind::InputFile, Some("-"), None),
        key("opt", Kind::Text, None, None),
        flag("f"),
    ],
    dependencies: no_dependencies,
    post_prefix: None,
    components: &[],
};

fn parse_words(line: &[&str]) -> Result<Args, ArgsError> {
    let words: Vec<String> = line.iter().map(|word| word.to_string()).collect();
    parse(&TOOL, &words)
}

fn messages(line: &[&str]) -> (&'static str, Vec<String>) {
    match parse_words(line) {
        Err(ArgsError::Invalid(ArgException { code, messages })) => (code, messages),
        other => panic!("{line:?}: {other:?}"),
    }
}

/// The texts as LOSAT writes them: the expected texts are NCBI's (`(m_Pos = <n>)` and the like),
/// in LOSAT's words (`crate::wording::losat_text`).
fn strings(texts: &[&str]) -> Vec<String> {
    texts
        .iter()
        .map(|text| crate::wording::losat_text(text).into_owned())
        .collect()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:115-129, 351-360, 383-390 (cited at
// `integer_value`).
#[test]
fn integer_keys_convert_as_ncbi_s_reader() {
    let value = |text: &str| parse_words(&["-n", text]).unwrap().integer("n");
    assert_eq!(value("+5"), Some(5));
    assert_eq!(value("-1"), Some(-1));
    assert_eq!(value("-2147483648"), Some(i32::MIN));
    assert_eq!(value("2147483647"), Some(i32::MAX));
    // Not octal: a leading zero is a decimal digit.
    assert_eq!(value("010"), Some(10));
    // Base 16 after base 10 fails, for a value that starts with `0x` (`big_win_0x10`).
    assert_eq!(value("0x10"), Some(16));
    assert_eq!(value("0XfF"), Some(255));
    let cases: &[(&str, &[&str])] = &[
        (
            "x",
            &[
                "Cannot convert string 'x' to Int8 (m_Pos = 0)",
                "Argument \"n\". Argument cannot be converted:  `x'",
            ],
        ),
        (
            "",
            &[
                "Cannot convert string '' to Int8 (m_Pos = 0)",
                "Argument \"n\". Argument cannot be converted",
            ],
        ),
        (
            "2.5",
            &[
                "Cannot convert string '2.5' to Int8 (m_Pos = 1)",
                "Argument \"n\". Argument cannot be converted:  `2.5'",
            ],
        ),
        (
            "5 ",
            &[
                "Cannot convert string '5 ' to Int8 (m_Pos = 1)",
                "Argument \"n\". Argument cannot be converted:  `5 '",
            ],
        ),
        (
            " 5",
            &[
                "Cannot convert string ' 5' to Int8 (m_Pos = 0)",
                "Argument \"n\". Argument cannot be converted:  ` 5'",
            ],
        ),
        (
            "9223372036854775808",
            &[
                "Cannot convert string '9223372036854775808' to Int8, overflow (m_Pos = 18)",
                "Argument \"n\". Argument cannot be converted:  `9223372036854775808'",
            ],
        ),
        // The exception of the conversion in base 16 is the one reported.
        (
            "0x1g",
            &[
                "Cannot convert string '0x1g' to Int8 (m_Pos = 3)",
                "Argument \"n\". Argument cannot be converted:  `0x1g'",
            ],
        ),
        (
            "2147483648",
            &["Argument \"n\". Integer value is out of range:  `2147483648'"],
        ),
        (
            "-2147483649",
            &["Argument \"n\". Integer value is out of range:  `-2147483649'"],
        ),
    ];
    for &(text, expected) in cases {
        assert_eq!(
            messages(&["-n", text]),
            ("eConvert", strings(expected)),
            "{text:?}"
        );
    }
    // A key whose constraint compares numbers refuses a hexadecimal value.
    assert_eq!(
        parse_words(&["-m", "0x10"]).unwrap_err(),
        ArgsError::Refused(
            "the hexadecimal value '0x10' of -m is not supported by LOSAT's T".to_string()
        )
    );
    assert_eq!(
        parse_words(&["-b", "0x1"]).unwrap_err(),
        ArgsError::Refused(
            "the hexadecimal value '0x1' of -b is not supported by LOSAT's T".to_string()
        )
    );
    assert_eq!(
        messages(&["-m", "0"]),
        (
            "eConstraint",
            strings(&["Argument \"m\". Illegal value, expected >=1:  `0'"])
        )
    );
}

#[test]
fn defaults_and_flags_are_ncbi_s_values() {
    let args = parse_words(&[]).unwrap();
    assert_eq!(args.integer("n"), Some(7));
    assert_eq!(args.value("n"), Some("7"));
    assert_eq!(args.integer("m"), None);
    assert_eq!(args.value("in"), Some("-"));
    assert!(!args.flag("f") && !args.given("in"));
    let args = parse_words(&["-f", "-m=3", "-in", "x.fa"]).unwrap();
    assert!(args.flag("f"));
    assert_eq!(args.integer("m"), Some(3));
    assert_eq!(args.value("in"), Some("x.fa"));
    assert!(args.given("in") && !args.given("n"));
    // A file key without a value has no stream (`CArg_NoValue`).
    assert_eq!(
        args.output_file("opt", &mut Vec::new()).err().unwrap(),
        ArgException::new(
            "eNoValue",
            "Argument \"opt\". The argument has no value".to_string()
        )
    );
}

#[test]
fn every_reader_knows_the_help_flags_and_the_file_arguments() {
    assert_eq!(parse_words(&["-h"]).unwrap_err(), ArgsError::Usage);
    assert_eq!(parse_words(&["-help-full"]).unwrap_err(), ArgsError::Help);
    assert_eq!(parse_words(&["-xmlhelp"]).unwrap_err(), ArgsError::XmlHelp);
    let args = parse_words(&["-logfile", "log.txt", "-conffile", "t.ini"]).unwrap();
    assert_eq!(args.value("logfile"), Some("log.txt"));
    assert_eq!(
        messages(&["-in", "a", "-in", "b"]),
        (
            "eSynopsis",
            strings(&["Argument with this name is defined already: in"])
        )
    );
    assert_eq!(
        messages(&["-foo"]),
        ("eInvalidArg", strings(&["Unknown argument: \"foo\""]))
    );
    assert_eq!(
        messages(&["-f", "true"]),
        (
            "eSynopsis",
            strings(&["Too many positional arguments (1), the offending value: true"])
        )
    );
}

// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:161-184
// (quoted at `Constraint::usage`). NCBI BLAST+ 2.17.0 blast_formatter `-sorthits 5` and
// `-sorthsps 5` (line RG item F3, probes P2.055 and P2.073): `Illegal value, expected (>=0 and
// =<4):  `5'`.
#[test]
fn a_key_between_two_numbers_takes_both_ends() {
    for (text, number) in [("0", 0), ("4", 4), ("+2", 2), ("004", 4)] {
        assert_eq!(
            parse_words(&["-b", text]).unwrap().integer("b"),
            Some(number),
            "{text}"
        );
    }
    for text in ["5", "-1", "2147483647"] {
        assert_eq!(
            messages(&["-b", text]),
            (
                "eConstraint",
                strings(&[&format!(
                    "Argument \"b\". Illegal value, expected (>=0 and =<4):  `{text}'"
                )])
            ),
            "{text}"
        );
    }
}
