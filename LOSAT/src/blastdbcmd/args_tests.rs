use super::*;
use crate::ncbi_app::args::{ArgException, Args, ArgsError};

fn words(line: &[&str]) -> Vec<std::ffi::OsString> {
    line.iter().map(std::ffi::OsString::from).collect()
}

/// The arguments of `line` as NCBI's reader reads blastdbcmd's.
fn parse(line: &[std::ffi::OsString]) -> Result<Args, ArgsError> {
    crate::ncbi_app::args::parse(&super::super::TOOL, line)
}

fn messages(line: &[&str]) -> (&'static str, Vec<String>) {
    match parse(&words(line)) {
        Err(ArgsError::Invalid(ArgException { code, messages })) => {
            let texts = messages
                .into_iter()
                .map(|message| String::from_utf8(message).unwrap());
            (code, texts.collect())
        }
        other => panic!("{line:?}: {other:?}"),
    }
}

#[test]
fn defaults_fill_the_keys_not_given() {
    let args = parse(&words(&["-db", "db/p5", "-info"])).unwrap();
    assert_eq!(args.value("db"), Some("db/p5"));
    assert_eq!(args.value("dbtype"), Some("guess"));
    assert_eq!(args.value("out"), Some("-"));
    assert!(args.given("info") && !args.given("exact_length"));
    // `-info` excludes `-outfmt` and `-strand`: no default for them.
    assert_eq!(args.value("outfmt"), None);
    assert_eq!(args.value("strand"), None);
    let args = parse(&words(&["-db=db/p5", "-info"])).unwrap();
    assert_eq!(args.value("db"), Some("db/p5"));
    assert_eq!(parse(&words(&["-info"])).unwrap().value("db"), Some("nr"));
}

// The texts of NCBI blastdbcmd 2.17.0 (probe B of line BDC, `manifest_b.tsv`; I33, E60, E62).
#[test]
fn refusals_have_ncbi_s_texts() {
    let cases: &[(&[&str], &str, &str)] = &[
        (
            &["-db", "db/p5", "-info", "-foo"],
            "eInvalidArg",
            "Unknown argument: \"foo\"",
        ),
        (&["-db"], "eNoArg", "Argument \"-db\". Value is missing"),
        (
            &["-db", "db/p5", "-info", "-dbtype", "protein"],
            "eConstraint",
            "Argument \"dbtype\". Illegal value, expected `guess', `nucl', `prot':  `protein'",
        ),
        (
            &["-db", "db/p5", "-info", "-dbtype", ""],
            "eConstraint",
            "Argument \"dbtype\". Illegal value, expected `guess', `nucl', `prot'",
        ),
        (
            &["-db", "db/p5", "-info", "-db", "db/n5"],
            "eSynopsis",
            "Argument with this name is defined already: db",
        ),
        (
            &["-db", "db/p5", "-exact_length"],
            "eConstraint",
            "Argument \"info\". Must be specified, as it is required by argument:  `exact_length'",
        ),
        (
            &["-db", "db/p5", "-info", "extra"],
            "eSynopsis",
            "Too many positional arguments (1), the offending value: extra",
        ),
        (
            &["-db", "db/p5", "-info", "--", "-x"],
            "eSynopsis",
            "Too many positional arguments (1), the offending value: -x",
        ),
        (
            &["-db", "db/p5", "-info", "-"],
            "eSynopsis",
            "Too many positional arguments (1), the offending value: -",
        ),
        (
            &["-db", "db/p5", "-info=1"],
            "eInvalidArg",
            "Invalid argument: -info=1",
        ),
        (
            &["-db", "db/p5", "-info", "-entry", "XP_000001.1"],
            "eConstraint",
            "Argument \"info\". Incompatible with argument:  `entry'",
        ),
        (
            &["-db", "x", "-entry", "a", "-entry_batch", "b"],
            "eConstraint",
            "Argument \"entry_batch\". Incompatible with argument:  `entry'",
        ),
        (
            &["-db", "x", "-entry", "a", "-get_dups", "-target_only"],
            "eConstraint",
            "Argument \"target_only\". Incompatible with argument:  `get_dups'",
        ),
        (
            &["-db", "db/p5", "-info", "-outfmt", "%f"],
            "eConstraint",
            "Argument \"outfmt\". Incompatible with argument:  `info'",
        ),
        (
            &["-db", "db/p5", "-info", "-list", "db"],
            "eConstraint",
            "Argument \"list\". Incompatible with argument:  `db'",
        ),
        (
            &["-info", "-show_blastdb_search_path"],
            "eConstraint",
            "Argument \"show_blastdb_search_path\". Incompatible with argument:  `info'",
        ),
        (
            &["-db", "db/p5", "-info", "-metadata_output_prefix", "x"],
            "eConstraint",
            "Argument \"metadata_output_prefix\". Incompatible with argument:  `info'",
        ),
        (
            &["-db", "db/p5", "-info", "-line_length", "0"],
            "eConstraint",
            "Argument \"line_length\". Illegal value, expected >=1:  `0'",
        ),
        (
            &["-db", "db/p5", "-info", "-ipg", "-1"],
            "eConstraint",
            "Argument \"ipg\". Illegal value, expected >=0:  `-1'",
        ),
        (
            &["-db", "db/p5", "-info", "-line_length", "2147483648"],
            "eConvert",
            "Argument \"line_length\". Integer value is out of range:  `2147483648'",
        ),
    ];
    for &(line, code, message) in cases {
        assert_eq!(
            messages(line),
            (code, vec![message.to_string()]),
            "{line:?}"
        );
    }
    assert_eq!(
        messages(&["-db", "db/p5", "-info", "-line_length", "abc"]),
        (
            "eConvert",
            vec![
                crate::wording::losat_text("Cannot convert string 'abc' to Int8 (m_Pos = 0)")
                    .into_owned(),
                "Argument \"line_length\". Argument cannot be converted:  `abc'".to_string(),
            ]
        )
    );
}

#[test]
fn help_words_stop_the_reader_where_they_stand() {
    assert_eq!(
        parse(&words(&["-db", "x", "-info", "-h"])).unwrap_err(),
        ArgsError::Usage
    );
    assert_eq!(
        parse(&words(&["-h", "-foo"])).unwrap_err(),
        ArgsError::Usage
    );
    assert_eq!(messages(&["-foo", "-h"]).1, ["Unknown argument: \"foo\""]);
    assert_eq!(
        parse(&words(&["-info", "-help"])).unwrap_err(),
        ArgsError::Help
    );
    assert_eq!(messages(&["-h=1"]).1, ["Invalid argument: -h=1"]);
    // Known to NCBI's reader: its help with the hidden arguments (none), its XML help.
    assert_eq!(
        parse(&words(&["-info", "-help-full"])).unwrap_err(),
        ArgsError::Help
    );
    assert_eq!(
        parse(&words(&["-info", "-xmlhelp"])).unwrap_err(),
        ArgsError::XmlHelp
    );
    // `-version` is taken before the reader; the reader does not know it.
    assert_eq!(
        messages(&["-version=1"]).1,
        ["Unknown argument: \"version\""]
    );
}

#[test]
fn names_are_checked_and_looked_up_as_ncbi_does() {
    let cases: &[(&[&str], &str)] = &[
        (&["-info", "-foo.bar"], "Invalid argument name: foo.bar"),
        (&["-info", "-db.x=1"], "Invalid argument name: db.x"),
        (&["-info", "---x"], "Invalid argument name: --x"),
        (&["-info", "-\u{e9}"], "Invalid argument name: \u{e9}"),
        (&["-info", "--x"], "Unknown argument: \"-x\""),
        (
            &["-info", "-=x"],
            "Unexpected extra argument, at position # 4294967295",
        ),
        (&["-info", "-x-file"], "Unknown argument: \"x-file\""),
    ];
    for &(line, message) in cases {
        assert_eq!(
            messages(line),
            ("eInvalidArg", vec![message.to_string()]),
            "{line:?}"
        );
    }
    // `-NAME-file` and `-NAME-verbatim` are the argument NAME.
    let args = parse(&words(&["-db-file", "db/p5", "-info-FILE"])).unwrap();
    assert_eq!(args.value("db"), Some("db/p5"));
    assert!(args.given("info"));
    let args = parse(&words(&["-db-verbatim=db/p5", "-info"])).unwrap();
    assert_eq!(args.value("db"), Some("db/p5"));
    assert_eq!(
        messages(&["-info", "-info-file"]).1,
        ["Argument with this name is defined already: info"]
    );
}

#[test]
fn integer_keys_convert_as_ncbi_s_reader() {
    // `s_StringToInt8` reads a hexadecimal value, and BLAST's constraint `>=` reads its text
    // with `NStr::StringToDouble`, whose `strtod` reads it too (NCBI 2.17.0: `-line_length 0x5`
    // writes lines of 5, probe K4i `edge.blastdbcmd.line_length.hex`).
    let args = parse(&words(&["-db", "x", "-info", "-line_length", "0x10"])).unwrap();
    assert_eq!(args.integer("line_length"), Some(16));
    let args = parse(&words(&["-db", "x", "-ipg", "0X0"])).unwrap();
    assert_eq!(args.integer("ipg"), Some(0));
    // `0x` alone reads 0 in base 16, but `strtod` reads only its `0` and the constraint fails.
    assert_eq!(
        messages(&["-db", "x", "-info", "-line_length", "0x"]),
        (
            "eConstraint",
            vec!["Argument \"line_length\". Illegal value, expected >=1:  `0x'".to_string()]
        )
    );
    assert_eq!(
        messages(&["-db", "x", "-info", "-line_length", "0x0"]),
        (
            "eConstraint",
            vec!["Argument \"line_length\". Illegal value, expected >=1:  `0x0'".to_string()]
        )
    );
    // The conversion in base 16 throws the exception that is reported.
    assert_eq!(
        messages(&["-db", "x", "-info", "-line_length", "0x1p3"]),
        (
            "eConvert",
            vec![
                crate::wording::losat_text("Cannot convert string '0x1p3' to Int8 (m_Pos = 3)")
                    .into_owned(),
                "Argument \"line_length\". Argument cannot be converted:  `0x1p3'".to_string(),
            ]
        )
    );
    assert_eq!(
        messages(&["-db", "x", "-info", "-line_length", "0x80000000"]),
        (
            "eConvert",
            vec![
                "Argument \"line_length\". Integer value is out of range:  `0x80000000'"
                    .to_string()
            ]
        )
    );
    // A value that starts with a character of several bytes.
    assert_eq!(
        messages(&["-db", "x", "-ipg", "\u{3042}"]),
        (
            "eConvert",
            vec![
                crate::wording::losat_text(
                    "Cannot convert string '\\343\\201\\202' to Int8 (m_Pos = 0)"
                )
                .into_owned(),
                "Argument \"ipg\". Argument cannot be converted:  `\u{3042}'".to_string(),
            ]
        )
    );
    // blastdbcmd's own conversions (`-range`, `-mask_sequence_with`, batch lines).
    assert_eq!(string_to_int8_message("80"), Ok(80));
    assert_eq!(string_to_int8_message("+7"), Ok(7));
    assert_eq!(string_to_int8_message("-9223372036854775808"), Ok(i64::MIN));
    assert_eq!(
        string_to_int8_message("9223372036854775808").unwrap_err(),
        "Cannot convert string '9223372036854775808' to Int8, overflow"
    );
    assert_eq!(
        string_to_int8_message("12ab").unwrap_err(),
        "Cannot convert string '12ab' to Int8"
    );
    assert_eq!(
        string_to_int8_message("-").unwrap_err(),
        "Cannot convert string '-' to Int8"
    );
    assert_eq!(
        string_to_int8_message("").unwrap_err(),
        "Cannot convert string '' to Int8"
    );
    // NCBI reads up to a NUL (its loop runs while the character at `str[pos]` is not zero,
    // `CHECK_ENDPTR` sees the NUL): a batch line `12\0x` is 12.
    assert_eq!(string_to_int8_message("12\u{0}x"), Ok(12));
    assert_eq!(
        string_to_int8_message("\u{0}1").unwrap_err(),
        "Cannot convert string '\\0001' to Int8"
    );
}

/// `CArgValue::AsInteger`: an integer key is read as the number the reader accepted, so a
/// value whose digits end at a NUL (a value the library can pass) is that number for its
/// consumer too, as in NCBI (`-line_length 40\0x` wraps at 40).
#[test]
fn an_integer_key_is_the_number_the_reader_accepted() {
    let args = parse(&words(&[
        "-db",
        "db/p5",
        "-entry",
        "all",
        "-line_length",
        "40\u{0}x",
    ]))
    .unwrap();
    assert_eq!(args.integer("line_length"), Some(40));
    let args = parse(&words(&["-db", "db/p5", "-ipg", "1001\u{0}x"])).unwrap();
    assert_eq!(args.integer("ipg"), Some(1001));
    let args = parse(&words(&["-db", "db/p5", "-entry", "all"])).unwrap();
    assert_eq!(args.integer("line_length"), Some(80));
}

// NCBI's `multimap` keeps the order of insertion among equal names: of two exclusions that a
// command line breaks, NCBI reports the first of the map (blastdbcmd.cpp:928-934: `-entry_batch`
// excludes `-entry` before `-ipg` does, so `entry`'s first entry is the reverse of the former).
#[test]
fn the_first_exclusion_of_ncbi_s_multimap_is_reported() {
    assert_eq!(
        messages(&["-db", "x", "-ipg", "1", "-entry", "a", "-entry_batch", "b"]),
        (
            "eConstraint",
            vec!["Argument \"entry_batch\". Incompatible with argument:  `entry'".to_string()]
        )
    );
    assert_eq!(
        messages(&["-db", "x", "-entry", "a", "-ipg", "1"]),
        (
            "eConstraint",
            vec!["Argument \"ipg\". Incompatible with argument:  `entry'".to_string()]
        )
    );
}
