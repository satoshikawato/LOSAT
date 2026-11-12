//! The arguments of `losat blastdbcmd`: the descriptions of blastdbcmd 2.17.0
//! (`CBlastDBCmdApp::Init`) that NCBI's argument reader reads (`ncbi_app::args`), so that every
//! command line that NCBI's reader refuses is refused with NCBI's text, also for the options that
//! LOSAT's BLASTDBCMD does not serve (the program refuses those after the reader accepts them).
//! With them blastdbcmd's copies of `NStr::StringToInt8` and `NStr::PrintableString` for the
//! values it converts itself (`-range`, `-mask_sequence_with`, the lines of its batch files), not
//! yet folded into `ncbi_str` (whose `StringToInt8` ends a number at a NUL byte).

use crate::ncbi_app::args::{flag, key, Constraint, Dep, Desc, Kind};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:901-919
// ```c++
// void CBlastDBCmdApp::Init()
// {
//     HideStdArgs(fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
//
//     unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);
//     ...
//     arg_desc->AddDefaultKey(kArgDb, "dbname", "BLAST database name",
//                             CArgDescriptions::eString, "nr");
//
//     arg_desc->AddDefaultKey(kArgDbType, "molecule_type",
//                             "Molecule type stored in BLAST database",
//                             CArgDescriptions::eString, "guess");
//     arg_desc->SetConstraint(kArgDbType, &(*new CArgAllow_Strings,
//                                         "nucl", "prot", "guess"));
// ```
/// The arguments of `CBlastDBCmdApp::Init` (blastdbcmd.cpp:901-1167) in the order it describes
/// them; NCBI's reader adds the help flags before them and `-logfile` and `-conffile` after them
/// (`ncbi_app::args`). The file keys are `eInputFile` (`-entry_batch`, `-ipg_batch`,
/// `-taxidlist`) and `eOutputFile` (`-out`; blastdbcmd.cpp:928-931, 945-948, 1029-1030,
/// 1042-1043).
pub(super) const DESCS: &[Desc] = &[
    key("db", Kind::Text, Some("nr"), None),
    key(
        "dbtype",
        Kind::Text,
        Some("guess"),
        Some(Constraint::Strings(&["guess", "nucl", "prot"])),
    ),
    key("entry", Kind::Text, None, None),
    key("entry_batch", Kind::InputFile, None, None),
    key("ipg", Kind::Integer, None, Some(Constraint::AtLeast(0))),
    key("ipg_batch", Kind::InputFile, None, None),
    flag("info"),
    flag("metadata"),
    key("metadata_output_prefix", Kind::Text, None, None),
    flag("tax_info"),
    key("range", Kind::Text, None, None),
    key(
        "strand",
        Kind::Text,
        Some("plus"),
        Some(Constraint::Strings(&["minus", "plus"])),
    ),
    key("mask_sequence_with", Kind::Text, None, None),
    key("taxids", Kind::Text, None, None),
    key("taxidlist", Kind::InputFile, None, None),
    flag("no_taxid_expansion"),
    key("out", Kind::OutputFile, Some("-"), None),
    key("outfmt", Kind::Text, Some("%f"), None),
    flag("target_only"),
    flag("get_dups"),
    key(
        "line_length",
        Kind::Integer,
        Some("80"),
        Some(Constraint::AtLeast(1)),
    ),
    flag("ctrl_a"),
    flag("show_blastdb_search_path"),
    key("list", Kind::Text, None, None),
    flag("remove_redundant_dbs"),
    flag("recursive"),
    key("list_outfmt", Kind::Text, Some("%f %p"), None),
    flag("exact_length"),
    flag("long_seqids"),
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:955-961
// ```c++
//     arg_desc->AddFlag("info", "Print BLAST database information", true);
//     // All other options to this program should be here
//     const char* exclusions[]  = { "entry", "entry_batch", "outfmt", "strand",
//         "target_only", "ctrl_a", "get_dups", "pig", "range",
//         "mask_sequence", "list", "remove_redundant_dbs", "recursive",
//         "list_outfmt", "metadata", "metadata_output_prefix", kArgTaxIdListFile.c_str(), kArgTaxIdList.c_str(),
//         kArgNoTaxIdExpansion.c_str()};
// ```
// `pig` and `mask_sequence` name no argument, so they never apply.
const INFO_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "metadata",
    "metadata_output_prefix",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:969-977
// ```c++
//     const char* exclusions_m[]  = { "entry", "entry_batch", "outfmt", "strand",
//         "target_only", "ctrl_a", "get_dups", "pig", "range",
//         "mask_sequence", "list", "remove_redundant_dbs", "recursive",
//         "list_outfmt", "info", kArgTaxIdListFile.c_str(), kArgTaxIdList.c_str(),
//         kArgNoTaxIdExpansion.c_str()};
// ```
const METADATA_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "info",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:994-998
// ```c++
//     const char* tax_info_exclusions[]  = { "info", "entry", "entry_batch", "strand",
//         "target_only", "ctrl_a", "get_dups", "pig", "range",
//         "mask_sequence", "list", "remove_redundant_dbs", "recursive",
//         "list_outfmt", kArgTaxIdListFile.c_str(), kArgTaxIdList.c_str(),
//         kArgNoTaxIdExpansion.c_str() };
// ```
const TAX_INFO_EXCLUSIONS: &[&str] = &[
    "info",
    "entry",
    "entry_batch",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1111-1113
// ```c++
//     const char* exclusions_discovery[]  = { "entry", "entry_batch", "outfmt",
//         "strand", "target_only", "ctrl_a", "get_dups", "pig", "range", kArgDb.c_str(),
//         "info", "mask_sequence", "line_length" };
// ```
const DISCOVERY_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "db",
    "info",
    "mask_sequence",
    "line_length",
];

/// The `SetDependency` calls of blastdbcmd (blastdbcmd.cpp:901-1167), in their order; NCBI's
/// reader holds them in its `multimap` (`ncbi_app::args::dependencies`).
pub(super) fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    let mut deps = Vec::new();
    let mut set = |first: &'static str, dep: Dep, second: &'static str| {
        deps.push((first, dep, second));
    };
    for second in ["entry", "range", "strand", "mask_sequence_with"] {
        set("entry_batch", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "target_only", "ipg_batch"] {
        set("ipg", Dep::Excludes, second);
    }
    for second in [
        "entry",
        "entry_batch",
        "range",
        "strand",
        "mask_sequence_with",
    ] {
        set("ipg_batch", Dep::Excludes, second);
    }
    for &second in INFO_EXCLUSIONS {
        set("info", Dep::Excludes, second);
    }
    for &second in METADATA_EXCLUSIONS {
        set("metadata", Dep::Excludes, second);
    }
    set("metadata_output_prefix", Dep::Requires, "metadata");
    for &second in TAX_INFO_EXCLUSIONS {
        set("tax_info", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig"] {
        set("taxids", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig", "taxids"] {
        set("taxidlist", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig"] {
        set("no_taxid_expansion", Dep::Excludes, second);
    }
    set("get_dups", Dep::Excludes, "target_only");
    for &second in DISCOVERY_EXCLUSIONS {
        for first in [
            "list",
            "recursive",
            "remove_redundant_dbs",
            "list_outfmt",
            "show_blastdb_search_path",
        ] {
            set(first, Dep::Excludes, second);
        }
    }
    for second in ["list", "recursive", "list_outfmt", "remove_redundant_dbs"] {
        set("show_blastdb_search_path", Dep::Excludes, second);
    }
    set("exact_length", Dep::Requires, "info");
    set("long_seqids", Dep::Excludes, "info");
    deps
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-871
// ```c++
// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
// {
//     ...
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     ...
//     }
//     SIZE_TYPE pos0 = pos;
//     ...
//     Int8 n = 0;
//     Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
//     Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);
//     ...
//     while (char ch = str[pos]) {
//         ...
//         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
//             break;
//         }
//         // Overflow check
//         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
//             S2N_CONVERT_ERROR_OVERFLOW(Int8);
//         }
//         n *= base;
//         n += delta;
//         pos++;
//     }
//
//     // Last checks
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     ...
//     n = sign ? -n : n;
//     CHECK_ENDPTR(Int8);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:546-557
// ```c++
//     s += "Cannot convert string '";
//     s += NStr::PrintableString(str);
//     s += "' to ";
//     s += to_type;
//     if ( !msg.empty() ) {
//         s += ", ";
//         s += msg;
//     }
// ```
// The `CStringException` adds its position (`m_Pos = N`, ncbistr.hpp:4438) to the message.
/// `NStr::StringToInt8(value)`: the number, or the message of NCBI's exception as
/// `GetMsg` gives it (without the position).
pub(super) fn string_to_int8_message(value: &str) -> Result<i64, String> {
    string_to_int8_at(value).map_err(|(message, _)| message)
}

/// `NStr::StringToInt8(value)`: the number, or the message of NCBI's exception and the
/// position the exception holds (`m_Pos`).
fn string_to_int8_at(value: &str) -> Result<i64, (String, usize)> {
    let bytes = value.as_bytes();
    let error = |extra: &str, pos: usize| {
        (
            format!(
                "Cannot convert string '{}' to Int8{extra}",
                printable_string(value)
            ),
            pos,
        )
    };
    let mut pos = 0;
    let sign = match bytes.first() {
        Some(b'-') => {
            pos = 1;
            true
        }
        Some(b'+') => {
            pos = 1;
            false
        }
        _ => false,
    };
    let pos0 = pos;
    let limdiv = i64::MAX / 10;
    let limoff = i64::MAX % 10 + i64::from(sign);
    let mut n: i64 = 0;
    while let Some(&byte) = bytes.get(pos).filter(|byte| byte.is_ascii_digit()) {
        let delta = i64::from(byte - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(", overflow", pos));
        }
        // `i64::MIN` passes the overflow check as its magnitude, one past `i64::MAX`.
        n = n.wrapping_mul(10).wrapping_add(delta);
        pos += 1;
    }
    if pos == pos0 || pos < bytes.len() {
        return Err(error("", pos));
    }
    Ok(if sign { n.wrapping_neg() } else { n })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3743-3810
// ```c++
//         case '\a':
//             if (lang == eLanguage_C)
//                 c = 'a';
//             ...
//         case '\n':
//             if (!(mode & NStr::fNewLine_Passthru))
//                 c = 'n';
//             /*FALLTHRU*/
//         case '\\':
//         case '\'':
//         case '"':
//             break;
//         case '&':
//             ...
//             continue;
//         case '?':
//             if (lang == eLanguage_C) {
//                 if (i  &&  str[i - 1] == '?')
//                     break;
//                 if (i < str.size() - 1  &&  str[i + 1] == '?')
//                     break;
//             }
//             continue;
//         default:
//             ...
//             if (!isprint((unsigned char) c)) {
//                 octal = true;
//                 break;
//             }
//             continue;
//         }
//         ...
//         out->put('\\');
//         if (c == '\n') {
//             out->write("n\\\n", 3);
//         } else if (octal) {
//             bool reduce;
//             if (!(mode & NStr::fPrintable_Full)) {
//                 reduce = (i == str.size() - 1  ||
//                           str[i + 1] < '0' || '7' < str[i + 1] ? true : false);
// ```
// `PrintableString(str)` has the mode 0 (`fNewLine_Quote`): a line feed becomes `\n`, the
// branch of `c == '\n'` is that of `fNewLine_Passthru`.
/// `NStr::PrintableString(str)`: C escapes for the quotes, the backslash, the control
/// characters and `??` (trigraphs), and octal for the other bytes that are not printable
/// ASCII, in as few digits as the next byte allows.
pub(super) fn printable_string(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(bytes.len());
    for (i, &byte) in bytes.iter().enumerate() {
        let escape = match byte {
            0x07 => Some("\\a".to_string()),
            0x08 => Some("\\b".to_string()),
            0x0c => Some("\\f".to_string()),
            b'\r' => Some("\\r".to_string()),
            b'\t' => Some("\\t".to_string()),
            0x0b => Some("\\v".to_string()),
            b'\n' => Some("\\n".to_string()),
            b'\\' | b'\'' | b'"' => Some(format!("\\{}", byte as char)),
            b'?' if (i > 0 && bytes[i - 1] == b'?') || bytes.get(i + 1) == Some(&b'?') => {
                Some("\\?".to_string())
            }
            0x20..=0x7e => None,
            _ => {
                // `char` is signed: a byte above 0x7f is below '0'.
                let reduce = bytes
                    .get(i + 1)
                    .is_none_or(|&next| !(b'0'..=b'7').contains(&next));
                let digits = [byte >> 6, (byte >> 3) & 7, byte & 7];
                // NCBI ncbistr.cpp:3817-3820: `if (v || !reduce)` writes the first digit (and
                // then every digit).
                let skip = if digits[0] != 0 || !reduce {
                    0
                } else if digits[1] != 0 {
                    1
                } else {
                    2
                };
                let octal: String = digits[skip..]
                    .iter()
                    .map(|digit| char::from(b'0' + digit))
                    .collect();
                Some(format!("\\{octal}"))
            }
        };
        match escape {
            Some(escape) => out.push_str(&escape),
            None => out.push(byte as char),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncbi_app::args::{ArgException, Args, ArgsError};

    fn words(line: &[&str]) -> Vec<String> {
        line.iter().map(|word| word.to_string()).collect()
    }

    /// The arguments of `line` as NCBI's reader reads blastdbcmd's.
    fn parse(line: &[String]) -> Result<Args, ArgsError> {
        crate::ncbi_app::args::parse(&super::super::TOOL, line)
    }

    fn messages(line: &[&str]) -> (&'static str, Vec<String>) {
        match parse(&words(line)) {
            Err(ArgsError::Invalid(ArgException { code, messages })) => (code, messages),
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
            (&["-db", "db/p5", "-info", "-foo"], "eInvalidArg", "Unknown argument: \"foo\""),
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
            (&["-db", "db/p5", "-info=1"], "eInvalidArg", "Invalid argument: -info=1"),
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
                    "Cannot convert string 'abc' to Int8 (m_Pos = 0)".to_string(),
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
        assert!(matches!(
            parse(&words(&["-info", "-xmlhelp"])),
            Err(ArgsError::Refused(_))
        ));
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
        // NCBI's constraint `>=` reads a hexadecimal value as a double: LOSAT's refusal.
        assert!(matches!(
            parse(&words(&["-db", "x", "-info", "-line_length", "0x10"])),
            Err(ArgsError::Refused(_))
        ));
        // A value that starts with a character of several bytes.
        assert_eq!(
            messages(&["-db", "x", "-ipg", "\u{3042}"]),
            (
                "eConvert",
                vec![
                    "Cannot convert string '\\343\\201\\202' to Int8 (m_Pos = 0)".to_string(),
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
            string_to_int8_at("-").unwrap_err(),
            ("Cannot convert string '-' to Int8".to_string(), 1)
        );
        assert_eq!(
            string_to_int8_at("").unwrap_err(),
            ("Cannot convert string '' to Int8".to_string(), 0)
        );
    }

    #[test]
    fn printable_string_escapes_as_ncbi() {
        assert_eq!(printable_string("a'b\"c\\d"), "a\\'b\\\"c\\\\d");
        assert_eq!(printable_string("x\ty"), "x\\ty");
        assert_eq!(printable_string("1\n2"), "1\\n2");
        assert_eq!(printable_string("a??b?c"), "a\\?\\?b?c");
        assert_eq!(printable_string("\u{1}9"), "\\19");
        assert_eq!(printable_string("\u{1}1"), "\\0011");
        assert_eq!(printable_string("é"), "\\303\\251");
    }

    #[test]
    fn the_dependency_order_is_that_of_ncbi_s_multimap() {
        let deps = crate::ncbi_app::args::dependencies(&super::super::TOOL);
        assert!(deps.windows(2).all(|pair| pair[0].0 <= pair[1].0));
        // Among the entries of `entry`, the reverse of `entry_batch excludes entry` comes first.
        let first_of_entry = deps.iter().find(|dep| dep.0 == "entry").unwrap();
        assert_eq!(*first_of_entry, ("entry", Dep::Excludes, "entry_batch"));
    }
}
