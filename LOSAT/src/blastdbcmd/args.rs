//! The arguments of `losat blastdbcmd`: the descriptions of blastdbcmd 2.17.0
//! (`CBlastDBCmdApp::Init`) that NCBI's argument reader reads (`ncbi_app::args`), so that every
//! command line that NCBI's reader refuses is refused with NCBI's text, also for the options that
//! LOSAT's BLASTDBCMD does not serve (the program refuses those after the reader accepts them).
//! With them `string_to_int8_message`, `NStr::StringToInt8` (`ncbi_str`) for the values
//! blastdbcmd converts itself (`-mask_sequence_with`, the lines of its batch files).

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

/// `NStr::StringToInt8(value)` (`ncbi_str::string_to_int8`): the number, or the message of
/// NCBI's exception as `GetMsg` gives it (without the position).
pub(super) fn string_to_int8_message(value: &str) -> Result<i64, String> {
    crate::ncbi_str::string_to_int8(value.as_bytes(), 0)
        .map_err(|error| error.get_msg().to_string())
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
            string_to_int8_message("-").unwrap_err(),
            "Cannot convert string '-' to Int8"
        );
        assert_eq!(
            string_to_int8_message("").unwrap_err(),
            "Cannot convert string '' to Int8"
        );
        // NCBI reads up to a NUL (`while (char ch = str[pos])`, `CHECK_ENDPTR` sees the NUL): a
        // batch line `12\0x` is 12.
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

    #[test]
    fn the_dependency_order_is_that_of_ncbi_s_multimap() {
        let deps = crate::ncbi_app::args::dependencies(&super::super::TOOL);
        assert!(deps.windows(2).all(|pair| pair[0].0 <= pair[1].0));
        // Among the entries of `entry`, the reverse of `entry_batch excludes entry` comes first.
        let first_of_entry = deps.iter().find(|dep| dep.0 == "entry").unwrap();
        assert_eq!(*first_of_entry, ("entry", Dep::Excludes, "entry_batch"));
    }
}
