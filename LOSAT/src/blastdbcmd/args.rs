//! The arguments of `losat blastdbcmd` as NCBI's argument reader (`CArgDescriptions`) reads
//! them, with the descriptions of blastdbcmd 2.17.0 (`CBlastDBCmdApp::Init`): every command line
//! that NCBI's reader refuses is refused with NCBI's text, also for the options that LOSAT's
//! BLASTDBCMD does not serve (the program refuses those after this reader accepts them).

use std::collections::{BTreeMap, BTreeSet};

/// What an option holds (`CArgDescriptions::EType` of a key, or a flag).
///
/// The file types of NCBI (`eInputFile`, `eOutputFile`) are strings here: blastdbcmd describes
/// them without `fPreOpen`, so the reader does not open them (ncbiargs.cpp:744-752).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Flag,
    Text,
    Integer,
}

/// A constraint on a key's value (`CArgAllow`).
#[derive(Debug, Clone, Copy)]
enum Constraint {
    /// `CArgAllow_Strings`: one of the values, compared with case; listed in the order of its
    /// `std::set` (ncbiargs.cpp:4605-4630).
    Strings(&'static [&'static str]),
    /// `CArgAllowValuesGreaterThanOrEqual(min)` (blast_input_aux.hpp:101-121).
    AtLeast(i64),
}

impl Constraint {
    fn usage(self) -> String {
        match self {
            Constraint::Strings(values) => values
                .iter()
                .map(|value| format!("`{value}'"))
                .collect::<Vec<_>>()
                .join(", "),
            Constraint::AtLeast(min) => format!(">={min}"),
        }
    }
}

/// One argument description (`CArgDesc`).
struct Desc {
    name: &'static str,
    kind: Kind,
    default: Option<&'static str>,
    constraint: Option<Constraint>,
}

const fn flag(name: &'static str) -> Desc {
    Desc {
        name,
        kind: Kind::Flag,
        default: None,
        constraint: None,
    }
}

const fn key(
    name: &'static str,
    kind: Kind,
    default: Option<&'static str>,
    constraint: Option<Constraint>,
) -> Desc {
    Desc {
        name,
        kind,
        default,
        constraint,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:901-918
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
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2207-2219
// ```c++
//     if ( m_AutoHelp ) {
//         AddFlag(s_AutoHelp,
//                 "Print USAGE and DESCRIPTION;  ignore all other parameters");
//     }
//     AddFlag(s_AutoHelpFull,
//             "Print USAGE, DESCRIPTION and ARGUMENTS;"
//             " ignore all other parameters");
// ```
/// The arguments in the order NCBI describes them (`m_KeyFlagArgs`): the two help flags of
/// `CArgDescriptions`, those of `CBlastDBCmdApp::Init` (blastdbcmd.cpp:901-1167), then the
/// standard `-logfile` and `-version` (`-help-full`, `-xmlhelp`, `-conffile`, `-version-full`
/// and `-dryrun` are hidden, so the reader does not know them).
const DESCS: &[Desc] = &[
    flag("h"),
    flag("help"),
    key("db", Kind::Text, Some("nr"), None),
    key(
        "dbtype",
        Kind::Text,
        Some("guess"),
        Some(Constraint::Strings(&["guess", "nucl", "prot"])),
    ),
    key("entry", Kind::Text, None, None),
    key("entry_batch", Kind::Text, None, None),
    key("ipg", Kind::Integer, None, Some(Constraint::AtLeast(0))),
    key("ipg_batch", Kind::Text, None, None),
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
    key("taxidlist", Kind::Text, None, None),
    flag("no_taxid_expansion"),
    key("out", Kind::Text, Some("-"), None),
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
    key("logfile", Kind::Text, None, None),
    flag("version"),
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Dep {
    Requires,
    Excludes,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:957-966
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

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:996-1004
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

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1118-1121
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2618-2629
// ```c++
// void CArgDescriptions::SetDependency(const string& arg1,
//                                      EDependency   dep,
//                                      const string& arg2)
// {
//     m_Dependencies.insert(TDependencies::value_type(arg1,
//         SArgDependency(arg2, dep)));
//     if (dep == eExcludes) {
//         // Exclusions must work in both directions
//         m_Dependencies.insert(TDependencies::value_type(arg2,
//             SArgDependency(arg1, dep)));
//     }
// }
// ```
/// The dependencies of blastdbcmd (`SetDependency` calls of blastdbcmd.cpp:901-1167, in their
/// order) as NCBI's `multimap` holds them: ordered by the first argument's name, and in the
/// order of insertion among equal names. Which of several conflicts NCBI reports depends on
/// this order.
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    let mut deps = Vec::new();
    let mut set = |first: &'static str, dep: Dep, second: &'static str| {
        deps.push((first, dep, second));
        if dep == Dep::Excludes {
            deps.push((second, dep, first));
        }
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
    // A stable sort keeps the insertion order among equal names, as `multimap::insert` does.
    deps.sort_by_key(|&(first, _, _)| first);
    deps
}

/// Why the arguments were not read.
#[derive(Debug, PartialEq, Eq)]
pub(super) enum ArgsError {
    /// `-h` (`CArgHelpException::eHelp`): the short usage on stdout, exit status 0.
    Usage,
    /// `-help` (`CArgHelpException::eHelpFull`): the full description on stdout, exit status 0.
    Help,
    /// A `CArgException` with the code name `code` (`eInvalidArg`, ...): `messages` are the
    /// messages of the exceptions it rethrows, oldest first, then its own.
    Invalid {
        code: &'static str,
        messages: Vec<String>,
    },
    /// A value that NCBI's reader takes in a way LOSAT does not reproduce: LOSAT's refusal.
    Refused(String),
}

fn invalid(code: &'static str, message: String) -> ArgsError {
    ArgsError::Invalid {
        code,
        messages: vec![message],
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:93-97
// ```c++
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
pub(super) fn arg_message(name: &str, what: &str, attr: &str) -> String {
    let name = if name.is_empty() { "...." } else { name };
    if attr.is_empty() {
        format!("Argument \"{name}\". {what}")
    } else {
        format!("Argument \"{name}\". {what}:  `{attr}'")
    }
}

/// The arguments of one command line: those given and the defaults NCBI adds.
#[derive(Debug)]
pub(super) struct Args {
    /// The value of each argument that has one: the value given (`None` for a flag), else the
    /// default of a key that no given argument excludes.
    values: BTreeMap<&'static str, Option<String>>,
    /// The names given on the command line, in their order.
    given: Vec<&'static str>,
}

impl Args {
    /// Whether the argument `name` was given on the command line.
    pub(super) fn given(&self, name: &str) -> bool {
        self.given.contains(&name)
    }

    /// The value of the key `name`: given, or its default.
    pub(super) fn value(&self, name: &str) -> Option<&str> {
        self.values.get(name).and_then(|value| value.as_deref())
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2856-2932
// ```c++
// bool CArgDescriptions::x_CreateArg(const string& arg1,
//                                    bool have_arg2, const string& arg2,
//                                    unsigned* n_plain, CArgs& args) const
// {
//     ...
//     if (*n_plain == kMax_UInt || m_PositionalMode == ePositionalMode_Loose) {
//         // Check for the s_ArgDelimiter delimiter
//         if (arg1.compare(s_ArgDelimiter) == 0) {
//             if (*n_plain == kMax_UInt) {
//                 *n_plain = 0;  // pos.args started
//             }
//             return false;
//         }
//         ...
//         if ((arg1.length() > 1)  &&  arg1[0] == '-') {
//             name = arg1.substr(1);
//             ...
//             // Check for '=' in the arg1
//             size_t eq = name.find('=');
//             if (eq != NPOS) {
//                 name = name.substr(0, eq);
//             }
//             ...
//         } else {
//             ...
//             if (*n_plain == kMax_UInt) {
//                 *n_plain = 0;  // pos.args started
//             }
//         }
//     }
//     ...
//     if (*n_plain != kMax_UInt && !is_keyflag) {
//         ...
//         (*n_plain)++;
//
//         // Check for too many positional arguments
//         if (kMax_UInt - m_nExtraOpt > m_nExtra + m_PosArgs.size()  &&
//             *n_plain > m_PosArgs.size() + m_nExtra + m_nExtraOpt) {
//             NCBI_THROW(CArgException,eSynopsis,
//                 "Too many positional arguments (" +
//                 NStr::UIntToString(*n_plain) +
//                 "), the offending value: "+ arg1);
//         }
//     }
// ```
/// Reads `words` (the command line after the program's name, without the words NCBI's
/// application takes before the reader: `args_for_reader` in `blastdbcmd::run_with`) as NCBI's
/// `CArgDescriptions::CreateArgs` reads blastdbcmd's: each word in its order, then the
/// dependencies and the defaults (`x_PostCheck`). blastdbcmd has no positional argument, so a
/// word that is not an option, and every word after `--`, is one too many.
pub(super) fn parse(words: &[String]) -> Result<Args, ArgsError> {
    let mut values: BTreeMap<&'static str, Option<String>> = BTreeMap::new();
    let mut given: Vec<&'static str> = Vec::new();
    // `n_plain`: `None` until the positional arguments start (NCBI's `kMax_UInt`).
    let mut n_plain: Option<u32> = None;
    let mut at = 0;
    while at < words.len() {
        let arg1 = words[at].as_str();
        let arg2 = words.get(at + 1);
        let mut name = "";
        if n_plain.is_none() {
            if arg1 == "--" {
                n_plain = Some(0);
                at += 1;
                continue;
            }
            if arg1.len() > 1 && arg1.starts_with('-') {
                name = &arg1[1..];
                if let Some(eq) = name.find('=') {
                    name = &name[..eq];
                }
            } else {
                n_plain = Some(0);
            }
        }
        if let Some(count) = n_plain {
            let count = count + 1;
            return Err(invalid(
                "eSynopsis",
                format!("Too many positional arguments ({count}), the offending value: {arg1}"),
            ));
        }
        let (desc, value, used) = create_arg(arg1, name, arg2)?;
        if values.contains_key(desc.name) {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1928-1941
            // ```c++
            //     CArgs::TArgsI arg_it = x_Find(arg->GetName());
            //     if ( arg_it !=  m_Args.end()) {
            //         ...
            //                 NCBI_THROW(CArgException,eSynopsis,
            //                    "Argument with this name is defined already: "
            //                    + arg->GetName());
            // ```
            return Err(invalid(
                "eSynopsis",
                format!("Argument with this name is defined already: {}", desc.name),
            ));
        }
        values.insert(desc.name, value);
        given.push(desc.name);
        at += if used { 2 } else { 1 };
    }
    post_check(&mut values, &given)?;
    Ok(Args { values, given })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3007-3081
// ```c++
//     if (it == m_Args.end()) {
//         if ( name.empty() ) {
//             NCBI_THROW(CArgException,eInvalidArg,
//                     "Unexpected extra argument, at position # " +
//                     NStr::UIntToString(n_plain));
//         } else {
//             NCBI_THROW(CArgException,eInvalidArg,
//                     "Unknown argument: \"" + name + "\"");
//         }
//     }
//     ...
//     if ( s_IsFlag(arg) ) {
//         x_CheckAutoHelp(arg1);
//         eq_separator = arg1.length() > name.length()  &&
//             (arg1[name.length() + 1] == '=');
//         if ( eq_separator ) {
//             NCBI_THROW(CArgException,eInvalidArg, "Invalid argument: " + arg1);
//         }
//     }
//
//     // Check value separated by '='
//     string arg_val;
//     if ( s_IsKey(arg) && !confidential) {
//         eq_separator = arg1.length() > name.length()  &&
//             (arg1[name.length() + 1] == '=');
//         ...
//     }
//
//     // Get argument value
//     string value;
//     if ( !eq_separator  &&  !no_separator ) {
//         if ( !s_IsKey(arg)  || (confidential && conf_method == 't')) {
//             value = arg1;
//         }
//         else {
//             // <key> <value> arg  -- advance from the arg.name to the arg.value
//             if ( !have_arg2 ) {
//                 ...
//                 NCBI_THROW(CArgException,eNoArg,s_ArgExptMsg(arg1,
//                     "Value is missing", kEmptyStr));
//             }
//             value = arg2;
//             arg2_used = true;
//         }
//     }
//     else {
//         ...
//             arg_val = arg1.substr(name.length() + 2);
//         ...
//         value = arg_val;
//     }
// ```
/// One option word `arg1` with the name `name` (and the next word `arg2`, the value of a key
/// written apart): its description, its value (`None` for a flag), and whether `arg2` was used.
fn create_arg(
    arg1: &str,
    name: &str,
    arg2: Option<&String>,
) -> Result<(&'static Desc, Option<String>, bool), ArgsError> {
    let Some(desc) = DESCS.iter().find(|desc| desc.name == name) else {
        if name.is_empty() {
            // `n_plain` is `kMax_UInt` here: the positional arguments have not started.
            return Err(invalid(
                "eInvalidArg",
                format!("Unexpected extra argument, at position # {}", u32::MAX),
            ));
        }
        return Err(invalid(
            "eInvalidArg",
            format!("Unknown argument: \"{name}\""),
        ));
    };
    let eq_separator = arg1.as_bytes().get(name.len() + 1) == Some(&b'=');
    if desc.kind == Kind::Flag {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2840-2847
        // ```c++
        // void CArgDescriptions::x_CheckAutoHelp(const string& arg) const
        // {
        //     if (arg.compare(string("-") + s_AutoHelp) == 0) {
        //         if (m_AutoHelp) {
        //             NCBI_THROW(CArgHelpException,eHelp,kEmptyStr);
        //         }
        //     } else if (arg.compare(string("-") + s_AutoHelpFull) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpFull,kEmptyStr);
        // ```
        match arg1 {
            "-h" => return Err(ArgsError::Usage),
            "-help" => return Err(ArgsError::Help),
            _ => {}
        }
        if eq_separator {
            return Err(invalid("eInvalidArg", format!("Invalid argument: {arg1}")));
        }
        return Ok((desc, None, false));
    }
    let (value, used) = if eq_separator {
        (arg1[name.len() + 2..].to_string(), false)
    } else {
        match arg2 {
            Some(value) => (value.clone(), true),
            None => return Err(invalid("eNoArg", arg_message(arg1, "Value is missing", ""))),
        }
    };
    process_value(desc, &value)?;
    Ok((desc, Some(value), used))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1228-1262
// ```c++
//     // Check against additional (user-defined) constraints, if any imposed
//     if ( m_Constraint ) {
//         bool err = false;
//         try {
//             bool check = m_Constraint->Verify(value);
//             ...
//                 err = !check;
//             ...
//         if (err) {
//             ...
//                 string err_msg;
//                 ...
//                     err_msg = "Illegal value, expected ";
//                 ...
//                 NCBI_THROW(CArgException, eConstraint, s_ArgExptMsg(GetName(),
//                            err_msg + m_Constraint->GetUsage(),value));
// ```
/// `CArgDesc::ProcessArgument`: the conversion of `value` to the key's type, then its
/// constraint.
fn process_value(desc: &Desc, value: &str) -> Result<(), ArgsError> {
    let number = match desc.kind {
        Kind::Integer => Some(integer_value(desc.name, value)?),
        Kind::Text | Kind::Flag => None,
    };
    let allowed = match desc.constraint {
        None => true,
        Some(Constraint::Strings(values)) => values.contains(&value),
        // `NStr::StringToDouble(value) >= m_MinValue`: the decimal integer `integer_value`
        // accepted is the same number as a double.
        Some(Constraint::AtLeast(min)) => number.is_some_and(|number| number >= min),
    };
    match desc.constraint {
        Some(constraint) if !allowed => Err(invalid(
            "eConstraint",
            arg_message(
                desc.name,
                &format!("Illegal value, expected {}", constraint.usage()),
                value,
            ),
        )),
        _ => Ok(()),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:118-129
// ```c++
// Int8 s_StringToInt8(const string& value)
// {
//     try {
//         return NStr::StringToInt8(value);
//     } catch (const CStringException&) {
//         if (NStr::StartsWith(value, "0x", NStr::eNocase)) {
//             return NStr::StringToInt8(value, 0, 16);
//         } else {
//             throw;
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:351-360, 383-390
// ```c++
// inline CArg_Int8::CArg_Int8(const string& name, const string& value)
//     : CArg_String(name, value)
// {
//     try {
//         m_Integer = s_StringToInt8(value);
//     } catch (const CException& e) {
//         NCBI_RETHROW(e, CArgException, eConvert, s_ArgExptMsg(GetName(),
//             "Argument cannot be converted", value));
//     }
// }
// ...
// inline CArg_Integer::CArg_Integer(const string& name, const string& value)
//     : CArg_Int8(name, value)
// {
//     if (m_Integer < kMin_Int  ||  m_Integer > kMax_Int) {
//         NCBI_THROW(CArgException, eConvert, s_ArgExptMsg(GetName(),
//             "Integer value is out of range", value));
//     }
// }
// ```
/// The value of an `eInteger` key (`CArg_Integer`). A hexadecimal value (`0x...`) is LOSAT's
/// refusal: NCBI's constraint check reads it as a double (`NStr::StringToDouble`), which LOSAT
/// does not reproduce.
fn integer_value(name: &str, value: &str) -> Result<i64, ArgsError> {
    let number = match string_to_int8(value) {
        Ok(number) => number,
        Err(_) if value.len() >= 2 && value[..2].eq_ignore_ascii_case("0x") => {
            return Err(ArgsError::Refused(format!(
                "the hexadecimal value '{value}' of -{name} is not supported by LOSAT's BLASTDBCMD"
            )));
        }
        Err(error) => {
            return Err(ArgsError::Invalid {
                code: "eConvert",
                messages: vec![
                    error,
                    arg_message(name, "Argument cannot be converted", value),
                ],
            });
        }
    };
    if number < i64::from(i32::MIN) || number > i64::from(i32::MAX) {
        return Err(invalid(
            "eConvert",
            arg_message(name, "Integer value is out of range", value),
        ));
    }
    Ok(number)
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
/// `NStr::StringToInt8(value)` in base 10 with no flags: the number, or the message of
/// NCBI's exception.
fn string_to_int8(value: &str) -> Result<i64, String> {
    let bytes = value.as_bytes();
    let error = |extra: &str, pos: usize| {
        format!(
            "Cannot convert string '{}' to Int8{extra} (m_Pos = {pos})",
            printable_string(value)
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp (s_PrintableString, eLanguage_C,
// mode 0)
// ```c++
//         switch (c) {
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
/// `NStr::PrintableString(str)`: C escapes for the quotes, the backslash, the control
/// characters and `??` (trigraphs), and octal for the other bytes that are not printable
/// ASCII, in as few digits as the next byte allows.
fn printable_string(text: &str) -> String {
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
            b'\n' => Some("\\n\\\n".to_string()),
            b'\\' | b'\'' | b'"' => Some(format!("\\{}", byte as char)),
            b'?' if (i > 0 && bytes[i - 1] == b'?') || bytes.get(i + 1) == Some(&b'?') => {
                Some("\\?".to_string())
            }
            0x20..=0x7e => None,
            _ => {
                // `char` is signed: a byte above 0x7f is below '0'.
                let reduce = bytes
                    .get(i + 1)
                    .map_or(true, |&next| !(b'0'..=b'7').contains(&next));
                let digits = [byte >> 6, (byte >> 3) & 7, byte & 7];
                let skip = if !reduce {
                    0
                } else if digits[0] != 0 {
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3162-3296
// ```c++
//     ITERATE(TDependencies, dep, m_Dependencies) {
//         // Skip missing and empty arguments
//         if (!args.Exist(dep->first)  ||  !args[dep->first]) {
//             continue;
//         }
//         switch ( dep->second.m_Dep ) {
//         case eRequires:
//             require.insert(make_pair(dep->second.m_Arg,dep->first));
//             ...
//             break;
//         case eExcludes:
//             // Excluded exists and is not empty?
//             if (args.Exist(dep->second.m_Arg)  &&  args[dep->second.m_Arg]) {
//                 NCBI_THROW(CArgException, eConstraint,
//                     s_ArgExptMsg(dep->second.m_Arg,
//                     "Incompatible with argument", dep->first));
//             }
//             exclude.insert(dep->second.m_Arg);
//             ...
//     // Set default values (if available) for the arguments not defined
//     // in the command line.
//     ITERATE (list<const CArgDesc*>, it, def_args) {
//         const CArgDesc& arg = **it;
//
//         // Nothing to do if defined in the command-line
//         if ( args.Exist(arg.GetName()) ) {
//             continue;
//         }
//
//         if (require.find(arg.GetName()) != require.end() ) {
//             string requester(require.find(arg.GetName())->second);
//             // Required argument must be present
//             NCBI_THROW(CArgException, eConstraint,
//                 s_ArgExptMsg(arg.GetName(),
//                 "Must be specified, as it is required by argument", requester));
//         }
//
//         if (exclude.find(arg.GetName()) != exclude.end()) {
//             CRef<CArg_ExcludedValue> arg_exvalue(
//                 new CArg_ExcludedValue(arg.GetName()));
//             // Add the excluded-value argument to "args"
//             args.Add(arg_exvalue);
//             continue;
//         }
//         // Use default argument value
//         try {
//             CRef<CArgValue> arg_value(arg.ProcessDefault());
// ```
/// `x_PostCheck`: the dependencies of the given arguments, then the defaults of the keys that
/// were not given and that no given argument excludes.
fn post_check(
    values: &mut BTreeMap<&'static str, Option<String>>,
    given: &[&'static str],
) -> Result<(), ArgsError> {
    let mut require: BTreeMap<&str, &str> = BTreeMap::new();
    let mut exclude: BTreeSet<&str> = BTreeSet::new();
    for (first, dep, second) in dependencies() {
        if !given.contains(&first) {
            continue;
        }
        match dep {
            Dep::Requires => {
                require.entry(second).or_insert(first);
            }
            Dep::Excludes => {
                if given.contains(&second) {
                    return Err(invalid(
                        "eConstraint",
                        arg_message(second, "Incompatible with argument", first),
                    ));
                }
                exclude.insert(second);
            }
        }
    }
    for desc in DESCS {
        if given.contains(&desc.name) {
            continue;
        }
        if let Some(requester) = require.get(desc.name) {
            return Err(invalid(
                "eConstraint",
                arg_message(
                    desc.name,
                    "Must be specified, as it is required by argument",
                    requester,
                ),
            ));
        }
        if exclude.contains(desc.name) {
            continue;
        }
        if let Some(default) = desc.default {
            values.insert(desc.name, Some(default.to_string()));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(line: &[&str]) -> Vec<String> {
        line.iter().map(|word| word.to_string()).collect()
    }

    fn messages(line: &[&str]) -> (&'static str, Vec<String>) {
        match parse(&words(line)) {
            Err(ArgsError::Invalid { code, messages }) => (code, messages),
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
    }

    #[test]
    fn int8_conversion_follows_ncbi() {
        assert_eq!(string_to_int8("80"), Ok(80));
        assert_eq!(string_to_int8("+7"), Ok(7));
        assert_eq!(string_to_int8("-9223372036854775808"), Ok(i64::MIN));
        assert_eq!(
            string_to_int8("9223372036854775808").unwrap_err(),
            "Cannot convert string '9223372036854775808' to Int8, overflow (m_Pos = 18)"
        );
        assert_eq!(
            string_to_int8("12ab").unwrap_err(),
            "Cannot convert string '12ab' to Int8 (m_Pos = 2)"
        );
        assert_eq!(
            string_to_int8("-").unwrap_err(),
            "Cannot convert string '-' to Int8 (m_Pos = 1)"
        );
        assert_eq!(
            string_to_int8("").unwrap_err(),
            "Cannot convert string '' to Int8 (m_Pos = 0)"
        );
        assert!(matches!(
            integer_value("line_length", "0x10"),
            Err(ArgsError::Refused(_))
        ));
    }

    #[test]
    fn printable_string_escapes_as_ncbi() {
        assert_eq!(printable_string("a'b\"c\\d"), "a\\'b\\\"c\\\\d");
        assert_eq!(printable_string("x\ty"), "x\\ty");
        assert_eq!(printable_string("a??b?c"), "a\\?\\?b?c");
        assert_eq!(printable_string("\u{1}9"), "\\19");
        assert_eq!(printable_string("\u{1}1"), "\\0011");
        assert_eq!(printable_string("é"), "\\303\\251");
    }

    #[test]
    fn the_dependency_order_is_that_of_ncbi_s_multimap() {
        let deps = dependencies();
        assert!(deps.windows(2).all(|pair| pair[0].0 <= pair[1].0));
        // Among the entries of `entry`, the reverse of `entry_batch excludes entry` comes first.
        let first_of_entry = deps.iter().find(|dep| dep.0 == "entry").unwrap();
        assert_eq!(*first_of_entry, ("entry", Dep::Excludes, "entry_batch"));
    }
}
