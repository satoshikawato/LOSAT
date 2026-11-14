//! NCBI's argument reader (`CArgDescriptions`, c++/src/corelib/ncbiargs.cpp) as an application of
//! `ncbi_app` reads its command line: the program gives the descriptions of its `Init`
//! (`Tool::descs`, `Tool::dependencies`), this module adds those that every `CArgDescriptions`
//! has (the help flags, `-logfile`, `-conffile`), reads the words, and refuses every command line
//! that NCBI's reader refuses with NCBI's text. The values are read as NCBI's `CArgs` gives them
//! (`Args`): the text of a key, a flag, the integer of an `eInteger` key, the number of an
//! `eDouble` key, and the streams of the file keys (`files::open_input_file`,
//! `Args::output_file`, module `files`).

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

use super::files::{self, OutputFile};
use super::Tool;

/// What an argument holds (`CArgDescriptions::EType` of a key, or a flag).
///
/// The file types (`eInputFile`, `eOutputFile`) are read as text: LOSAT's programs describe
/// them without `fPreOpen`, so the reader opens nothing (ncbiargs.cpp:686-693, 744-752); the
/// program opens them where its `Run` asks for the stream (`files::open_input_file`,
/// `Args::output_file`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Flag,
    Text,
    Integer,
    /// `eDouble`: `CArg_Double` (`NStr::StringToDouble` with `fDecimalPosixOrLocal`).
    Double,
    InputFile,
    OutputFile,
}

impl Kind {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2387-2400
    // ```c++
    //     static const char* s_TypeName[k_EType_Size] = {
    //         "String",
    //         "Boolean",
    //         "Int8",
    //         "Integer",
    //         "IntId",
    //         "Real",
    //         "File_In",
    //         "File_Out",
    // ```
    /// `CArgDescriptions::GetTypeName` of a key's type.
    fn type_name(self) -> &'static str {
        match self {
            Kind::Flag => "Boolean",
            Kind::Text => "String",
            Kind::Integer => "Integer",
            Kind::Double => "Real",
            Kind::InputFile => "File_In",
            Kind::OutputFile => "File_Out",
        }
    }
}

/// The number of an `eInteger` or `eDouble` key.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Number {
    Integer(i32),
    Double(f64),
}

/// A constraint on a key's value (`CArgAllow`).
#[derive(Debug, Clone, Copy)]
pub(crate) enum Constraint {
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
pub(crate) struct Desc {
    name: &'static str,
    kind: Kind,
    default: Option<&'static str>,
    constraint: Option<Constraint>,
    /// `AddKey`: the key must be given (`CArgDescMandatory`).
    mandatory: bool,
}

/// `AddFlag(name, comment, true)`: a flag that its word sets (every flag of LOSAT's programs).
pub(crate) const fn flag(name: &'static str) -> Desc {
    Desc {
        name,
        kind: Kind::Flag,
        default: None,
        constraint: None,
        mandatory: false,
    }
}

/// `AddKey`, `AddOptionalKey` or `AddDefaultKey(name, ..., kind, default)`, and the constraint
/// that `SetConstraint` gives it.
pub(crate) const fn key(
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
        mandatory: false,
    }
}

/// `AddKey(name, ..., kind)`: a key without a default that the command line must give.
pub(crate) const fn mandatory_key(name: &'static str, kind: Kind) -> Desc {
    Desc {
        name,
        kind,
        default: None,
        constraint: None,
        mandatory: true,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2207-2219
// ```c++
//     if ( m_AutoHelp ) {
//         AddFlag(s_AutoHelp,
//                 "Print USAGE and DESCRIPTION;  ignore all other parameters");
//     }
//     AddFlag(s_AutoHelpFull,
//             "Print USAGE, DESCRIPTION and ARGUMENTS;"
//             " ignore all other parameters");
//     AddFlag(s_AutoHelpShowAll,
//             "Print USAGE, DESCRIPTION and ARGUMENTS, including hidden ones;"
//             " ignore all other parameters");
//     AddFlag(s_AutoHelpXml,
//             "Print USAGE, DESCRIPTION and ARGUMENTS in XML format;"
//             " ignore all other parameters");
// ```
/// The arguments of every `CArgDescriptions` before the program's own: the four help flags of
/// its constructor (`s_AutoHelp`, `s_AutoHelpFull`, `s_AutoHelpShowAll`, `s_AutoHelpXml`).
/// `HideStdArgs` takes effect later (`x_AddDefaultArgs`, after a help or an argument error,
/// before the usage is printed): `-help-full` and `-xmlhelp` are known to the reader (probes
/// C05, C06 of line BDC), and `-version` (which NCBI's application takes before the reader,
/// `ncbi_environment::preparse`) is not.
const HELP_FLAGS: &[Desc] = &[flag("h"), flag("help"), flag("help-full"), flag("xmlhelp")];

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1213-1218
// ```c++
//         if ( !m_DisableArgDesc ) {
//             for(CArgDescriptions* desc : m_ArgDesc->GetAllDescriptions()) {
//                 desc->AddDefaultFileArguments(m_DefaultConfig);
//             }
//         }
//         m_Args.reset(arg_desc->CreateArgs(GetArguments()));
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2229-2242
// ```c++
// void CArgDescriptions::AddDefaultFileArguments(const string& default_config)
// {
//     if (!Exist(s_ArgLogFile + 1) ) {
//         AddOptionalKey
//             (s_ArgLogFile+1, "File_Name",
//                 "File to which the program log should be redirected",
//                 CArgDescriptions::eOutputFile);
//     }
//     if (!Exist(s_ArgCfgFile + 1) ) {
//         if (default_config.empty()) {
//             AddOptionalKey
//                 (s_ArgCfgFile + 1, "File_Name",
//                     "Program's configuration (registry) data file",
//                     CArgDescriptions::eInputFile);
// ```
/// The arguments of every application after the program's own (`AddDefaultFileArguments`):
/// the application opens `-logfile` itself (`ncbi_environment::diag_target`), and LOSAT refuses
/// `-conffile` before the reader (`ncbi_environment::preparse`).
const FILE_ARGUMENTS: &[Desc] = &[
    key("logfile", Kind::OutputFile, None, None),
    key("conffile", Kind::InputFile, None, None),
];

/// The arguments in the order NCBI describes them (`m_KeyFlagArgs`) when it reads a command
/// line: the help flags, the program's (`tool.descs`), then the file arguments.
fn descs(tool: &Tool) -> impl Iterator<Item = &'static Desc> {
    HELP_FLAGS
        .iter()
        .chain(tool.descs)
        .chain(FILE_ARGUMENTS.iter())
}

/// A dependency between two arguments (`CArgDescriptions::EDependency`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Dep {
    Requires,
    Excludes,
}

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
/// The dependencies of `tool` (its `SetDependency` calls, `Tool::dependencies`) as NCBI's
/// `multimap` holds them: each exclusion in both directions, ordered by the first argument's
/// name, and in the order of insertion among equal names. Which of several conflicts NCBI
/// reports depends on this order.
pub(crate) fn dependencies(tool: &Tool) -> Vec<(&'static str, Dep, &'static str)> {
    let mut deps = Vec::new();
    for (first, dep, second) in (tool.dependencies)() {
        deps.push((first, dep, second));
        if dep == Dep::Excludes {
            deps.push((second, dep, first));
        }
    }
    // A stable sort keeps the insertion order among equal names, as `multimap::insert` does.
    deps.sort_by_key(|&(first, _, _)| first);
    deps
}

/// A `CArgException` with the code name `code` (`eInvalidArg`, ...): `messages` are the
/// messages of the exceptions it was thrown for (`NCBI_RETHROW`), oldest first, then its own.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct ArgException {
    pub(crate) code: &'static str,
    pub(crate) messages: Vec<String>,
}

impl ArgException {
    /// `NCBI_THROW(CArgException, code, message)`.
    pub(crate) fn new(code: &'static str, message: String) -> Self {
        ArgException {
            code,
            messages: vec![message],
        }
    }

    /// `GetMsg()`: the message of the exception itself.
    pub(crate) fn message(&self) -> &str {
        self.messages.last().map_or("", String::as_str)
    }
}

/// Why the arguments were not read.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum ArgsError {
    /// `-h` (`CArgHelpException::eHelp`): the short usage on stdout, exit status 0.
    Usage,
    /// `-help` (`CArgHelpException::eHelpFull`): the full description on stdout, exit status 0.
    Help,
    /// A `CArgException` of the reader.
    Invalid(ArgException),
    /// A value that NCBI's reader takes in a way LOSAT does not reproduce: LOSAT's refusal.
    Refused(String),
}

fn invalid(code: &'static str, message: String) -> ArgsError {
    ArgsError::Invalid(ArgException::new(code, message))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
// ```c++
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
pub(crate) fn arg_message(name: &str, what: &str, attr: &str) -> String {
    let name = if name.is_empty() { "...." } else { name };
    if attr.is_empty() {
        format!("Argument \"{name}\". {what}")
    } else {
        format!("Argument \"{name}\". {what}:  `{attr}'")
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:207-209, 223-224
// ```c++
// #define THROW_CArg_NoValue \
//     NCBI_THROW(CArgException,eNoValue, s_ArgExptMsg(GetName(), \
//         "The argument has no value", ""));
// ...
// CNcbiIstream& CArg_NoValue::AsInputFile (CArgValue::TFileFlags) const { THROW_CArg_NoValue; }
// CNcbiOstream& CArg_NoValue::AsOutputFile(CArgValue::TFileFlags) const { THROW_CArg_NoValue; }
// ```
/// The exception of a stream asked of the argument `name` that has no value.
fn no_value(name: &str) -> ArgException {
    ArgException::new(
        "eNoValue",
        arg_message(name, "The argument has no value", ""),
    )
}

/// The arguments of one command line: those given and the defaults NCBI adds.
#[derive(Debug)]
pub(crate) struct Args {
    /// The value of each argument that has one: the value given (`None` for a flag), else the
    /// default of a key that no given argument excludes.
    values: BTreeMap<&'static str, Option<String>>,
    /// The number of each `eInteger` and `eDouble` key that has a value (`CArg_Integer`,
    /// `CArg_Double`).
    numbers: BTreeMap<&'static str, Number>,
    /// The names given on the command line, in their order.
    given: Vec<&'static str>,
}

impl Args {
    /// Whether the argument `name` was given on the command line.
    pub(crate) fn given(&self, name: &str) -> bool {
        self.given.contains(&name)
    }

    /// The value of the key `name` as text (`AsString`): given, or its default.
    pub(crate) fn value(&self, name: &str) -> Option<&str> {
        self.values.get(name).and_then(|value| value.as_deref())
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:509-516
    // ```c++
    // CArg_Flag::CArg_Flag(const string& name, bool value)
    //     : CArg_Boolean(name, value)
    // {
    // }
    // bool CArg_Flag::HasValue(void) const
    // {
    //     return AsBoolean();
    // }
    // ```
    /// `args[name]` (and `AsBoolean`) of the flag `name`: whether its word was given.
    pub(crate) fn flag(&self, name: &str) -> bool {
        self.given(name)
    }

    /// `AsInteger()` of the `eInteger` key `name`: given, or its default.
    pub(crate) fn integer(&self, name: &str) -> Option<i32> {
        match self.numbers.get(name) {
            Some(&Number::Integer(number)) => Some(number),
            _ => None,
        }
    }

    /// `AsDouble()` of the `eDouble` key `name`: given, or its default.
    pub(crate) fn double(&self, name: &str) -> Option<f64> {
        match self.numbers.get(name) {
            Some(&Number::Double(number)) => Some(number),
            _ => None,
        }
    }

    /// `args[name].AsOutputFile()`: the stream of the `eOutputFile` key `name`, `stdout` for
    /// `-`, or NCBI's exception when it does not open. A `CArgException` thrown in `Run` is
    /// reported by the application as one of the reader's (`Diag::run_arg_exception`). The
    /// stream of an `eInputFile` key is `files::open_input_file` of its value, which the
    /// programs open where NCBI's `Run` reads it.
    pub(crate) fn output_file<'a>(
        &self,
        name: &str,
        stdout: &'a mut dyn Write,
    ) -> Result<OutputFile<'a>, ArgException> {
        match self.value(name) {
            Some(value) => files::open_output_file(name, value, stdout),
            None => Err(no_value(name)),
        }
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
/// application takes before the reader: `ncbi_environment::preparse`) as NCBI's
/// `CArgDescriptions::CreateArgs` reads `tool`'s: each word in its order, then the dependencies
/// and the defaults (`x_PostCheck`). LOSAT's programs have no positional argument, so a word that
/// is not an option, and every word after `--`, is one too many.
pub(crate) fn parse(tool: &Tool, words: &[String]) -> Result<Args, ArgsError> {
    let mut values: BTreeMap<&'static str, Option<String>> = BTreeMap::new();
    let mut numbers: BTreeMap<&'static str, Number> = BTreeMap::new();
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
        let Created {
            desc,
            value,
            number,
            used,
        } = create_arg(tool, arg1, name, arg2)?;
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
        if let Some(number) = number {
            numbers.insert(desc.name, number);
        }
        given.push(desc.name);
        at += if used { 2 } else { 1 };
    }
    post_check(tool, &mut values, &mut numbers, &given)?;
    Ok(Args {
        values,
        numbers,
        given,
    })
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
/// written apart).
fn create_arg(
    tool: &Tool,
    arg1: &str,
    name: &str,
    arg2: Option<&String>,
) -> Result<Created, ArgsError> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:921-929
    // ```c++
    // CArgDesc::CArgDesc(const string& name, const string& comment, CArgDescriptions::TFlags flags)
    //     : m_Name(name), m_Comment(comment), m_Flags(flags)
    // {
    //     if ( !CArgDescriptions::VerifyName(m_Name) ) {
    //         NCBI_THROW(CArgException,eInvalidArg,
    //             "Invalid argument name: " + m_Name);
    //     }
    // ```
    // The lookup of `name` (`x_Find`) makes a `CArgDesc_NameOnly` of it, so a name that is no
    // argument name fails there.
    if !is_argument_name(name) {
        return Err(invalid(
            "eInvalidArg",
            format!("Invalid argument name: {name}"),
        ));
    }
    let find = |name: &str| descs(tool).find(|desc| desc.name == name);
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2975-2989
    // ```c++
    //     size_t dash = name.rfind('-');
    //     if (it == m_Args.end() && dash != NPOS && dash != 0) {
    //         string test(name.substr(0, dash));
    //         string suffix(name.substr(dash+1));
    //         if (NStr::strcasecmp(suffix.c_str(), "file") == 0 ||
    //             NStr::strcasecmp(suffix.c_str(), "verbatim") == 0)
    //         {
    //             try {
    //                 it = x_Find(test);
    //             } catch (const CArgException&) {
    //                 it = m_Args.end();
    //             }
    //             if (it != m_Args.end()) {
    // ```
    // The confidential arguments read their value from a file (`-name-file`); LOSAT's programs
    // have none, but the lookup of the name before the dash stays: `-info-file` is `-info`, and
    // `-db-file NAME` is `-db NAME` (the separator is looked for after the whole name).
    let desc = find(name).or_else(|| {
        let dash = name.rfind('-').filter(|&dash| dash != 0)?;
        let suffix = &name[dash + 1..];
        if suffix.eq_ignore_ascii_case("file") || suffix.eq_ignore_ascii_case("verbatim") {
            find(&name[..dash])
        } else {
            None
        }
    });
    let Some(desc) = desc else {
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
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2840-2852
        // ```c++
        // void CArgDescriptions::x_CheckAutoHelp(const string& arg) const
        // {
        //     if (arg.compare(string("-") + s_AutoHelp) == 0) {
        //         if (m_AutoHelp) {
        //             NCBI_THROW(CArgHelpException,eHelp,kEmptyStr);
        //         }
        //     } else if (arg.compare(string("-") + s_AutoHelpFull) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpFull,kEmptyStr);
        //     } else if (arg.compare(string("-") + s_AutoHelpXml) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpXml,kEmptyStr);
        //     } else if (arg.compare(string("-") + s_AutoHelpShowAll) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpShowAll,kEmptyStr);
        //     }
        // ```
        match arg1 {
            "-h" => return Err(ArgsError::Usage),
            // `-help-full` adds the hidden arguments, of which LOSAT's programs have none: the
            // text is that of `-help` (probe C05 of line BDC).
            "-help" | "-help-full" => return Err(ArgsError::Help),
            // The help in XML.
            "-xmlhelp" => {
                return Err(ArgsError::Refused(tool.unsupported_option(arg1)));
            }
            _ => {}
        }
        if eq_separator {
            return Err(invalid("eInvalidArg", format!("Invalid argument: {arg1}")));
        }
        return Ok(Created {
            desc,
            value: None,
            number: None,
            used: false,
        });
    }
    let (value, used) = if eq_separator {
        (arg1[name.len() + 2..].to_string(), false)
    } else {
        match arg2 {
            Some(value) => (value.clone(), true),
            None => return Err(invalid("eNoArg", arg_message(arg1, "Value is missing", ""))),
        }
    };
    let number = process_value(tool, desc, &value)?;
    Ok(Created {
        desc,
        value: Some(value),
        number,
        used,
    })
}

/// The argument that one option word gives (`create_arg`).
struct Created {
    /// Its description.
    desc: &'static Desc,
    /// Its value (`None` for a flag).
    value: Option<String>,
    /// The number of an `eInteger` or `eDouble` key.
    number: Option<Number>,
    /// Whether the next word was its value.
    used: bool,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3338-3365
// ```c++
// bool CArgDescriptions::VerifyName(const string& name, bool extended)
// {
//     if ( name.empty() )
//         return true;
//
//     string::const_iterator it = name.begin();
//     if (extended  &&  *it == '#') {
//         ...
//     } else {
//         if (name[0] == '-') {
//             // Prohibit names like '-' or '--foo'.
//             // The second char must be present and may not be '-'.
//             if (name.length() == 1  ||  name[1] == '-') {
//                 return false;
//             }
//         }
//         for ( ;  it != name.end();  ++it) {
//             if ( !s_IsArgNameChar((unsigned char)(*it)) )
//                 return false;
//         }
//     }
//
//     return true;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1786-1789
// ```c++
// inline bool s_IsArgNameChar(char c)
// {
//     return isalnum(c)  ||  c == '_'  ||  c == '-';
// }
// ```
/// `CArgDescriptions::VerifyName(name)`: whether `name` can be an argument's name (ASCII letters
/// and digits, `_` and `-`, not `-` alone or starting with `--`; the empty name can).
fn is_argument_name(name: &str) -> bool {
    let bytes = name.as_bytes();
    if bytes.first() == Some(&b'-') && matches!(bytes.get(1), None | Some(b'-')) {
        return false;
    }
    bytes
        .iter()
        .all(|&byte| byte.is_ascii_alphanumeric() || byte == b'_' || byte == b'-')
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
/// `CArgDesc::ProcessArgument` (and `ProcessDefault`): the conversion of `value` to the key's
/// type, then its constraint; the number of an `eInteger` or `eDouble` key.
fn process_value(tool: &Tool, desc: &Desc, value: &str) -> Result<Option<Number>, ArgsError> {
    let number = match desc.kind {
        Kind::Integer => Some(Number::Integer(integer_value(tool, desc, value)?)),
        Kind::Double => Some(double_value(desc, value)?),
        Kind::Text | Kind::Flag | Kind::InputFile | Kind::OutputFile => None,
    };
    let allowed = match desc.constraint {
        None => true,
        Some(Constraint::Strings(values)) => values.contains(&value),
        // `NStr::StringToDouble(value) >= m_MinValue`: the decimal integer `integer_value`
        // accepted is the same number as a double.
        Some(Constraint::AtLeast(min)) => match number {
            Some(Number::Integer(number)) => i64::from(number) >= min,
            Some(Number::Double(number)) => number >= min as f64,
            None => false,
        },
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
        _ => Ok(number),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:115-129
// ```c++
// // Allow autodetection among decimal and hex, but NOT octal, in case
// // anyone's been relying on leading zeros being meaningless.
// inline
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
// The exception of the conversion is the predecessor of the reader's: its text (with the
// position, `NcbiException::text`) is the first message (`Cannot convert string 'x' to Int8
// (m_Pos = 0)`). A hexadecimal value is converted in base 16, whose exception is then the one
// thrown.
/// The value of an `eInteger` key (`CArg_Integer`). A hexadecimal value (`0x...`) of a key with
/// a constraint that compares numbers is LOSAT's refusal: NCBI's constraint check reads it as a
/// double (`NStr::StringToDouble`), which LOSAT does not reproduce.
fn integer_value(tool: &Tool, desc: &Desc, value: &str) -> Result<i32, ArgsError> {
    use crate::ncbi_str::{string_to_int8, string_to_int8_base};
    let bytes = value.as_bytes();
    let hexadecimal = bytes
        .get(..2)
        .is_some_and(|head| head.eq_ignore_ascii_case(b"0x"));
    let converted = match string_to_int8(bytes, 0) {
        Ok(number) => Ok(number),
        Err(_) if hexadecimal && matches!(desc.constraint, Some(Constraint::AtLeast(_))) => {
            return Err(ArgsError::Refused(format!(
                "the hexadecimal value '{value}' of -{} is not supported by LOSAT's {}",
                desc.name, tool.refusal_name
            )));
        }
        Err(_) if hexadecimal => string_to_int8_base(bytes, 0, 16),
        Err(error) => Err(error),
    };
    let number = converted.map_err(|error| {
        ArgsError::Invalid(ArgException {
            code: "eConvert",
            messages: vec![
                error.text(),
                arg_message(desc.name, "Argument cannot be converted", value),
            ],
        })
    })?;
    i32::try_from(number).map_err(|_| {
        invalid(
            "eConvert",
            arg_message(desc.name, "Integer value is out of range", value),
        )
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
// ```c++
// inline CArg_Double::CArg_Double(const string& name, const string& value)
//     : CArg_String(name, value)
// {
//     try {
//         m_Double = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
//     } catch (const CException& e) {
//         NCBI_RETHROW(e,CArgException,eConvert,
//             s_ArgExptMsg(GetName(),"Argument cannot be converted",value));
//     }
// }
// ```
// As for `CArg_Int8`, the conversion's exception is the predecessor (its text is the first
// message: `Cannot convert string 'abc' to double (m_Pos = 0)`).
/// The number of an `eDouble` key (`CArg_Double`).
fn double_value(desc: &Desc, value: &str) -> Result<Number, ArgsError> {
    match crate::ncbi_str::string_to_double_posix_or_local(value.as_bytes()) {
        Ok(number) => Ok(Number::Double(number)),
        Err(error) => Err(ArgsError::Invalid(ArgException {
            code: "eConvert",
            messages: vec![
                error.text(),
                arg_message(desc.name, "Argument cannot be converted", value),
            ],
        })),
    }
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
    tool: &Tool,
    values: &mut BTreeMap<&'static str, Option<String>>,
    numbers: &mut BTreeMap<&'static str, Number>,
    given: &[&'static str],
) -> Result<(), ArgsError> {
    let mut require: BTreeMap<&str, &str> = BTreeMap::new();
    let mut exclude: BTreeSet<&str> = BTreeSet::new();
    for (first, dep, second) in dependencies(tool) {
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
    for desc in descs(tool) {
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
            if let Some(number) = process_value(tool, desc, default)? {
                numbers.insert(desc.name, number);
            }
            values.insert(desc.name, Some(default.to_string()));
        } else if desc.mandatory {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1158-1173,1270-1275
            // ```c++
            // string CArgDescMandatory::GetUsageCommentAttr(void) const
            // {
            //     CArgDescriptions::EType type = GetType();
            //     // Print type name
            //     string str = CArgDescriptions::GetTypeName(type);
            //     ...
            //     // Print constraint info, if any
            //     string constr = GetUsageConstraint();
            //     if ( !constr.empty() ) {
            //         str += ", ";
            //         str += constr;
            //     }
            //     return str;
            // }
            // ...
            // CArgValue* CArgDescMandatory::ProcessDefault(void) const
            // {
            //     NCBI_THROW(CArgException, eNoArg,
            //                s_ArgExptMsg(GetName(), "Mandatory value is missing",
            //                             GetUsageCommentAttr()));
            // }
            // ```
            // (`x_PostCheck`, quoted above, asks each argument that is not given for its
            // default in the order of the descriptions.)
            let mut attr = desc.kind.type_name().to_string();
            if let Some(constraint) = desc.constraint {
                attr += ", ";
                attr += &constraint.usage();
            }
            return Err(invalid(
                "eNoArg",
                arg_message(desc.name, "Mandatory value is missing", &attr),
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod args_tests;
