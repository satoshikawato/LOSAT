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
    /// `eBoolean`: `CArg_Boolean` (`NStr::StringToBool`).
    Boolean,
    Integer,
    /// `eDouble`: `CArg_Double` (`NStr::StringToDouble` with `fDecimalPosixOrLocal`).
    Double,
    InputFile,
    OutputFile,
}

impl Kind {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2387-2400
    /// `CArgDescriptions::GetTypeName` of a key's type.
    fn type_name(self) -> &'static str {
        match self {
            Kind::Flag => "Boolean",
            Kind::Text => "String",
            Kind::Boolean => "Boolean",
            Kind::Integer => "Integer",
            Kind::Double => "Real",
            Kind::InputFile => "File_In",
            Kind::OutputFile => "File_Out",
        }
    }
}

/// The converted value of an `eBoolean`, `eInteger` or `eDouble` key.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Number {
    Boolean(bool),
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
    /// `CArgAllowValuesBetween(min, max, true)`: from `min` to `max`, both included.
    Between(i64, i64),
    /// `CArgAllow_Integers(min, max)`: an integer from `min` to `max` (`CArgAllow_Int8s::Verify`
    /// reads the value with `s_StringToInt8`, as `CArg_Integer` does).
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4686-4695
    Integers(i64, i64),
    /// `CArgAllow_Doubles(min, max)`: a number from `min` to `max`, and its usage text
    /// (`NStr::DoubleToString` of each end: `0..100`).
    Doubles(f64, f64, &'static str),
}

impl Constraint {
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:161-184
    // ```c++
    //     virtual bool Verify(const string& value) const {
    //         double val = NStr::StringToDouble(value);
    //         bool retval = false;
    //         if ( !m_Inclusive ) {
    //             retval = (val > m_MinValue && val < m_MaxValue);
    //         } else {
    //             retval = (val >= m_MinValue && val <= m_MaxValue);
    //         }
    //         return retval;
    //     }
    //
    //     /// Overloaded method from CArgAllow
    //     virtual string GetUsage(void) const {
    //         string retval;
    //         if ( !m_Inclusive ) {
    //             retval = "(>" + NStr::DoubleToString(m_MinValue) + " and <"
    //                 + NStr::DoubleToString(m_MaxValue) + ")";
    //         } else {
    //             retval = "(>=" + NStr::DoubleToString(m_MinValue) + " and =<"
    //                 + NStr::DoubleToString(m_MaxValue) + ")";
    //         }
    //         return retval;
    //     }
    // ```
    // (`NStr::DoubleToString` of a whole number writes it without a decimal point.)
    fn usage(self) -> String {
        match self {
            Constraint::Doubles(_, _, usage) => usage.to_string(),
            Constraint::Strings(values) => values
                .iter()
                .map(|value| format!("`{value}'"))
                .collect::<Vec<_>>()
                .join(", "),
            Constraint::AtLeast(min) => format!(">={min}"),
            Constraint::Between(min, max) => format!("(>={min} and =<{max})"),
            Constraint::Integers(min, max) => integers_usage(min, max),
        }
    }

    /// Whether the constraint compares the value as a number (`NStr::StringToDouble`).
    fn numeric(self) -> bool {
        matches!(self, Constraint::AtLeast(_) | Constraint::Between(..))
    }
}

/// `CArgAllow_Integers::GetUsage` of one range: the number when `min == max`, a bound when the
/// other is `int`'s limit, else `min..max` (`CArgAllow_Int8s::GetUsage`).
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4778-4794
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4698-4724
fn integers_usage(min: i64, max: i64) -> String {
    let (int_min, int_max) = (i64::from(i32::MIN), i64::from(i32::MAX));
    if min == max {
        min.to_string()
    } else if min == int_min && max != int_max {
        format!("less or equal to {max}")
    } else if min != int_min && max == int_max {
        format!("greater or equal to {min}")
    } else if min == int_min && max == int_max {
        String::new()
    } else {
        format!("{min}..{max}")
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

#[cfg(test)]
impl Desc {
    /// The argument's name.
    pub(crate) fn name(&self) -> &'static str {
        self.name
    }

    /// The usage of its constraint (`CArgAllow::GetUsage`), as `-help` prints it.
    pub(crate) fn constraint_usage(&self) -> Option<String> {
        self.constraint.map(Constraint::usage)
    }
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
/// The arguments of every `CArgDescriptions` before the program's own: the four help flags of
/// its constructor (`s_AutoHelp`, `s_AutoHelpFull`, `s_AutoHelpShowAll`, `s_AutoHelpXml`).
/// `HideStdArgs` takes effect later (`x_AddDefaultArgs`, after a help or an argument error,
/// before the usage is printed): `-help-full` and `-xmlhelp` are known to the reader (probes
/// C05, C06 of line BDC), and `-version` (which NCBI's application takes before the reader,
/// `ncbi_environment::preparse`) is not.
const HELP_FLAGS: &[Desc] = &[flag("h"), flag("help"), flag("help-full"), flag("xmlhelp")];

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1213-1218
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2229-2242
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
    /// `-xmlhelp` (`CArgHelpException::eHelpXml`): the description in XML on stdout, exit status
    /// 0 (`ncbi_xmlhelp`).
    XmlHelp,
    /// A `CArgException` of the reader.
    Invalid(ArgException),
    /// A value that NCBI's reader takes in a way LOSAT does not reproduce: LOSAT's refusal.
    Refused(String),
}

fn invalid(code: &'static str, message: String) -> ArgsError {
    ArgsError::Invalid(ArgException::new(code, message))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
pub(crate) fn arg_message(name: &str, what: &str, attr: &str) -> String {
    let name = if name.is_empty() { "...." } else { name };
    if attr.is_empty() {
        format!("Argument \"{name}\". {what}")
    } else {
        format!("Argument \"{name}\". {what}:  `{attr}'")
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:207-209, 223-224
/// The exception of a stream asked of the argument `name` that has no value.
fn no_value(name: &str) -> ArgException {
    ArgException::new("eNoValue", no_value_message(name))
}

/// The message of `THROW_CArg_NoValue` (quoted above) for the argument `name`: the value of an
/// argument that has none asked for (`AsString()` too, ncbiargs.cpp:213), which a program's `Run`
/// reports (blast_formatter's `-rid`).
pub(crate) fn no_value_message(name: &str) -> String {
    arg_message(name, "The argument has no value", "")
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
    /// `args[name]` (and `AsBoolean`) of the flag `name`: whether its word was given.
    pub(crate) fn flag(&self, name: &str) -> bool {
        self.given(name)
    }

    /// `args[name].AsString()` inside `Run`: the value of the key `name`, or NCBI's exception
    /// of an argument without a value (`CArg_NoValue`), which the application reports as one of
    /// the reader's (`Diag::run_arg_exception`).
    pub(crate) fn as_string(&self, name: &str) -> Result<&str, ArgException> {
        self.value(name).ok_or_else(|| no_value(name))
    }

    /// `AsBoolean()` of the `eBoolean` key `name`: given, or its default.
    pub(crate) fn boolean(&self, name: &str) -> Option<bool> {
        match self.numbers.get(name) {
            Some(&Number::Boolean(value)) => Some(value),
            _ => None,
        }
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
/// One option word `arg1` with the name `name` (and the next word `arg2`, the value of a key
/// written apart).
fn create_arg(
    tool: &Tool,
    arg1: &str,
    name: &str,
    arg2: Option<&String>,
) -> Result<Created, ArgsError> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:921-929
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
        match arg1 {
            "-h" => return Err(ArgsError::Usage),
            // `-help-full` adds the hidden arguments, of which LOSAT's programs have none: the
            // text is that of `-help` (probe C05 of line BDC).
            "-help" | "-help-full" => return Err(ArgsError::Help),
            // The help in XML (`ncbi_xmlhelp`).
            "-xmlhelp" => return Err(ArgsError::XmlHelp),
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
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1786-1789
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
/// `CArgDesc::ProcessArgument` (and `ProcessDefault`): the conversion of `value` to the key's
/// type, then its constraint; the number of an `eInteger` or `eDouble` key.
fn process_value(tool: &Tool, desc: &Desc, value: &str) -> Result<Option<Number>, ArgsError> {
    let number = match desc.kind {
        Kind::Boolean => Some(boolean_value(desc, value)?),
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
            Some(Number::Boolean(_)) | None => false,
        },
        // `s_StringToInt8(value)` within the range: the integer `integer_value` converted
        // (hexadecimal included).
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4686-4695
        Some(Constraint::Integers(min, max)) => match number {
            Some(Number::Integer(number)) => (min..=max).contains(&i64::from(number)),
            _ => false,
        },
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
        // (the conversion of `double_value`; a NaN is in no range.)
        Some(Constraint::Doubles(min, max, _)) => match number {
            Some(Number::Double(number)) => min <= number && number <= max,
            _ => false,
        },
        Some(Constraint::Between(min, max)) => match number {
            Some(Number::Integer(number)) => (min..=max).contains(&i64::from(number)),
            Some(Number::Double(number)) => number >= min as f64 && number <= max as f64,
            _ => false,
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
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:351-360, 383-390
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
        Err(_) if hexadecimal && desc.constraint.is_some_and(Constraint::numeric) => {
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
                String::from_utf8_lossy(&error.text()).into_owned(),
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

/// The value of an `eBoolean` key (`CArg_Boolean`): `NStr::StringToBool`, whose exception
/// (`String cannot be converted to bool`, position 0) is the predecessor of the reader's.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:489-498
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2817-2839
fn boolean_value(desc: &Desc, value: &str) -> Result<Number, ArgsError> {
    match crate::ncbi_str::string_to_bool(value) {
        Some(value) => Ok(Number::Boolean(value)),
        None => Err(ArgsError::Invalid(ArgException {
            code: "eConvert",
            messages: vec![
                crate::wording::losat_text("String cannot be converted to bool (m_Pos = 0)")
                    .into_owned(),
                arg_message(desc.name, "Argument cannot be converted", value),
            ],
        })),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
// As for `CArg_Int8`, the conversion's exception is the predecessor (its text is the first
// message: `Cannot convert string 'abc' to double (m_Pos = 0)`).
/// The number of an `eDouble` key (`CArg_Double`).
fn double_value(desc: &Desc, value: &str) -> Result<Number, ArgsError> {
    match crate::ncbi_str::string_to_double_posix_or_local(value.as_bytes()) {
        Ok(number) => Ok(Number::Double(number)),
        Err(error) => Err(ArgsError::Invalid(ArgException {
            code: "eConvert",
            messages: vec![
                String::from_utf8_lossy(&error.text()).into_owned(),
                arg_message(desc.name, "Argument cannot be converted", value),
            ],
        })),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3162-3296
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
            // (`x_PostCheck`, cited above, asks each argument that is not given for its
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
