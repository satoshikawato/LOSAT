//! NCBI's argument reader (`CArgDescriptions`, c++/src/corelib/ncbiargs.cpp) as an application of
//! `ncbi_app` reads its command line: the descriptions of the program's `Init` are those of its
//! NCBI `-xmlhelp` XML (`Arguments`, made from `ncbi_xmlhelp/<program>.xml` when LOSAT is built:
//! `LOSAT/build/ncbi_args.rs`; line RG DECISIONS K4i-1), this module adds those that every
//! `CArgDescriptions` has (the help flags, `-logfile`, `-conffile`), reads the words, and refuses
//! every command line that NCBI's reader refuses with NCBI's text: NCBI's conversion of each type
//! (`Kind`) and NCBI's constraints (`Constraint`). NCBI reads the bytes of each word, and so does
//! this reader (`parse_in`): an error quotes the bytes of the word, and a value keeps its bytes. The
//! values are read as NCBI's `CArgs` gives them (`Args`): the text of a key, a flag, the integer of
//! an `eInteger` key, the number of an `eDouble` key, and the streams of the file keys
//! (`files::open_input_file`, `Args::output_file`, module `files`).

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::{OsStr, OsString};
use std::io::Write;

use super::files::{self, OutputFile};
use super::Tool;
use crate::report::run_files::OpenFiles;

#[path = "constraint.rs"]
mod constraint;

pub(crate) use constraint::{is_genetic_code, string_to_int8, Constraint};
// The file name rule of LOSAT Web ABI v1's own output path parser (`web_v1_argv`).
#[cfg(any(losat_web_abi_v1, test))]
pub(crate) use constraint::file_name_length;

/// What an argument holds (`CArgDescriptions::EType` of a key, or a flag).
///
/// The file types (`eInputFile`, `eOutputFile`) are read as text: NCBI's programs describe
/// them without `fPreOpen` (their XML has no `<preOpen/>`), so the reader opens nothing
/// (ncbiargs.cpp:686-693, 744-752); the program opens them where its `Run` asks for the stream
/// (`files::open_input_file`, `Args::output_file`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Flag,
    /// `eString`.
    Text,
    /// `eBoolean`: `CArg_Boolean` (`NStr::StringToBool`).
    Boolean,
    /// `eInt8`: `CArg_Int8` (`s_StringToInt8`).
    Int8,
    /// `eInteger`: `CArg_Integer` (`s_StringToInt8`, then the range of `int`).
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
            Kind::Flag | Kind::Boolean => "Boolean",
            Kind::Text => "String",
            Kind::Int8 => "Int8",
            Kind::Integer => "Integer",
            Kind::Double => "Real",
            Kind::InputFile => "File_In",
            Kind::OutputFile => "File_Out",
        }
    }
}

/// The converted value of an `eBoolean`, `eInt8`, `eInteger` or `eDouble` key.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Number {
    Boolean(bool),
    Int8(i64),
    Integer(i32),
    Double(f64),
}

/// A text of the tables (a name, a default, an allowed value): `len` bytes of `TEXT` from `at`.
/// The tables hold texts and lists as offsets, not as `&str` and slices, so that no address in
/// them is relocated when LOSAT is loaded (`build/ncbi_args.rs`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Text {
    at: u16,
    len: u16,
}

impl Text {
    /// The text.
    fn get(self) -> &'static str {
        let at = usize::from(self.at);
        &TEXT[at..at + usize::from(self.len)]
    }

    /// Whether the text is `value`.
    fn is(self, value: &str) -> bool {
        usize::from(self.len) == value.len() && self.get() == value
    }
}

/// A list of the tables: `count` rows of an array from `first`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Span {
    first: u16,
    count: u16,
}

impl Span {
    /// The rows of `table` in the span.
    fn of<T>(self, table: &'static [T]) -> &'static [T] {
        let first = usize::from(self.first);
        &table[first..first + usize::from(self.count)]
    }
}

/// One argument description (`CArgDesc`), as the program's XML gives it (`<key>` or `<flag>`).
#[derive(Debug, PartialEq)]
pub(crate) struct Desc {
    name: Text,
    kind: Kind,
    /// `<default>`: the value of a key of `AddDefaultKey`.
    default: Option<Text>,
    /// `<constraint>`: what `SetConstraint` gave the key.
    constraint: Option<Constraint>,
    /// `AddKey` (a key without `optional="true"`): the key must be given (`CArgDescMandatory`).
    mandatory: bool,
    /// `<description>`: the argument's comment, which NCBI's `-help` prints (`described`).
    description: Text,
    /// `fPreOpen` of an `eOutputFile` key (blastdb_aliastool's `-gi_file_out`): the reader opens
    /// the file as it reads the word (`files::pre_open_output`).
    pre_open: bool,
}

impl Desc {
    /// The argument's name.
    fn name(&self) -> &'static str {
        self.name.get()
    }

    /// The default of a key of `AddDefaultKey`.
    fn default(&self) -> Option<&'static str> {
        self.default.map(Text::get)
    }
}

/// The argument descriptions of one program as its NCBI `-xmlhelp` XML gives them
/// (`CArgDescriptions::PrintUsageXml`, cited at `ncbi_xmlhelp::print`), without the arguments of
/// every `CArgDescriptions` (the help flags, `-logfile`, `-conffile`, `-version`, ...; the reader
/// adds the ones it knows: `HELP_FLAGS`, `FILE_ARGUMENTS`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Arguments {
    /// The keys, then the flags, each in the order of NCBI's `m_Args` (by name; rows of `DESCS`).
    /// The order in which NCBI's `Init` described them (`m_KeyFlagArgs`), which `x_PostCheck`
    /// follows, is the one of the usage (`usage_position`).
    descs: Span,
    /// `AddAlias` of an alias and a name (`CArgDesc_Alias`, `<alias>`): another name of the
    /// argument (rows of `ALIASES`).
    aliases: Span,
    /// NCBI's `m_Dependencies` (`<dependencies>`, rows of `DEPENDENCIES`): the `multimap` that
    /// `SetDependency` fills, an exclusion in both directions; its `eRequires` entries, then its
    /// `eExcludes` ones, each in the map's order (ordered by the first name, in the order of
    /// insertion among equal names). `x_PostCheck` (`post_check`) reports the first exclusion of
    /// the map that the command line breaks, and names, for each argument required and missing,
    /// the first argument of the map that requires it; neither depends on where the `eRequires`
    /// entries stand among the others.
    dependencies: Span,
    /// The arguments that LOSAT refuses (rows of `REFUSALS`): the name and how (`Refusal`).
    refusals: Span,
}

impl Arguments {
    /// The program's descriptions.
    fn descs(self) -> &'static [Desc] {
        self.descs.of(&DESCS)
    }

    /// The aliases: (alias, name).
    fn aliases(self) -> impl Iterator<Item = (&'static str, &'static str)> {
        let aliases = self.aliases.of(&ALIASES).iter();
        aliases.map(|&(alias, name)| (alias.get(), name.get()))
    }

    /// The entries of `m_Dependencies`: (first name, dependency, second name).
    pub(crate) fn dependencies(self) -> impl Iterator<Item = (&'static str, Dep, &'static str)> {
        let dependencies = self.dependencies.of(&DEPENDENCIES).iter();
        dependencies.map(|&(first, dep, second)| (first.get(), dep, second.get()))
    }

    /// The arguments that LOSAT refuses: (name, how).
    fn refusals(self) -> impl Iterator<Item = (&'static str, Refusal)> {
        let refusals = self.refusals.of(&REFUSALS).iter();
        refusals.map(|&(name, refused)| (name.get(), refused.view()))
    }
}

/// A row of `REFUSALS` as the tables hold it (`Refusal`, its texts as spans).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refused {
    Read(Option<Text>),
    Run(Text),
}

impl Refused {
    fn view(self) -> Refusal {
        match self {
            Refused::Read(why) => Refusal::Read(why.map(Text::get)),
            Refused::Run(what) => Refusal::Run(what.get()),
        }
    }
}

/// How LOSAT refuses one of NCBI's arguments that it does not run (AGENTS.md rule 2: explicit
/// unsupported errors): the one list of them, the rows of `REFUSED` in `build/ncbi_args.rs` (line
/// RG DECISIONS K4i-2), which the readers' and the runs' refusals (`cli::reader::named_refusal`,
/// `cli::run_refusal`) and `-xmlhelp` (`cli::xml_help_text`) read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// By name, where the reader meets it: why, if the refusal says it after the program's name.
    Read(Option<&'static str>),
    /// Where NCBI starts the search, whatever the value (NCBI's reader reads it, so that NCBI's
    /// own errors about it and the other arguments come first): what the argument is.
    Run(&'static str),
}

/// How LOSAT refuses `program`'s argument `name`; `None` for an argument it runs, or one that
/// `program` does not have.
pub(crate) fn refusal(program: &str, name: &str) -> Option<Refusal> {
    let mut refusals = arguments(program)?.refusals();
    refusals
        .find(|&(refused, _)| refused == name)
        .map(|(_, refusal)| refusal)
}

/// Every argument that LOSAT refuses: (program, name, how).
#[cfg(test)]
pub(crate) fn refusals() -> Vec<(&'static str, &'static str, Refusal)> {
    let programs = PROGRAMS.iter();
    programs
        .flat_map(|&(program, arguments)| {
            let refusals = arguments.refusals();
            refusals.map(move |(name, refusal)| (program.get(), name, refusal))
        })
        .collect()
}

/// The genetic codes that `program`'s key `name` allows (`CArgAllowGeneticCodeInteger`, with
/// LOSAT's code 32 where PD-TLOSAN-LOCAL-GENCODE-32 adds it), in increasing order; `None` for a
/// key of no genetic code. What the library API's and LOSAT Web's descriptions of a program list
/// (`api::describe`, `web_adapter`).
pub fn genetic_codes(program: &str, name: &str) -> Option<Vec<u32>> {
    match find_desc(arguments(program)?, name)?.constraint? {
        Constraint::GeneticCode { with_32 } => Some(constraint::genetic_codes(with_32)),
        _ => None,
    }
}

/// The values that `program`'s key `name` allows where its constraint is a set of strings
/// (`CArgAllow_Strings`, `CArgAllowStringSet`), in the order of NCBI's set; `None` for a key of
/// another constraint or none. LOSAT Web's defaults are resolved for each of `-task`'s
/// (`web_adapter`).
pub fn string_values(program: &str, name: &str) -> Option<Vec<&'static str>> {
    match find_desc(arguments(program)?, name)?.constraint? {
        Constraint::Strings(values) | Constraint::StringSet(values) => {
            Some(constraint::strings(values).collect())
        }
        _ => None,
    }
}

/// The values that `program`'s key `name` allows where its constraint is a set (`CArgAllow_Strings`,
/// `CArgAllowStringSet`, `CArgAllowIntegerSet`), as LOSAT's descriptions of the program list them:
/// in the order of `DISPLAY_ORDERS` where it has a row (line RG DECISIONS K4i-2b-1), else of NCBI's
/// set; `None` for a key of another constraint or none, and for a flag. The library API's and LOSAT
/// Web's descriptions (`api::describe`, `web_adapter`) and the check of a task that a caller set
/// (`value_parsers::TaskParser::check`) read it.
pub fn choices(program: &str, name: &str) -> Option<Vec<String>> {
    let display = DISPLAY_ORDERS
        .iter()
        .find(|(of, key, _)| of.is(program) && key.is(name));
    if let Some((_, _, values)) = display {
        return Some(constraint::strings(*values).map(str::to_string).collect());
    }
    match find_desc(arguments(program)?, name)?.constraint? {
        Constraint::Strings(values) | Constraint::StringSet(values) => {
            Some(constraint::strings(values).map(str::to_string).collect())
        }
        Constraint::IntegerSet(values) => Some(
            constraint::integers(values)
                .iter()
                .map(i32::to_string)
                .collect(),
        ),
        _ => None,
    }
}

/// One argument that LOSAT's descriptions of a program list (`described`), as NCBI's tables give
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Described {
    /// The name, without its dash.
    pub(crate) name: &'static str,
    pub(crate) kind: Kind,
    /// `<default>`: the default of a key of `AddDefaultKey`, which the library API's and LOSAT
    /// Web's descriptions give as the option's default (line RG DECISIONS K4i-2b-1 (3)).
    pub(crate) default: Option<&'static str>,
    /// The values of its set (`choices`); empty for a key of another constraint or none, and for
    /// a flag.
    pub(crate) choices: Vec<String>,
    /// `<description>`: NCBI's comment of the argument, which its `-help` prints (NCBI's words).
    pub(crate) description: &'static str,
}

// No NCBI counterpart: line RG DECISIONS K4i-3-1. NCBI's `-help` prints the usage line first, which
// names the arguments in the order in which the program described them (`usage_position`), then the
// arguments with their comments, by group.
/// The arguments of `tool`'s program that LOSAT's descriptions of it list (the library API's
/// `describe`, LOSAT Web's): every argument of its XML (the arguments of every NCBI application,
/// the help flags, `-logfile`, `-conffile`, `-version`, ..., are not in its tables) but those that
/// LOSAT refuses by name (`Refusal::Read`), in the order of the usage line of NCBI's `-help`
/// (`usage_position` in `Tool::help`; blastn: `-task`, `-db`, `-dbsize`, ...), not of the list of
/// arguments under it (blastn: `-query`, `-query_loc`, `-strand`, `-task`, ...).
pub(crate) fn described(tool: &Tool) -> impl Iterator<Item = Described> {
    let mut listed: Vec<(usize, Described)> = tool_arguments(tool)
        .descs()
        .iter()
        .filter(|desc| !matches!(refusal(tool.program, desc.name()), Some(Refusal::Read(_))))
        .map(|desc| {
            let name = desc.name();
            let described = Described {
                name,
                kind: desc.kind,
                default: desc.default(),
                choices: choices(tool.program, name).unwrap_or_default(),
                description: desc.description.get(),
            };
            (usage_position(tool.help, name), described)
        })
        .collect();
    listed.sort_by_key(|(position, _)| *position);
    listed.into_iter().map(|(_, described)| described)
}

// The tables (`TEXT`, `HELP_FLAGS`, `FILE_ARGUMENTS`, `DESCS`, `ALIASES`, `DEPENDENCIES`,
// `STRING_VALUES`, `INTEGER_VALUES`, `REFUSALS`, `DISPLAY_ORDERS`, `PROGRAMS`), written by
// `LOSAT/build.rs` (`build/ncbi_args.rs`) from the XML of every program LOSAT embeds. `HELP_FLAGS`
// are the arguments of every `CArgDescriptions` before the program's own (the help flags of its
// constructor), `FILE_ARGUMENTS` those of every application after them (`AddDefaultFileArguments`:
// `-logfile`, `-conffile`); the build script cites NCBI's lines.
include!(concat!(env!("OUT_DIR"), "/ncbi_args.rs"));

/// Whether `name` is one of the arguments that every NCBI application has after its own
/// (`FILE_ARGUMENTS`: `-logfile`, `-conffile`), which the application start reads.
pub(crate) fn is_file_argument(name: &str) -> bool {
    FILE_ARGUMENTS.iter().any(|desc| desc.name.is(name))
}

/// The argument descriptions of `program` (`PROGRAMS`), if LOSAT has its XML.
pub(crate) fn arguments(program: &str) -> Option<Arguments> {
    PROGRAMS
        .iter()
        .find(|(name, _)| name.is(program))
        .map(|&(_, arguments)| arguments)
}

/// The argument descriptions of `tool`'s program: every program on NCBI's reader has its XML
/// (`ncbi_xmlhelp`).
fn tool_arguments(tool: &Tool) -> Arguments {
    arguments(tool.program).expect("the -xmlhelp XML of a program on NCBI's reader")
}

/// The usage of the constraint of `program`'s argument `name` (`CArgAllow::GetUsage`), as its
/// `-help` prints it.
#[cfg(test)]
pub(crate) fn constraint_usage(program: &str, name: &str) -> Option<String> {
    find_desc(arguments(program)?, name)?
        .constraint
        .map(Constraint::usage)
}

/// NCBI's `m_Dependencies` of `program` (`Arguments::dependencies`).
#[cfg(test)]
pub(crate) fn dependencies(program: &str) -> Vec<(&'static str, Dep, &'static str)> {
    arguments(program).map_or_else(Vec::new, |arguments| arguments.dependencies().collect())
}

/// The arguments that NCBI's reader of `arguments` knows (its `m_Args`): the help flags, the
/// program's, then the file arguments.
fn descs(arguments: Arguments) -> impl Iterator<Item = &'static Desc> {
    HELP_FLAGS
        .iter()
        .chain(arguments.descs())
        .chain(FILE_ARGUMENTS.iter())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2718-2734
/// `x_Find`: the description of the argument `name`, or of the argument that `name` is an
/// alias of.
fn find_desc(arguments: Arguments, name: &str) -> Option<&'static Desc> {
    let name = arguments
        .aliases()
        .find(|(alias, _)| *alias == name)
        .map_or(name, |(_, aliased)| aliased);
    descs(arguments).find(|desc| desc.name.is(name))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2618-2629
/// A dependency between two arguments (`CArgDescriptions::EDependency`), an entry of NCBI's
/// `m_Dependencies` (`Arguments::dependencies`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Dep {
    Requires,
    Excludes,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3181-3204
/// The first walk of `x_PostCheck` over the dependencies `deps`, the entries of NCBI's
/// `m_Dependencies` in the order of its `multimap` (`Arguments::dependencies`), of the arguments
/// that `given` says were given: the first exclusion between two given arguments is the error
/// (`Err`, its message, bytes as `arg_message` gives them), else the arguments that the given
/// ones require, each with its first requester, and those they exclude: the walk of the reader's
/// `post_check`.
fn walk_dependencies<'a>(
    deps: impl IntoIterator<Item = (&'a str, Dep, &'a str)>,
    given: &dyn Fn(&str) -> bool,
) -> Result<DependencyWalk<'a>, Vec<u8>> {
    let mut walk = DependencyWalk::default();
    for (first, dep, second) in deps {
        if !given(first) {
            continue;
        }
        match dep {
            Dep::Requires => {
                walk.require.entry(second).or_insert(first);
            }
            Dep::Excludes => {
                if given(second) {
                    return Err(arg_message(
                        second,
                        "Incompatible with argument",
                        first.as_bytes(),
                    ));
                }
                walk.exclude.insert(second);
            }
        }
    }
    Ok(walk)
}

/// What [`walk_dependencies`] found: the second walk of `x_PostCheck` goes over the arguments
/// in the order of their descriptions and asks each one that was not given.
#[derive(Debug, Default)]
pub(crate) struct DependencyWalk<'a> {
    /// The arguments that given ones require, with the first argument that requires each.
    require: BTreeMap<&'a str, &'a str>,
    /// The arguments that given ones exclude (they take no default).
    exclude: BTreeSet<&'a str>,
}

impl DependencyWalk<'_> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3270-3276
    /// The error of the argument `name`, which was not given, when a given argument requires it.
    pub(crate) fn missing(&self, name: &str) -> Option<Vec<u8>> {
        let requester = self.require.get(name)?;
        Some(arg_message(
            name,
            "Must be specified, as it is required by argument",
            requester.as_bytes(),
        ))
    }

    /// Whether a given argument excludes the argument `name`.
    pub(crate) fn excludes(&self, name: &str) -> bool {
        self.exclude.contains(name)
    }
}

/// A `CArgException` with the code name `code` (`eInvalidArg`, ...): `messages` are the
/// messages of the exceptions it was thrown for (`NCBI_RETHROW`), oldest first, then its own.
/// A message is bytes: NCBI's quote the bytes of the words, which need not be UTF-8.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct ArgException {
    pub(crate) code: &'static str,
    pub(crate) messages: Vec<Vec<u8>>,
}

impl ArgException {
    /// `NCBI_THROW` of a `CArgException` with `code` and `message`.
    pub(crate) fn new(code: &'static str, message: impl Into<Vec<u8>>) -> Self {
        ArgException {
            code,
            messages: vec![message.into()],
        }
    }

    /// `GetMsg()`: the message of the exception itself.
    pub(crate) fn message(&self) -> &[u8] {
        self.messages.last().map_or(&[], Vec::as_slice)
    }

    /// The message as text, for a program's own message about the exception (a byte that is not
    /// UTF-8 replaced).
    pub(crate) fn message_text(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(self.message())
    }
}

/// Why the arguments were not read.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum ArgsError {
    /// `-h` (`CArgHelpException::eHelp`): the short usage on stdout, exit status 0.
    Usage,
    /// `-help` (`CArgHelpException::eHelpFull`): the full description on stdout, exit status 0.
    Help,
    /// `-help-full` (`CArgHelpException::eHelpShowAll`): the full description with the hidden
    /// arguments on stdout, exit status 0.
    HelpFull,
    /// `-xmlhelp` (`CArgHelpException::eHelpXml`): the description in XML on stdout, exit status
    /// 0 (`ncbi_xmlhelp`).
    XmlHelp,
    /// A `CArgException` of the reader.
    Invalid(ArgException),
}

fn invalid(code: &'static str, message: impl Into<Vec<u8>>) -> ArgsError {
    ArgsError::Invalid(ArgException::new(code, message))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
/// The message of an exception about the argument `name`: what is wrong, and the bytes `attr`
/// (a value of the command line), if any.
pub(crate) fn arg_message(name: &str, what: &str, attr: &[u8]) -> Vec<u8> {
    let name = if name.is_empty() { "...." } else { name };
    let mut message = format!("Argument \"{name}\". {what}").into_bytes();
    if !attr.is_empty() {
        message.extend_from_slice(b":  `");
        message.extend_from_slice(attr);
        message.push(b'\'');
    }
    message
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:207-209, 223-224
/// The exception of a stream or a value asked of the argument `name` that has no value
/// (`CArg_NoValue`: `args[name].AsString()` inside `Run` of a key without a value, which the
/// application reports as one of the reader's, `Diag::run_arg_exception`).
pub(crate) fn no_value(name: &str) -> ArgException {
    ArgException::new(
        "eNoValue",
        arg_message(name, "The argument has no value", b""),
    )
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:239-249
/// The exception of a value asked of the argument `name` that a given argument excludes
/// (`CArg_ExcludedValue`): blastdb_aliastool's default `-dbtype` beside `-seqid_file_out`.
fn excluded_value(name: &str) -> ArgException {
    ArgException::new(
        "eExcludedValue",
        arg_message(name, "The value is excluded by other arguments.", b""),
    )
}

/// A key's value as the command line gives it: its text, or its bytes where they are not UTF-8
/// (NCBI's argument holds the bytes, which the programs of `ncbi_app` read: `Args::bytes`). The
/// word is read as UTF-8 once, when the reader takes it.
#[derive(Debug, Clone)]
enum Value {
    Text(String),
    Bytes(OsString),
}

impl Value {
    /// The value of the word `word`.
    fn of(word: OsString) -> Self {
        match word.into_string() {
            Ok(text) => Value::Text(text),
            Err(bytes) => Value::Bytes(bytes),
        }
    }
}

/// The arguments of one command line: those given and the defaults NCBI adds.
#[derive(Debug, Clone)]
pub(crate) struct Args {
    /// The descriptions that the arguments were read with (the program's), which the check of a
    /// read of a value as text looks the key up in (`check_text_read`) and which the readers of
    /// the options of the programs on their own option types check their names against and take
    /// a key's table default from (`option_or_default`).
    arguments: Arguments,
    /// The value of each argument that has one: the value given (`None` for a flag), else the
    /// default of a key that no given argument excludes.
    values: BTreeMap<&'static str, Option<Value>>,
    /// The converted value of each `eBoolean`, `eInt8`, `eInteger` and `eDouble` key that has a
    /// value (`CArg_Boolean`, `CArg_Int8`, `CArg_Integer`, `CArg_Double`).
    numbers: BTreeMap<&'static str, Number>,
    /// The names given on the command line, in their order.
    given: Vec<&'static str>,
    /// The arguments not given that a given argument excludes (`CArg_ExcludedValue`), whose
    /// values `as_string` refuses with NCBI's exception.
    excluded: Vec<&'static str>,
    /// Where the files that the program opens by name for writing go (`output_file`): the file
    /// system, or LOSAT Web's files (`OpenFiles`), which the application gives the arguments.
    files: OpenFiles,
    /// The names that the options' readers read (`given_option`, `option_or_default`, `flag`),
    /// for the test that every argument of a program is read.
    #[cfg(test)]
    read: std::sync::Arc<std::sync::Mutex<BTreeSet<&'static str>>>,
}

impl Args {
    /// Whether the argument `name` was given on the command line.
    pub(crate) fn given(&self, name: &str) -> bool {
        self.given.contains(&name)
    }

    /// The value of the key `name` as text (`AsString`): given, or its default. Only a key whose
    /// value the reader has made ASCII is read so (`check_text_read`); the others are read by
    /// their bytes (`bytes`), as NCBI's `std::string` holds them.
    pub(crate) fn value(&self, name: &str) -> Option<&str> {
        self.check_text_read(name);
        match self.values.get(name) {
            Some(Some(Value::Text(text))) => Some(text),
            _ => None,
        }
    }

    /// The value of the key `name` as its bytes (`AsString`, NCBI's `std::string`): given, or
    /// its default, UTF-8 or not (on a system whose names are not bytes, the text with each
    /// byte that is not UTF-8 replaced: `utils::os_bytes`).
    pub(crate) fn bytes(&self, name: &str) -> Option<std::borrow::Cow<'_, [u8]>> {
        self.os_value(name)
            .map(crate::utils::os_bytes::os_str_bytes)
    }

    // No NCBI counterpart: NCBI's `AsString` is the value's bytes.
    /// The check, in a debug build, that the key `name` read as text (`value`) is one whose value
    /// the reader has made ASCII: a flag, a number (`eBoolean`, `eInt8`, `eInteger`, `eDouble`),
    /// or a key of a constraint other than a file name's length (a list of choices, a range of
    /// numbers, a genetic code). Every other key's value may be any bytes (`bytes`).
    fn check_text_read(&self, name: &str) {
        debug_assert!(
            find_desc(self.arguments, name).is_some_and(|desc| {
                matches!(
                    desc.kind,
                    Kind::Flag | Kind::Boolean | Kind::Int8 | Kind::Integer | Kind::Double
                ) || desc.constraint.as_ref().is_some_and(|constraint| {
                    !matches!(constraint, Constraint::MaxFileNameLength(_))
                })
            }),
            "-{name} is not a key whose value is ASCII: read it with `Args::bytes`"
        );
    }

    /// The arguments given on the command line, in their order, each with its value's bytes
    /// (`None` for a flag), which LOSAT's refusals of the programs on `cli::reader` read
    /// (`cli::reader::refusal_after_read`).
    pub(crate) fn given_values(&self) -> impl Iterator<Item = (&'static str, Option<&OsStr>)> + '_ {
        self.given.iter().map(|&name| (name, self.os_value(name)))
    }

    /// Whether the value of the key `name` is not UTF-8.
    pub(crate) fn is_bytes(&self, name: &str) -> bool {
        matches!(self.values.get(name), Some(Some(Value::Bytes(_))))
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:509-516
    /// `args[name]` (and `AsBoolean`) of the flag `name`: whether its word was given.
    pub(crate) fn flag(&self, name: &'static str) -> bool {
        self.check_read(name, &[Kind::Flag]);
        self.given(name)
    }

    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each
    // value; LOSAT's ten programs on their own option types (`BlastnArgs`, ...) read them once,
    // in each type's `from_args` (line RG DECISIONS K4i-3-1).
    /// The option `T` of the key `name` when it was given on the command line, `None` when it was
    /// not, whatever its default: an option whose default the program resolves itself (`-out`'s
    /// standard output, `-mt_mode`'s 0, a task's preset). `NotConverted` where `T` does not hold
    /// the value (`OptionValue::of`).
    pub(crate) fn given_option<T: OptionValue>(
        &self,
        name: &'static str,
    ) -> Result<Option<T>, NotConverted> {
        self.check_read(name, T::KINDS);
        if !self.given(name) {
            return Ok(None);
        }
        let value = self.os_value(name).ok_or(NotConverted { name })?;
        T::of(value, self.numbers.get(name).copied())
            .map(Some)
            .ok_or(NotConverted { name })
    }

    /// The option `T` of the key `name`: given, else its default in the program's tables, which
    /// NCBI's reader added (`post_check`) unless a given argument excludes the key; there it is
    /// still the option's value (tblastn's and psiblast's `-query` under `-in_pssm`, `-in_msa`:
    /// line RG DECISIONS K4i-3-1 Q8). `NotConverted` for a key without a value and without a
    /// default (a mandatory key that a given argument excludes), or where `T` does not hold the
    /// value.
    pub(crate) fn option_or_default<T: OptionValue>(
        &self,
        name: &'static str,
    ) -> Result<T, NotConverted> {
        self.check_read(name, T::KINDS);
        let converted = match self.os_value(name) {
            Some(value) => T::of(value, self.numbers.get(name).copied()),
            None => {
                let desc = find_desc(self.arguments, name).ok_or(NotConverted { name })?;
                let default = desc.default().ok_or(NotConverted { name })?;
                let number = process_value(desc, default.as_bytes()).ok().flatten();
                T::of(OsStr::new(default), number)
            }
        };
        converted.ok_or(NotConverted { name })
    }

    /// The check, in a debug build, that `name` is an argument of the program, by its own name
    /// (not an alias), of one of the kinds `kinds`; in a test build, the record of the read.
    fn check_read(&self, name: &'static str, kinds: &[Kind]) {
        debug_assert!(
            find_desc(self.arguments, name)
                .is_some_and(|desc| desc.name.is(name) && kinds.contains(&desc.kind)),
            "-{name} is not an argument of one of the kinds {kinds:?}"
        );
        #[cfg(test)]
        self.read.lock().unwrap().insert(name);
    }

    /// The names that the options' readers have read (`check_read`).
    #[cfg(test)]
    pub(crate) fn names_read(&self) -> BTreeSet<&'static str> {
        self.read.lock().unwrap().clone()
    }

    /// The arguments with the value of each key that `lossy` takes and that is not UTF-8 read as
    /// text, each byte that is not UTF-8 replaced (U+FFFD): what a dry run's caller, which ends
    /// before LOSAT's refusal of such a value (`cli::reader::refusal_after_read`), reads.
    pub(crate) fn with_lossy_text(&self, lossy: impl Fn(&str) -> bool) -> Args {
        let mut args = self.clone();
        for (name, value) in &mut args.values {
            if let Some(Value::Bytes(bytes)) = value {
                if lossy(name) {
                    *value = Some(Value::Text(bytes.to_string_lossy().into_owned()));
                }
            }
        }
        args
    }

    /// `args[name].AsString()` inside `Run`: the value of the key `name`, or NCBI's exception
    /// of an argument without a value (`CArg_NoValue`), which the application reports as one of
    /// the reader's (`Diag::run_arg_exception`).
    /// An argument that a given one excludes throws `eExcludedValue` instead
    /// (`CArg_ExcludedValue`). Read as text as [`Args::value`] reads it.
    pub(crate) fn as_string(&self, name: &str) -> Result<&str, ArgException> {
        if self.excluded.contains(&name) {
            return Err(excluded_value(name));
        }
        self.value(name).ok_or_else(|| no_value(name))
    }

    /// The value of the key `name` as text, in a program that reads its values as text and has
    /// refused those that are not UTF-8 where its body starts (`text_values_refusal`:
    /// blastdb_aliastool): given, or its default.
    pub(crate) fn text(&self, name: &str) -> Option<&str> {
        let value = self.values.get(name);
        debug_assert!(
            !matches!(value, Some(Some(Value::Bytes(_)))),
            "-{name} is not UTF-8: refuse it first (`text_values_refusal`)"
        );
        match value {
            Some(Some(Value::Text(text))) => Some(text),
            _ => None,
        }
    }

    // No NCBI counterpart: NCBI's `AsString` is the value's bytes; blastdb_aliastool reads its
    // values as text.
    /// LOSAT's refusal of a value that is not UTF-8 of an argument that a program reads as text
    /// (`text`: blastdb_aliastool), which the program writes where its body starts, so that NCBI's
    /// errors of the reader and `-dryrun` come first: NCBI's `Run` reads the bytes. It names the
    /// first such argument in the command line's order, in the words of the programs on
    /// `cli::reader` (`cli::non_utf8_value_refusal`). The values read by their bytes are not refused: an
    /// `eOutputFile` key, which the programs open by its name's bytes (`output_file`,
    /// `os_value`), `-logfile`, which the application opens before the reader, and `-conffile`,
    /// whose file the application loads before the reader (`ncbi_environment`).
    pub(crate) fn text_values_refusal(&self, tool: &Tool) -> Option<String> {
        let read_as_bytes = |name: &str| {
            FILE_ARGUMENTS.iter().any(|desc| desc.name.is(name))
                || find_desc(tool_arguments(tool), name)
                    .is_some_and(|desc| desc.kind == Kind::OutputFile)
        };
        let refused = self.given.iter().find(|&&name| {
            matches!(self.values.get(name), Some(Some(Value::Bytes(_)))) && !read_as_bytes(name)
        })?;
        Some(crate::cli::non_utf8_value_refusal(
            tool.refusal_name,
            refused,
        ))
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

    /// The run's files, where the files that the program opens by name for writing go.
    pub(crate) fn files(&self) -> &OpenFiles {
        &self.files
    }

    /// The arguments of a run whose files go to `files` (the application's: `ncbi_app::run_in`).
    pub(crate) fn with_files(self, files: &OpenFiles) -> Self {
        Args {
            files: files.clone(),
            ..self
        }
    }

    /// `args[name].AsOutputFile()`: the stream of the `eOutputFile` key `name`, the file of its
    /// value's bytes as NCBI opens it in the run's files (`files`), `stdout` for `-`, or NCBI's
    /// exception when it does not open. A `CArgException` thrown in `Run` is reported by the
    /// application as one of the reader's (`Diag::run_arg_exception`). The stream of an
    /// `eInputFile` key is `files::open_input_file` of its value, which the programs open where
    /// NCBI's `Run` reads it.
    pub(crate) fn output_file<'a>(
        &self,
        name: &str,
        stdout: &'a mut dyn Write,
    ) -> Result<OutputFile<'a>, ArgException> {
        let value = self.os_value(name).ok_or_else(|| no_value(name))?;
        files::open_output_file(name, value, stdout, &self.files)
    }

    /// The value of the key `name` as its bytes (`AsString` of a name that a program opens):
    /// given, or its default, UTF-8 or not. The programs open a file by it (`output_file`,
    /// `files::open_input_file`).
    pub(crate) fn os_value(&self, name: &str) -> Option<&OsStr> {
        match self.values.get(name) {
            Some(Some(Value::Text(text))) => Some(OsStr::new(text.as_str())),
            Some(Some(Value::Bytes(bytes))) => Some(bytes.as_os_str()),
            _ => None,
        }
    }
}

/// The key that an option's reader (`Args::given_option`, `Args::option_or_default`) could not
/// read into the option's type: LOSAT's error (`cli::reader::conversion_error`), which no value
/// that NCBI's reader took reaches (the constraint of each key that a narrower type reads keeps
/// it in that type's range).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NotConverted {
    /// The key.
    pub(crate) name: &'static str,
}

/// A type of a program's option that `Args` reads from a key's value: the kinds of the keys it
/// reads, and the option of a value (`None` where the type does not hold it). The value's bytes
/// are those of the command line, or of the table's default.
pub(crate) trait OptionValue: Sized {
    /// The kinds of the keys whose values it holds.
    const KINDS: &'static [Kind];
    /// The option of a key's value `value` and its number (`Number`, for the kinds that have
    /// one).
    fn of(value: &OsStr, number: Option<Number>) -> Option<Self>;
}

/// The kinds of the keys whose value is text.
const TEXT_KINDS: &[Kind] = &[Kind::Text, Kind::InputFile, Kind::OutputFile];

/// The text of a key's value (`AsString`): a value that is not UTF-8 is refused before the options
/// are read (`cli::reader::refusal_after_read`), or read lossily by a dry run's caller
/// (`Args::with_lossy_text`).
impl OptionValue for String {
    const KINDS: &'static [Kind] = TEXT_KINDS;
    fn of(value: &OsStr, _: Option<Number>) -> Option<Self> {
        value.to_str().map(str::to_string)
    }
}

/// The bytes of a key's value, by which the program opens its file.
impl OptionValue for std::path::PathBuf {
    const KINDS: &'static [Kind] = TEXT_KINDS;
    fn of(value: &OsStr, _: Option<Number>) -> Option<Self> {
        Some(value.into())
    }
}

/// The bytes of a key's value.
impl OptionValue for OsString {
    const KINDS: &'static [Kind] = TEXT_KINDS;
    fn of(value: &OsStr, _: Option<Number>) -> Option<Self> {
        Some(value.to_os_string())
    }
}

/// `AsBoolean()` of an `eBoolean` key.
impl OptionValue for bool {
    const KINDS: &'static [Kind] = &[Kind::Boolean];
    fn of(_: &OsStr, number: Option<Number>) -> Option<Self> {
        match number? {
            Number::Boolean(value) => Some(value),
            _ => None,
        }
    }
}

/// `AsInt8()` of an `eInt8` key.
impl OptionValue for i64 {
    const KINDS: &'static [Kind] = &[Kind::Int8];
    fn of(_: &OsStr, number: Option<Number>) -> Option<Self> {
        match number? {
            Number::Int8(value) => Some(value),
            _ => None,
        }
    }
}

/// `AsInteger()` of an `eInteger` key.
impl OptionValue for i32 {
    const KINDS: &'static [Kind] = &[Kind::Integer];
    fn of(_: &OsStr, number: Option<Number>) -> Option<Self> {
        match number? {
            Number::Integer(value) => Some(value),
            _ => None,
        }
    }
}

/// `AsInteger()` of an `eInteger` key whose constraint keeps it at or above 0 (a count, a number
/// of threads, a length).
impl OptionValue for usize {
    const KINDS: &'static [Kind] = &[Kind::Integer];
    fn of(value: &OsStr, number: Option<Number>) -> Option<Self> {
        usize::try_from(i32::of(value, number)?).ok()
    }
}

/// `AsInteger()` of an `eInteger` key whose constraint keeps it in 0 to 255 (a genetic code, a
/// mode, a template length).
impl OptionValue for u8 {
    const KINDS: &'static [Kind] = &[Kind::Integer];
    fn of(value: &OsStr, number: Option<Number>) -> Option<Self> {
        u8::try_from(i32::of(value, number)?).ok()
    }
}

/// `AsDouble()` of an `eDouble` key.
impl OptionValue for f64 {
    const KINDS: &'static [Kind] = &[Kind::Double];
    fn of(_: &OsStr, number: Option<Number>) -> Option<Self> {
        match number? {
            Number::Double(value) => Some(value),
            _ => None,
        }
    }
}

/// [`parse_in`] of a run on the file system (the tests' reader).
#[cfg(test)]
pub(crate) fn parse(tool: &Tool, words: &[OsString]) -> Result<Args, ArgsError> {
    parse_in(tool, words, &OpenFiles::disk())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2856-2932
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:744-752
/// Reads `words` (the command line after the program's name, without the words NCBI's
/// application takes before the reader: `ncbi_environment::preparse`) as NCBI's
/// `CArgDescriptions::CreateArgs` reads `tool`'s: each word in its order, then the dependencies
/// and the defaults (`x_PostCheck`). NCBI's programs on this reader have no positional argument,
/// so a word that is not an option, and every word after `--`, is one too many. The words are
/// read as their bytes, as NCBI reads `argv`: a name is the bytes after the dash up to `=`, and
/// a name that is not ASCII is not an argument's name (`is_argument_name`).
///
/// The run's files go to `files`: the file of an `fPreOpen` key (blastdb_aliastool's
/// `-gi_file_out`) is made or emptied there as its word is read, before the reader looks at the
/// next word and at the dependencies (NCBI 2.17.0: `-gi_file_out x.bgl` alone leaves an empty
/// `x.bgl` and the error that it requires `-gi_file_in`, probe G10 of row RG-D2a).
pub(crate) fn parse_in(
    tool: &Tool,
    words: &[OsString],
    files: &OpenFiles,
) -> Result<Args, ArgsError> {
    let arguments = tool_arguments(tool);
    let mut values: BTreeMap<&'static str, Option<Value>> = BTreeMap::new();
    let mut numbers: BTreeMap<&'static str, Number> = BTreeMap::new();
    let mut given: Vec<&'static str> = Vec::new();
    // `n_plain`: `None` until the positional arguments start (NCBI's `kMax_UInt`).
    let mut n_plain: Option<u32> = None;
    let mut at = 0;
    while at < words.len() {
        let arg1 = words[at].as_os_str();
        let bytes = arg1.as_encoded_bytes();
        let arg2 = words.get(at + 1);
        let mut name: &[u8] = b"";
        if n_plain.is_none() {
            if bytes == b"--" {
                n_plain = Some(0);
                at += 1;
                continue;
            }
            if bytes.len() > 1 && bytes[0] == b'-' {
                name = &bytes[1..];
                if let Some(eq) = name.iter().position(|&byte| byte == b'=') {
                    name = &name[..eq];
                }
            } else {
                n_plain = Some(0);
            }
        }
        if let Some(count) = n_plain {
            let count = count + 1;
            let mut message =
                format!("Too many positional arguments ({count}), the offending value: ")
                    .into_bytes();
            message.extend_from_slice(bytes);
            return Err(invalid("eSynopsis", message));
        }
        let Created {
            desc,
            value,
            number,
            used,
        } = create_arg(arguments, arg1, name, arg2)?;
        if desc.pre_open {
            if let Some(value) = &value {
                files::pre_open_output(desc.name(), value, files).map_err(ArgsError::Invalid)?;
            }
        }
        if values.contains_key(desc.name()) {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1928-1941
            return Err(invalid(
                "eSynopsis",
                format!(
                    "Argument with this name is defined already: {}",
                    desc.name()
                ),
            ));
        }
        values.insert(desc.name(), value.map(Value::of));
        if let Some(number) = number {
            numbers.insert(desc.name(), number);
        }
        given.push(desc.name());
        at += if used { 2 } else { 1 };
    }
    let excluded = post_check(tool, arguments, &mut values, &mut numbers, &given)?;
    Ok(Args {
        arguments,
        values,
        numbers,
        given,
        excluded,
        files: files.clone(),
        #[cfg(test)]
        read: Default::default(),
    })
}

// No NCBI counterpart: NCBI has no Web engine.
/// The arguments of the words `given` of `tool` (each name and its value's bytes, `None` for a
/// flag) that another reader has checked: LOSAT Web's ABI v1, whose frozen reader checks its
/// words with its own texts (`web_v1_argv`, line RG DECISIONS K4i-3-1 Q4). Each value is
/// converted by its key's kind, without its constraint (`kind_value`), and no default is added:
/// the options' readers take a key's table default (`Args::option_or_default`). `None` where a
/// name is not one of `tool`'s or a value does not convert.
#[cfg(any(losat_web_abi_v1, test))]
pub(crate) fn checked(tool: &Tool, given: &[(String, Option<OsString>)]) -> Option<Args> {
    let arguments = tool_arguments(tool);
    let mut values: BTreeMap<&'static str, Option<Value>> = BTreeMap::new();
    let mut numbers: BTreeMap<&'static str, Number> = BTreeMap::new();
    let mut names: Vec<&'static str> = Vec::new();
    for (name, value) in given {
        let desc = find_desc(arguments, name)?;
        if let Some(value) = value {
            if let Some(number) = kind_value(desc, value.as_encoded_bytes()).ok()? {
                numbers.insert(desc.name(), number);
            }
        }
        values.insert(desc.name(), value.clone().map(Value::of));
        names.push(desc.name());
    }
    Some(Args {
        arguments,
        values,
        numbers,
        given: names,
        excluded: Vec::new(),
        files: OpenFiles::disk(),
        #[cfg(test)]
        read: Default::default(),
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3007-3081
/// One option word `arg1` with the name `name` (and the next word `arg2`, the value of a key
/// written apart).
fn create_arg(
    arguments: Arguments,
    arg1: &OsStr,
    name: &[u8],
    arg2: Option<&OsString>,
) -> Result<Created, ArgsError> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:921-929
    // The lookup of `name` (`x_Find`) makes a `CArgDesc_NameOnly` of it, so a name that is no
    // argument name fails there.
    let text = is_argument_name(name)
        .then(|| std::str::from_utf8(name).ok())
        .flatten();
    let Some(name) = text else {
        return Err(invalid(
            "eInvalidArg",
            [b"Invalid argument name: ".as_slice(), name].concat(),
        ));
    };
    let bytes = arg1.as_encoded_bytes();
    let find = |name: &str| find_desc(arguments, name);
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
    let eq_separator = bytes.get(name.len() + 1) == Some(&b'=');
    if desc.kind == Kind::Flag {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2840-2852
        match bytes {
            b"-h" => return Err(ArgsError::Usage),
            b"-help" => return Err(ArgsError::Help),
            // `-help-full` adds the hidden arguments: makeblastdb's `-randomize`
            // (`Tool::help_full`); the other programs have none, and their text is that of
            // `-help` (probe C05 of line BDC).
            b"-help-full" => return Err(ArgsError::HelpFull),
            // The help in XML (`ncbi_xmlhelp`).
            b"-xmlhelp" => return Err(ArgsError::XmlHelp),
            _ => {}
        }
        if eq_separator {
            return Err(invalid(
                "eInvalidArg",
                [b"Invalid argument: ".as_slice(), bytes].concat(),
            ));
        }
        return Ok(Created {
            desc,
            value: None,
            number: None,
            used: false,
        });
    }
    let (value, used) = if eq_separator {
        // SAFETY: the value is the rest of the word after the ASCII `=`, so it starts at a
        // boundary of the platform's encoding and ends where the word ends.
        let value = unsafe { OsStr::from_encoded_bytes_unchecked(&bytes[name.len() + 2..]) };
        (value.to_os_string(), false)
    } else {
        match arg2 {
            Some(value) => (value.clone(), true),
            // Without `=` the word is the dash and the name, which is ASCII.
            None => {
                let word = format!("-{name}");
                return Err(invalid(
                    "eNoArg",
                    arg_message(&word, "Value is missing", b""),
                ));
            }
        }
    };
    let number = process_value(desc, value.as_encoded_bytes())?;
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
    value: Option<OsString>,
    /// The number of an `eInteger` or `eDouble` key.
    number: Option<Number>,
    /// Whether the next word was its value.
    used: bool,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3338-3365
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1786-1789
/// `CArgDescriptions::VerifyName`: whether `name` can be an argument's name (ASCII letters
/// and digits, `_` and `-`, not `-` alone or starting with `--`; the empty name can).
fn is_argument_name(bytes: &[u8]) -> bool {
    if bytes.first() == Some(&b'-') && matches!(bytes.get(1), None | Some(b'-')) {
        return false;
    }
    bytes
        .iter()
        .all(|&byte| byte.is_ascii_alphanumeric() || byte == b'_' || byte == b'-')
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1228-1262
/// `CArgDesc::ProcessArgument` (and `ProcessDefault`): the conversion of `value` to the key's
/// type, then its constraint, which reads the text of the value (`Constraint::verify`); the
/// converted value of an `eBoolean`, `eInteger` or `eDouble` key. `value` is the bytes of the
/// value.
fn process_value(desc: &Desc, value: &[u8]) -> Result<Option<Number>, ArgsError> {
    let number = kind_value(desc, value)?;
    match desc.constraint {
        Some(constraint) if !constraint.verify(value) => Err(invalid(
            "eConstraint",
            arg_message(
                desc.name(),
                &format!("Illegal value, expected {}", constraint.usage()),
                value,
            ),
        )),
        _ => Ok(number),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1228-1262
/// The conversion of `value` to the type of the key `desc` (`process_value` before the
/// constraint): the converted value of an `eBoolean`, `eInt8`, `eInteger` or `eDouble` key.
fn kind_value(desc: &Desc, value: &[u8]) -> Result<Option<Number>, ArgsError> {
    Ok(match desc.kind {
        Kind::Integer => Some(Number::Integer(integer_value(desc, value)?)),
        Kind::Int8 => Some(Number::Int8(int8_value(desc, value)?)),
        Kind::Double => Some(double_value(desc, value)?),
        Kind::Boolean => Some(Number::Boolean(boolean_value(desc, value)?)),
        Kind::Text | Kind::Flag | Kind::InputFile | Kind::OutputFile => None,
    })
}

/// The `eConvert` exception of a value of the argument `desc` that its conversion did not read:
/// the conversion's exception is its predecessor, whose text (with the position,
/// `NcbiException::text`, in LOSAT's words) is the first message (`Cannot convert string 'x' to
/// Int8` and the position).
fn not_converted(desc: &Desc, value: &[u8], predecessor: Vec<u8>) -> ArgsError {
    ArgsError::Invalid(ArgException {
        code: "eConvert",
        messages: vec![
            predecessor,
            arg_message(desc.name(), "Argument cannot be converted", value),
        ],
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:351-360
/// The value of an `eInt8` key (`CArg_Int8`: `s_StringToInt8`, `string_to_int8`).
fn int8_value(desc: &Desc, value: &[u8]) -> Result<i64, ArgsError> {
    string_to_int8(value).map_err(|error| not_converted(desc, value, error.text()))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:383-390
/// The value of an `eInteger` key (`CArg_Integer`): that of `CArg_Int8`, in the range of `int`.
fn integer_value(desc: &Desc, value: &[u8]) -> Result<i32, ArgsError> {
    let number = int8_value(desc, value)?;
    i32::try_from(number).map_err(|_| {
        invalid(
            "eConvert",
            arg_message(desc.name(), "Integer value is out of range", value),
        )
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
/// The number of an `eDouble` key (`CArg_Double`).
fn double_value(desc: &Desc, value: &[u8]) -> Result<Number, ArgsError> {
    crate::ncbi_str::string_to_double_posix_or_local(value)
        .map(Number::Double)
        .map_err(|error| not_converted(desc, value, error.text()))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:489-497
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2837-2838
/// The value of an `eBoolean` key (`CArg_Boolean`: `NStr::StringToBool`, whose exception, at
/// position 0, is the predecessor of the reader's; in LOSAT's words, `crate::wording`).
fn boolean_value(desc: &Desc, value: &[u8]) -> Result<bool, ArgsError> {
    let text = std::str::from_utf8(value).ok();
    text.and_then(crate::ncbi_str::string_to_bool)
        .ok_or_else(|| {
            let predecessor =
                crate::wording::losat_text("String cannot be converted to bool (m_Pos = 0)");
            not_converted(desc, value, predecessor.into_owned().into_bytes())
        })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3162-3296
/// `x_PostCheck`: the dependencies of the given arguments, then the defaults of the keys that
/// were not given and that no given argument excludes.
fn post_check(
    tool: &Tool,
    arguments: Arguments,
    values: &mut BTreeMap<&'static str, Option<Value>>,
    numbers: &mut BTreeMap<&'static str, Number>,
    given: &[&'static str],
) -> Result<Vec<&'static str>, ArgsError> {
    let walk = walk_dependencies(arguments.dependencies(), &|name| given.contains(&name))
        .map_err(|message| invalid("eConstraint", message))?;
    // NCBI asks the arguments in the order of its `m_KeyFlagArgs` and stops at the first one that
    // throws; of the arguments that throw, the one first in the usage is the one NCBI names.
    let mut first: Option<(usize, ArgsError)> = None;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3279-3283
    let mut excluded = Vec::new();
    for desc in descs(arguments) {
        if given.contains(&desc.name()) {
            continue;
        }
        let error = if let Some(message) = walk.missing(desc.name()) {
            invalid("eConstraint", message)
        } else if walk.excludes(desc.name()) {
            excluded.push(desc.name());
            continue;
        } else if let Some(default) = desc.default() {
            match process_value(desc, default.as_bytes()) {
                Ok(number) => {
                    if let Some(number) = number {
                        numbers.insert(desc.name(), number);
                    }
                    values.insert(desc.name(), Some(Value::Text(default.to_string())));
                    continue;
                }
                Err(error) => error,
            }
        } else if desc.mandatory {
            mandatory_missing(desc)
        } else {
            continue;
        };
        let position = usage_position(tool.help, desc.name());
        if first
            .as_ref()
            .is_none_or(|(earlier, _)| position < *earlier)
        {
            first = Some((position, error));
        }
    }
    match first {
        Some((_, error)) => Err(error),
        None => Ok(excluded),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1158-1173,1270-1275
/// The exception of the mandatory key `desc` that the command line does not give.
fn mandatory_missing(desc: &Desc) -> ArgsError {
    let mut attr = desc.kind.type_name().to_string();
    if let Some(constraint) = desc.constraint {
        attr += ", ";
        attr += &constraint.usage();
    }
    invalid(
        "eNoArg",
        arg_message(desc.name(), "Mandatory value is missing", attr.as_bytes()),
    )
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3245-3250
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3609-3620
/// Where the argument `name` stands in NCBI's `m_KeyFlagArgs` (the order in which the program's
/// `Init` described its arguments, after the help flags): its place in the synopsis of the usage
/// `help` (a `Tool::help`), which lists the arguments in that order (no program of these sets
/// `fUsageSortArgs`). A hidden argument (makeblastdb's `-randomize`) is not in the usage; it comes
/// after the others, and is a flag that no argument requires, so `x_PostCheck` never stops at it
/// (`args_tests`).
fn usage_position(help: &str, name: &str) -> usize {
    let synopsis = help.split("\n\n").next().unwrap_or("");
    synopsis
        .split_ascii_whitespace()
        .filter_map(|word| word.trim_start_matches('[').strip_prefix('-'))
        .position(|word| word.trim_end_matches(']') == name)
        .unwrap_or(usize::MAX)
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod args_tests;
