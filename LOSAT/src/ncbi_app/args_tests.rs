//! NCBI's argument reader (`args`) on the descriptions of NCBI's programs (their XML,
//! `Arguments`): the conversion of each type, the constraints, the defaults, the flags, the
//! dependencies, the order of the errors and the streams of the file keys. The texts are NCBI
//! BLAST+ 2.17.0's: the dustmasker probe of line RG item K1 (`LOSAT/tests/dustmasker/ncbi`),
//! probe B of line BDC for blastdbcmd, and the verdicts of NCBI's reader on the probe of item
//! K4i (`LOSAT/tests/unit/ncbi_args_verdicts.tsv`).

use super::*;

/// `program` as these tests read its command line. The synopsis of its usage (`Tool::help`),
/// which orders the errors of `x_PostCheck`, is that of the programs on this reader
/// (`ON_THE_READER`); the others' have no two such errors in these tests.
fn tool(program: &'static str) -> Tool {
    let help = match program {
        "blastdbcmd" => include_str!("../blastdbcmd/help.txt"),
        "dustmasker" => include_str!("../dustmasker/help.txt"),
        "segmasker" => include_str!("../segmasker/help.txt"),
        "convert2blastmask" => include_str!("../convert2blastmask/help.txt"),
        "makembindex" => include_str!("../makembindex/help.txt"),
        "windowmasker" => include_str!("../windowmasker/help.txt"),
        _ => "",
    };
    Tool {
        program,
        refusal_name: "T",
        help,
        post_prefix: None,
        components: &[],
        help_full: None,
        short_description: None,
        write_failure: None,
    }
}

/// The programs that run on this reader (`ncbi_app`).
const ON_THE_READER: [&str; 6] = [
    "blastdbcmd",
    "dustmasker",
    "segmasker",
    "convert2blastmask",
    "makembindex",
    "windowmasker",
];

fn parse_words(program: &'static str, line: &[&str]) -> Result<Args, ArgsError> {
    let words: Vec<OsString> = line.iter().map(OsString::from).collect();
    parse(&tool(program), &words)
}

/// The messages of an exception as text (these tests' are UTF-8).
fn texts(messages: &[Vec<u8>]) -> Vec<String> {
    let texts = messages
        .iter()
        .map(|message| String::from_utf8(message.clone()));
    texts.collect::<Result<_, _>>().unwrap()
}

fn messages(program: &'static str, line: &[&str]) -> (&'static str, Vec<String>) {
    match parse_words(program, line) {
        Err(ArgsError::Invalid(ArgException { code, messages })) => (code, texts(&messages)),
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

/// Every program's table.
fn programs() -> impl Iterator<Item = (&'static str, Arguments)> {
    PROGRAMS
        .iter()
        .map(|&(program, arguments)| (program.get(), arguments))
}

/// The text of the XML of `program`, unchanged (`ncbi_xmlhelp::print` refusing nothing).
fn xml(program: &str) -> String {
    crate::ncbi_xmlhelp::print(program, |_| false).unwrap()
}

/// `NStr::XmlEncode` undone (`&amp;`, `&lt;`, `&gt;`, `&apos;`, `&quot;`).
fn xml_decode(text: &str) -> String {
    text.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&apos;", "'")
        .replace("&quot;", "\"")
        .replace("&amp;", "&")
}

/// The elements `<key name=...` and `<flag name=...` of `xml`: (tag, name, the element's text).
fn xml_arguments(xml: &str) -> Vec<(&str, &str, &str)> {
    let mut out = Vec::new();
    for (start, _) in xml.match_indices("\n<") {
        let element = &xml[start + 2..];
        let Some((tag, rest)) = element.split_once(" name=\"") else {
            continue;
        };
        if tag != "key" && tag != "flag" {
            continue;
        }
        let name = rest.split_once('"').unwrap().0;
        let end = element.find(&format!("</{tag}>")).unwrap();
        out.push((tag, name, &element[..end]));
    }
    out
}

/// The `<alias>` lines of the elements `elements` (alias, name) and the entries of the
/// `<dependencies>` element of `xml` (first name, dependency, second name), each in the XML's
/// order: read here apart from the build's reader (`build/ncbi_args.rs`).
fn xml_aliases_and_dependencies(
    xml: &str,
    elements: &[&(&str, &str, &str)],
) -> (Vec<(String, String)>, Vec<(String, Dep, String)>) {
    let text = |within: &str, tag: &str| {
        let (_, rest) = within.split_once(&format!("<{tag}>"))?;
        Some(xml_decode(rest.split_once(&format!("</{tag}>")).unwrap().0))
    };
    let mut aliases = Vec::new();
    for &&(_, name, element) in elements {
        assert!(!element.contains("<negated_alias>"), "-{name}");
        let mut rest = element;
        while let Some(alias) = text(rest, "alias") {
            aliases.push((alias, name.to_string()));
            rest = rest.split_once("</alias>").unwrap().1;
        }
    }
    let mut dependencies = Vec::new();
    if let Some((_, block)) = xml.split_once("\n<dependencies>\n") {
        let block = block.split_once("\n</dependencies>\n").unwrap().0;
        for entry in block.split("<first_").skip(1) {
            let dep = match entry.split_once("_second>").unwrap().0 {
                "requires" => Dep::Requires,
                "excludes" => Dep::Excludes,
                other => panic!("{other}"),
            };
            let (first, second) = (text(entry, "arg1"), text(entry, "arg2"));
            dependencies.push((first.unwrap(), dep, second.unwrap()));
        }
    }
    (aliases, dependencies)
}

#[test]
fn every_program_of_the_xml_has_its_table() {
    let tables: Vec<&str> = programs().map(|(program, _)| program).collect();
    let xmls: Vec<&str> = crate::ncbi_xmlhelp::programs().collect();
    assert_eq!(tables.len(), xmls.len());
    for program in xmls {
        assert!(tables.contains(&program), "{program}");
    }
}

/// The XML read again on its own, apart from the build's reader: each argument of the program's
/// own (not those of every `CArgDescriptions`) is in the table, with its type, its default, its
/// mandatory state, and its constraint, whose usage is the XML's text.
#[test]
fn every_argument_of_the_xml_is_in_its_table() {
    let standard = [
        "h",
        "help",
        "help-full",
        "xmlhelp",
        "logfile",
        "conffile",
        "version",
        "version-full",
        "version-full-xml",
        "version-full-json",
        "dryrun",
    ];
    let mut constraints = 0;
    let (mut aliases, mut dependencies) = (0, 0);
    for (program, arguments) in programs() {
        let xml = xml(program);
        let elements = xml_arguments(&xml);
        let own: Vec<_> = elements
            .iter()
            .filter(|(_, name, _)| !standard.contains(name))
            .collect();
        assert_eq!(own.len(), arguments.descs().len(), "{program}");
        // The aliases and the dependencies, every entry in the XML's order: NCBI writes the
        // `eRequires` entries of its multimap, then the `eExcludes` ones.
        let (xml_aliases, xml_dependencies) = xml_aliases_and_dependencies(&xml, &own);
        let table_aliases: Vec<_> = arguments
            .aliases()
            .map(|(alias, name)| (alias.to_string(), name.to_string()))
            .collect();
        assert_eq!(table_aliases, xml_aliases, "{program}");
        let table_dependencies: Vec<_> = arguments
            .dependencies()
            .map(|(first, dep, second)| (first.to_string(), dep, second.to_string()))
            .collect();
        assert_eq!(table_dependencies, xml_dependencies, "{program}");
        let excludes = xml_dependencies
            .iter()
            .position(|(_, dep, _)| *dep == Dep::Excludes)
            .unwrap_or(xml_dependencies.len());
        assert!(
            xml_dependencies[excludes..]
                .iter()
                .all(|(_, dep, _)| *dep == Dep::Excludes),
            "{program}"
        );
        aliases += xml_aliases.len();
        dependencies += xml_dependencies.len();
        for &&(tag, name, element) in &own {
            let desc = arguments.descs().iter().find(|desc| desc.name() == name);
            let desc = desc.unwrap_or_else(|| panic!("{program} -{name}"));
            if tag == "flag" {
                assert_eq!(desc.kind, Kind::Flag, "{program} -{name}");
            } else {
                let kind = element
                    .split_once(" type=\"")
                    .and_then(|(_, rest)| rest.split_once('"'))
                    .unwrap()
                    .0;
                assert_eq!(desc.kind.type_name(), kind, "{program} -{name}");
                let optional = element
                    .lines()
                    .next()
                    .unwrap()
                    .contains("optional=\"true\"");
                assert_eq!(desc.mandatory, !optional, "{program} -{name}");
            }
            let default = element
                .split_once("<default>")
                .map(|(_, rest)| xml_decode(rest.split_once("</default>").unwrap().0));
            assert_eq!(
                desc.default().map(str::to_string),
                default,
                "{program} -{name}"
            );
            let usage = element.split_once("<constraint").map(|(_, rest)| {
                let description = rest.split_once("<description>").unwrap().1;
                xml_decode(description.split_once("</description>").unwrap().0)
            });
            assert_eq!(
                desc.constraint.map(Constraint::usage),
                usage,
                "{program} -{name}"
            );
            constraints += usize::from(usage.is_some());
        }
    }
    // The 169 constraints of the 14 programs of line RG item K4i (SURVEY), and makembindex's 6
    // and windowmasker's 16.
    assert_eq!(constraints, 169 + 6 + 16);
    // rpsbproc's 8 aliases; the 1,294 entries of the dependencies of the 14 programs,
    // windowmasker's 76 (its 38 exclusions, each in both directions) and makembindex's 3.
    assert_eq!(aliases, 8);
    assert_eq!(dependencies, 1294 + 76 + 3);
}

/// Which class a constraint is.
fn class(constraint: Constraint) -> &'static str {
    match constraint {
        Constraint::Strings(_) => "Strings",
        Constraint::AtLeast(_) => "AtLeast",
        Constraint::AtMost(_) => "AtMost",
        Constraint::Between {
            inclusive: true, ..
        } => "Between inclusive",
        Constraint::Between {
            inclusive: false, ..
        } => "Between exclusive",
        Constraint::Doubles { .. } => "Doubles",
        Constraint::Integers { .. } => "Integers",
        Constraint::StringSet(_) => "StringSet",
        Constraint::IntegerSet(_) => "IntegerSet",
        Constraint::GeneticCode { .. } => "GeneticCode",
        Constraint::MaxFileNameLength(_) => "MaxFileNameLength",
    }
}

const CLASSES: [&str; 11] = [
    "Strings",
    "AtLeast",
    "AtMost",
    "Between inclusive",
    "Between exclusive",
    "Doubles",
    "Integers",
    "StringSet",
    "IntegerSet",
    "GeneticCode",
    "MaxFileNameLength",
];

#[test]
fn every_class_is_used() {
    let used: BTreeSet<&str> = programs()
        .flat_map(|(_, arguments)| arguments.descs().iter())
        .filter_map(|desc| desc.constraint.map(class))
        .collect();
    assert_eq!(used, CLASSES.into_iter().collect());
}

#[test]
fn every_default_passes_its_own_checks() {
    for (program, arguments) in programs() {
        for desc in arguments.descs() {
            if let Some(default) = desc.default() {
                assert!(
                    process_value(desc, default.as_bytes()).is_ok(),
                    "{program} -{} {default:?}",
                    desc.name()
                );
            }
        }
    }
}

/// One case of NCBI's verdicts (`ncbi_args_verdicts.tsv`).
struct Verdict {
    case: &'static str,
    program: &'static str,
    /// `ok`, or the code of the `CArgException`.
    code: &'static str,
    messages: Vec<&'static str>,
    words: Vec<&'static str>,
}

fn verdicts() -> Vec<Verdict> {
    let table = include_str!("../../tests/unit/ncbi_args_verdicts.tsv");
    table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let count: usize = fields[3].parse().unwrap();
            Verdict {
                case: fields[0],
                program: fields[1],
                code: fields[2],
                messages: fields[4..4 + count].to_vec(),
                words: fields[4 + count..].to_vec(),
            }
        })
        .collect()
}

/// The program of a case: the 14 programs have one name each.
fn program_name(program: &str) -> &'static str {
    programs()
        .map(|(name, _)| name)
        .find(|name| *name == program)
        .unwrap()
}

/// Every case of the probe, on every program's table, as NCBI's reader read it: the same
/// acceptance, or the same exception (its code and its messages).
#[test]
fn the_reader_gives_ncbi_s_verdicts() {
    let cases = verdicts();
    assert_eq!(cases.len(), 478);
    let mut differing = Vec::new();
    for case in &cases {
        let got = parse_words(program_name(case.program), &case.words);
        let same = match (&got, case.code) {
            (Ok(_), "ok") => true,
            (Err(ArgsError::Invalid(error)), code) => {
                error.code == code && texts(&error.messages) == strings(&case.messages)
            }
            _ => false,
        };
        if !same {
            differing.push(format!("{}: {got:?}", case.case));
        }
    }
    assert!(differing.is_empty(), "{}", differing.join("\n"));
}

/// The cases reach every constraint: each constrained argument of each program of the probe
/// refuses a value outside it, and each constraint (its usage) takes a value inside it in some
/// program. (makembindex and windowmasker, whose XML came after the probe, are not in it; their
/// constraints are of the classes the probe reaches, and their refusals are NCBI's in their own
/// fixtures, `tests/makembindex`, `tests/windowmasker`.)
#[test]
fn the_verdicts_reach_every_constraint() {
    let cases = verdicts();
    let mut inside: BTreeSet<String> = BTreeSet::new();
    let mut classes: BTreeSet<&str> = BTreeSet::new();
    for case in cases.iter().filter(|case| case.code == "ok") {
        let arguments = arguments(case.program).unwrap();
        for word in &case.words {
            let Some(name) = word.strip_prefix('-') else {
                continue;
            };
            let name = name.split('=').next().unwrap();
            if let Some(constraint) = find_desc(arguments, name).and_then(|desc| desc.constraint) {
                inside.insert(constraint.usage());
                classes.insert(class(constraint));
            }
        }
    }
    assert_eq!(classes, CLASSES.into_iter().collect());
    let probed: BTreeSet<&str> = cases.iter().map(|case| case.program).collect();
    assert_eq!(probed.len(), 14);
    for (program, arguments) in programs().filter(|(program, _)| probed.contains(program)) {
        for desc in arguments.descs() {
            let Some(constraint) = desc.constraint else {
                continue;
            };
            let refusal = format!(
                "Argument \"{}\". Illegal value, expected {}",
                desc.name(),
                constraint.usage()
            );
            let outside = cases.iter().any(|case| {
                case.program == program
                    && case.code == "eConstraint"
                    && case.messages.last().unwrap().starts_with(&refusal)
            });
            assert!(outside, "{program} -{}", desc.name());
            assert!(
                inside.contains(&constraint.usage()),
                "{program} -{}",
                desc.name()
            );
        }
    }
}

/// `x_PostCheck` stops at the first argument of `m_KeyFlagArgs` that throws, which is not in the
/// XML: `usage_position` reads it from the synopsis. Every argument of the programs on this
/// reader is in their synopsis, and the hidden arguments of the others are flags that no
/// argument requires, so none of them is an error of `x_PostCheck`.
#[test]
fn the_synopsis_orders_the_errors_of_the_post_check() {
    for program in ON_THE_READER {
        let help = tool(program).help;
        for desc in arguments(program).unwrap().descs() {
            assert_ne!(
                usage_position(help, desc.name()),
                usize::MAX,
                "{program} -{}",
                desc.name()
            );
        }
    }
    let mut hidden = 0;
    for (program, arguments) in programs() {
        let xml = xml(program);
        for (tag, name, element) in xml_arguments(&xml) {
            if !element.contains("<hidden/>") || arguments.descs().iter().all(|d| d.name() != name)
            {
                continue;
            }
            hidden += 1;
            assert_eq!(tag, "flag", "{program} -{name}");
            let required = arguments
                .dependencies()
                .any(|(_, dep, second)| dep == Dep::Requires && second == name);
            assert!(!required, "{program} -{name}");
        }
    }
    assert_eq!(hidden, 1, "makeblastdb -randomize");
    // With the synopsis turned around, the first of two mandatory keys missing is the other one.
    let turned = Tool {
        help: "USAGE\n  convert2blastmask -masking_options x -masking_algorithm y\n\n",
        ..tool("convert2blastmask")
    };
    assert_eq!(
        parse(&turned, &[]).unwrap_err(),
        ArgsError::Invalid(ArgException::new(
            "eNoArg",
            "Argument \"masking_options\". Mandatory value is missing:  `String'".to_string()
        ))
    );
    assert_eq!(
        messages("convert2blastmask", &[]),
        (
            "eNoArg",
            strings(&["Argument \"masking_algorithm\". Mandatory value is missing:  `String'"])
        )
    );
}

/// An alias names its argument (`x_Find`): rpsbproc's `--evalue` is `-e`.
#[test]
fn an_alias_is_its_argument() {
    let args = parse_words("rpsbproc", &["--evalue", "0.5", "--data-mode=full"]).unwrap();
    assert_eq!(args.value("e"), Some("0.5"));
    assert_eq!(args.double("e"), Some(0.5));
    assert_eq!(args.value("m"), Some("full"));
    assert!(args.given("e") && args.given("m"));
    assert_eq!(
        messages("rpsbproc", &["--evalue", "0.5", "-e", "1"]),
        (
            "eSynopsis",
            strings(&["Argument with this name is defined already: e"])
        )
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:115-129, 351-360, 383-390 (cited at
// `constraint::string_to_int8`, `int8_value` and `integer_value`).
#[test]
fn integer_keys_convert_as_ncbi_s_reader() {
    let value = |text: &str| {
        parse_words("dustmasker", &["-window", text])
            .unwrap()
            .integer("window")
    };
    assert_eq!(value("+5"), Some(5));
    assert_eq!(value("-1"), Some(-1));
    assert_eq!(value("-2147483648"), Some(i32::MIN));
    assert_eq!(value("2147483647"), Some(i32::MAX));
    // Not octal: a leading zero is a decimal digit.
    assert_eq!(value("010"), Some(10));
    // Base 16 after base 10 fails, for a value that starts with `0x` (`big_win_0x10`).
    assert_eq!(value("0x10"), Some(16));
    assert_eq!(value("0XfF"), Some(255));
    assert_eq!(value("0x"), Some(0));
    let cases: &[(&str, &[&str])] = &[
        (
            "x",
            &[
                "Cannot convert string 'x' to Int8 (m_Pos = 0)",
                "Argument \"window\". Argument cannot be converted:  `x'",
            ],
        ),
        (
            "",
            &[
                "Cannot convert string '' to Int8 (m_Pos = 0)",
                "Argument \"window\". Argument cannot be converted",
            ],
        ),
        (
            "2.5",
            &[
                "Cannot convert string '2.5' to Int8 (m_Pos = 1)",
                "Argument \"window\". Argument cannot be converted:  `2.5'",
            ],
        ),
        // No space is skipped, after the number or before it (NCBI on `64 ` and ` 64`: probes
        // dust.window.08 and 07 of the line RG K4i review).
        (
            "5 ",
            &[
                "Cannot convert string '5 ' to Int8 (m_Pos = 1)",
                "Argument \"window\". Argument cannot be converted:  `5 '",
            ],
        ),
        (
            " 5",
            &[
                "Cannot convert string ' 5' to Int8 (m_Pos = 0)",
                "Argument \"window\". Argument cannot be converted:  ` 5'",
            ],
        ),
        (
            "9223372036854775808",
            &[
                "Cannot convert string '9223372036854775808' to Int8, overflow (m_Pos = 18)",
                "Argument \"window\". Argument cannot be converted:  `9223372036854775808'",
            ],
        ),
        // The exception of the conversion in base 16 is the one reported.
        (
            "0x1g",
            &[
                "Cannot convert string '0x1g' to Int8 (m_Pos = 3)",
                "Argument \"window\". Argument cannot be converted:  `0x1g'",
            ],
        ),
        (
            "2147483648",
            &["Argument \"window\". Integer value is out of range:  `2147483648'"],
        ),
        (
            "-2147483649",
            &["Argument \"window\". Integer value is out of range:  `-2147483649'"],
        ),
    ];
    for &(text, expected) in cases {
        assert_eq!(
            messages("dustmasker", &["-window", text]),
            ("eConvert", strings(expected)),
            "{text:?}"
        );
    }
}

/// A constraint reads the text of the value, after the conversion of its type: a hexadecimal
/// value is a number to both (blastdbcmd `-line_length`, `>=1`).
#[test]
fn a_constraint_reads_the_text_of_the_value() {
    let line_length = |text: &str| {
        parse_words("blastdbcmd", &["-line_length", text]).map(|args| args.integer("line_length"))
    };
    assert_eq!(line_length("0x5"), Ok(Some(5)));
    assert_eq!(line_length("+5"), Ok(Some(5)));
    assert_eq!(
        messages("blastdbcmd", &["-line_length", "0x"]),
        (
            "eConstraint",
            strings(&["Argument \"line_length\". Illegal value, expected >=1:  `0x'"])
        )
    );
}

#[test]
fn defaults_and_flags_are_ncbi_s_values() {
    let args = parse_words("dustmasker", &[]).unwrap();
    assert_eq!(args.integer("window"), Some(64));
    assert_eq!(args.value("window"), Some("64"));
    assert_eq!(args.value("in"), Some("-"));
    assert!(!args.flag("parse_seqids") && !args.given("in"));
    let args = parse_words("dustmasker", &["-parse_seqids", "-level=3", "-in", "x.fa"]).unwrap();
    assert!(args.flag("parse_seqids"));
    assert_eq!(args.integer("level"), Some(3));
    assert_eq!(args.value("in"), Some("x.fa"));
    assert!(args.given("in") && !args.given("window"));
    let args = parse_words("segmasker", &["-locut", "1e"]).unwrap();
    assert_eq!(args.double("locut"), Some(1.0));
    assert_eq!(args.double("hicut"), Some(2.5));
    // A file key without a value has no stream (`CArg_NoValue`).
    assert_eq!(
        args.output_file("logfile", &mut Vec::new()).err().unwrap(),
        ArgException::new(
            "eNoValue",
            "Argument \"logfile\". The argument has no value".to_string()
        )
    );
}

#[test]
fn every_reader_knows_the_help_flags_and_the_file_arguments() {
    let parse_words = |line: &[&str]| parse_words("dustmasker", line);
    let messages = |line: &[&str]| messages("dustmasker", line);
    assert_eq!(parse_words(&["-h"]).unwrap_err(), ArgsError::Usage);
    assert_eq!(
        parse_words(&["-help-full"]).unwrap_err(),
        ArgsError::HelpFull
    );
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
        messages(&["-parse_seqids", "true"]),
        (
            "eSynopsis",
            strings(&["Too many positional arguments (1), the offending value: true"])
        )
    );
}

/// The generated tables hold no pointer: one text (`TEXT`) and rows of integers, so that loading
/// LOSAT relocates nothing of them (a `&str` or a slice in a static is an address that the dynamic
/// loader relocates in every run of every program; line RG item K4i-1 fix 1). Every line of the
/// generated file but its comments and `TEXT` is made of the row forms alone: no reference, no
/// string, no other punctuation; names in capitals only where a table is declared (a use of another
/// static could copy a reference); lower-case names only for the rows' fields, `static`, the
/// booleans and the integer type of a table (a function or a variable could be a pointer).
#[test]
fn the_tables_hold_no_pointer() {
    const LOWER: [&str; 22] = [
        "name",
        "kind",
        "default",
        "constraint",
        "mandatory",
        "first",
        "count",
        "at",
        "len",
        "min",
        "max",
        "inclusive",
        "descs",
        "aliases",
        "dependencies",
        "refusals",
        "with_32",
        "static",
        "true",
        "false",
        "i32",
        "u32",
    ];
    let generated = include_str!(concat!(env!("OUT_DIR"), "/ncbi_args.rs"));
    let mut texts = 0;
    let mut tables = 0;
    for line in generated.lines() {
        if line.starts_with("//") {
            continue;
        }
        if line.starts_with("const TEXT: &str = \"") {
            texts += 1;
            continue;
        }
        assert!(
            line.chars()
                .all(|c| c.is_ascii_alphanumeric() || " _()[]{},.:;=-".contains(c)),
            "{line}"
        );
        let declares = line.starts_with("static ");
        tables += usize::from(declares);
        let bytes = line.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let c = bytes[i];
            let end = i + bytes[i..]
                .iter()
                .position(|b| {
                    !(b.is_ascii_alphanumeric() || *b == b'_' || (c.is_ascii_digit() && *b == b'.'))
                })
                .unwrap_or(bytes.len() - i);
            if end == i {
                i += 1;
                continue;
            }
            let word = &line[i..end];
            if c.is_ascii_lowercase() {
                assert!(LOWER.contains(&word), "{word} in {line}");
            } else if c.is_ascii_uppercase()
                && word.len() > 1
                && word.bytes().all(|b| !b.is_ascii_lowercase())
            {
                assert!(declares, "{word} in {line}");
            }
            i = end;
        }
    }
    assert_eq!(texts, 1);
    assert!(tables > 0);
}

// No NCBI counterpart: PD-TLOSAN-LOCAL-GENCODE-32. tblastn's `-db_gencode` allows LOSAT's code 32
// beside NCBI's codes, with NCBI's usage text; tblastx's keys and the other programs' keep NCBI's
// codes alone.
#[test]
fn tblastn_s_db_gencode_alone_allows_code_32() {
    let constraint = |program, name| find_desc(arguments(program).unwrap(), name)?.constraint;
    let ncbi: Vec<u32> = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16]
        .into_iter()
        .chain(21..=31)
        .chain([33])
        .collect();
    let mut with_32 = ncbi.clone();
    with_32.insert(with_32.len() - 1, 32);
    assert_eq!(genetic_codes("tblastn", "db_gencode"), Some(with_32));
    let tblastn = constraint("tblastn", "db_gencode").unwrap();
    assert!(tblastn.verify(b"32") && tblastn.verify(b"33") && !tblastn.verify(b"7"));
    assert_eq!(tblastn.usage(), "values between: 1-6, 9-16, 21-31, 33");
    for (program, name) in [
        ("tblastx", "db_gencode"),
        ("tblastx", "query_gencode"),
        ("blastx", "query_gencode"),
        ("rpstblastn", "query_gencode"),
    ] {
        assert_eq!(
            genetic_codes(program, name),
            Some(ncbi.clone()),
            "{program} -{name}"
        );
        let constraint = constraint(program, name).unwrap();
        assert!(
            !constraint.verify(b"32") && constraint.verify(b"1"),
            "{program} -{name}"
        );
    }
    assert_eq!(genetic_codes("tblastn", "evalue"), None);
    assert_eq!(genetic_codes("blastp", "db_gencode"), None);
}

// The one list of LOSAT's refusals of NCBI's arguments (`REFUSED` of `build/ncbi_args.rs`): how
// each is refused, and the values of a set of strings that describe reads.
#[test]
fn the_refusals_and_the_sets_come_from_the_tables() {
    assert_eq!(refusal("blastn", "remote"), Some(Refusal::Read(None)));
    assert_eq!(refusal("blastn", "db"), None);
    assert_eq!(refusal("nosuchprogram", "remote"), None);
    assert_eq!(
        refusal("psiblast", "use_sw_tback"),
        Some(Refusal::Run("the Smith-Waterman traceback"))
    );
    // tblastn's `-ungapped` runs (line NP item O4c).
    assert_eq!(refusal("tblastn", "ungapped"), None);
    // makeblastdb's `-gi_mask` and `-gi_mask_name` run (line NP item M4c-1): it refuses none.
    assert_eq!(refusal("makeblastdb", "gi_mask"), None);
    assert_eq!(refusal("makeblastdb", "gi_mask_name"), None);
    // The programs on NCBI's reader refuse none of their arguments.
    // (No count of the rows: un-refusing an argument is deleting its row of `REFUSED`.)
    let rows = refusals();
    for program in [
        "blastdbcmd",
        "convert2blastmask",
        "dustmasker",
        "makembindex",
        "segmasker",
        "windowmasker",
        "rpsbproc",
    ] {
        assert!(!rows.iter().any(|row| row.0 == program), "{program}");
    }
    assert_eq!(
        string_values("blastn", "task"),
        Some(vec![
            "blastn",
            "blastn-short",
            "dc-megablast",
            "megablast",
            "rmblastn"
        ])
    );
    assert_eq!(string_values("tblastx", "task"), None);
    assert_eq!(string_values("blastn", "evalue"), None);
}
