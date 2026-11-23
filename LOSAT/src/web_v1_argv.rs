//! The words of a request of LOSAT Web's ABI v1 (frozen, plan TD-1), as its BLASTX and TBLASTN
//! read them: with the reader that LOSAT's command line had before it read its arguments with
//! NCBI's argument reader (`cli::reader`). v1 keeps its own copy of that reader and of its texts
//! (line RG DECISIONS K4i-2 D5), so that its frozen messages stay
//! (`tests/fixtures/web_v1_tblastn`); no other caller reads an argv here.
//!
//! No NCBI counterpart: NCBI has no Web engine.

use std::ffi::OsString;
use std::path::PathBuf;

use clap::builder::{PossibleValuesParser, ValueParser};
use clap::error::ErrorKind;
use clap::{CommandFactory, FromArgMatches};

use crate::blastinput::value_parsers::{genetic_code_number, ncbi_double, ncbi_int8, ncbi_integer};
use crate::cli::{Cli, Commands};
use crate::utils::genetic_code::GeneticCode;

/// The command of v1's `program` (BLASTX or TBLASTN): the command line's (`cli::command_of`), with
/// the value parsers and the exclusions of the arguments that v1's reader had when it was frozen
/// (`v1_parsers`, `V1_EXCLUSIONS`): the command line's arguments are NCBI's reader's since line RG
/// item K4i-2b, clap only converts them.
pub(crate) fn command(program: &str) -> clap::Command {
    crate::cli::command_of(program).mut_subcommand(program, |mut scope| {
        for (name, parser) in v1_parsers(program) {
            scope = scope.mut_arg(name, |arg| arg.value_parser(parser));
        }
        for (name, excluded) in V1_EXCLUSIONS {
            if scope.get_arguments().any(|arg| arg.get_id() == *name) {
                scope = scope.mut_arg(*name, |arg| arg.conflicts_with_all(*excluded));
            }
        }
        if program == "tblastn" {
            scope = scope.mut_arg("in_pssm", |arg| {
                arg.conflicts_with_all(["query", "query_loc"])
            });
        }
        scope.mut_arg("mt_mode", |arg| arg.requires("num_threads"))
    })
}

/// The exclusions of v1's BLASTX and TBLASTN (the clap attributes of their arguments when v1 was
/// frozen): an argument and the arguments it excludes.
const V1_EXCLUSIONS: &[(&str, &[&str])] = &[
    (
        "subject_loc",
        &crate::blastinput::db_filter::DB_FILTER_OPTIONS,
    ),
    ("db", &["subject", "subject_loc"]),
    ("best_hit_overhang", &["culling_limit"]),
    ("best_hit_score_edge", &["culling_limit"]),
    ("max_target_seqs", &["num_descriptions", "num_alignments"]),
    (
        "seqidlist",
        &[
            "negative_seqidlist",
            "taxids",
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "no_taxid_expansion",
            "subject",
        ],
    ),
    (
        "negative_seqidlist",
        &[
            "taxids",
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "no_taxid_expansion",
            "subject",
        ],
    ),
    (
        "taxids",
        &[
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "subject",
        ],
    ),
    (
        "taxidlist",
        &["negative_taxids", "negative_taxidlist", "subject"],
    ),
    ("negative_taxids", &["negative_taxidlist", "subject"]),
    ("negative_taxidlist", &["subject"]),
    ("no_taxid_expansion", &["subject"]),
    ("db_soft_mask", &["subject", "subject_loc", "db_hard_mask"]),
    ("db_hard_mask", &["subject", "subject_loc"]),
];

/// The value parsers of v1's `program` (BLASTX or TBLASTN) where they are not the command line's
/// conversion: NCBI's constraints, checked with v1's own messages.
fn v1_parsers(program: &str) -> Vec<(&'static str, ValueParser)> {
    let mut parsers = vec![
        ("searchsp", ValueParser::new(nonnegative_ncbi_int8)),
        ("out", ValueParser::new(output_path())),
        ("export_search_strategy", ValueParser::new(output_path())),
        ("word_size", ValueParser::new(protein_word_size)),
        ("qcov_hsp_perc", ValueParser::new(percentage)),
        (
            "max_intron_length",
            ValueParser::new(nonnegative_ncbi_integer),
        ),
        ("window_size", ValueParser::new(nonnegative_ncbi_integer)),
        ("culling_limit", ValueParser::new(nonnegative_ncbi_integer)),
        ("threshold", ValueParser::new(ncbi_threshold)),
        ("best_hit_overhang", ValueParser::new(best_hit_value)),
        ("best_hit_score_edge", ValueParser::new(best_hit_value)),
        ("num_descriptions", ValueParser::new(report_count)),
        ("num_alignments", ValueParser::new(report_count)),
        ("mt_mode", ValueParser::new(ncbi_mt_mode)),
        ("line_length", ValueParser::new(formatting_line_length)),
        ("sort_hits", ValueParser::new(formatting_sort_option)),
        ("sort_hsps", ValueParser::new(formatting_sort_option)),
    ];
    if program == "blastx" {
        parsers.extend([
            (
                "task",
                ValueParser::new(PossibleValuesParser::new(["blastx", "blastx-fast"])),
            ),
            ("max_hsps", ValueParser::new(ncbi_positive_integer)),
            (
                "strand",
                ValueParser::new(PossibleValuesParser::new(["both", "plus", "minus"])),
            ),
            ("query_gencode", ValueParser::new(genetic_code)),
            ("max_target_seqs", ValueParser::new(ncbi_positive_integer)),
            ("num_threads", ValueParser::new(ncbi_positive_integer)),
        ]);
    } else {
        parsers.extend([
            (
                "task",
                ValueParser::new(PossibleValuesParser::new(["tblastn", "tblastn-fast"])),
            ),
            ("db_gencode", ValueParser::new(tblastn_genetic_code)),
            ("max_target_seqs", ValueParser::new(blastn_count)),
            ("max_hsps", ValueParser::new(blastn_count)),
            ("num_threads", ValueParser::new(blastn_count)),
        ]);
    }
    parsers
}

/// The arguments of the request `argv` (`losat <program> <words>`, BLASTX or TBLASTN) read with
/// `command` (the program's command, `cli::command_of`, with v1's own value parsers), or the
/// message of their error as v1 gives it (`render`). ABI v1 runs no NCBI application start: the
/// versions and the toolkit's words that the start reads (`-logfile`, `-conffile`, `-dryrun`) are
/// refused; its callers refuse the help words before.
pub(crate) fn parse(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, String> {
    read(command, argv).map_err(|error| render(&error))
}

/// `parse` without the rendering of the error.
fn read(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, clap::Error> {
    use crate::blastinput::ncbi_environment::{preparse, Preparsed};
    let program = argv
        .get(1)
        .and_then(|word| word.to_str())
        .unwrap_or_default()
        .to_string();
    let words = argv.get(2..).unwrap_or(&[]);
    // The refusals name the program as its other refusals do (`Tool::refusal_name`).
    let tool =
        crate::cli::reader::tool(&program).expect("v1 reads the words of one of its programs");
    let words = match preparse(words) {
        Preparsed::Version(request) => {
            return Err(raw(crate::cli::version_word_refusal(tool, request.word())))
        }
        Preparsed::Words { words: read, .. } => {
            if let Some(word) = crate::cli::toolkit_file_word(words, &program) {
                return Err(raw(crate::cli::toolkit_file_word_refusal(tool, word)));
            }
            read
        }
    };
    let scope = command
        .find_subcommand(&program)
        .expect("v1 reads the words of one of its programs");
    let mut translated: Vec<OsString> = vec!["losat".into(), program.clone().into()];
    let mut input = words.into_iter();
    while let Some(token) = input.next() {
        let text = token.to_str().unwrap_or("").to_string();
        // A last `--` changes nothing in TBLASTN, as NCBI's reader; BLASTX's reader of v1 did not
        // take it.
        if text == "--" && input.len() == 0 && program == "tblastn" {
            continue;
        }
        let key = text.strip_prefix('-').filter(|key| !key.starts_with('-'));
        let (name, inline) = key
            .unwrap_or("")
            .split_once('=')
            .map_or((key.unwrap_or(""), None), |(k, v)| (k, Some(v)));
        let Some(arg) = scope
            .get_arguments()
            .find(|arg| arg.get_long() == Some(name))
        else {
            return Err(unknown_option(&program, name, &text));
        };
        if arg.get_action().takes_values() {
            let value = match inline {
                Some(value) => OsString::from(value),
                None => input.next().ok_or_else(|| {
                    clap::Error::raw(
                        ErrorKind::InvalidValue,
                        format!("-{name} requires one value"),
                    )
                })?,
            };
            let mut word = OsString::from(format!("--{name}="));
            word.push(value);
            translated.push(word);
        } else {
            if inline.is_some() {
                return Err(clap::Error::raw(
                    ErrorKind::TooManyValues,
                    format!("-{name} is a flag and does not take a value"),
                ));
            }
            translated.push(format!("--{name}").into());
        }
    }
    let mut matches = command.clone().try_get_matches_from(&translated)?;
    Cli::from_arg_matches_mut(&mut matches)
        .map(|cli| cli.command)
        .map_err(|error| error.format(&mut command.clone()))
}

/// A refusal of v1's reader, as its reader wrote them (a clap error of LOSAT's text).
fn raw(message: String) -> clap::Error {
    clap::Error::raw(ErrorKind::InvalidValue, message)
}

/// v1's rejection of an option that the program does not have: NCBI's options that LOSAT
/// refuses by name, and the rest.
fn unknown_option(program: &str, name: &str, text: &str) -> clap::Error {
    let tool = crate::cli::reader::tool(program).expect("v1's programs are on the reader");
    match crate::cli::reader::named_refusal(tool, name) {
        Some(refusal) => raw(refusal),
        None => clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!("unknown option or argument '{text}'; use -help to list the options"),
        ),
    }
}

/// The text of v1's `error`: clap's, with each program's long options written with one dash, as
/// v1's messages have them.
fn render(error: &clap::Error) -> String {
    let mut message = error
        .to_string()
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let command = Cli::command();
    let toolkit = ["logfile", "conffile"].into_iter();
    let longs = std::iter::once(&command)
        .chain(command.get_subcommands())
        .flat_map(|scope| scope.get_arguments().filter_map(clap::Arg::get_long))
        .map(str::to_string)
        .chain(toolkit.map(str::to_string));
    for name in longs.filter(|name| !matches!(name.as_str(), "help" | "version") && name.len() > 1)
    {
        message = message.replace(&format!("--{name}"), &format!("-{name}"));
    }
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
        .replace("see a summary with '-help'", "see a summary with '-h'")
        .replace("see more with '--help'", "see more with '-help'")
}

// v1's copies of the value parsers of the command line before line RG item K4i-2b (D5).

/// An integer argument with a constraint: after the conversion, NCBI's constraint reads the
/// value again with `NStr::StringToDouble` (default flags, `ncbi_str::string_to_double`), which
/// does not read a `0x` prefix without digits (the only form that `ncbi_integer` reads and
/// `strtod` does not end at), so the value is illegal.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_constrained_integer(value: &str) -> Result<i32, String> {
    let n = ncbi_integer(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    Ok(n)
}

/// Whether the constraint's `NStr::StringToDouble` of `value` (default flags) reads `value`
/// (`ncbi_constrained_integer`).
fn reads_as_double(value: &str) -> bool {
    crate::ncbi_str::string_to_double(value.as_bytes()).is_ok()
}

/// An integer argument of `CArg_Integer` with NCBI's lower bound (`at_least`,
/// `CArgAllowValuesGreaterThanOrEqual`).
fn ncbi_integer_at_least(value: &str, at_least: i32) -> Result<i32, String> {
    let n = ncbi_constrained_integer(value)?;
    if n < at_least {
        return Err(format!("expected an integer >= {at_least}"));
    }
    Ok(n)
}

/// A non-negative `Int8` argument (`CArg_Int8` with `CArgAllowValuesGreaterThanOrEqual` of 0):
/// `-searchsp`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
/// The constraint reads the value as `ncbi_constrained_integer`'s does.
pub(crate) fn nonnegative_ncbi_int8(value: &str) -> Result<i64, String> {
    let n = ncbi_int8(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    if n < 0 {
        return Err("expected an integer >= 0".into());
    }
    Ok(n)
}

/// A non-negative integer argument (`CArg_Integer` with `CArgAllowValuesGreaterThanOrEqual` of 0),
/// such as `-culling_limit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3302
fn nonnegative_ncbi_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 0)
}

/// A `-num_descriptions` or `-num_alignments` value: 0 or more, an `Int4` in NCBI's
/// integer forms (the sign and hexadecimal forms of `nonnegative_ncbi_integer`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685
fn report_count(value: &str) -> Result<usize, String> {
    nonnegative_ncbi_integer(value).map(|n| n as usize)
}

/// A BLASTN count of 1 or more (`-num_threads`, `-max_target_seqs`, `-max_hsps`), as
/// NCBI's arguments (c++/src/algo/blast/blastinput/blast_args.cpp:203-207,2731-2732,3162-3163).
fn blastn_count(value: &str) -> Result<usize, String> {
    ncbi_positive_integer(value).map(|n| n as usize)
}

/// `blastn_count` as an `int`, the type of BLASTX's fields.
fn ncbi_positive_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 1)
}

/// A percentage, in 0..=100 (NaN is outside): `-perc_identity`, and
/// `-qcov_hsp_perc` of every program (the `*_query_coverage` parsers).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201
/// The constraint reads the value with the argument's own conversion (`ncbi_double`), so it
/// tests the argument's number (`1e` is 1, `0x10` is 16).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
fn percentage(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(0.0..=100.0).contains(&n) {
        return Err("expected a percentage in 0..=100".into());
    }
    Ok(n)
}

/// A `-best_hit_overhang` or `-best_hit_score_edge` value: above 0 and below 0.5 (NaN is
/// outside), of every program that has the options.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/hspfilter_besthit.h:80-91
/// The constraint reads the value with the argument's own conversion (`ncbi_double`); NCBI
/// BLAST+ 2.17.0 refuses 0 and 0.5 with "Illegal value, expected (>0 and <0.5)" (line NP O6a's
/// probe).
fn best_hit_value(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(n > 0.0 && n < 0.5) {
        return Err("expected a number above 0 and below 0.5".into());
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:995-1006
fn genetic_code(value: &str) -> Result<u8, String> {
    let n = genetic_code_number(value).ok_or("invalid genetic code")?;
    if !crate::ncbi_app::args::is_genetic_code(i32::from(n)) {
        return Err("unsupported genetic code".into());
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-51
/// A BLASTN output file: a name shorter than 256 bytes, as NCBI's constraint; `-` is
/// standard output, and a file that cannot be created is reported when it is opened
/// (`blastn/blast_engine/run.rs`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
fn output_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().try_map(|value| {
        if crate::ncbi_app::args::file_name_length(value.as_encoded_bytes()) >= 256 {
            return Err("Illegal value, expected file name length < 256".to_string());
        }
        Ok(PathBuf::from(value))
    })
}

/// A `-threshold` value: NCBI's constraint reads it again with `NStr::StringToDouble` and its
/// default flags (`ncbi_str::string_to_double`: glibc's `strtod` alone, so a value that only
/// `StringToDoublePosix` reads, such as `1e`, makes it throw), and requires at least 0 (a NaN
/// fails it); a value that the constraint reads is the same number as the argument's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_threshold(value: &str) -> Result<f64, String> {
    let threshold = ncbi_double(value)?;
    // The check that `value` is at least `m_MinValue` (a NaN fails it).
    let allowed =
        crate::ncbi_str::string_to_double(value.as_bytes()).is_ok_and(|checked| checked >= 0.0);
    if !allowed {
        return Err("Illegal value, expected greater or equal to 0".into());
    }
    Ok(threshold)
}

/// A `-mt_mode` value: an integer from 0 to 2 (`CArgAllowValuesBetween` from 0 to 2, with `true`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3177
fn ncbi_mt_mode(value: &str) -> Result<u8, String> {
    let mode = ncbi_integer_at_least(value, 0)?;
    u8::try_from(mode)
        .ok()
        .filter(|mode| *mode <= 2)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 2".to_string())
}

/// A protein `-word_size` (`CArgAllowValuesGreaterThanOrEqual` of 2 for a protein query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170
fn protein_word_size(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 2)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
/// `NStr::StringToInt` of `value` in `CArgAllowGeneticCodeInteger` is
/// `value_parsers::genetic_code_number` (LOSAT's "invalid" message for the values outside 0..=255).
fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = genetic_code_number(value).ok_or("invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

/// A `-sorthits` or `-sorthsps` value: an integer from 0 to 4 (`CArgAllowValuesBetween` from 0 to 4,
/// with `true`, cited at `FormattingArgs`).
fn formatting_sort_option(value: &str) -> Result<u8, String> {
    let option = ncbi_integer_at_least(value, 0)?;
    u8::try_from(option)
        .ok()
        .filter(|option| *option <= 4)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 4".to_string())
}

/// A `-line_length` value: an integer of 1 or more (`CArgAllowValuesGreaterThanOrEqual` of 1,
/// cited at `FormattingArgs`).
fn formatting_line_length(value: &str) -> Result<usize, String> {
    let length = ncbi_integer_at_least(value, 1)?;
    usize::try_from(length).map_err(|_| "Illegal value: expected an integer >= 1".to_string())
}
