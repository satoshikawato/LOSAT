//! The options of a program, from the engine's own argument definitions.
//!
//! NCBI builds each program's options once, as `CArgDescriptions` (`SetArgumentDescriptions`),
//! and prints them from there (`PrintUsage`); LOSAT's clap definitions are that table, and
//! `describe` reads it, so a caller shows exactly the options, defaults and allowed values of the
//! command line and holds no table of its own (LOSAT Web ABI v2 *describe*, docs/web/abi_v2.md
//! §9, reads the same definitions).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
//! ```c++
//! string& CArgDescriptions::PrintUsage(string& str, bool detailed) const
//! {
//!     CPrintUsage x(*this);
//!     list<string> arr;
//!
//!     // SYNOPSIS
//!     arr.push_back("USAGE");
//!     x.AddSynopsis(arr, m_UsageName, "    ");
//!
//!     // DESCRIPTION
//!     arr.push_back(kEmptyStr);
//!     x.AddDescription(arr, detailed);
//! ```
#![warn(missing_docs, unused)]

use clap::CommandFactory;

use super::program::Program;
use super::search::API_OWNED;
use crate::cli::Cli;

// NCBI reference (598d8ae6): c++/include/corelib/ncbiargs.hpp:588-598
// ```c++
//     enum EType {
//         eString = 0, ///< An arbitrary string
//         eBoolean,    ///< {'true', 't', 'false', 'f'},  case-insensitive
//         eInt8,       ///< Convertible into an integer number (Int8 only)
//         eInteger,    ///< Convertible into an integer number (int or Int8)
//         eIntId,      ///< Convertible to TIntId (int or Int8 depending on NCBI_INT8_GI
//         eDouble,     ///< Convertible into a floating point number (double)
//         eInputFile,  ///< Name of file (must exist and be readable)
//         eOutputFile, ///< Name of file (must be writable)
//         eIOFile,     ///< Name of file (must be writable)
//         eDirectory,  ///< Name of file directory
// ```
/// What an option's value is, from the type that the command line's parser makes of it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ValueKind {
    /// No value: the flag stands alone.
    Flag,
    /// `true` or `false`, as the command line spells them.
    Boolean,
    /// An integer.
    Integer,
    /// A floating-point number.
    Real,
    /// A file's name.
    Path,
    /// Any other text; `Parameter::choices` lists the values of an option that has a list.
    Text,
}

/// One option of a program, as the command line defines it.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Parameter {
    /// The flag with its single dash, as on the command line (`-evalue`).
    pub flag: String,
    /// The help text of the command line.
    pub help: String,
    /// Whether the flag takes a value (`-evalue 1e-5`) or stands alone (`-ungapped`).
    pub takes_value: bool,
    /// What the value is.
    pub value_kind: ValueKind,
    /// The default value that the parser fills in when the flag is absent, if it has one. An
    /// option without a parser default (an `Option` of the arguments, "given or not" to NCBI)
    /// has none here even when its help text names NCBI's default.
    pub default: Option<String>,
    /// The values the parser accepts, when it lists them (`-task` of blastp); empty otherwise.
    pub choices: Vec<String>,
}

/// The options of a program and the genetic codes its parser accepts.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Description {
    /// The described program.
    pub program: Program,
    /// The options that an argv of the library API may hold, in the order of the command
    /// line's help: every option of the program but those the API sets itself (`-query`,
    /// `-subject`, `-db`, `-out`, `-outfmt`; `search::API_OWNED`).
    pub parameters: Vec<Parameter>,
    /// The values that `-query_gencode` accepts; empty for a program without the option.
    pub query_gencodes: Vec<u32>,
    /// The values that `-db_gencode` accepts; empty for a program without the option.
    pub subject_gencodes: Vec<u32>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1021-1024
// ```c++
//         arg_desc.AddDefaultKey(kArgQueryGeneticCode, "int_value",
//                                "Genetic code to use to translate query (see https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes for details)\n",
//                                CArgDescriptions::eInteger,
//                                NStr::IntToString(BLAST_GENETIC_CODE));
// ```
// NCBI's argument takes any integer and the genetic-code table decides later; LOSAT's parser
// decides at once, so the accepted values are those the parser takes.
/// Genetic codes are small integers; every accepted one is below this bound.
const GENETIC_CODE_PROBE: u32 = 64;

/// The genetic codes that the parser of `program` accepts for `flag`.
fn accepted_gencodes(program: Program, flag: &str) -> Vec<u32> {
    (0..GENETIC_CODE_PROBE)
        .filter(|code| {
            let code = code.to_string();
            let argv = [
                "losat",
                program.name(),
                "-query",
                "q",
                "-subject",
                "s",
                flag,
                &code,
            ];
            crate::cli::try_parse_from::<Cli, _, _>(argv).is_ok()
        })
        .collect()
}

/// The kind of the value of `arg`, from the type its parser makes.
fn value_kind(arg: &clap::Arg) -> ValueKind {
    use std::any::TypeId;
    if !arg.get_action().takes_values() {
        return ValueKind::Flag;
    }
    let parser = arg.get_value_parser();
    let made = parser.type_id();
    let is_one_of = |types: &[TypeId]| types.iter().any(|candidate| made == *candidate);
    if is_one_of(&[TypeId::of::<bool>()]) {
        ValueKind::Boolean
    } else if is_one_of(&[
        TypeId::of::<u8>(),
        TypeId::of::<u32>(),
        TypeId::of::<u64>(),
        TypeId::of::<usize>(),
        TypeId::of::<i32>(),
        TypeId::of::<i64>(),
    ]) {
        ValueKind::Integer
    } else if is_one_of(&[TypeId::of::<f32>(), TypeId::of::<f64>()]) {
        ValueKind::Real
    } else if is_one_of(&[TypeId::of::<std::path::PathBuf>()]) {
        ValueKind::Path
    } else {
        ValueKind::Text
    }
}

/// The options of `program`, read from the engine's argument definitions.
pub fn describe(program: Program) -> Description {
    let root = Cli::command();
    let command = root
        .find_subcommand(program.name())
        .expect("every API program is a command of the CLI");
    let mut parameters = Vec::new();
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else { continue };
        if matches!(long, "help" | "version") || API_OWNED.contains(&long) {
            continue;
        }
        let defaults: Vec<String> = arg
            .get_default_values()
            .iter()
            .map(|value| value.to_string_lossy().into_owned())
            .collect();
        parameters.push(Parameter {
            flag: format!("-{long}"),
            help: arg.get_help().map(ToString::to_string).unwrap_or_default(),
            takes_value: arg.get_action().takes_values(),
            value_kind: value_kind(arg),
            default: (!defaults.is_empty()).then(|| defaults.join(" ")),
            choices: arg
                .get_possible_values()
                .iter()
                .map(|value| value.get_name().to_string())
                .collect(),
        });
    }
    let has = |long: &str| {
        command
            .get_arguments()
            .any(|arg| arg.get_long() == Some(long))
    };
    Description {
        program,
        parameters,
        query_gencodes: if has("query_gencode") {
            accepted_gencodes(program, "-query_gencode")
        } else {
            Vec::new()
        },
        subject_gencodes: if has("db_gencode") {
            accepted_gencodes(program, "-db_gencode")
        } else {
            Vec::new()
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn api_owned_options_are_not_described_and_search_options_are() {
        for program in Program::ALL {
            let description = describe(program);
            assert_eq!(description.program, program);
            let flags: Vec<&str> = description
                .parameters
                .iter()
                .map(|parameter| parameter.flag.as_str())
                .collect();
            for owned in API_OWNED {
                assert!(
                    !flags.contains(&format!("-{owned}").as_str()),
                    "{program} lists -{owned}"
                );
            }
            for expected in ["-evalue", "-num_threads", "-max_target_seqs"] {
                assert!(flags.contains(&expected), "{program} lacks {expected}");
            }
            let evalue = description
                .parameters
                .iter()
                .find(|parameter| parameter.flag == "-evalue")
                .unwrap();
            assert!(evalue.takes_value);
            // The help texts are the command line's; some options have none.
            assert!(description
                .parameters
                .iter()
                .any(|parameter| !parameter.help.is_empty()));
        }
    }

    #[test]
    fn genetic_codes_are_listed_for_the_translated_programs_only() {
        let tblastx = describe(Program::Tblastx);
        assert!(tblastx.query_gencodes.contains(&1));
        assert!(tblastx.query_gencodes.contains(&11));
        assert!(!tblastx.query_gencodes.contains(&0));
        assert!(!tblastx.subject_gencodes.is_empty());
        let blastx = describe(Program::Blastx);
        assert!(!blastx.query_gencodes.is_empty());
        assert!(blastx.subject_gencodes.is_empty());
        let tblastn = describe(Program::Tblastn);
        assert!(tblastn.query_gencodes.is_empty());
        assert!(!tblastn.subject_gencodes.is_empty());
        for program in [Program::Blastn, Program::Blastp] {
            let description = describe(program);
            assert!(description.query_gencodes.is_empty(), "{program}");
            assert!(description.subject_gencodes.is_empty(), "{program}");
        }
    }
}
