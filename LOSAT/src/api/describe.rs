//! The options of a program, from the engine's own argument definitions.
//!
//! NCBI builds each program's options once, as `CArgDescriptions` (`SetArgumentDescriptions`),
//! and prints them from there (`PrintUsage`); LOSAT's argument tables are NCBI's (generated from
//! each program's `-xmlhelp` XML, `ncbi_app::args`), which the command line's reader reads, and
//! `describe` takes every part of an option from them (line RG DECISIONS K4i-3-1): the options
//! and their order (the `-help`'s), their value kinds, defaults and allowed values, and their help
//! (NCBI's comment, in LOSAT's words), so a caller shows exactly what the command line reads and
//! holds no table of its own (LOSAT Web ABI v2 *describe*, docs/web/abi_v2.md §9, reads the same
//! rows).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
#![warn(missing_docs, unused)]

use super::program::Program;
use super::search::{alignment_owned, file_options, API_OWNED};
use crate::ncbi_app::args::Kind;

// NCBI reference (598d8ae6): c++/include/corelib/ncbiargs.hpp:588-598
/// What an option's value is: its type in NCBI's argument tables.
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
    /// The help text of the command line: NCBI's comment of the option, which its `-help`
    /// prints, in LOSAT's words.
    pub help: String,
    /// Whether the flag takes a value (`-evalue 1e-5`) or stands alone (`-ungapped`).
    pub takes_value: bool,
    /// What the value is.
    pub value_kind: ValueKind,
    /// The default of NCBI's argument reader (`AddDefaultKey`, the XML's `<default>`), if the
    /// option has one. An option without one (a flag, or a key whose default NCBI's option
    /// handlers set, "given or not" to the reader) has none here even when its help text names
    /// NCBI's default.
    pub default: Option<String>,
    /// The values that NCBI's argument reader accepts, when its constraint is a set (`-task`,
    /// `-strand`, `-template_type`, `-template_length`); empty otherwise, a flag's included.
    pub choices: Vec<String>,
}

/// The options of a program and the genetic codes its parser accepts.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Description {
    /// The described program.
    pub program: Program,
    /// The options that an argv of the library API may hold, in the order of the command
    /// line's `-help`: every option of the program but those the API sets itself (`-query`,
    /// `-subject`, `-db`, `-in_pssm`, `-out`, `-outfmt`; `search::API_OWNED`; psiblast's
    /// `-in_msa`, `search::alignment_owned`) and the search strategy file, which it refuses
    /// (`search::file_options`). psiblast's `-msa_master_idx` and `-ignore_msa_master` are
    /// listed: they choose the master of the alignment that `Inputs::msa` gives.
    pub parameters: Vec<Parameter>,
    /// The values that `-query_gencode` accepts; empty for a program without the option.
    pub query_gencodes: Vec<u32>,
    /// The values that `-db_gencode` accepts; empty for a program without the option.
    pub subject_gencodes: Vec<u32>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:998-1012 (cited at
// `ncbi_app::constraint`)
// NCBI's key of a genetic code allows the codes of `CArgAllowGeneticCodeInteger`.
/// The genetic codes that `program`'s key `name` allows: its constraint in the argument tables
/// (`ncbi_app::args::genetic_codes`, NCBI's codes and LOSAT's 32 of tblastn's `-db_gencode`).
fn accepted_gencodes(program: Program, name: &str) -> Vec<u32> {
    crate::ncbi_app::args::genetic_codes(program.name(), name).unwrap_or_default()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2387-2400
/// The kind of a value of NCBI's type `kind` (`CArgDescriptions::GetTypeName`: an `Int8` or an
/// `Integer` is an integer, a `File_In` or a `File_Out` a file's name).
fn value_kind(kind: Kind) -> ValueKind {
    match kind {
        Kind::Flag => ValueKind::Flag,
        Kind::Boolean => ValueKind::Boolean,
        Kind::Int8 | Kind::Integer => ValueKind::Integer,
        Kind::Double => ValueKind::Real,
        Kind::InputFile | Kind::OutputFile => ValueKind::Path,
        Kind::Text => ValueKind::Text,
    }
}

/// Every option of `program` that LOSAT's descriptions list, in the order of the command line's
/// help (`ncbi_app::args::described`): `describe`'s, and those the library API sets itself. LOSAT
/// Web's *describe* lists these rows but those that its adapter sets (`web_adapter`).
pub(crate) fn parameters(program: Program) -> impl Iterator<Item = Parameter> {
    let tool =
        crate::cli::reader::tool(program.name()).expect("every API program is on NCBI's reader");
    crate::ncbi_app::args::described(tool).map(|argument| Parameter {
        flag: format!("-{}", argument.name),
        help: crate::wording::losat_text(argument.description).into_owned(),
        takes_value: argument.kind != Kind::Flag,
        value_kind: value_kind(argument.kind),
        default: argument.default.map(str::to_string),
        choices: argument.choices,
    })
}

/// The options of `program`, read from the engine's argument definitions.
pub fn describe(program: Program) -> Description {
    let owned = |flag: &str| {
        let name = flag.trim_start_matches('-');
        API_OWNED.contains(&name)
            || file_options(program).contains(&name)
            || alignment_owned(program).contains(&name)
    };
    Description {
        program,
        parameters: parameters(program)
            .filter(|parameter| !owned(&parameter.flag))
            .collect(),
        query_gencodes: accepted_gencodes(program, "query_gencode"),
        subject_gencodes: accepted_gencodes(program, "db_gencode"),
    }
}

#[cfg(test)]
#[path = "describe_tests.rs"]
mod tests;
