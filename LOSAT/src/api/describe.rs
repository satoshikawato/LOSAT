//! The options of a program, from the engine's own argument definitions.
//!
//! NCBI builds each program's options once, as `CArgDescriptions` (`SetArgumentDescriptions`),
//! and prints them from there (`PrintUsage`); LOSAT's argument tables are NCBI's (generated from
//! each program's `-xmlhelp` XML, `ncbi_app::args`), and `describe` takes the defaults and the
//! allowed values from them, the options, their help and value kinds from the clap definitions
//! that convert them, so a caller shows exactly what the command line reads and holds no table
//! of its own (LOSAT Web ABI v2 *describe*, docs/web/abi_v2.md §9, reads the same).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
#![warn(missing_docs, unused)]

use clap::CommandFactory;

use super::program::Program;
use super::search::{alignment_owned, file_options, API_OWNED};
use crate::cli::Cli;

// NCBI reference (598d8ae6): c++/include/corelib/ncbiargs.hpp:588-598
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
    /// line's help: every option of the program but those the API sets itself (`-query`,
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
        if matches!(long, "help" | "version")
            || API_OWNED.contains(&long)
            || file_options(program).contains(&long)
            || alignment_owned(program).contains(&long)
        {
            continue;
        }
        // The default and the values that NCBI's argument reader allows (its tables, line RG
        // DECISIONS K4i-2b-1): none for a flag.
        let name = program.name();
        parameters.push(Parameter {
            flag: format!("-{long}"),
            help: arg.get_help().map(ToString::to_string).unwrap_or_default(),
            takes_value: arg.get_action().takes_values(),
            value_kind: value_kind(arg),
            default: crate::ncbi_app::args::default_value(name, long).map(str::to_string),
            choices: crate::ncbi_app::args::choices(name, long).unwrap_or_default(),
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
            accepted_gencodes(program, "query_gencode")
        } else {
            Vec::new()
        },
        subject_gencodes: if has("db_gencode") {
            accepted_gencodes(program, "db_gencode")
        } else {
            Vec::new()
        },
    }
}

#[cfg(test)]
#[path = "describe_tests.rs"]
mod tests;
