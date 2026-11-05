//! The search programs of the library API.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
//! ```c++
//!     static const string kProgram("blastx");
//!     arg.Reset(new CProgramDescriptionArgs(kProgram,
//!                                   "Translated Query-Protein Subject BLAST"));
//! ```
//! Each NCBI application names its program once, in lower case; the API takes the same names.
#![warn(missing_docs, unused)]

use std::fmt;

use super::error::Error;

/// A search program of the library API: the five that have an in-memory entry (`run_local`).
/// `psiblast`, `makeblastdb` and `blastdbcmd` are command-line programs only
/// (docs/losat_library_surface_plan.md §8).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Program {
    /// Nucleotide query against nucleotide subjects (`megablast`, `dc-megablast`, `blastn`,
    /// `blastn-short`).
    Blastn,
    /// Protein query against protein subjects.
    Blastp,
    /// Translated nucleotide query against protein subjects.
    Blastx,
    /// Protein query against translated nucleotide subjects.
    Tblastn,
    /// Translated nucleotide query against translated nucleotide subjects.
    Tblastx,
}

impl Program {
    /// Every program, in the order of the command line's help.
    pub const ALL: [Program; 5] = [
        Program::Blastx,
        Program::Blastn,
        Program::Blastp,
        Program::Tblastx,
        Program::Tblastn,
    ];

    /// The command-line name of the program (`blastn`, ...): the first argv word of the
    /// command line and of LOSAT Web ABI v2.
    pub fn name(self) -> &'static str {
        match self {
            Program::Blastn => "blastn",
            Program::Blastp => "blastp",
            Program::Blastx => "blastx",
            Program::Tblastn => "tblastn",
            Program::Tblastx => "tblastx",
        }
    }

    /// The program of a command-line name. A program that the command line runs and the API
    /// does not (`psiblast`, `makeblastdb`, `makedb`, `blastdbcmd`) is named in the error.
    pub fn parse(name: &str) -> Result<Self, Error> {
        match name {
            "blastn" => Ok(Program::Blastn),
            "blastp" => Ok(Program::Blastp),
            "blastx" => Ok(Program::Blastx),
            "tblastn" => Ok(Program::Tblastn),
            "tblastx" => Ok(Program::Tblastx),
            "psiblast" | "makeblastdb" | "makedb" | "blastdbcmd" => Err(Error::Argv(format!(
                "{name} is available on the command line only, not through the library API"
            ))),
            other => Err(Error::Argv(format!("unknown program '{other}'"))),
        }
    }

    /// The engine's program, for what depends on the program in the reading of `-outfmt`.
    pub(crate) fn engine(self) -> crate::blastinput::app::Program {
        use crate::blastinput::app::Program as Engine;
        match self {
            Program::Blastn => Engine::Blastn,
            Program::Blastp => Engine::Blastp,
            Program::Blastx => Engine::Blastx,
            Program::Tblastn => Engine::Tblastn,
            Program::Tblastx => Engine::Tblastx,
        }
    }

    /// The upper-case name that the engine's messages and readers use (`BLASTN`, ...).
    pub(crate) fn reader_name(self) -> &'static str {
        self.engine().name()
    }
}

impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl std::str::FromStr for Program {
    type Err = Error;

    fn from_str(name: &str) -> Result<Self, Error> {
        Program::parse(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip_and_cli_only_programs_are_named() {
        for program in Program::ALL {
            assert_eq!(Program::parse(program.name()), Ok(program));
            assert_eq!(program.to_string(), program.name());
            assert_eq!(program.name().parse::<Program>(), Ok(program));
        }
        for cli_only in ["psiblast", "makeblastdb", "makedb", "blastdbcmd"] {
            let error = Program::parse(cli_only).unwrap_err();
            assert!(error.message().contains("command line only"), "{error}");
        }
        assert_eq!(
            Program::parse("nosuch"),
            Err(Error::Argv("unknown program 'nosuch'".to_string()))
        );
    }
}
