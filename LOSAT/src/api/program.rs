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

/// A search program of the library API. `makeblastdb`, `blastdbcmd` and `dustmasker` are not
/// searches: they have functions of their own (`run_makeblastdb`, `run_blastdbcmd`,
/// `run_dustmasker`).
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
    /// Protein query against protein subjects: LOSAT runs psiblast's first round, which searches
    /// each query with blastp's engine, before any PSSM exists.
    Psiblast,
}

impl Program {
    /// Every program, in the order of the command line's help.
    pub const ALL: [Program; 6] = [
        Program::Blastx,
        Program::Blastn,
        Program::Blastp,
        Program::Tblastx,
        Program::Tblastn,
        Program::Psiblast,
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
            Program::Psiblast => "psiblast",
        }
    }

    /// The program of a command-line name. A program that the command line runs and the API
    /// does not search with (`makeblastdb`, `makedb`, `blastdbcmd`, `dustmasker`) is named in
    /// the error.
    pub fn parse(name: &str) -> Result<Self, Error> {
        match name {
            "blastn" => Ok(Program::Blastn),
            "blastp" => Ok(Program::Blastp),
            "blastx" => Ok(Program::Blastx),
            "tblastn" => Ok(Program::Tblastn),
            "tblastx" => Ok(Program::Tblastx),
            "psiblast" => Ok(Program::Psiblast),
            "makeblastdb" | "makedb" => Err(Error::Argv(format!(
                "{name} is not a search program: run_makeblastdb runs it"
            ))),
            "blastdbcmd" => Err(Error::Argv(
                "blastdbcmd is not a search program: run_blastdbcmd runs it".to_string(),
            )),
            "dustmasker" => Err(Error::Argv(
                "dustmasker is not a search program: run_dustmasker runs it".to_string(),
            )),
            // rpsblast (series PR, item PS11a) runs on the command line; the API has no program
            // for it yet (its results and options are those of a later item).
            "rpsblast" => Err(Error::Argv(
                "rpsblast is not available through the library API".to_string(),
            )),
            // rpstblastn (series PR, item PS12a) likewise.
            "rpstblastn" => Err(Error::Argv(
                "rpstblastn is not available through the library API".to_string(),
            )),
            other => Err(Error::Argv(format!("unknown program '{other}'"))),
        }
    }

    /// The engine's program, for what depends on the program in the reading of `-outfmt`.
    pub(crate) fn engine(self) -> crate::blastinput::app::Program {
        use crate::blastinput::app::Program as Engine;
        match self {
            Program::Blastn => Engine::Blastn,
            // psiblast's first round is blastp's search (`algorithm::psiblast`).
            Program::Blastp | Program::Psiblast => Engine::Blastp,
            Program::Blastx => Engine::Blastx,
            Program::Tblastn => Engine::Tblastn,
            Program::Tblastx => Engine::Tblastx,
        }
    }

    /// The upper-case name that the engine's messages and readers use (`BLASTN`, ...).
    pub(crate) fn reader_name(self) -> &'static str {
        match self {
            Program::Psiblast => "PSIBLAST",
            program => program.engine().name(),
        }
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
    fn names_round_trip_and_the_tools_are_named() {
        for program in Program::ALL {
            assert_eq!(Program::parse(program.name()), Ok(program));
            assert_eq!(program.to_string(), program.name());
            assert_eq!(program.name().parse::<Program>(), Ok(program));
        }
        for tool in ["makeblastdb", "makedb", "blastdbcmd"] {
            let error = Program::parse(tool).unwrap_err();
            assert!(error.message().contains("not a search program"), "{error}");
        }
        // rpsblast and rpstblastn run on the command line only (line PS, PS11a, PS12a).
        for program in ["rpsblast", "rpstblastn"] {
            let error = Program::parse(program).unwrap_err();
            assert!(
                error.message().contains(&format!(
                    "{program} is not available through the library API"
                )),
                "{error}"
            );
        }
        assert_eq!(
            Program::parse("nosuch"),
            Err(Error::Argv("unknown program 'nosuch'".to_string()))
        );
    }
}
