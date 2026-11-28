//! The program of an archive's search as blast_formatter's formatter names it.

use crate::cli::NativeError;
use crate::report::archive::read::options::{program_name_to_enum, EProgram};

use super::refusals::Refusal;

/// The program type of a search's options handle (the factory handle of the program that the
/// option builder computes and adjusts), as the formatter reads it: its name
/// (`Blast_ProgramNameFromType`), whether its subjects are nucleotide (`Blast_SubjectIsNucleotide`)
/// and whether it is `eBlastTypeBlastn` (the check of SAM).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ProgramType {
    pub(crate) name: &'static str,
    pub(crate) subject_nucleotide: bool,
    pub(crate) blastn: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:124-200
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:707-751
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:319-357
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:72-89
// The factory makes a nucleotide handle of type blastn for blastn, megablast, dc-megablast and
// vecscreen; the RPS, RPS-tblastn, PSI-BLAST, PSI-tblastn, PHI-BLAST and mapping handles for those
// programs; and a DELTA-BLAST handle whose type is PSI-BLAST's. The type's name is the one the
// core gives it. The subjects' kind is the program type's: the PSSM types are blastp's and
// tblastn's with a PSSM query, blastp's and blastx's with a PSSM subject; PHI-BLAST's are its
// query's kind; the mapper's are nucleotide.
/// The program type of the options handle made for `program`.
pub(crate) fn program_type(program: EProgram) -> ProgramType {
    let (name, subject_nucleotide) = match program {
        EProgram::Blastn | EProgram::Megablast | EProgram::DiscMegablast | EProgram::VecScreen => {
            ("blastn", true)
        }
        EProgram::Blastp => ("blastp", false),
        EProgram::Blastx => ("blastx", false),
        EProgram::Tblastn => ("tblastn", true),
        EProgram::Tblastx => ("tblastx", true),
        EProgram::RpsBlast => ("rpsblast", false),
        EProgram::RpsTblastn => ("rpstblastn", false),
        EProgram::PsiBlast | EProgram::DeltaBlast => ("psiblast", false),
        EProgram::PsiTblastn => ("psitblastn", true),
        EProgram::PhiBlastp => ("phiblastp", false),
        EProgram::PhiBlastn => ("phiblastn", true),
        EProgram::Mapper => ("mapper", true),
    };
    ProgramType {
        name,
        subject_nucleotide,
        blastn: name == "blastn",
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:226-231
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:185-187
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:666-673
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:325-400
// The formatter keeps the handle type's name as its program, DELTA-BLAST's as `deltablast`, and
// maps it back to a program for each Seq-annot it writes: the Blast Type is that program's task
// name and number. The PHI-BLAST names are no program of the map (NCBI's `CBlastException`, which
// `CATCH_ALL` reports as an engine error, exit status 3). The mapper's name maps back to the
// mapper, whose archives NCBI's programs do not write: refused.
/// The label and value of the Seq-annot's "Blast Type" of an archive of `program`; the error
/// the run ends with where NCBI's formatter has none.
pub(crate) fn blast_type(program: EProgram) -> anyhow::Result<(&'static str, i32)> {
    match program_name_to_enum(formatter_program(program).as_bytes()) {
        Ok(EProgram::Mapper) => Err(Refusal::Program(EProgram::Mapper.task_name()).into()),
        Ok(program) => Ok((program.task_name(), program.number())),
        Err(exception) => Err(NativeError::exception(&exception).into()),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:185-187
/// The formatter's program (`m_Program`) of an archive of `program`: its handle type's name,
/// `deltablast` for DELTA-BLAST (whose handle type is PSI-BLAST's). The Seq-annot's Blast Type
/// (`blast_type`) and the tabular formats' version line (`# DELTABLAST 2.17.0+`) are of it.
pub(crate) fn formatter_program(program: EProgram) -> &'static str {
    match program {
        EProgram::DeltaBlast => "deltablast",
        _ => program_type(program).name,
    }
}

#[cfg(test)]
#[path = "search_program_tests.rs"]
mod search_program_tests;
