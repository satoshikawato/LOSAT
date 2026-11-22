//! The program of an archive's search as blast_formatter's formatter names it.

use crate::report::archive::read::options::EProgram;

/// The program of a search's options handle (the factory handle of the program that the option
/// builder computes and adjusts), as blast_formatter's formatter names it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SearchProgram {
    /// A program whose type the formatter names (`Blast_ProgramNameFromType`, then
    /// `ProgramNameToEnum`).
    Named(ProgramType),
    /// A program whose Seq-annot NCBI's formatter does not name as the others (PHI-BLAST, whose
    /// type's name `ProgramNameToEnum` does not know; the mapper of magicblast): its task name
    /// (`EProgramToTaskName`). LOSAT's blast_formatter refuses its archives.
    Other(&'static str),
}

/// The program type of a search as blast_formatter's formatter writes it: the label and value of
/// the Seq-annot's "Blast Type" (the task name and number of the program that
/// `ProgramNameToEnum` gives for the handle's type name), whether its subjects are nucleotide
/// (`Blast_SubjectIsNucleotide`) and whether it is `eBlastTypeBlastn` (the check of SAM).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ProgramType {
    pub(crate) task: &'static str,
    pub(crate) number: i32,
    pub(crate) subject_nucleotide: bool,
    pub(crate) blastn: bool,
}

impl ProgramType {
    const fn new(task: &'static str, number: i32, subject_nucleotide: bool) -> Self {
        ProgramType {
            task,
            number,
            subject_nucleotide,
            blastn: number == 1,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:124-200
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:676-699
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:707-751
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:757-811
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:319-345
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:56-76
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:72-89
// The factory makes a nucleotide handle of type blastn for blastn, megablast, dc-megablast and
// vecscreen; the RPS, RPS-tblastn, PSI-BLAST and PSI-tblastn handles for those programs; and a
// DELTA-BLAST handle whose type is PSI-BLAST's. The formatter names the handle's type, maps that
// name back to a program and writes the program's task name and number (the numbers of
// `EProgram`). The subjects' kind is the program type's: the PSSM types are blastp's and
// tblastn's with a PSSM query, blastp's and blastx's with a PSSM subject.
/// The program of the options handle made for `program`, as the formatter names it.
pub(crate) fn search_program(program: EProgram) -> SearchProgram {
    SearchProgram::Named(match program {
        EProgram::Blastn | EProgram::Megablast | EProgram::DiscMegablast | EProgram::VecScreen => {
            ProgramType::new("blastn", 1, true)
        }
        EProgram::Blastp => ProgramType::new("blastp", 2, false),
        EProgram::Blastx => ProgramType::new("blastx", 3, false),
        EProgram::Tblastn => ProgramType::new("tblastn", 4, true),
        EProgram::Tblastx => ProgramType::new("tblastx", 5, true),
        EProgram::RpsBlast => ProgramType::new("rpsblast", 6, false),
        EProgram::RpsTblastn => ProgramType::new("rpstblastn", 7, false),
        EProgram::PsiBlast | EProgram::DeltaBlast => ProgramType::new("psiblast", 10, false),
        EProgram::PsiTblastn => ProgramType::new("psitblastn", 11, true),
        EProgram::PhiBlastp | EProgram::PhiBlastn | EProgram::Mapper => {
            return SearchProgram::Other(program.task_name())
        }
    })
}
