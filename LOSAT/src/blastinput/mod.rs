//! BLAST Input Handling
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/
//!
//! This module provides input handling for BLAST searches
//! including argument parsing and FASTA input processing.
//! The arguments of each program are `algorithm::<program>::args`.

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp and c++/src/algo/blast/blastinput/blast_args.cpp: the application layer.
pub mod app;
// NCBI objtools/readers/aln_reader.cpp, alnread.cpp, aln_formatguess.cpp: the alignment reader
// of psiblast's `-in_msa`.
pub(crate) mod aln_reader;
pub mod compressed;
pub(crate) mod title_parser;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp (s_GetTaxIDList): the taxid options of a database search.
pub mod db_filter;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2342-2356, c++/src/algo/blast/api/uniform_search.cpp: the database masks of a
// database search (`-db_soft_mask`, `-db_hard_mask`).
pub mod db_mask;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp, c++/src/objtools/readers/fasta.cpp, c++/src/util/line_reader.cpp:
// the FASTA input of BLASTN, TBLASTX, TBLASTN and BLASTP.
pub mod fasta_reader;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp: the input files of the command line as NCBI opens them.
pub mod input_files;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp, c++/src/corelib/metareg.cpp: the application layer's environment and registry.
pub mod ncbi_environment;
// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp (IRWRegistry::x_Read): the reader of the registry files.
pub(crate) mod ncbi_registry;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c BLAST_ValidateOptions: the protein option checks of blastp, tblastn
// and psiblast.
pub mod protein_options_check;
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp: psiblast's PSSM in ASN.1 text (`-out_pssm`, `-in_pssm`).
pub mod pssm_text;
pub mod query_batch;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978: the hit list size and the numbers of descriptions
// and alignments of a run.
pub mod report_counts;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:145-179 and c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-460: the
// `-query_loc` and `-subject_loc` ranges.
pub mod seq_range;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
pub mod value_parsers;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520 (`GetCompressedPlusStrand`): EXPERIMENT (LOSAT_X_TBNSTORE,
// LOSAT_X_TBXSTORE), the records of `-subject` packed as read.
pub(crate) mod x_packed_subjects;

#[cfg(test)]
mod damaged_input_tests;
