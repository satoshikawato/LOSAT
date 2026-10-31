//! BLAST Input Handling
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/
//!
//! This module provides input handling for BLAST searches
//! including argument parsing and FASTA input processing.
//! The arguments of each program are `algorithm::<program>::args`.

// NCBI app/blast/blast_app_util.hpp and blastinput/blast_args.cpp: the application layer.
pub mod app;
// NCBI blastinput/blast_args.cpp (s_GetTaxIDList): the taxid options of a database search.
pub mod db_filter;
// NCBI blastinput/blast_fasta_input.cpp, objtools/readers/fasta.cpp, util/line_reader.cpp:
// the FASTA input of BLASTN, TBLASTX, TBLASTN and BLASTP.
pub mod fasta_reader;
// NCBI corelib/ncbiargs.cpp: the input files of the command line as NCBI opens them.
pub mod input_files;
// NCBI corelib/ncbiapp.cpp, metareg.cpp: the application layer's environment and registry.
pub mod ncbi_environment;
pub mod query_batch;
// NCBI blastinput/blast_args.cpp:2894-2978: the hit list size and the numbers of descriptions
// and alignments of a run.
pub mod report_counts;
// NCBI blastinput/blast_input_aux.cpp:145-179 and blast_fasta_input.cpp:433-460: the
// `-query_loc` and `-subject_loc` ranges.
pub mod seq_range;

// NCBI blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
pub mod value_parsers;
// NCBI api/blast_objmgr_tools.cpp:515-520 (`GetCompressedPlusStrand`): EXPERIMENT (LOSAT_X_TBNSTORE,
// LOSAT_X_TBXSTORE), the records of `-subject` packed as read.
pub(crate) mod x_packed_subjects;
