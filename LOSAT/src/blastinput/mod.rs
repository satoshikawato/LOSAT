//! BLAST Input Handling
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/
//!
//! This module provides input handling for BLAST searches
//! including argument parsing and FASTA input processing.
//!
//! # Structure
//!
//! - `blast_args` - Common argument handling
//! - `blastn_args` - BLASTN-specific arguments
//! - `blastp_args` - BLASTP-specific arguments
//! - `tblastx_args` - TBLASTX-specific arguments

// NCBI app/blast/blast_app_util.hpp and blastinput/blast_args.cpp: the application layer.
pub mod app;
// ABI v1's and the Web adapter's (until S10) checks of `bio` FASTA records.
pub mod bio_checks;
pub mod blast_args;
pub mod blastn_args;
pub mod blastp_args;
// NCBI blastinput/blast_fasta_input.cpp, objtools/readers/fasta.cpp, util/line_reader.cpp:
// the FASTA input of BLASTN, TBLASTX, TBLASTN and BLASTP.
pub mod fasta_reader;
// NCBI corelib/ncbiargs.cpp: the input files of the command line as NCBI opens them.
pub mod input_files;
// NCBI corelib/ncbiapp.cpp, metareg.cpp: the application layer's environment and registry.
pub mod ncbi_environment;
pub mod query_batch;
// NCBI blastinput/blast_input_aux.cpp:145-179 and blast_fasta_input.cpp:433-460: the
// `-query_loc` and `-subject_loc` ranges.
pub mod seq_range;
pub mod tblastx_args;

// NCBI blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
pub mod value_parsers;
