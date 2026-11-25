//! The search programs (blastn, blastp, blastx, phi_search, psiblast, pssm_search, rps_search, rpsblast,
//! rpstblastn, tblastn, tblastx) and
//! `common`, the code they share.
//!
//! The programs are layered, bottom to top:
//! `common < blastn < tblastx < blastp < tblastn < blastx`. A program uses only the programs
//! below it, and `common` uses none; psiblast and pssm_search join the order above blastp once
//! their cycle with blastp is cut. `LOSAT/tests/check_program_edges.py` (CI) fails on a reference
//! up the order and on a reference between two programs that `LOSAT/tests/program_edges.tsv` does
//! not list. Rows are added only by review, for an edge down the layering; remove a row when its
//! edge goes (fewer edges are welcome).

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
pub mod blastn;
pub mod blastp;
pub mod common;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:388-688
// PHI-BLAST (psiblast `-phi_pattern`): the pattern, its scans, statistics, alignments and engine.
pub mod phi_search;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:47-51
pub mod psiblast;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_program.c:46-47
// The search of a query given as a PSSM (psiblast `-in_pssm`): its kernels and stages.
pub mod pssm_search;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_program.c:73-74
// The search of a profile database (rpsblast): its scan, word finders, stages and engine.
pub mod rps_search;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:77-79
pub mod rpsblast;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:48-50
pub mod rpstblastn;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:45-62
pub mod tblastn;
pub mod tblastx;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-55
pub mod blastx;

#[cfg(test)]
mod mapped_row_range_tests;
