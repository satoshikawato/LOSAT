//! PHI-BLAST (psiblast `-phi_pattern`): the search of the alignments that contain a pattern's
//! occurrences in the query and in the subject. Cold, scalar code, as NCBI's PHI code: the pattern
//! compiler, the pattern scans, the statistics, the gapped alignments around the pattern, the
//! engine of round 1 and the pattern file reader.

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1382-1411
pub(crate) mod file;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:141-680
pub(crate) mod pattern;
