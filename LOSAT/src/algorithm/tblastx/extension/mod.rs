//! Extension algorithms for TBLASTX
//!
//! This module implements the two-hit extension for protein sequences, following NCBI BLAST's
//! approach with X-drop termination.
//!
//! Sequences are encoded in NCBISTDAA (0-27), and scoring uses BLOSUM62 via
//! the ncbistdaa_to_blosum62 conversion table.
//!
//! # Module Structure
//!
//! - `two_hit` - Two-hit extension algorithm (NCBI style)

// The one-hit extension (`-window_size 0`).
mod one_hit;
// The query batch the walks go on into past a context's ends.
mod query_batch;
mod two_hit;
// The table walks of the two-hit extension.
mod x_two_hit_fast;

pub(crate) use one_hit::{
    extend_hit_one_hit_alone, extend_hit_one_hit_outside, extend_hit_one_hit_resume, OneHit,
    OneHitAlone,
};
pub(crate) use query_batch::{Again, InBatch, QueryBatch};
pub(crate) use two_hit::extend_hit_two_hit_fitting;
pub(crate) use two_hit::extend_hit_two_hit_outside;
pub(crate) use two_hit::extend_hit_two_hit_with;
// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the walks the sentinels stop, and the checks that choose them.
pub(crate) use two_hit::extend_hit_two_hit_sent;
pub(crate) use x_two_hit_fast::{
    matrix_walk_table, sent_contexts as x_tbx_sent_contexts, sent_frame as x_tbx_sent_frame,
    sent_subject as x_tbx_sent_subject, two_hit_in, two_hit_in2, two_hit_in4, two_hit_sent_in,
    MatrixWalkTable, TwoHit,
};
// Diagnostic (LOSAT_X_STATS): the cells of one extension (tblastx/x_unit_probe.rs).
pub(crate) use two_hit::x_two_hit_cells;

use crate::utils::matrix::blosum62_score;

/// Get the substitution matrix score for two amino acids in NCBISTDAA encoding.
///
/// NCBI aa_ungapped uses the blast core 28x28 score matrix directly; NULLB is
/// not special-cased in the extension loops.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-849
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1579
#[inline(always)]
pub fn get_score(a: u8, b: u8) -> i32 {
    blosum62_score(a, b)
}

/// Convert amino acid coordinates to nucleotide coordinates.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1093-1106 + output adjustment (c++/src/objtools/align_format/tabular.cpp:1037-1058):
/// The C core produces 0-indexed coordinates, then C++ output layer adds +1.
/// We combine both steps here for 1-indexed output.
///
/// For positive frames:
///   Output +1: start = 3*offset + frame - 1 + 1 = 3*offset + shift + 1
///   Output +1: end = 3*end + frame - 2 + 1 = 3*end + shift
///
/// For negative frames:
///   Output +1: start = len - (3*offset + shift + 1) + 1 = len - 3*offset - shift
///   Output +1: end = len - (3*end + shift) + 1
pub fn convert_coords(aa_start: usize, aa_end: usize, frame: i8, dna_len: usize) -> (usize, usize) {
    let f_abs = frame.unsigned_abs() as usize;
    let shift = f_abs - 1;

    if frame > 0 {
        // NCBI core: start = 3*offset + frame - 1
        // NCBI output (c++/src/objtools/align_format/tabular.cpp:1040): +1 for 1-indexed
        // Combined: start = 3*offset + shift + 1
        let start_bp = aa_start * 3 + shift + 1;
        // NCBI core: end = 3*end + frame - 2
        // NCBI output (c++/src/objtools/align_format/tabular.cpp:1041): +1 for 1-indexed
        // Combined: end = 3*end + shift - 1 + 1 = 3*end + shift
        let end_bp = aa_end * 3 + shift;
        (start_bp, end_bp)
    } else {
        // NCBI core: start = len - 3*offset + frame
        // NCBI output (c++/src/objtools/align_format/tabular.cpp:1055): +1 for 1-indexed
        // Combined: start = len - (3*offset + shift + 1) + 1 = len - 3*offset - shift
        let start_bp = dna_len - aa_start * 3 - shift;
        // NCBI core: end = len - 3*end + frame + 1
        // NCBI output (c++/src/objtools/align_format/tabular.cpp:1054): +1 for 1-indexed
        // Combined: end = len - (3*end + shift) + 1
        let end_bp = dna_len - aa_end * 3 - shift + 1;
        (start_bp, end_bp)
    }
}

#[cfg(test)]
mod extension_tests;
