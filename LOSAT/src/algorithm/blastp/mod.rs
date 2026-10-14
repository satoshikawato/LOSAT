//! BLASTP algorithm module
//!
//! This module implements BLASTP (protein vs protein search) in pure Rust,
//! using NCBI BLAST as the behavioral reference.

mod alignment;
pub mod args;
pub mod blast_engine;
// NCBI c++/src/algo/blast/core/blast_setup.c:382-385:
// if (Blast_QueryIsProtein(program_number)) BLAST_ScoreSetAmbigRes(sbp, 'X');
// BLASTP and TBLASTN both encode an untranslated protein query in NCBISTDAA.
pub(crate) mod encoding;
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:217-228
// status = s_BlastAaWordFinder_TwoHit(subject, query, ..., matrix, ...);
// The same protein two-hit extension primitive is used by BLASTP and TBLASTN.
pub(crate) mod extension;
// EXPERIMENT (LOSAT_X_BPUNGAP / LOSAT_X_BPUNGAPSHADOW, off unless set): the ungapped X-drop
// walks of aa_ungapped.c:804-883 with a byte score table, same results (`pub` for the summary line
// of `main`).
pub mod x_ungapped_fast;
// NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// max_offset = BlastGetStartForGappedAlignment(...);
// status = s_BlastProtGappedAlignment(...);
// BLASTP and TBLASTN use the same protein gapped-score path.
pub(crate) mod gapalign;
mod hsp;
mod query_split;
// EXPERIMENT (LOSAT_X_BPQBLOCK): the queries searched in blocks of whole query batches.
mod x_query_blocks;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-3691
// ```c
// static void s_HSPListNormalizeScores(...);
// static void s_HitlistReapContained(...);
// static int s_HitlistEvaluateAndPurge(...);
// static void s_ComputeNumIdentities(...);
// ```
mod kappa;
// EXPERIMENT (LOSAT_X_BPSEED): the subject scan and the two-hit test in one loop.
mod x_seed_fused;
// EXPERIMENT (LOSAT_X_BPPFREE): preliminary hit lists without spare room, each freed after its
// query's redo.
mod x_prelim_free;
// EXPERIMENT (LOSAT_X_BPPCOMPACT): preliminary hit lists with a 64-byte element until the redo.
mod x_prelim_compact;

/// EXPERIMENT (LOSAT_X_BPSEEDSHADOW): the shadow mode's line at exit.
pub fn x_bp_seed_print_summary() {
    x_seed_fused::x_bp_seed_print_summary();
}

pub use args::BlastpArgs;
pub use blast_engine::run;
// EXPERIMENT (LOSAT_X_BPGAP / LOSAT_X_BPGAPSHADOW): the summary lines `main` prints at exit.
// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub use gapalign::x_restricted_gap_print_summary;
// EXPERIMENT (LOSAT_X_BPPCOMPACT=shadow): the summary line `main` prints at exit.
pub use x_prelim_compact::x_bp_prelim_compact_print_summary;
